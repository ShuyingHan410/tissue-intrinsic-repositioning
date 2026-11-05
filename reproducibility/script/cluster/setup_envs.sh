#!/bin/bash
#SBATCH -J setup_envs
#SBATCH -p shared
#SBATCH -c 4
#SBATCH --mem 16G
#SBATCH -t 0-03:00
#SBATCH -o setup_envs_%j.out
#SBATCH -e setup_envs_%j.err
#
# One-time conda env setup on FASRC. Two ways to run:
#   sbatch setup_envs.sh                                   # batch (uses the proxy)
#   salloc -p test -c 4 --mem 16G -t 0-02:00 ; bash setup_envs.sh   # interactive
#
# Builds the envs (STalign needs old numpy, so it can't share with moscot/paste2):
#   otenv       moscot + PASTE2 + scanpy   (CPU)
#   stalignenv  STalign                    (CPU, numpy 1.23)
set -e
export http_proxy=http://rcproxy.rc.fas.harvard.edu:3128
export https_proxy=$http_proxy
# initialize conda + mamba in a non-interactive (batch) shell
# (~/.bashrc returns early for non-interactive shells, so mamba's function is never
#  defined -> init conda's hook explicitly, and fall back to conda if no mamba.)
CB="$(conda info --base 2>/dev/null)"; [ -d "$CB" ] || CB="${CONDA_BASE:?set CONDA_BASE to your conda installation}"
source "$CB/etc/profile.d/conda.sh"
[ -f "$CB/etc/profile.d/mamba.sh" ] && source "$CB/etc/profile.d/mamba.sh"
command -v mamba >/dev/null 2>&1 || mamba() { conda "$@"; }

# FASRC's proxy 403s repo.anaconda.com (the Anaconda 'defaults' channel). Switch to
# conda-forge, which the proxy allows, and silence the notices/ToS fetch to defaults.
conda config --set number_channel_notices 0   2>/dev/null || true
conda config --remove channels defaults        2>/dev/null || true
conda config --add channels conda-forge        2>/dev/null || true
conda config --set channel_priority flexible   2>/dev/null || true

# The compute-node proxy 403s github.com, and the comparator packages are github-only
# (not on PyPI). CLONE THEM ON THE LOGIN NODE FIRST (login nodes have internet):
#   unset http_proxy https_proxy; mkdir -p ~/repos
#   git clone https://github.com/raphael-group/paste2.git  ~/repos/paste2
#   git clone https://github.com/JEFworks-Lab/STalign.git  ~/repos/STalign
# Then this script installs them from the local copies (no github access needed).
PASTE2_DIR="$HOME/repos/paste2"
CAST_DIR="$HOME/repos/CAST"
STALIGN_DIR="$HOME/repos/STalign"
[ -d "$PASTE2_DIR" ]  || { echo "ERROR: clone paste2 on the LOGIN node first: git clone https://github.com/raphael-group/paste2.git $PASTE2_DIR"; exit 1; }
[ -d "$CAST_DIR" ]    || { echo "ERROR: clone CAST on the LOGIN node first: git clone https://github.com/wanglab-broad/CAST.git $CAST_DIR"; exit 1; }
[ -d "$STALIGN_DIR" ] || { echo "ERROR: clone STalign on the LOGIN node first: git clone https://github.com/JEFworks-Lab/STalign.git $STALIGN_DIR"; exit 1; }

# create env only if it doesn't already exist (so re-runs resume, not error out)
mkenv() {
  if conda env list | awk '{print $1}' | grep -qx "$1"; then echo "env $1 exists - skip create"
  else mamba create -y -c conda-forge --override-channels -n "$@"; fi
}

# --- otenv: moscot + PASTE2 (CPU) -------------------------------------------
mkenv otenv python=3.11
mamba run -n otenv pip install "moscot[spatial]" scanpy anndata pot h5py
mamba run -n otenv pip install "$PASTE2_DIR"

# --- stalignenv: STalign (CPU, needs numpy 1.23) ----------------------------
mkenv stalignenv python=3.11 "numpy=1.23"
mamba run -n stalignenv pip install "$STALIGN_DIR" scanpy h5py

# --- castenv (GPU) -----------------------------------------------------------
# cu121 wheels run on the 12.x GPU nodes (module cuda/12.4.1). torch and DGL MUST be
# a matched pair - pin torch to a version DGL ships wheels for (DGL lags torch, so
# don't let pip grab the newest torch). If DGL's URL 404s, try the plain
# https://data.dgl.ai/wheels/cu121/repo.html or bump torch-2.x to match a DGL build.
mkenv castenv python=3.10
mamba run -n castenv pip install scanpy anndata h5py scikit-learn
mamba run -n castenv pip install "torch==2.3.1" --index-url https://download.pytorch.org/whl/cu121
mamba run -n castenv pip install dgl -f https://data.dgl.ai/wheels/torch-2.3/cu121/repo.html
mamba run -n castenv pip install "$CAST_DIR"


# --- slatenv: SLAT (scSLAT, CPU is enough for one pair) ----------------------
# scSLAT is on PyPI; its PyG dependencies are fetched by the console script
# install_pyg_dependencies from data.pyg.org (needs the proxy). If that URL is blocked,
# run the two pip lines on the LOGIN node instead (login nodes have internet).
mkenv slatenv python=3.11
mamba run -n slatenv pip install scSLAT scanpy anndata h5py faiss-cpu
mamba run -n slatenv install_pyg_dependencies || echo "install_pyg_dependencies failed: run it on the login node"

echo "envs built: otenv, stalignenv, castenv, slatenv"
