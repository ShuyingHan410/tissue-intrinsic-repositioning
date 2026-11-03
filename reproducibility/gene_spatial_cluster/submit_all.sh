#!/bin/bash
# One-shot submission with dependencies:  bash submit_all.sh
set -e
mkdir -p logs results chunks
# refuse to submit if inputs are missing / coords not harmonised
bash 99_check_paths.sh | tee logs/precheck.txt
grep -q "ALL CHECKS PASSED" logs/precheck.txt || { echo "precheck failed - not submitting"; exit 1; }
j0=$(sbatch --parsable 00_extract.sbatch)
# NOTE: 01_sweep.sbatch requests --array=1-33 for 990 columns at CHUNK_SIZE=30.
# If you change CHUNK_SIZE, update the array range: expected chunks = ceil(990 / CHUNK_SIZE)
j1=$(sbatch --parsable --dependency=afterok:$j0 01_sweep.sbatch)
sbatch --dependency=afterok:$j1 02_merge.sbatch
echo "submitted: extract=$j0  sweep(array)=$j1  merge queued after sweep"
