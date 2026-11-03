# FitCurve - final MorphoGAM curve/coordinate fitter

`curve_finder.R` is the single, self-contained implementation. `morphogam_final.html`
demonstrates it (how it works + audit table + gallery vs the manual plots).

## The method (thin wrapper around MorphoGAM's own CurveFinder)
1. **Connectivity guard (auto-k):** raise kNN k until the graph is one component (bridge
   leftovers) - the one fix that lets MorphoGAM run automatically on every slice
   (default knn=5 disconnected on sparse SM/ME).
2. **Open vs closed = MorphoGAM's native `loop="auto"`** (Mapper-graph cycle test, PC1
   lens); fed the connected graph it decides and fits in one step. Not reinvented.
3. **H1 persistence (ripserr) cross-check:** `agree=FALSE` flags where it disagrees with
   loop="auto" - the genuinely ambiguous (thick/fat/fragmented) shapes to eyeball. loop
   mostly affects `t`, not the `r_for_jsd` the benchmark uses.
4. **Curve smoothness = REML.** MorphoGAM gives the coordinate t (from its geodesic
   embedding); we refit the curve f(t) with `mgcv` splines using **REML + select=TRUE**,
   which auto-selects wiggliness and is more stable than MorphoGAM's default GCV (kills the
   pathological spikes). r = signed orthogonal distance to this curve.
5. **Fallback (tries BOTH modes):** if the REML curve still spikes or CurveFinder errors,
   fit a robust binned centerline in BOTH closed and open mode and keep the lower-spike one.
   This auto-resolves ambiguous fat rings (e.g. petal MU: detectors call it open, but it
   fits cleanly only as closed).

## Parameters - how knn and knot.fx are determined (formal)
- **knn** (the geodesic graph -> the coordinate t): use the **smallest k that connects**.
  A residual-variance sweep shows a sharp short-circuit jump just above the connectivity
  threshold (e.g. res.var 0.007 at knn<=12 -> 0.25 at knn>=20), so smallest-connecting-k is
  the *formally optimal* choice - most faithful geodesics that still connect. Raising knn
  smooths t but short-circuits the shape; lowering it disconnects.
- **knot.fx** (curve wiggliness): replaced by **REML** auto-selection (above). Note there is
  no fully-automatic "correct" smoothness - on a noisy ruffled band the medial curve
  genuinely wiggles, so going smoother than REML means deliberately trading away petal
  detail (a scale choice). REML is the stable default; lower a k cap only to force smoother.

## Use
```r
source("curve_finder.R")
f  <- FitCurve(as.matrix(cells[, c("x","y")]))         # attrs: mode, method, h1, agree, quality, spike
xyt <- fit_layer_xyt(slice_df, "MU")                   # three-step-ready: r_for_jsd[0,1] + metadata
# override closed/open for a hard/ambiguous layer (your eye wins):
xyt <- fit_layer_xyt(slice_df, "ME", force_loop = FALSE)   # force OPEN
```
`fit_layer_xyt` output is a drop-in for the `*_xyt_mapped.csv` the three-step engine reads.

## Honest limitations (hard shapes)
On CLEAR shapes the fit is good. On messy colon shapes the result is inherently limited and
no knn/knot/detector setting fixes it - the 1-D-curve assumption itself is weak there:
- **diffuse blobs** (e.g. 062921 MU): cells fill the interior -> a 1-D curve can't "follow
  petals"; r_for_jsd is approximate. Flagged `quality="diffuse"`.
- **fat rings** (e.g. petal MU): detectors call them open; fit cleanly only as closed
  (auto-resolved by the both-modes fallback).
- **fragmented tissue** (e.g. 082421 ME): scattered arcs -> no clean 1-D fit either way.
Review the flagged layers (`agree=FALSE`, `quality="diffuse"`, high `spike`) and set
`force_loop` where your eye disagrees. Since the three-step uses r_for_jsd DISTRIBUTIONS
(JSD between cell types), imperfect fits may still separate types adequately - Phase 3
measures whether these limitations actually matter downstream.

## Superseded (deleted) approaches
Harmonic ring fit + seal/bridge, binned-primary, geometric closed/open detectors, and a
phi1-lens Mapper were all tested and dropped (see git / PLAN.md history). Key lessons:
Fourier curves spike in sparse gaps; geometric detectors fail on non-convex colon shapes;
MorphoGAM's own loop="auto" + our connectivity guard is the reliable, native choice.

## Next: Phase 3
Wire `fit_layer_xyt` into `three_layer.Rmd` to regenerate `output_threelayer/
per_slice_morphogam/*_xyt_mapped.csv` (into a copy), then re-run the three-step benchmark.
