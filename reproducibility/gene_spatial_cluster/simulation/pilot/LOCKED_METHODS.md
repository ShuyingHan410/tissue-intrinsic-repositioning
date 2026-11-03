# FROZEN COMPARATOR SPECIFICATIONS — do not change except to fix a coding bug
Locked 2026-08-28, after Check 0 passed. Any change invalidates gate2.rds and every
downstream benchmark.

## Environment
- R 4.6.1 (2026-06-24) "Happy Hop", aarch64-apple-darwin23
- Bioconductor 3.23
- DESpace 2.4.0 · edgeR 4.10.3 · mgcv 1.9.4 · metafor 5.0.1   (local pilot, k = 5 results)
- FASRC v16 runs (from 2026-09-13): conda env `rsim`, R 4.5.3, DESpace 2.2.2 (dsp_test present),
  edgeR 4.8.2, metafor 5.0.1. DESpace is permitted as any 2.x (the multi-sample dsp_test interface);
  1.x is NOT acceptable. Every batch log prints the versions it ran with.

## 1. Proposed method (`m_proposed`)
delta_s = sum(n_i r_i)/sum(n_i) - mean(r_i)   per slice
SE0_s   = sqrt( var(r) * (1/N_eff - 1/N) ),  N_eff = (sum n_i)^2 / sum(n_i^2)
gate    : SE0 < 0.05  AND  cells >= 30
          [FROZEN 2026-09-13, plan v16 §0.6 option (a): this is now ALSO the headline real-data gate.
           delta_core.R and inclusion_criteria.R edited the same day; the molecule/rate floor
           (MIN_MOL | MIN_RATE_PER_1K) no longer gates. Real D1 calls at q<0.05: 227 -> 199;
           attrition table in simulation/v16/gate_dryrun/. Pre-gate outputs backed up in
           gene_spatial_cluster/pre_SEgate_backup_20260913/.]
pool    : metafor::rma(yi=delta, sei=se, method="REML", test="knha")   per arm
contrast: D = mu9 - mu0 ; SE = sqrt(se0^2+se9^2) ; Welch-Satterthwaite df
p-value : two-sided t on that df
lambda  : max(1, median(z_blank^2)/qchisq(.5,1)), n<20 usable blanks => lambda=1

## 2. DESpace2 — REAL PACKAGE (`m_despace2_panel`)
DESpace::dsp_test(sce, cluster_col="dom", sample_col="slice", condition_col="C",
                  test="QLF", verbose=FALSE)          # QLF = package default
- pseudobulk sample x cluster, NB via edgeR, cluster x condition interaction
- K = 5 equal-width domains on r
- RUN ON A GENE PANEL (50 genes): edgeR needs cross-gene dispersion.
  Per-gene testing is NOT how the method works and inflates Type I.
- p-value extraction: returned data.frame column `PValue`, matched on `gene_id`.

## 3. GAM — hierarchical factor-smooth (`m_gam(., "fs")`)  [the real comparator]
Co <- as.ordered(C)          # MUST be ordered
bam(Y ~ C + s(r,k=8) + s(r,by=Co,k=8) + s(r,slice,bs="fs",k=5) + s(mouse,bs="re"),
    family=nb(), method="fREML", discrete=TRUE)
p-value extraction: summary(m)$s.table, rows matching ":Co" (the DIFFERENCE smooth).
  With an UNORDERED factor, s(r,by=C) fits one smooth per level, so the test asks
  "is there a gradient?" not "did it change?" -> 93.5% Type I under Null 2.

## 4. GAM — naive cell-level (`m_gam(., "naive")`)  [NEGATIVE CONTROL ONLY]
bam(Y ~ C + s(r,k=8) + s(r,by=Co,k=8), family=nb(), ...)
Appears in Type I panels only. EXCLUDED from power panels: power at an
uncontrolled false-positive rate is not power, and its matched-threshold repair
requires a known null, unavailable on real data.

## 5. Pseudobulk DE (`m_pseudobulk_de`)
glm(X_s ~ offset(log(ncell_s)) + C, family=quasipoisson); Wald p on C.
For DE-specificity controls only.

## Substrate rules
- Pseudo-arm substrate REDRAWN EVERY REPLICATE. A substrate fixed across replicates
  freezes that draw's domain imbalance (realised max share diff 0.002-0.036) and any
  domain-sensitive test rejects it consistently.
- Randomisation unit is the MOUSE.
- Background panel genes are redrawn every replicate. Do NOT cache one panel across
  replicates: it removes a real source of Monte Carlo variability and understates
  uncertainty. Reuse within a replicate across methods is fine.
- Background panel genes are NULL (theta = 0, no DE) — one gene changes, backdrop
  stable. [Bug fix 2026-08-28: run_b1.R gave backgrounds the target's theta_c, so
  in power cells the whole panel shifted together and edgeR normalisation absorbed
  the common interaction, suppressing DESpace2 power (12.5% vs 79.0% on the
  identical cell). B1 null cells unaffected (theta = 0). DESpace2 power re-run
  with null backdrop: b1_ds2fix_results.rds.]
