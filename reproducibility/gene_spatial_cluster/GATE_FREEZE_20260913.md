# Per-slice eligibility gate — frozen 2026-09-13

**Rule.** A gene x cell-type x slice observation is eligible iff `cells >= 30` and `0 < SE0 < 0.05`,
where `SE0 = sd(r) * sqrt(sum(n^2)/sum(n)^2 - 1/N)` is the exact permutation-null SE of delta.
Replaces `mol >= 20 | mol >= 8 per 1000 cells`. Source: simulation plan v16 §0.6, option (a).

**Code.** `delta_core.R` (eligibility moved after `se_delta`; `SE_MAX` constant),
`inclusion_criteria.R::coverage_of` (same rule for coverage/routing).

**Data.** `results/slices_*.csv` re-gated in place; original flag kept as `eligible_old`.
`03_meta.R` rerun -> `arms.csv`, `contrasts.csv`, `blank_diagnostics.csv` (log: `03_meta_SEgate.log`).
Pre-gate copies of all of these, plus the edited scripts' previous versions, in `pre_SEgate_backup_20260913/`.

**Effect.** D1 calls at q < 0.05: 227 -> 199 (195 kept, 32 lost — all had q in (0.01, 0.05) — 4 gained).
Median k: D0 parent 10 -> 10; D9 daughter 20 -> 17; D9 parent 9 -> 7. Blank lambda essentially unchanged.
Full attrition tables: `simulation/v16/gate_dryrun/`.

**Outstanding.** `gene_spatial_lineage_results.csv` (routing) was produced by `01_run_chunk.R` under the old
coverage rule. Rerun the chunk array on the cluster with the edited `inclusion_criteria.R` to regenerate
it; the dry run predicts 0 routing changes.
