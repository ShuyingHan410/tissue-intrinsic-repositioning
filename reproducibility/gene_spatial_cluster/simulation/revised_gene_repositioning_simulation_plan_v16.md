# Revised Simulation Plan for Benchmarking Directional Gene Repositioning

*Version 16 changes from v15 — a restructure, not a redesign, prompted by a second external review. (1) **§0 is rewritten as the single operative protocol**; §0a–§0c and the v12 core-benchmark text (now §0-hist) are the change log and must not be used for new runs. (2) **Experiment 1 runs first and alone**; a comparator non-nominal at a τ_θ level is reported only at its matched threshold there. (3) **Null 1 at the new reference τ_θ = 0.05 (size 1) is calibrated**: 4.44% [4.01, 4.87], n_eff 9,014, realised I² 15.2% — all four prespecified levels now pass. (4) **Experiment 0 completed as specified**: all four quantities with p10/p90 at every stratum, size 1 added (`v16/ppc_neff.csv`, `ppc_neff.png`). (5) τ_θ = 0.05 propagated: Experiments 7 and 8 take it as primary; Experiment 3 stays at 0.15 with the reason stated. (6) τ_θ is the *imposed* parameter; on real-geometry substrates total heterogeneity is summarised only by realised I². (7) The bimodal shape is written as the one-parameter family the engine implements. The gate decision (§0.6) was resolved the same day: option (a) frozen, real pipeline re-gated, 227 → 199 calls. **In-place edits after a third review (same day):** matched thresholds are design-specific \(c_{m,k,\tau}\) and a comparator Null 1 at k = 20 (Experiment 1b) is added; the FDR repair is pre-specified as BH on empirically calibrated p-values \(\tilde p=\hat F_0(p)\) alongside native BH; method roles are stated once in §0.3; the Experiment 3 rationale is reworded as a design intent, not a predicted result.*

*Version 15 changes from v14, written after the local calibration stage (`simulation/v16/CALIBRATION_STAGE.md`) and before any cluster run (see §0c). (1) **k = 10 confirmed calibrated**: Type I 4.4–4.7% at every τ_θ on ~10,000 effective null tests per cell; the earlier 8% was noise. (2) **NB size = 1** for all new runs: real Fibro 2 genes have N_eff/N_event of 0.944 / 0.889 / 0.611 at 27 / 64 / 331 molecules, size 1 matches within 1% at every stratum, size 2 overstates by 4–10%. (3) **τ_θ = 0.05 is the real-median heterogeneity** (realised I² 14% at k = 10; the estimator floor at k = 10 is 0, not 11%), so Experiments 1–2 gain a fourth τ_θ level and the reference becomes 0.05. (4) **σ_mouse ≈ 0** for the typical real gene (three-level REML on the Day-9 daughter arms, 8 mice; 75th percentile 0.04–0.07, 90th indistinguishable from blanks), so Experiment 9 runs σ_mouse ∈ {0.07, 0.20} and the reference is 0. (5) Engine extensions for Experiments 6b, 9, 11 smoke-tested. One item remains open: the gate decision (§0b B1).*

*Version 14 changes from v13, written the same day and before any v13 run (see §0b), prompted by an external review and two same-day measurements. (1) **Gate-induced selection is part of the estimand**: every cell reports P(testable), k₀, k₉, power | testable and P(correct discovery), and the benchmark and the headline real analysis must use the same per-slice gate — the real main branch currently uses a transcript/rate floor while the simulation and the stratified branch use SE₀ < 0.05. (2) **Generator fidelity check** (Experiment 0): real versus simulated N_event, N_eff, N_eff/N_event and max-share per abundance stratum, with a pass criterion, before the abundance experiments run. (3) **Substrate realism stated and tested**: the ten real Day-0 slices differ in spread enough to inject τ_θ ≈ 0.06 of geometric heterogeneity into an effect gene, which pooled resampling erases; a within-slice bootstrap cell (Type I 6.5% vs pooled 4.0%; power 97.5% vs 97.0%) is added to Experiment 5, and pseudo-mice are described as grouping labels. (4) Four robustness cells: awkward directional shapes (Exp 11), heavy-tailed and single-outlier mouse effects (in Exp 9), FDR under correlated gene blocks (Exp 6b), and the Fibro 6 substrate promoted to Exp 10. (5) Required QQ diagnostic for λ down to p = 10⁻³. (6) A wording rule scoping every comparative claim to directional first-moment repositioning on a matched substrate at a matched operating point. No method, comparator specification, or completed result is altered.*

*Version 13 changes from v12, all recorded before any run at the new design (see §0a): (1) the **primary matched pseudo-arm design becomes 10 slices per arm** (4 pseudo-mice, 2/3/3/2), which is the real Day-0 Fibro 2 design; the v12 value of 5 was inherited from the v9 literal split and never matched the real parent arm, so the v12 main experiments were harsher than the analysis they validate. The completed k = 5 results are retained verbatim as a labelled half-design stress test; a k = 20 secondary design matches the daughter arm. (2) **Null cells score the 50-gene background panel**: the v12 scripts scored one gene per replicate (200 null tests per cell, not the 2,000 specified), so the calibration band was ±3 points. (3) **§8 abundance strata are run** (M = 27, 331), plus a per-arm abundance scenario (M_B/M_A ∈ {1/3, 1, 3}) that replaces amplified/dampened and reflects the real Day 0 → Day 9 abundance change. (4) **σ_mouse is estimated from the Day-9 arms** by three-level REML rather than assumed at 0.20, and a new Experiment 9 runs the proposed method with λ on — the only cell that tests the pipeline as deployed; v12's main results are the core without λ at σ_mouse = 0. (5) Heterogeneity is labelled by **realised I²**, and τ_θ = 0 is identified as the typical real case (real median I² ≈ 14%), so the v12 reference 'I² = 15%' is met by τ_θ = 0, not 0.15. (6) The **real pipeline adopts the simulation's SE₀ < 0.05 eligibility gate**; the rate-floor alternative admitted slices with 1–2 molecules. Nothing in the method, the comparator specifications or the completed results is altered.*

*Version 12. Restructures the plan around a **core benchmark of six experiments** (new §0) —
one clean experiment per scientific claim, varying one factor at a time around a reference
setting of M=64, I²=15%, σ_mouse=0.20 — with everything from §1 onward retained as the master
validation checklist rather than a run list. Cuts the primary benchmark from ~12–17 core-hours to
~35 min plus a ~1 h FDR run. Three refinements over a pure one-factor-at-a-time design, all
forced by pilot evidence: (1) **τ_θ is crossed with method**, not varied separately, because the
pilot measured a ranking reversal — DESpace2-style overtakes the proposed method at τ=0.30
(0.299 vs 0.276) while losing at τ=0 and 0.15; (2) **σ_mouse is crossed with Ablation C**, which
has nothing to show at σ_mouse=0; (3) **pattern-classification accuracy is retained** in reduced
form, since demoting lost/amplified/dampened would otherwise remove the only metric no comparator
produces. Records two preconditions the restructure does NOT remove: comparator check 0 currently
fails for all three comparators, and the reference σ_mouse=0.20 has never been run.*

*Version 11 changes from v10, both driven by the pilot run: (1) new §1b subsection on
**implementation traps** — mgcv's `by=` factor must be ORDERED or `s(r, by = C)` fits one smooth
per level instead of the difference, which gave the GAM comparators **93.5% Type I error under
Null 2** in the pilot; generalised to the rule that **Null 2 is a comparator-validity check**
and any method rejecting there is testing "is there a gradient?" rather than "did it change?" —
a failure that *inflates* headline power rather than depressing it; (2) the no-random-effects
GAM is **retained as an explicit negative control** rather than deleted or treated as a
competitor: it appears in Panel B (6.0% at \(\tau_\theta=0\) versus **43.5%** at
\(\tau_\theta=0.30\), against 9.0% for the factor-smooth version) and is **excluded from the
power panels**, because power at an uncontrolled false-positive rate is not power and the
matched-threshold repair is a simulation-only device unavailable on real data.*

*Version 10 changes from v9: (1) the "fit the profile model per slice and meta-analyze the
interaction" option is **deleted as non-identifiable** — condition is constant within a slice,
so \(C_i s_\Delta(r_i)\) cannot be estimated there; the GAM comparator is now a single
hierarchical model with a factor-smooth and mouse effect, kept only if it passes the gate;
(2) Benchmark B **split into B1 (matched-geometry, primary) and B2 (real-geometry robustness)**,
because a 5% threshold learned on the matched null does not transfer to mismatched geometry —
the very fact the truth audit establishes; (3) all 500 Benchmark C genes pinned to one
substrate; (4) the pseudo-arm randomization unit fixed as the **mouse**, with two committed
constructions (synthetic pseudo-mice for B1, literal 2v2 split as realism check), replacing the
self-contradictory "assign slices respecting mouse nesting"; (5) §11's "shared null" sentence
qualified to matched-substrate scenarios; (6) DE-only expected results split by variant.*

*Version 9 changes from v8: (1) **shared nulls now require matched cell populations.** Under
Null 1 a completely flat gene (\(h\equiv1\)) already violates \(H_{0,\mathrm{DESpace2}}\) in
Simulation A, because Fibro 2 (D0) and IAF 4 (D9) occupy different parts of the axis — the
expected domain share equals the cell share, and the real coordinates differ by up to 6.4
percentage points per domain (\(\chi^2\ p=2.4\times10^{-97}\)). \(\delta\) is immune since it
subtracts \(\bar r_{\text{cells}}\), so \(\delta\) and the profile methods **cannot** share a
null on differing populations. The §1b calibration gate, Benchmark C's 400 common nulls and the
matched-geometry DE-only variant are therefore respecified onto a **matched pseudo-arm
substrate** — Day-0 slices randomly split into two pseudo-conditions, preserving slice as the
replicate unit and respecting mouse nesting. Scenarios that deliberately keep the real Day 0 /
Day 9 populations are unchanged but are scored per method via the truth audit. Each reported
number must record which substrate produced it.*

*Version 8 changes from v7: (1) DESpace2 now specified as it actually works — pseudobulk
\(X_{gsd}=\sum_{i\in(s,d)}Y_{gsi}\) within sample × domain, then the multi-sample DSP test with
slice as the replicate unit; (2) new **comparator-truth audit**: solving \(\theta\) per arm (the
v4 fix) makes \(h_0\ne h_9\) whenever \(SD(r_0)\ne SD(r_9)\), so **Null 2 and DE-only are not
automatically shared nulls** — in Fibro 6 → IAF 3 the \(h_9/h_0\) ratio swings **41%** across
the axis while \(D_\delta=0\); nulls are now classified per method per scenario before
simulating, and **DE-only is split** into realistic-geometry and matched-geometry
(\(h_0\equiv h_9\)) variants, the latter being the shared null; (3) Ablation B fully specified
as fixed-effect inverse-variance pooling with \(\tau^2=0\), replacing "such as", since
cell-pooling would confound replicate unit, weighting and uncertainty at once; (4) C2's quantile
convention pinned to the order statistic \(k=\lceil0.95(B+1)\rceil\).*

*Version 7 changes from v6: (1) new §1b **comparator fairness precondition** — the §1
cell-level GAM has no slice/mouse term and rejects a true null 54% of the time under
slice-varying gradients (a random intercept makes it 61%, intercept+slope 47%, while the
slice-level \(\delta\) contrast is 5%); comparators must match the replicate structure, and
Panel C now compares power at **matched empirical Type I error**; (2) DE-only added to the
method-specific null table as **not null for pseudobulk DE**; (3) \(\sigma_{\text{mouse}}\)
grid specified {0, 0.10, 0.20, 0.40} and the residual "matched to the observed \(\lambda\)"
circularity removed; (4) abundance strata quantified as 27 / 64 / 331 molecules per slice
(25th/median/90th, post-gate); (5) \(I^2\) targets quantified as {0, 15, 45, 60}% against an
observed median of ~14%; (6) blank median corrected from ~10 to **16** molecules; (7) Simulation
B sizes fixed and a new §12b reproducibility section (per-replicate seeds, achieved-truth
storage, comparator versions).*

*Version 6 changes from v5: (1) C1 now uses the **floored** estimator
\(\hat\lambda_{\text{applied}}=\max(1,\hat\lambda_{\text{raw}})\) with the real pipeline's
\(n<20\Rightarrow\lambda=1\) fallback — v5 benchmarked an unfloored variant that would apply
\(\hat\lambda<1\) in ~51% of null runs and *inflate* \(z\); the floor fires ~50% of the time and
biases \(E[\hat\lambda_{\text{applied}}]\) to ~1.07 at 205 blanks, close to the observed D1
value of 1.08; (2) Benchmark A reframed as **null for the population-relative target** — Null 3
is not null for GAM/DESpace2, whose rejection rates there must be labelled rejection rate, not
Type I error; (3) primary Benchmark C restricted to a **directional-only** mixture
(400 nulls + 100 directional alternatives) so a single \(\pi_0\) is well defined for all
methods, with the non-directional scenario kept in Benchmark B; (4) the BH family specified as
one lineage x layer x route family of 500; (5) amplified/dampened **anchored** at
\(\delta=0.025\) so \(b-a=|D|\) no longer admits infinitely many designs; (6) the pilot's
"simulator is at fault" attribution softened to a diagnostic trigger.*

*Version 5 changes from v4: (1) the benchmark is split into **three sized benchmarks** —
calibration, comparative power, mixed-gene FDR — with \(G\), \(\pi_0=0.8\), \(R\) and a
cluster-robust Monte Carlo SE now specified, since FDR was previously undefined and \(R=100\)
could not separate 5% from 8.8%; (2) Ablation C gains a **third arm**: median-\(\lambda\)
repairs tails only under pure scale inflation \(Z\sim N(0,\lambda)\), and the real blanks look
instead like a tail-shape problem, so a direct tail calibration \(c_{0.05}=Q_{0.95}(|Z_{\text{blank}}|)\)
is added as a *diagnostic* — deliberately not folded into the primary method label;
(3) blank calibration now targets **usable statistics (~200 per statistic), not 50 barcodes**,
because \(\hat\lambda\) from 50 blanks spans [0.50, 1.85] by chance, plus a \(B_{\text{usable}}\)
sensitivity sweep; (4) Ablation A gains an explicit **Dirichlet-multinomial concentration
parameter** so the \(N_{\mathrm{eff}}\) contrast is constructed rather than hoped for;
(5) mouse-nesting hedge removed; (6) pilot resized per arm with explicit acceptance checks.*

*Version 4 changes from v3: (1) simulated blanks now carry a zero-mean **mouse effect** in the
correlated regime — with \(h=1\) they stay independent across slices and \(\hat\lambda\) would
read ~1 even when gene tests are inflated, which would make Ablation C measure nothing;
(2) Null 3 resampling upgraded to **moment-balanced** two-parameter tilting, since a
one-parameter tilt controls the mean only; (3) §10 no longer claims mouse nesting **causes**
the observed \(\lambda\), and no longer treats reproducing \(\lambda\) as independent
evidence; (4) DE-only truth defined on the \(\delta\) scale with \(\theta\) solved per arm —
the same \(\theta\) gives an 18% spurious \(D\) in the Fibro6→IAF3 geometry; (5) new §7
table mapping \(a\) to \(|D|\) per pattern — `reversed` needs \(a=0.025\), not 0.05;
(6) \(I^2\) wording corrected.*

*Version 3 changes from v2: (1) Ablation C now specifies how \(\lambda\) is obtained — a
simulated 50-blank set, so the ablation tests the estimation procedure and not just the
arithmetic; (2) Null 3 resampling switched from truncation to exponential tilting, because
truncation shrinks \(\operatorname{var}(r)\) and therefore the positional SE, confounding the
comparison; (3) Null 3's expected-results block now separates a genuine false positive from
the GAM/DESpace2 correctly answering a different question; (4) the \(I^2\) heterogeneity
target is now tied to a named arm.*

---

# 0. THE CORE BENCHMARK — OPERATIVE PROTOCOL (v16)

**This section is the only run specification.** §0a–§0c are the change log from v13–v15 and
the v12 core-benchmark text is kept after them under "§0-hist"; all of that is retained for
auditability and pre-registration. **Superseded settings in those sections must not be used for
new runs.** Where §0 and any later section disagree, §0 governs.

## 0.1 The four claims

1. the method is **calibrated** under a true null;
2. it has **power** for directional repositioning;
3. it measures something **different** from ordinary DE and from general spatial-profile change;
4. its distinctive components — \(N_{\mathrm{eff}}\), random effects, calibration — **matter**.

## 0.2 Reference setting

\[
\boxed{k=10\ \text{slices per arm (4 pseudo-mice, 2/3/3/2)},\qquad M=64,\qquad \text{NB size}=1,\qquad
\tau_\theta=0.05\ (\text{realised } I^2\approx14\%),\qquad \sigma_{\text{mouse}}=0,\qquad
\text{matched pseudo-arm substrate, redrawn per replicate},\qquad R=200}
\]

Each value was fixed by a measurement, recorded in §0c: k = 10 is the real Day-0 Fibro 2 design
(C1, calibrated at 4.4–4.7%); size 1 matches real \(N_{\mathrm{eff}}/N_{\text{event}}\) within
1% at M = 27 / 64 / 331 (C2); \(\tau_\theta = 0.05\) realises the real median \(I^2\) (C3);
\(\sigma_{\text{mouse}}\) is indistinguishable from 0 for the typical real gene (C4).
Stress levels are crossed where the experiment's question needs them and are labelled:
\(\tau_\theta = 0.15\) "upper quartile (\(I^2\approx54\%\))", \(0.30\) "beyond real 90th pct
(\(\approx84\%\))"; \(\sigma_{\text{mouse}} = 0.07\) "real upper quartile", \(0.20\) "stress".

**Heterogeneity language.** \(\tau_\theta\) is the *imposed* heterogeneity parameter. On the
pooled pseudo-arm substrate it is also the total, because pooling removes the slices' own
geometric differences. On real-geometry substrates (Experiment 5) total heterogeneity includes
geometry (≈ \(\tau_\theta\) 0.06 equivalent for an effect gene, §0b B3) and is summarised
**only by realised \(I^2\)**. Every figure reports realised \(I^2\) beside \(\tau_\theta\).

## 0.3 Methods

Proposed: \(\delta_s\), \(SE_0\) from \(N_{\mathrm{eff}}\), per-slice gate (§0.6), REML +
Knapp–Hartung per arm, D1 t-contrast, floored median \(\hat\lambda\) where blanks are simulated.
Comparators at frozen defaults (`LOCKED_METHODS.md`): hierarchical GAM-fs (slice factor-smooth +
mouse RE), DESpace2 on a 50-gene panel with null backdrop, naive GAM (negative control only),
pseudobulk DE (non-spatial control). Every comparative sentence follows §0b B8.

**Method roles, stated once.** *Primary spatial comparison:* proposed vs GAM-fs vs DESpace2 — these three
appear in power, matched-power and FDR panels. *Pseudobulk DE:* orthogonal abundance control; appears in
calibration and specificity panels, never in power panels. *Naive GAM:* calibration negative control only;
appears in Type I panels, **excluded from every headline power and FDR comparison** because its
false-positive rate is uncontrolled. "all 5" in the table below means *run*, not *compared on equal
footing*; inferential status is as stated here.

## 0.4 Experiments — the operative table

All at the reference setting unless a cell says otherwise; NB size 1 throughout except the
retained k = 5 results and Ablation 4A (size 2, labelled). Every cell reports
\(P(\text{testable}), k_0, k_9\), power | testable, \(P(\text{correct discovery})\) and realised
\(I^2\) (§0b B1, §0a A5).

| # | experiment | k | scenarios | grid | methods | where |
|---|---|---:|---|---|---|---|
| 0 | Generator fidelity | — | real vs simulated slices | M ∈ {27,64,331}, blanks 16; size ∈ {0.7,1,2}; N_event, N_eff, ratio, max share; median [p10,p90] | — | **done** (`v16/ppc_neff.csv`, `.png`) |
| **1** | **Calibration** | 10 | Null 1, Null 2 | \(\tau_\theta\in\{0,.05,.15,.30\}\); 50 null genes/rep | all 5 | cluster — **runs first, alone** |
| **1b** | **Calibration at k = 20** | 20 | Null 1 | \(\tau_\theta\in\{0,.05,.15,.30\}\); 50 null genes/rep | all 5 | cluster, with 1 — supplies \(c_{m,20,\tau}\) and \(\hat F_{0,m,20,\tau}\) for Experiments 6 and 8 |
| 2 | Directional power | 10 | gained, reversed | \(|D|\in\{.025,.05,.075\}\times\tau_\theta\in\{0,.05,.15,.30\}\) | all 5 | **done** |
| 3 | Specificity | 10 | DE-only, repositioning-only, non-directional | **\(\tau_\theta = 0.15\), deliberately**: specificity is evaluated under upper-quartile heterogeneity to test whether estimand separation persists under a challenging but empirically relevant level of slice-to-slice variation | all 5 + pseudobulk | **done** |
| 4A | \(N_{\mathrm{eff}}\) ablation | 5 | as v12 | as v12 | proposed | keep (size 2) |
| 4B | REML+KH vs fixed | 10 | null, gained | \(\tau_\theta\in\{0, .18\}\) (realised \(I^2\approx0, 60\%\)) | proposed | **done** (`v16/ablB_k10.csv`; null re-run at 50 genes/rep: 4.8% vs 21.9%) |
| 5 | Real geometry | real 10 v 18; bootstrap 5 v 5 | Null 3A/3B, gained, reversed; Null 1, gained | \(|D|=.05\); imposed \(\tau_\theta\in\{0,.15\}\); report realised \(I^2\) | all 5 | **done** |
| 6 | FDR | 10, 20 | 400 null + 100 gained | \(\tau_\theta = .15\) | proposed, GAM-fs, DESpace2, pseudobulk | **done** (k = 10, 20) |
| 6b | FDR, correlated blocks | 10 | as 6; ρ = .7, blocks of 20 | \(\tau_\theta = .15\) | proposed, GAM-fs, DESpace2 | **done** |
| 7 | Abundance | 10 | Null 1, gained {.025,.05} | **primary: \(M\in\{27,64,331\}\) at \(\tau_\theta=0.05\)**; secondary: \(M\in\{27,331\}\times\tau_\theta\in\{0,.15,.30\}\) | all 5 | **done** |
| 7b | Per-arm abundance | 10 + real | gained .05, matching null | \(M_B/M_A\in\{1/3,1,3\}\times\tau_\theta\in\{.05,.15\}\) | all 5 | **done** |
| 8 | Slice count | 10 vs 20 | Null 1, gained .05 | **primary: \(\tau_\theta=0.05\)**; secondary {0,.15,.30} | all 5 | **done** |
| 9 | Deployed pipeline | 10 | Null 1, gained .05 | \(\sigma_{\text{mouse}}\in\{.07,.20\}\times\tau_\theta\in\{.05,.15\}\times\{\text{normal},t_3,\text{outlier}\}\); λ on; QQ of blank z | proposed full, GAM-fs, DESpace2 | **done** |
| 10 | Second lineage | 10 | Null 1, gained .05 | Fibro 6 pool; **M = 43** (the all-panel measurement that gives Fibro 2 = 42 gives Fibro 6 = 28; 64 × 28/42); \(\tau_\theta\in\{.05,.15\}\); **NB size 2** (Exp 0 on Fibro 6 done: size 2 within 1% / 1% / 7% at 27 / 64 / 331, size 1 −19% at 331 — half the cells per slice make over-dispersion bite harder) | proposed, GAM-fs | **done** |
| 11 | Awkward shapes | 10 | localised (π = .155), bimodal (κ = 19.0) | \(|D|=.05\), \(\tau_\theta=.15\) | all 5 | **done** |
| — | Null 1, proposed, k ∈ {5,10,20}, \(\tau_\theta\) ∈ {0,.05,.15,.30} | | | | | **done** (§0c C1) |
| — | \(\tau_\theta\to I^2\) map, k ∈ {5,10} | | | | | **done** (§0c C3) |
| — | \(\sigma_{\text{mouse}}\) estimate, real D9 | | | | | **done** (§0c C4) |
| — | Exp 0 on Fibro 6 | | | | | **done** (`v16/ppc_neff_Fibro6.csv`): size 2 |
| — | Proposed-only previews of Exp 2, 8, 4B at the reference design | | | | | **done** (`v16/power_preview.csv`, `ablB_k10.csv`); previews, not results — superseded by the cluster runs |
| — | k = 5 results, all v12 cells, size 2 | | | | | retained, labelled |

**Shape definitions for Experiment 11** (matches `engine_v16.R`):
localised: \(h_i = (1-\pi)/n + \pi\,w_i/n_w\), \(w_i = \mathbb 1[r_i\in[0.75,0.85]]\),
\(\pi = \delta^\ast/(\bar r_w - \bar r)\) closed-form;
bimodal: \(h(r) = 1 + \tfrac{\kappa}{3}\,\tilde\phi(r;0.25,0.05) + \kappa\,\tilde\phi(r;0.75,0.05)\)
with \(\tilde\phi\) the normal density scaled to peak 1 — **one free parameter**, the 1:3 bump
ratio is fixed, \(\kappa\) solved by `uniroot` on expected \(\delta\). Total mass is preserved by
the \(\mu_i = M h_i/\sum h\) normalisation for every \(h\); it is not a second constraint.

## 0.5 Sequencing and the comparator gate

1. **Experiment 1 runs first and alone.** Inspect Null 1 and Null 2 for every comparator at all
   four \(\tau_\theta\) levels with the 50-gene scoring (≈ 10,000 null tests per cell).
2. **Matched thresholds are design-specific.** For method \(m\) at design \((k, \tau_\theta,
   \text{substrate}, \text{gate})\),
   \[ c_{m,k,\tau} = Q_{0.05}\!\left(p_{m,\text{Null 1}}\mid k,\tau_\theta,\text{substrate},\text{gate}\right), \]
   the 5th percentile of that method's own Null 1 p-values from Experiment 1 (k = 10) or 1b (k = 20)
   at the same \(\tau_\theta\). A threshold learned at one k is never applied at another; the k = 20
   cells of Experiments 6 and 8 use \(c_{m,20,\tau}\). A comparator whose Type I error lies outside
   [3%, 7%] at a design is **non-nominal there**, and its power at that design is reported **only at
   \(c_{m,k,\tau}\)**, never at nominal 0.05. Real-geometry cells (Experiment 5) take their thresholds
   from the null cells run on that substrate, per the truth audit. The historical expectation, from
   the k = 5 runs, is that GAM-fs and DESpace2 are nominal at \(\tau_\theta\le0.15\) and non-nominal
   at 0.30; Experiments 1 and 1b decide.
**Experiment 1 / 1b — DONE 2026-09-14** (`v16/EXP1_RESULTS.md`). Proposed 4.7–5.1% at all eight designs. GAM-fs and DESpace2
nominal at τ_θ ≤ 0.05 at both k; matched at 0.15 and 0.30 (GAM-fs at k = 20, τ = 0.15 sits on the 7.0% boundary — both
thresholds are reported everywhere). Null 2 tracks Null 1 for all methods (no gradient trap). Doubling k does not reduce the
comparators' inflation. Naive GAM 19% / 44%. Batch 2 cleared.

3. **FDR repair, pre-specified.** A single 5% cutoff does not calibrate the whole p-value distribution
   that BH uses. For Experiments 6 and 6b, each comparator's p-values are first transformed by the
   empirical CDF of its own Null 1 p-values at the same design,
   \[ \tilde p = \hat F_{0,m,k,\tau}(p), \]
   estimated from the ≈ 10,000 null tests of Experiment 1 / 1b (independent of the FDR replicates).
   Under a valid null \(\tilde p\sim\mathrm{U}(0,1)\). BH is applied to \(\tilde p\) and realised FDR / TPR
   are reported at \(q\in\{0.05, 0.10\}\) as the **calibrated** result. The **native** BH result (on raw
   \(p\)) is reported beside it, because native miscalibration is itself a finding. The proposed method is
   reported both ways too, so the transform is applied symmetrically. The v14 oracle-matched-FDR curve
   (q\* at realised FDR = 5%) is retained as a third, descriptive view.
4. Only then launch Experiments 2, 3, 6. Experiments 5, 7–11 follow.

**Experiments 2, 3, 6 — DONE 2026-09-14, FROZEN** (`v16/BATCH2_RESULTS.md`). Power at matched operating points: parity among
proposed, GAM-fs, DESpace2 at \(|D|\ge0.05\) at every \(\tau_\theta\); no GAM-fs edge and no \(\tau=0.30\) reversal at k = 10;
DESpace2 behind at \(|D|=0.025\). Specificity table as designed (proposed 6.0% on DE-only, 3.5–5.5% on non-directional; comparators
84–100% on non-directional; pseudobulk 100% on DE-only). FDR: only the proposed method holds BH FDR natively (3.0% / 4.3% at
q = 0.05, TPR 85% / 100% at k = 10 / 20) — **the k = 5 TPR of 5% was the small-df t**. GAM-fs 7.7–7.9% native → 3.5% calibrated;
DESpace2 19–27% native, 12–19% calibrated because its null is panel-dependent (13–19% of the mixture's nulls fall below 0.05
vs 7% in the all-null calibration panel). Experiment 4B at k = 10 done locally (fixed-effect Type I 21.9% vs REML+KH 4.8% at
realised \(I^2\approx60\%\)).

5. Re-freeze `FROZEN_CONCLUSIONS.md` after Experiments 1, 2, 3, 4B, 6.

**Batch 3 (Experiments 5, 7, 7b, 8, 9, 10, 11, 6b) — scripted 2026-09-14** as one 85-cell array (`v16/run_batch3.R`,
`batch3.sbatch`, `summarise_batch3.R`; `LIST=1 Rscript run_batch3.R` prints the cell table). Thresholds: pseudo-arm cells at
M = 64 use Experiment 1's \(c_{m,k,\tau}\); every other design uses its own in-batch Null 1 partner (same substrate, τ_θ, M, k,
mouse setting), so nominal and matched are both reported everywhere. Exp 9 reports the proposed method with and without λ,
the blank \(\hat\lambda\) distribution, floor rate, and blank-z tail quantiles. Exp 7b's null cells are the shared-null check
for per-arm abundance.

**Batch 3 — DONE 2026-09-14, FROZEN** (`v16/BATCH3_RESULTS.md`; 85/85 cells, no failed replicates). Real geometry: proposed 3.0%
under Null 3A, DESpace2 71% (rejection rate, the truth-audit case); bootstrap reproduces the pooled substrate. Abundance: calibrated
at M = 27 and 331; DESpace2's inflation grows with M (19.5% at 331, τ = 0.30); power floored by heterogeneity above 64; a 3× abundance
change in either direction produces no call. k = 20: parity at ceiling. Mouse correlation (σ = 0.07, real upper quartile, 10,000-test
recheck): proposed 7.9% / 7.4% uncorrected → **5.0% / 4.9% with λ**; at σ = 0.20 all methods 10–50% uncorrected, λ over-corrects to
0–3% because blanks share the genes' mouse draw (adaptive-conservative; power 77% → 60%). Fibro 6 (size 2, M = 43): 5.2% / 4.9% on
10,000 tests, power 100% / 98.5%. Shapes: 99.5% / 98.5%. Correlated-gene FDR: proposed 3.6% native, GAM-fs 10%, DESpace2 17.7%.
Pseudobulk is uninformative on the real substrate (dispersion ≈ 130 from 33–2,567-cell slices) and is not scored there.

**All cluster experiments complete.** Report regenerated 2026-09-14: `v16/publication_simulation_report_v16.Rmd` → `.html`
(every number computed at render time from the frozen result files; every figure states k, size, substrate, gate and realised I²;
every comparative sentence scoped per §0b B8). **The v16 programme is closed.** One pipeline loose end outside the simulation:
the routing file needs a chunk-array rerun under the frozen gate (predicted no-op).

## 0.6 The per-slice gate — DECIDED 2026-09-13: option (a), ≥ 30 cells and \(SE_0 < 0.05\)

The benchmark and the headline real analysis must use the same gate (§0b B1). Options:
(a) the real main branch adopts ≥ 30 cells **and** \(SE_0 < 0.05\), replacing the
molecule/rate floor; (b) the engine adopts the real floor. **Recommended: (a)** — it measures
positional information directly, rejects 1–2-molecule slices automatically, removes the rate
constant, and is already the primary rule of the stratified branch. On choosing, re-run the real
eligibility branch and `03_meta.R`, produce the attrition table (rows changing eligibility;
\(k\) before/after; D1 calls changing at q = 0.05), record the gate in `LOCKED_METHODS.md` §1.
**The gate is recorded; cluster scripts use \(SE_0 < 0.05\) and ≥ 30 cells, fixed.**

**Dry run of option (a) on the real data (2026-09-13, `v16/05_gate_dryrun.R`, read-only copy in
`v16/gate_dryrun/`).** Eligibility recomputed as ≥ 30 cells and \(0 < SE_0 < 0.05\); unmodified `03_meta.R`
rerun on the copy.

| arm | rows dropped | rows newly admitted | median k, old → new |
|---|---:|---:|---|
| Day 0 parent | 310 / 11,629 | 71 | 10 → 10 |
| Day 9 daughter | 3,657 / 22,953 | 250 | 20 → 17 |
| Day 9 parent | 4,079 / 11,256 | 0 | 9 → 7 |

Dropped rows have median 4 molecules in ~100 cells; a few dozen rows with 40–143 molecules also fail, all from
two small Day-9 slices whose positional spread pushes \(SE_0\) just over 0.05. D1 calls at q < 0.05:
**227 → 199** (195 kept, 32 lost, 4 gained, 0 routing changes). All 32 lost calls had q ∈ (0.01, 0.05) under
the current gate and q ∈ (0.05, 0.10) under the SE gate; none of the 141 calls with q < 0.01 is lost. Option (a)
removes the marginal tier and nothing else. Changed genes: `gate_dryrun/attrition_calls_changed.csv`.
**Decided and frozen 2026-09-13 (option a).** `delta_core.R` and `inclusion_criteria.R` now gate on
≥ 30 cells and \(0 < SE_0 < 0.05\); the per-slice result files were re-gated in place (column
`eligible_old` retained) and `03_meta.R` rerun: 199 D1 calls at q < 0.05. Pre-gate outputs are in
`gene_spatial_cluster/pre_SEgate_backup_20260913/`. Blank \(\hat\lambda\) is essentially unchanged by the
gate (D1 1.26 → 1.29; parent-D0 arm 1.21 → 1.53; daughter 1.60 → 1.42), so the loss of 28 calls is the
gate, not the correction. Recorded in `LOCKED_METHODS.md` §1. *Caveat:* gene routing
(`gene_spatial_lineage_results.csv`) still reflects the old coverage rule until `01_run_chunk.R` is rerun
on the cluster with the edited `inclusion_criteria.R`; the dry run showed 0 routing changes, so this is
expected to be a no-op, but it must be confirmed on that rerun.

## 0.7 Status of the v12 preconditions

| v12 blocker | status |
|---|---|
| comparator calibration fails (check 0) | re-tested at k = 10 by Experiment 1 (§0.5); matched-threshold repair specified |
| \(\sigma_{\text{mouse}} = 0.20\) untested as reference | resolved: reference is 0 (measured); 0.20 is a stress level in Exp 9 |
| Null 3 resampler, non-directional solver, matched DE-only not implemented | implemented and run in the k = 5 benchmark (`null3_resampler.R`, `nd_generator.R`, `spec_results.rds`) |
| — new — gate mismatch between simulation and real main branch | **resolved 2026-09-13**: option (a) frozen (§0.6) |

## 0.8 Cost

Experiment 1 ≈ 1.5 h single-core (4 τ levels × 2 nulls × 50 genes × GAM-fs); 1b ≈ 1.5 h (k = 20, one null); Experiments 2–3
≈ 1.5 h; 6 ≈ 2 h parallelised per k; 7 + 7b ≈ 3 h; 9 ≈ 1.5 h; 5, 8, 10, 11, 6b ≈ 4 h together.
Total ≈ 15 h single-core; an afternoon on FASRC with replicates parallelised.

---

> **Historical material follows.** §0a (v13), §0b (v14), §0c (v15) and §0-hist (the v12 core
> benchmark text) are retained for auditability. Settings stated there that differ from §0 —
> k = 5, NB size 2, \(\tau_\theta\) reference 0 or 0.15 labelled "15% heterogeneity",
> \(\sigma_{\text{mouse}} = 0.20\) as reference, three \(\tau_\theta\) levels, the per-stratum
> size rule — are **superseded and must not be used for new runs**.

---

# 0a. VERSION 13 AMENDMENTS — the slice-count correction and what follows from it

**Status.** Written 2026-09-13, before any benchmark run at the new primary design. The six
v12 experiments are complete at 5 slices per arm (`FROZEN_CONCLUSIONS.md`, 2026-08-28;
`publication_simulation_report.Rmd`). Nothing in this section deletes or re-labels those results.
It changes the primary design going forward, adds three experiments the conversation around the
step-2/step-3 explainer pages showed were missing, and corrects four places where the plan text
and the run scripts disagree. Every amendment names the file it touches.

## A0. Why v13 exists

The v12 matched pseudo-arm substrate uses **5 slices per arm** (4 pseudo-mice contributing
1, 1, 1, 2 slices). That number was inherited from the v9 idea of splitting the 10 real Day-0
slices into two literal halves. Once slices became synthetic resamples (v10 §1), the halving had
no remaining justification, but the count was carried forward and the v10 text describes
"4 mice / ~5 slices per arm" as matching the real design. It does not. The real design is:

| real arm | slices pooled for a typical gene | mice | slices per mouse |
|---|---:|---:|---|
| parent, Day 0 Fibro 2 | **10** (1,005 of 1,026 genes) | 4 | 2, 3, 3, 2 |
| daughter, Day 9 IAF 4 | 18–21 (most genes) | 7 | 1–3 |
| parent, Day 9 Fibro 2 | 9 or 19 | 7 | 1–3 |

(`arms.csv`, column `k`; slice/mouse inventory from `gene_spatial_lineage_slices.csv`.) The
real-geometry experiment (B2) already runs at 10 vs 18 real slices. The main experiments at 5
were therefore **harsher than the analysis they validate**, which is not conservative, only
inaccurate. Three v12 findings are plausibly consequences of k = 5 rather than of the method:

1. FDR experiment: proposed-method TPR 5.4% at q = 0.05 despite ~70% single-test power on the
   same genes — attributed in the report to the Knapp–Hartung t on 4–8 df.
2. Ranking reversal at \(\tau_\theta = 0.30\) (DESpace2-style 0.299 vs proposed 0.276, pilot).
3. Matched power 73% at \(|D| = 0.05,\ \tau_\theta = 0.15\).

A 200-replicate check on 2026-09-13 (seed 20260913, proposed method only, \(M = 64\),
\(\tau_\theta = 0.15\), `make_pseudo_arms()` with the patterns below) gave:

| slices per arm | pseudo-mice | Type I (Null 1) | power, gained \(|D|=0.05\) |
|---:|---:|---:|---:|
| 5 | 4 (1,1,1,2) | 4.5% [1.6, 7.4] | 69.5% |
| 10 | 4 (2,3,3,2) | **8.0% [4.2, 11.8]** | 97.0% |
| 20 | 8 (2,3,3,2 ×2) | 4.5% [1.6, 7.4] | 100% |

The k = 5 row reproduces B1 (73%) within Monte Carlo error. The k = 10 Type I point is the one
value outside the band; its interval includes 5% and it is **not** evidence of miscalibration,
but it must be resolved (A1, precondition) before k = 10 power is reported.

**Pre-registration statement.** The switch to k = 10 is recorded here, with its reason, before
the first k = 10 benchmark run. The k = 5 results are retained verbatim as a labelled
supplementary stress test. Both are reported; neither is chosen after the fact.

## A1. Primary design: 10 slices per arm

**Substrate.** `make_pseudo_arms(pool, mice_per_arm = 4, slices_per_mouse = c(2, 3, 3, 2))`.
Pooled Day-0 Fibro 2 coordinates (10,211 cells), slice sizes drawn from the observed ten,
substrate **redrawn every replicate**, randomisation unit the mouse — all unchanged from v10/v12.
Only the slice pattern changes, and it now equals the real Day-0 pattern exactly.

**Secondary design, k = 20.** `mice_per_arm = 8, slices_per_mouse = rep(c(2, 3, 3, 2), 2)`.
Matches the daughter-arm slice count. Run for the FDR experiment and the \(|D| = 0.05\) power
cells only (Experiment 8), not the full grid.

**Retained design, k = 5.** No new runs. Existing `b1_results.rds`, `b1_ds2fix_results.rds`,
`spec_results.rds`, `fdr_results.rds`, `abl*_results.rds` are reported as
"**half-design stress test (k = 5)**" alongside the k = 10 primary. The v12 conclusions that
depend on k — FDR shortfall, \(\tau_\theta = 0.30\) reversal — are restated as k = 5 findings
until the k = 10 runs confirm or retire them.

**Precondition before any k = 10 power cell is reported.** Run Null 1 at k = 10,
\(\tau_\theta \in \{0, 0.15, 0.30\}\), with the A2 background scoring (≥ 1,000 effective null
tests per cell), proposed method and both calibrated comparators. Acceptance: proposed-method
Type I inside [3%, 7%] at every \(\tau_\theta\). If the 8.0% reproduces, diagnose before
proceeding — candidate causes are the 3-slice mice dominating Knapp–Hartung's observed-scatter
rescaling, or a seed artefact; neither is mouse correlation, which the simulator did not generate.

**Which experiments re-run at k = 10.**

| experiment | k = 10 re-run? | note |
|---|---|---|
| 1 Calibration | yes | with A2 background scoring |
| 2 Directional power | yes | gained + reversed, full \(|D| \times \tau_\theta\) grid |
| 3 Specificity | yes | DE-only, non-directional, repositioning-only |
| 4A \(N_{\mathrm{eff}}\) ablation | no | asks about within-slice information; k-independent; keep k = 5 result |
| 4B REML+KH vs fixed | yes | proposed only, minutes |
| 4C blank λ | folded into Experiment 9 | see A4 |
| 5 Real geometry | no | already real 10 vs 18; unchanged |
| 6 FDR | yes, at k = 10 **and** k = 20 | the k-sensitive result |

**Reporting rule.** Every table and figure states k. Power is compared across k only at matched
realised Type I; nominal-power comparisons across k are not made.

**Cost.** GAM-fs run time scales with cells: ×2 at k = 10, ×4 at k = 20 relative to v12.
Proposed method and DESpace2 are negligible. Estimate: Experiments 1–3 at k = 10 ≈ 1.2 h
single-core; FDR ≈ 2 h parallelised per k; 4B ≈ 5 min.

## A2. Null-cell sample size: score the background panel

**The discrepancy.** §0 and the "Benchmark A — calibration" sizing rule specify
\(R = 200 \times G_0 = 10 = 2{,}000\) raw null tests per cell. The run scripts (`run_b1.R`,
`run_spec.R`, `run_b2.R`) simulate a 50-gene null background panel per replicate for DESpace2's
dispersion estimate but **score gene 1 only**, so every null cell in the frozen results rests on
200 tests. Its binomial 95% interval is ±3 percentage points — the width of the grey band in
Figure B — which separates 5% from 10% but not 5% from 7%, and falls short of the
\(n_{\text{eff}} \ge 1{,}000\) target the plan itself sets.

**Amendment.** In every **null** cell, score all 50 background genes for every method. They are
genuine nulls by construction (\(\theta = 0\) in both arms, \(M = 64\), the same \(\tau_\theta\)
and substrate as the target), so this is 10,000 raw null tests per cell for no additional
simulation. Monte Carlo SE stays cluster-robust with the replicate as the cluster (§0 formula);
\(n_{\text{eff}}\) is then reported and the ≥ 1,000 target becomes attainable.

**Power cells unchanged.** One target gene per replicate; the 50-gene backdrop stays null
(`LOCKED_METHODS.md`, bug fix 2026-08-28) and is not scored.

**Cost control.** DESpace2 already returns all 50 p-values from one call. The proposed method is
milliseconds per gene. GAM-fs is the cost: 50 fits per replicate in null cells. Either run null
cells in parallel across replicates, or score 10 of the 50 for GAM-fs and all 50 for the other
methods — state which in the results, and compute each method's \(n_{\text{eff}}\) from the
genes it actually scored.

**Script change.** `run_b1.R` line ~37 and the equivalent lines in `run_spec.R`, `run_b2.R`:
in cells where `pattern == "null"`, loop `m_proposed`, `m_gam` over `cnt[g, ]` for
\(g = 1 \dots 50\) and keep the full `m_despace2_panel(cnt, d)` vector instead of `[1]`.

## A3. Abundance: run §8 as written, and add the per-arm case

**What v12 ran.** §8 specifies three empirical strata (27 / 64 / 331 molecules per slice =
25th / 50th / 90th percentile post-gate). Only 64 was run; 27 appears inside Ablation A only;
331 was never run. Real genes of interest sit in the high stratum (Ackr3: 204–388 molecules per
Day-0 Fibro 2 slice).

**Experiment 7 — abundance strata.** \(M \in \{27, 64, 331\}\) × \(\tau_\theta \in \{0, 0.15, 0.30\}\);
Null 1 and gained at \(|D| \in \{0.025, 0.05\}\); all methods; k = 10; R = 200; A2 scoring in
null cells. The \(M = 64\) column is Experiments 1–2 at k = 10 and is not re-run.

*Pre-stated expectation, so the result is a test rather than a finding.* Per-slice
\(SE_0 \propto N_{\mathrm{eff}}^{-1/2}\): ≈ 0.016 at 64, ≈ 0.007 at 331. \(\tau^2\) does not
shrink with \(M\). At 331 the pooled SE is therefore floored by heterogeneity and power gains
saturate; methods that ignore slice structure become *more* confident as within-slice noise
shrinks, so the naive GAM's inflation should grow with \(M\), and whether GAM-fs and DESpace2
hold at 331 under \(\tau_\theta \ge 0.15\) is open. At 27 the proposed method will decline to test
some slices (\(SE_0 \ge 0.05\)); report testability rate separately, as Ablation A does.

**Experiment 7b — per-arm abundance (replaces amplified/dampened).** The real comparison
changes abundance *and* position together: Ackr3's per-cell rate in Day-9 IAF 4 is roughly
one-third of Day-0 Fibro 2, and IAF 4 slices have a median of 426 cells against 998. v12 tests
abundance change only in isolation (DE-only). Add:

- gained at \(|D| = 0.05\) with \(M_A = 64\) and \(M_B / M_A \in \{1/3, 1, 3\}\)
  (\(M_B = 21, 64, 192\)); \(\tau_\theta \in \{0, 0.15\}\);
- the matching null (\(\theta = 0\) both arms, same \(M_B/M_A\)) — this is DE-only generalised
  to both directions and is the shared-null check for this experiment;
- both substrates: matched pseudo-arms (k = 10) and real Day 0 / Day 9 geometry;
- all methods; R = 200.

`sim_design()` already accepts `M = c(M0, M9)` (engine extension log, 2026-08-28); no engine
change is needed. For the proposed method report gate attrition (slices failing
\(SE_0 < 0.05\); genes that would route `parent_only` / `daughter_only` under the real
pipeline's rules) separately from Type I and power.

**§6 consequence.** `amplified` and `dampened` are retired as scenario names; 7b covers them
with one parameterisation. The pattern-classification extra keeps `lost` only.

## A4. Mouse correlation: test the pipeline as it is actually run

**The gap.** The real pipeline applies the blank-derived \(\hat\lambda\) (`03_meta.R`, step 4a).
In the simulation \(\hat\lambda\) is applied only inside Ablation C (`run_ablC.R`); `m_proposed()`
has no λ argument, so every main result is the core method **without** λ, at
\(\sigma_{\text{mouse}} = 0\). The v12 reference setting names \(\sigma_{\text{mouse}} = 0.20\),
but no main experiment used it; the only evidence at 0.20 is Ablation C's C0 arm (Type I 6.3%
[4.8, 7.8], λ_raw 1.15, C1 → 4.5%). Precondition 2 of §0 was therefore met only in the narrow
sense.

**Amendment (i) — estimate \(\sigma_{\text{mouse}}\) from the data instead of assuming it.**
Fit a three-level REML model (`metafor::rma.mv`, slices nested in mice) to the real per-slice
\(\delta\) values in the arm where it is identifiable: the Day-9 daughter arm (IAF 4: 18 slices
from 7 mice; also Fibro 6 → IAF 2/3: 20 slices). Report, across the ~1,000 genes per lineage,
the distribution of the mouse-level variance share and the median implied
\(\sigma_{\text{mouse}}\) on the θ scale (via \(\tau_\delta \approx \tau_\theta\,SD(r)\), §9).
The Day-0 arm (4 mice) is not used for estimation.

**Amendment (ii) — Experiment 9, the deployed pipeline.** At k = 10, Null 1 and gained
\(|D| = 0.05\), \(\tau_\theta \in \{0, 0.15\}\), \(\sigma_{\text{mouse}}\) at the estimated
median **and** at 0.20 (the v12 stress value): run the proposed method **with λ on** (300 blanks
per replicate at \(M = 16\), mouse-correlated regime of §3), GAM-fs (has a mouse random effect)
and DESpace2 (does not). This is the only cell in the benchmark that tests the method exactly
as deployed, on the dependence structure the data actually have. Ablation C's C0/C1/C2 arms are
reported from this run; `run_ablC.R` at k = 5 is retained as the supplementary stress test.

**Reference-setting consequence.** If the estimated \(\sigma_{\text{mouse}}\) is near 0, the
reference becomes \(\sigma_{\text{mouse}} = \hat\sigma\) (≈ 0) and 0.20 is re-labelled a stress
level. The §0 reference box is corrected below to say so.

## A5. Heterogeneity labelling

**Promote the frozen rule to plan text.** Every figure and table reports the **realised** \(I^2\)
(median of the REML fits in that cell) beside \(\tau_\theta\). Never "τ = 0.15 (15%
heterogeneity)". Mapping on this substrate at \(M = 64\), k = 5 (`ablB_tau_cal.rds`):

| \(\tau_\theta\) | realised median \(I^2\) |
|---:|---:|
| 0 | ~11% (estimator floor at k = 5) |
| 0.15 | 44–50% |
| 0.18 | 60% |
| 0.30 | 79% |

Re-derive at k = 10 before the first k = 10 figure: the \(I^2\) estimator's floor falls as k
rises, so the τ-to-\(I^2\) map changes with the design.

**Correct the "typical" label.** Real Day-0 parent-arm median \(I^2\) is 13.5%, 90th percentile
60% (§9 table). Hence \(\tau_\theta = 0\) **is** the typical real case, \(\tau_\theta = 0.15\)
is an upper-quartile stress case, and 0.30 lies beyond the real 90th percentile. The v12
reference "\(I^2 = 15\%\)" is satisfied by \(\tau_\theta = 0\), not by \(\tau_\theta = 0.15\).
Experiments 1–2 keep all three levels; only the labels change.

**v15 correction (§0c C3).** At k = 10 the estimator floor is 0 and the real median \(I^2\) is realised by
\(\tau_\theta\approx0.05\), which becomes the reference level; the statements above about \(\tau_\theta=0\)
being the typical case hold only at k = 5.

**v16 correction.** \(\tau_\theta\) is the *imposed* heterogeneity; it equals the total only on the pooled
substrate. On real-geometry substrates total heterogeneity is summarised only by realised \(I^2\) (§0.2).

**v14 addition (§0b B3).** The \(\tau_\theta\) levels are *total* heterogeneity (on the pooled substrate). On the real
Day-0 slices about 0.06 of it is geometric — present even for a biologically homogeneous effect
gene — and the pooled substrate erases it.

## A6. Eligibility gate: align the real pipeline with the simulated one

**The mismatch.** The simulated proposed method gates each slice on \(SE_0 < 0.05\) and
≥ 30 cells (`LOCKED_METHODS.md` §1). The real pipeline (`delta_core.R` ~line 104) gates on
≥ 30 cells **and** (≥ 20 molecules **or** molecules ≥ `MIN_RATE_PER_1K` × cells / 1000). The
rate alternative admits very sparse slices: for Ackr3 in the Day-9 daughter arm,
`100221_D9_m3_2_slice_2` (50 cells, **1** molecule) and in the parent-D9 arm
`062921_D9_m5_2_slice_2` (77 cells, **2** molecules) are both `eligible = TRUE`. A \(\delta\)
from one molecule is one cell's position. Its \(SE_0\) is finite and large, so it enters REML
with low weight, but it still counts toward \(k\), toward the Knapp–Hartung df, and toward the
\(I^2\) denominator.

**Amendment.** The real pipeline adopts the simulation's gate: ≥ 30 cells **and**
\(SE_0 < 0.05\), replacing the molecule/rate floor. \(SE_0\) is self-calibrating (a sparse
slice automatically fails), it is the criterion every simulated result assumes, and it removes a
tunable constant. Re-run `01_run_chunk.R` → `03_meta.R`; report the number of gene × slice rows
whose eligibility changes, the distribution of \(k\) before and after, and any D1 calls that
change at q = 0.05. Until this is done the benchmark validates a gate the real analysis does not
use.

## A7. Experiment table, v13 — superseded by §0b B9

| # | experiment | k | scenarios | grid | methods | status |
|---|---|---:|---|---|---|---|
| 1 | Calibration | 10 | Null 1, Null 2 | \(\tau_\theta \in \{0, .15, .30\}\); A2 scoring | all 5 | re-run |
| 2 | Directional power | 10 | gained, reversed | \(|D| \in \{.025,.05,.075\} \times \tau_\theta\) | all 5 | re-run |
| 3 | Specificity | 10 | DE-only, repositioning-only, non-directional | reference | all 5 + pseudobulk | re-run |
| 4A | \(N_{\mathrm{eff}}\) ablation | 5 | null, gained .05 (DM) | \(N_{\mathrm{eff}} \in \{5,20,50,80\}\) | proposed | keep |
| 4B | REML+KH vs fixed | 10 | null, gained | realised \(I^2 \approx \{11\%, 60\%\}\) | proposed | re-run |
| 5 | Real geometry | real (10 v 18) | Null 3A/3B, gained, reversed | \(|D| = .05\) | all 5 | keep |
| 6 | FDR | 10, 20 | 400 null + 100 gained | \(\tau_\theta = .15\) | proposed, GAM-fs, DESpace2, pseudobulk | re-run ×2 |
| 7 | Abundance strata | 10 | Null 1, gained \(\{.025,.05\}\) | \(M \in \{27, 331\} \times \tau_\theta\) | all 5 | new |
| 7b | Per-arm abundance | 10 + real | gained .05, matching null | \(M_B/M_A \in \{1/3, 1, 3\} \times \tau_\theta \in \{0,.15\}\) | all 5 | new |
| 8 | Slice count | 20 | Null 1, gained .05 | \(\tau_\theta\) | all 5 | new |
| 9 | Deployed pipeline | 10 | Null 1, gained .05 | \(\sigma_{\text{mouse}} \in \{\hat\sigma, .20\} \times \tau_\theta \in \{0,.15\}\), λ on | proposed (full), GAM-fs, DESpace2 | new; replaces 4C |
| — | k = 5 results | 5 | all v12 cells | as run | as run | retained, labelled |

**Cost, single-core unless stated.** Experiments 1–3 at k = 10 ≈ 1.2 h. 4B ≈ 5 min. 6 ≈ 2 h
parallelised per k (×2). 7 ≈ 1.5 h (two new M columns). 7b ≈ 1.5 h (two substrates). 8 ≈ 1 h
(GAM at ×4). 9 ≈ 1 h (blanks are cheap; GAM-fs dominates). Total ≈ 10–12 h, or an afternoon on
FASRC with null cells parallelised across replicates.

## A8. Order of operations — superseded by §0b B10

1. **A6** — adopt the \(SE_0\) gate in the real pipeline; re-run `03_meta.R`; record what changed.
2. **A2** — patch the three run scripts to score the background panel in null cells.
3. **A1 precondition** — Null 1 at k = 10 with A2 scoring; resolve the 8.0% point.
4. **A5** — re-derive the τ-to-\(I^2\) map at k = 10.
5. Re-run Experiments 1, 2, 3, 4B, 6 at k = 10. Re-freeze `FROZEN_CONCLUSIONS.md` as v13.
6. **A4(i)** — estimate \(\sigma_{\text{mouse}}\) from the Day-9 arms; fix the reference value.
7. Experiments 7, 7b, 9, then 8.
8. Regenerate `publication_simulation_report.Rmd` with k on every figure.

## A9. What v13 does not change

\(\delta_s\), \(SE_0\), the \(SE_0 < 0.05\) gate, REML + Knapp–Hartung, the D1 contrast and its
Welch–Satterthwaite df, the floored median \(\hat\lambda\) — all as in `LOCKED_METHODS.md` §1.
Comparator specifications (`LOCKED_METHODS.md` §2–5). θ solved per arm per substrate by
`uniroot`. \(M = 64\) as the reference stratum, \(M = 16\) for blanks. NB size 2 as the
reference; **one addition**: at \(M = 331\) the over-dispersion term \(\mu^2/\text{size}\) is no
longer negligible, so the Experiment 7 null cell at \(M = 331\) is also run at size 0.7 (the
real fitted median), proposed method only, as a sensitivity row. Substrate redrawn every
replicate; randomisation unit the mouse; cluster-robust Monte Carlo SE. Deferred, recorded, not
run: Fibro 6 substrate replication; \(\pi_0 \in \{0.5, 0.9\}\).

---

# 0b. VERSION 14 AMENDMENTS — selection, generator fidelity, and robustness

**Status.** Written 2026-09-13, same day as v13, before any v13 run. Prompted by an external
review of v13 and by two measurements made the same day (geometry heterogeneity of the real
slices; pooled versus within-slice-bootstrap substrates). v14 does not change the v13 primary
design (k = 10, §0a A1) or any method. It adds reporting rules that make selection visible,
a generator-fidelity check that must pass before the abundance experiments run, four robustness
cells, and one wording rule for the paper. The v13 experiment table (§0a A7) is superseded by
the table in B9 below.

## B1. Gate-induced selection is part of the estimand

**The issue.** Any per-slice gate decides which slices enter the meta-analysis. If the
probability of passing depends on the effect, the pooled estimate is computed on an
outcome-dependent subset. This is true of the transcript-floor gate in the main real pipeline
and of the \(SE_0 < 0.05\) gate in the simulation and in the stratified sensitivity branch
(`01c_stratified_chunk.R`, `SE_MAX`). A directional lean concentrates transcripts into fewer
cells, lowering \(N_{\mathrm{eff}}\) and raising \(SE_0\), so an effect gene is marginally
*less* testable than a null gene at the same \(M\). At \(M = 64\) with ~1,000 cells the gate
essentially never fires (\(SE_0 \approx 0.016\)); it fires at \(M = 27\), in Ablation A, and in
the per-arm abundance scenario (\(M_B = 21\)) — exactly the cells v13 added.

**Reporting rule, every cell, every method that gates.** Alongside Type I error and power:

\[
P(\text{testable}),\qquad k_0,\qquad k_9,\qquad
\text{power}\mid\text{testable},\qquad
P(\text{correct discovery}) = P(\text{testable})\times\text{power}\mid\text{testable}.
\]

"Testable" means the arm-level pool received at least `MIN_K` slices in both arms and the
D1 contrast was emitted. Ablation A already reports the first quantity; this extends the
convention to all cells. Under the null, report \(P(\text{testable})\) and Type I conditional
on testable; a null gene that is not tested is neither an error nor a success.

**One gate, not two.** v13 A6 proposes the main real pipeline adopt \(SE_0 < 0.05\). v14 makes
the requirement symmetric and leaves the choice open: **the benchmark and the headline real
analysis must use the same per-slice gate.** Either (a) the real main branch adopts
\(SE_0 < 0.05\) (A6 as written), or (b) the simulation adopts the real transcript/rate floor
(`MIN_MOL`, `MIN_RATE_PER_1K`) in `engine.R` `per_slice_table()` and the results are re-run.
Option (a) is preferred: the rule is self-calibrating, it is already the primary rule of the
stratified branch, and it has no tunable constant. Whichever is chosen, record it in
`LOCKED_METHODS.md` §1 and report the gate-attrition table for the real data (rows changing
eligibility; \(k\) distribution before and after; D1 calls changing at q = 0.05).

## B2. Generator fidelity: a posterior-predictive check on molecule concentration — done, see §0c C2 (size 1)

**The issue.** The primary generator is NB with size 2. The method's SE depends on
\(N_{\mathrm{eff}} = (\sum n_i)^2 / \sum n_i^2\), so if the NB does not reproduce how real
transcripts concentrate across cells, the simulated SE is systematically wrong in a direction
that could favour or penalise the method. `FROZEN_CONCLUSIONS.md` records one number: simulated
\(N_{\mathrm{eff}}\) is ~4% above real genes at \(M = 64\) (real fitted size ≈ 0.7). That is one
abundance and one summary.

**Required check, before Experiment 7.** For real Day-0 Fibro 2 gene × slice observations,
stratified by transcript total (bins around 27, 64, 331) and for simulated slices at the same
\(M\) under size ∈ {2, 0.7}, compare the distributions of

\[
N_{\text{event}}=\sum_i n_i,\qquad
N_{\mathrm{eff}},\qquad
N_{\mathrm{eff}}/N_{\text{event}},\qquad
\max_i n_i \big/ \sum_i n_i .
\]

Report medians and 10th/90th percentiles per stratum, real versus simulated, as one table
(`ppc_neff.csv`) and one figure. **Acceptance:** simulated median \(N_{\mathrm{eff}}/N_{\text{event}}\)
within 10% of real in every stratum. If size 2 fails at any stratum, that stratum is run at the
real fitted size, and the sensitivity row already planned for \(M = 331\) (§0a A9) becomes the
primary at that abundance. The blank distribution (\(M = 16\)) is checked the same way against
the 48 real blanks.

## B3. Substrate realism: state it, and test it once

**Geometry heterogeneity measured (2026-09-13).** The ten real Day-0 Fibro 2 slices have
mean \(r\) 0.43–0.50 and SD 0.09–0.15. Painting the gained gene (\(\theta = 0.356\)) onto each
real slice gives expected \(\delta\) from 0.035 to 0.057, SD 0.0078, equivalent to
\(\tau_\theta \approx 0.06\). The pooled pseudo-slice substrate erases this: every pseudo-slice has
SD 0.128. Under the null it is irrelevant (\(\theta = 0\) has nothing to scale). Under an
effect, real geometry adds slice-to-slice spread the simulation currently attributes entirely to
\(\tau_\theta\).

**Measured consequence, 200 replicates, \(\tau_\theta = 0\), k = 5, proposed method:**

| substrate | Type I (Null 1) | median realised \(I^2\) | power, gained 0.05 |
|---|---:|---:|---:|
| pooled pseudo-slices | 4.0% | 6.5% | 97.0% |
| within-slice bootstrap of the real slices, arms by mouse (2 v 2) | 6.5% | 4.2% | 97.5% |

Calibration and power agree within Monte Carlo error.

**Amendments.**
- §0a A5 gains the sentence: *the \(\tau_\theta\) levels are total heterogeneity; on real
  slices about 0.06 of it is geometric and present even for a biologically homogeneous gene.*
- The plan states plainly that the matched pseudo-arm substrate is the **internal-validity**
  experiment (exact null, independent arms, any k) and that pseudo-mice are grouping labels,
  not animals: they do not differ in geometry, abundance, variance or composition unless
  \(\sigma_{\text{mouse}}\) or B5 says so. Experiment 5 and the within-slice bootstrap are the
  **realism** experiments. Results from the two are never pooled in one figure.
- The within-slice bootstrap (k = 5, 2 v 2 by mouse) is added to Experiment 5 as one cell:
  Null 1 and gained 0.05 at \(\tau_\theta \in \{0, 0.15\}\), all methods. The numbers above are
  its proposed-method pilot.

## B4. Directional shapes outside the exponential family

**The issue.** All directional alternatives use \(h(r) = e^{\theta z}\), a log-linear tilt. That
is the comparators' own model family, so it is not generous to \(\delta\); but a first-moment
statistic could still be favoured if every alternative moves the whole distribution smoothly.
Two shapes with the same \(D\) but different form test whether power depends on the family.

**Experiment 11 — awkward directional shapes.** Gained at \(|D| = 0.05\), \(\tau_\theta = 0.15\),
k = 10, all methods, R = 200:

- **Localised shift.** A fraction \(\pi\) of the gene's expected mass is moved to a narrow
  window at the outer end (\(r \in [0.75, 0.85]\)); the remainder keeps \(h = 1\).
  \(\pi\) solved so the expected \(\delta\) equals 0.05 on the substrate.
- **Bimodal redistribution.** \(h(r) = 1 + \tfrac{\kappa}{3}\,\tilde\phi(r; 0.25, 0.05) + \kappa\,\tilde\phi(r; 0.75, 0.05)\),
  \(\tilde\phi\) the normal density scaled to peak 1; **one free parameter** \(\kappa\) (bump ratio fixed at 1:3),
  solved by `uniroot` so expected \(\delta = 0.05\). Total mass is preserved by the
  \(\mu_i = M h_i / \sum h\) normalisation for any \(h\) and is not a second constraint.
  *(v16 wording; matches `engine_v16.R::h_bimodal`. The v14 text "solved jointly so mass is preserved" was
  imprecise.)*

Both keep \(M\) fixed. Engine change: `sim_slice()` accepts an arbitrary weight vector `h` in
place of `theta`; the solver generalises to any one-parameter family via `uniroot` on
`delta_exp`. Expectation stated in advance: \(\delta\) has the same power as for the exponential
tilt (it sees only the first moment); GAM-fs gains on the bimodal shape; DESpace2 gains on the
localised shift if the window aligns with a domain boundary and loses if it straddles one.

## B5. Non-normal mouse effects

**The issue.** Mouse and slice heterogeneity enter as additive normal effects. With four mice per
arm, REML + Knapp–Hartung is sensitive to one unusual animal, and real between-mouse
differences are not guaranteed to be normal.

**Experiment 9, two added cells** (proposed with λ on, GAM-fs, DESpace2; k = 10; Null 1 and
gained 0.05; \(\tau_\theta = 0.15\); R = 200):

- **Heavy-tailed:** \(b_m \sim \sigma_{\text{mouse}} \cdot t_3 / \sqrt{3}\) (variance matched to
  the normal case).
- **One outlying mouse:** \(b_m = 0\) for three mice and \(b_m = 3\sigma_{\text{mouse}}\) for one,
  the outlier's arm randomised per replicate.

Report Type I and power, and for the proposed method the realised \(\hat\lambda\) and the
Knapp–Hartung inflation factor \(q\), so the reader can see whether the correction responds.
Blanks receive the same mouse effect distribution (§3 rule), otherwise \(\hat\lambda\) cannot
see it.

## B6. FDR under correlated gene effects

**The issue.** The 500 genes in Experiment 6 share a substrate but have independent
heterogeneity nudges, so there are no gene modules. Real tissue has co-regulated blocks whose
slice-level deviations move together; BH's guarantee holds under positive regression
dependence, but the realised FDR of a miscalibrated method can be much worse under
correlation than under independence.

**Experiment 6b — one sensitivity condition.** As Experiment 6 (400 null + 100 gained,
\(\tau_\theta = 0.15\), k = 10, 50 replicates) but with the genes partitioned into 25 blocks of
20; within a block the slice nudges are \(u_{gs} = \rho\, w_{bs} + \sqrt{1-\rho^2}\, e_{gs}\)
with \(w_{bs}\) shared by the block, \(\rho = 0.7\), marginal variance \(\tau_\theta^2\)
unchanged. Report realised FDR and TPR at q ∈ {0.05, 0.10} beside the independent-gene
result. Proposed, GAM-fs, DESpace2.

## B7. Required diagnostics for λ

The median \(\hat\lambda\) corrects scale, not tail shape (§3, `FROZEN_CONCLUSIONS.md`
claim 3). A 5% rejection rate after correction can conceal poor behaviour at \(10^{-3}\).
**Required figure in Experiment 9 and in the real-data results:** QQ plot of blank
\(z\)-statistics against \(N(0,1)\), before and after division by \(\sqrt{\hat\lambda}\), with
the 0.05, 0.01 and 0.001 tail quantiles tabulated. The C2 tail statistic
\(c_{0.05} = Q_{0.95}(|Z_{\text{blank}}|)\) is reported on the same figure. A method-level claim of
calibration requires the post-λ QQ plot to lie inside the pointwise 95% band down to
\(p = 10^{-3}\); below that, state that calibration is unverified.

## B8. Wording rule for cross-method claims

Every comparative statement in the paper is scoped. The template is:

> *For directional first-moment repositioning, on a matched substrate, at a matched
> false-positive operating point, method X had power P₁ against P₂.*

The unscoped forms "more powerful", "more sensitive to spatial change", and "best spatial
method" are not used. The non-directional scenario (Experiment 3) is cited wherever the
comparators' sensitivity to changes \(\delta\) ignores is relevant, and it is described as the
comparators answering a different question, not as a failure of either side.

## B9. Experiment table, v14 — superseded by §0c C7

| # | experiment | k | scenarios | grid | methods | status |
|---|---|---:|---|---|---|---|
| 0 | Generator fidelity (B2) | — | real vs simulated slices | M ∈ {27, 64, 331}; size ∈ {2, 0.7}; blanks | — | **new; gate for 7** |
| 1 | Calibration | 10 | Null 1, Null 2 | \(\tau_\theta \in \{0,.15,.30\}\); 50-gene scoring | all 5 | re-run |
| 2 | Directional power | 10 | gained, reversed | \(|D| \in \{.025,.05,.075\} \times \tau_\theta\) | all 5 | re-run |
| 3 | Specificity | 10 | DE-only, repositioning-only, non-directional | reference | all 5 + pseudobulk | re-run |
| 4A | \(N_{\mathrm{eff}}\) ablation | 5 | null, gained .05 (DM) | \(N_{\mathrm{eff}} \in \{5,20,50,80\}\) | proposed | keep |
| 4B | REML+KH vs fixed | 10 | null, gained | realised \(I^2 \approx \{11\%,60\%\}\) | proposed | re-run |
| 5 | Real geometry | real 10 v 18; **+ bootstrap 5 v 5** | Null 3A/3B, gained, reversed; **+ Null 1, gained (B3)** | \(|D| = .05\) | all 5 | keep + 1 cell |
| 6 | FDR | 10, 20 | 400 null + 100 gained | \(\tau_\theta = .15\) | proposed, GAM-fs, DESpace2, pseudobulk | re-run ×2 |
| 6b | FDR, correlated genes (B6) | 10 | as 6, 25 blocks, \(\rho = .7\) | \(\tau_\theta = .15\) | proposed, GAM-fs, DESpace2 | new |
| 7 | Abundance strata | 10 | Null 1, gained {.025,.05} | \(M \in \{27,331\} \times \tau_\theta\); size per B2 | all 5 | new |
| 7b | Per-arm abundance | 10 + real | gained .05, matching null | \(M_B/M_A \in \{1/3,1,3\} \times \tau_\theta \in \{0,.15\}\) | all 5 | new |
| 8 | Slice count | 20 | Null 1, gained .05 | \(\tau_\theta\) | all 5 | new |
| 9 | Deployed pipeline | 10 | Null 1, gained .05 | \(\sigma_{\text{mouse}} \in \{\hat\sigma,.20\} \times \tau_\theta \in \{0,.15\}\); λ on; **+ \(t_3\) and one-outlier mouse (B5)**; QQ (B7) | proposed full, GAM-fs, DESpace2 | new; replaces 4C |
| 10 | Second lineage | 10 | Null 1, gained .05 | Fibro 6 substrate, \(M\) = Fibro 6 median, \(\tau_\theta \in \{0,.15\}\) | proposed, GAM-fs | promoted from deferred |
| 11 | Awkward directional shapes (B4) | 10 | localised shift, bimodal | \(|D| = .05\), \(\tau_\theta = .15\) | all 5 | new |
| — | k = 5 results | 5 | all v12 cells | as run | as run | retained, labelled |

Every cell reports \(P(\text{testable})\), \(k_0\), \(k_9\), power | testable and
\(P(\text{correct discovery})\) (B1), and realised \(I^2\) beside \(\tau_\theta\) (§0a A5).

**Added cost.** B2 minutes. 5-bootstrap cell ≈ 15 min. 6b ≈ 2 h parallelised. 9 additions
≈ 40 min. 10 ≈ 45 min (two methods). 11 ≈ 1 h. v14 total over v13 ≈ 5 h; whole plan ≈ 15–17 h
single-core, an afternoon on FASRC.

## B10. Order of operations — superseded by §0c C8

1. **B1 gate decision** — choose (a) or (b); apply to the real pipeline or the engine; record in
   `LOCKED_METHODS.md`; re-run `03_meta.R` and report the attrition table.
2. **§0a A2** — patch run scripts to score the 50-gene background in null cells; add the B1
   testability columns to every output.
3. **§0a A1 precondition** — Null 1 at k = 10, ≥ 1,000 effective null tests; resolve the 8.0%.
4. **B2** — generator fidelity table; fix size per stratum if needed.
5. **§0a A5** — τ-to-\(I^2\) map at k = 10.
6. Experiments 1, 2, 3, 4B, 6 at k = 10; re-freeze `FROZEN_CONCLUSIONS.md` as v14.
7. **§0a A4(i)** — estimate \(\sigma_{\text{mouse}}\) from the Day-9 arms.
8. Experiments 9 (with B5, B7), 7, 7b, 5-bootstrap, 11.
9. Experiments 8, 6b, 10.
10. Regenerate `publication_simulation_report.Rmd`; every figure states k, substrate,
    gate, and realised \(I^2\); every comparative sentence follows B8.

## B11. What v14 does not change

Everything in §0a A9. The k = 10 primary. The exponential tilt as the *reference* alternative
(B4 adds shapes; it does not replace the family used for θ-calibration). The NB generator as
primary, subject to B2's per-stratum size choice. Normal random effects as the reference
(B5 adds robustness cells; it does not change the generative default). BH within one
lineage × layer × route family. Deferred and recorded: \(\pi_0 \in \{0.5, 0.9\}\); a substrate from
a tissue other than colon.

---

# 0c. VERSION 15 — calibration-stage outcomes and the settings they fix

**Status.** Written 2026-09-13 after the five local calibration checks in
`simulation/v16/` (`CALIBRATION_STAGE.md`; scripts `00_fidelity.R`, `01_null_k10.R`,
`02_tau_i2_map_k10.R`, `03_sigma_mouse.R`, `04_smoke.R`). These were the first four steps of
§0b B10 plus the engine smoke tests. No cluster run has started. v15 records what the checks
found, fixes the four generative settings they determine, and marks the corresponding v13/v14
items resolved. Every number below is from those outputs.

## C1. The k = 10 precondition is met — §0a A1 resolved

Null 1, proposed method, 50 null genes × 200 replicates per cell, cluster-robust MC SE:

| k | τ_θ | Type I | 95% CI | n_eff |
|---:|---:|---:|---|---:|
| 10 | 0 | 4.71% | [4.29, 5.13] | 9,546 |
| 10 | 0.15 | 4.46% | [4.08, 4.84] | 11,174 |
| 10 | 0.30 | 4.42% | [3.99, 4.85] | 8,736 |
| **10** | **0.05, size 1** | **4.44%** | **[4.01, 4.87]** | **9,014** |
| 5 | 0 / 0.15 | 4.47% / 4.34% | | ~10,000 |
| 20 | 0 / 0.15 | 4.96% / 4.99% | | ~10,000 |

The 8.0% observed with 200 tests on 2026-09-13 was Monte Carlo noise. The τ_θ = 0.05 row was added in v16
after C3 made it the reference (realised I² 15.2%). **k = 10 is confirmed as primary, calibrated at all
four prespecified τ_θ levels.** The A2 scoring rule (50 background genes per null replicate) delivers the
n_eff ≥ 1,000 target with a factor of ~10 to spare; it stays.

## C2. NB size = 1 — supersedes §0b B2's per-stratum rule and §0a A9's size-0.7 row

Real Day-0 Fibro 2, 9,900 gene × slice rows, median \(N_{\mathrm{eff}}/N_{\text{event}}\)
versus simulated pseudo-slices at the same abundance:

| stratum | real | size 0.7 | **size 1** | size 2 |
|---|---:|---:|---:|---:|
| 27 | 0.944 | 0.941 | 0.949 | 1.000 |
| 64 | 0.889 | 0.875 | 0.893 | 0.921 |
| 331 | 0.611 | 0.560 | 0.604 | 0.673 |
| blanks, 16 | 1.000 | 1.000 | 1.000 | 1.000 |

Size 2 overstates the effective count by 4% at the median abundance and 10% at the high
stratum (the B2 acceptance boundary). Size 0.7 undershoots at 331. **Size 1 is within 1% of
real at every stratum** and matches the real maximum-cell-share distribution (median 0.036 at
M = 64: real 0.036, size 1 0.034, size 2 0.031).

**Setting.** `size = 1` in every new run, all experiments, genes and blanks. B2's "run the
stratum at the real fitted size" clause is moot; B2 itself is retained as the check that was
run, and is re-run only if the substrate changes (Experiment 10, Fibro 6). The frozen k = 5
results stay at size 2 and carry the label "size 2; simulated SE ≈ 2% optimistic relative to
real genes at M = 64".

*Interpretation of the tail (v16).* Real genes have a wider lower tail in \(N_{\mathrm{eff}}/N_{\text{event}}\)
at M = 64 (p10 0.758 vs 0.815 for size 1). One NB size reproduces the central information concentration,
not the full gene-to-gene burstiness distribution. The pre-specified criterion was the median, size 1
passes it, and it removes size 2's systematic optimism; no further generator layer is added. The caveat
is recorded and Ablation 4A remains the stress test for extreme concentration.

## C3. τ_θ levels — corrects §0a A5 and §0b B3

Realised median \(I^2\) on null genes, 100 replicates:

| τ_θ | k = 10 | k = 5 |
|---:|---:|---:|
| 0 | 0 | 0 |
| **0.05** | **14.1** | 0 |
| 0.10 | 35.1 | 27.6 |
| 0.15 | 53.9 | 55.2 |
| 0.20 | 68.3 | 63.2 |
| 0.30 | 83.8 | 79.2 |

At k = 10 the estimator's floor is 0, not the ~11% recorded at k = 5, and the real Day-0 median
\(I^2\) of 13.5% is realised by **τ_θ ≈ 0.05**, not by τ_θ = 0. The geometric heterogeneity of
the ten real slices measured in B3 (τ_θ ≈ 0.06 equivalent for an effect gene) is the same size,
which supports reading the real median heterogeneity as largely geometric.

**Setting.** Experiments 1 and 2 run **four** τ_θ levels, {0, 0.05, 0.15, 0.30}, labelled
"none", "real median (I² ≈ 14%)", "upper quartile (I² ≈ 54%)", "beyond real 90th percentile
(I² ≈ 84%)". Every other experiment that fixes τ_θ = 0.15 keeps it and carries the
"upper quartile" label. The §0 reference setting is τ_θ = 0.05.

## C4. σ_mouse — resolves §0a A4(i); fixes Experiment 9's levels

Three-level REML (slices in mice) on the real Day-9 daughter arms, 8 mice, 18–21 slices,
~1,000 genes across four lineages:

| | genes | blanks |
|---|---:|---:|
| median mouse share of between-slice variance | 0–11% | — |
| median σ_mouse, θ scale | 0.000 | 0.000 |
| 75th percentile | 0.04–0.07 | — |
| 90th percentile | 0.08–0.12 | 0.14 |

Blanks carry no mouse biology, so their upper tail is the estimation floor with 8 mice; the
genes' tail is not distinguishable from it. \(\hat\sigma_{\text{mouse}} \approx 0\) for the
typical gene. The v2 decision not to model mouse and the observed real \(\hat\lambda \approx 1\)
are consistent with this. v12's 0.20 is about twice the real 90th percentile.

**Setting.** Experiment 9 runs \(\sigma_{\text{mouse}} \in \{0.07, 0.20\}\), labelled
"real upper quartile" and "stress"; the near-zero estimate itself is not run because it
duplicates Experiment 1. The §0 reference setting is \(\sigma_{\text{mouse}} = 0\).

## C5. Engine extensions verified — §0b B4, B5, B6 implementable as written

`engine_v16.R` adds arbitrary-shape generation, block-correlated nudges, t₃ and one-outlier
mouse effects, and per-arm M, without redefining any `engine.R` function. Smoke tests pass:
the localised-shift and bimodal solvers hit δ = 0.05 exactly; totals stay at M for every
shape; block nudges give marginal sd 0.147 and within-block correlation 0.46 (target 0.49);
`sim_design_u` reproduces `sim_design` bit-for-bit given identical draws. On the pool,
δ = 0.05 requires π = 0.155 (localised, window [0.75, 0.85]) and κ = 19.0 (bimodal).

## C6. Reference setting, v15 (supersedes the §0 box and the v13 correction under it)

\[
\boxed{k=10,\qquad M=64,\qquad \text{NB size}=1,\qquad \tau_\theta=0.05\ (I^2\approx14\%),\qquad
\sigma_{\text{mouse}}=0,\qquad \text{matched geometry},\qquad R=200}
\]

Crossed stress levels: τ_θ ∈ {0, 0.15, 0.30} in Experiments 1–2; σ_mouse ∈ {0.07, 0.20} in
Experiment 9; k = 20 in Experiments 6 and 8; M ∈ {27, 331} in Experiment 7.

## C7. Experiment table, v15 (supersedes §0b B9)

| # | experiment | k | scenarios | grid | methods | status |
|---|---|---:|---|---|---|---|
| 0 | Generator fidelity | — | real vs simulated | M ∈ {27,64,331}, blanks; size ∈ {0.7,1,2} | — | **done → size 1** |
| 1 | Calibration | 10 | Null 1, Null 2 | τ_θ ∈ {0,.05,.15,.30}; 50-gene scoring | all 5 | cluster |
| 2 | Directional power | 10 | gained, reversed | \(|D|\) ∈ {.025,.05,.075} × τ_θ (4 levels) | all 5 | cluster |
| 3 | Specificity | 10 | DE-only, repositioning-only, non-directional | τ_θ = .15 | all 5 + pseudobulk | cluster |
| 4A | \(N_{\mathrm{eff}}\) ablation | 5 | as v12 | as v12 | proposed | keep (size 2, labelled) |
| 4B | REML+KH vs fixed | 10 | null, gained | τ_θ ∈ {0, .18} | proposed | local, minutes |
| 5 | Real geometry + bootstrap | real 10 v 18; 5 v 5 | Null 3A/3B, gained, reversed; Null 1, gained | \(|D|\) = .05 | all 5 | **done** |
| 6 | FDR | 10, 20 | 400 null + 100 gained | τ_θ = .15 | proposed, GAM-fs, DESpace2, pseudobulk | cluster ×2 |
| 6b | FDR, correlated blocks | 10 | as 6, ρ = .7, blocks of 20 | τ_θ = .15 | proposed, GAM-fs, DESpace2 | **done** |
| 7 | Abundance strata | 10 | Null 1, gained {.025,.05} | M ∈ {27,331} × τ_θ ∈ {0,.15,.30} | all 5 | **done** |
| 7b | Per-arm abundance | 10 + real | gained .05, matching null | \(M_B/M_A\) ∈ {1/3,1,3} × τ_θ ∈ {0,.15} | all 5 | **done** |
| 8 | Slice count | 20 | Null 1, gained .05 | τ_θ ∈ {0,.15,.30} | all 5 | **done** |
| 9 | Deployed pipeline | 10 | Null 1, gained .05 | σ_mouse ∈ {.07,.20} × τ_θ ∈ {0,.15} × {normal, t₃, outlier}; λ on; QQ | proposed full, GAM-fs, DESpace2 | **done** |
| 10 | Second lineage | 10 | Null 1, gained .05 | Fibro 6 pool, M = Fibro 6 median, τ_θ ∈ {0,.15}; fidelity re-checked | proposed, GAM-fs | **done** |
| 11 | Awkward shapes | 10 | localised (π = .155), bimodal (κ = 19.0) | \(|D|\) = .05, τ_θ = .15 | all 5 | **done** |
| — | Null 1, proposed, k ∈ {5,10,20} | | | | | **done (C1)** |
| — | τ_θ → I² map | 10, 5 | | | | **done (C3)** |
| — | σ_mouse estimate | real D9 | | | | **done (C4)** |
| — | k = 5 results | 5 | all v12 cells | as run, size 2 | as run | retained, labelled |

All runs at NB size 1 on the Fibro 2 substrate; **size is substrate-specific** and is set by Experiment 0 on
each substrate (Fibro 6: size 2). Exceptions: 4A and the retained k = 5 results (size 2, Fibro 2, labelled). Every cell reports
P(testable), k₀, k₉, power | testable, P(correct discovery) (B1) and realised I² (A5).

## C8. Order of operations, v15 (supersedes §0b B10)

1. **§0b B1 gate decision** — the one open item: (a) real main branch adopts SE₀ < 0.05, or
   (b) the engine adopts the transcript/rate floor. Record in `LOCKED_METHODS.md`; produce the
   real-data attrition table. *Awaiting the analyst's choice.*
2. ~~A2 scoring patch~~ implemented in `v16/01_null_k10.R`; port to the cluster scripts.
3. ~~A1 precondition~~ done (C1).
4. ~~B2 fidelity~~ done (C2): size 1.
5. ~~A5 τ map~~ done (C3): add τ_θ = 0.05.
6. ~~A4(i) σ_mouse~~ done (C4): 0.07 / 0.20.
7. Write cluster scripts against C6; gate as a parameter until step 1 is decided.
8. Experiments 1, 2, 3, 4B, 6 → re-freeze `FROZEN_CONCLUSIONS.md` as v15.
9. Experiments 9, 7, 7b, 5, 11, then 8, 6b, 10.
10. Regenerate `publication_simulation_report.Rmd`: every figure states k, size, substrate,
    gate, realised I²; every comparative sentence follows §0b B8.

## C9. What v15 does not change

Everything in §0a A9 and §0b B11 except the two settings they fixed differently: NB size is
now 1 (was 2 with a 0.7 sensitivity), and the reference τ_θ is 0.05 (was 0). The method,
comparator specifications, k = 10, the exponential tilt as the θ-calibration family, normal
random effects as the generative default, and all completed k = 5 results are unchanged.

---

# 0-hist. v12 CORE BENCHMARK TEXT — HISTORICAL, SUPERSEDED BY §0

> Retained verbatim (with the v13/v15 in-place notes) for pre-registration. Do not run from this section.

**Everything after §1 is the master validation specification. It is a checklist, not a run
list.** Implementing every branch of it as the primary benchmark would take 12–17 core-hours and
produce a figure no reader can follow. This section defines the **six experiments that are
actually run**, on the principle:

\[
\boxed{\text{one clean experiment per scientific claim}}
\qquad\text{not}\qquad
\boxed{\text{every scenario crossed with every parameter}}
\]

The manuscript should say plainly that the full specification exists and that the primary
benchmark is a principled subset of it. That answers *"why didn't you test X?"* without running X.

## The four claims the paper must establish

1. the method is **calibrated** under a true null;
2. it has **power** for directional repositioning;
3. it measures something **different** from ordinary DE and from general spatial-profile change;
4. its distinctive components — \(N_{\mathrm{eff}}\), random effects, calibration — **matter**.

\[
\boxed{\text{Calibration}}\rightarrow\boxed{\text{Directional power}}
\rightarrow\boxed{\text{Specificity}}\rightarrow\boxed{\text{Ablation}}
\]

## Reference setting

Vary **one factor at a time** around a fixed reference, rather than estimating
\(\text{power}=f(D,\ \text{abundance},\ I^2,\ \sigma_{\text{mouse}},\ \text{pattern},\dots)\)
over a full Cartesian product:

\[
\boxed{k=10,\qquad M=64,\qquad \text{NB size}=1,\qquad \tau_\theta=0.05\ (I^2\approx14\%),\qquad
\sigma_{\text{mouse}}=0,\qquad \text{matched geometry},\qquad R=200}
\]

**v15.** Fixed by the calibration stage (§0c C6). The v13 box below it is retained as the record.

**v13.** The v12 box read \(M=64,\ I^2=15\%,\ \sigma_{\text{mouse}}=0.20\) with an implicit
k = 5. Three corrections (§0a): k = 10 is the real Day-0 design; \(I^2=15\%\) is realised by
\(\tau_\theta=0\), not 0.15 (A5); and \(\sigma_{\text{mouse}}=0.20\) was never run in a main
experiment — it is now a stress level, with the reference estimated from data (A4).

**Two exceptions, both forced by pilot evidence — see "Factors that must be crossed" below.**

## The six experiments

| # | experiment | scenarios | methods | grid |
|---|---|---|---|---|
| **1** | **Calibration** | Null 1, Null 2 | proposed, hierarchical GAM, DESpace2 (+ naive GAM once, as negative control) | \(\tau_\theta\in\{0,0.15,0.30\}\) |
| **2** | **Directional power** | gained, reversed | proposed, hierarchical GAM, DESpace2 | \(\|D\|\in\{0,\ 0.025,\ 0.05,\ 0.075\}\times\tau_\theta\in\{0,0.15,0.30\}\) |
| **3** | **Specificity / boundary** | DE-only (matched), non-directional | + pseudobulk DE | reference only |
| **4** | **Ablations** | A: \(N_{\mathrm{eff}}=5\) vs \(80\) at \(N_{\mathrm{event}}=100\); B: \(I^2=0\) vs \(60\%\); C: no \(\lambda\) vs \(\lambda\) | proposed only | C crossed with \(\sigma_{\text{mouse}}\in\{0,0.20\}\) |
| **5** | **Realism check** | Null 3, gained, reversed | all | real Day 0 / Day 9 geometry, \(\|D\|=0.05\) |
| **6** | **FDR** | one mixture: \(G=500\), \(\pi_0=0.8\) | all | reference only |

**v13.** All rows now at k = 10 slices per arm (§0a A1); null cells score the 50-gene
background panel (A2); Experiments 7, 7b, 8, 9 added and 4C folded into 9 — see the v13 table
in §0a A7. The k = 5 results for rows 1–6 are retained as a labelled stress test.

Plus one cheap extra: **pattern-classification accuracy** — all five patterns at \(\|D\|=0.05\),
one \(\tau\), proposed only. It is the one metric no comparator produces, and dropping
`lost`/`amplified`/`dampened` from the primary power figure would otherwise remove the material
for it.

**Demoted to supplementary:** `lost`, `amplified`, `dampened` power curves; \(\pi_0\in\{0.5,0.9\}\);
the C2 tail-calibration diagnostic; the full \(\sigma_{\text{mouse}}\) and \(I^2\) grids; B2 as a
headline ranking claim.

## Factors that must be crossed, not varied one at a time

**One-factor-at-a-time is right for abundance, \(\|D\|\) and \(\pi_0\). It is wrong for
\(\tau_\theta\), and wrong for \(\sigma_{\text{mouse}}\) inside Ablation C.**

The pilot measured a **ranking reversal** across \(\tau_\theta\). Power at matched empirical
threshold, mean over five patterns:

| \(\tau_\theta\) | proposed | DESpace2-style |
|---|---:|---:|
| 0 | **0.957** | 0.756 |
| 0.15 | **0.760** | 0.698 |
| 0.30 | 0.276 | **0.299** |

DESpace2 **overtakes the proposed method at \(\tau_\theta=0.30\)**. Reporting Experiment 2 at a
single \(\tau=15\%\) would state an ordering that reverses just outside the reference setting,
and any reviewer who re-ran it at higher heterogeneity would find the reversal. \(\tau_\theta\)
is therefore crossed with method in Experiments 1 and 2.

Likewise \(\sigma_{\text{mouse}}\) cannot be held fixed inside **Ablation C**: inflation is the
thing that ablation exists to correct, so at \(\sigma_{\text{mouse}}=0\) it has nothing to show
(the pilot confirmed this — C0 4.5%, C1 4.25%, both already nominal).

## Two preconditions before any of this runs

**These are blockers, and the restructure does not remove them.**

1. **Comparator calibration (§1b check 0) currently FAILS.** Pilot Type I error under a true
   null: hierarchical GAM 3.0% → **9.0%**, DESpace2-style **10.5% → 11.0%**, naive GAM
   6.0% → **43.5%** (at \(\tau_\theta=0\) → \(0.30\)). The proposed method passes (6.0% → 3.5%).
   Comparing power against three miscalibrated comparators is the strawman §1b exists to
   prevent, and simplifying the design does not make them comparable. Fix or formally scope
   them, and re-run the gate, before Experiments 1–3.
2. **The reference setting contains an untested value.** \(\sigma_{\text{mouse}}=0.20\) has never
   been run — the pilot used \(\sigma_{\text{mouse}}=0\) throughout. Verify the proposed method
   still calibrates there before anchoring six experiments to it.

Additionally, three §6 scenarios have **no implementation yet** and are required by Experiments 3
and 5: Null 3's moment-balanced \(\gamma_1,\gamma_2\) resampler, the non-directional \(\kappa\)/\(c\)
solver, and the matched-geometry DE-only construction. Smoke-test each before scaling.

## Final figure

| panel | content |
|---|---|
| **A** | simulation concept — null, directional, non-directional |
| **B** | Type I calibration (Exp 1), across \(\tau_\theta\) |
| **C** | power versus \(\|D\|\) (Exp 2) |
| **D** | specificity: DE-only vs directional vs non-directional (Exp 3) |
| **E** | ablations: \(N_{\mathrm{eff}}\), random effects, calibration (Exp 4) |
| **F** | real-geometry robustness (Exp 5) + FDR (Exp 6) |

The expected Experiment-3 result is the paper's clearest single statement:

| scenario | proposed | GAM / DESpace2 | pseudobulk DE |
|---|---|---|---|
| DE-only | − | − | **+** |
| directional repositioning | **+** | possibly + | − |
| non-directional spatial change | − | **+** | − |

## Cost

| | single-core |
|---|---:|
| master specification (§1 onward, every branch) | ~12–17 h |
| Experiments 1–5 | **~35 min** |
| Experiment 6 (FDR) | ~7.6 h, or **~1 h** parallelised |

## What must not be simplified away

\[
\delta_s=\frac{\sum_i n_i r_i}{\sum_i n_i}-\bar r_{\text{cells}},
\qquad
SE_0=\sqrt{s_r^2\!\left(\tfrac{1}{N_{\mathrm{eff}}}-\tfrac1N\right)},
\qquad
\text{REML}+\text{Knapp–Hartung}.
\]

These are what make the method statistically distinctive. What is simplified away is **repeated
demonstration of the same point** — four directional patterns all showing directional detection,
four \(I^2\) levels all showing heterogeneity sensitivity, three calibration arms where the
headline needs two.

---


# Objective

## Objective

Benchmark the proposed repositioning method against the most relevant alternatives for **cross-condition spatial change**, while separately testing whether each component of the proposed framework contributes to performance.

The benchmark should answer four questions:

1. **Detection:** Can the method detect true directional repositioning while controlling Type I error and FDR?
2. **Comparison:** How does it perform relative to a continuous profile model and a discrete-domain model?
3. **Specificity:** Can it distinguish spatial repositioning from ordinary differential expression?
4. **Ablation:** Which parts of the framework—effective positional SE, random-effects pooling, and blank calibration—actually improve performance?

## Benchmark structure — three benchmarks, not one grid

Running every scenario x effect size x abundance x method x ablation at the sample size needed
for *calibration* claims is both unnecessary and unaffordable. Calibration questions need many
null tests and few methods; power questions need many methods and few replicates. The benchmark
is therefore split into three, each sized for the question it answers.

### Benchmark A — calibration

Scenarios that are **null for the population-relative repositioning target** \(D_\delta=0\)
(Null 1, Null 2, Null 3, DE-only). Proposed method plus the three ablations; GAM and DESpace2
are run for reference but are not the target. This is where "is Type I error 5%, 7% or 9%?" is
settled, so it is sized for **effective null tests**, not for replicates.

**"Null" is method-specific, and Null 3 is the case where it diverges.** Null 3 holds
\(D_\delta=0\) while deliberately changing the marginal spatial profile, so the null hypothesis
is true for the repositioning estimand and *false* for the profile-based estimands:

| scenario | proposed (\(\delta\)) | absolute-position variant | GAM | DESpace2 | pseudobulk DE |
|---|---|---|---|---|---|
| Null 1 unpositioned | null | null | null | null | null |
| Null 2 stable positioned | null | null | null | null | null |
| **DE-only** | null | null | null | null | **NOT null** |
| **Null 3 population movement** | **null** | **null** (its stated claim) | **NOT null** | **NOT null** | null |

**This table is a prediction, not a fact — it must be verified per scenario and per stratum.**
See the comparator-truth audit below, which shows that Null 2 and DE-only are *not* reliably
shared nulls once θ is solved per arm.

Both bolded rows are cases where a method is *supposed* to reject. Pseudobulk DE rejecting under
DE-only is the intended result of §2, not a Type I error — it is the whole point of including
it. Only \(M\) changes there, so the profile methods stay null (their condition main effect
absorbs the abundance shift) while the abundance test fires.

Consequently, **a GAM rejection rate of 70% under Null 3 is not a 70% Type I error.** The GAM is
answering a different question correctly: \(f_0(r)\neq f_9(r)\) genuinely holds. Report
GAM/DESpace2 values under Null 3 as a **rejection rate**, never as Type I error, or omit them
from that bar group entirely.

This distinction strengthens rather than weakens the comparison: it shows the methods target
different estimands, instead of implying the comparators are miscalibrated. The one genuine
false positive under Null 3 remains an absolute-position statistic that *claims* to measure
population-relative displacement — see §6.

\[
R=200\ \text{replicates}\ \times\ G_0=50\ \text{null genes}=10{,}000\ \text{raw null tests per cell.}
\]

**v13.** \(G_0\) was 10 in v12, and the run scripts actually scored \(G_0=1\) (the target gene
only), leaving 200 tests per cell. The 50-gene null background panel is already simulated in
every replicate for DESpace2; v13 scores all of it in null cells (§0a A2).

**Report the cluster-robust Monte Carlo SE, not the naive one.** Genes within a replicate share
the same coordinates, the same mouse effects and the same \(\hat\lambda\), so they are not
independent and \(\sqrt{\hat\alpha(1-\hat\alpha)/(RG_0)}\) understates the error. Let \(p_r\) be
the rejection proportion within replicate \(r\). Then

\[
\hat\alpha=\frac{1}{R}\sum_r p_r,
\qquad
\widehat{SE}=\frac{\mathrm{sd}(p_r)}{\sqrt{R}},
\qquad
n_{\text{eff}}=\frac{\hat\alpha(1-\hat\alpha)}{\widehat{SE}^{2}}.
\]

Target \(n_{\text{eff}}\ge1{,}000\), which resolves 5% against 7% at ~80% power. **If
\(n_{\text{eff}}\) falls short, raise \(R\), not \(G_0\)** — replicates are independent by
construction, extra genes inside a replicate are not.

### Shared nulls require matched cell populations — and Simulation A cannot provide them

**A shared null cannot be built on genuinely different Day 0 and Day 9 populations.** This is
forced, not incidental. Under Null 1 the gene is flat (\(h\equiv1\)), so
\(\mu_i=M/N\) for every cell and the expected expression share of domain \(d\) equals the
*cell* share of domain \(d\). If the two arms' coordinate distributions differ, the domain
profile differs, and \(H_{0,\mathrm{DESpace2}}\) is false for a completely flat gene. Measured
on the real Fibro 2 (D0) versus IAF 4 (D9) coordinates with \(K=5\):

| domain | Day 0 share | Day 9 share | difference |
|---|---:|---:|---:|
| D1 | 0.0074 | 0.0326 | +0.0251 |
| D2 | 0.3222 | 0.3358 | +0.0136 |
| D3 | 0.5212 | 0.5525 | +0.0313 |
| D4 | 0.1385 | 0.0745 | **−0.0639** |
| D5 | 0.0107 | 0.0046 | −0.0061 |

\(\chi^2\) test of equal domain profile: \(p=2.4\times10^{-97}\).

The proposed method is unaffected — \(\delta\) subtracts \(\bar r_{\text{cells}}\), so the
population difference cancels, which is exactly the design property Null 3 tests. That is the
whole difficulty: **\(\delta\) and the profile methods cannot share a null on differing
populations, by construction.**

#### Consequence

Three things in this plan require a genuine shared null and therefore **must not be run on
mismatched populations**:

1. the §1b comparator-calibration gate (pilot check 0);
2. Benchmark C's \(G_0=400\) common nulls;
3. the matched-geometry DE-only variant.

#### The construction to use

Build shared-null scenarios by **splitting one condition's cells into two pseudo-arms**:

> Take Day-0 slices only. Simulate the gene with the same \(h(r)\) in both pseudo-arms. Every
> method is then genuinely null: same population, same geometry, same profile, no condition
> effect of any kind.

**The randomization unit must be the mouse, not the slice.** "Assign slices while respecting
mouse nesting" is self-contradictory: if two slices from mouse A land in opposite pseudo-arms,
the arms are not independent and the nesting is destroyed. Assign **whole mice**, carrying all
their slices:

\[
\text{mouse}\ \longrightarrow\ \text{pseudo-condition},
\qquad
\text{all slices of that mouse follow.}
\]

There are only 4 Day-0 mice, so a literal split gives ~2 mice per arm — enough for a
calibration null, likely too small for the B1 power comparison. Commit to **both**, labelled:

1. **Synthetic pseudo-mice (primary for B1).** Build pseudo-slices by resampling cells from the
   pooled Day-0 coordinate distribution, with slice sizes drawn from the observed distribution;
   group pseudo-slices into pseudo-mice; assign pseudo-mice to arms. This preserves real
   geometry while letting the mouse and slice counts match the real design (4 mice / ~5 slices
   per arm) instead of halving them.
2. **Literal mouse split (realism check).** The real 4 Day-0 mice split 2 v 2, carrying their
   slices. Small, but uses the untouched empirical design.

Report which construction produced each number. Do not leave the assignment rule implicit.

This keeps the real colon geometry, unlike a fully synthetic substitute, while removing the
population mismatch. It is the standard negative-control construction and is what the §1b gate
and the Benchmark C nulls should use. Simulation B (§4), where both arms draw from the same
coordinate distribution by definition, is an acceptable alternative substrate.

**Scenarios that deliberately keep the real Day 0 / Day 9 populations** — Null 3, the
realistic-geometry DE-only, and all power scenarios — remain as specified, but their comparator
rejections are scored per method through the truth audit below, never as shared Type I error.

Record which substrate each reported number came from. A Type I error from a
mismatched-population scenario and one from a matched pseudo-arm scenario are not the same
quantity and must not appear in the same bar group.

### The comparator-truth audit — run before any simulation

Declaring hypotheses shared in advance is not rigorous; compute them. The v4 fix that solves
\(\theta_0\) and \(\theta_9\) **separately** (so that \(\delta_0=\delta_9\) when
\(SD(r_0)\neq SD(r_9)\)) has a consequence that reaches the comparators: if \(\theta_0\ne\theta_9\)
then \(h_0\) and \(h_9\) are genuinely different functions, so a profile-based method can
legitimately detect a difference even though \(D_\delta=0\).

The GAM null is \(s_\Delta(r)=0\), i.e. \(h_9(r)/h_0(r)\) constant in \(r\). With
\(h_c(r)=\exp(\theta_c z)\) the ratio is \(\exp\!\left[(\theta_9-\theta_0)z\right]\), which is
constant **only** if \(\theta_9=\theta_0\). Measured in the real geometries at a
\(\delta=0.05\) target:

| stratum | \(SD(r_0)\) | \(SD(r_9)\) | \(\theta_0\) | \(\theta_9\) | \(h_9/h_0\) swing over \(z\in[-2.5,2.5]\) | GAM null |
|---|---:|---:|---:|---:|---:|---|
| Fibro 6 → IAF 3 | 0.1178 | 0.1406 | 0.424 | 0.356 | **41%** | **FALSE** |
| Fibro 2 → IAF 4 | 0.1276 | 0.1280 | 0.392 | 0.391 | 1% | effectively true |

So in Fibro 6 → IAF 3, **Null 2 and DE-only are not shared nulls**: a 41% swing in the
condition ratio across the axis is a real profile-shape difference that a GAM fitted to
thousands of cells will detect, correctly. Counting that as a Type I error would repeat the
Null-3 mistake in a subtler place.

**Required audit.** For every scenario × stratum, before drawing counts, compute the expected
condition-specific profiles and domain means from \(\mu_{sci}\) and classify each null
separately:

\[
H_{0,\delta}:\ D_\delta=0,
\qquad
H_{0,\mathrm{GAM}}:\ \frac{h_9(r)}{h_0(r)}\ \text{constant},
\qquad
H_{0,\mathrm{DESpace2}}:\ \frac{\sum_{i\in d}\mu_{9i}}{\sum_i\mu_{9i}}
=\frac{\sum_{i\in d}\mu_{0i}}{\sum_i\mu_{0i}}\ \ \forall d .
\]

Record the achieved ratio swing and the maximum domain-proportion difference alongside the
achieved \(\delta\) and \(D\) (§5). **Only score a comparator rejection as a Type I error when
that comparator's own generative null is actually true.**

### DE-only comes in two variants

Because of the above, DE-only is split:

1. **DE-only (realistic geometry)** — different cell populations per arm, \(\theta\) solved
   separately, \(\delta_0=\delta_9=0.05\). Null for \(\delta\); *possibly not* null for
   GAM/DESpace2, per the audit. Use this for the \(\delta\)-versus-pseudobulk specificity claim
   of §2 and Panel E.
2. **DE-only (matched geometry)** — the **shared** null, built on the matched pseudo-arm
   substrate above. Use the *same* coordinate distribution and the *same* \(h(r)\) in both
   pseudo-conditions, so \(h_0\equiv h_9\) exactly, and change only \(M_9=2M_0\). Every spatial method is then genuinely null, and this is the
   variant that belongs in the shared Type I error comparison.

The matched-geometry variant sacrifices the realism of distinct Day 0 / Day 9 populations; that
is the deliberate trade for a clean common null. Report both and say which supports which claim.

### Benchmark B — comparative power

Directional patterns and the non-directional control across the \(|D|\) grid, all three spatial
methods. \(R=100\) for the pilot, \(R=200\) for the final curves. Power differences between
methods are large relative to their Monte Carlo error; calibration differences are not, which is
why this benchmark can be smaller than A.

### Benchmark C — mixed-gene FDR

The **only** benchmark containing null and alternative genes simultaneously, and therefore the
only one in which FDR is defined at all.

**Primary Benchmark C is directional only.** This is a deliberate restriction. If the mixture
also contained non-directional profile changes or population-movement genes, there would be
**no single \(\pi_0\)**: a non-directional gene has \(D_\delta=0\) but a changed profile, so it
is null for the proposed method and alternative for GAM and DESpace2, and the three methods
would be scored against three different truths. Restricting the mixture to a directional
contrast gives all methods **one common scientific truth**:

\[
G=500,
\qquad
G_0=400\ \text{common nulls}\ (\text{matched pseudo-arms; }D=0,\ \text{no profile change}),
\qquad
G_1=100\ \text{directional alternatives}\ (D\neq0)\ \text{on the same substrate},
\]
\[
\pi_0=0.8,
\qquad
R=100,
\]

\[
\boxed{\text{All 500 genes in Benchmark C are generated on the same matched pseudo-arm substrate}}
\]

— the same slices, the same pseudo-conditions, the same geometry. The **only** difference
between a null and an alternative gene is its simulated localization effect. Otherwise one BH
family would span genes simulated under different experimental designs.

With \(\pi_0\in\{0.5,\ 0.9\}\) as sensitivity analyses. The 100 alternatives are a mixture of
gained, lost, reversed, amplified and dampened, with effect sizes drawn from
\(|D|\in\{0.01,\ 0.025,\ 0.05,\ 0.075\}\), so FDR is evaluated against a realistic blend of
strong and marginal signals rather than one magnitude.

The non-directional profile scenario stays in **Benchmark B** as a boundary/specificity
experiment, where its differing estimands are the point rather than a confound. A secondary
method-specific FDR experiment could be added later; it is not needed for the main result.

### Multiple-testing family

BH is applied within a testing family, so the family must be stated or the reported \(q\)
values are not reproducible. For the primary simulation:

\[
\boxed{\text{all }500\text{ genes form ONE lineage}\times\text{layer}\times\text{route family}}
\qquad
q_g=\mathrm{BH}(p_1,\dots,p_{500}).
\]

Every method's own correction is applied over the same 500 genes, so the FDR comparison is
like-for-like. If a variant splits the genes over several families, say so explicitly and
report the family sizes.

**These are design choices, not quantities implied by the colon data.** State them as such.

---

# 1. Main methods to compare

## Proposed repositioning method

For each slice,

\[
\delta_s
=
\frac{\sum_i n_i r_i}{\sum_i n_i}
-
\frac{1}{N}\sum_i r_i.
\]

This measures the displacement of the transcript-weighted mean position from the mean position of the corresponding cell population.

Condition-level effects are estimated by random-effects pooling, and the main contrast is

\[
D=\hat\mu_9-\hat\mu_0.
\]

The method is specifically designed to detect **directional positional shifts along a biologically defined continuous axis**.

---

## Comparator 1: continuous condition-interaction GAM

Fit a count model with a condition-specific smooth spatial effect:

\[
Y_i\sim NB(\mu_i,\phi),
\]

\[
\log\mu_i
=
\beta_0
+
\beta_1C_i
+
s_0(r_i)
+
C_i\,s_{\Delta}(r_i).
\]

Test

\[
H_0:s_{\Delta}(r)=0.
\]

Interpretation:

- proposed method: tests whether the **first spatial moment moves**;
- GAM: tests whether the **entire continuous spatial expression profile changes**.

This comparator is especially important because it directly asks whether modeling the full profile provides an advantage over the lower-dimensional \(\delta\) statistic.

---

## Comparator 2: DESpace2

Discretize the continuous axis into spatial domains and fit a condition × domain interaction model.

Conceptually,

\[
\log \mu
=
\text{domain}
+
\text{condition}
+
\text{condition}\times\text{domain}.
\]

**Implement it as the framework actually works, not as a cell-level regression.** DESpace2
pseudobulks within each biological sample × spatial domain and fits the negative-binomial
domain × condition model *across samples*. So aggregate first,

\[
X_{gsd}=\sum_{i\in(s,d)}Y_{gsi},
\]

where \(s\) indexes slice/sample and \(d\) indexes domain, then apply the package's multi-sample
DSP test with the slice as the unit of replication. A cell-level
`domain + condition + interaction` fit is **not** DESpace2 and would reproduce the
pseudoreplication failure quantified in §1b.

Test

\[
H_0:
\beta_{\text{condition}\times\text{domain}}=0.
\]

Use \(K=5\) domains as the primary analysis, with \(K=3\) and \(K=10\) as sensitivity analyses.

Interpretation:

- proposed method: continuous directional displacement;
- GAM: continuous full-profile change;
- DESpace2: discrete domain-level change.

---

---

# 1b. Comparator specification — the fairness precondition

**A comparator specified without the hierarchical structure the proposed method has is a
strawman, and its Type I error will be catastrophic for reasons that have nothing to do with
its estimand.** This section must be settled before any result is generated.

## The problem, measured

The comparator in §1 is written as a cell-level model,

\[
\log\mu_i=\beta_0+\beta_1C_i+s_0(r_i)+C_i\,s_\Delta(r_i),
\]

with no slice or mouse term. Fitted to pooled cells it performs inference on ~3,000 cells while
there are only ~10 independent biological units. Simulated under a **true null** — no condition
effect whatsoever, only slice-to-slice variation in the gradient (\(\tau_\theta=0.30\),
5 slices per arm, 300 cells per slice):

| specification | Type I error |
|---|---:|
| GAM, no random effects (§1 as written) | **54%** |
| GAM + random intercept `s(slice, bs="re")` | **61%** |
| GAM + random intercept and random slope | **47%** |
| proposed \(\delta\), slice-level contrast | **5%** |
| nominal | 5% |

Two things follow. First, publishing Panel B from the §1 specification would report a ~54% Type
I error for the GAM that is an artifact of **our** specification choice. Second — and less
obvious — **a random intercept does not fix it.** The heterogeneity is in the *slope*: slices
differ in their spatial gradient, not their mean expression. Even a random intercept plus random
slope leaves ~47%, because with 5 slices per arm the condition-level smooth is barely
distinguishable from between-slice variation, and the smooth's *p*-value is computed as though
cells were the replicate unit.

This is the same pseudoreplication the proposed method's REML slice-level pooling exists to
avoid. It is not evidence that GAMs are bad; it is evidence that **the replicate unit must be
the slice for every method in the comparison.**

## Required specification

Pre-specify, before seeing results, and report the exact model formula for each comparator:

1. **Preferred — one hierarchical GAM across all slices.** Condition is **constant within a
   slice**, so there is no within-slice variation in \(C\) and the interaction
   \(C_i s_\Delta(r_i)\) is **not identifiable inside a single slice**. A "per-slice condition
   interaction" does not exist and cannot be meta-analyzed. Fit instead

   \[
   \log\mu_{si}=\beta_0+\beta_1C_s+s_0(r_{si})+C_s\,s_\Delta(r_{si})
   +\underbrace{f_s(r_{si})}_{\text{slice-specific random smooth}}+\underbrace{b_{m(s)}}_{\text{mouse effect}},
   \]

   with the slice-specific term as a factor-smooth (`s(r, slice, bs="fs")`) and a mouse random
   effect, then **keep the specification only if it passes the matched-null gate below.**
2. **Or pseudobulk to the replicate unit.** Aggregate counts to slice x domain and fit with
   slice as the unit of replication. This is the standard multi-sample spatial-DE design and is
   the natural mode for a domain-based method; use DESpace2's replicate-aware/multi-sample
   interface if available rather than pooling cells.
3. **Or a genuine mixed model** with a factor-smooth, e.g. `s(r, slice, bs="fs")`, and inference
   that respects the number of slices — verified, not assumed, by the gate below.

## Implementation traps that silently invalidate the benchmark

**The pilot hit one of these, and it produced a publishable-looking result that was pure
mis-specification.** Recording it because it will recur in any reimplementation.

### mgcv: the `by=` factor must be ORDERED

With an **unordered** factor `C`, `s(r, by = C)` fits **one smooth per level** — a separate
\(s_0(r)\) and \(s_9(r)\), not their difference. Testing those terms answers *"is there a
spatial gradient?"*, not *"does the gradient differ between conditions?"* Any stably positioned
gene then rejects almost always, because a gradient genuinely exists in both arms.

In the pilot this gave the GAM comparators **93.5% and 92.0% Type I error under Null 2**, a
scenario the comparator-truth audit certifies as a genuine shared null:

```
unordered C:   s(r):C0   F = 41.1   p < 1e-4        <- both arms have a gradient
               s(r):C1   F =  4.6   p = 0.024
ordered  Co:   s(r):Co1  F =  1.7   p = 0.147       <- the difference smooth, correctly null
```

**Required:** make the condition factor ordered, so `s(r, by = Co)` *is* the difference smooth
\(s_\Delta(r)\), and extract the p-value from that term specifically — never `min()` over
whatever rows happen to match a name pattern.

### Null 2 is a comparator-validity check, not only a scenario

The general lesson is larger than mgcv. **Any method that rejects under Null 2 is testing the
wrong hypothesis** — it is detecting the presence of spatial structure rather than a
*between-condition change* in it. Null 2's stated purpose in §6 ("verify that spatial
localization itself is not mistaken for cross-condition repositioning") therefore applies to the
comparators' *implementations*, not just to the estimands.

**Run Null 2 on every comparator before trusting any of its numbers, and treat a high rejection
rate there as an implementation bug until proven otherwise.** A method that fails Null 2 will
also show spuriously high "power" in every directional scenario, because those scenarios also
contain a gradient — so the failure inflates the headline result rather than depressing it,
which is the direction that does not get noticed.

## The gate: comparator calibration is a precondition, not a result

**Run every method under Null 1 first and check its Type I error before running any power
scenario** — on the **matched pseudo-arm substrate** defined in §1's Benchmark A block, not on
mismatched Day 0 / Day 9 populations, where a flat gene already violates
\(H_{0,\mathrm{DESpace2}}\) at \(p\approx10^{-97}\). Power measured at an uncontrolled false-positive rate is meaningless: a method
rejecting 54% of nulls will also show high "power", and the comparison is vacuous.

If a comparator cannot be brought near nominal by any specification above, do **not** simply
report its inflated power. Instead:

> **Compare power at matched empirical Type I error.** For each method, calibrate its rejection
> threshold on the Null 1 simulations so that its empirical false-positive rate is exactly 5%,
> then apply that threshold in the power scenarios.

### The no-random-effects GAM is a negative control, not a competitor

Keep the cell-level GAM in the benchmark — **do not delete it** — but assign it an explicit role:

| | role | appears in |
|---|---|---|
| **GAM, no random effects** | **negative control**: demonstrates why replicate-awareness is required | Panel B (Type I) only |
| **GAM, factor-smooth + mouse RE** | the real profile-based competitor | Panels B, C, D |

It earns its place by *failing*: the pilot measured its Type I error at 6.0% when slices are
homogeneous (\(\tau_\theta=0\)) and **43.5%** at \(\tau_\theta=0.30\), while the factor-smooth
version stays at 9.0%. That contrast **is** the §1b argument, made with numbers rather than
assertion, and it belongs in the paper.

**But exclude it from the power panels.** Two reasons:

1. Power at an uncontrolled false-positive rate is not power. At \(\tau_\theta>0\) this
   specification has no valid operating point to compare against.
2. The matched-threshold repair is a **simulation-only device**. It requires calibrating against
   a *known* null — which is precisely what real data does not provide, and precisely why the
   blank-derived \(\lambda\) machinery of §3 exists. A method that becomes usable only when
   handed the true null distribution is not usable in practice, however well it scores once
   handed it.

Report its matched-threshold power in a supplement if desired, labelled as **not achievable on
real data**. The head-to-head power claim in Panel C is between the proposed method, the
**factor-smooth** GAM, and DESpace2.

**The threshold and the power scenario must share a substrate.** Calibrating \(c_m\) so that
\(P_{\text{matched null}}(T_m>c_m)=0.05\) and then applying \(c_m\) to power under real
Day 0 / Day 9 geometry is only valid if

\[
\mathcal{L}\!\left(T_m\mid H_0,\text{matched}\right)\approx
\mathcal{L}\!\left(T_m\mid H_0,\text{real geometry}\right),
\]

and the comparator-truth audit exists precisely because that does **not** hold for GAM and
DESpace2 — geometry changes their null behaviour materially. A threshold learned on pseudo-arms
is therefore not a 5% threshold under mismatched populations. Benchmark B is split accordingly:

\[
\boxed{\textbf{B1: matched-geometry method comparison (primary)}}
\qquad
\boxed{\textbf{B2: real-geometry robustness (secondary)}}
\]

**B1** — threshold *and* power both on the matched pseudo-arm substrate, with gained / lost /
reversed / amplified / dampened introduced purely through \(h(r)\). Every method then sees the
same population geometry, the same replicate structure, the same false-positive operating point,
and the only thing that changes is the gene's localization. This is the clean answer to "which
method has more power for a directional repositioning signal?"

**B2** — real Day 0 / Day 9 geometry, reported as a realism check. **No shared Type I error is
claimed here**, because the methods genuinely target different quantities under mismatched
populations; results are interpreted through the truth audit, not against a common 5% level.

This is the standard remedy and makes Panel C interpretable regardless of each method's native
calibration. Report both the nominal-threshold and matched-threshold power curves, and state
which is primary. Panel B keeps the native, uncalibrated Type I errors — that is its purpose —
but Panel C must not compare power across methods operating at different real error rates.

---

# 2. Orthogonal comparator: pseudobulk differential expression

Pseudobulk DE is not a competing spatial method. Its purpose is to demonstrate that

\[
\text{abundance change}
\neq
\text{spatial repositioning}.
\]

It should be included in specific simulation scenarios to show:

### DE-only

\[
M_9\neq M_0,
\qquad
\delta_9=\delta_0.
\]

Expected:

- pseudobulk DE: significant;
- repositioning method: not significant.

### Repositioning-only

\[
M_9=M_0,
\qquad
\delta_9\neq\delta_0.
\]

Expected:

- pseudobulk DE: not significant;
- repositioning method: significant.

This directly demonstrates the scale-invariant nature of \(\delta\).

---

# 3. Internal ablation analyses

Three ablations should be compared against the full repositioning framework.

## Ablation A: naive molecule-level SE

Replace the effective-position uncertainty

\[
SE_{0,s}
=
\sqrt{
s_r^2
\left(
\frac{1}{N_{\mathrm{eff}}}
-
\frac{1}{N}
\right)
}
\]

with the explicit naive estimate

\[
SE_{\mathrm{naive},s}
=
\frac{s_r}{\sqrt{N_{\mathrm{event}}}},
\qquad
N_{\mathrm{event}}=\sum_i n_i.
\]

This treats every transcript molecule as an independent positional observation and therefore
ignores the fact that many molecules may come from the same cell. It is the direct
pseudoreplication error that the \(N_{\mathrm{eff}}\)-based uncertainty is designed to avoid.

Purpose:

> Test whether \(N_{\mathrm{eff}}\)-aware uncertainty improves calibration when the same number of transcripts is distributed across very different numbers of cells.

### The contrast must be generated, not hoped for

The stress test this ablation needs is

- 100 molecules spread across many cells;
- 100 molecules concentrated in only a few cells,

i.e. two genes with the **same** \(N_{\mathrm{event}}=\sum_i n_i\) but very different

\[
N_{\mathrm{eff}}=\frac{\left(\sum_i n_i\right)^{2}}{\sum_i n_i^{2}}.
\]

**The generative model of §5 does not produce this contrast on its own.** Varying gene
abundance \(M\) changes \(N_{\mathrm{event}}\) and \(N_{\mathrm{eff}}\) together, so it cannot
isolate concentration; and \(\phi_g\) appears in \(Y\sim NB(\mu,\phi_g)\) without any defined
sweep. Left as written, the ablation would depend on whichever \(N_{\mathrm{eff}}\) spread the
NB draws happened to deliver.

Add an explicit **concentration (burstiness) parameter** that is independent of abundance. With
baseline cell probabilities \(p_i\propto L_{si}h_c(r_{si})\), draw

\[
q\sim\operatorname{Dirichlet}(\alpha\,p),
\qquad
(n_1,\dots,n_N)\sim\operatorname{Multinomial}\!\left(N_{\mathrm{event}},\,q\right).
\]

Large \(\alpha\) spreads molecules diffusely across cells; small \(\alpha\) concentrates them
in a few. Holding \(N_{\mathrm{event}}=100\) fixed, tune \(\alpha\) to hit targets such as

\[
N_{\mathrm{event}}=100,\ N_{\mathrm{eff}}\approx80
\qquad\text{versus}\qquad
N_{\mathrm{event}}=100,\ N_{\mathrm{eff}}\approx5 .
\]

The naive estimator reads both as \(SE_{\mathrm{naive}}=s_r/10\), while the
\(N_{\mathrm{eff}}\)-based estimator separates them by roughly a factor of four. That is the
comparison the ablation exists to make.

Controlling \(N_{\mathrm{eff}}\) directly is preferred over sweeping \(\phi_g\), because
\(N_{\mathrm{eff}}\) *is* the quantity the ablation is meant to justify. Record the achieved
\(N_{\mathrm{event}}\) and \(N_{\mathrm{eff}}\) per gene and slice, since they are the true
independent variable.

**Dirichlet-multinomial is used only for this ablation** unless stated otherwise; the main
benchmark retains the NB model of §5. Note that DM preserves \(\sum_i n_i\) exactly, which is
what makes the abundance-matched contrast clean.

---

## Ablation B: no random-effects pooling

Compare the full random-effects procedure

\[
\delta_s
\rightarrow
\hat\mu_c
\]

with **one** specified alternative — fixed-effect inverse-variance pooling, i.e. the same
estimator with the between-slice variance component forced to zero:

\[
\hat\mu_c^{\text{fixed}}
=
\frac{\sum_s \delta_s/s_s^{2}}{\sum_s 1/s_s^{2}},
\qquad
\tau^{2}=0,
\]

with the corresponding fixed-effect standard error and a normal reference distribution.

\[
\boxed{\text{Full: REML }\hat\tau^{2}+\text{Knapp–Hartung}}
\qquad\text{versus}\qquad
\boxed{\text{Ablation B: fixed-effect pooling, }\tau^{2}=0}
\]

**This is the clean contrast and the ablation should not be left as "such as".** The two
estimators share weighting scheme, replicate unit and data; they differ *only* in whether
between-slice heterogeneity is estimated and propagated. Pooling all cells across slices — the
other candidate — would change the replicate unit, the weighting and the uncertainty model
simultaneously, so it cannot attribute any difference to the random-effects component. Keep
cell-pooling as an optional secondary contrast if desired, clearly labelled as confounded.

Purpose:

> Test whether REML + Knapp–Hartung improves inference when slice effects are heterogeneous,
> holding everything except the variance component fixed.

---

## Ablation C: no blank calibration

Compare **three** calibration arms:

\[
\text{C0: no calibration } (z)
\qquad
\text{C1: floored median-}\lambda\ \left(z_{\mathrm{cal}}=z/\sqrt{\hat\lambda_{\text{applied}}}\right)
\qquad
\text{C2: tail calibration (diagnostic)}
\]

where \(\hat\lambda_{\text{applied}}=\max(1,\hat\lambda_{\text{raw}})\) — the floor is part of
the method under test, not a detail (see below).

C1 is the current pipeline and remains the **primary** proposed method. C2 is a *diagnostic
sensitivity analysis*, not a redefinition of the method — see below.

Purpose:

> Test whether empirical-null calibration improves Type I error and FDR when null statistics are
> inflated, and diagnose *what kind* of inflation is present.

### Why a third arm: the estimator is centred, the error is in the tail

The calibration estimates

\[
\lambda=\frac{\operatorname{median}(Z^{2})}{\operatorname{median}(\chi_1^{2})},
\]

which asks *is the typical null statistic over-dispersed?* Type I error asks
\(P(|Z|>1.96)\) — *is the tail too heavy?* These coincide **only when the null departs from
\(N(0,1)\) by a scale factor**. If the inflation really is pure scale,

\[
Z_{\text{null}}\sim N(0,\lambda),
\]

then estimating the central scale and dividing by \(\sqrt\lambda\) repairs the tail as well, and
C1 is sufficient. **Median-\(\lambda\) is not inherently incapable of fixing tails.**

But a null can be near-normal in the middle and heavy in the tail, e.g.

\[
0.95\,N(0,1)+0.05\,N(0,4),
\]

for which \(\hat\lambda\approx1\) while \(P(|Z|>1.96)\) is badly inflated. The real blanks look
like this case rather than the scale case: the median-based \(\lambda\) for D1 is 1.08, which is
**not distinguishable from 1** at 205 blanks (the null 97.5th percentile of \(\hat\lambda\) is
about 1.35), yet the D1 blank false-positive rate is 9.27% against a 5% target
(\(p\approx0.005\)). Calibration accordingly moves it only to 8.78%. The same holds for the
parent arm (\(\lambda=1.22\), null upper 1.34, FP 8.5%, \(p=0.020\)); only the daughter arm has
a \(\lambda\) that clearly exceeds chance (1.43 against a null upper of 1.29).

### C2, concretely

From the usable blank statistics, estimate the rejection threshold directly:

\[
c_{0.05}=Q_{0.95}\!\left(|Z_{\text{blank}}|\right),
\qquad
\text{reject when } |Z_{\text{gene}}|>c_{0.05}.
\]

**Fix the quantile convention explicitly.** At \(B\approx200\) different interpolation rules
move the threshold enough to matter, so specify the order statistic rather than relying on a
default:

\[
c_{0.05}=|Z|_{(k)},
\qquad
k=\left\lceil 0.95\,(B+1)\right\rceil,
\]

where \(|Z|_{(1)}\le\dots\le|Z|_{(B)}\) are the sorted absolute blank statistics. Equivalently,
name the software convention used (e.g. R's `quantile(..., type = 1)`) and keep it constant
across every arm and replicate. This is reproducibility, not statistics — but an unstated
convention makes C2 unrepeatable.

This estimates the cut-point at the place calibration is actually required. Ablation C then
answers a sharper question than "does calibration help?":

> Is the residual inflation a **scale** problem, which median-\(\lambda\) repairs, or a
> **tail-shape** problem, which requires direct tail calibration?

A result of the form

| arm | Type I error |
|---|---:|
| C0 no calibration | 9% |
| C1 median-\(\lambda\) | 8.5% |
| C2 tail calibration | 5.2% |

would diagnose exactly why the real-data calibration only partly works.

**Keep C2 out of the primary method label.** Reporting it as a diagnostic — *"we additionally
evaluated direct tail calibration to determine whether residual inflation reflected
non-Gaussian empirical-null tails rather than simple variance inflation"* — resolves the
conceptual question without changing the method being benchmarked partway through its own
evaluation.

Note C2 is specified here for **Type I error**. Full empirical-null *p*-values,

\[
p_{\text{emp}}=\frac{1+\sum_b I\!\left(|Z_b|\ge|Z|\right)}{B+1},
\]

have resolution no finer than \(1/(B+1)\), so with a few hundred blanks they are too coarse to
drive FDR. Diagnosing the 5% Type I problem does not require solving that.

### Where \(\lambda\) comes from in the simulation

The real pipeline estimates \(\lambda\) from 50 blank barcodes run through the identical
analysis path — which yields roughly 200 *usable statistics* per calibration target, not 50
(see "Match the calibration information content" below; the barcode count is not the
calibration sample size). The simulation has no blanks, so one must be built; otherwise this
ablation has no \(\lambda\) to switch on and off.

**Simulate a blank set — in two regimes.** For each replicate, generate \(B_{\text{raw}}=50\)
additional "blank" barcodes, spread across enough strata/routes/layers to reach the usable-count
target below. Their abundance \(M_{gs}\) is drawn from the
**observed blank-barcode abundance distribution** rather than the gene distribution (real
blanks are sparser: median **16** molecules per slice against 64 for genes, measured after the
same \(SE_\delta<0.05\) and \(\ge30\)-cell gate — see §8).

*Ideal-null regime.* \(h_{bmsc}(r)=1\) in every slice and condition.

*Mouse-correlated regime.* **A blank with \(h=1\) in every slice is independent across slices
and therefore cannot inherit the dependence §10 is trying to detect.** If benchmark genes
follow \(\theta_{gmsc}=\theta_{gc}+b_{gm}+u_{gms}\) while blanks are fixed at
\(\theta_{\text{blank}}=0\), mouse nesting can inflate the gene statistics while the blanks
stay clean, giving \(\hat\lambda\approx1\) even when the tests are inflated. **This would make
Ablation C measure nothing.** So give each blank barcode a zero-mean mouse effect and no
condition effect:

\[
\theta_{bms}=b_{bm}+u_{bms},
\qquad
b_{bm}\sim N(0,\sigma_{\text{mouse}}^2),
\qquad
h_{bms}(r)=\exp(\theta_{bms}z).
\]

The blank remains a null *with respect to the condition contrast* — \(E(D_{\text{blank}})=0\) —
but slices from the same mouse are now correlated, so the blank-derived \(\lambda\) can
actually measure the inflation that mouse nesting produces.

**Blanks must traverse the identical analysis path**: the same \(SE_0<0.05\) precision gate,
the same minimum-slice requirement, the same REML/Knapp–Hartung pooling, and the same D1 and
arm-level tests. Computing 50 unrestricted p-values instead would estimate \(\lambda\) on a
different filtering regime from the one the genes experience.

Then estimate, **exactly as the real pipeline does** (`02d_stratified_merge.R:51-60`),
separately per statistic (D1, daughter arm, parent arm):

\[
\hat\lambda_{\text{raw}}
=
\frac{\operatorname{median}\!\left(z_{\text{blank}}^{2}\right)}
     {\operatorname{median}\!\left(\chi^{2}_{1}\right)}
=
\frac{\operatorname{median}\!\left(z_{\text{blank}}^{2}\right)}{0.4549364},
\qquad
\boxed{\ \hat\lambda_{\text{applied}}=\max\!\left(1,\ \hat\lambda_{\text{raw}}\right)\ }
\]

with the additional real-pipeline fallback that **fewer than 20 usable blank statistics yields
\(\hat\lambda=1\)** (no correction) rather than an estimate.

**The floor is not optional and not cosmetic.** Omitting it would benchmark a different
procedure from the published one, and in the anti-conservative direction: when
\(\hat\lambda_{\text{raw}}<1\), dividing by \(\sqrt{\hat\lambda}\) *inflates* \(z\), making the
test more liberal — the opposite of calibration's purpose. That is not a rare edge case. When
the truth is exactly \(\lambda=1\):

| usable blanks | P(floor fires) | \(\hat\lambda_{\text{raw}}\) | \(\hat\lambda_{\text{applied}}\) |
|---:|---:|---:|---:|
| 50 | 50% | 1.04 ± 0.34 | 1.15 ± 0.24 |
| 100 | 49% | 1.02 ± 0.23 | 1.10 ± 0.16 |
| 205 | 50% | 1.01 ± 0.16 | 1.07 ± 0.11 |
| 294 | 50% | 1.01 ± 0.14 | 1.06 ± 0.09 |

Two consequences the simulation must be able to see:

1. Unfloored, ~51% of runs would apply \(\hat\lambda<1\), inflating \(z\) by a median of 5.8%
   and by up to 40% in the worst draw.
2. Floored, the estimator is **biased upward by construction**:
   \(E[\hat\lambda_{\text{applied}}]\approx1.07\) at 205 blanks *even when the true
   \(\lambda\) is 1*. This is close to the observed D1 value of 1.08 — so D1's calibration is
   barely distinguishable from what pure sampling noise plus the floor would produce on their
   own. Ablation C should report this comparison explicitly.

**Always report both \(\hat\lambda_{\text{raw}}\) and \(\hat\lambda_{\text{applied}}\), plus the
floor activation rate,** in every replicate and at every \(B_{\text{usable}}\) in the
sensitivity sweep. The activation rate is what tells you whether the procedure is estimating
inflation or merely refusing to de-inflate.

### Match the calibration *information content*, not the barcode count

"\(B=50\) blanks" is the wrong target. The real analysis did not estimate \(\lambda\) from 50
numbers: 50 blank barcodes generate an eligible null statistic in **each** stratum, route and
layer they survive in, yielding **205 (D1), 212 (parent), 294 (daughter)** usable statistics
after the \(SE_0<0.05\) gate. A one-stratum pilot with 50 sparse blanks (median ~10 molecules)
would retain far fewer, and

\[
\hat\lambda=\frac{\operatorname{median}(z_b^{2})}{0.4549}
\]

is then dominated by its own sampling noise. When the truth is exactly \(\lambda=1\):

| usable blanks | \(\hat\lambda\) | 95% range |
|---:|---:|---:|
| 20 | 1.11 ± 0.56 | [0.33, 2.48] |
| 50 | 1.05 ± 0.35 | [0.50, 1.85] |
| 100 | 1.02 ± 0.24 | [0.60, 1.53] |
| 205 | 1.01 ± 0.17 | [0.71, 1.37] |

At 50 blanks, \(\hat\lambda\) ranges over [0.50, 1.85] by chance alone. Ablation C would then be
answering *"does randomly over- or under-estimating \(\lambda\) from a tiny sample help?"* rather
than *"does empirical-null calibration help?"*

**Requirement.** Simulate the same 50 conceptual blank barcodes, but across enough
strata/routes/layers that the number of usable statistics after filtering approximately matches
the real analysis (~200 per statistic). Record and report

\[
B_{\text{raw}},
\qquad
B_{\text{usable,D1}},
\qquad
B_{\text{usable,parent}},
\qquad
B_{\text{usable,daughter}}
\]

for every replicate — these, not the barcode count, are the calibration sample size.

**Calibration sample-size sensitivity.** Sweep

\[
B_{\text{usable}}=50,\ 100,\ 200,\ 300
\]

and report Type I error after calibration at each. This answers a question worth knowing on its
own: *how many null statistics does the median-\(\lambda\) procedure need before its own
estimation noise stops dominating?*

This matters because it makes Ablation C a test of **the whole calibration procedure** —
estimation included — rather than of the correction arithmetic applied to a \(\lambda\) handed
down by the simulator. It is also the only way to answer question 2 of §10: does \(\lambda\)
*recover* nominal error, or only partially?

Evaluate under both:

- **ideal null**, where slices are independent and \(\lambda\approx1\) is expected — here
  calibration should be nearly inert, and a large estimated \(\lambda\) would itself be a
  finding;
- **correlated null** (\(\sigma_{\text{mouse}}>0\), §10), where \(\lambda>1\) is expected —
  this is the regime the real data occupies.

A useful auxiliary output is the estimated \(\lambda\) itself, reported as a *function of*
\(\sigma_{\text{mouse}}\). Note this is **not** independent evidence that the simulation
matches reality — if \(\sigma_{\text{mouse}}\) were tuned to reproduce the observed
\(\lambda\), recovering that \(\lambda\) would be the calibration target, not a validation of
it. The defensible statement is that a particular \(\sigma_{\text{mouse}}\) places the
simulation in the same inflation regime as the data.

---

# 4. Simulation architecture

Use a **hybrid framework**.

## Simulation A: colon-anchored semi-simulation

Keep the observed colon spatial structure and simulate only gene expression.

Retain:

- real MorphoGAM \(r_i\) coordinates;
- real cell numbers;
- real slice-to-slice variation;
- real cell distributions along the tissue axis;
- observed Day 0 / Day 9 sample structure;
- realistic abundance and sparsity distributions;
- the observed mouse-to-slice nesting (31 slices from 12 mice).

This gives biological realism while preserving known simulated truth.

---

## Simulation B: simple fully synthetic validation

Use at least two coordinate distributions:

### Symmetric reference

\[
r_i\sim U(0,1).
\]

### Off-center / skewed stress test

Use a distribution with mean away from 0.5 and non-zero skew so that assumptions based on
axis symmetry are explicitly challenged. A convenient simple example is

\[
r_i\sim \mathrm{Beta}(2,2.5),
\]

which has mean approximately 0.444 and provides an off-center skewed geometry. If the goal is
to reproduce the observed IAF 4 moments more closely, resample the empirical IAF 4
coordinates rather than assuming a symmetric parametric distribution.

Purpose:

> Verify that observed performance is not specific to the geometry of the colon dataset, and explicitly test the skewed-coordinate setting in which a nominally symmetric spatial perturbation can acquire a non-zero first moment.

### Sizes

To be runnable, fix the synthetic design explicitly:

\[
10\ \text{slices}\ (5\ \text{per condition}),
\qquad
300\ \text{cells per slice},
\qquad
M=64\ \text{molecules per slice},
\]

matching the colon medians of §8, with 4 mice per arm to retain the nesting structure of §10.
This synthetic benchmark can otherwise be smaller than the colon-anchored benchmark: use the
same scenario set but only the primary effect size \(|D|=0.05\) unless a discrepancy with
Simulation A appears, in which case expand it to the full grid to localize the cause.

---

# 5. Gene-expression generative model

For gene \(g\), slice \(s\), cell \(i\), condition \(c\):

\[
Y_{gsi}\sim NB(\mu_{gsi},\phi_g).
\]

Standardize position within each slice:

\[
z_{si}
=
\frac{r_{si}-\bar r_s}{SD(r_s)}.
\]

Place transcripts with an arbitrary **positional weight function** \(h_c(r)\), normalized
within slice:

\[
\mu_{gsi}
=
M_{gs}L_{si}
\frac{h_c(r_{si})}
{\sum_j L_{sj}h_c(r_{sj})}.
\]

Here:

- \(M_{gs}\): expected total abundance of the gene in slice \(s\);
- \(L_{si}\): optional exposure / library-size term;
- \(\phi_g\): gene-specific dispersion;
- \(h_c(\cdot)\): the positional weight in condition \(c\).

**The normalization is the essential feature.** Summing over cells gives
\(\sum_i\mu_{gsi}=M_{gs}\) for *any* \(h_c\), so changing the spatial weight moves transcripts
**without changing total abundance**. Abundance is controlled solely by \(M_{gs}\), and
localization solely by \(h_c\). This is what makes the DE-only null clean.

### Choices of \(h_c\)

**Directional (monotone).** Standardize position within slice,

\[
z_{si}=\frac{r_{si}-\bar r_s}{SD(r_s)},
\qquad
h_c(r)=\exp(\theta_c z),
\]

which recovers the earlier model as a special case.

\[
\theta_c=0
\]

means no directional localization,

\[
\theta_c>0
\]

means bias toward increasing \(r\),

\[
\theta_c<0
\]

means bias toward decreasing \(r\).

**Non-directional first-moment-neutral change.** For the boundary scenario in §6, use

\[
h_9(r)
\propto
1+\kappa\,(r-c)^2,
\]

where the centering constant \(c\) is **not fixed at 0.5**. Instead, solve \(c\) numerically
using the actual coordinates in each simulated slice so that the expected
repositioning statistic is exactly neutral:

\[
\delta^{\mathrm{exp}}_9(c)
=
\frac{\sum_i \mu_{9i}(c)r_i}{\sum_i \mu_{9i}(c)}
-
\bar r
=
0.
\]

This is necessary because an off-center or skewed cell-position distribution can make
\((r-0.5)^2\) carry a non-zero first moment. Centering at the observed mean is a useful
initial value, but it is not guaranteed to make \(\delta=0\) exactly.

The \(\exp(\theta z)\) form is **monotone in \(r\) by construction and cannot express this
kind of non-directional profile change**, so the general \(h_c\) formulation is required.

**Other forms** (mid-axis peak, sigmoid, step) follow by substituting \(h_c\); nothing else in
the simulator changes.

Simulation truth should ultimately be summarized using

\[
\delta
=
\bar r_{\mathrm{transcripts}}
-
\bar r_{\mathrm{cells}}
\]

and

\[
D=\delta_9-\delta_0.
\]

## Truth calibration and validation

Do **not** assume that a chosen \(\theta\), \(\kappa\), or centering constant automatically
produces the desired \(\delta\). Before drawing negative-binomial counts, compute the
expected positional effect directly from the simulated mean counts:

\[
\delta^{\mathrm{exp}}_{sc}
=
\frac{\sum_i \mu_{sci}r_{si}}
     {\sum_i \mu_{sci}}
-
\bar r_s.
\]

### Directional scenarios

For a target displacement \(D^*\), use the first-order approximation

\[
\theta_{\mathrm{init}}
\approx
\frac{D^*}{SD(r)}
\]

as the starting value for a numerical root finder, and then solve using the exact expected
\(\delta\) above. In the observed IAF 4 geometry this approximation is already very close
(e.g. \(\theta\approx0.4\) is a useful starting value for \(|D|\approx0.05\)), but the
numerical solve should remain the final calibration step.

For the colon-anchored simulation, calibrate at the level at which truth is defined
(slice or stratum) and verify the pooled target after calibration.

### Non-directional scenario

For each chosen \(\kappa\), solve the centering constant \(c\) in

\[
h_9(r)=1+\kappa(r-c)^2
\]

so that

\[
\delta^{\mathrm{exp}}_9=0
\]

within numerical tolerance. For the main benchmark, solve this **per slice**, not only on a
pooled stratum, so that no residual slice-level first-moment signal remains. Do this using the
actual coordinate distribution rather than assuming symmetry around 0.5.

### Exposure safeguard

If the optional exposure term \(L_{si}\) is used, include it in the exact truth calculation.
Because \(\mu_{gsi}\propto L_{si}h(r_{si})\), even \(h(r)=1\) need not imply
\(\delta=0\) if exposure is associated with position. The cleanest primary benchmark is
therefore to set \(L_{si}=1\); empirical exposure can be added as a sensitivity analysis
with truth recalibrated using the full \(\mu_{gsi}\).

### Required truth checks

Before accepting a simulated parameter setting:

1. verify \(\delta_0^{\mathrm{exp}}\) matches its target;
2. verify \(\delta_9^{\mathrm{exp}}\) matches its target;
3. verify \(D^{\mathrm{exp}}=\delta_9^{\mathrm{exp}}-\delta_0^{\mathrm{exp}}\) matches the requested effect;
4. for the non-directional scenario, require \(|D^{\mathrm{exp}}|\) to be essentially zero;
5. store the achieved rather than nominal \(\delta\) and \(D\) as the simulation truth;
6. run the **comparator-truth audit** (§1's Benchmark A block): record the achieved
   \(h_9/h_0\) ratio swing across the observed \(z\) range and the maximum between-condition
   domain-proportion difference, and classify \(H_{0,\delta}\), \(H_{0,\mathrm{GAM}}\) and
   \(H_{0,\mathrm{DESpace2}}\) separately. A scenario intended as a shared null but showing a
   large ratio swing is **not** a shared null, and its comparator rejections must not be scored
   as Type I error.

A practical tolerance is

\[
|\delta^{\mathrm{exp}}-\delta^{\mathrm{target}}|<10^{-4}.
\]

---

# 6. Core simulation scenarios

## Null 1: unpositioned

\[
\delta_0=0,
\qquad
\delta_9=0.
\]

Purpose: basic Type I error.

---

## Null 2: stable positioned

\[
\delta_0=+a,
\qquad
\delta_9=+a.
\]

Purpose: verify that spatial localization itself is not mistaken for cross-condition repositioning.

---

## DE-only

**Explicit truth — defined on the \(\delta\) scale, not the \(\theta\) scale.**

\[
M_9=2M_0,
\qquad
\delta_0=\delta_9=+0.05,
\qquad
D=0,
\]

with \(\theta_0\) and \(\theta_9\) **calibrated separately** by the exact expected-\(\delta\)
solver of §5.

**Do not set \(\theta_0=\theta_9\).** Because \(\delta\approx\theta\,SD(r)\) and the two arms are
different cell populations with different coordinate distributions, an identical \(\theta\)
does not give an identical \(\delta\). This is the same trap fixed in Null 3, and it is live:
in the real data \(SD(r)=0.118\) for Fibro 6 at Day 0 against \(0.141\) for IAF 3 at Day 9, a
ratio of 1.19, so \(\theta=0.4\) in both arms would produce \(\delta_0=0.047\),
\(\delta_9=0.056\) and a **spurious \(D=+0.009\) — 18% of a nominal 0.05 effect, from geometry
alone**. (For Fibro 2 → IAF 4 the SDs happen to match to 0.3%, so the artifact is negligible
there; it is stratum-dependent and must be checked, not assumed.)

**Run both DE-only variants** — realistic-geometry and matched-geometry — as defined in §1's
Benchmark A block. The realistic variant supports the \(\delta\)-versus-pseudobulk specificity
claim; only the matched-geometry variant (\(h_0\equiv h_9\), identical coordinate distribution,
only \(M\) doubling) is a shared null suitable for the cross-method Type I comparison. The
truth defined above is the realistic variant.

A non-zero \(\delta\) is deliberate. With \(\delta_0=\delta_9=0\) this scenario would be
indistinguishable from Null 1 with a different \(M\), and would not test what it is meant to.
Holding a *real* gradient fixed while doubling abundance is the case that separates the two
questions.

Purpose: demonstrate separation of abundance and spatial localization.

Expected — **split by variant**:

*Matched-geometry DE-only* (the shared null):

- pseudobulk DE: significant;
- proposed method, GAM, DESpace2: not significant.

*Realistic-geometry DE-only* (\(\theta\) solved per arm, populations differ):

- pseudobulk DE: significant;
- proposed method: not significant (\(D_\delta=0\) by construction);
- GAM, DESpace2: **determined by the truth audit, not assumed.** Where
  \(SD(r_0)\neq SD(r_9)\) forces \(\theta_0\neq\theta_9\), \(h_0\not\equiv h_9\) and their nulls
  are false, so rejections there are correct behaviour rather than Type I error.

---

## Null 3: population movement without transcript repositioning

The **cell population** shifts along the axis between conditions while the gene's transcript
bias *relative to its own cells* stays fixed:

\[
E[r_{\text{cells},0}]\approx0.35,
\qquad
E[r_{\text{cells},9}]\approx0.55,
\]

but define the null by

\[
\delta_0=\delta_9=a
\]

rather than by \(\theta_0=\theta_9\).

This distinction matters because

\[
\delta\approx\theta\,SD(r),
\]

so resampling the Day-9 cell population can change \(SD(r)\), and the same \(\theta\) can
therefore produce a different \(\delta\).

Implementation:

1. choose a baseline target such as \(a=0.05\);
2. resample Day-9 cells from the upper part of the observed \(r\) distribution so the
   population genuinely occupies a different part of the same axis;
3. calibrate \(\theta_0\) and \(\theta_9\) separately using the exact expected-\(\delta\)
   solver so that both conditions have the same population-relative displacement \(a\);
4. verify \(D^{\mathrm{exp}}=0\) before count sampling.

Implementation note: \(r\) is fitted per slice and normalized to [0,1], so one cannot simply
"add a constant" because the coordinate would renormalize.

### The resampling must shift the mean without shrinking the variance

Hard truncation — "take cells with \(r\) above some cut" — moves the mean but also **truncates
the distribution**, reducing \(\operatorname{var}(r)\). That is a confound, because the
positional standard error depends on it directly:

\[
SE_{0}
=
\sqrt{
\operatorname{var}(r)
\left(
\frac{1}{N_{\mathrm{eff}}}-\frac{1}{N}
\right)}.
\]

A truncated Day-9 arm would therefore be systematically **more precise** than Day 0 for a
reason that has nothing to do with repositioning, and any method could look better or worse at
Day 9 as an artifact of the resampling scheme rather than of the biology.

Use **moment-balanced weighted resampling**. A one-parameter exponential tilt
\(\Pr(i)\propto\exp(\gamma r_i)\) controls only the mean and does **not** independently control
\(\operatorname{Var}(r)\), so it cannot be claimed to preserve it. Use two parameters,

\[
\Pr(i)\propto\exp(\gamma_1 r_i+\gamma_2 r_i^2),
\]

and solve \(\gamma_1,\gamma_2\) jointly against both moments:

\[
\frac{\sum_i w_i r_i}{\sum_i w_i}=\mu_9^{*},
\qquad
\frac{\sum_i w_i r_i^{2}}{\sum_i w_i}=\sigma_0^{2}+(\mu_9^{*})^{2}.
\]

This deliberately changes \(\bar r\) while holding \(SD(r)\) fixed. A one-parameter tilt may be
used **only** if the achieved SD passes the prespecified tolerance below.

If the empirical coordinate support makes the requested \(0.35\to0.55\) shift impossible while
preserving the variance, **reduce the target shift** rather than accepting a large variance
change.

### Required checks for this scenario

Report, per condition, the achieved

\[
\bar r_{\text{cells}},
\qquad
SD(r),
\qquad
N,
\qquad
\delta^{\mathrm{exp}},
\]

and require that \(SD(r)\) match between conditions to within a stated tolerance
(a few percent). If it cannot be matched, report the residual difference so that any
performance gap can be attributed correctly rather than read as a repositioning effect.

Purpose:

> This is the **only scenario that tests the rationale for subtracting** \(\bar r_{\text{cells}}\).

Expected:

- proposed method (\(\delta=\bar r_{\text{transcripts}}-\bar r_{\text{cells}}\)): **not significant**;
- any statistic using **absolute transcript position** while claiming to measure
  within-population displacement: **significant — and that is a genuine false positive**,
  because the quantity it claims to estimate has not changed;
- GAM and DESpace2: likely significant — **but this is not a false positive for them.** Their
  target is the marginal spatial profile, and the marginal profile genuinely did change. They
  are answering a different question correctly.

**This distinction must be stated wherever the panel is shown.** Reporting "the comparators
fire and our method does not" without it would misrepresent DESpace2 and the GAM, and a
reviewer would be right to object. The scenario demonstrates that \(\delta\) targets a
*population-relative* quantity, not that the other methods are miscalibrated.

Without this null there is no test defending the population-relative design, and the obvious
reviewer question — *why not just use transcript position?* — has no figure to answer it.

---

## Gained

\[
\delta_0=0,
\qquad
\delta_9=+a.
\]

Purpose: test acquisition of directional localization.

---

## Lost

\[
\delta_0=+a,
\qquad
\delta_9=0.
\]

Purpose: test disappearance of directional localization.

---

## Reversed

\[
\delta_0=-a,
\qquad
\delta_9=+a.
\]

Purpose: test switching from one side of the population center to the other.

---

## Amplified

\[
0<\delta_0<\delta_9.
\]

Anchored form (§7), shown at \(|D|=0.05\):

\[
\delta_0=0.025,
\qquad
\delta_9=0.025+|D|=0.075.
\]

The baseline \(\delta_0=0.025\) is held **fixed across the whole \(|D|\) grid**; only
\(\delta_9\) moves. See §7 for why the anchor is required.

Purpose: test strengthening of a same-direction spatial bias.

---

## Dampened

\[
0<\delta_9<\delta_0.
\]

Anchored form (§7), shown at \(|D|=0.05\):

\[
\delta_0=0.025+|D|=0.075,
\qquad
\delta_9=0.025.
\]

The endpoint \(\delta_9=0.025\) is held **fixed across the whole \(|D|\) grid**; only
\(\delta_0\) moves. See §7 for why the anchor is required.

Purpose: test weakening of a same-direction spatial bias.

---

## Non-directional spatial-pattern change

Use

\[
h_0(r)=1,
\]

and

\[
h_9(r)
=
1+\kappa(r-c)^2,
\]

where \(c\) is solved numerically for each simulated slice so that

\[
\delta_9^{\mathrm{exp}}=0.
\]

Under the primary clean simulation with \(L_i=1\),

\[
\delta_0^{\mathrm{exp}}=0,
\qquad
\delta_9^{\mathrm{exp}}=0,
\qquad
D^{\mathrm{exp}}=0
\]

by construction, even when the coordinate distribution is off-center or skewed.

Use more than one \(\kappa\) value so that the profile-shape change can become progressively
stronger while the first moment remains neutral. The achieved \(\delta\) must be checked
numerically for every \(\kappa\).

Purpose:

> Establish the boundary of the repositioning method using a profile change that is genuinely first-moment neutral.

Expected:

- GAM may have high power;
- DESpace2 may detect the domain-level change;
- the \(\delta\)-based method should remain near the null rejection rate if the first moment is truly unchanged.

This is an important positive control for the broader spatial methods.

---

# 7. Effect sizes

## \(a\) is not \(|D|\) — set \(a\) per pattern

The scenarios in §6 are written in terms of a per-condition displacement \(a\), but the power
comparison is made at a common \(|D|\). These are not the same number, and the mapping differs
by pattern:

| pattern | \(\delta_0\) | \(\delta_9\) | \(\|D\|\) | \(a\) needed for \(\|D\|=0.05\) |
|---|---:|---:|---:|---:|
| gained | \(0\) | \(+a\) | \(a\) | 0.05 |
| lost | \(+a\) | \(0\) | \(a\) | 0.05 |
| **reversed** | \(-a\) | \(+a\) | \(2a\) | **0.025** |
| amplified | \(+0.025\) | \(+(0.025+\|D\|)\) | \(\|D\|\) | anchored, see below |
| dampened | \(+(0.025+\|D\|)\) | \(+0.025\) | \(\|D\|\) | anchored, see below |

**Amplified and dampened need a fixed baseline anchor, not just \(b-a=\|D\|\).** That condition
alone admits infinitely many simulations — \(0.01\to0.06\) and \(0.10\to0.15\) both give
\(\|D\|=0.05\) but differ in power and in classification difficulty, because the pattern label
depends on whether each endpoint CI excludes zero, not on \(D\) alone. Fix the anchor at
\(\delta_0=0.025\) (amplified) and \(\delta_9=0.025\) (dampened) across the whole effect-size
grid:

\[
\text{amplified: }\ \delta_0=0.025,\ \ \delta_9=0.025+|D|,
\qquad
\text{dampened: }\ \delta_0=0.025+|D|,\ \ \delta_9=0.025 .
\]

**Using \(a=0.05\) for `reversed` would test it at \(|D|=0.10\), double every other pattern**,
and Panel D would then compare power across patterns at unequal effect sizes. Always set \(a\)
(or \(b-a\)) from the target \(|D|\), and record the achieved \(|D|\) per scenario as required
by §5.


Use several true repositioning magnitudes:

\[
|D|
=
0.01,\ 0.025,\ 0.05,\ 0.075.
\]

Interpretation:

- 0.01: very weak;
- 0.025: weak;
- 0.05: moderate;
- 0.075: strong.

The primary power figure should show

\[
\text{Power}
\quad\text{versus}\quad
|D|.
\]

For each target \(|D|\), initialize the directional parameter with

\[
\theta_{\mathrm{init}}
=
\frac{D}{SD(r)}
\]

and then solve numerically for the exact target using \(\delta^{\mathrm{exp}}\). The shortcut
is used only to initialize the solver; the achieved \(D\), not the nominal \(\theta\), is
the effect-size truth stored for benchmarking.

---

# 8. Gene abundance

**v13.** All three strata are run in Experiment 7, and a per-arm abundance scenario
(\(M_B/M_A\in\{1/3,1,3\}\)) is added as Experiment 7b — see §0a A3. In v12 only the median
stratum was run.

Use empirical colon abundance strata, taken from the real per-slice molecule counts **after**
the analysis gate (\(SE_\delta<0.05\), \(\ge30\) cells; n = 34,126 gene slice-observations):

| stratum | molecules per slice | source |
|---|---:|---|
| low | 27 | 25th percentile |
| medium | 64 | median |
| high | 331 | 90th percentile |

For reference, blank barcodes under the same gate have a median of **16** molecules per slice
(n = 3,987), which is the distribution the simulated blanks of §3 should be drawn from.

Purpose:

> Evaluate how sparsity affects power and whether performance differences are strongest for low-information genes.

---

# 9. Between-slice heterogeneity

Simulate slice-specific directional effects on the \(\theta\) scale:

\[
\theta_{sc}
=
\theta_c+u_s,
\]

\[
u_s\sim N(0,\tau_\theta^2).
\]

Because

\[
\delta\approx\theta\,SD(r),
\]

heterogeneity on the \(\theta\) scale approximately maps to heterogeneity on the
repositioning scale as

\[
\tau_\delta
\approx
\tau_\theta\,SD(r).
\]

Therefore \(\tau_\theta\) should not be chosen arbitrarily. Calibrate it so that the
resulting \(\tau_\delta\) produces heterogeneity comparable to the observed range of
\(I^2\) in the real analysis.

Using a representative within-slice variance \(\sigma_\delta^2\),

\[
I^2
\approx
\frac{\tau_\delta^2}
{\tau_\delta^2+\sigma_\delta^2},
\]

so a target \(I^2\) gives the approximate starting value

\[
\tau_\delta
\approx
\sigma_\delta
\sqrt{\frac{I^2}{1-I^2}},
\qquad
\tau_\theta
\approx
\frac{\tau_\delta}{SD(r)}.
\]

### Targets, specified

Observed \(I^2\) in the real analysis (n = 1,984 parent / 1,936 daughter contrasts):

| arm | median | IQR | 90th pct | max |
|---|---:|---:|---:|---:|
| parent (D0) | 13.5% | [0, 44.0] | 59.7% | 93.9% |
| daughter (D9) | 14.9% | [0, 43.1] | 61.6% | 93.5% |

Use four levels:

\[
I^2\in\{0\%,\ 15\%,\ 45\%,\ 60\%\}
\]

**v13.** On the matched substrate the REML \(I^2\) estimator has a floor of ~11% at
\(\tau_\theta=0\) (k = 5), so the 0% and 15% targets are both realised by \(\tau_\theta=0\),
which is therefore the *typical* real case; \(\tau_\theta=0.15\) realises ~45–50% and 0.30
~79%. Every figure reports realised \(I^2\) beside \(\tau_\theta\), and the map is re-derived
at k = 10 (§0a A5).

corresponding to none, **typical** (the observed median), upper-quartile, and high (the observed
90th percentile). Note the real median is only ~14%, so 45% and 60% are stress tests rather than
the common case — do not present them as typical. The two arms are close enough here that a
single set of targets serves both; state that this was checked rather than assumed.

**State which arm the target is matched to.** \(I^2\) is estimated per arm in the real
analysis and the two differ systematically — observed values span roughly 2--68%, with the
daughter arm generally showing higher estimated \(I^2\) than the parent arm. Calibrating
\(\tau_\theta\) against a single pooled \(I^2\) would set heterogeneity too high for the
parent arm and too low for the daughter arm. Either calibrate the two arms separately, or
state explicitly which arm the quoted target refers to.

As with the directional effect-size calibration, use these formulas to initialize the
simulation and verify the achieved heterogeneity on the \(\delta\) scale.

Purpose:

> Evaluate when the random-effects component improves performance under heterogeneity matched to the scale actually estimated by the method.

---

# 10. Mouse-level correlation — MAIN colon robustness analysis

**This is not an optional add-on.** In the real colon data 31 slices come from only 12 mice,
and the blank-derived inflation factors are \(\lambda=1.08\!-\!1.43\) (D1, daughter arm,
parent arm) — i.e. **the empirical null is measurably inflated**. Mouse-level dependence is one
plausible contributor, but the observed \(\lambda\) **does not identify its cause**; detection
artifacts, residual slice effects and model misspecification are equally consistent with it.

The simulation therefore runs the causal arrow in the *estimable* direction: vary
\(\sigma_{\text{mouse}}\) across a prespecified grid and measure the resulting blank-derived
\(\hat\lambda\),

\[
\sigma_{\text{mouse}}
\;\longrightarrow\;
\text{simulated dependence}
\;\longrightarrow\;
\hat\lambda,
\]

rather than tuning \(\sigma_{\text{mouse}}\) to reproduce the observed \(\lambda\) and then
citing that agreement as evidence. Values of \(\sigma_{\text{mouse}}\) producing \(\hat\lambda\)
in the observed range identify a **biologically relevant simulation regime** — which is a
weaker and defensible claim.

Simulating only \(\sigma_{\text{mouse}}=0\) would place every method in a regime the real data
does not occupy and would make all of them look well calibrated.

The mouse-nested condition is therefore run for **every core scenario in §6**, not as a
side analysis. Simulate

\[
\theta_{msc}
=
\theta_c+b_m+u_{ms},
\]

where

\[
b_m\sim N(0,\sigma_{\mathrm{mouse}}^2).
\]

Compare:

\[
\sigma_{\mathrm{mouse}}=0
\]

with

\[
\sigma_{\mathrm{mouse}}>0.
\]

Purpose:

> Determine whether treating slices as independent inflates Type I error and whether a hierarchical animal-level model is required.

### The grid, specified

\(\sigma_{\text{mouse}}\) is on the \(\theta\) scale; the induced displacement-scale mouse SD is
\(\sigma_{\text{mouse}}\cdot SD(r)\), with \(SD(r)\approx0.13\) in the colon strata:

| \(\sigma_{\text{mouse}}\) (\(\theta\) scale) | induced \(\delta\)-scale mouse SD | regime |
|---:|---:|---|
| 0 | 0 | independent slices; reference |
| 0.10 | 0.013 | mild |
| 0.20 | 0.026 | moderate |
| 0.40 | 0.052 | strong — comparable to a \(\|D\|=0.05\) effect |

Run the full grid for Null 1 and Null 2; run \(\{0,\ 0.20\}\) for the remaining §6 scenarios to
contain cost. **Extend the grid upward if \(\hat\lambda\) has not yet bracketed the observed
1.08–1.43 range** — the point is to locate that regime, not to assume 0.40 reaches it.

Note this is a **prespecified grid, not a value matched to the observed \(\lambda\)**: tuning
\(\sigma_{\text{mouse}}\) to reproduce \(\lambda\) and then citing the agreement as evidence is
the circularity removed in v5. Report \(\hat\lambda_{\text{raw}}\),
\(\hat\lambda_{\text{applied}}\), the floor activation rate and Type I error at every grid
point.

Three specific questions this answers:

1. Does treating slices as independent inflate Type I error, and by how much?
2. Does the \(\lambda\) recalibration (Ablation C) actually recover nominal error, or only
   partially — the real blanks residual is 8.8%, not 5% — and **is the shortfall a scale
   problem or a tail-shape problem?** Arm C2 is what distinguishes these: if median-\(\lambda\)
   and tail calibration both reach ~5%, the inflation was pure scale; if only C2 does, the
   empirical null is non-Gaussian in shape and the current estimator is measuring the wrong
   feature of it.
3. Would an explicit mouse random effect \(\delta_{ms}=\mu+b_m+\epsilon_{ms}\) do better, and
   is it estimable at 4 and 8 mice per arm?

Question 3 is the acknowledged weak point of the current framework; the simulation is the
only place it can be answered, since the real data has too few mice to fit the model reliably.

---

# 11. Performance metrics

## Type I error

\[
\widehat{\alpha}
=
\frac{FP}{FP+TN}.
\]

Evaluate especially under:

- null unpositioned;
- null stable positioned;
- null population-movement (Null 3);
- DE-only.

Type I error is defined **per method against that method's own null hypothesis.** For the
**matched-substrate versions** of Null 1, Null 2 and DE-only the spatial null is shared across
methods. Under **real Day 0 / Day 9 geometry**, null status is not assumed — it is determined
per scenario and per stratum by the comparator-truth audit (§1), and Null 2 and DE-only are
frequently *not* shared nulls there. For Null 3 only
\(H_{0,\delta}\) is true; \(H_{0,\mathrm{GAM}}\) and \(H_{0,\mathrm{DESpace2}}\) are false, so
those methods' Null-3 rejection rates are not Type I errors and must not be pooled into a Type I
error figure.

Always accompany \(\hat\alpha\) with the **cluster-robust Monte Carlo SE and the achieved
\(n_{\text{eff}}\)** defined in the Objective's Benchmark A block. A Type I error reported
without its Monte Carlo error cannot support any of the calibration claims in §10, since at
\(R=100\) with one null gene the 95% interval is [0.007, 0.093] — wide enough to contain both
the nominal 5% and the observed 8.8% residual.

---

## Power

\[
\mathrm{Power}
=
\frac{TP}{TP+FN}.
\]

This is the main sensitivity metric.

---

## Type II error

\[
\beta
=
1-\mathrm{Power}.
\]

---

## False discovery rate

For simulations containing both null and alternative genes,

\[
FDR
=
E
\left[
\frac{FP}{TP+FP}
\right].
\]

Evaluate after each method's multiple-testing correction.

**FDR is defined only in Benchmark C**, because it depends on the null fraction \(\pi_0=G_0/G\)
and is meaningless without it. Type I error above is *not* affected: it is measured under
pure-null scenarios in Benchmark A, where every gene is null by construction and no \(\pi_0\)
is required. Do not report an FDR for Benchmark A or B.

---

## Pattern-classification accuracy

For the proposed method only, evaluate whether detected genes are correctly labeled as:

- gained;
- lost;
- reversed;
- amplified;
- dampened.

Do not compare this directly against methods that do not produce these labels.

---

# 12. Recommended pilot

Before running the full benchmark:

### Use one colon stratum

Choose one well-populated mucosal fibroblast transition.

### Use six scenarios

- null unpositioned;
- stable positioned;
- gained;
- lost;
- reversed;
- amplified.

### Use one effect size

\[
|D|=0.05.
\]

Set the per-condition displacement \(a\) from this target using the mapping table in §7 —
it is **not** \(a=0.05\) for every pattern. In particular `reversed` uses \(a=0.025\), since
\(\delta_0=-a,\ \delta_9=+a\) gives \(|D|=2a\). Otherwise Panel D would compare power across
patterns at unequal effect sizes.

### Use one abundance level

Median abundance.

### Size each pilot arm for its own question

The pilot is for **debugging the simulator**, not for calibration claims. Size accordingly:

| pilot arm | scenarios | replicates | genes/replicate | purpose |
|---|---|---:|---:|---|
| A-pilot (calibration) | null unpositioned, stable positioned | \(R=200\) | \(G_0=10\) | check \(\hat\alpha\), \(n_{\text{eff}}\), \(\hat\lambda\) stability |
| B-pilot (power) | gained, lost, reversed, amplified | \(R=100\) | 1 | check power ordering is sensible |
| C-pilot (FDR) | mixed | \(R=50\) | \(G=500,\ \pi_0=0.8\) | check the FDR machinery runs end to end |

\(R=100\) with one gene is adequate for the B-pilot, where effects are large, and **not**
adequate for any statement about Type I error — see §11.

### Blank set — the pilot cannot reach the target the way the full benchmark does

There is a real tension between "use one colon stratum" and
"\(B_{\text{usable}}\approx200\) per statistic" (§3). In the **full** benchmark the target is
met for free: 50 barcodes surviving across many strata, routes and layers generate ~200 usable
statistics, exactly as in the real pipeline. In a **one-stratum pilot** 50 barcodes yield at
most 50 raw statistics and fewer after the \(SE_0<0.05\) gate, which is squarely inside the
regime where \(\hat\lambda\) spans [0.50, 1.85] by chance.

Resolve it one of two ways, and state which was used:

1. **Raise \(B_{\text{raw}}\) in the pilot** to roughly 500 barcodes within the single stratum,
   so ~200 survive filtering. Blanks are cheap — they carry no effect to calibrate — so this
   costs far less than adding strata. *This is the recommended pilot route.*
2. **Simulate blanks across all strata** even while genes are simulated in one, reproducing the
   real pipeline's pooling structure.

Do not simply accept whatever count 50 barcodes happens to give: confirm the achieved
\(B_{\text{usable}}\) per statistic before trusting any Ablation C output, and report it
alongside every \(\hat\lambda\).

### Run

1. full repositioning method;
2. GAM interaction model;
3. DESpace2;
4. pseudobulk DE for the relevant abundance controls;
5. three ablations — including all three calibration arms C0 / C1 / C2.

### Pilot acceptance checks

Do not proceed to the full grid until:

0. **every comparator has passed the §1b calibration gate under Null 1, run on the matched
   pseudo-arm substrate** — Type I error near nominal under the chosen specification, or an explicit decision to report matched-threshold
   power. This is check zero because a comparator that fails it invalidates Panel C regardless
   of anything else in the pilot;
1. every \(\delta^{\mathrm{exp}}\) and \(D^{\mathrm{exp}}\) matches its target to \(10^{-4}\) (§5),
   **and the comparator-truth audit has been run and its per-method null classification
   recorded** for every scenario × stratum in the pilot;
2. Null 3 achieves matched \(SD(r)\) between conditions within tolerance (§6);
3. \(B_{\text{usable}}\) reaches ~200 per statistic, and \(\hat\lambda\) is stable across replicates;
4. \(n_{\text{eff}}\ge1{,}000\) in the A-pilot, or \(R\) is raised until it is;
5. achieved \(N_{\mathrm{event}}\) and \(N_{\mathrm{eff}}\) hit the Ablation A targets (~80 vs ~5);
6. Under the ideal null (\(\sigma_{\text{mouse}}=0\)), the simulator reproduces the intended
   data-generating truth exactly (checks 1-2), and the uncalibrated procedure is assessed for
   approximate nominality. **Unexpected inflation is a diagnostic trigger, not an attribution.**
   It should prompt checks of both the simulator *and* the inferential procedure, since even at
   \(\sigma_{\text{mouse}}=0\) the method still involves a permutation \(SE_0\) used as a
   sampling SE, Knapp-Hartung inference at small \(k\), REML variance estimation, the
   \(SE_0<0.05\) selection step, and NB count generation. An uncalibrated 6-7% here could be a
   genuine inferential limitation worth reporting, not a coding error.

---

# 12b. Reproducibility requirements

Every reported number must be regenerable from the recorded inputs alone.

- **Seeds.** Set one master seed per replicate, derived deterministically as
  `seed = f(scenario_id, sigma_mouse, |D|, abundance, replicate_index)`, and record it in the
  output. Do not rely on a single global seed set once at the top of a run — that makes any
  single replicate impossible to reproduce in isolation for debugging.
- **Store achieved truth, not nominal truth.** Per §5 the achieved \(\delta_0^{\exp}\),
  \(\delta_9^{\exp}\), \(D^{\exp}\), and per §3 the achieved \(N_{\mathrm{event}}\),
  \(N_{\mathrm{eff}}\), \(B_{\text{usable}}\), \(\hat\lambda_{\text{raw}}\),
  \(\hat\lambda_{\text{applied}}\). Benchmarking against nominal parameters when the achieved
  values differ is the most common way a simulation quietly misreports its own effect sizes.
- **Record the exact comparator formulas and package versions** (mgcv, DESpace2, edgeR/DESeq2),
  since §1b makes the comparator specification a substantive result, not an implementation
  detail.
- **Log every filtered-out unit** — slices dropped by the \(SE_0<0.05\) gate or the cell
  minimum — so that selection effects are visible rather than silent.

---

# 13. Final benchmark figure

## Panel A — Simulation design

Show example Day 0 / Day 9 spatial profiles for:

- null;
- gained;
- lost;
- reversed;
- amplified;
- non-directional change.

---

## Panel B — Type I error

From Benchmark A. Compare:

- full method;
- GAM;
- DESpace2;
- ablations.

Reference line:

\[
\alpha=0.05.
\]

Plot the cluster-robust Monte Carlo interval on every bar. No \(\pi_0\) is involved: rejection
rates here are measured under scenarios null for the repositioning target.

**Label Null 3 separately.** GAM and DESpace2 are not null under Null 3 (the marginal profile
genuinely changes), so their bars there are **rejection rates, not Type I error**. Either split
Null 3 into its own annotated bar group or drop the two profile methods from it. Mislabelling
those bars would assert that DESpace2 and the GAM are miscalibrated when they are answering a
different question correctly.

---

## Panel C — Power versus true repositioning magnitude

\[
x=|D|,
\qquad
y=\mathrm{Power}.
\]

Compare:

- proposed method;
- GAM;
- DESpace2.

**Plot power at matched empirical Type I error** (thresholds calibrated on Null 1, §1b), with
the nominal-threshold curves shown as a secondary panel or supplement. Annotate each method's
native Null-1 error rate on the figure, so a reader can see what the matching corrected for.
Comparing power across methods running at different real false-positive rates is not a
comparison.

---

## Panel D — Power by repositioning pattern

Compare performance for:

- gained;
- lost;
- reversed;
- amplified;
- dampened.

---

## Panel E — DE versus repositioning specificity

Show the two complementary scenarios:

1. DE-only;
2. repositioning-only.

Compare:

- pseudobulk DE;
- proposed repositioning method.

---

## Panel F — Ablation study

Show how performance changes for:

- naive SE (Ablation A), plotted against achieved \(N_{\mathrm{eff}}\) at fixed
  \(N_{\mathrm{event}}=100\), since that is the contrast the ablation constructs;
- no random-effects pooling (Ablation B);
- calibration arms C0 / C1 / C2 (Ablation C), with C2 labelled a diagnostic rather than a
  competing method.

Panel F shows Type I error and power from Benchmark A/B. Add FDR here **only** if it is drawn
from Benchmark C, where \(\pi_0\) is defined.

---

## Panel G — Non-directional spatial change

Show that the GAM and/or DESpace2 can outperform the repositioning statistic when the profile changes without a directional center shift.

---

# 14. Main statistical hypotheses

## Primary hypothesis

> A statistic explicitly targeting transcript displacement along a continuous biological axis should provide higher power for directional repositioning than more general full-profile or domain-based spatial models, while maintaining calibrated Type I error and FDR.

## Secondary hypothesis

> The continuous GAM should outperform the repositioning statistic for non-directional spatial-profile changes that leave the first spatial moment unchanged.

## Ablation hypothesis

> Effective positional uncertainty, random-effects pooling, and empirical-null calibration each contribute measurably to either power or error control under the conditions they were designed to address.

---

# 15. Final methodological framing

The central comparison is:

| Method | Spatial representation | Main inferential target |
|---|---|---|
| Proposed method | Continuous first moment \(\delta\) | Directional positional shift |
| GAM interaction | Continuous smooth \(s(r)\) | Any continuous profile change |
| DESpace2 | Discrete spatial domains | Any domain-level pattern change |
| Pseudobulk DE | No spatial information | Total abundance change |

The goal is not to claim universal superiority.

The intended conclusion is:

> **The repositioning framework is a targeted method for directional shifts along an ordered tissue axis. Its reduced spatial representation may improve sensitivity for gain, loss, reversal, amplification, and dampening of positional bias, while full-profile approaches remain preferable for complex non-directional spatial changes.**
