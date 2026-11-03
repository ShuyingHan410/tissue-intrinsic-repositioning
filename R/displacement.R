#' Population-relative displacement of a gene in one section
#'
#' Computes \eqn{\delta = \sum_i n_i r_i / \sum_i n_i - N^{-1} \sum_i r_i}: the
#' transcript-weighted mean position of a gene minus the mean position of all cells of
#' the host population in one section and layer. Positive values mean the transcripts sit
#' deeper (towards \eqn{r = 1}) than the population itself. The analytic null standard
#' error assumes position–expression exchangeability:
#' \eqn{SE_0 = \widehat{sd}(r)\,\sqrt{\sum_i n_i^2 / (\sum_i n_i)^2 - 1/N}}.
#'
#' @param r numeric vector, tissue-intrinsic position of every cell of the population
#'   in the section (0 = inner / surface, 1 = outer / deep). Cells with zero counts must
#'   be included.
#' @param n non-negative numeric vector of transcript counts of the gene, same length as `r`.
#' @param min_cells,se_max eligibility gate: a section is eligible when it has at least
#'   `min_cells` cells and `0 < SE_0 < se_max` (defaults are the values used in the paper).
#' @return a one-row data.frame with `cells`, `molecules`, `n_eff` (effective number of
#'   transcripts, \eqn{(\sum n)^2/\sum n^2}), `delta`, `se0`, `z` and `eligible`.
#' @examples
#' set.seed(1); r <- runif(500); n <- rpois(500, 0.3 * exp(2 * r))
#' section_delta(r, n)
#' @export
section_delta <- function(r, n, min_cells = 30, se_max = 0.05) {
  ok <- is.finite(r) & is.finite(n)
  r <- r[ok]; n <- n[ok]
  N <- length(r); S <- sum(n)
  delta <- if (S > 0) sum(n * r) / S - mean(r) else NA_real_
  sw2   <- if (S > 0) sum(n^2) / S^2 - 1 / N else NA_real_
  se0   <- if (is.finite(sw2) && sw2 > 0 && N > 1) stats::sd(r) * sqrt(sw2) else NA_real_
  n_eff <- if (S > 0) S^2 / sum(n^2) else NA_real_
  eligible <- N >= min_cells && is.finite(se0) && se0 > 0 && se0 < se_max
  data.frame(cells = N, molecules = S, n_eff = n_eff, delta = delta, se0 = se0,
             z = if (is.finite(se0) && se0 > 0) delta / se0 else NA_real_,
             eligible = eligible)
}

#' Per-section displacement table for one gene
#'
#' Applies [section_delta()] to every section of a cell-level table.
#'
#' @param data data.frame with one row per cell of the host population.
#' @param gene name of the column holding the transcript counts of the gene.
#' @param r name of the position column (default `"r"`).
#' @param section name of the section identifier column (default `"section"`).
#' @param keep names of additional columns that are constant within a section and should
#'   be carried over (for example `"condition"` or `"mouse"`).
#' @inheritParams section_delta
#' @return a data.frame with one row per section.
#' @examples
#' f <- system.file("extdata", "igfbp5_fibro2_iaf4.csv", package = "TICR")
#' d <- read.csv(f)
#' head(delta_by_section(d, "Igfbp5", keep = c("condition", "mouse", "population")))
#' @export
delta_by_section <- function(data, gene, r = "r", section = "section",
                             keep = c("condition"), min_cells = 30, se_max = 0.05) {
  stopifnot(all(c(gene, r, section) %in% names(data)))
  keep <- intersect(keep, names(data))
  out <- lapply(split(data, data[[section]]), function(d) {
    s <- section_delta(d[[r]], d[[gene]], min_cells = min_cells, se_max = se_max)
    meta <- lapply(keep, function(k) d[[k]][1]); names(meta) <- keep
    cbind(data.frame(section = d[[section]][1], stringsAsFactors = FALSE),
          as.data.frame(meta, stringsAsFactors = FALSE), s)
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  out$gene <- gene
  out
}

#' Pool per-section displacement estimates within one condition
#'
#' Random-effects meta-analysis of the eligible per-section estimates by REML with the
#' Knapp–Hartung adjustment (`metafor::rma(..., method = "REML", test = "knha")`), falling
#' back to DerSimonian–Laird if REML fails to converge.
#'
#' @param delta,se numeric vectors of per-section estimates and their standard errors.
#' @param min_k minimum number of sections required (default 3).
#' @return a one-row data.frame with `k`, `estimate`, `se`, `ci_lo`, `ci_hi`, `tau2`,
#'   `I2`, `p` (test of estimate = 0), or `NULL` when fewer than `min_k` sections are
#'   available.
#' @examples
#' pool_sections(c(0.02, 0.03, 0.01, 0.04), c(0.01, 0.012, 0.009, 0.011))
#' @export
pool_sections <- function(delta, se, min_k = 3) {
  ok <- is.finite(delta) & is.finite(se) & se > 0
  y <- delta[ok]; s <- se[ok]; k <- length(y)
  if (k < min_k) return(NULL)
  fit <- tryCatch(metafor::rma(yi = y, sei = s, method = "REML", test = "knha"),
                  error = function(e) tryCatch(metafor::rma(yi = y, sei = s, method = "DL", test = "knha"),
                                               error = function(e2) NULL))
  if (is.null(fit)) return(NULL)
  data.frame(k = k, estimate = as.numeric(fit$beta), se = as.numeric(fit$se),
             ci_lo = fit$ci.lb, ci_hi = fit$ci.ub, tau2 = fit$tau2, I2 = fit$I2,
             p = fit$pval)
}

#' Contrast two pooled conditions
#'
#' \eqn{\Delta = \hat\mu_B - \hat\mu_A} with standard errors added in quadrature and a
#' Welch–Satterthwaite t test, each arm contributing `k - 1` degrees of freedom.
#'
#' @param pool_a,pool_b one-row data.frames returned by [pool_sections()] for the
#'   reference (A) and the contrasted (B) condition.
#' @return a one-row data.frame with `estimate`, `se`, `df`, `ci_lo`, `ci_hi`, `p`.
#' @export
contrast_conditions <- function(pool_a, pool_b) {
  if (is.null(pool_a) || is.null(pool_b)) return(NULL)
  est <- pool_b$estimate - pool_a$estimate
  se  <- sqrt(pool_a$se^2 + pool_b$se^2)
  df  <- (pool_a$se^2 + pool_b$se^2)^2 /
         (pool_a$se^4 / (pool_a$k - 1) + pool_b$se^4 / (pool_b$k - 1))
  tc  <- stats::qt(0.975, df)
  data.frame(estimate = est, se = se, df = df, ci_lo = est - tc * se, ci_hi = est + tc * se,
             p = 2 * stats::pt(-abs(est / se), df))
}

#' Test population-relative repositioning of one gene between two conditions
#'
#' Pools the eligible sections of each condition with [pool_sections()], contrasts the
#' two conditions with [contrast_conditions()] and classifies the change from the two
#' condition-level estimates and intervals.
#'
#' @param sections data.frame from [delta_by_section()] (one row per section) with a
#'   condition column.
#' @param condition name of the condition column.
#' @param levels character vector of length 2: the reference condition first, the
#'   contrasted condition second.
#' @inheritParams pool_sections
#' @return a list with `arms` (one pooled row per condition), `contrast` (the test of
#'   \eqn{\Delta = 0}) and `pattern` (see [classify_pattern()]).
#' @examples
#' f <- system.file("extdata", "igfbp5_fibro2_iaf4.csv", package = "TICR")
#' d <- read.csv(f)
#' s <- delta_by_section(d, "Igfbp5", keep = c("condition", "population"))
#' test_displacement(s, levels = c("Day0", "Day9"))
#' @export
test_displacement <- function(sections, condition = "condition", levels = NULL, min_k = 3) {
  stopifnot(condition %in% names(sections))
  if (is.null(levels)) levels <- unique(sections[[condition]])[1:2]
  el <- sections[sections$eligible %in% TRUE, , drop = FALSE]
  arms <- lapply(levels, function(l) {
    d <- el[el[[condition]] == l, , drop = FALSE]
    p <- pool_sections(d$delta, d$se0, min_k = min_k)
    if (is.null(p)) return(NULL)
    cbind(data.frame(condition = l, stringsAsFactors = FALSE), p)
  })
  names(arms) <- levels
  ct <- contrast_conditions(arms[[1]], arms[[2]])
  list(arms = do.call(rbind, arms), contrast = ct,
       pattern = if (is.null(ct)) NA_character_ else classify_pattern(arms[[1]], arms[[2]]))
}

#' Classify a repositioning event from the two condition-level estimates
#'
#' A condition is "positioned" when its pooled 95\% interval excludes zero. The event is
#' `gained` (positioned only in B), `lost` (only in A), `reversed` (opposite signs),
#' `amplified` or `dampened` (same sign, larger or smaller magnitude in B), or
#' `unpositioned` when neither interval excludes zero.
#'
#' @inheritParams contrast_conditions
#' @return a character scalar.
#' @export
classify_pattern <- function(pool_a, pool_b) {
  pos <- function(p) is.finite(p$ci_lo) && !(p$ci_lo < 0 && p$ci_hi > 0)
  a <- pos(pool_a); b <- pos(pool_b)
  if (!a && !b) return("unpositioned")
  if (!a &&  b) return("gained")
  if ( a && !b) return("lost")
  if (sign(pool_a$estimate) != sign(pool_b$estimate)) return("reversed")
  if (abs(pool_b$estimate) > abs(pool_a$estimate)) "amplified" else "dampened"
}
