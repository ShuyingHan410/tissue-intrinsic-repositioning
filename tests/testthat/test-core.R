test_that("section_delta matches the closed form and gates sparse sections", {
  set.seed(1)
  r <- runif(400); n <- rpois(400, 0.5)
  s <- section_delta(r, n)
  expect_equal(s$delta, sum(n * r) / sum(n) - mean(r))
  expect_equal(s$se0, sd(r) * sqrt(sum(n^2) / sum(n)^2 - 1 / length(r)))
  expect_equal(s$n_eff, sum(n)^2 / sum(n^2))
  # uniform counts: transcript-weighted mean equals the cell mean
  expect_equal(section_delta(r, rep(2, 400))$delta, 0)
  # too few cells -> ineligible
  expect_false(section_delta(r[1:20], n[1:20])$eligible)
})

test_that("pooling and contrast reproduce metafor and Welch degrees of freedom", {
  y <- c(0.02, 0.03, 0.01, 0.04); s <- c(0.010, 0.012, 0.009, 0.011)
  p <- pool_sections(y, s)
  fit <- metafor::rma(yi = y, sei = s, method = "REML", test = "knha")
  expect_equal(p$estimate, as.numeric(fit$beta))
  expect_equal(p$se, as.numeric(fit$se))
  expect_null(pool_sections(y[1:2], s[1:2]))
  q <- pool_sections(-y, s)
  ct <- contrast_conditions(p, q)
  expect_equal(ct$estimate, q$estimate - p$estimate)
  expect_equal(ct$df, (p$se^2 + q$se^2)^2 / (p$se^4 / 3 + q$se^4 / 3))
  expect_equal(classify_pattern(p, q), "reversed")
})

test_that("test_displacement recovers the Igfbp5 example", {
  f <- system.file("extdata", "igfbp5_fibro2_iaf4.csv", package = "TICR")
  d <- read.csv(f)
  sec <- delta_by_section(d, "Igfbp5", keep = c("condition", "population"))
  res <- test_displacement(sec, levels = c("Day0", "Day9"))
  expect_equal(nrow(res$arms), 2)
  expect_gt(res$contrast$estimate, 0.05)
  expect_lt(res$contrast$p, 1e-4)
  expect_equal(res$pattern, "reversed")
})

test_that("population_similarity is symmetric under identical groups and ranks the atlas match", {
  f <- system.file("extdata", "colon_fibroblast_two_sections.csv", package = "TICR")
  d <- read.csv(f)
  S <- population_similarity(d, groups = c("Day0", "Day9"))
  expect_true(is.matrix(S)); expect_true(all(S[is.finite(S)] >= 0))
  rk <- rank_candidates(S, "Fibro 2", target = "IAF 4")
  expect_lte(attr(rk, "rank"), 3)
  # a group compared with itself scores its own populations highest on the diagonal
  d0 <- d[d$condition == "Day0", ]; d0b <- d0; d0b$condition <- "copy"
  S2 <- population_similarity(rbind(d0, d0b), groups = c("Day0", "copy"))
  for (p in rownames(S2)) expect_equal(names(which.max(S2[p, ])), p)
})

test_that("js_distance is zero for identical samples and increases with separation", {
  set.seed(2); x <- runif(300)
  expect_equal(js_distance(x, x, bw = 0.05), 0)
  expect_gt(js_distance(x, pmin(1, x + 0.3), bw = 0.05), js_distance(x, pmin(1, x + 0.05), bw = 0.05))
})
