# v3 merge -- STRATIFIED inference, POOLED calibration.
#   Rscript 02d_stratified_merge.R [se]      ("se" = the SE<0.05 sensitivity run)
#
# Inference is stratified: each (lineage x layer) is a self-contained analysis and
# BH runs within (lineage x layer x route). Layers are never pooled.
#
# CALIBRATION is NOT stratified by layer. Per-layer lambda was tested and rejected:
# blank z-scores per layer gave lambda 1.27 (MU, n=1461), 0.93 (SM, n=197), 1.42
# (ME, n=115), and a Kruskal-Wallis on z^2 across layers gave p = 0.11 -- the spread
# is sampling noise, not a real layer effect. Estimating nine lambdas would fit that
# noise, and ME/SM would often fall below the 20-null minimum and silently revert to
# lambda = 1 (no correction at all). The dominant inflation mechanism -- 31 slices
# from 12 mice -- is layer-independent. So lambda is estimated per STATISTIC, pooled
# across layers: three values, each from a 3x larger null sample.
suppressMessages({library(dplyr);library(readr);library(tidyr);library(purrr)})
SENS <- identical(commandArgs(trailingOnly=TRUE)[1],"mol")
IN   <- if (SENS) "resultsC_mol" else "resultsC"
SUF  <- if (SENS) "_mol" else ""
arms <- bind_rows(lapply(list.files(IN,"^part_.*csv$",full.names=TRUE),
                         function(f) suppressMessages(read_csv(f,show_col_types=FALSE))))
write_csv(arms, sprintf("armsC%s.csv",SUF))

wide <- arms %>%
  select(gene,lineage,layer,side,route_s2,delta,se,ci_lo,ci_hi,p,k_slices,cells,mol,rate_tx,I2,tau2) %>%
  group_by(gene,lineage,layer) %>% mutate(route_s2=first(na.omit(route_s2))) %>% ungroup() %>%
  pivot_wider(names_from=side,
              values_from=c(delta,se,ci_lo,ci_hi,p,k_slices,cells,mol,rate_tx,I2,tau2))
for (nm in c("delta_parent_D0","delta_daughter_D9","se_parent_D0","se_daughter_D9",
             "ci_lo_parent_D0","ci_hi_parent_D0","ci_lo_daughter_D9","ci_hi_daughter_D9",
             "p_parent_D0","p_daughter_D9","k_slices_parent_D0","k_slices_daughter_D9",
             "rate_tx_parent_D0","rate_tx_daughter_D9"))
  if (!nm %in% names(wide)) wide[[nm]] <- NA_real_

con <- wide %>% mutate(
  route = case_when(
    route_s2=="D1_full"       & !is.na(delta_parent_D0) & !is.na(delta_daughter_D9) ~ "D1_full",
    route_s2=="daughter_only" & !is.na(delta_daughter_D9) ~ "daughter_only",
    route_s2=="parent_only"   & !is.na(delta_parent_D0)   ~ "parent_only",
    TRUE ~ "excluded"),
  d1     = delta_daughter_D9 - delta_parent_D0,
  d1_se  = sqrt(se_daughter_D9^2 + se_parent_D0^2),
  d1_df  = d1_se^4/(se_daughter_D9^4/(k_slices_daughter_D9-1)+se_parent_D0^4/(k_slices_parent_D0-1)),
  d1_p_raw = 2*pt(-abs(d1/d1_se), d1_df),
  d1_ci_lo = d1 - qt(.975,d1_df)*d1_se,
  d1_ci_hi = d1 + qt(.975,d1_df)*d1_se,
  log2fc_expr = log2(rate_tx_daughter_D9/rate_tx_parent_D0),
  eff   = ifelse(route=="D1_full", d1, ifelse(route=="daughter_only", delta_daughter_D9, delta_parent_D0)),
  p_raw = ifelse(route=="D1_full", d1_p_raw, ifelse(route=="daughter_only", p_daughter_D9, p_parent_D0)))

# ---- lambda per STATISTIC, pooled across layers -----------------------------------
lam_of <- function(p){ p<-p[is.finite(p)&p>0&p<1]
  if(length(p)<20) return(NA_real_); median(qnorm(p/2)^2)/qchisq(.5,1) }
bl <- con %>% filter(grepl("^Blank",gene))
LAM <- tibble(
  statistic = c("D1","daughter_arm","parent_arm"),
  n_blank   = c(sum(is.finite(bl$d1_p_raw)), sum(is.finite(bl$p_daughter_D9)), sum(is.finite(bl$p_parent_D0))),
  fp_before = c(mean(bl$d1_p_raw<.05,na.rm=TRUE), mean(bl$p_daughter_D9<.05,na.rm=TRUE),
                mean(bl$p_parent_D0<.05,na.rm=TRUE)),
  lambda    = c(lam_of(bl$d1_p_raw), lam_of(bl$p_daughter_D9), lam_of(bl$p_parent_D0)))
LAM$lambda[!is.finite(LAM$lambda) | LAM$lambda < 1] <- 1
apply_lam <- function(p,l){ z<-abs(qnorm(p/2)); 2*pnorm(-z/sqrt(l)) }
lam_for <- function(rt) case_when(rt=="D1_full" ~ LAM$lambda[1],
                                  rt=="daughter_only" ~ LAM$lambda[2],
                                  rt=="parent_only" ~ LAM$lambda[3], TRUE ~ 1)
con <- con %>% mutate(lambda = lam_for(route), p_cal = apply_lam(p_raw, lambda))

# ---- calibrate the ENDPOINT tests too -------------------------------------------
# The endpoint CIs and p-values come straight from rma() and are UNcalibrated, yet
# they decide the `pattern` label via the positioned / not-positioned call. Because
# the blanks show the arm-level nulls are inflated (lambda_parent, lambda_daughter),
# raw intervals are ~11-19% too narrow, i.e. too ready to call a gene positioned.
# Widen each endpoint interval about its own estimate by sqrt(lambda_arm), and apply
# the same factor to the endpoint p-values, so pattern classification and the
# displayed endpoint significance use the SAME calibrated bar as the D1 result.
widen <- function(est, lo, hi, lam) list(lo = est - (est - lo)*sqrt(lam),
                                         hi = est + (hi - est)*sqrt(lam))
lamP <- LAM$lambda[LAM$statistic == "parent_arm"]
lamD <- LAM$lambda[LAM$statistic == "daughter_arm"]
.p0 <- widen(con$delta_parent_D0,   con$ci_lo_parent_D0,   con$ci_hi_parent_D0,   lamP)
.p9 <- widen(con$delta_daughter_D9, con$ci_lo_daughter_D9, con$ci_hi_daughter_D9, lamD)
con <- con %>% mutate(
  ci_lo_parent_D0_raw = ci_lo_parent_D0, ci_hi_parent_D0_raw = ci_hi_parent_D0,
  ci_lo_daughter_D9_raw = ci_lo_daughter_D9, ci_hi_daughter_D9_raw = ci_hi_daughter_D9,
  ci_lo_parent_D0 = .p0$lo, ci_hi_parent_D0 = .p0$hi,
  ci_lo_daughter_D9 = .p9$lo, ci_hi_daughter_D9 = .p9$hi,
  p_parent_D0_raw = p_parent_D0, p_daughter_D9_raw = p_daughter_D9,
  p_parent_D0   = apply_lam(p_parent_D0,   lamP),
  p_daughter_D9 = apply_lam(p_daughter_D9, lamD))
LAM$fp_after <- c(mean(apply_lam(bl$d1_p_raw,LAM$lambda[1])<.05,na.rm=TRUE),
                  mean(apply_lam(bl$p_daughter_D9,LAM$lambda[2])<.05,na.rm=TRUE),
                  mean(apply_lam(bl$p_parent_D0,LAM$lambda[3])<.05,na.rm=TRUE))
write_csv(LAM, sprintf("lambdaC%s.csv",SUF)); print(as.data.frame(LAM))
# blank nulls per layer, recorded so the pooling choice can be audited
blk_layer <- bl %>% group_by(layer) %>% summarise(
  D1=sum(is.finite(d1_p_raw)), daughter=sum(is.finite(p_daughter_D9)),
  parent=sum(is.finite(p_parent_D0)), .groups="drop")
write_csv(blk_layer, sprintf("blankC_by_layer%s.csv",SUF)); print(as.data.frame(blk_layer))

# ---- BH within lineage x LAYER x route (inference IS stratified) -------------------
real <- con %>% filter(!grepl("^Blank",gene), route %in% c("D1_full","daughter_only","parent_only")) %>%
  group_by(lineage, layer, route) %>% mutate(q = p.adjust(p_cal,"BH"), m_family = n()) %>% ungroup()
pos <- function(lo,hi) is.finite(lo)&is.finite(hi)&(lo>0|hi<0)
real <- real %>% mutate(
  d0_positioned = pos(ci_lo_parent_D0, ci_hi_parent_D0),
  d9_positioned = pos(ci_lo_daughter_D9, ci_hi_daughter_D9),
  pattern = case_when(route!="D1_full" ~ NA_character_,
    !d0_positioned &  d9_positioned ~ "gained",
     d0_positioned & !d9_positioned ~ "lost",
     d0_positioned &  d9_positioned & sign(delta_parent_D0)!=sign(delta_daughter_D9) ~ "reversed",
     d0_positioned &  d9_positioned & abs(delta_daughter_D9)>abs(delta_parent_D0) ~ "amplified",
     d0_positioned &  d9_positioned ~ "dampened", TRUE ~ "weak-both-ends"),
  direction = ifelse(eff<0,"luminal","crypt-base"),
  stratum   = paste(lineage, layer, sep=" / "))
write_csv(real, sprintf("contrastsC%s.csv",SUF))
cat(sprintf("\nv3%s: %d rows | D1_full %d | significant %d | strata %d\n", SUF,
  nrow(real), sum(real$route=="D1_full"),
  sum(real$route=="D1_full" & real$q<0.05, na.rm=TRUE), n_distinct(real$stratum)))
