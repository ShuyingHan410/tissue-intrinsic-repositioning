## v16 §0a A4(i): estimate the mouse-level variance from the real Day-9 arms by
## three-level REML (slices nested in mice), per gene x lineage x arm.
## Reports the mouse variance share and the implied sigma_mouse on the theta scale.
## Output: sigma_mouse.csv (per fit), sigma_mouse_summary.csv
suppressMessages({source("engine_v16.R"); library(parallel)})
CORES <- max(1L, min(6L, detectCores() - 1L))
GSC <- "../.."
slc <- read_csv(file.path(GSC, "gene_spatial_lineage_slices.csv"), show_col_types = FALSE, progress = FALSE)
gi  <- read_csv(file.path(GSC, "genes_index.csv"), show_col_types = FALSE)
dat <- slc %>% filter(eligible, is.finite(delta), is.finite(se_delta), se_delta > 0,
                      side %in% c("daughter_D9", "parent_D9")) %>%
  left_join(gi %>% select(gene, is_blank), by = "gene") %>%
  distinct(gene, lineage, side, Slice_ID, .keep_all = TRUE)

keys <- dat %>% group_by(gene, lineage, side, is_blank) %>%
  summarise(k = n(), n_mice = n_distinct(mouse), sd_r = median(sd_r), .groups = "drop") %>%
  filter(k >= 6, n_mice >= 3)
cat("fits to run:", nrow(keys), "\n")

fit_one <- function(i) {
  z <- keys[i, ]
  d <- dat %>% filter(gene == z$gene, lineage == z$lineage, side == z$side)
  fit <- try(rma.mv(yi = delta, V = se_delta^2, random = ~ 1 | mouse / Slice_ID,
                    data = d, method = "REML"), silent = TRUE)
  if (inherits(fit, "try-error")) return(data.frame(z, s2_mouse = NA, s2_slice = NA, mu = NA))
  data.frame(z, s2_mouse = fit$sigma2[1], s2_slice = fit$sigma2[2], mu = as.numeric(fit$beta))
}
res <- bind_rows(mclapply(seq_len(nrow(keys)), fit_one, mc.cores = CORES))
res <- res %>% mutate(share_mouse = s2_mouse / (s2_mouse + s2_slice),
                      sigma_mouse_delta = sqrt(s2_mouse),
                      sigma_mouse_theta = sqrt(s2_mouse) / sd_r,       # tau_delta ~ tau_theta * SD(r)
                      tau_slice_theta   = sqrt(s2_slice) / sd_r)
write_csv(res, "sigma_mouse.csv")

summ <- res %>% filter(!is_blank, is.finite(s2_mouse)) %>% group_by(lineage, side) %>%
  summarise(n_genes = n(), k_med = median(k), mice_med = median(n_mice),
            frac_mouse_var_pos = mean(s2_mouse > 1e-8),
            share_mouse_med = median(share_mouse, na.rm = TRUE),
            share_mouse_p75 = quantile(share_mouse, .75, na.rm = TRUE),
            sigma_mouse_theta_med = median(sigma_mouse_theta), sigma_mouse_theta_p75 = quantile(sigma_mouse_theta, .75),
            sigma_mouse_theta_p90 = quantile(sigma_mouse_theta, .90),
            tau_slice_theta_med = median(tau_slice_theta), .groups = "drop")
write_csv(summ, "sigma_mouse_summary.csv")
options(width = 200); print(as.data.frame(summ %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
cat("\nblanks (should show ~no mouse variance):\n")
print(as.data.frame(res %>% filter(is_blank, is.finite(s2_mouse)) %>% group_by(side) %>%
  summarise(n = n(), frac_pos = mean(s2_mouse > 1e-8), sigma_mouse_theta_med = median(sigma_mouse_theta),
            sigma_mouse_theta_p90 = quantile(sigma_mouse_theta, .9)) %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
cat("DONE 03_sigma_mouse\n")
