## R packages for the tissue-intrinsic repositioning framework (R >= 4.5).
## Run from the repository root:  Rscript environment/install_R_packages.R
cran <- c("dplyr", "readr", "tidyr", "purrr", "stringr", "tibble", "data.table",
          "ggplot2", "patchwork", "knitr", "rmarkdown", "pheatmap", "scales",
          "FNN", "Matrix", "igraph", "RSpectra", "mgcv", "metafor", "ripserr",
          "SeuratObject", "remotes", "BiocManager")
new <- cran[!cran %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")

## MorphoGAM (Nicol et al. 2024), at the commit used for the paper (v1.2.0)
if (!requireNamespace("MorphoGAM", quietly = TRUE))
  remotes::install_github("phillipnicol/MorphoGAM", ref = "fc01ea98d7a95c97842a73e802804f94846170cd")

## Simulation comparators only (DESpace2 >= 2.2 and edgeR; Bioconductor)
bioc <- c("DESpace", "SingleCellExperiment", "edgeR")
newb <- bioc[!bioc %in% rownames(installed.packages())]
if (length(newb)) BiocManager::install(newb, update = FALSE, ask = FALSE)

for (p in c(cran, "MorphoGAM", bioc))
  cat(sprintf("%-22s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")))
