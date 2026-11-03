Colon coordinate stage. `per_slice_morphogam/` (the coordinate files the paper uses) is not tracked in git:
extract it with `tar -xzf data/colon_tissue_intrinsic_coordinates.tar.gz` from the repository root.
`flip_orientation.R` harmonises orientation (lumen = 0 -> outer = 1); `orientation_harmonisation_audit.csv`
and `need_flip.csv` are its audit; `per_slice_morphogam_fit_summary.csv` and `per_slice_morphogam_failure_log.csv`
record which fits needed a larger kNN or a hand-traced curve (`three_layer.Rmd`).
