Example data (mouse colon MERFISH, Cadinu et al. 2024, Cell; coordinates from the paper's pipeline).

* `colon_fibroblast_two_sections.csv`: fibroblast-lineage cells of one healthy (Day 0, `062921_D0_m3a_1_slice_1`) and one DSS-colitis (Day 9, `100221_D9_m5_2_slice_1`) section: `section, condition, layer (MU/SM/ME), population, cell_id, r`.
* `igfbp5_fibro2_iaf4.csv`: every Fibro 2 cell of the ten Day 0 sections and every IAF 4 cell of the Day 9 sections (mucosa): `section, condition, mouse, population, layer, cell_id, r, Igfbp5` (transcript count, zeros included).

`r` is the normalized within-layer position, 0 = lumen, 1 = outer boundary.
