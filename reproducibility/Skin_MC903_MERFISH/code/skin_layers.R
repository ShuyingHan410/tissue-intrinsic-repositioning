# =====================================================================
# skin_layers.R - the ONE place the skin layer labels are defined. Sourced by overview.Rmd and export_slices.R.
#
#   assign_skin_layers(d)  ->  d with depth, anchor, b_sub, b_fas, layer_depth, layer_nn, name
#
# A. DEPTH-DEFINED COMPARTMENTS (the layer the method uses). Label-blind except for anatomical landmarks:
#    the interfollicular keratinocyte sheet marks the skin surface (as the lumen marks the colon), adipocytes
#    mark the fat onset, muscle marks the panniculus carnosus. No fibroblast, immune or neighbourhood label enters.
#    depth = distance (um) from a cell to the nearest keratinocyte of the CONNECTED epidermal sheet in its strip.
#      (keratinocytes < 40 um apart are linked; only components of >= 300 cells count; this drops hair-follicle
#       bulbs deep in lesional skin that carry "stressed KC" labels and would act as false surface points)
#    Boundaries, superficial -> deep:
#      EPI      depth <= 10 um                 the keratinocyte sheet and cells embedded in it
#      PAP      10 um < depth <= 0.4 * b_sub   papillary dermis = upper 40% of the dermis (epidermis -> fat onset).
#                                               Mouse dermis is only ~80-280 um thick, so a human 100 um cut is meaningless;
#                                               the RELATIVE cut is scale-free, absorbs oedema, and is not fitted to the labels
#      RET      0.4 * b_sub < depth <= b_sub   reticular dermis, down to the first fat
#      SUB      b_sub < depth <= b_fas         subcutis: b_sub = 10th percentile of adipocyte depth (dermal white adipose
#                                               tissue), clamped to 80-350 um; 150 um if the strip has < 30 adipocytes
#      SUBFAS   depth > b_fas                  subfascia: below the panniculus; b_fas = median depth of the strip's muscle
#                                               (SM) cells when >= 100 of them lie deeper than b_sub + 50 um; else not captured
#    Ear cartilage splits the ear into two faces: ear depths beyond ~100 um are unreliable (chondrocytes shown as reference).
#
# B. AUTHORS' NEIGHBOURHOOD CLUSTERS (CellCharter k = 11) named and collapsed to layers, mimicking MU/SM/ME in the colon.
#    Named from: cell-type composition, median depth under A, prevalence in normal vs disease, and the paper's text
#    (NN1 epidermis, NN2 upper pilosebaceous unit, NN3 dermis, NN4 subcutis, NN5 ear cartilage; IAN1-6 inflammation).
#    Collapsed by depth to EPI / DER / SUB / CART; cluster 3 (neutrophil) spans crust and dermis and is kept apart.
#    Used ONLY as the authors'-label comparator: it puts most papillary fibroblasts into "EPI".
# =====================================================================
suppressMessages({ library(data.table); library(FNN) })

PAP_FRAC <- 0.4
LAY5 <- c("EPI", "PAP", "RET", "SUB", "SUBFAS")
PAL5 <- c(EPI = "#d62728", PAP = "#ff7f0e", RET = "#1f77b4", SUB = "#f7e300", SUBFAS = "#984ea3", `CART (ref)` = "#7f7f7f")

SKIN_CC <- data.table(
  cc11 = c(1, 9, 0, 3, 5, 4, 2, 10, 7, 6, 8),
  # cluster -> paper name. The paper (Results, Fig 3C-D, S3D-F) places: IAN1 epidermis; IAN2 dermo-epidermal junction ("slightly deeper");
  # IAN3 and IAN4 dermis (IAN2/IAN3 rise early, IAN4/IAN5 later and more in OXA); IAN5 epicutaneous serous crust; IAN6 lower anagen hair
  # follicles of OXA tissue (not truly inflammation-associated). Matched to clusters by composition, depth and time course (see overview sec. 4):
  #   9 = IAN1 (stressed KC, 0 um)   0 = IAN2 (KC/SG/DC at 11 um)   2 = IAN3 (FIB/MAC/DC/mast, early, MC903 d6 peak)
  #   10 = IAN4 (FIB/MAC, late, OXA d16)   3 = IAN5 (neutrophil crust)   7 = IAN6 (anagen HF, OXA)
  name = c("NN1 epidermis", "IAN1 inflamed epidermis (stressed KC)", "IAN2 dermo-epidermal junction (KC/SG/DC)", "IAN5 epicutaneous crust (neutrophil)", "NN2 upper pilosebaceous unit",
           "NN3 dermis", "IAN3 inflamed dermis, early (FIB/MAC/DC/mast)", "IAN4 inflamed dermis, late (FIB/MAC)", "IAN6 lower anagen hair follicle (OXA)", "NN4 subcutis", "NN5 cartilage (ear)"),
  layer_nn = c("EPI", "EPI", "EPI", "EPI/infiltrate", "DER", "DER", "DER", "DER", "DER", "SUB", "CART"),
  reading = c("normal interfollicular epidermis; vanishes by MC903 d6",
              "hyperplastic epidermis of stressed basal / cycling / differentiated KC; 1% of normal cells, 30% at MC903 d3 (NN1 -> IAN1 switch)",
              "dermo-epidermal interface in disease: stressed KC with sebaceous and dendritic cells",
              "neutrophil state: surface crusts plus dermal infiltrate; oxazolone only",
              "pilosebaceous niche (sebaceous gland + outer root sheath); present in all conditions",
              "normal dermis: fibroblasts with Schwann cells and vessels; replaced by clusters 2 and 10 in disease",
              "inflamed dermis rich in macrophages, DC and mast cells (MC903 d6 peak; NN3 -> IAN2/3)",
              "inflamed dermis of fibroblasts + macrophages, emerging later and more in OXA; holds 50% of pro-inflammatory fibroblasts; same depth as IAN3 (not deeper)",
              "anagen hair follicles, almost specific to oxazolone (inflammation-induced hair cycle)",
              "subcutis: fibroblasts, muscle, adipocytes; stable anatomy", "ear cartilage with perichondrial fibroblasts; stable anatomy"))
PAL_NN  <- setNames(c("#d62728","#e377c2","#ffbb78","#2ca02c","#ff7f0e","#1f77b4","#17becf","#9467bd","#bcbd22","#f7e300","#7f7f7f"), SKIN_CC$name)
PAL_NN3 <- c(EPI = "#d62728", `EPI/infiltrate` = "#2ca02c", DER = "#1f77b4", SUB = "#f7e300", CART = "#7f7f7f")

.sheet <- function(xy) {
  nn <- get.knn(xy, k = 10); i <- rep(seq_len(nrow(xy)), 10); j <- as.vector(nn$nn.index); keep <- as.vector(nn$nn.dist) < 40
  g <- igraph::graph_from_edgelist(cbind(i[keep], j[keep]), directed = FALSE); g <- igraph::add_vertices(g, max(0, nrow(xy) - igraph::vcount(g)))
  cm <- igraph::components(g)$membership; cm %in% which(tabulate(cm) >= 300)
}

assign_skin_layers <- function(d) {
  d <- copy(d)
  d[, kc := celltype %in% c("KC_IFE", "KC_Stressed_Basal", "KC_Stressed_Cycling", "KC_Stressed_Diff")]
  d[, anchor := FALSE]; d[kc == TRUE, anchor := .sheet(cbind(x, y)), by = sample]
  d[, depth := { k <- which(anchor); if (length(k) < 50) NA_real_ else get.knnx(cbind(x[k], y[k]), cbind(x, y), k = 1)$nn.dist[, 1] }, by = sample]
  B <- d[, {
    ad <- depth[celltype == "Adipocyte"]; b_sub <- if (length(ad) >= 30) min(350, max(80, quantile(ad, .10))) else 150
    sm <- depth[celltype == "SM" & depth > b_sub + 50]; b_fas <- if (length(sm) >= 100) median(sm) else Inf
    .(n_adipo = length(ad), b_sub = round(b_sub), n_sm_deep = length(sm), b_fas = round(b_fas))
  }, by = sample]
  d <- merge(d, B[, .(sample, b_sub, b_fas)], by = "sample")
  d[, layer_depth := factor(fcase(depth <= 10, "EPI", depth <= PAP_FRAC * b_sub, "PAP", depth <= b_sub, "RET", depth <= b_fas, "SUB", default = "SUBFAS"), levels = LAY5)]
  d <- merge(d, SKIN_CC[, .(cc11, name, layer_nn)], by = "cc11", all.x = TRUE)
  d[, name := factor(name, levels = SKIN_CC$name)]; d[, layer_nn := factor(layer_nn, levels = names(PAL_NN3))]
  attr(d, "bounds") <- B
  d
}
