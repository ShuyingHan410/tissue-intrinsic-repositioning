# v3 STRATIFIED sweep: each (lineage x LAYER) is an independent analysis.
#   Rscript 01c_stratified_chunk.R <chunk_id>  -> resultsC/part_<id>.csv
#
# Differs from v2 (Method B) in that layers are NEVER pooled. delta is computed
# per (gene x type x layer x slice) and pooled over SLICES WITHIN THAT LAYER ONLY.
# No coordinate is ever averaged across layers, and no data-dependent decision is
# made about whether to pool -- stratification is pre-specified.
#
# FILTERS, in order, per (layer x slice) observation:
#   1. delta and SE must be finite and SE > 0        (degenerate cases)
#   2. >= 30 cells of the cell type                  (MIN_CELLS)
#   3. PRIMARY: SE < 0.05                            (SE_MAX)  -- precision rule
#      ALTERNATIVE: >= 20 molecules                  (MIN_MOL) -- set SENS_MOL=1
# Then a stratum needs >= 3 surviving slices (MIN_SLICES) to be pooled at all.
# Each stage's counts are written to audit_NN.csv so the filtering is auditable.
suppressMessages({library(dplyr); library(readr); library(tidyr); library(purrr); library(metafor)})
args <- commandArgs(trailingOnly=TRUE); cid <- as.integer(args[1]); stopifnot(is.finite(cid))
OUT<-Sys.getenv("OUT","."); COORD<-Sys.getenv("COORD_DIR","../output_threelayer/per_slice_morphogam")   # config.sh exports both
LAYERS<-c("MU","SM","ME")
MIN_CELLS <- 30L          # >= 30 cells   (unchanged from v1/v2)
MIN_MOL   <- 20L          # >= 20 molecules per (layer x slice).
# Chosen on power, not convention: observations below 20 molecules are 30% of the
# data but carry only 2.5% of the inverse-variance weight, so excluding them barely
# moves the estimate while removing genes that could never reach significance from
# the BH denominator. Their median per-slice SE (0.035-0.085) is comparable to or
# larger than the real effects (0.01-0.09).
SE_MAX    <- 0.05         # sensitivity alternative: threshold on PRECISION instead of
# counts. Adapts to cell number and var(r), which a flat count cannot. Applied in the
# SENS_SE run; the primary run uses MIN_MOL.
SENS_MOL  <- identical(Sys.getenv("SENS_MOL"), "1")   # set to run the count rule instead
MIN_SLICES<- 3L
Sys.setenv(OUT=OUT); source("inclusion_criteria.R")     # Table S2 routing, as v1/v2
LIN<-list(list(name="Fibro6->IAF2",parent="Fibro 6",daughter="IAF 2"),
          list(name="Fibro6->IAF3",parent="Fibro 6",daughter="IAF 3"),
          list(name="Fibro2->IAF4",parent="Fibro 2",daughter="IAF 4"),
          list(name="Fibro7->Fibro4",parent="Fibro 7",daughter="Fibro 4"))
`%||%`<-function(a,b) if (is.null(a)) b else a

meta<-read_csv(file.path(OUT,"meta.csv.gz"),show_col_types=FALSE,progress=FALSE)%>%
  mutate(Cell_ID=as.character(Cell_ID))
cf<-file.path(OUT,"chunks",sprintf("chunk_%02d.csv.gz",cid)); stopifnot(file.exists(cf))
counts<-read_csv(cf,show_col_types=FALSE,progress=FALSE)%>%mutate(Cell_ID=as.character(Cell_ID))
genes<-setdiff(names(counts),c("Slice_ID","Cell_ID"))
co<-map_dfr(LAYERS,function(ly){fs<-list.files(COORD,sprintf("_%s_xyt_mapped[.]csv$",ly),full.names=TRUE)
 map_dfr(fs,function(f)suppressMessages(read_csv(f,show_col_types=FALSE,progress=FALSE))%>%
  transmute(Slice_ID=as.character(Slice_ID),Cell_ID=as.character(Cell_ID),
            r=as.numeric(r_for_jsd),layer=ly))})%>%filter(is.finite(r),r>=0,r<=1)
J<-counts%>%inner_join(meta%>%select(Slice_ID,Cell_ID,Tier2,day,mouse,total_counts),
                       by=c("Slice_ID","Cell_ID"))%>%inner_join(co,by=c("Slice_ID","Cell_ID"))
cat(sprintf("chunk %d: %d genes, %s cells\n",cid,length(genes),format(nrow(J),big.mark=",")))

# routing is a DETECTION question -> computed on all layers pooled, exactly as v1/v2
ROUTE<-map_dfr(genes,function(g) map_dfr(LIN,function(Ln){
  rt<-tryCatch(route_gene(J,g,list(name=Ln$name,parent=Ln$parent,daughter=Ln$daughter)),
               error=function(e) NULL)
  if(is.null(rt)) return(tibble())
  tibble(gene=g,lineage=Ln$name,route_s2=rt$route)}))

pool<-function(d,s){ if(length(d)<MIN_SLICES) return(NULL)
  m<-tryCatch(rma(yi=d,vi=s^2,method="REML",test="knha"),error=function(e)NULL)
  if(is.null(m)) return(NULL)
  list(est=as.numeric(m$beta[1]),se=as.numeric(m$se),lo=m$ci.lb[1],hi=m$ci.ub[1],
       p=m$pval[1],tau2=m$tau2,I2=m$I2,k=length(d))}

res<-list(); AUDIT<-list()
for(g in genes){ gv<-J[[g]]
 for(Ln in LIN){
  rt<-ROUTE$route_s2[ROUTE$gene==g & ROUTE$lineage==Ln$name]; rt<-if(length(rt)) rt[1] else NA
  for(ly in LAYERS){
   for(side in c("parent_D0","daughter_D9")){
    ct<-if(side=="parent_D0") Ln$parent else Ln$daughter
    dy<-if(side=="parent_D0") "Day 0" else "Day 9"
    k<-J$Tier2==ct & J$day==dy & J$layer==ly
    if(!any(k)) next
    z<-tibble(Slice_ID=J$Slice_ID[k],r=J$r[k],n=gv[k],tot=J$total_counts[k])
    sl<-z%>%group_by(Slice_ID)%>%
      summarise(cells=n(),mol=sum(n),tot=sum(tot),
                delta=ifelse(sum(n)>0,sum(n*r)/sum(n)-mean(r),NA_real_),
                se=sqrt(var(r)*(sum(n^2)/sum(n)^2-1/n())),.groups="drop")%>%
      filter(is.finite(delta), is.finite(se), se>0)
    n_all <- nrow(sl); n_cells_ok <- sum(sl$cells >= MIN_CELLS)
    sl <- sl %>% filter(cells >= MIN_CELLS) %>%
      { if (SENS_MOL) filter(., mol >= MIN_MOL) else filter(., se < SE_MAX) }
    AUDIT[[length(AUDIT)+1]] <- tibble(gene=g, lineage=Ln$name, layer=ly, side=side,
      slices_total=n_all, pass_cells=n_cells_ok, pass_filter=nrow(sl))
    if(nrow(sl)<MIN_SLICES) next
    pl<-pool(sl$delta,sl$se); if(is.null(pl)) next
    res[[length(res)+1]]<-tibble(gene=g,lineage=Ln$name,layer=ly,side=side,cell_type=ct,
      route_s2=rt, k_slices=pl$k, cells=sum(sl$cells), mol=sum(sl$mol),
      rate_tx=1e4*sum(sl$mol)/sum(sl$tot),
      delta=pl$est,se=pl$se,ci_lo=pl$lo,ci_hi=pl$hi,p=pl$p,tau2=pl$tau2,I2=pl$I2)
   }}}}
OUTDIR <- if (SENS_MOL) "resultsC_mol" else "resultsC"
dir.create(OUTDIR,showWarnings=FALSE)
out<-bind_rows(res); write_csv(out,file.path(OUTDIR,sprintf("part_%02d.csv",cid)))
write_csv(bind_rows(AUDIT), file.path(OUTDIR,sprintf("audit_%02d.csv",cid)))
cat(sprintf("chunk %d done: %d (gene x lineage x layer x side) rows\n",cid,nrow(out)))
