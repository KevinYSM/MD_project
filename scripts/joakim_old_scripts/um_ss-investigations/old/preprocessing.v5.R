library(Seurat)
library(future)
library(tidyverse)
library(djvdj)
library(plyranges)
library(rtracklayer)
library(DoubletFinder)
library(pheatmap)
library(SeuratWrappers)
library(rliger)

extract_raw <- function(dat.cca, project.name = NULL){
  
  if (!is.null(project.name)){
    snames <- unique(dat.cca@meta.data$Sample.ID[dat.cca@meta.data$project.name %in% project.name])
  } else {
    snames <- unique(dat.cca@meta.data$Sample.ID)
  }
  
  dat.processed <- list()
  for (sname in snames){
    dat <- subset(dat.cca, Sample.ID == sname)
    meta <- dat.cca@meta.data[dat.cca$Sample.ID == sname,]
    dat <- GetAssayData(dat, slot = "counts", assay = "RNA")
    dat <- CreateSeuratObject(counts = dat, meta.data = meta)
    dat.processed[[sname]] <- dat
  } 
  
  return(dat.processed)
}

outdir <- "~/proj/um_ss/Investigations/seurat/results/liger_all/"
dir.create(outdir,showWarnings = F,recursive = T)

dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")

dat.tumor <- extract_raw(dat.cca=dat.cca, project.name = "GWA-JN-388")

# LIGER ------------------------------------------------------------------------
lig.tumor <- seuratToLiger(objects = dat.tumor,
                           combined.seurat = F,
                           names = "use-meta",
                           meta.var = "Sample.ID",
                           renormalize=F,
                           use.idents=T,
                           use.seurat.genes=F,
                           remove.missing=F)

lig.tumor <- normalize(lig.tumor)
lig.tumor <- selectGenes(lig.tumor, combine = "union", do.plot = T)
lig.tumor <- scaleNotCenter(lig.tumor)
lig.tumor <- optimizeALS(lig.tumor, k = 50, lambda = 5, verbose = T)

saveRDS(lig.tumor,paste0(outdir,"lig.tumor.rda"))

# Tumor only -------------------------------------------------------------------
outdir <- "~/proj/um_ss/Investigations/seurat/results/liger_all/"

lig.tumor <- readRDS(paste0(outdir,"lig.tumor.rda"))
gsea_output <- rliger::runGSEA(lig.tumor)

dims.remove <- c(6,23,29) # Cell cycle
dims.remove <- c(dims.remove,28,31,32,35,36,42,43) # Translation and DNA repair

dims.use <- setdiff(1:50,dims.remove)
lig.tumor <- quantile_norm(lig.tumor,ref_dataset = "C8-GEX",dims.use=dims.use,
                           quantiles = 100, eps = 0.01)
lig.tumor <- runUMAP(lig.tumor, distance = "cosine", dims.use=dims.use, 
                     n_neighbors = 20, min_dist = 0.05)

lig.tumor.louvain <- louvainCluster(lig.tumor, resolution = 1)
all.plots.louvain  <- plotByDatasetAndCluster(lig.tumor.louvain, 
                                              axis.labels = c('UMAP 1', 'UMAP 2'), 
                                              return.plots = T)

all.plots.reannotated <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                                 return.plots = T, clusters = lig.reannotated)

#all.plots.clusters <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
#                                     return.plots = T)
vasu_phenotype <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/vasu_phenotype.rda")
all.plots.vasu <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T, clusters = vasu_phenotype,do.shuffle = F,
                                     reorder.idents = T,new.order = c("None","Unidentified","MART1","CD3+41bb+"))

all.plots.reannotated[[2]] + all.plots.vasu[[2]] + all.plots.louvain[[2]]

# Investigate components
gene_loadings <- plotGeneLoadings(lig.tumor.louvain, do.spec.plot = T, return.plots = TRUE,
                                  factor.share.thresh = Inf)
setdiff(which(unlist(lapply(gene_loadings,length))>0),dims.remove)

# 16: Senescence
# 28,31,32,36,42: Translation, NMD
# 35: DNA repair
# 43: apoptosis, DNA damage
# 40: Tregs
# 12: The outlier melanoma
i <- 12
gene_loadings[[i]]
g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10)


all.plots.reannotated[[2]]
all.plots.reannotated[[2]] + gene_loadings[[i]]

# Plot specific genes
rliger::plotGene(lig.tumor,"TRBV30",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"TRGV10",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"TNFRSF9",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"HAVCR2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"PDCD1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"ENTPD1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"LAG3",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"TRGV9",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"TRDV2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"XIST",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"CD3D",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD4",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD8A",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD8B",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"FOXP3",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"NCAM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"NCR1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"PMEL",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"MLANA",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"TYR",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"CPXM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"APOE",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"HBA1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"HBA2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"HBB",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"JCHAIN",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD19",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"MS4A1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD14",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

lig.reannotated <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/lig.reannotated.rda")
lig.reannotated <- setNames(as.character(lig.reannotated),names(lig.reannotated))

all.plots.reannotated <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T, clusters = lig.reannotated)


all.plots[[2]] + all.plots.clusters[[2]] + all.plots.reannotated[[2]]

all.plots[[1]] + all.plots.reannotated[[2]]


# Reannotate clusters

seu <- ligerToSeurat(lig.tumor.louvain)
Idents(seu) <- lig.tumor.louvain@clusters
seu@meta.data$clusters <- Idents(seu) 

DimPlot(seu,label = T)
DotPlot(seu,features = c(
  "CD3D","CD8A","CD4","FOXP3","NCAM1","NCR1","PMEL","MLANA","TYR",
  "HBA1","HBA2","HBB","JCHAIN","CD19","MS4A1","CD14","TNFRSF21",
  "ITGAX","MZB1","TNFRSF9",
  "HAVCR2","ENTPD1","PDCD1","TRGV9","TRDV2"
),assay = "RNA",group.by = "clusters") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

seu <- RenameIdents(seu,
                    `0` = "Melanoma",
                    `1` = "Melanoma",
                    `2` = "Melanoma",
                    `3` = "Melanoma",
                    `4` = "Melanoma",
                    `5` = "CD8 T cell",
                    `6` = "CD4 T cell",
                    `7` = "CD8 T cell (activated/exhausted)",
                    `8` = "CD8 T cell (activated/exhausted)",
                    `9` = "CD8 T cell",
                    `10` = "NK cell",
                    `11` = "CD8 T cell",
                    `12` = "CD8 T cell (activated/exhausted)",
                    `13` = "CD8 T cell (activated/exhausted)",
                    `14` = "CD8 T cell",
                    `15` = "Melanoma",
                    `16` = "Monocyte / DC",
                    `17` = "NK cell",
                    `18` = "CD8 T cell (activated/exhausted)",
                    `19` = "Melanoma",
                    `20` = "CD4 T cell (Treg)",
                    `21` = "Melanoma (CPXM1+)",
                    `22` = "Melanoma",
                    `23` = "pDC (CD4+)",
                    `24` = "B cell",
                    `25` = "pDC (CD4-)",
                    `26` = "Mixed",
                    `27` = "Mixed",
                    `28` = "CD8 T cell (activated/exhausted)",
                    `29` = "Melanoma"
                    )
seu@meta.data$clusters.annotated <- Idents(seu)

DimPlot(seu,label = T)
DimPlot(seu,label = T,group.by = "clusters")

s <- subset(seu,clusters.annotated != "Mixed")
DimPlot(s,label = T)
DimPlot(s,label = T,group.by = "clusters")

DotPlot(s,features = c(
  "CD3D","CD8A","CD4","FOXP3","NCAM1","NCR1","PMEL","MLANA","TYR","CPXM1",
  "HBA1","HBA2","HBB","JCHAIN","CD19","MS4A1","CD14","TNFRSF21",
  "ITGAX","MZB1","TNFRSF9",
  "HAVCR2","ENTPD1","PDCD1"
),assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(s,features = c(
  "CD3D","CD4","CD8A","CD8B","TRAC","TRBC2","CTLA4","FOXP3",
  "CCL5","NKG7","GZMA","TNFRSF4","IL7R","TCF7","IFNG","TIGIT",
  "LAG3","PDCD1","KIR2DL4","ISG15","ITM2A","GZMB","CCL3","CCL4",
  "CXCL13","GZMK","GZMH","KLRB1","CCR5","PRF1","IL32","FGFBP2",
  "FCGR3A","CX3CR1","CD5","GNLY"
),assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

s.tnk <- subset(s,clusters.annotated %in% 
                  c("CD8 T cell",
                    "CD8 T cell (activated/exhausted)",
                    "CD4 T cell",
                    "CD4 T cell (Treg)",
                    "NK cell"))

cells_remove <- names(which(s.tnk@assays$RNA@counts["MLANA",] > 0 | 
  s.tnk@assays$RNA@counts["PMEL",] > 0 | 
  s.tnk@assays$RNA@counts["TYR",] > 0 | 
  s.tnk@assays$RNA@counts["CPXM1",] > 0 | 
  s.tnk@assays$RNA@counts["CD19",] > 0 | 
  s.tnk@assays$RNA@counts["CD14",] > 0))
cells_keep <- setdiff(colnames(s.tnk),cells_remove)


s.tnk <- subset(s.tnk, cells = cells_keep)

#AOAH
#GZMH
#CD160

cells_keep_liger <- str_split_fixed(cells_keep,"GEX_",2)[,2]

lig.tnk <- subsetLiger(lig.tumor.louvain,cells.use = cells_keep_liger)

annot <- setNames(as.character(s.tnk@meta.data$clusters.annotated),
                  rownames(s.tnk@meta.data))
names(annot) <- str_split_fixed(names(annot),"GEX_",2)[,2]
lig.tnk <- runUMAP(lig.tnk, 
                           distance = "cosine", dims.use=dims.use, 
                           n_neighbors = 20, min_dist = 0.05)

lig.tnk.plots <- plotByDatasetAndCluster(
  lig.tnk, axis.labels = c('UMAP 1', 'UMAP 2'), 
  return.plots = T, clusters = annot)
lig.tnk.plots[[2]]

lig.tnk.plots.cluster <- plotByDatasetAndCluster(
  lig.tnk, axis.labels = c('UMAP 1', 'UMAP 2'), 
  return.plots = T)
lig.tnk.plots.cluster[[2]]


s.tnk.reumap <- ligerToSeurat(lig.tnk)
stopifnot(identical(names(Idents(s.tnk.reumap)),names(Idents(s.tnk))))
Idents(s.tnk.reumap) <- Idents(s.tnk)
stopifnot(identical(rownames(s.tnk.reumap@meta.data),
                    rownames(s.tnk@meta.data)))
s.tnk.reumap@meta.data$clusters <- s.tnk@meta.data$clusters

DimPlot(s.tnk.reumap)
DimPlot(s.tnk.reumap, group.by = "clusters", label = T)

s.tnk.reumap@meta.data$clusters.annotated <- Idents(s.tnk.reumap)
Idents(s.tnk.reumap) <- s.tnk.reumap@meta.data$clusters
markers.all <- FindAllMarkers(s.tnk.reumap)


cells_keep <- rownames(s.tnk.reumap@meta.data)[which(s.tnk.reumap@meta.data$clusters.annotated %in% 
  c("CD8 T cell (activated/exhausted)","CD8 T cell"))]

s.t <- subset(s.tnk,cells = cells_keep)

cells_keep <- str_split_fixed(cells_keep,"GEX_",2)[,2]

lig.t <- subsetLiger(lig.tnk,cells.use = cells_keep)
lig.t <- runUMAP(lig.t, distance = "cosine", dims.use=dims.use, 
                  n_neighbors = 20, min_dist = 0.05)

lig.t.plots <- plotByDatasetAndCluster(
  lig.t, axis.labels = c('UMAP 1', 'UMAP 2'), 
  return.plots = T, clusters = annot)
lig.t.plots[[2]]


lig.t.clusters.plots <- plotByDatasetAndCluster(
  lig.t, axis.labels = c('UMAP 1', 'UMAP 2'), 
  return.plots = T)
lig.t.clusters.plots[[2]]



s.t.reumap <- ligerToSeurat(lig.t)
stopifnot(identical(names(Idents(s.t.reumap)),names(Idents(s.t))))
Idents(s.t.reumap) <- Idents(s.t)
stopifnot(identical(rownames(s.t.reumap@meta.data),
                    rownames(s.t@meta.data)))
s.t.reumap@meta.data$clusters <- s.t@meta.data$clusters

s.t.reumap@meta.data$clusters.annotated <- Idents(s.t.reumap)
Idents(s.t.reumap) <- s.t.reumap@meta.data$clusters
s.t <- s.t.reumap


markers.t <- FindAllMarkers(s.t,assay = "RNA",slot = "data",only.pos = T)
markers.t.roc <- FindAllMarkers(s.t,assay = "RNA",slot = "data",
                                test.use = "roc",only.pos = T)


lig.t.clusters.plots[[2]] + DotPlot(s.t,features = 
                        c("CD27","FAS","FASLG","ITGAL",
                         "ITGB2","ITGB1","IL2","SELL",
                         "ITGAM","ITGA4","PRF1","TNF",
                         "CCR7","ITGAE","CD69","CXCR3","CXCR6",
                         "CCR5","CCR6","CX3CR1","ITGAV","CD44",
                         "ADGRE5","VCAM1",
                         "CCL20","XCL1","IL7R","LITAF","NR4A2",
                         "AOAH","FCRL6","CCL5","CD8A","CD8B",
                         "TRAC","TRBC2","NKG7","IFNG","PDCD1",
                         "ITM2A","KIR2DL4","GZMB","CCL4","CCL3",
                         "CXCL13",
                         "CD160","GZMH","GZMK",
                         "LAG3","TIGIT","TNFRSF9","GNLY",
                         "FGFBP2","LTB","KLRF1","TCF7",
                         "FCGR3A"),
        cluster.idents = T,assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

lig.t.clusters.plots[[2]] + DotPlot(s.t,features = 
    c("CD3D",
      "GNLY","FGFBP2","FCGR3A","CX3CR1",
      "GZMH","PRF1","NKG7","GZMA","GZMB",
      "CCL5","CD8A","CD8B","TRAC","TRBC2","IL32",
      "IFNG","TIGIT","LAG3","PDCD1","KIR2DL4",
      "ITM2A","CCL3","CCL4","CXCL13","CTLA4",
      "CCR5","GZMK","CD5","ISG15","TCF7","IL7R","KLRB1",
      "CD4","FOXP3","TNFRSF4",
      
      "AOAH","FCRL6","CD160","JUNB","KLRD1","JUN",
      "LTB","TPT1","TXNIP","SLC4A10","KLRG1","LTK","NCR3",
      "XIST","DDX17","PARP8","SMCHD1","MACF1","CELF2","SYNE2","IKZF1","ETS1","ZFP36L2",
      "TNFRSF9","EGR2","ZBED2","TBC1D4","AHI1","RNF19A","CD200","LYST","GNG4","ENTPD1","CD7",
      "STMN1","MKI67","MCM2",
      "CD27","SIRPG","CD74","PTMS","CST7","HAVCR2","CD81",
      "GINS2","TYMS","E2F1"
      ),
  cluster.idents = T,assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


lig.t.clusters.plots[[2]] + DotPlot(s.t,features = 
                                      c("CD3D",
                                        "GNLY","FGFBP2","FCGR3A","CX3CR1",
                                        "GZMH","PRF1","NKG7","GZMA","GZMB",
                                        "CCL5","CD8A","CD8B","TRAC","TRBC2","IL32",
                                        "IFNG","TIGIT","LAG3","PDCD1","KIR2DL4",
                                        "ITM2A","CCL3","CCL4","CXCL13","CTLA4",
                                        "CCR5","GZMK","CD5","ISG15","TCF7","IL7R","KLRB1",
                                        "CD4","FOXP3","TNFRSF4",
                                        "AOAH","FCRL6","CD160",
                                        "EGR2","XCL2","TNFRSF9","NFKBID","NR4A1","XCL1","GNG4","TNFSF14","GOLIM4",
                                        "CD27","HAVCR2",
                                        "KLRC2","ZNF683","CD38","CD244",
                                        "TMSB4X","LTB","EEF1A1","RPL34","KLF2",
                                        "CD69","HLA-DRA",
                                        "MKI67"
                                      ),
                                    cluster.idents = T,assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(s.t,features = c("CD69","PDCD1","LAG3"),ncol = 3)


#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5438859/:
# HLA-DR ... appears at the late stages of activation on T and NK cells; thus, it is considered to be a very late activation marker
# The earliest activation marker is CD69, which is an inducible cell surface glycoprotein expressed upon activation via the TCR or the IL-2 receptor (CD25)

markers.t_13_vs_7_8 <- FindMarkers(
  s.t,ident.1 = c("13"),ident.2 = c("7","8"),
  assay = "RNA",slot = "data",only.pos = F)

markers.t_5_vs_9 <- FindMarkers(
  s.t,ident.1 = c("5"),ident.2 = c("9"),
  assay = "RNA",slot = "data",only.pos = F)


# XLC2, GNG4, TNFSF14, GOLIM4

#LTB
#EEF1A1
#RPL34
#KLF2

#GZMH
#TMSB4X

#KLRC2
#ZNF683
#CD38
#CD244


# 14: "Cytotoxic" (Li et al). Antigen-inexperienced?
# 5: Naive-like
# 11: Naive-like
# 9: GZMH high, KLRC2+, ZNF683+, CD38+, CD244+. Early activated
# 13: Antigen-experienced, activated/exhausted/dysfunctional (early)
# 12: Antigen-experienced, activated/exhausted/dysfunctional (early, proliferative)
# 18: Antigen-experienced, activated/exhausted/dysfunctional (early, proliferative)
# 7: Antigen-experienced, activated/exhausted/dysfunctional (late)
# 8: Antigen-experienced, activated/exhausted/dysfunctional (late)
# 28: Antigen-experienced, activated/exhausted/dysfunctional (late, proliferative)



s.t <- RenameIdents(s.t, 
                    `9` = "Early activated (KLRC2+, CD38+)",
                    `5` = "Naive-like / early activated",
                    `11` = "Naive-like", 
                    `14` = "Activated, cytotoxic (FCGR3A+)",
                    `13` = "Late activated",
                    `12` = "Late activated, proliferative",
                    `18` = "Late activated, proliferative",
                    `7` = "Exhausted / dysfunctional",
                    `8` = "Exhausted / dysfunctional",
                    `28` = "Exhausted / dysfunctional, proliferative")






gene_loadings <- plotGeneLoadings(lig.tnk, do.spec.plot = T, return.plots = TRUE,
                                  factor.share.thresh = Inf)

i <- 47
gene_loadings[[i]]
g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10)

# 8: first NK cluster
# 26: GZMK, etc. Separates upper half
# 33: IL2RB, lower CD8 cluster and second NK
# 47: GZMH, lower CD8 cluster and lower CD4
# 50: Upper, exhausted CD8 cluster. CD27, GZMK, NKG7, LAG3



