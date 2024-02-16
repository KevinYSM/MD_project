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

mark_9_cd8 <- FindMarkers(s.tnk,group.by = "clusters",
                        ident.1 = 9,ident.2 = c(5,11,13,12,18,7,8,28,14))
mark_9_cd8[mark_9_cd8$p_val_adj < 0.05 & mark_9_cd8$avg_log2FC > 0 & 
             mark_9_cd8$pct.1 > (mark_9_cd8$pct.2+0.1),]

m <- mark_9_cd8[mark_9_cd8$p_val_adj < 0.05 & mark_9_cd8$avg_log2FC > 0 & 
             mark_9_cd8$pct.1 > (mark_9_cd8$pct.2+0.1),]

m <- mark_9_5[mark_9_5$p_val_adj < 0.05 & mark_9_5$avg_log2FC > 0 & 
                  mark_9_5$pct.1 > (mark_9_5$pct.2+0.1),]

m[order(m$avg_log2FC,decreasing = T),]

DimPlot(s.tnk)
FeaturePlot(s.tnk,
            features = c("PMEL","MLANA","TYR","CD3D","CD8A","CD4",
                         "CPXM1","CD19","CD14"))


FeaturePlot(s.tnk,
            features = c("PMEL","MLANA","TYR","CD3D","CD8A","CD4",
                         "CPXM1","CD19","CD14"))

AOAH
GZMH
CD160

cells_keep_liger <- str_split_fixed(cells_keep,"GEX_",2)[,2]

lig.tnk <- subsetLiger(lig.tumor.louvain,cells.use = cells_keep_liger)

lig.tnk.plots.reannotated <- plotByDatasetAndCluster(
  lig.tnk, axis.labels = c('UMAP 1', 'UMAP 2'), 
  return.plots = T, clusters = lig.reannotated)

lig.tnk.louvain <- louvainCluster(lig.tnk, resolution = 2,
                                  eps = 0.001)

lig.tnk.plots.louvain <- plotByDatasetAndCluster(
  lig.tnk.louvain, axis.labels = c('UMAP 1', 'UMAP 2'), 
  return.plots = T)
lig.tnk.plots.louvain[[2]]

lig.tnk.louvain <- runUMAP(lig.tnk.louvain, 
                                 distance = "cosine", dims.use=dims.use, 
                     n_neighbors = 20, min_dist = 0.05)

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

lig.tnk.plots.cluster[[1]] + lig.tnk.plots.cluster[[2]]

#saveRDS(lig.tumor,paste0(outdir,"lig.tumor.reannotated.rda"))
#lig.tumor.bak <- lig.tumor

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




gene_loadings[[i]]
rliger::plotGene(lig.tnk,"CD69",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")










clus <- c(50,2,43,34,26,30)
cells <- names(lig.tumor@clusters)[lig.tumor@clusters %in% clus]
saveRDS(cells,paste0(outdir,"cells.rda"))

saveRDS(lig.tumor@clusters,paste0(outdir,"lig.tumor_clusters.rda"))

stopifnot(identical(names(lig.tumor@clusters),names(celltypes3)))
saveRDS(celltypes3,paste0(outdir,"celltypes3.rda"))
