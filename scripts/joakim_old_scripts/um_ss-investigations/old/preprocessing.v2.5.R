# Reload session ---------------------------------------------------------------
library(Seurat)
library(future)
library(tidyverse)
library(djvdj)
library(plyranges)
library(rtracklayer)
library(DoubletFinder)
library(pheatmap)

#plan("multiprocess", workers = 63)
#plan("multiprocess", workers = 15)
#options(future.globals.maxSize = 8000 * 1024^2)
options(future.globals.maxSize = 80000 * 1024^2)

outdir <- "~/proj/um_ss/Investigations/seurat/results/"

#dat.integrated <- readRDS("~/nimbus/data/proj/um_ss/Investigions/seurat/results/dat.integrated.rda")
#saveRDS(dat.integrated,file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.rda")

dat.integrated <- readRDS("~/proj/um_ss/Investigations/seurat/results/dat.integrated.rda")
#dat.integrated <- readRDS("/data/proj/um_ss/Investigations/seurat/results/dat.integrated.rda")

dat.integrated_RNA <- dat.integrated
DefaultAssay(dat.integrated_RNA) <- "RNA"

# Markers to plot --------------------------------------------------------------
dat.integrated <- readRDS("~/proj/um_ss/Investigations/seurat/results/dat.integrated.rda")

canonical_markers <- c("ITGAL","ITGB1","ITGB2","ITGA4","ITGAE","CXCR3","SELL","IL2RA","IL2RB","CD2","CD5","CD69","CD47","CXCR5","PECAM1")
other_markers <- c("CDKN2A","BRAF","KRAS","NRAS", "VDR", "MITF", "TYR", "PDCD1", "CD274", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "ENTPD1")
markers.seurat_vignette <- union(c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                                   "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                                   "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57"),
                                 c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                                   "CCL2", "PPBP"))
m.to_plot <- unique(c(markers.seurat_vignette,
                      canonical_markers,other_markers, "HBA1" ,"CD4", "NCAM1", "CD3E", "CD14", "CD19", "PMEL","MLANA", "FOXP3"))
markers_plot <- setdiff(c(m.to_plot,"CD8B",
                          grep("^TRDV",rownames(dat.integrated),value=T),
                          grep("^TRGV",rownames(dat.integrated),value=T),
                          "CLEC4C","THBD","NRP1","NRP2","IL12A","IL6","TLR2","TLR4","TLR7","TLR9","CD1C","CCR7",
                          "IRF4","KLF4","ZEB2","MAFB","ITGAX","IRF8","CD83","TRAV10","TRAJ18","CD207","PTPRC",
                          "JCHAIN","CD68","CD1A","AXL","CHGA","CHGB","KIT","SOX10","DCT","EDNRB","WNT5A","FGFR3",
                          "EGFR","NGFR","SOX9","SNAI1","ZEB1","TWIST1","SNAI2","SLC7A8","SLC12A7","DLX5","CD163","XIST"),"IGJ")

saveRDS(markers_plot,file=paste0(outdir,"markers_plot.rda"))

# Basic preprocessing ----------------------------------------------------------
library(Seurat)
library(future)
library(tidyverse)
library(djvdj)
library(plyranges)
library(rtracklayer)

plan("multiprocess", workers = 63)
#plan("multiprocess", workers = 15)

options(future.globals.maxSize = 8000 * 1024^2)

pths <- list()
pths[["G18-023"]] <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/gex/G18-023_*/outs/filtered_feature_bc_matrix")
pths[["G18-049"]] <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/gex/G18-049_*/outs/filtered_feature_bc_matrix")
pths[["GWA-JN-388"]] <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/gex/GWA-JN-388_*/outs/filtered_feature_bc_matrix")

dat <- lapply(pths,function(x){Read10X(data.dir=x, strip.suffix = T)})
dat <- sapply(names(dat),function(nm){CreateSeuratObject(counts = dat[[nm]], project = nm)})
dat <- lapply(dat,function(x){SCTransform(x, vst.flavor = "v2", verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)})

features <- SelectIntegrationFeatures(object.list = dat, nfeatures = 3000)

dat.bak <- dat
dat <- PrepSCTIntegration(object.list = dat, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = dat, 
                                  normalization.method = "SCT",
                                  anchor.features = features)

dat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
dat.integrated <- RunPCA(dat.integrated, verbose = FALSE)
dat.integrated <- RunUMAP(dat.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
dat.integrated <- FindNeighbors(dat.integrated, reduction = "pca", dims = 1:30)
dat.integrated <- FindClusters(dat.integrated, resolution = 0.3)

project <- stringr::str_split_fixed(rownames(dat.integrated@meta.data),"_",3)[,3]
names(project) <- colnames(dat.integrated)
stopifnot(identical(colnames(dat.integrated),rownames(dat.integrated@meta.data)))
project.name <- project
project.name[project.name=="1"] <- "G18-023"
project.name[project.name=="2"] <- "G18-049"
project.name[project.name=="3"] <- "GWA-JN-388"

dat.integrated <- AddMetaData(object = dat.integrated, metadata = project, col.name = "project")
dat.integrated <- AddMetaData(object = dat.integrated, metadata = project.name, col.name = "project.name")

# Add external metadata
dat.integrated.bak <- dat.integrated
df <- data.frame(sample_nr=apply(str_split_fixed(rownames(dat.integrated@meta.data),"_",3)[,c(1,3)],1,
                                 function(x){paste0(x,collapse="_")}))
rownames(df) <- rownames(dat.integrated@meta.data)
df$cell_id <- rownames(df)

library(readxl)
meta <- read_excel("/data/proj/um_ss/Investigations/samples_10x.xlsx")
meta[["Sample ID"]]
meta[["sample_nr"]] <- unique(df$sample_nr)

df <- merge(df,meta,by="sample_nr",all.x=T)
rownames(df) <- df$cell_id
df <- df[rownames(dat.integrated@meta.data),]
stopifnot(identical(rownames(df),rownames(dat.integrated@meta.data)))

if (all(is.na(df$...20))){
  df$...20 <- NULL
}
colnames(df) <- gsub(" ",".",colnames(df))

dat.integrated <- AddMetaData(object = dat.integrated, metadata = df)

dat.integrated <- PercentageFeatureSet(dat.integrated, pattern = "^MT-", col.name = "percent_mito", assay = "RNA")
dat.integrated <- PercentageFeatureSet(dat.integrated, pattern = "^RP[SL]", col.name = "percent_ribo", assay = "RNA")
dat.integrated <- PercentageFeatureSet(dat.integrated, "^HB[^(P)]", col.name = "percent_hb", assay = "RNA")
dat.integrated <- PercentageFeatureSet(dat.integrated, "PECAM1|PF4", col.name = "percent_plat", assay = "RNA")

dat.integrated <- AddMetaData(dat.integrated,
                              apply(dat.integrated@meta.data[,c("seurat_clusters","UM.ID","project.name")],1,
                                    function(x){paste0(x,collapse = ":")}),col.name = "seurat_clusters__UM.ID__project.name")
dat.integrated <- AddMetaData(dat.integrated,
                              apply(dat.integrated@meta.data[,c("UM.ID","project.name")],1,
                                    function(x){paste0(x,collapse = ":")}),col.name = "UM.ID__project.name")
dat.integrated <- AddMetaData(dat.integrated,
                              apply(dat.integrated@meta.data[,c("active.ident","UM.ID","project.name")],1,
                                    function(x){paste0(x,collapse = ":")}),col.name = "active.ident__UM.ID__project.name")

dat.integrated_RNA <- dat.integrated
DefaultAssay(dat.integrated_RNA) <- "RNA"

#saveRDS(dat.integrated,file = paste0(outdir,"/dat.integrated.rda"))

## Random plots -----------------------------------------------------------------

#outdir <- "/data/proj/um_ss/Investigations/seurat/results/"
outdir <- "~/proj/um_ss/Investigations/seurat/results/"
dir.create(outdir,showWarnings = F,recursive = T)

pdf(file = paste0(outdir,"dimplot.integrated.projname.pdf"),width = 20,height = 10)
DimPlot(dat.integrated, reduction = "umap", split.by = "project.name")
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.clusters.pdf"),width = 10,height = 10)
DimPlot(dat.integrated, reduction = "umap", group.by = "seurat_clusters", 
        label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.samples.pdf"),width = 20,height = 10)
DimPlot(dat.integrated, reduction = "umap", group.by = "Sample.ID", split.by = "project.name",
        label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.samples2.pdf"),width = 80,height = 7)
DimPlot(dat.integrated, reduction = "umap", group.by = "seurat_clusters", split.by = "UM.ID",
        label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.patient.pdf"),width = 20,height = 10)
DimPlot(dat.integrated, reduction = "umap", group.by = "UM.ID", split.by = "project.name",
        label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.patient2.pdf"),width = 10,height = 10)
DimPlot(dat.integrated, reduction = "umap", group.by = "UM.ID", label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.sex.pdf"),width = 20,height = 10)
DimPlot(dat.integrated, reduction = "umap", group.by = "project.name", split.by = "Sex",
        label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.project.pdf"),width = 10,height = 10)
DimPlot(dat.integrated, reduction = "umap", group.by = "project.name",
        label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.sex2.pdf"),width = 10,height = 10)
DimPlot(dat.integrated, reduction = "umap", group.by = "Sex",
        label = TRUE, repel = TRUE)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.sex_expr.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("XIST"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.CD8A.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("CD8A"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.CD4.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("CD4"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.CD14.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("CD14"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.CD19.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("CD19"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.NCAM1.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("NCAM1"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.PMEL.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("PMEL"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.MLANA.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("MLANA"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.ALB.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("ALB"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.ARG1.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("ARG1"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.PCK1.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("PCK1"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.AFP.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("AFP"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.BCHE.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("BCHE"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.MARCO.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("MARCO"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.CD5L.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("CD5L"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.VSIG4.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("VSIG4"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.CD68.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("CD68"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.MALAT1.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("MALAT1"), cols = c("white", "red"),max.cutoff = Inf)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.percent_mito.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("percent_mito"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.percent_ribo.pdf"),width = 6,height = 6)
FeaturePlot(dat.integrated_RNA, slot = "counts",features = c("percent_ribo"), cols = c("white", "red"),max.cutoff = 10)
dev.off()

pdf(file = paste0(outdir,"dimplot.integrated.qc.pdf"),width = 15,height = 6)
#feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
#feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
feats <- c("nFeature_RNA", "percent_mito")
VlnPlot(dat.integrated_RNA, group.by = "seurat_clusters", features = feats, pt.size = 0.1, ncol = 3,raster = T) +
  NoLegend()
dev.off()

canonical_markers <- c("ITGAL","ITGB1","ITGB2","ITGA4","ITGAE","CXCR3","SELL","IL2RA","IL2RB","CD2","CD5","CD69","CD47","CXCR5","PECAM1")
for (m in canonical_markers){
  pdf(file = paste0(outdir,"dimplot.integrated.",m,".pdf"),width = 6,height = 6)
  print(FeaturePlot(dat.integrated_RNA, slot = "counts",features = m, cols = c("white", "red"),max.cutoff = 10))
  dev.off()
}

other_markers <- c("CDKN2A","BRAF","KRAS","NRAS", "VDR", "MITF", "TYR", "PDCD1", "CD274", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "ENTPD1")
for (m in other_markers){
  pdf(file = paste0(outdir,"dimplot.integrated.",m,".pdf"),width = 6,height = 6)
  print(FeaturePlot(dat.integrated_RNA, slot = "counts",features = m, cols = c("white", "red"),max.cutoff = 10))
  dev.off()
}

markers <- FindConservedMarkers(dat.integrated_RNA, ident.1 = 0, 
                                grouping.var = c("project.name","Sex"), 
                                verbose = F)
head(markers)

markers <- list()
for (cl in levels(dat.integrated_RNA@meta.data$seurat_clusters)){
  markers[[cl]] <- FindConservedMarkers(dat.integrated_RNA, ident.1 = cl, 
                                        grouping.var = c("project.name","Sex"),
                                        verbose = T)
}

saveRDS(markers,file = paste0(outdir,"/markers.rda"))

filtered_markers <- list()
for (cl in levels(dat.integrated_RNA@meta.data$seurat_clusters)){
  idx_filter <- markers[[cl]]$`G18-049_avg_log2FC` > 0 & 
    markers[[cl]]$`G18-023_avg_log2FC` > 0 & markers[[cl]]$`GWA-JN-388_avg_log2FC` > 0
  filtered_markers[[cl]] <- rownames(markers[[cl]][idx_filter,])
}

lapply(markers,function(x){intersect(rownames(x),canonical_markers)})

lapply(markers,function(x){x[intersect(rownames(x)[x$max_pval < 0.05 & 
  rowSums(x[,grep("_avg_log2FC",colnames(x))] > 0) >= 2,],canonical_markers),]})

markers.findAllMarkers <- FindAllMarkers(dat.integrated_RNA, only.pos = T)
saveRDS(markers.findAllMarkers,file = paste0(outdir,"/markers.findAllMarkers.rda"))

## Cluster annotation (exploration) ---------------------------------------------

# Cluster designations (compare also: https://satijalab.org/seurat/articles/integration_introduction.html):

# 0: CD4 T cell, CD4+, CD2+, CD5+, CD69+, CXCR3+, SELL+, PDCD1 (PD1)+
# 1: Exhausted/dysfunctional CD8 T cell, CD8A+, CD2+, CD5+, CD69+, CXCR3+, LAG3+, TIGIT+, HAVCR2+, PDCD1+, CTLA4+, ENTPD1+
# 2: Melanocytic, MLANA+, PMEL+, MITF+ (high)
# 3: CD8 T cell (GNLY+, FCGR3A+), CD8A+, CD2+, CD5+, CD69+, CXCR3+
# 4: Melanocytic, MLANA+, PMEL+, TYR+, MITF+
# 5: Melanocytic/Lymphocyte doublets, MLANA+, PMEL+, CD3E+, NKG7+, high mitochondrial content, low detected gene counts
# 6: NK cells, NCAM1+, CD2+, CD69+, ITGB2+ (high)
# 7: Monocytes/macrophages/dendritic cells, CD14+, CD4+, ITGB2+ (high), PECAM1+, SELL+, CD4, S100A9+ (high)
# 8: B cell/lymphocyte doublet, CD19+, CD69+, CXCR5+, SELL+, CD3D+, MS4A1+, CD79A+, ITGA4+, HLA-DQA1+, CD3E+
# 9: Exhausted/dysfunctional CD8 T cell, CD8A+, CD2+, CD5+, CD69+, CXCR3+, LAG3+, TIGIT+, HAVCR2+, PDCD1+, CTLA4+, ENTPD1+
# 10: Melanocytic/failed cells, PMEL+, MLANA+, high mitochondrial content, low detected gene counts
# 11: CD8 T cell, CD8A+, CXCR3+
# 12: Melanocytic, MLANA+, PMEL+, CD274 (PD-L1)+, TYR+ (high), MITF+ (high), CDKN2A+ (high)
# 13: CD8 T cell, CD3D, NKG7, CCL5, ITGA4, CD3E, CD8A
# 14: Erythrocytes, HBB, HBA1, HBA2
# 15: pDC: IL3RA+, TSPAN13+, CXCR3+ (high), SELL+, HLA-DQA1+, GPR183+, CD4+
# 16: B cells, SELL+, MS4A1 (CD20)+, CD79A+, HLA-DQA1+, ITGA4+, CD69+, CD19+
# 17: Macrophages, CCL5, PPGP, GNG11, ITGB1

canonical_markers <- c("ITGAL","ITGB1","ITGB2","ITGA4","ITGAE","CXCR3","SELL","IL2RA","IL2RB","CD2","CD5","CD69","CD47","CXCR5","PECAM1")
other_markers <- c("CDKN2A","BRAF","KRAS","NRAS", "VDR", "MITF", "TYR", "PDCD1", "CD274", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "ENTPD1")
markers.seurat_vignette <- union(c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57"),
                     c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                       "CCL2", "PPBP"))

m.to_plot <- unique(c(markers.seurat_vignette,
                   canonical_markers,other_markers, "HBA1" ,"CD4", "NCAM1", "CD3E", "CD14", "CD19", "PMEL","MLANA", "FOXP3"))

pdf(file=paste0(outdir,"cell_type_markers_annotation.pdf"),width = 20,height = 10)
DotPlot(dat.integrated_RNA, features = m.to_plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
dev.off()

# Cluster annotation (actions) -------------------------------------------------

dat.integrated <- RenameIdents(dat.integrated, 
                                `0` = "CD4 T cell", 
                                `1` = "CD8 T cell (Exhausted/dysfunctional)", 
                                `2` = "Melanocytic (MITF high)",
                                `3` = "CD8 T cell (GNLY+, FCGR3A+)", 
                                `4` = "Melanocytic", 
                                `5` = "Melanocytic/Lymphocyte doublets, poor quality", 
                                `6` = "NK cells", 
                                `7` = "Monocytes/macrophages/dendritic cells (ITGB2, S100A9 high)", 
                                `8` = "B cell/lymphocyte doublet", 
                                `9` = "CD8 T cell (Exhausted/dysfunctional)",
                                `10` = "Melanocytic (poor quality)", 
                                `11` = "CD8 T cell", 
                                `12` = "Melanocytic (TYR, MITF, CDKN2A high)", 
                                `13` = "CD8 T cell", 
                                `14` = "Erythrocytes",
                                `15` = "pDC",
                                `16` = "B cells",
                                `17` = "Macrophages")

## Cluster annotation (further exploration) -------------------------------------

DimPlot(dat.integrated, label = TRUE, split.by = "project.name")

dat <- subset(x = dat.integrated, subset = project.name == "GWA-JN-388")

pdf(file=paste0(outdir,"GWA-JN-388_vs_UM_ID.pdf"),width = 50,height = 10)
DimPlot(dat, label = TRUE, split.by = "UM.ID", repel=T, label.size=3)
dev.off()

pdf(file=paste0(outdir,"all_by_UM_ID.pdf"),width = 90,height = 10)
DimPlot(dat.integrated, label = TRUE, split.by = "UM.ID", repel=T, label.size=3)
dev.off()

ggplot(dat.integrated_RNA@meta.data, aes(x=nCount_RNA, y=percent_mito, color = active.ident)) + geom_point()
ggplot(dat.integrated_RNA@meta.data, aes(x=log2(nCount_RNA), y=percent_mito, color = active.ident)) + geom_point()
ggplot(dat.integrated_RNA@meta.data, aes(x=log2(nCount_RNA), y=percent_mito, color = seurat_clusters)) + geom_point()
ggplot(dat.integrated_RNA@meta.data, aes(x=nCount_RNA, y=percent_mito, color = seurat_clusters, fill = seurat_clusters)) + geom_density2d()
ggplot(dat.integrated_RNA@meta.data, aes(x=nFeature_RNA, y=percent_mito, color = active.ident)) + geom_point()

# Investigate sex (actions) ----------------------------------------------------
gtf <- import("~/nimbus/data/local/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
#gtf <- import("/data/local/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

chrX.gene <- unique((gtf |> filter(seqnames == "chrX"))$gene_name)
chrY.gene <- unique((gtf |> filter(seqnames == "chrY"))$gene_name)

pct_chrX = colSums(dat.integrated_RNA@assays$RNA@counts[chrX.gene, ])/colSums(dat.integrated_RNA@assays$RNA@counts)
pct_chrY = colSums(dat.integrated_RNA@assays$RNA@counts[chrY.gene, ])/colSums(dat.integrated_RNA@assays$RNA@counts)

dat.integrated <- AddMetaData(dat.integrated, metadata = pct_chrX, col.name = "pct_chrX")
dat.integrated <- AddMetaData(dat.integrated, metadata = pct_chrY, col.name = "pct_chrY")

## Investigate sex (exploration) ------------------------------------------------

FeatureScatter(dat.integrated, feature1 = "XIST", feature2 = "pct_chrY", group.by = "UM.ID")
FeatureScatter(dat.integrated, feature1 = "pct_chrX", feature2 = "pct_chrY", group.by = "UM.ID")

FeatureScatter(dat.integrated, feature1 = "XIST", feature2 = "pct_chrY", group.by = "UM.ID__project.name")
FeatureScatter(dat.integrated, feature1 = "pct_chrX", feature2 = "pct_chrY", group.by = "UM.ID__project.name")

VlnPlot(dat.integrated, features = c("XIST", "pct_chrY"), group.by = "seurat_clusters")

# UM1 looks suspicious, but shares TCRs

pdf(file = paste0(outdir,"X_vs_Y.pdf"), width = 6, height = 6)
DotPlot(dat.integrated, features = c("XIST","pct_chrY"), cols = c("blue", "red"),group.by = "UM.ID__project.name") + RotatedAxis()
dev.off()
pdf(file = paste0(outdir,"X_vs_Y_2.pdf"), width = 6, height = 6)
DotPlot(dat.integrated, features = c("pct_chrX","pct_chrY"), cols = c("blue", "red"),group.by = "UM.ID__project.name") + RotatedAxis()
dev.off()

dat.integrated.male <- subset(x = dat.integrated, subset = Sex == "M")
dat.integrated.female <- subset(x = dat.integrated, subset = Sex == "F")

DotPlot(dat.integrated.male, features = c("pct_chrX","pct_chrY"), cols = c("blue", "red"),
        group.by = "UM.ID__project.name") + RotatedAxis()

DotPlot(dat.integrated.female, features = c("pct_chrX","pct_chrY"), cols = c("blue", "red"),
        group.by = "UM.ID__project.name") + RotatedAxis()

# Add VDJ data (actions) -------------------------------------------------------
# library(djvdj)
# fix(import_vdj)
# Add the following, after "contigs <- .load_vdj_data(vdj_dir)":
# contigs <- lapply(contigs,function(x){
#   x[,c("v_gene", "d_gene", "j_gene", "c_gene")][is.na(x[,c("v_gene", "d_gene", "j_gene", "c_gene")])] <- "None"
#   return(x)
# })

pths <- list()
pths[["G18-023"]] <- Sys.glob("~/nimbus/data/proj/um_ss/Pipelines/10x/results/vdj/G18-023_*/outs")
pths[["G18-049"]] <- Sys.glob("~/nimbus/data/proj/um_ss/Pipelines/10x/results/vdj/G18-049_*/outs")
pths[["GWA-JN-388"]] <- Sys.glob("~/nimbus/data/proj/um_ss/Pipelines/10x/results/vdj/GWA-JN-388_*/outs")

names(pths[["G18-023"]]) <- paste0(1:length(pths[["G18-023"]]),"_")
names(pths[["G18-049"]]) <- paste0(1:length(pths[["G18-049"]]),"_")
names(pths[["GWA-JN-388"]]) <- paste0(1:length(pths[["GWA-JN-388"]]),"_")

vdj <- list()
for (nm in names(pths)) {
  vdj[[nm]] <- import_vdj(input = NULL, vdj_dir = pths[[nm]],
                          filter_chains = T, define_clonotypes = "cdr3_gene", 
                          include_indels = T)
}

rownames(vdj[["G18-023"]]) <- gsub("-1$","_1",rownames(vdj[["G18-023"]]))
rownames(vdj[["G18-049"]]) <- gsub("-1$","_2",rownames(vdj[["G18-049"]]))
rownames(vdj[["GWA-JN-388"]]) <- gsub("-1$","_3",rownames(vdj[["GWA-JN-388"]]))

vdj <- do.call("rbind", vdj)
rownames(vdj) <- str_split_fixed(rownames(vdj),"\\.",2)[,2]

dat.integrated <- AddMetaData(dat.integrated, metadata = vdj)
dat.integrated <- AddMetaData(dat.integrated, metadata = dat.integrated@active.ident, col.name = "active.ident")

pdf(paste0(outdir,"umap.cell_types.pdf"),width = 20,height = 14)
DimPlot(dat.integrated, label = T, group.by = "active.ident")
dev.off()

## Investigate VDJ sequences (exploration) --------------------------------------
dat.integrated@meta.data$cdr3_top100 <- dat.integrated@meta.data$cdr3
dat.integrated@meta.data$cdr3_top100[which(! dat.integrated@meta.data$cdr3_top100 %in% 
                                             names(sort(table(dat.integrated@meta.data$cdr3),decreasing = T)[1:100]))] <- NA

dat.integrated@meta.data$cdr3_top1000 <- dat.integrated@meta.data$cdr3
dat.integrated@meta.data$cdr3_top1000[which(! dat.integrated@meta.data$cdr3_top1000 %in% 
                                              names(sort(table(dat.integrated@meta.data$cdr3),decreasing = T)[1:1000]))] <- NA

dat.integrated <- calc_abundance(dat.integrated, cluster_col = "seurat_clusters", clonotype_col = "cdr3", prefix = "djvdj.")
dat.integrated <- calc_diversity(dat.integrated, cluster_col = "seurat_clusters", clonotype_col = "cdr3", prefix = "djvdj.")
dat.integrated <- calc_similarity(dat.integrated, cluster_col = "seurat_clusters", clonotype_col = "cdr3", prefix = "djvdj.similarity.")
#saveRDS(dat.integrated,file = paste0(outdir,"/dat.integrated.rda"))

vj_gene_usage <- calc_gene_usage(dat.integrated, cluster_col = "seurat_clusters", gene_cols = c("v_gene", "j_gene"), chain = "TRB")
cluster_clonotype_similarity <- calc_similarity(dat.integrated, cluster_col = "seurat_clusters", clonotype_col = "cdr3", prefix = "djvdj.similarity.", return_mat = T)

t_size <- 5
l_size <- 0.1
theme_text_small <- theme(text = element_text(size = t_size),
                          axis.text = element_text(size = t_size),
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                          axis.title = element_text(size = t_size),
                          legend.text = element_text(size = t_size),
                          title = element_text(size = t_size),
                          legend.title = element_blank(),
                          plot.title = element_text(size = t_size),
                          plot.subtitle = element_text(size = t_size),
                          strip.text = element_text(size = t_size, angle = 45, vjust = 0,hjust = 0),
                          legend.key.size = unit(l_size,"line"))
  
#plot_abundance(dat.integrated, cluster_col = "seurat_clusters", yaxis = "frequency") + theme_text_small
#plot_abundance(dat.integrated, cluster_col = "project.name", yaxis = "frequency") + theme_text_small
pdf(file = paste0(outdir,"clonotypes_vs_cell_type.pdf"), width = 15, height = 5)
plot_abundance(dat.integrated, cluster_col = "active.ident", yaxis = "frequency", clonotype_col = "cdr3") + theme_text_small
dev.off()

#plot_diversity(dat.integrated, cluster_col = "seurat_clusters__UM.ID__project.name", clonotype_col = "clonotype_id") + 
#  theme_text_small + facet_wrap(~ str_split_fixed(seurat_clusters__UM.ID__project.name,":",3)[,2], scales = "free_x")

dat.integrated.t <- subset(dat.integrated,
                           cells = rownames(dat.integrated@meta.data[dat.integrated@meta.data$active.ident %in% 
                            c("CD4 T cell", "CD8 T cell", "CD8 T cell (Exhausted/dysfunctional)", 
                              "CD8 T cell (GNLY+, FCGR3A+)"),]))

pdf(file = paste0(outdir,"diversity.samples_vs_cluster_vs_project_tcellclusters_only.pdf"), width = 8, height = 6)
plot_diversity(dat.integrated.t, cluster_col = "seurat_clusters__UM.ID__project.name", clonotype_col = "cdr3") + 
  theme_text_small + facet_wrap(~ str_split_fixed(seurat_clusters__UM.ID__project.name,":",3)[,2], scales = "free_x",nrow = 2)
dev.off()

plot_diversity(dat.integrated.t, cluster_col = "active.ident", clonotype_col = "cdr3") + theme_text_small

pdf(file = paste0(outdir,"diversity.active_ident.pdf"), width = 8, height = 6)
plot_diversity(dat.integrated.t, cluster_col = "active.ident__UM.ID__project.name", clonotype_col = "cdr3") + 
  theme_text_small + facet_wrap(~ str_split_fixed(active.ident__UM.ID__project.name,":",3)[,2], scales = "free_x",nrow = 2)
dev.off()


#plot_diversity(dat.integrated, cluster_col = "seurat_clusters", clonotype_col = "clonotype_id")
#plot_diversity(dat.integrated, cluster_col = "UM.ID", clonotype_col = "clonotype_id")
#plot_diversity(dat.integrated, cluster_col = "project.name", clonotype_col = "clonotype_id")

plot_similarity(dat.integrated, cluster_col = "UM.ID", clonotype_col = "cdr3")
plot_similarity(dat.integrated, cluster_col = "project.name", clonotype_col = "cdr3")
plot_similarity(dat.integrated, cluster_col = "UM.ID__project.name", clonotype_col = "cdr3")

cluster_clonotype_similarity.UM.ID__project.name <- calc_similarity(
  dat.integrated, cluster_col = "UM.ID__project.name", clonotype_col = "cdr3", 
  prefix = "djvdj.similarity.", return_mat = T)

#pheatmap(cluster_clonotype_similarity.UM.ID__project.name,border_color = NA)

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

pheatmap(makeSymm(cluster_clonotype_similarity.UM.ID__project.name),border_color = NA)

pdf(file = paste0(outdir,"similarity.jaccard.hclust.pdf"), width = 6, height = 6)
plot(hclust(as.dist(cluster_clonotype_similarity.UM.ID__project.name)))
dev.off()

pdf(file = paste0(outdir,"similarity.clusters.jaccard.hclust.pdf"), width = 6, height = 6)
plot(hclust(as.dist(cluster_clonotype_similarity)))
dev.off()

cluster_clonotype_similarity.seurat_clusters__UM.ID__project.name <- calc_similarity(
  dat.integrated, cluster_col = "seurat_clusters__UM.ID__project.name", clonotype_col = "cdr3", 
  prefix = "djvdj.similarity.", return_mat = T)

pdf(file = paste0(outdir,"similarity.clusters.um_id_project.jaccard.hclust.pdf"), width = 15, height = 6)
plot(hclust(as.dist(cluster_clonotype_similarity.seurat_clusters__UM.ID__project.name)),cex=0.35)
dev.off()

DimPlot(dat.integrated, label = F, group.by = "cdr3_top100") + theme(legend.position="none")

pdf(paste0(outdir,"umap.cdr3.pdf"),width = 10,height = 10)
DimPlot(dat.integrated, label = F, group.by = "cdr3") + theme(legend.position="none")
dev.off()

table(dat.integrated@meta.data$cdr3,dat.integrated@meta.data$seurat_clusters)

pheatmap(table(dat.integrated@meta.data$cdr3_top100,dat.integrated@meta.data$seurat_clusters),show_rownames = F,border_color = NA)
pheatmap(table(dat.integrated@meta.data$cdr3_top1000,dat.integrated@meta.data$seurat_clusters),show_rownames = F,border_color = NA)
pheatmap(table(dat.integrated@meta.data$cdr3_top100,dat.integrated@meta.data$UM.ID__project.name),show_rownames = F,border_color = NA)

x <- table(dat.integrated@meta.data$cdr3,dat.integrated@meta.data$UM.ID__project.name)
x.shared <- x[which(rowSums(x>0)>1),]
dim(x.shared)
pheatmap(x.shared, show_rownames = F)

x.shared.binary <- x.shared
x.shared.binary[x.shared.binary>1] <- 1
pdf(paste0(outdir,"shared.cdr3.binary.pdf"),width = 10,height = 10)
pheatmap(x.shared.binary, show_rownames = F)
dev.off()


# Data cleaning as prescribed by DoubletFinder instructions (see github) (actions) -------
poor_quality_cluster <- setNames(object = dat.integrated@meta.data$active.ident %in% 
                                   c("Melanocytic/Lymphocyte doublets, poor quality",
                                     "Melanocytic (poor quality)"), 
                                 rownames(dat.integrated@meta.data))
dat.integrated <- AddMetaData(dat.integrated, metadata = poor_quality_cluster, 
                              col.name = "poor_quality_cluster")

umi_threshold <- 900
too_low_umi <- setNames(dat.integrated@meta.data$nCount_RNA < umi_threshold, 
                        nm = rownames(dat.integrated@meta.data))
dat.integrated <- AddMetaData(dat.integrated, metadata = too_low_umi, col.name = "too_low_umi")

mito_threshold <- 30
too_high_mito <- setNames(dat.integrated@meta.data$percent_mito > mito_threshold, 
                          nm = rownames(dat.integrated@meta.data))
dat.integrated <- AddMetaData(dat.integrated, metadata = too_high_mito, col.name = "too_high_mito")

#saveRDS(dat.integrated,file="/data/proj/um_ss/Investigations/seurat/results/dat.integrated.rda")

dat.integrated.filtered <- subset(dat.integrated, too_low_umi==F & too_high_mito==F)
rm(dat.integrated)
gc()

## Data cleaning as prescribed by DoubletFinder instructions (see github) (exploration) -------

# Identifying a "nUMI" threshold (now called "nCount_RNA")

#table(dat.integrated@meta.data$poor_quality_cluster,
#      dat.integrated@meta.data$Sample.ID)

VlnPlot(object = dat.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"))

x <- subset(dat.integrated, nCount_RNA < 2000)
g <- FeatureScatter(x,feature1 = "nCount_RNA", feature2 = "percent_mito")
g + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  )

ggplot(x@meta.data, aes(x=nCount_RNA, y=percent_mito) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlim(c(0,2000)) + geom_vline(xintercept = 900) + 
  theme(
    legend.position='none'
  )

ggplot(subset(x,poor_quality_cluster==T)@meta.data, aes(x=nCount_RNA, y=percent_mito) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlim(c(0,2000)) + geom_vline(xintercept = 900) + 
  theme(
    legend.position='none'
  )

ggplot(subset(x,poor_quality_cluster==F)@meta.data, aes(x=nCount_RNA, y=percent_mito) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlim(c(0,2000)) + geom_vline(xintercept = 900) + 
  theme(
    legend.position='none'
  )

g <- FeatureScatter(subset(x,poor_quality_cluster==T),feature1 = "nCount_RNA", feature2 = "percent_mito")
g + xlim(c(0,2000)) + geom_vline(xintercept = 900)

g <- FeatureScatter(subset(x,poor_quality_cluster==F),feature1 = "nCount_RNA", feature2 = "percent_mito")
g + xlim(c(0,2000)) + geom_vline(xintercept = 900)

table(poor_quality_cluster = dat.integrated@meta.data$poor_quality_cluster, 
      low_umi = dat.integrated@meta.data$nCount_RNA < 900)

hist(log10(dat.integrated$nFeature_RNA) / log10(dat.integrated$nCount_RNA))

FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_quantile(quantiles = c(0.01, 0.05, 0.25 ,0.5, 0.75, 0.95, 0.99), method = "rqss", lambda = 0.1)

x <- subset(dat.integrated, nCount_RNA > 900 & poor_quality_cluster==F)

FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_quantile(quantiles = c(0.01, 0.05, 0.25 ,0.5, 0.75, 0.95, 0.99), method = "rqss", lambda = 0.1)


FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & poor_quality_cluster==F) + 
  geom_quantile(quantiles = c(0.01, 0.05, 0.25 ,0.5, 0.75, 0.95, 0.99), method = "rqss", lambda = 0.1)

FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample.ID") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & poor_quality_cluster==F) + 
  geom_quantile(quantiles = c(0.01, 0.05, 0.25 ,0.5, 0.75, 0.95, 0.99), method = "rqss", lambda = 0.1)


FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "seurat_clusters") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & poor_quality_cluster==F) + 
  geom_quantile(quantiles = c(0.01, 0.05, 0.25 ,0.5, 0.75, 0.95, 0.99), method = "rqss", lambda = 0.1) + 
  geom_vline(xintercept = 900)


FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "percent_mito", group.by = "seurat_clusters") + 
  scale_x_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & poor_quality_cluster==F) + 
  geom_vline(xintercept = 900) + geom_hline(yintercept = 30)

x <- subset(dat.integrated, nCount_RNA <= 900 | poor_quality_cluster==T)

FeatureScatter(x, feature1 = "nCount_RNA", feature2 = "percent_mito") + 
  scale_x_log10() + 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  geom_vline(xintercept = 900) + geom_hline(yintercept = 30)

hist(x@meta.data$percent_mito[x@meta.data$nCount_RNA > 1000 & x@meta.data$nCount_RNA < 3000])

FeatureScatter(subset(x, nCount_RNA > 1000 & nCount_RNA < 3000), feature1 = "nCount_RNA", feature2 = "percent_mito") + 
  scale_x_log10() + 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=1) +
  geom_vline(xintercept = 900) + geom_hline(yintercept = 30)

FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "percent_ribo", 
               group.by = "seurat_clusters", plot.cor = F) + 
  scale_x_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & poor_quality_cluster==F) + 
  geom_vline(xintercept = 900) + geom_hline(yintercept = 2)

# Mature RBCs are characterized by a lack of a nucleus and other cell organelles. In addition to being enucleate, RBCs are devoid of RNA, ribosomes, and the capacity for de novo protein synthesis. There is also a loss of mitochondria, leaving RBC dependent on glycolysis for ATP and energy production

FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "percent_ribo", 
               group.by = "active.ident", plot.cor = F) + 
  scale_x_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & poor_quality_cluster==F) + 
  geom_vline(xintercept = 900) + geom_hline(yintercept = 2)


FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "percent_ribo", 
               group.by = "active.ident", plot.cor = F) + 
  scale_x_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & poor_quality_cluster==F) + 
  geom_vline(xintercept = 900) + geom_hline(yintercept = 2)

ggplot(dat.integrated@meta.data,aes(x=active.ident,y=percent_ribo)) + geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept = 5)

ggplot(dat.integrated@meta.data,aes(x=active.ident,y=percent_mito)) + geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept = 30)

ggplot(dat.integrated@meta.data,aes(x=active.ident,y=log10(nCount_RNA))) + geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept = log10(900))

ggplot(subset(dat.integrated@meta.data, percent_mito < 30 & nCount_RNA > 900), 
       aes(x=active.ident,y=percent_ribo)) + geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept = 5)

ggplot(dat.integrated@meta.data,aes(x=active.ident,y=log10(nFeature_RNA))) + geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept = log10(900))



FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "seurat_clusters") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & poor_quality_cluster==F & dat.integrated@meta.data$percent_mito < 30) + 
  geom_quantile(quantiles = c(0.01, 0.05, 0.25 ,0.5, 0.75, 0.95, 0.99), method = "rqss", lambda = 0.1) + 
  geom_vline(xintercept = 900)

FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "seurat_clusters") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ nCount_RNA > 900 & dat.integrated@meta.data$percent_mito < 30) + 
  geom_quantile(quantiles = c(0.01, 0.05, 0.25 ,0.5, 0.75, 0.95, 0.99), method = "rqss", lambda = 0.1) + 
  geom_vline(xintercept = 900)

FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ dat.integrated@meta.data$seurat_clusters)

FeatureScatter(dat.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ dat.integrated@meta.data$active.ident)

DimPlot(dat.integrated.filtered, label = TRUE)

remaining_cells_per_cell_type <- table(dat.integrated@meta.data$active.ident,
      (dat.integrated@meta.data$too_high_mito==F & dat.integrated@meta.data$too_low_umi==F))

remaining_cells_per_cluster <- table(dat.integrated@meta.data$seurat_clusters,
                                       (dat.integrated@meta.data$too_high_mito==F & dat.integrated@meta.data$too_low_umi==F))

remaining_cells_per_cell_type
remaining_cells_per_cluster

# Preprocess again, without these cells ----------------------------------------

DefaultAssay(dat.integrated.filtered) <- "RNA"

dat.integrated.filtered.raw <- GetAssayData(dat.integrated.filtered, slot = "counts")
dat.integrated.filtered.raw <- CreateSeuratObject(dat.integrated.filtered.raw, 
                                                  meta.data = dat.integrated.filtered@meta.data)
dat.integrated.filtered.split <- SplitObject(dat.integrated.filtered.raw, split.by = "project.name")

rm(dat.integrated.filtered)
gc()

dat.integrated.filtered.split <- lapply(dat.integrated.filtered.split,
                                        function(x){
                                          SCTransform(x, vst.flavor = "v2", verbose = T) |> 
                                            RunPCA(npcs = 30, verbose = T)
                                          }
                                        )

features <- SelectIntegrationFeatures(object.list = dat.integrated.filtered.split, nfeatures = 3000)
 
dat.integrated.filtered.split <- PrepSCTIntegration(object.list = dat.integrated.filtered.split, 
                                                    anchor.features = features, verbose = T)

anchors <- FindIntegrationAnchors(object.list = dat.integrated.filtered.split,
                                  normalization.method = "SCT",
                                  anchor.features = features, verbose = T)
saveRDS(anchors,file="/data/proj/um_ss/Investigations/seurat/results/anchors.reprocessing.rda")

rm(dat.integrated.filtered.split)
rm(dat.integrated.filtered.raw)
gc()

dat.integrated.filtered.reprocecessed <- IntegrateData(anchorset = anchors, 
                                                       normalization.method = "SCT",verbose = T)
dat.integrated.filtered.reprocecessed <- RunPCA(dat.integrated.filtered.reprocecessed, 
                                                verbose = T)
dat.integrated.filtered.reprocecessed <- RunUMAP(dat.integrated.filtered.reprocecessed, 
                                                 reduction = "pca", dims = 1:30, verbose = T)
dat.integrated.filtered.reprocecessed <- FindNeighbors(dat.integrated.filtered.reprocecessed, 
                                                       reduction = "pca", dims = 1:30, verbose = T)
dat.integrated.filtered.reprocecessed <- FindClusters(dat.integrated.filtered.reprocecessed, 
                                                      resolution = 0.3, verbose = T)

#saveRDS(dat.integrated.filtered.reprocecessed,
#        file="/data/proj/um_ss/Investigations/seurat/results/dat.integrated.filtered.reprocecessed.rda")

## Eploration of filtered and reprocessed dataset (exploration) -----------------

dat.integrated.filtered.reprocecessed <- readRDS(
  "~/proj/um_ss/Investigations/seurat/results/dat.integrated.filtered.reprocecessed.rda")

pdf(paste0(outdir,"reprocessed.umap.cell_types.pdf"),width = 20,height = 14)
DimPlot(dat.integrated.filtered.reprocecessed, label = TRUE,group.by = "active.ident")
dev.off()

table(dat.integrated.filtered.reprocecessed$active.ident,
      dat.integrated.filtered.reprocecessed$seurat_clusters)

VlnPlot(object = dat.integrated.filtered.reprocecessed, 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito"))
VlnPlot(object = dat.integrated.filtered.reprocecessed, 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito"),
        group.by = "active.ident", pt.size=0) + xlab(NULL)

FeatureScatter(dat.integrated.filtered.reprocecessed, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ dat.integrated.filtered.reprocecessed@meta.data$seurat_clusters)

FeatureScatter(dat.integrated.filtered.reprocecessed, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_x_log10() + 
  scale_y_log10() + 
  facet_wrap( ~ dat.integrated.filtered.reprocecessed@meta.data$active.ident)

FeatureScatter(dat.integrated.filtered.reprocecessed, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_x_log10() + 
  scale_y_log10()

FeatureScatter(dat.integrated.filtered.reprocecessed, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "active.ident") + 
  scale_x_log10() + 
  scale_y_log10()

# DoubletFinder ----------------------------------------------------------------

dat.integrated.filtered.reprocecessed <- readRDS(
  "/data/proj/um_ss/Investigations/seurat/results/dat.integrated.filtered.reprocecessed.rda")

dat.integrated.filtered.reprocecessed_SCT <- dat.integrated.filtered.reprocecessed
DefaultAssay(dat.integrated.filtered.reprocecessed_SCT) <- "SCT"

dat.integrated.filtered.reprocecessed_SCT.split <- SplitObject(dat.integrated.filtered.reprocecessed_SCT, 
                                                               split.by = "Sample.ID")

# elbowdata <- function (object, ndims = 50, reduction = "pca") 
# {
#   data.use <- Stdev(object = object, reduction = reduction)
#   if (length(x = data.use) == 0) {
#     stop(paste("No standard deviation info stored for", reduction))
#   }
#   if (ndims > length(x = data.use)) {
#     warning("The object only has information for ", length(x = data.use), 
#             " reductions")
#     ndims <- length(x = data.use)
#   }
#   stdev <- "Standard Deviation"
#   data = data.frame(dims = 1:ndims, stdev = data.use[1:ndims])
#   return(data)
# }
# 
# # Preprocessing specifically for DuplicateFinder
# 
# elbowstrength <- function(ebd){
#   # Based on: https://community.ibm.com/community/user/datascience/blogs/moloy-de1/2020/07/02/points-to-ponder
#   ebd$delta_1 <- c(NA,diff(ebd$stdev))
#   ebd$delta_2 <- c(NA,diff(ebd$delta_1))
#   ebd$strength <- NA
#   for (k in 1:nrow(ebd)){
#     ebd$strength[k] <- ebd$delta_2[k+1] - ebd$delta_1[k+1]
#   }
#   return(ebd)
# }
# 
# plot_elbowstrength <- function(ebd){
#   g <- ggplot(ebd, aes(x = dims, y = stdev)) + geom_point() + 
#     geom_point(data=ebd,aes(x=dims,y=strength),colour="blue")
#   return(g)
# }

dat.duplicatefinder <- list()
for (nm in names(dat.integrated.filtered.reprocecessed_SCT.split)){
  dat.duplicatefinder[[nm]] <- dat.integrated.filtered.reprocecessed_SCT.split[[nm]]
  dat.duplicatefinder[[nm]] <- CreateSeuratObject(GetAssayData(dat.duplicatefinder[[nm]], 
                                slot = "counts"), meta.data = dat.duplicatefinder[[nm]]@meta.data)
  dat.duplicatefinder[[nm]] <- SCTransform(dat.duplicatefinder[[nm]])
  dat.duplicatefinder[[nm]] <- RunPCA(dat.duplicatefinder[[nm]])
}

# ndim_min <- 10
# ebd <- list()
# pca_ndim <- list()
# for (nm in names(dat.duplicatefinder)){
#   ebd[[nm]] <- elbowdata(dat.duplicatefinder[[nm]])
#   ebd[[nm]] <- elbowstrength(ebd[[nm]])
#   
#   pca_ndim[[nm]] <- ebd[[nm]][which(ebd[[nm]]$strength == 
#                       max(ebd[[nm]][ebd[[nm]]$dims >= ndim_min,]$strength, 
#                           na.rm = T)),]$dims
# }
# 
# for (nm in names(ebd)){
#   print(nm)
#   print(pca_ndim[[nm]])
#   plot(ebd[[nm]][,c("dims","stdev")])
#   abline(v = pca_ndim[[nm]])
#   readline()
# }

# SampleID_1_11june18: s
# SampleID_2_11june18: s
# SampleID_3_11june18: 16
# SampleID_4_11june18: 18
# SampleID_6_11june18: 18
# SampleID_7_11june18: 16
# SampleID_8_11june18: 16
# SampleID_9_11june18: 16
# 1A: s
# 1B: s
# 2A: s
# 2B: s
# 6A: s
# 6B: s
# 7A: 13
# 7B: 13
# A2-GEX: 16
# B3-GEX: 17
# C5-GEX: s
# C7-GEX: 16
# C8-GEX: 14
# D2-GEX: 16
# D4-GEX: s
# D6-GEX: 14

# Conclusion: the automatically determined inflection point can somethimes be debated,
# and DublicateFinder seems to benefint from more PCs (clearer peaks in bcmvn chart, as well
# as better results for one investigated case when using more PCs, with regards to nFeature_RNA violing plots
# comparing predicted duplicates to singlets). 15 is close to some average of number of typical
# inflection points for each sample.

ndim <- 15 # Use for all (see also https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/91 for further justification)

for (nm in names(dat.duplicatefinder)){
  dat.duplicatefinder[[nm]] <- RunUMAP(dat.duplicatefinder[[nm]], dims = 1:ndim)
}

saveRDS(dat.duplicatefinder,"/data/proj/um_ss/Investigations/seurat/results/dat.duplicatefinder.rda")

#dat.integrated <- readRDS("~/proj/um_ss/Investigations/seurat/results/dat.integrated.rda")
#dat.integrated <- readRDS("/data/proj/um_ss/Investigations/seurat/results/dat.integrated.rda")

# cells_recovered_per_sample <- table(dat.integrated@meta.data$Sample.ID)
# 
# assign_multiplet_rate <- function(cells_recovered_sample){
#   #CG000183_ChromiumSingleCell3__v3_UG_Rev_C.pdf
#   
#   mutliplet_rates_10x <- data.frame(
#     multiplet_rate = c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6)/100,
#     cells_loaded = c(800,1600,3200,4800,6400,8000,9600,11200,12800,14400,16000),
#     cells_recovered = c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000
#     )
#   )
#   
#   idx_closest <- which.min(abs(mutliplet_rates_10x$cells_recovered - cells_recovered_sample))
#   mutliplet_rate <- mutliplet_rates_10x[idx_closest,]$multiplet_rate 
#   
#   return(mutliplet_rate)
# }
# 
# sample_multiplet_rates <- sapply(cells_recovered_per_sample,assign_multiplet_rate)
# saveRDS(sample_multiplet_rates,"/data/proj/um_ss/Investigations/seurat/results/sample_multiplet_rates.rda")
sample_multiplet_rates <- readRDS("/data/proj/um_ss/Investigations/seurat/results/sample_multiplet_rates.rda")

outdir <- "/data/proj/um_ss/Investigations/seurat/results/DoubletFinder"
dir.create(outdir,recursive = T,showWarnings = F)

dat.duplicatefinder.res <- list()
sweep.res.list <- list()
sweep.stats <- list()
pK <- list()
nExp_poi <- list()
for (nm in setdiff(names(dat.duplicatefinder),names(dat.duplicatefinder.res))){
  print(nm)
  x <- dat.duplicatefinder[[nm]]

  sweep.res.list[[nm]] <- paramSweep_v3(x, PCs = 1:ndim, sct = T,num.cores = 4)
  sweep.stats[[nm]] <- summarizeSweep(sweep.res.list[[nm]], GT = F)
  bcmvn_x <- find.pK(sweep.stats[[nm]])
  
  bcmvn_x$pK <- as.numeric(as.character(bcmvn_x$pK))
  pK[[nm]] <- bcmvn_x$pK[which.max(bcmvn_x$BCmetric)]
  
  multiplet_rate <- sample_multiplet_rates[[nm]]
  nExp_poi[[nm]] <- round(multiplet_rate*nrow(x@meta.data))
  
  x <- doubletFinder_v3(x, PCs = 1:ndim, pN = 0.25, pK = pK[[nm]],
                        nExp = nExp_poi[[nm]], reuse.pANN = F, sct = T)
  print(table(x@meta.data[,grep("DF.classifications_",colnames(x@meta.data))]))

  #annotations <- x@meta.data$seurat_clusters
  #homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  #nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  #x <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
  
  dat.duplicatefinder.res[[nm]] <- x
}

saveRDS(dat.duplicatefinder.res,"/data/proj/um_ss/Investigations/seurat/results/dat.duplicatefinder.res.rda")

# g <- ggplot(bcmvn_x,aes(x=pK,y=BCmetric)) + geom_bar(stat="identity")
# ggsave(file=paste0(outdir,"/",nm,".bcmvn.pdf"))
# print(g)
# dev.off()

# Inspect DoubletFinder results ------------------------------------------------
dat.duplicatefinder.res <- readRDS("~/nimbus/data/proj/um_ss/Investigations/seurat/results/dat.duplicatefinder.res.rda")

proportion_doublets <- lapply(dat.duplicatefinder.res,function(x){
  DF.name = colnames(x@meta.data)[
    grepl("DF.classification", colnames(x@meta.data))]
  res <- table(x@meta.data[,DF.name])
  print(DF.name)
  return(res)
})

lapply(proportion_doublets,function(x){
  x[1]/sum(x)
})

outdir <- "/data/proj/um_ss/Investigations/seurat/results/DoubletFinder"
dir.create(outdir, recursive = T)

for (nm in names(dat.duplicatefinder.res)){
  DF.name = colnames(dat.duplicatefinder.res[[nm]]@meta.data)[
    grepl("DF.classification", colnames(dat.duplicatefinder.res[[nm]]@meta.data))]
  
  pdf(file = paste0(outdir,"/",nm,".pdf"),width = 15,height = 7)
  print(cowplot::plot_grid(
    ncol = 2, 
    DimPlot(dat.duplicatefinder.res[[nm]], group.by = "orig.ident") + NoAxes(),
    DimPlot(dat.duplicatefinder.res[[nm]], group.by = DF.name) + NoAxes()))
  dev.off()
  
  pdf(file = paste0(outdir,"/",nm,".nFeature_RNA.pdf"),width = 7,height = 7)
  print(VlnPlot(dat.duplicatefinder.res[[nm]], features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1))
  dev.off()
  
  pdf(file = paste0(outdir,"/",nm,".nCount_RNA.pdf"),width = 7,height = 7)
  print(VlnPlot(dat.duplicatefinder.res[[nm]], features = "nCount_RNA", group.by = DF.name, pt.size = 0.1))
  dev.off()
}

# Investigation of possible suboptimal DoubletFinder results (exploration) -----

# Possibly suboptimal sample results:
# B3-GEX, C5-GEX, C7-GEX, D4-GEX, D6-GEX, SampleID_8_11june18, SampleID_9_11june18

nms_rerun <- c("C5-GEX", "C7-GEX", "D6-GEX", "SampleID_8_11june18")
all(nms_rerun %in% names(dat.duplicatefinder))

sweep.res.list <- list()
sweep.stats <- list()
for (nm in nms_rerun){
  sweep.res.list[[nm]] <- paramSweep_v3(dat.duplicatefinder[[nm]], PCs = 1:ndim, sct = T)
  sweep.stats[[nm]] <- summarizeSweep(sweep.res.list[[nm]], GT = F)
}

reselected_pK <- setNames(c(0.03),
                          c("B3-GEX"))

nm <- names(sweep.stats)[1]

x <- dat.duplicatefinder[[nm]]

#sweep.res.list.redo <- paramSweep_v3(x, PCs = 1:15, sct = T)
#sweep.stats.redo <- summarizeSweep(sweep.res.list.redo, GT = F)

#bcmvn_x <- find.pK(sweep.stats[[nm]])
bcmvn_x <- find.pK(sweep.stats.redo)
bcmvn_x$pK <- as.numeric(as.character(bcmvn_x$pK))
bcmvn_x$pK[which.max(bcmvn_x$BCmetric)]

bcmvn_x <- bcmvn_x[setdiff(1:nrow(bcmvn_x),which.max(bcmvn_x$BCmetric)),]
pK <- bcmvn_x$pK[which.max(bcmvn_x$BCmetric)]

multiplet_rate <- sample_multiplet_rates[[nm]]
nExp_poi <- round(multiplet_rate*nrow(x@meta.data))

x <- doubletFinder_v3(x, PCs = 1:15, pN = 0.25, pK = 0.13,
                      nExp = nExp_poi, reuse.pANN = F, sct = T)
print(table(x@meta.data[,grep("DF.classifications_",colnames(x@meta.data))]))

grep("DF.classifications_",colnames(x@meta.data),value = T)

DF.name <- paste0("DF.classifications_","0.25_",0.03,"_",nExp_poi)
VlnPlot(x, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

VlnPlot(x, features = "nCount_RNA", group.by = DF.name, pt.size = 0.1)

# Filtering and re-clustering after DoubletFinder ------------------------------
dat.duplicatefinder.res <- readRDS("/data/proj/um_ss/Investigations/seurat/results/dat.duplicatefinder.res.rda")

dat.integrated.filtered.reprocecessed <- readRDS(
  "/data/proj/um_ss/Investigations/seurat/results/dat.integrated.filtered.reprocecessed.rda")

DefaultAssay(dat.integrated.filtered.reprocecessed) <- "SCT"


cells_remove <- list()
for (nm in names(dat.duplicatefinder.res)){
  DF.name = colnames(dat.duplicatefinder.res[[nm]]@meta.data)[
    grepl("DF.classification", colnames(dat.duplicatefinder.res[[nm]]@meta.data))]
  
  cells_remove[[nm]] <- rownames(dat.duplicatefinder.res[[nm]]@meta.data)[
    dat.duplicatefinder.res[[nm]]@meta.data[,DF.name] == "Doublet"]
}
cells_remove <- as.character(unlist(cells_remove))

dat.integrated.filtered.reprocecessed <- subset(dat.integrated.filtered.reprocecessed, 
       cells = setdiff(colnames(dat.integrated.filtered.reprocecessed),cells_remove))

dat.integrated.doubletfiltered <- SplitObject(dat.integrated.filtered.reprocecessed, 
                                                               split.by = "project.name")

dat.integrated.doubletfiltered.split <- list()
for (nm in names(dat.integrated.doubletfiltered)){
  DefaultAssay(dat.integrated.doubletfiltered[[nm]]) <- "RNA"
  dat.integrated.doubletfiltered.split[[nm]] <- CreateSeuratObject(
    counts = GetAssayData(dat.integrated.doubletfiltered[[nm]], slot = "counts"),
    meta.data = dat.integrated.doubletfiltered[[nm]]@meta.data)
}
rm(tmp)
rm(dat.integrated.doubletfiltered)
rm(dat.duplicatefinder.res)
rm(dat.integrated.filtered.reprocecessed)
gc()

dat.integrated.doubletfiltered.split <- lapply(dat.integrated.doubletfiltered.split,
                                        function(x){
                                          SCTransform(x, vst.flavor = "v2", verbose = T) |> 
                                            RunPCA(npcs = 50, verbose = T)
                                        }
)

features <- SelectIntegrationFeatures(object.list = dat.integrated.doubletfiltered.split, nfeatures = 3000)

dat.integrated.doubletfiltered.split <- PrepSCTIntegration(
  object.list = dat.integrated.doubletfiltered.split, 
  anchor.features = features, verbose = T)

saveRDS(dat.integrated.doubletfiltered.split,file="/data/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.split.rda")

anchors <- FindIntegrationAnchors(object.list = dat.integrated.doubletfiltered.split,
                                  normalization.method = "SCT",
                                  anchor.features = features, verbose = T, dims = 1:50, n.trees = 100)

saveRDS(anchors,
        file=paste0("/data/proj/um_ss/Investigations/seurat/results/anchors.doubletfiltered.rda"))

rm(dat.integrated.filtered.split)
rm(dat.integrated.filtered.raw)
gc()



Tool(object = dat.integrated.doubletfiltered.reprocecessed, slot = "Integration")@sample.tree

dat.integrated.doubletfiltered.reprocecessed <- IntegrateData(
  anchorset = anchors, normalization.method = "SCT",verbose = T, dims = 1:50, 
  sample.tree = matrix(c(-3, 1, -2, -1),ncol=2))
dat.integrated.doubletfiltered.reprocecessed <- RunPCA(
  dat.integrated.doubletfiltered.reprocecessed, verbose = T,npcs = 50)
dat.integrated.doubletfiltered.reprocecessed <- RunUMAP(
  dat.integrated.doubletfiltered.reprocecessed, reduction = "pca", dims = 1:50, 
  verbose = T)
dat.integrated.doubletfiltered.reprocecessed <- FindNeighbors(
  dat.integrated.doubletfiltered.reprocecessed, reduction = "pca", dims = 1:50, n.trees = 100,
  verbose = T)
dat.integrated.doubletfiltered.reprocecessed <- FindClusters(
  dat.integrated.doubletfiltered.reprocecessed, resolution = 0.8, verbose = T)

#saveRDS(dat.integrated.doubletfiltered.reprocecessed,
#        file="/data/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.rda")
saveRDS(dat.integrated.doubletfiltered.reprocecessed,
        file="/data/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.v2.rda")

outdir <- "/data/proj/um_ss/Investigations/seurat/results/"

pdf(file=paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.v2.annotated.dimplot.pdf"),
    width = 11,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = F)
dev.off()

pdf(file=paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.v2.annotated.dimplot.by_project.pdf"),
    width = 35,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = F, split.by = "project.name")
dev.off()

pdf(file=paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.v2.annotated.dimplot.tcr.pdf"),
    width = 11,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = F, group.by = "cdr3") + 
  theme(legend.position="none") + ggtitle(element_blank())
dev.off()

dat.integrated.doubletfiltered.reprocecessed <- RunTSNE(
  dat.integrated.doubletfiltered.reprocecessed, dims = 1:50)

saveRDS(dat.integrated.doubletfiltered.reprocecessed,
        file="/data/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.v2.rda")

pdf(file=paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.v2.annotated.dimplot.tsne.pdf"),
    width = 11,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = F, reduction = "tsne")
dev.off()

pdf(file=paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.v2.annotated.dimplot.tsne.by_project.pdf"),
    width = 35,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = F, split.by = "project.name", reduction = "tsne")
dev.off()

pdf(file=paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.v2.annotated.dimplot.tsne.tcr.pdf"),
    width = 11,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = F, group.by = "cdr3", reduction = "tsne") + 
  theme(legend.position="none") + ggtitle(element_blank())
dev.off()


# Visualization ----------------------------------------------------------------

#dat.integrated.doubletfiltered.reprocecessed <- readRDS("~/nimbus/data/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.rda")
#dat.integrated.doubletfiltered.reprocecessed <- readRDS("~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.rda")
dat.integrated.doubletfiltered.reprocecessed <- readRDS("~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.v2.rda")

dat.integrated.doubletfiltered.reprocecessed.old <- dat.integrated.doubletfiltered.reprocecessed
# 
# dat.integrated.doubletfiltered.reprocecessed <- RunPCA(
#   dat.integrated.doubletfiltered.reprocecessed, verbose = T,npcs = 50)
# 
# dat.integrated.doubletfiltered.reprocecessed <- RunUMAP(
#   dat.integrated.doubletfiltered.reprocecessed, reduction = "pca", dims = 1:50, 
#   verbose = T)
# 
DimPlot(dat.integrated.doubletfiltered.reprocecessed, group.by = "active.ident", reduction = "tsne")
DimPlot(dat.integrated.doubletfiltered.reprocecessed, group.by = "seurat_clusters", reduction = "tsne")

DimPlot(dat.integrated.doubletfiltered.reprocecessed, group.by = "active.ident", reduction = "tsne", split.by = "project.name")

dat.integrated.doubletfiltered.reprocecessed.tst <- subset(dat.integrated.doubletfiltered.reprocecessed, active.ident != "Melanocytic/Lymphocyte doublets, poor quality")

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_1.pdf"),width = 33,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, group.by = "active.ident", 
        split.by = "project.name", reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_2.pdf"),width = 33,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.tst, group.by = "active.ident", 
        split.by = "project.name", reduction = "tsne")
dev.off()

dat.integrated.doubletfiltered.reprocecessed.tst_2 <- subset(dat.integrated.doubletfiltered.reprocecessed, active.ident != "Melanocytic (poor quality)")
pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_3.pdf"),width = 33,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.tst_2, group.by = "active.ident", 
        split.by = "project.name", reduction = "tsne")
dev.off()

dat.integrated.doubletfiltered.reprocecessed.tst_4 <- subset(dat.integrated.doubletfiltered.reprocecessed, active.ident != "Melanocytic (MITF high)")
pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_4.pdf"),width = 33,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.tst_4, group.by = "active.ident", 
        split.by = "project.name", reduction = "tsne")
dev.off()


pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_1_1.pdf"),width = 15,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, group.by = "active.ident", 
        reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_4_1.pdf"),width = 15,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.tst_4, group.by = "active.ident", 
        reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_2_1.pdf"),width = 10,height = 8)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, group.by = "seurat_clusters", 
        reduction = "tsne")
dev.off()


dat.integrated.doubletfiltered.reprocecessed_RNA <- dat.integrated.doubletfiltered.reprocecessed
DefaultAssay(dat.integrated.doubletfiltered.reprocecessed_RNA) <- "RNA"

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.seurat_clusters.nFeature_RNA.pdf"),width = 20,height = 8)
ggplot(dat.integrated.doubletfiltered.reprocecessed_RNA@meta.data,aes(x=seurat_clusters,y=log10(nFeature_RNA))) + 
  geom_violin() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.seurat_clusters.percent_mito.pdf"),width = 20,height = 8)
ggplot(dat.integrated.doubletfiltered.reprocecessed_RNA@meta.data,aes(x=seurat_clusters,y=percent_mito)) + 
  geom_violin() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.seurat_clusters.nCount_RNA.pdf"),width = 20,height = 8)
ggplot(dat.integrated.doubletfiltered.reprocecessed_RNA@meta.data,aes(x=seurat_clusters,y=log10(nCount_RNA))) + 
  geom_violin() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

m.to_plot <- c(m.to_plot,"CCR4","CCR6")

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.seurat_clusters.dotplot.pdf"),width = 33,height = 12)
DotPlot(dat.integrated.doubletfiltered.reprocecessed_RNA, features = m.to_plot, cols = c("blue", "red"), 
        dot.scale = 8) +
  RotatedAxis()
dev.off()

tab <- table(dat.integrated.doubletfiltered.reprocecessed@meta.data$seurat_clusters,
             dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident)

tab[tab==0] <- NA

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.seurat_clusters.pheatmap.pdf"),
    width = 10,height = 12)
pheatmap(tab,border_color = NA,cluster_rows = F,na_col = "#FFFFFF",cluster_cols = F)
dev.off()

dat.integrated.doubletfiltered.reprocecessed.bak <- dat.integrated.doubletfiltered.reprocecessed

dat.integrated.doubletfiltered.reprocecessed <- AddMetaData(dat.integrated.doubletfiltered.reprocecessed,
                                                            metadata = setNames(dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident,rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)),
                                                            col.name = "active.ident.old")


# 1: Melanocytic (I)
# 5: Melanocytic (II)
# 14: Melanocytic (III)
# 21: Melanocytic (IV)
# 24: Melanocytic (V)
# 25: Melanocytic (VI)
# 29: Melanocytic (VII)
# 32: Melanocytic (VIII)
# 33: Melanocytic (IX)

# 2: CD4 T cells (I)
# 3: CD4 T cells (II)
# 6: CD4 T cells (III)
# 13: CD4 T cells (IV)
# 30: CD4 T cells (V)

# 0: CD8 T cells (I)
# 4: CD8 T cells (II)
# 7: CD8 T cells (III)
# 9: CD8 T cells (IV)
# 10: CD8 T cells (V)
# 12: CD8 T cells (VI)
# 16: CD8 T cells (VII)
# 19: CD8 T cells (VIII)
# 20: CD8 T cells (IX)
# 22: CD8 T cells (X)
# 23: CD8 T cells (XI)
# 24: CD8 T cells (XII)
# 26: CD8 T cells (XIII)
# 28: CD8 T cells (XIX)

# 8: NK cells (I)
# 17: NK cells (II)
# 35: NK cells (III)

# 11: Monocytes / Macrophages / DCs
# 18: Monocytes / Macrophages / DCs
# 27: Monocytes / Macrophages / DCs

# 31: pDC

# 15: B cells

# 34: Erythrocytes

dat.integrated.doubletfiltered.reprocecessed <- RenameIdents(
  dat.integrated.doubletfiltered.reprocecessed,
  `1` = "Melanocytic (I)",
  `5` = "Melanocytic (II)",
  `14` = "Melanocytic (III)",
  `21` = "Melanocytic (IV)",
  `24` = "Melanocytic (V)",
  `25` = "Melanocytic (VI)",
  `29` = "Melanocytic (VII)",
  `32` = "Melanocytic (VIII)",
  `33` = "Melanocytic (IX)",
  `2` = "CD4 T cells (I)",
  `3` = "CD4 T cells (II)",
  `6` = "CD4 T cells (III)",
  `13` = "CD4 T cells (IV)",
  `30` = "CD4 T cells (V)",
  `0` = "CD8 T cells (I)",
  `4` = "CD8 T cells (II)",
  `7` = "CD8 T cells (III)",
  `9` = "CD8 T cells (IV)",
  `10` = "CD8 T cells (V)",
  `12` = "CD8 T cells (VI)",
  `16` = "CD8 T cells (VII)",
  `19` = "CD8 T cells (VIII)",
  `20` = "CD8 T cells (IX)",
  `23` = "CD8 T cells (X)",
  `24` = "CD8 T cells (XI)",
  `28` = "CD8 T cells (XII)",
  `22` = "Gamma-delta T cells (mixed)",
  `26` = "Gamma-delta T cells (TRDV2, TRGV9)",
  `8` = "NK cells (I)",
  `17` = "NK cells (II)",
  `35` = "NK cells (III)",
  `11` = "CD14 Monocytes / DCs (I)",
  `18` = "CD14 Monocytes / DCs (II)",
  `27` = "CD16 Monocytes",
  `31` = "pDC",
  `15` = "B cells",
  `34` = "Erythrocytes"
  )

# Human dendritic cell subsets: an update (Matthew Collin, article)
# See also: https://ashpublications.org/blood/article/116/16/e74/27833/Nomenclature-of-monocytes-and-dendritic-cells-in:
# "In fact, it is very difficult at present to identify a single marker that can be used to clearly assign a cell to either the monocyte or the DC lineage."

# CD83 is a "marker of mature DCs". Might indicate that cluster 18 are DCs matured from
# monocytes in cluster 11. (but according to wikipedica monocytes also express CD83, 
# but as a soluble rather than membrane bound form)

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.pdf"),width = 33,height = 12)
DotPlot(dat.integrated.doubletfiltered.reprocecessed, 
        features = c(m.to_plot,"CD8B", 
                     grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                     grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                     "CLEC4C","THBD","NRP1","NRP2","IL12A","IL6","TLR2","TLR4","TLR7","TLR9","CD1C","CCR7",
                     "IRF4","KLF4","ZEB2","MAFB","ITGAX","IRF8","CD83","TRAV10","TRAJ18"),
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()
 
# grep("^TRAV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T), 
# grep("^TRBV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T)),

x <- dat.integrated.doubletfiltered.reprocecessed@assays$RNA@counts[rownames(dat.integrated.doubletfiltered.reprocecessed@assays$RNA@counts) %in% c(m.to_plot,"CD8B", 
                                                          grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                                                          grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                                                          "CLEC4C","THBD","NRP1","NRP2","IL12A","IL6","TLR2","TLR4","TLR7","TLR9","CD1C","CCR7",
                                                          "IRF4","KLF4","ZEB2","MAFB","ITGAX","IRF8","CD83","TRAV10","TRAJ18"),]
x <- as.matrix(x)
stopifnot(identical(colnames(dat.integrated.doubletfiltered.reprocecessed),names(dat.integrated.doubletfiltered.reprocecessed@active.ident)))
y <- aggregate(t(x), list(dat.integrated.doubletfiltered.reprocecessed@active.ident), mean)
rownames(y) <- y[,1]
y <- y[,2:ncol(y)]

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.markers_hclust.pdf"),width = 15,height = 10)
plot(hclust(dist(log2(y+1))))
dev.off()

dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident <- NULL

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_5.pdf"),width = 15,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = T, reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_5_1.pdf"),width = 35,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = T, reduction = "tsne",split.by = "project.name")
dev.off()


pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.G18-023.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "G18-023"), 
        features = c(m.to_plot,"CD8B", 
                     grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                     grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                     "CLEC4C","THBD","NRP1","NRP2","IL12A","IL6","TLR2","TLR4","TLR7","TLR9","CD1C","CCR7",
                     "IRF4","KLF4","ZEB2","MAFB","ITGAX","IRF8","CD83","TRAV10","TRAJ18"),
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.G18-049.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "G18-049"), 
        features = c(m.to_plot,"CD8B", 
                     grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                     grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                     "CLEC4C","THBD","NRP1","NRP2","IL12A","IL6","TLR2","TLR4","TLR7","TLR9","CD1C","CCR7",
                     "IRF4","KLF4","ZEB2","MAFB","ITGAX","IRF8","CD83","TRAV10","TRAJ18"),
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.GWA-JN-388.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "GWA-JN-388"), 
        features = c(m.to_plot,"CD8B", 
                     grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                     grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                     "CLEC4C","THBD","NRP1","NRP2","IL12A","IL6","TLR2","TLR4","TLR7","TLR9","CD1C","CCR7",
                     "IRF4","KLF4","ZEB2","MAFB","ITGAX","IRF8","CD83","TRAV10","TRAJ18"),
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

dat.integrated.doubletfiltered.reprocecessed <- AddMetaData(dat.integrated.doubletfiltered.reprocecessed,
            metadata = dat.integrated.doubletfiltered.reprocecessed@active.ident,
            col.name = "active.ident")

x <- dat.integrated.doubletfiltered.reprocecessed@assays$RNA@counts
dat.integrated.doubletfiltered.reprocecessed <- StashIdent(object = dat.integrated.doubletfiltered.reprocecessed, 
                                                           save.name = "old.ident")
dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident <- as.character(
  dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident)

idx_should_not_have_tcr <- dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident %in% c(
  "B cells",
  "CD14 Monocytes / DCs (I)",
  "CD14 Monocytes / DCs (II)",
  "CD16 Monocytes",
  "Erythrocytes",
  "Melanocytic (I)",
  "Melanocytic (II)",
  "Melanocytic (III)",
  "Melanocytic (IV)",
  "Melanocytic (V)",
  "Melanocytic (VI)",
  "Melanocytic (VII)",
  "Melanocytic (VIII)",
  "Melanocytic (IX)",
  "NK cells (I)",
  "NK cells (II)",
  "NK cells (III)",
  "pDC"
) & !is.na(dat.integrated.doubletfiltered.reprocecessed@meta.data$cdr3)

table(dat.integrated.doubletfiltered.reprocecessed@meta.data[idx_should_not_have_tcr
  ,]$project.name)

cells_unknown <- rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[idx_should_not_have_tcr]

# False melanocytes
cells <- rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[
  grepl("Melanocytic",dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident)]
#cells <- setdiff(cells,cells_unknown)
#"ITGB2","CCL5",
banned_genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD19",
                  "TIGIT","CTLA4","PDCD1","FCGR3A","PECAM1","PTPRC",
                  grep("^TRAV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRBV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRAJ",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T))
cells_unknown <- union(names(which(colSums(x[banned_genes,cells]) > 0)),cells_unknown)

cells_non_tumor_projects <- intersect(rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name %in% c("G18-023","G18-049")],cells)

cells_unknown <- union(cells_unknown,
                       names(which(colSums(x[c("MLANA","PMEL","TYR"),cells_non_tumor_projects])==0)))

y <- x[c(banned_genes,c("PMEL","MLANA","TYR","MITF")),intersect(cells_unknown,
             rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[
               dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "GWA-JN-388"])]
y <- as.matrix(y)
y[y > 0] <- 1
pheatmap(rbind(y[grep("^TR",rownames(y),invert = T),],colSums(y[grep("^TR",rownames(y)),])),show_colnames = F)

table(as.character(dat.integrated.doubletfiltered.reprocecessed@active.ident[cells_unknown]),
      dat.integrated.doubletfiltered.reprocecessed@meta.data[cells_unknown,]$project.name)
remaining <- setdiff(cells,cells_unknown)
table(as.character(dat.integrated.doubletfiltered.reprocecessed@active.ident[remaining]),
      dat.integrated.doubletfiltered.reprocecessed@meta.data[remaining,]$project.name)

# False erythrocytes
cells <- rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[
  grepl("Erythrocytes",dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident)]

# banned_genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD14","CD19",
#                   "HAVCR2","LAG3","TIGIT","CTLA4","PDCD1","FCGR3A","NCAM1",
#                   grep("^TRAV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
#                   grep("^TRBV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
#                   grep("^TRAJ",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
#                   grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
#                   grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T))

#banned_genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD19",
#                  "HAVCR2","LAG3","TIGIT","CTLA4","PDCD1",
banned_genes <- c(grep("^TRAV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRBV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRAJ",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T))


cells_unknown <- union(names(which(colSums(x[banned_genes,cells]) > 0)),cells_unknown)

table(as.character(dat.integrated.doubletfiltered.reprocecessed@active.ident[cells_unknown]),
      dat.integrated.doubletfiltered.reprocecessed@meta.data[cells_unknown,]$project.name)
remaining <- setdiff(cells,cells_unknown)
table(as.character(dat.integrated.doubletfiltered.reprocecessed@active.ident[remaining]),
      dat.integrated.doubletfiltered.reprocecessed@meta.data[remaining,]$project.name)

y <- x[c(banned_genes,c("HBB","HBA1","HBA2")),intersect(cells_unknown,
                                                                rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[
                                                                  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "GWA-JN-388"])]
y <- as.matrix(y)
y[y > 0] <- 1
pheatmap(rbind(y[grep("^TR",rownames(y),invert = T),],colSums(y[grep("^TR",rownames(y)),])),show_colnames = F)

length(cells)
length(cells_unknown)

# False monocytes / DCs / B cells
cells <- rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident %in% c(
    "B cells",
    "CD16 Monocytes",
    "CD14 Monocytes / DCs (I)",
    "CD14 Monocytes / DCs (II)",
    "pDC"
  )]

banned_genes <- c(grep("^TRAV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRBV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRAJ",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T))

#cells_unknown <- names(which(x[banned_genes,cells] > 0))
cells_unknown <- union(names(which(colSums(x[banned_genes,cells]) > 0)),cells_unknown)

table(as.character(dat.integrated.doubletfiltered.reprocecessed@active.ident[cells_unknown]),
      dat.integrated.doubletfiltered.reprocecessed@meta.data[cells_unknown,]$project.name)
remaining <- setdiff(cells,cells_unknown)
table(as.character(dat.integrated.doubletfiltered.reprocecessed@active.ident[remaining]),
      dat.integrated.doubletfiltered.reprocecessed@meta.data[remaining,]$project.name)

length(cells)
length(cells_unknown)

# #CD3G
# genes <- names(which(rowSums(x[,intersect(remaining,rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[
#   dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-023"
# ])])>0))
# 
# intersect(genes,m.to_plot)

# False B cells

cells <- rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident %in% c(
    "B cells")]
banned_genes <- c("CD3D","CD3E","CD3G","CD247")
cells_unknown <- union(names(which(colSums(x[banned_genes,cells]) > 0)),cells_unknown)

table(as.character(dat.integrated.doubletfiltered.reprocecessed@active.ident[cells_unknown]),
      dat.integrated.doubletfiltered.reprocecessed@meta.data[cells_unknown,]$project.name)
remaining <- setdiff(cells,cells_unknown)
table(as.character(dat.integrated.doubletfiltered.reprocecessed@active.ident[remaining]),
      dat.integrated.doubletfiltered.reprocecessed@meta.data[remaining,]$project.name)

length(cells)
length(cells_unknown)

active.ident <- as.character(dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident)
names(active.ident) <- rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)
active.ident[cells_unknown] <- "Unknown"

dat.integrated.doubletfiltered.reprocecessed <- AddMetaData(dat.integrated.doubletfiltered.reprocecessed,
            metadata = active.ident,col.name = "active.ident.updated")

dat.integrated.doubletfiltered.reprocecessed$active.ident.updated
Idents(dat.integrated.doubletfiltered.reprocecessed) <- dat.integrated.doubletfiltered.reprocecessed$active.ident.updated

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.G18-023.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "G18-023"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.G18-049.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "G18-049"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.GWA-JN-388.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "GWA-JN-388"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()



pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_6.pdf"),width = 15,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = T, reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_6_1.pdf"),width = 35,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = T, reduction = "tsne",split.by = "project.name")
dev.off()

identical(Idents(dat.integrated.doubletfiltered.reprocecessed),dat.integrated.doubletfiltered.reprocecessed@active.ident)

dat.melanocytes <- subset(dat.integrated.doubletfiltered.reprocecessed, 
                          idents = sort(unique(grep("Melanocytic",dat.integrated.doubletfiltered.reprocecessed@active.ident,value = T))))

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_7.pdf"),width = 10,height = 10)
DimPlot(dat.melanocytes, label = T, reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_7_1.pdf"),width = 10,height = 10)
DimPlot(dat.melanocytes, label = T, reduction = "tsne",split.by = "project.name")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.tmp_7.dotplot.G18-023.pdf"),width = 33,height = 12)
DotPlot(subset(dat.melanocytes, project.name == "G18-023"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.tmp_7.dotplot.G18-049.pdf"),width = 33,height = 12)
DotPlot(subset(dat.melanocytes, project.name == "G18-049"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.tmp_7.dotplot.GWA-JN-388.pdf"),width = 33,height = 12)
DotPlot(subset(dat.melanocytes, project.name == "GWA-JN-388"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()



dat.unknown <- subset(dat.integrated.doubletfiltered.reprocecessed, 
                          idents = sort(unique(grep("Unknown",dat.integrated.doubletfiltered.reprocecessed@active.ident,value = T))))

dat.unknown <- AddMetaData(dat.unknown, 
                           metadata = setNames(ifelse(grepl("Melanocytic",dat.unknown@meta.data$active.ident),"Melanocytic","Unknown"),
                                               rownames(dat.unknown@meta.data)),col.name = "is_melanocytic")

dat.unknown.melanocytic <- subset(dat.unknown,is_melanocytic == "Melanocytic")

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_8.pdf"),width = 10,height = 10)
DimPlot(dat.unknown.melanocytic, label = T, reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_8_1.pdf"),width = 33,height = 10)
DimPlot(dat.unknown.melanocytic, label = T, reduction = "tsne",split.by = "project.name")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.tmp_8.dotplot.G18-023.pdf"),width = 33,height = 12)
DotPlot(subset(dat.unknown.melanocytic, project.name == "G18-023"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.tmp_8.dotplot.G18-049.pdf"),width = 33,height = 12)
DotPlot(subset(dat.unknown.melanocytic, project.name == "G18-049"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.tmp_8.dotplot.GWA-JN-388.pdf"),width = 33,height = 12)
DotPlot(subset(dat.unknown.melanocytic, project.name == "GWA-JN-388"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.tmp_8.dotplot.pdf"),width = 33,height = 12)
DotPlot(dat.unknown.melanocytic, 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

dat.unknown.melanocytic_RNA <- dat.unknown.melanocytic
DefaultAssay(dat.unknown.melanocytic_RNA) <- "RNA"

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.unknown_melanocytic_CD8A.pdf"),width = 33,height = 12)
FeaturePlot(dat.unknown.melanocytic_RNA, reduction = "tsne",
            split.by = "project.name",features = "CD8A", cols = c("white","blue"),order = T)
dev.off()
pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.unknown_melanocytic_PTPCR.pdf"),width = 33,height = 12)
FeaturePlot(dat.unknown.melanocytic_RNA, reduction = "tsne",
            split.by = "project.name",features = "PTPRC", cols = c("white","blue"),order = T)
dev.off()
pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.unknown_melanocytic_PMEL.pdf"),width = 33,height = 12)
FeaturePlot(dat.unknown.melanocytic_RNA, reduction = "tsne",
            split.by = "project.name",features = "PMEL", cols = c("white","blue"),order = T)
dev.off()
pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.unknown_melanocytic_MLANA.pdf"),width = 33,height = 12)
FeaturePlot(dat.unknown.melanocytic_RNA, reduction = "tsne",
            split.by = "project.name",features = "MLANA", cols = c("white","blue"),order = T)
dev.off()


dat.unknown.melanocytic <- subset(dat.unknown,is_melanocytic == "Melanocytic")
dat.unknown.melanocytic.tils <- subset(dat.unknown.melanocytic,project.name == "G18-023")
dat.unknown.melanocytic.blood <- subset(dat.unknown.melanocytic,project.name == "G18-049")


banned_genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD19",
                  "TIGIT","CTLA4","PDCD1","FCGR3A","MS4A1","PECAM1","PTPRC",
                  grep("^TRAV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRBV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRAJ",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRDV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T),
                  grep("^TRGV",rownames(dat.integrated.doubletfiltered.reprocecessed),value=T))

m <- intersect(markers_plot,rownames(dat.unknown.melanocytic.tils@assays$RNA@counts))
x <- dat.unknown.melanocytic.tils@assays$RNA@counts
y <- x[m,]
y <- as.matrix(y)
y[y > 0] <- 1
y <- y[rowSums(y)>0,]
y <- y[c("CD3D","CD4","CD8A","CD8B"),]

cells_tils_cd8 <- names(which(y["CD3D",] > 0 & y["CD4",] > 0 & ! (y["CD8A",] > 0 | y["CD8B",] > 0 ) ))
cells_tils_cd4 <- names(which(y["CD3D",] > 0 & (y["CD8A",] > 0 | y["CD8B",] > 0 ) & ! (y["CD4",] > 0 )))

pheatmap(rbind(y[grep("^TR",rownames(y),invert = T),],colSums(y[grep("^TR",rownames(y)),])),
         show_colnames = F,treeheight_row = 0,treeheight_col = 0,fontsize_row = 5)

x <- dat.unknown.melanocytic.blood@assays$RNA@counts
y <- x[m,]
y <- as.matrix(y)
y[y > 0] <- 1
y <- y[rowSums(y)>0,]
y <- y[c("CD3D","CD4","CD8A","CD8B"),]
pheatmap(rbind(y[grep("^TR",rownames(y),invert = T),],colSums(y[grep("^TR",rownames(y)),])),
         show_colnames = F,treeheight_row = 0,treeheight_col = 0,fontsize_row = 5)

cells_blood_cd8 <- names(which(y["CD3D",] > 0 & y["CD4",] > 0 & ! (y["CD8A",] > 0 | y["CD8B",] > 0 ) ))
cells_blood_cd4 <- names(which(y["CD3D",] > 0 & (y["CD8A",] > 0 | y["CD8B",] > 0 ) & ! (y["CD4",] > 0 )))

dat.integrated.doubletfiltered.reprocecessed@meta.data[c(cells_tils_cd8,cells_blood_cd8),]$active.ident.updated <- "CD8 T cells (reclassified I)"
dat.integrated.doubletfiltered.reprocecessed@meta.data[c(cells_tils_cd4,cells_blood_cd4),]$active.ident.updated <- "CD4 T cells (reclassified I)"

dat.integrated.doubletfiltered.reprocecessed <- AddMetaData(
  dat.integrated.doubletfiltered.reprocecessed,
  metadata = dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident.updated,col.name = "active.ident.updated_2")

dat.integrated.doubletfiltered.reprocecessed$active.ident.updated_2
Idents(dat.integrated.doubletfiltered.reprocecessed) <- dat.integrated.doubletfiltered.reprocecessed$active.ident.updated_2

table(Idents(dat.integrated.doubletfiltered.reprocecessed))




library(ggplot2)
ggplot(dat.unknown.melanocytic.blood@meta.data,aes(x=active.ident,y=log10(nFeature_RNA))) + geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(dat.unknown.melanocytic.tils@meta.data,aes(x=active.ident,y=log10(nFeature_RNA))) + geom_violin() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



dat.integrated.doubletfiltered.reprocecessed <- AddMetaData(dat.integrated.doubletfiltered.reprocecessed,
                                                            metadata = setNames(ifelse(grepl("Unknown",dat.integrated.doubletfiltered.reprocecessed@active.ident),"Unknown","Retained"),rownames(dat.integrated.doubletfiltered.reprocecessed@meta.data)),
                                                            col.name = "is_unknown")


ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "GWA-JN-388" & 
    grepl("Melanocytic",dat.integrated.doubletfiltered.reprocecessed@meta.data$active.ident)
,],
       aes(x=active.ident,y=log10(nFeature_RNA), color = is_unknown)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#G18-023
#G18-049
ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-023",],
  aes(x=is_unknown,y=log10(nFeature_RNA), color = is_unknown)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-049",],
  aes(x=is_unknown,y=log10(nFeature_RNA), color = is_unknown)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "GWA-JN-388",],
  aes(x=is_unknown,y=log10(nFeature_RNA), color = is_unknown)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-023",],
  aes(x=is_unknown,y=nFeature_RNA, color = percent_mito)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-049",],
  aes(x=is_unknown,y=log10(nFeature_RNA), color = percent_mito)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "GWA-JN-388",],
  aes(x=is_unknown,y=log10(nFeature_RNA), color = percent_mito)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


table(is.na(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-023" & 
    dat.integrated.doubletfiltered.reprocecessed@meta.data$is_unknown == "Unknown"
  ,]$cdr3))

table(is.na(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-049" & 
    dat.integrated.doubletfiltered.reprocecessed@meta.data$is_unknown == "Unknown"
    ,]$cdr3))

table(is.na(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "GWA-JN-388" & 
    dat.integrated.doubletfiltered.reprocecessed@meta.data$is_unknown == "Unknown"
  ,]$cdr3))

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.tmp_9.G18-023.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "G18-023"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.tmp_9.G18-049.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "G18-049"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.renamed_clusters.dotplot.tmp_9.GWA-JN-388.pdf"),width = 33,height = 12)
DotPlot(subset(dat.integrated.doubletfiltered.reprocecessed, project.name == "GWA-JN-388"), 
        features = markers_plot,
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA") +
  RotatedAxis()
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_9.pdf"),width = 15,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = T, reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_9_1.pdf"),width = 35,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed, label = T, reduction = "tsne",split.by = "project.name")
dev.off()

dat.integrated.doubletfiltered.reprocecessed <- AddMetaData(dat.integrated.doubletfiltered.reprocecessed,
                                                            metadata = setNames(ifelse(grepl("reclassified",dat.integrated.doubletfiltered.reprocecessed@active.ident),"reclassified","same"),
                                                                                names(dat.integrated.doubletfiltered.reprocecessed@active.ident)),col.name = "is_reclassified")

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.tmp_9_2.pdf"),width = 35,height = 10)
DimPlot(subset(dat.integrated.doubletfiltered.reprocecessed, is_reclassified == "reclassified"), label = T, reduction = "tsne", 
        split.by = "project.name")
dev.off()

dat.integrated.doubletfiltered.reprocecessed.filtered <- subset(dat.integrated.doubletfiltered.reprocecessed,idents = setdiff(dat.integrated.doubletfiltered.reprocecessed@active.ident,"Unknown"))

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.pdf"),width = 15,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.filtered, label = T, reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered_split.pdf"),width = 35,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.filtered, label = T, reduction = "tsne",split.by = "project.name")
dev.off()



pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.t_cell_nFeatures.G18-023.pdf"),width = 10,height = 7)
ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-023" & 
    grepl("T cell",dat.integrated.doubletfiltered.reprocecessed@active.ident),],
  aes(x=active.ident.updated_2,y=nFeature_RNA)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.t_cell_nFeatures.G18-049.pdf"),width = 10,height = 7)
ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-049" & 
    grepl("T cell",dat.integrated.doubletfiltered.reprocecessed@active.ident),],
  aes(x=active.ident.updated_2,y=nFeature_RNA)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.t_cell_nFeatures.GWA-JN-388.pdf"),width = 10,height = 7)
ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "GWA-JN-388" & 
    grepl("T cell",dat.integrated.doubletfiltered.reprocecessed@active.ident),],
  aes(x=active.ident.updated_2,y=nFeature_RNA)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.t_cell_mito.G18-023.pdf"),width = 10,height = 7)
ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-023" & 
    grepl("T cell",dat.integrated.doubletfiltered.reprocecessed@active.ident),],
  aes(x=active.ident.updated_2,y=percent_mito)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.t_cell_mito.G18-049.pdf"),width = 10,height = 7)
ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "G18-049" & 
    grepl("T cell",dat.integrated.doubletfiltered.reprocecessed@active.ident),],
  aes(x=active.ident.updated_2,y=percent_mito)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.t_cell_mito.GWA-JN-388.pdf"),width = 10,height = 7)
ggplot(dat.integrated.doubletfiltered.reprocecessed@meta.data[
  dat.integrated.doubletfiltered.reprocecessed@meta.data$project.name == "GWA-JN-388" & 
    grepl("T cell",dat.integrated.doubletfiltered.reprocecessed@active.ident),],
  aes(x=active.ident.updated_2,y=percent_mito)) + geom_violin(draw_quantiles = 0.5) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# The reclassified T cells have significantly higher mitochondiral percentage than other T cells from
# the same datasets, and fewer discovered genes. Remove them.

dat.integrated.doubletfiltered.reprocecessed.filtered <- subset(dat.integrated.doubletfiltered.reprocecessed.filtered,
                                                                idents = setdiff(dat.integrated.doubletfiltered.reprocecessed.filtered@active.ident,
                                                                                 c("CD4 T cells (reclassified I)","CD8 T cells (reclassified I)")))

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.pdf"),width = 15,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.filtered, label = T, reduction = "tsne")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered_split.pdf"),width = 35,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.filtered, label = T, reduction = "tsne", split.by = "project.name")
dev.off()

pdf(file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered_split_samples.pdf"),width = 200,height = 10)
DimPlot(dat.integrated.doubletfiltered.reprocecessed.filtered, label = T, reduction = "tsne", split.by = "Sample.ID")
dev.off()

saveRDS(dat.integrated.doubletfiltered.reprocecessed.filtered,
        file = paste0(outdir,"dat.integrated.doubletfiltered.reprocecessed.filtered.rda"))

# Reintegrate after further filtering ------------------------------------------
dat.integrated.doubletfiltered.reprocecessed.filtered <- readRDS("/data/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.rda")


dat.split <- SplitObject(dat.integrated.doubletfiltered.reprocecessed.filtered, 
                                              split.by = "project.name")

for (nm in names(dat.split)){
  DefaultAssay(dat.split[[nm]]) <- "RNA"
  dat.split[[nm]] <- CreateSeuratObject(
    counts = GetAssayData(dat.split[[nm]], slot = "counts"),
    meta.data = dat.split[[nm]]@meta.data)
}

dat.split <- lapply(dat.split,
                    function(x){
                      SCTransform(x, vst.flavor = "v2", verbose = T) |> 
                        RunPCA(npcs = 50, verbose = T)
                      }
)

features <- SelectIntegrationFeatures(object.list = dat.split, nfeatures = 3000)

dat.split <- PrepSCTIntegration(
  object.list = dat.split, 
  anchor.features = features, verbose = T)

saveRDS(dat.split,file="/data/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.split.rda")

anchors <- FindIntegrationAnchors(object.list = dat.split,
                                  normalization.method = "SCT",
                                  anchor.features = features, verbose = T, dims = 1:50, n.trees = 100)

saveRDS(anchors,
        file=paste0("/data/proj/um_ss/Investigations/seurat/results/anchors.dat.integrated.doubletfiltered.reprocecessed.filtered.split.rda"))

dat.split <- IntegrateData(
  anchorset = anchors, normalization.method = "SCT",verbose = T, dims = 1:50, 
  sample.tree = matrix(c(-3, 1, -2, -1),ncol=2))
dat.split <- RunPCA(
  dat.split, verbose = T,npcs = 50)
dat.split <- RunUMAP(
  dat.split, reduction = "pca", dims = 1:50, 
  verbose = T)
dat.split <- FindNeighbors(
  dat.split, reduction = "pca", dims = 1:50, n.trees = 100,
  verbose = T)
dat.split <- FindClusters(
  dat.split, resolution = 0.8, verbose = T)

saveRDS(dat.split,file="/data/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")

