library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(WriteXLS)
library(scater)
library(djvdj)
library(scDblFinder)
library(Matrix)

vertical_xlabels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                     hjust=1))

outdir <- "/data/proj/um_ss/Investigations/pdx_human/results/v1/"
dir.create(outdir,recursive=T,showWarnings=F)

outdir <- "~/proj/um_ss/Investigations/pdx_human/results/v1/"

# Read cellbender data ---------------------------------------------------------
Read10X_h5.custom <- function (filename, use.names = TRUE, unique.features = TRUE)
{
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  blacklist <- c("cell.data","cell_sums","gene_means","gene_sum_sq","gene_vars","norm.data","scale.data")
  genomes <- names(x = infile)
  genomes <- setdiff(genomes,blacklist)
  output <- list()
  if (hdf5r::existsGroup(infile, "matrix")) {
    if (use.names) {
      feature_slot <- "features/name"
    }
    else {
      feature_slot <- "features/id"
    }
  }
  else {
    if (use.names) {
      feature_slot <- "gene_names"
    }
    else {
      feature_slot <- "genes"
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    features <- infile[[paste0(genome, "/", feature_slot)]][]
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[],
                               x = as.numeric(x = counts[]), dims = shp[], repr = "T")
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as.sparse(x = sparse.mat)
    if (infile$exists(name = paste0(genome, "/features"))) {
      types <- infile[[paste0(genome, "/features/feature_type")]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(X = types.unique, FUN = function(x) {
          return(sparse.mat[which(x = types == x), ])
        }, simplify = FALSE, USE.NAMES = TRUE)
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  }
  else {
    return(output)
  }
}

fnames <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/multi/PDX*_graft/outs/multi/count/cellbender_feature_bc_matrix/cellbender_feature_bc_matrix_filtered.h5")

snames <- basename(dirname(dirname(dirname(dirname(dirname(fnames))))))
snames <- gsub("_","",snames)

names(fnames) <- snames
fnames <- as.list(fnames)

seu <- list()
for (sname in names(fnames)){
  seu[[sname]] <- Read10X_h5.custom(fnames[[sname]])
  colnames(seu[[sname]]) <- paste0(sname,"_",colnames(seu[[sname]]))
}

seu <- lapply(seu,Seurat::CreateSeuratObject)
seu <- lapply(seu,Seurat::NormalizeData)

seu <- RunFastMNN(object.list = seu,
                  features = 4000, auto.merge = T, d = 150, k = 15)

seu <- RunUMAP(seu, reduction = "mnn", dims = 1:150)
#saveRDS(seu,file = paste0(outdir,"seu.rda"))

pdf(file=paste0(outdir,"markers_pre_doublet_filtering.pdf"),width = 40,height = 40)
FeaturePlot(object = seu, reduction = "umap", pt.size = .1,
            features = c("MLANA","CD3D","NCR1","CD4","CD8A","CD14","HBA1",
                         "JCHAIN","MS4A1","FOXP3","CD34","MKI67",
                         "CD1C","THBD","FCGR3A","FCGR3B","CLEC4C",
                         "IL3RA","ITGAX","SIRPA","CLEC10A","SERPINF1",
                         "LILRA4"))
dev.off()

pdf(file=paste0(outdir,"umap_samples.pre_doublet_filtering.pdf"),width = 10,height = 10)
DimPlot(seu, group.by="orig.ident")
dev.off()

# Add VDJ data -----------------------------------------------------------------
add_vdj_to_seurat <- function(lig.tumor, vdj_paths){
  library(djvdj)
  fnames <- vdj_paths
  vdj <- list()
  for (nm in names(fnames)) {
    vdj[[nm]] <- djvdj::import_vdj(input = NULL, vdj_dir = fnames[[nm]],
                                   filter_chains = T,
                                   define_clonotypes = "cdr3_gene",prefix = )
  }
  
  for (nm in names(fnames)) {
    rownames(vdj[[nm]]) <- paste0(nm,"_",rownames(vdj[[nm]]))
  }
  
  vdj <- do.call("rbind",vdj)
  rownames(vdj) <- str_split_fixed(rownames(vdj),"\\.",2)[,2]
  
  orig.rownames <- rownames(lig.tumor@meta.data)
  meta <- merge(lig.tumor@meta.data,vdj,by="row.names",all.x=T)
  rownames(meta) <- meta$Row.names
  meta$Row.names <- NULL
  meta <- meta[orig.rownames,]
  stopifnot(identical(rownames(meta),rownames(lig.tumor@meta.data)))
  lig.tumor@meta.data <- meta
  return(lig.tumor)
}

fnames <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/multi/PDX*_graft/outs/per_sample_outs/*/vdj_t/")
fnames <- fnames[gsub("_","",basename(dirname(fnames))) %in% unique(seu@meta.data$orig.ident)]
names(fnames) <- gsub("_","",basename(dirname(fnames)))
fnames <- fnames[names(fnames) %in% seu@meta.data$orig.ident]

seu <- add_vdj_to_seurat(seu, vdj_paths = fnames)

saveRDS(seu,file = paste0(outdir,"/seu.rda"))

# Doublets ---------------------------------------------------------------------
library(scDblFinder)

# The authors do not recommend further QC filtering than this before doublet detection
idx_keep <- seu@meta.data[["nCount_RNA"]] > 200

cells_keep <- rownames(seu@meta.data[idx_keep,])

seu <- subset(seu, subset = nCount_RNA > 200)

sce <- as.SingleCellExperiment(seu)

library(BiocParallel)
sce.standard <- scDblFinder(sce, samples = "orig.ident",
                            BPPARAM = MulticoreParam(20))

# Investigate outcome of doublet detection -------------------------------------
mat <- counts(sce)

stopifnot(identical(colnames(mat),rownames(seu@meta.data)))
idx <- (!is.na(seu@meta.data$cdr3) |
          colSums(mat[c("CD3D","CD4","CD8A","CD8B"),]>0) > 0) &
  colSums(mat[c("MLANA","PMEL","TYR","HBA1","HBA2","HBB"),]>0) > 0
idx <- idx | ((!is.na(seu@meta.data$cdr3) &
                 colSums(mat[c("CD14","CD19","MS4A1","JCHAIN"),]>0) > 0))
idx <- idx | ((colSums(mat[c("MLANA","PMEL","TYR"),]>0) > 0) &
                (colSums(mat[c("CD14","CD19","MS4A1","JCHAIN","NCR1"),]>0) > 0))
knownDoublets <- names(which(idx))
seu@meta.data$knownDoublets <- rownames(seu@meta.data) %in% knownDoublets

knownDoublets <- rownames(colData(sce)) %in% knownDoublets
colData(sce) <- cbind(colData(sce),knownDoublets)

table(colData(sce.standard)[["scDblFinder.class"]])

stopifnot(identical(rownames(colData(sce.standard)),
                    rownames(seu@meta.data)))

seu@meta.data$scDblFinder.standard.score <-
  colData(sce.standard)[["scDblFinder.score"]]
seu@meta.data$scDblFinder.standard.class <-
  colData(sce.standard)[["scDblFinder.class"]]

pdf(file=paste0(outdir,"scDblFinder.standard.score.pdf"),width = 10,height = 10)
FeaturePlot(seu, reduction = "umap",features = "scDblFinder.standard.score")
dev.off()

pdf(file=paste0(outdir,"scDblFinder.standard.class.pdf"),width = 10,height = 10)
DimPlot(seu, reduction = "umap",group.by = "scDblFinder.standard.class")
dev.off()

pdf(file=paste0(outdir,"knownDoublets.pdf"),width = 10,height = 10)
FeaturePlot(seu, reduction = "umap",features = "knownDoublets")
dev.off()

table(doublet_agreement=seu@meta.data$scDblFinder.standard.class)["doublet"] /
  nrow(seu@meta.data)

table(seu@meta.data$knownDoublets,
      seu@meta.data$scDblFinder.standard.class)

seu@meta.data$is_doublet <- seu@meta.data$scDblFinder.standard.class=="doublet" |
  seu@meta.data$knownDoublets


pdf(file=paste0(outdir,"is_doublet_incl_cluster.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "is_doublet")
dev.off()

# Add some QC parameters -------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

pdf(file=paste0(outdir,"percent.mt.pdf"),width = 10,height = 10)
FeaturePlot(object = seu, features = "percent.mt")
dev.off()

pdf(file=paste0(outdir,"percent.ribo.pdf"),width = 10,height = 10)
FeaturePlot(object = seu, features = "percent.ribo")
dev.off()

#https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html

sce <- as.SingleCellExperiment(seu)

mt.genes <- grep(pattern = "^MT-", x = rownames(seu@assays$RNA@counts), value = T)

sce <- addPerCellQC(sce, subsets = list(mito = mt.genes),
                    percent.top=100, use.altexps=F)

stopifnot(identical(rownames(colData(sce)),rownames(seu@meta.data)))
seu@meta.data$percent.top_100 <- colData(sce)[,"percent.top_100"]

# Plots to establish thresholds ------------------------------------------------
seu@meta.data$percent.ribo.40 <- seu@meta.data$percent.ribo > 40
pdf(file=paste0(outdir,"percent.ribo.40.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "percent.ribo.40")
dev.off()

seu@meta.data$percent.mt.30 <- seu@meta.data$percent.mt > 30
pdf(file=paste0(outdir,"percent.mt.30.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "percent.mt.30")
dev.off()

hist(seu@meta.data$percent_mito)
hist(seu@meta.data$percent_ribo)

HBA1 <- seu@assays$RNA@counts["HBA1",]
CD3D <- seu@assays$RNA@counts["CD3D",]
MLANA <- seu@assays$RNA@counts["MLANA",]
PMEL <- seu@assays$RNA@counts["PMEL",]
NCAM1 <- seu@assays$RNA@counts["NCAM1",]
stopifnot(identical(names(HBA1),rownames(seu@meta.data)))
stopifnot(identical(names(CD3D),rownames(seu@meta.data)))
stopifnot(identical(names(MLANA),rownames(seu@meta.data)))
seu@meta.data <- cbind(seu@meta.data, HBA1)
seu@meta.data <- cbind(seu@meta.data, CD3D)
seu@meta.data <- cbind(seu@meta.data, MLANA)
seu@meta.data <- cbind(seu@meta.data, PMEL)
seu@meta.data <- cbind(seu@meta.data, NCAM1)

### Distinguish erythrocytes
stopifnot(identical(colnames(seu@assays$RNA@data),colnames(seu@assays$RNA@counts)))
mat <- seu@assays$RNA@data[,seu@assays$RNA@counts["HBA1",] > 0]
mat <- mat[which(rowSums(mat) > 0),]
dim(mat)

library(pheatmap)
annot <- as.data.frame(t(mat[c("MLANA","PMEL","CD3D","NCAM1","HBA1"),]))
annot <- merge(annot, seu@meta.data[,c("nFeature_RNA","nCount_RNA","is_doublet","percent.top_100")],
               all.x=T, all.y=F, by="row.names")
annot$is_doublet <- ifelse(annot$is_doublet,1,0)
rownames(annot) <- annot$Row.names
annot$Row.names <- NULL
pdf(file=paste0(outdir,"pheatmap_HBA1.pdf"))
p <- pheatmap(cor(as.matrix(mat),method = "spearman"),
              show_colnames = F,show_rownames = F,annotation_col = annot)
dev.off()
clus <- cutree(p$tree_col,2)
table(clus)

# Inspection of the heatmap suggests that non of the HBA1 cells are non-doublets
# Blacklist all of them.

#cl <- which.max(c(sum(mat["MLANA",names(which(clus==1))]>0),
#                  sum(mat["MLANA",names(which(clus==2))]>0)))

#hba1_blacklist <- names(which(clus==cl))
hba1_blacklist <- names(clus)

seu@meta.data$hba1_blacklist <- F
seu@meta.data[hba1_blacklist,]$hba1_blacklist <- T
seu@meta.data[hba1_blacklist,]$is_doublet <- T

pdf(file=paste0(outdir,"HBA1_blacklist.pdf"))
ggplot(seu@meta.data[seu@meta.data$HBA1 > 0,],
       aes(x=nFeature_RNA,y=nCount_RNA, color = hba1_blacklist)) +
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)
dev.off()

#min_nGene_erythrocytes <- summary(annot[names(which(clus!=cl)),]$nFeature_RNA)[["1st Qu."]]
#min_nUMI_erythrocytes <- summary(annot[names(which(clus!=cl)),]$nCount_RNA)[["1st Qu."]]

### Distinguish plasma cells
stopifnot(identical(colnames(seu@assays$RNA@data),colnames(seu@assays$RNA@counts)))
mat <- seu@assays$RNA@data[,seu@assays$RNA@counts["JCHAIN",] > 0 &
                             seu@assays$RNA@counts["IL3RA",] == 0 &
                             seu@assays$RNA@counts["LILRA4",] == 0]
mat <- mat[which(rowSums(mat) > 0),]
dim(mat)

annot <- as.data.frame(t(mat[intersect(c("MLANA","PMEL","CD3D","NCAM1","JCHAIN"),
                                       rownames(mat)),]))
annot <- merge(annot, seu@meta.data[,c("nFeature_RNA","nCount_RNA","is_doublet",
                                       "percent.top_100")],
               all.x=T, all.y=F, by="row.names")
annot$is_doublet <- ifelse(annot$is_doublet,1,0)
rownames(annot) <- annot$Row.names
annot$Row.names <- NULL

pdf(file=paste0(outdir,"pheatmap_JCHAIN.pdf"))
p <- pheatmap(cor(as.matrix(mat),method = "spearman"),
              show_colnames = F,show_rownames = F,annotation_col = annot,
              cutree_cols = min(3,ncol(mat)),cutree_rows = min(3,ncol(mat)))
dev.off()
clus <- cutree(p$tree_col,min(3,ncol(mat)))
table(clus)

# All appear to be doublets. Blacklist all

#JCHAIN_blacklist <- names(which(clus==1))
JCHAIN_blacklist <- names(clus)

seu@meta.data$JCHAIN_blacklist <- F
seu@meta.data[JCHAIN_blacklist,]$JCHAIN_blacklist <- T
seu@meta.data[JCHAIN_blacklist,]$is_doublet <- T

pdf(file=paste0(outdir,"JCHAIN_blacklist.pdf"))
ggplot(seu@meta.data[!is.na(seu@meta.data$JCHAIN_clusters),],
       aes(x=nFeature_RNA,y=nCount_RNA, color = JCHAIN_blacklist)) +
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)
dev.off()

#min_nGene_JCHAIN <- summary(annot[names(which(clus!=1)),]$nFeature_RNA)[["1st Qu."]]
#min_nUMI_JCHAIN <- summary(annot[names(which(clus!=1)),]$nCount_RNA)[["1st Qu."]]

# Final filtering --------------------------------------------------------------
library(miQC)
library(flexmix)

sce.list <- list()
for (sname in unique(colData(sce)$orig.ident)){
  sce.list[[sname]] <- sce[,colData(sce)$orig.ident==sname]
  
  pdf(file=paste0(outdir,"miQC_plotMetrics.",sname,".pdf"),width = 10,height = 10)
  print(plotMetrics(sce.list[[sname]]))
  dev.off()
}

sce.filtered.list <- list()
model.list <- list()
for (sname in unique(colData(sce)$orig.ident)){
  if (sname == "PDX11graft"){
    model.list[[sname]] <- mixtureModel(sce.list[[sname]],model_type = "one_dimensional")
  } else {
    model.list[[sname]] <- mixtureModel(sce.list[[sname]])
  }

  pdf(file=paste0(outdir,"miQC_plotFiltering.",sname,".pdf"),width = 10,height = 10)
  print(plotFiltering(sce.list[[sname]], model.list[[sname]],
                      posterior_cutoff = 0.9,
                      keep_all_below_boundary = TRUE))
  dev.off()
  
  pdf(file=paste0(outdir,"miQC_plotModel.",sname,".pdf"),width = 10,height = 10)
  print(plotModel(sce.list[[sname]], model.list[[sname]]))
  dev.off()
  
  
  sce.filtered.list[[sname]] <- filterCells(sce.list[[sname]],
                                            model.list[[sname]],
                                            posterior_cutoff = 0.9,
                                            keep_all_below_boundary = TRUE)
}

miqc.list <- list()
for (sname in names(sce.filtered.list)){
  miqc.list[[sname]] <- rownames(colData(sce.filtered.list[[sname]]))
}
miqc <- as.character(unlist(miqc.list))

seu@meta.data$miqc_kept <- rownames(seu@meta.data) %in% miqc

pdf(file=paste0(outdir,"miQC_kept.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "miqc_kept")
dev.off()

pdf(file=paste0(outdir,"percent_mito_miQC_kept.pdf"),width = 10,height = 10)
ggplot(seu@meta.data,aes(x=nFeature_RNA,y=percent.mt,color=miqc_kept)) +
  geom_point() + theme_classic() + vertical_xlabels +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 0.3) +
  facet_wrap(~ orig.ident)
dev.off()

seu@meta.data$nFeature_RNA_le_500 <- seu@meta.data$nFeature_RNA < 500
pdf(file=paste0(outdir,"nFeature_RNA_le_500.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "nFeature_RNA_le_500")
dev.off()

seu@meta.data$nFeature_RNA_le_200 <- seu@meta.data$nFeature_RNA < 200
pdf(file=paste0(outdir,"nFeature_RNA_le_200.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "nFeature_RNA_le_200")
dev.off()

seu@meta.data$nFeature_RNA_le_100 <- seu@meta.data$nFeature_RNA < 100
pdf(file=paste0(outdir,"nFeature_RNA_le_100.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "nFeature_RNA_le_100")
dev.off()

seu@meta.data$nFeature_RNA_le_50 <- seu@meta.data$nFeature_RNA < 50
pdf(file=paste0(outdir,"nFeature_RNA_le_50.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "nFeature_RNA_le_50")
dev.off()

pdf(file=paste0(outdir,"percent.top_100.pdf"),width = 10,height = 10)
FeaturePlot(object = seu, features = "percent.top_100")
dev.off()

seu@meta.data$percent.top_100_gr_80 <- seu@meta.data$percent.top_100 > 80
pdf(file=paste0(outdir,"percent.top_100_gr_80.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "percent.top_100_gr_80")
dev.off()

# idx_erythrocyte <- seu@meta.data$HBA1 > 0 &
#   !seu@meta.data$hba1_blacklist &
#   seu@meta.data[["nFeature_RNA"]] > min_nGene_erythrocytes &
#   seu@meta.data[["nCount_RNA"]] > min_nUMI_erythrocytes
# 
# idx_plasma_cell <- !is.na(seu@meta.data$JCHAIN_clusters) &
#   !seu@meta.data$JCHAIN_blacklist &
#   seu@meta.data[["nFeature_RNA"]] > min_nGene_JCHAIN &
#   seu@meta.data[["nCount_RNA"]] > min_nUMI_JCHAIN

idx_keep <- seu@meta.data[["nFeature_RNA"]] > 200 & 
               seu@meta.data[["nCount_RNA"]] > 500 & 
               seu@meta.data[["percent.top_100"]] < 80 & 
               seu@meta.data$miqc_kept & 
               seu@meta.data[["percent.ribo"]] < 40 & 
               !seu@meta.data$is_doublet & 
               !seu@meta.data$hba1_blacklist & 
               !seu@meta.data$JCHAIN_blacklist

table(seu@meta.data$miqc_kept)
table(idx_keep)
table(idx_keep,seu@meta.data$orig.ident)

cells_keep <- rownames(seu@meta.data[idx_keep,])

seu@meta.data$cells_filtered <- ! rownames(seu@meta.data) %in%
  cells_keep
table(seu@meta.data$cells_filtered)
table(seu@meta.data$cells_filtered)["TRUE"]/nrow(seu@meta.data)
table(seu@meta.data$cells_filtered,
      seu@meta.data$orig.ident)

pdf(file=paste0(outdir,"cells_filtered.pdf"),width = 25,height = 10)
DimPlot(object = seu, group.by = "cells_filtered", split.by = "cells_filtered")
dev.off()


seu.filt <- subset(seu,cells_filtered==F)
saveRDS(seu.filt, paste0(outdir, "seu.filt.rda"))

# FastMNN ----------------------------------------------------------------------

# boxplot(seu.filt@meta.data$nFeature_RNA[seu.filt@assays$RNA@counts["NCR1",] > 0 & 
#                                           seu.filt@assays$RNA@counts["CD3D",] > 0],
#         seu.filt@meta.data$nFeature_RNA[seu.filt@assays$RNA@counts["NCR1",] == 0 & 
#                                           seu.filt@assays$RNA@counts["CD3D",] > 0])

idx_remove <- seu.filt@assays$RNA@counts["NCR1",] > 0 & 
  (seu.filt@assays$RNA@counts["CD3D",] > 0 | !is.na(seu.filt@meta.data$cdr3))
seu.filt$cells_filtered[idx_remove] <- T
seu.filt <- subset(seu.filt, cells_filtered==F)

seu.list <- SplitObject(seu.filt, split.by = "orig.ident")
seu.filt <- RunFastMNN(object.list = seu.list,
                       features = 4000, auto.merge = T, d = 150, k = 15)
seu.filt <- RunUMAP(seu.filt, reduction = "mnn", dims = 1:150)
#saveRDS(seu.filt, paste0(outdir, "seu.filt.rda"))

pdf(file=paste0(outdir,"markers_post_filtering_mnn.pdf"),width = 40,height = 40)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("MLANA","KIT","CD3D","NCR1","CD4","CD8A","CD14","HBA1",
                         "JCHAIN","MS4A1","FOXP3","CD34","MKI67",
                         "CD1C","THBD","FCGR3A","FCGR3B","CLEC4C",
                         "IL3RA","ITGAX","SIRPA","CLEC10A","SERPINF1",
                         "LILRA4"))
dev.off()

pdf(file=paste0(outdir,"t_cell_markers_post_filtering_mnn.pdf"),width = 40,height = 40)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("IL7R","TCF7","LAG3","HAVCR2","TIGIT","CD68",
                         "HLA-DRA","PDCD1","GZMA","GZMH","TNFRSF9"))
dev.off()

pdf(file=paste0(outdir,"melanoma_markers_post_filtering_mnn.pdf"),width = 30,height = 30)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("rna_MLANA","rna_PMEL","TYR","CPXM1","CDKN2A","BAP1",
                         "CHGA","CHGB"))
dev.off()

pdf(file = paste0(outdir,"percent_mt_ribo_top100_post_filtering_mnn.pdf"), width = 30, height = 30)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("percent.mt","percent.ribo","percent.top_100"))
dev.off()

seu.filt@meta.data$is_double_pos <- rownames(seu.filt@meta.data) %in% names(which(seu.filt@assays$RNA@counts["CD4",] > 0 & (seu.filt@assays$RNA@counts["CD8A",] > 0 | seu.filt@assays$RNA@counts["CD8B",] > 0)))

pdf(file=paste0(outdir,"double_pot_post_filtering_mnn.pdf"),width = 10,height = 10)
DimPlot(object = seu.filt, reduction = "umap", pt.size = 1,order = T,
        group.by = "is_double_pos")
dev.off()

pdf(file=paste0(outdir,"double_pos_vln.pdf"),width = 10,height = 10)
VlnPlot(object = seu.filt, features = "nFeature_RNA",group.by = "is_double_pos")
dev.off()

pdf(file=paste0(outdir,"double_pos_vln.nCount.pdf"),width = 10,height = 10)
VlnPlot(object = seu.filt, features = "nCount_RNA",group.by = "is_double_pos")
dev.off()

seu.filt@meta.data$CD4_pos <- rownames(seu.filt@meta.data) %in% names(which(seu.filt@assays$RNA@counts["CD4",] > 0 & ! (seu.filt@assays$RNA@counts["CD8A",] > 0 | seu.filt@assays$RNA@counts["CD8B",] > 0)))

pdf(file=paste0(outdir,"CD4_pos_vln.nFeature_RNA.pdf"),width = 10,height = 10)
VlnPlot(object = seu.filt, features = "nFeature_RNA",group.by = "CD4_pos")
dev.off()

pdf(file=paste0(outdir,"CD4_pos__vln.nCount.pdf"),width = 10,height = 10)
VlnPlot(object = seu.filt, features = "nCount_RNA",group.by = "CD4_pos")
dev.off()

# CD4+CD8+ double positive cells have substatially higher counts / features and are likely doublets
seu.filt <- subset(seu.filt,is_double_pos==F)

#seu.filt <- readRDS(paste0(outdir, "seu.filt.rda"))
seu.list <- SplitObject(seu.filt, split.by = "orig.ident")
seu.filt <- RunFastMNN(object.list = seu.list,
                       features = 4000, auto.merge = T, d = 150, k = 10)
seu.filt <- RunUMAP(seu.filt, reduction = "mnn", dims = 1:150, n.neighbors = 10)

seu.filt <- FindNeighbors(seu.filt, reduction = "mnn", dims = 1:150, k.param = 10)
#seu.filt.bak <- seu.filt
seu.filt <- FindClusters(seu.filt, resolution = 2)
#saveRDS(seu.filt,paste0(outdir, "seu.filt.rda"))

pdf(file = paste0(outdir,"clusters.pdf"), width = 20, height = 20)
DimPlot(seu.filt, label = T)
dev.off()

markers.all <- FindAllMarkers(seu.filt)

library(WriteXLS)
WriteXLS(markers.all, ExcelFileName = paste0(outdir,"markers.all.res_2.xlsx"),
         BoldHeaderRow = T, AutoFilter = T, AdjWidth = T)

# Annotate clusters ------------------------------------------------------------
#outdir <- "~/proj/um_ss/Investigations/pdx_human/results/v1/"
#seu.filt <- readRDS("~/proj/um_ss/Investigations/pdx_human/results/v1/seu.filt.rda")

seu.filt@meta.data$CD3D <- NULL
seu.filt@meta.data$NCAM1 <- NULL
seu.filt@meta.data$HBA1 <- NULL
seu.filt@meta.data$PMEL <- NULL
seu.filt@meta.data$MLANA <- NULL

pdf(file=paste0(outdir,"markers_I.pdf"),width = 20,height = 10)
DotPlot(seu.filt,
        features = c("CD3D","CD3G","CD8A","CD4","FOXP3",
                     "NCAM1","NCR1","TYROBP","FCGR3A",
                     "PDCD1","LAG3","HAVCR2","TIGIT","CD68","HLA-DRA",
                     "CD19","MS4A1",
                     "CD1C","THBD","TNFRSF21","JCHAIN",
                     "ITGAX","LILRA4","SERPINF1",
                     "CD14","LYZ", "S100A8", "FCN1","S100A9", "CSF3R","CDKN1C",
                     "LILRB2","ITGAL","ITGAE","ITGAM",
                     "MZB1","ITM2C","PLD4","CD74","FCGR3B",
                     "CD86","CSF1R","C1QA","C1QC","C1QB","SPI1","FCER1A",
                     "CLEC10A","CLEC9A","XCR1","SIRPA","CCR5",
                     "FCGR2A","MERTK","FCGR1A","ANPEP","MRC1",
                     "HBA1","HBA2","HBB",
                     "PF4","PPBP",
                     "CD34","PECAM1","ACTA1",
                     "PMEL","MLANA","TYR","MKI67","ALB"),
        cluster.idents = F, assay = "RNA") + vertical_xlabels
dev.off()

pdf(file=paste0(outdir,"markers_II.pdf"),width = 20,height = 10)
DotPlot(seu.filt,features=c(
  "CD34","CD38","KIT","PTPRC","PROM1",
  "NT5E","THY1","ENG","MCAM","ITGB1","CD44","NGFR",
  "NCAM1","MYF5","DES","CDH15","VCAM1","KDR","FLT1","MYOD1",
  "SELL","KERA","ALDH1A1","ALDH2",
  "VIM","GJA1","PDGFRB",
  "CD80","CD86","HLA-A","HLA-B","HLA-DRA",
  "ITGA6","MME","TFRC","S100A4","DKK3",
  "CDH5","CD14","PECAM1","VWF","FLT4",
  "PDGFRA","ACTA2","SERPINH1","COL1A1")
  ,cluster.idents = F, assay = "RNA") + vertical_xlabels
dev.off()

pdf(file=paste0(outdir,"markers_III.pdf"),width = 20,height = 10)
DotPlot(seu.filt,features=c(
  "CD34","CD38","KIT","PTPRC","PROM1",
  "NT5E","THY1","ENG","MCAM","ITGB1","CD44","NGFR",
  "NCAM1","MYF5","DES","CDH15","VCAM1","KDR","FLT1","MYOD1",
  "SELL","KERA","ALDH1A1","ALDH2",
  "VIM","GJA1","PDGFRB",
  "CD80","CD86","HLA-A","HLA-B","HLA-DRA",
  "ITGA6","MME","TFRC","S100A4","DKK3",
  "CDH5","CD14","PECAM1","VWF","FLT4",
  "PDGFRA","ACTA2","SERPINH1","COL1A1","MLANA","PMEL","TYR","CPXM1","MKI67",
  "VGF")
  ,cluster.idents = F, assay = "RNA") + vertical_xlabels
dev.off()

pdf(file=paste0(outdir,"main_markers_post_filtering_mnn.pdf"),width = 30,height = 30)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("CD3D","CD8A","CD8B","CD4","FOXP3","TYR","MLANA",
                         "PMEL","NCR1","NCAM1","THBD","MAML3",
                         "HMCN1","DIP2C","MKI67"),order=T)
dev.off()


pdf(file=paste0(outdir,"um_classes_post_filtering_mnn.pdf"),width = 30,height = 30)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("CDH1","ECM1","EIF1B","FXR1","HTR2B","ID2","LMCD1",
                         "LTA4H","MTUS1","RAB31","ROBO1","SAT1","PRAME"),order=T)
dev.off()

pdf(file=paste0(outdir,"rna_content_post_filtering_mnn.pdf"),width = 30,height = 30)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("nFeature_RNA","nCount_RNA","percent.mt",
                         "percent.ribo","percent.top_100"),order=T)
dev.off()

pdf(file=paste0(outdir,"nCount_RNA_vln.nCount.pdf"),width = 10,height = 10)
VlnPlot(object = seu.filt, features = "nCount_RNA",group.by = "seurat_clusters")
dev.off()

pdf(file=paste0(outdir,"nFeature_RNA_vln.nCount.pdf"),width = 10,height = 10)
VlnPlot(object = seu.filt, features = "nCount_RNA",group.by = "seurat_clusters")
dev.off()

seu.filt <- RenameIdents(seu.filt,
                    `0` = "UM",
                    `1` = "UM",
                    `2` = "UM",
                    `3` = "CD8 T cell",
                    `4` = "UM",
                    `5` = "UM (TYR low)",
                    `6` = "UM",
                    `7` = "UM",
                    `8` = "UM",
                    `9` = "UM",
                    `10` = "UM",
                    `11` = "UM",
                    `12` = "UM",
                    `13` = "UM",
                    `14` = "UM",
                    `15` = "CD8 T cell",
                    `16` = "UM",
                    `17` = "UM",
                    `18` = "UM, proliferative",
                    `19` = "UM",
                    `20` = "UM",
                    `21` = "UM",
                    `22` = "UM",
                    `23` = "T, proliferative",
                    `24` = "CD8 T cell",
                    `25` = "CD4 T cell",
                    `26` = "UM (TYR low)",
                    `27` = "UM, proliferative",
                    `28` = "UM",
                    `29` = "UM"
)
seu.filt@meta.data$clusters.annotated <- Idents(seu.filt)

saveRDS(seu.filt, paste0(outdir,"seu.filt.rda"))

pdf(file = paste0(outdir,"clusters.annotated.pdf"), width = 20, height = 20)
DimPlot(seu.filt, label = T)
dev.off()

pdf(file=paste0(outdir,"main_markers_post_filtering_mnn.annotated.pdf"),width = 10,height = 6)
DotPlot(object = seu.filt,
            features = c("CD3D","CD8A","CD8B","CD4","FOXP3","TYR","MLANA",
                         "PMEL","NCR1","NCAM1","THBD","MAML3",
                         "HMCN1","DIP2C","KIT","VGF","MKI67"),
        cluster.idents = F, assay = "RNA") + vertical_xlabels
dev.off()


# Add metadata -----------------------------------------------------------------
library(readxl)

#seu.filt <- readRDS(paste0(outdir,"seu.filt.rda"))

samples_10x <- read_excel(path="~/proj/um_ss/Investigations/samples_10x.xlsx")
samples_10x <- data.frame(samples_10x)

seu.filt@meta.data$Sample.ID <- gsub("graft","",seu.filt@meta.data$orig.ident)

meta <- seu.filt@meta.data
meta$Rownames <- rownames(meta)
meta <- merge(meta, samples_10x, by="Sample.ID",all.x=T,all.y=F)
rownames(meta) <- meta$Rownames
meta$Rownames <- NULL
meta <- meta[rownames(seu.filt@meta.data),]
stopifnot(identical(rownames(meta),rownames(seu.filt@meta.data)))
seu.filt@meta.data <- meta

seu.filt@meta.data$Mouse_biopsy_site_and_TIL_selection <- paste0(
  seu.filt@meta.data$Mouse_biopsy_site,"_",seu.filt@meta.data$TIL_selection)

#saveRDS(seu.filt, paste0(outdir,"seu.filt.rda"))

DimPlot(seu.filt, split.by = "Mouse_biopsy_site")
DimPlot(seu.filt, split.by = "TIL_selection")
DimPlot(seu.filt, split.by = "Sample.ID")

pdf(file = paste0(outdir, "UMAP_Mouse_biopsy_site_and_TIL_selection.pdf"),
    width = 12, height = 3.5)
DimPlot(seu.filt, split.by = "Mouse_biopsy_site_and_TIL_selection")
dev.off()

# markers_plot <- c("PTPRC","CD3D","CD8A","CD8B",
#                   "CD69","NR4A1","PDCD1","ENTPD1","TNFRSF9","CTLA4","HAVCR2","LAG3","TIGIT","KIR2DL4",
#                   "CD38","HLA-DRA",
#                   "NKG7","CCL5",
#                   "PRF1","GZMA","GZMB","GZMH","GZMK",
#                   "IL32","CCL4",
#                   "CD27","ITM2A","TRAC","TRBC2",
#                   "ISG15","CD5","GNG4",
#                   "IFNG","TNF","XCL2","XCL1",
#                   "CXCL13","CCL3","CX3CR1",
#                   "GNLY","FCRL6","FGFBP2","FCGR3A","KLRG1","KLRB1","NCR3","RORC",
#                   "IL7R","TCF7","SATB1","LEF1","CCR7","CCR5",
#                   "KLF2","TOX","EOMES","TBX21","ZNF683",
#                   "TNFSF14","CD160",
#                   "MKI67","CDKN1A",
#                   "ITGA1","CXCR6","CXCR3","ITGAE","ITGAL","PRDM1",
#                   "CD44","RUNX3","S1PR1","SELL")
# 
# DotPlot(seu.filt,
#         features = markers_plot,
#         cluster.idents = T, assay = "RNA",scale = T,
#         group.by="Mouse_biopsy_site_and_TIL_selection") + vertical_xlabels

FeaturePlot(seu.filt, features = c("CX3CR1","XCL1","XCL2","KLF2"))

FeaturePlot(seu.filt, features = c("ALB","MUC1","KRT1","CD163"))

FeaturePlot(seu.filt, features = c("KRT7","KRT18","CXCL1","KRT19","SCGB3A1",
                                   "PIGR","APOC1","FGA","RGS5","NGFR","CXCL12",
                                   "PDGFRA","IGFBP3","DCN","LUM","IL6","THY1",
                                   "SOD2","MSLN","MYH11","CNN1"))

cells <- rownames(seu.filt@meta.data)[!is.na(seu.filt@meta.data$CD8_detailed)]
seu.filt@meta.data$cell_cd8 <- rownames(seu.filt@meta.data) %in% cells

seu.filt.cd8_t <- subset(seu.filt,cell_cd8==T)

# Liver MART1 vs liver bulk
cells <- rownames(seu.filt.cd8_t@meta.data)[
  seu.filt.cd8_t@meta.data$Mouse_biopsy_site_and_TIL_selection %in% c("Liver_MART1","Liver_Bulk")]

seu.filt.cd8_t@meta.data$is_liver <- rownames(seu.filt.cd8_t@meta.data) %in% cells
seu.filt.cd8_t.liver <- subset(seu.filt.cd8_t, is_liver==T)

cells <- rownames(seu.filt.cd8_t.liver@meta.data)[
  ! seu.filt.cd8_t.liver@meta.data$CD8_detailed %in% c("Mixed, Proliferative")]
seu.filt.cd8_t.liver@meta.data$is_nonmixed <- rownames(seu.filt.cd8_t.liver@meta.data) %in% cells

seu.filt.cd8_t.liver.nonmixed <- subset(seu.filt.cd8_t.liver, is_nonmixed==T)

Idents(seu.filt.cd8_t.liver.nonmixed) <- seu.filt.cd8_t.liver.nonmixed@meta.data$Mouse_biopsy_site_and_TIL_selection

#seu.filt.cd8_t.liver.nonmixed.male <- subset(seu.filt.cd8_t.liver.nonmixed,Sex=="M"): All are female

#"if you are comparing two groups ident.2 will be your reference ( I assume control for you)": https://github.com/satijalab/seurat/issues/5127
markers.liver.mart1_vs_bulk <- FindMarkers(seu.filt.cd8_t.liver.nonmixed, 
                                            ident.1 = "Liver_MART1",
                                            ident.2 = "Liver_Bulk")

WriteXLS(markers.liver.mart1_vs_bulk,
         ExcelFileName = paste0(outdir,"markers.liver.mart1_vs_bulk.xlsx"),
         AdjWidth = T, AutoFilter = T,row.names = T,BoldHeaderRow = T)

markers.liver.mart1_vs_bulk.filtered <- markers.liver.mart1_vs_bulk[!grepl("^MT",rownames(markersl.liver.mart1_vs_bulk)) & 
                              !grepl("^TR[AB][VCJ]",rownames(markersl.liver.mart1_vs_bulk)),]


markers.liver.mart1_vs_bulk.filtered$gene <- rownames(markers.liver.mart1_vs_bulk.filtered )

library(ggrepel)
ggplot(markers.liver.mart1_vs_bulk.filtered,aes(x=avg_log2FC,y=-log10(p_val_adj))) + 
  geom_point() + geom_text_repel(aes(label=gene),
                           data = markers.liver.mart1_vs_bulk.filtered[
                             abs(markers.liver.mart1_vs_bulk.filtered$avg_log2FC) > 0.4 & 
                               -log10(markers.liver.mart1_vs_bulk.filtered$p_val_adj) > 20,
                           ],max.overlaps = Inf,min.segment.length = 0) + 
  xlab("Average Log2 fold change") + ylab("-log10 q") + theme_minimal()
ggsave(file=paste0(outdir,"pdx_liver_mart1_vs_bulk.volcano.pdf"), width = 6, height = 6)

ggplot(markers.liver.mart1_vs_bulk.filtered,aes(x=avg_log2FC,y=-log10(p_val_adj))) + 
  geom_point() + geom_text_repel(aes(label=gene),
                                 data = markers.liver.mart1_vs_bulk.filtered[
                                   markers.liver.mart1_vs_bulk.filtered$gene %in% 
                                     c("ENTPD1","ICOS","TNFRSF9","KIR2DL4") | 
                                     abs(markers.liver.mart1_vs_bulk.filtered$avg_log2FC) > 0.4 & 
                                     -log10(markers.liver.mart1_vs_bulk.filtered$p_val_adj) > 30,
                                 ],max.overlaps = Inf,min.segment.length = 0) + 
  xlab("Average Log2 fold change") + ylab("-log10 q") + theme_minimal()
ggsave(file=paste0(outdir,"pdx_liver_mart1_vs_bulk.volcan.highlight_2.pdf"), width = 6, height = 6)

# Spleen MART1 vs spleen bulk
cells <- rownames(seu.filt.cd8_t@meta.data)[
  seu.filt.cd8_t@meta.data$Mouse_biopsy_site_and_TIL_selection %in% c("Spleen_MART1","Spleen_Bulk")]

seu.filt.cd8_t@meta.data$is_spleen <- rownames(seu.filt.cd8_t@meta.data) %in% cells
seu.filt.cd8_t.spleen <- subset(seu.filt.cd8_t, is_spleen==T)

cells <- rownames(seu.filt.cd8_t.spleen@meta.data)[
  ! seu.filt.cd8_t.spleen@meta.data$CD8_detailed %in% c("Mixed, Proliferative")]
seu.filt.cd8_t.spleen@meta.data$is_nonmixed <- rownames(seu.filt.cd8_t.spleen@meta.data) %in% cells

seu.filt.cd8_t.spleen.nonmixed <- subset(seu.filt.cd8_t.spleen, is_nonmixed==T)

Idents(seu.filt.cd8_t.spleen.nonmixed) <- seu.filt.cd8_t.spleen.nonmixed@meta.data$Mouse_biopsy_site_and_TIL_selection

markers.liver_m1_spleen_m1 <- FindMarkers(seu.filt.cd8_t.spleen.nonmixed, 
                                           ident.1 = "Spleen_MART1",
                                           ident.2 = "Spleen_Bulk")

WriteXLS(markers.liver_m1_spleen_m1,
         ExcelFileName = paste0(outdir,"markers.liver_m1_spleen_m1.xlsx"),
         AdjWidth = T, AutoFilter = T,row.names = T,BoldHeaderRow = T)

markers.liver_m1_spleen_m1$gene <- rownames(markers.liver_m1_spleen_m1)
ggplot(markers.liver_m1_spleen_m1,aes(x=avg_log2FC,y=-log10(p_val_adj))) + 
  geom_point() + geom_text_repel(aes(label=gene),
                                 data = markers.liver_m1_spleen_m1[
                                   abs(markers.liver_m1_spleen_m1$avg_log2FC) > 0.4 & 
                                     -log10(markers.liver_m1_spleen_m1$p_val_adj) > 20,
                                 ],max.overlaps = Inf,min.segment.length = 0) + 
  xlab("Average Log2 fold change") + ylab("-log10 q") + theme_minimal()
ggsave(file=paste0(outdir,"pdx_spleen_mart1_vs_bulk.volcano.pdf"), width = 6, height = 6)

# Liver MART1 vs spleen MART1
cells <- rownames(seu.filt.cd8_t@meta.data)[
  seu.filt.cd8_t@meta.data$Mouse_biopsy_site_and_TIL_selection %in% c("Liver_MART1","Spleen_MART1")]

seu.filt.cd8_t@meta.data$is_liver_m1_spleen_m1 <- rownames(seu.filt.cd8_t@meta.data) %in% cells
seu.filt.cd8_t.liver_m1_spleen_m1 <- subset(seu.filt.cd8_t, is_liver_m1_spleen_m1==T)

cells <- rownames(seu.filt.cd8_t.liver_m1_spleen_m1@meta.data)[
  ! seu.filt.cd8_t.liver_m1_spleen_m1@meta.data$CD8_detailed %in% c("Mixed, Proliferative")]
seu.filt.cd8_t.liver_m1_spleen_m1@meta.data$is_nonmixed <- rownames(seu.filt.cd8_t.liver_m1_spleen_m1@meta.data) %in% cells

seu.filt.cd8_t.liver_m1_spleen_m1.nonmixed <- subset(seu.filt.cd8_t.liver_m1_spleen_m1, is_nonmixed==T)

Idents(seu.filt.cd8_t.liver_m1_spleen_m1.nonmixed) <- seu.filt.cd8_t.liver_m1_spleen_m1.nonmixed@meta.data$Mouse_biopsy_site_and_TIL_selection

markers.liver_m1_spleen_m1 <- FindMarkers(seu.filt.cd8_t.liver_m1_spleen_m1.nonmixed, 
                                            ident.1 = "Liver_MART1",
                                            ident.2 = "Spleen_MART1")

WriteXLS(markers.liver_m1_spleen_m1,
         ExcelFileName = paste0(outdir,"markers.liver_m1_spleen_m1.xlsx"),
         AdjWidth = T, AutoFilter = T,row.names = T,BoldHeaderRow = T)

markers.liver_m1_spleen_m1$gene <- rownames(markers.liver_m1_spleen_m1)
ggplot(markers.liver_m1_spleen_m1,aes(x=avg_log2FC,y=-log10(p_val_adj))) + 
  geom_point() + geom_text_repel(aes(label=gene),
                                 data = markers.liver_m1_spleen_m1[
                                   abs(markers.liver_m1_spleen_m1$avg_log2FC) > 1 & 
                                     -log10(markers.liver_m1_spleen_m1$p_val_adj) > 100,
                                 ],max.overlaps = Inf,min.segment.length = 0) + 
  xlab("Average Log2 fold change") + ylab("-log10 q") + theme_minimal()
ggsave(file=paste0(outdir,"pdx_liver_m1_spleen_m1.volcano.pdf"), width = 6, height = 6)

# Liver bulk vs spleen bulk
cells <- rownames(seu.filt.cd8_t@meta.data)[
  seu.filt.cd8_t@meta.data$Mouse_biopsy_site_and_TIL_selection %in% c("Liver_Bulk","Spleen_Bulk")]

seu.filt.cd8_t@meta.data$is_liver_bulk_spleen_bulk <- rownames(seu.filt.cd8_t@meta.data) %in% cells
seu.filt.cd8_t.liver_bulk_spleen_bulk <- subset(seu.filt.cd8_t, is_liver_bulk_spleen_bulk==T)

cells <- rownames(seu.filt.cd8_t.liver_bulk_spleen_bulk@meta.data)[
  ! seu.filt.cd8_t.liver_bulk_spleen_bulk@meta.data$CD8_detailed %in% c("Mixed, Proliferative")]
seu.filt.cd8_t.liver_bulk_spleen_bulk@meta.data$is_nonmixed <- rownames(seu.filt.cd8_t.liver_bulk_spleen_bulk@meta.data) %in% cells

seu.filt.cd8_t.liver_bulk_spleen_bulk.nonmixed <- subset(seu.filt.cd8_t.liver_bulk_spleen_bulk, is_nonmixed==T)

Idents(seu.filt.cd8_t.liver_bulk_spleen_bulk.nonmixed) <- seu.filt.cd8_t.liver_bulk_spleen_bulk.nonmixed@meta.data$Mouse_biopsy_site_and_TIL_selection

markers.is_liver_bulk_spleen_bulk <- FindMarkers(seu.filt.cd8_t.liver_bulk_spleen_bulk.nonmixed, 
                                          ident.1 = "Liver_Bulk",
                                          ident.2 = "Spleen_Bulk")

WriteXLS(markers.is_liver_bulk_spleen_bulk,
         ExcelFileName = paste0(outdir,"markers.is_liver_bulk_spleen_bulk.xlsx"),
         AdjWidth = T, AutoFilter = T, row.names = T, BoldHeaderRow = T, SheetNames = "Null")

markers.is_liver_bulk_spleen_bulk$gene <- rownames(markers.is_liver_bulk_spleen_bulk)
ggplot(markers.is_liver_bulk_spleen_bulk, aes(x=avg_log2FC,y=-log10(p_val_adj))) + 
  geom_point() + geom_text_repel(aes(label=gene),
                                 data = markers.is_liver_bulk_spleen_bulk[
                                   abs(markers.is_liver_bulk_spleen_bulk$avg_log2FC) > 1 & 
                                     -log10(markers.is_liver_bulk_spleen_bulk$p_val_adj) > 50,
                                 ],max.overlaps = Inf,min.segment.length = 0) + 
  xlab("Average Log2 fold change") + ylab("-log10 q") + theme_minimal()
ggsave(file=paste0(outdir,"pdx_liver_bulk_spleen_bulk.volcano.pdf"), width = 6, height = 6)

# Liver MART1 vs spleen bulk



# Sub-clustering, CD8 T --------------------------------------------------------

#seu.filt <- readRDS(paste0(outdir, "seu.filt.rda"))

cd8_t_clusters <- sort(as.numeric(as.character(unique(seu.filt@meta.data$seurat_clusters[
  grepl("^CD8",seu.filt@active.ident) | grepl("^T, proliferative",seu.filt@active.ident)]))))
seu.filt.cd8_t <- subset(seu.filt,seurat_clusters %in% cd8_t_clusters)
seu.filt.cd8_t <- subset(seu.filt.cd8_t,CD4 == 0 & NCR1 == 0 &
                           (CD3D > 0 | CD3G > 0 | CD8A > 0 | CD8B >0) &
                           MLANA == 0 & CD19 == 0 & CD14 == 0)
seu.filt.cd8_t <- FindNeighbors(seu.filt.cd8_t,reduction = "mnn", k.param = 5, dims=1:150)
seu.filt.cd8_t <- FindClusters(seu.filt.cd8_t, algorithm = 1, resolution = 1)
seu.filt.cd8_t <- RunUMAP(seu.filt.cd8_t, reduction = "mnn", dims=1:150,
                          min.dist = 0.05, n.neighbors = 5)

#https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00202-4/fulltext
#https://www.cell.com/trends/immunology/fulltext/S1471-4906(16)30123-5
#https://pubmed.ncbi.nlm.nih.gov/30595452/
#https://www.nature.com/articles/s41423-019-0344-8/tables/1

DimPlot(seu.filt.cd8_t, label = T)

markers.cd8_t <- FindAllMarkers(seu.filt.cd8_t)
WriteXLS(markers.cd8_t,ExcelFileName = paste0(outdir,"markers.cd8_t.xlsx"),
         BoldHeaderRow = T,AutoFilter = T,AdjWidth = T)

DotPlot(seu.filt.cd8_t,
        features =markers_plot,
        cluster.idents = T, assay = "RNA",scale = T, group.by="seurat_clusters") + vertical_xlabels


# Early activated and activated defined based on NR4A1+TNF
# (see discussion of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522381/; Nur77)
# Essentially, most liver T cells express CD69, irregardless of whether recent activation or not;
# Nur77 is supposedly a marker for recent activation for these instead.
# "That the expression of CD69 was not promoted by recent TCR cross-linking"
# "in the liver immune activation is not the primary mechanism of HLA-DR/CD38 upregulation"

markers_plot <- c("PTPRC","CD3D","CD8A","CD8B",
                  "CD69","NR4A1","PDCD1","ENTPD1","TNFRSF9","CTLA4","HAVCR2","LAG3","TIGIT","KIR2DL4",
                  "CD38","HLA-DRA",
                  "NKG7","CCL5",
                  "PRF1","GZMA","GZMB","GZMH","GZMK",
                  "IL32","CCL4",
                  "CD27","ITM2A","TRAC","TRBC2",
                  "ISG15","CD5","GNG4",
                  "IFNG","TNF","XCL2","XCL1",
                  "CXCL13","CCL3","CX3CR1",
                  "GNLY","FCRL6","FGFBP2","FCGR3A","KLRG1","KLRB1","NCR3","RORC",
                  "IL7R","TCF7","SATB1","LEF1","CCR7","CCR5",
                  "KLF2","TOX","EOMES","TBX21","ZNF683",
                  "TNFSF14","CD160",
                  "MKI67","CDKN1A",
                  "ITGA1","CXCR6","CXCR3","ITGAE","ITGAL","PRDM1",
                  "CD44","RUNX3","S1PR1","SELL")#,
                  #"SLCO3A1","RPTOR","SSBP3"
                  #)#,"NCAM1","C1QBP","PRDX1","NME1","NME2","KLRK1","ATXN2")

# Naive-like: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7767439/
# TRM: https://rupress.org/jem/article/219/4/e20212059/213066/Modulation-of-tissue-resident-memory-T-cells-by
# TRM: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7696520/
# Liver T cells: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522381/

pdf(file=paste0(outdir,"dotplot_cd8_markers.clusters.pdf"), width = 13, height = 5)
DotPlot(seu.filt.cd8_t,
        features =markers_plot,
        cluster.idents = T, assay = "RNA",scale = T, group.by="seurat_clusters") + vertical_xlabels
dev.off()

pdf(file=paste0(outdir,"featureplot_cd8_markers.pdf"), width = 15, height = 10)
FeaturePlot(seu.filt.cd8_t,
            features = c("PDCD1","HAVCR2","TNFRSF9","ICOS",
                         "HLA-DRA","IL7R","TCF7","NR4A1",
                         "PRF1","GZMB","KLRB1","CCR7"),ncol = 4)
dev.off()

#Idents(seu.filt.cd8_t) <- seu.filt.cd8_t@meta.data$seurat_clusters

seu.filt.cd8_t <- RenameIdents(seu.filt.cd8_t,
                               `0` = "XCL1+XLC2+",
                               `1` = "CX3CR1+TCF7+",
                               `2` = "XCL1+XLC2+",
                               `3` = "XCL1+XLC2+",
                               `4` = "CX3CR1+TCF7+",
                               `5` = "Mixed, Proliferative",
                               `6` = "XCL1+XLC2+",
                               `7` = "CX3CR1+TCF7+",
                               `8` = "XCL1+XLC2+")

pdf(file=paste0(outdir,"dotplot_cd8_markers.labels.pdf"), width = 15, height = 4)
DotPlot(seu.filt.cd8_t,
        features =markers_plot,
        cluster.idents = T, assay = "RNA",scale = F) + vertical_xlabels
dev.off()

pdf(file=paste0(outdir,"cd8_t_clusters.pdf"),width = 11,height = 10)
DimPlot(seu.filt.cd8_t, label = T, group.by="seurat_clusters")
dev.off()

seu.filt.cd8_t@active.ident <- factor(seu.filt.cd8_t@active.ident,
                                      levels=sort(levels(seu.filt.cd8_t@active.ident)))

pdf(file=paste0(outdir,"cd8_t_clusters.labels.pdf"),width = 13,height = 10)
DimPlot(seu.filt.cd8_t, label = T)
dev.off()

saveRDS(seu.filt.cd8_t,file=paste0(outdir,"seu.filt.cd8_t.rda"))

seu.filt <- AddMetaData(seu.filt, metadata = seu.filt.cd8_t@active.ident, 
                        col.name = "CD8_detailed")
saveRDS(seu.filt,file=paste0(outdir,"seu.filt.rda"))
DimPlot(seu.filt,group.by="CD8_detailed",label = T)

# infercnv ---------------------------------------------------------------------
library(infercnv)

#seu.filt <- readRDS(paste0(outdir, "seu.filt.rda"))

annotations_file <- data.frame(sname=seu.filt@meta.data$orig.ident,
                               ct=as.character(seu.filt@active.ident))
idx <- grep("^UM",as.character(annotations_file$ct))
annotations_file$ct[idx] <- paste0("malignant_",annotations_file$sname[idx])
annotations_file$cname <- colnames(seu.filt)
annotations_file <- annotations_file[,c("cname","ct")]
rownames(annotations_file) <- annotations_file$cname
annotations_file$cname <- NULL

infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix=seu.filt@assays$RNA@counts,
  annotations_file=annotations_file,
  delim="\t",
  gene_order_file="/data/local/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes_pos.gtf",
  ref_group_names=grep("^malignant",unique(annotations_file$ct),value=T,invert = T))

dir.create(paste0(outdir,"infercnv"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir=paste0(outdir,"infercnv"),
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             HMM=FALSE,
                             num_threads = 1,
                             analysis_mode="samples",
                             output_format = "pdf",
                             png_res = 600)
