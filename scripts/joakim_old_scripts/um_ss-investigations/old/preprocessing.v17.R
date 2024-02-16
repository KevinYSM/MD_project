# The most important reference for liver CD8 T cells: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522381/

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

outdir <- "/data/proj/um_ss/Investigations/seurat/results/v16/"
dir.create(outdir, showWarnings = F, recursive = T)

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/"
dir.create(outdir, showWarnings = F, recursive = T)

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

fnames <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/multi/*/outs/multi/count/cellbender_feature_bc_matrix/cellbender_feature_bc_matrix_filtered.h5")

snames <- basename(dirname(dirname(dirname(dirname(dirname(fnames))))))
idx <- ! snames %in% c("1A","1B","2A","2B","6A","6B","7A","7B") & 
  ! grepl("^SampleID",snames) & ! grepl("^PDX",snames)
fnames <- fnames[idx]
snames <- snames[idx]

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

#seu <- readRDS(file = paste0(outdir,"seu.rda"))
seu <- RunUMAP(seu, reduction = "mnn", dims = 1:150)

#JCHAIN: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8114807/
pdf(file=paste0(outdir,"markers_pre_doublet_filtering.pdf"),width = 40,height = 40)
FeaturePlot(object = seu, reduction = "umap", pt.size = .1,
            features = c("MLANA","CD3D","NCR1","CD4","CD8A","CD14","HBA1",
                         "JCHAIN","MS4A1","FOXP3","CD34","MKI67",
                         "CD1C","THBD","FCGR3A","FCGR3B","CLEC4C",
                         "IL3RA","ITGAX","SIRPA","CLEC10A","SERPINF1",
                         "LILRA4"))
dev.off()

# Add VDJ data -----------------------------------------------------------------
add_vdj_to_seurat <- function(lig.tumor, vdj_paths){
  library(djvdj)
  fnames <- vdj_paths
  vdj <- list()
  for (nm in names(fnames)) {
    vdj[[nm]] <- djvdj::import_vdj(input = NULL, vdj_dir = fnames[[nm]],
                                   filter_chains = T, 
                                   define_clonotypes = "cdr3_gene")
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

fnames <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/multi/*/outs/per_sample_outs/*/vdj_t/")
fnames <- fnames[basename(dirname(fnames)) %in% unique(seu@meta.data$orig.ident)]
names(fnames) <- basename(dirname(fnames))
fnames <- fnames[names(fnames) %in% seu@meta.data$orig.ident]

seu <- add_vdj_to_seurat(seu, vdj_paths = fnames)

saveRDS(seu,file = paste0(outdir,"/seu.rda"))

# Doublets ---------------------------------------------------------------------
library(scDblFinder)

# The authors do not recommend further QC filtering than this before doublet detection
idx_keep <- seu@meta.data[["nCount_RNA"]] > 200

cells_keep <- rownames(seu@meta.data[idx_keep,])

lig.tumor.filt <- rliger::subsetLiger(lig.tumor, cells.use = cells_keep)
library(Seurat)
seu.subset <- subset(seu, subset = nCount_RNA > 200)

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

seu <- FindNeighbors(seu, reduction = "mnn", dims = 1:150, k.param = 5)
seu <- FindClusters(seu, resolution = 5)

pdf(file=paste0(outdir,"louvain_pre_doublet_cluster_identification.pdf"),width = 10,height = 10)
DimPlot(seu, reduction = "umap",label = T)
dev.off()

pdf(file=paste0(outdir,"louvain_pre_doublet_cluster_identification_samples.pdf"),width = 10,height = 10)
DimPlot(seu, reduction = "umap",label = T,group.by = "orig.ident")
dev.off()


doublet_cluster_fisher <- list()
for (clus in unique(seu@active.ident)){
  cells_clus <- names(seu@active.ident)[
    which(seu@active.ident == clus)]
  
  n_is_doublet_clus <- sum(seu@meta.data[cells_clus,
                                               "is_doublet"])
  n_is_not_doublet_clus <- sum(!seu@meta.data[cells_clus,
                                                    "is_doublet"])
  
  n_is_doublet_rest <- sum(seu@meta.data[
    ! rownames(seu@meta.data) %in% cells_clus,
    "is_doublet"])
  n_is_not_doublet_rest <- sum(!seu@meta.data[
    ! rownames(seu@meta.data) %in% cells_clus,
    "is_doublet"])
  
  tab <- rbind(c(n_is_not_doublet_clus,n_is_doublet_clus),
        c(n_is_not_doublet_rest,n_is_doublet_rest))
  
  ft <- fisher.test(tab)
  
  doublet_cluster_fisher[[clus]] <- data.frame(
    clus=clus, 
    tab=paste0(as.character(tab),collapse="_"), 
    p=ft$p.value, 
    odds_ratio=ft$estimate,
    percent_doublets=n_is_doublet_clus/(n_is_doublet_clus+n_is_not_doublet_clus))
}

doublet_cluster_fisher <- do.call("rbind",doublet_cluster_fisher)
doublet_cluster_fisher$q <- p.adjust(doublet_cluster_fisher$p,method = "BH")

doublet_cluster_fisher$median_nUMI <- NA
doublet_cluster_fisher$median_nGene <- NA
for (clus in unique(seu@active.ident)){
  idx <- seu@active.ident == clus
  median_nUMI <- median(seu@meta.data$nCount_RNA[idx])
  doublet_cluster_fisher[doublet_cluster_fisher$clus == clus,]$median_nUMI <- median_nUMI
  median_nGene <- median(seu@meta.data$nFeature_RNA[idx])
  doublet_cluster_fisher[doublet_cluster_fisher$clus == clus,]$median_nGene <- median_nGene
}

median_nUMI_per_sample <- list()
upper_90 <- list()
upper_95 <- list()
for (sname in unique(seu@meta.data$orig.ident)){
  median_nUMI_per_sample[[sname]] <- median(seu@meta.data$nCount_RNA[
    seu@meta.data$orig.ident==sname])
  upper_90[[sname]] <- quantile(seu@meta.data$nCount_RNA[
    seu@meta.data$orig.ident==sname],seq(0,1,0.05))["90%"]
  upper_95[[sname]] <- quantile(seu@meta.data$nCount_RNA[
    seu@meta.data$orig.ident==sname],seq(0,1,0.05))["95%"]
}

perc_above_thresh <- list()
for (clus in unique(seu@active.ident)){
  idx <- seu@active.ident == clus
  snames <- unique(seu@meta.data[idx,]$orig.ident)
  ntot <- nrow(seu@meta.data[idx,])
  n_above_thresh_sample <- list()
  for (sname in snames){
    n_above_thresh_sample[[sname]] <- sum(seu@meta.data[idx,][
      seu@meta.data[idx,]$orig.ident == sname,]$nCount_RNA > upper_95[[sname]])
  }
  n_above_thresh_tot <- sum(unlist(n_above_thresh_sample))
  perc_above_thresh[[clus]] <- n_above_thresh_tot/ntot
}
perc_above_thresh <- unlist(perc_above_thresh)
perc_above_thresh <- as.data.frame(perc_above_thresh)

doublet_cluster_fisher <- merge(doublet_cluster_fisher,perc_above_thresh,by.x="clus",by.y="row.names")

thresh_doublet_perc <- quantile(doublet_cluster_fisher$percent_doublets,seq(0,1,0.05))["95%"]

thresh_perc_above_thresh <- quantile(doublet_cluster_fisher$perc_above_thresh,seq(0,1,0.05))["95%"]

pdf(file=paste0(outdir,"doublet_cluster_fisher_per_sample_UMI.pdf"),width = 6,height = 6)
plot(doublet_cluster_fisher$percent_doublets,doublet_cluster_fisher$perc_above_thresh)
abline(h=thresh_perc_above_thresh,v=thresh_doublet_perc)
dev.off()

doublet_clusters <- as.numeric(doublet_cluster_fisher[which(
  doublet_cluster_fisher$q < 0.05 & 
    doublet_cluster_fisher$odds_ratio < 1 & 
    (doublet_cluster_fisher$percent_doublets >= thresh_doublet_perc) | 
    (doublet_cluster_fisher$perc_above_thresh > thresh_perc_above_thresh)),]$clus)

cells_doublet_clusters <- names(seu@active.ident[
  seu@active.ident %in% doublet_clusters])

seu@meta.data$in_doublet_cluster <- 
  rownames(seu@meta.data) %in% cells_doublet_clusters

meta <- seu@meta.data
rownames(meta) <- rownames(seu@meta.data)
meta$barcode <- rownames(meta)
meta <- merge(meta,doublet_cluster_fisher[,
          c("clus","percent_doublets","median_nUMI","median_nGene",
            "perc_above_thresh")],
          by.x="seurat_clusters",by.y="clus",all.x=T)
rownames(meta) <- meta$barcode
meta$barcode <- NULL
meta <- meta[rownames(seu@meta.data),]
stopifnot(identical(rownames(meta),rownames(seu@meta.data)))
seu@meta.data <- meta

pdf(file=paste0(outdir,"in_doublet_cluster.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "in_doublet_cluster")
dev.off()

pdf(file=paste0(outdir,"is_doublet.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "is_doublet")
dev.off()

pdf(file=paste0(outdir,"percent_doublets.pdf"),width = 10,height = 10)
FeaturePlot(object = seu, features = "percent_doublets")
dev.off()

pdf(file=paste0(outdir,"in_doublet_cluster.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "in_doublet_cluster")
dev.off()

pdf(file=paste0(outdir,"nUMI.pdf"),width = 10,height = 10)
FeaturePlot(object = seu, features = "nCount_RNA")
dev.off()

pdf(file=paste0(outdir,"perc_above_thresh.pdf"),width = 10,height = 10)
FeaturePlot(object = seu, features = "perc_above_thresh")
dev.off()

seu@meta.data$is_doublet <- seu@meta.data$scDblFinder.standard.class=="doublet" | 
  seu@meta.data$knownDoublets | seu@meta.data$in_doublet_cluster

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

cl <- which.max(c(sum(mat["MLANA",names(which(clus==1))]>0),
                  sum(mat["MLANA",names(which(clus==2))]>0)))

hba1_blacklist <- names(which(clus==cl))

seu@meta.data$hba1_blacklist <- F
seu@meta.data[names(which(clus==cl)),]$hba1_blacklist <- T
seu@meta.data[names(which(clus==cl)),]$is_doublet <- T

pdf(file=paste0(outdir,"HBA1_blacklist.pdf"))
ggplot(seu@meta.data[seu@meta.data$HBA1 > 0,],
       aes(x=nFeature_RNA,y=nCount_RNA, color = hba1_blacklist)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)
dev.off()

min_nGene_erythrocytes <- summary(annot[names(which(clus!=cl)),]$nFeature_RNA)[["1st Qu."]]
min_nUMI_erythrocytes <- summary(annot[names(which(clus!=cl)),]$nCount_RNA)[["1st Qu."]]
summary(annot[names(which(clus==2)),]$percent.top_100)

### Distinguish plasma cells
stopifnot(identical(colnames(seu@assays$RNA@data),colnames(seu@assays$RNA@counts)))
mat <- seu@assays$RNA@data[,seu@assays$RNA@counts["JCHAIN",] > 0 &
                             seu@assays$RNA@counts["IL3RA",] == 0 & 
                             seu@assays$RNA@counts["LILRA4",] == 0]
mat <- mat[which(rowSums(mat) > 0),]
dim(mat)

annot <- as.data.frame(t(mat[c("MLANA","PMEL","CD3D","NCAM1","JCHAIN"),]))
annot <- merge(annot, seu@meta.data[,c("nFeature_RNA","nCount_RNA","is_doublet",
                                       "percent.top_100")],
               all.x=T, all.y=F, by="row.names")
annot$is_doublet <- ifelse(annot$is_doublet,1,0)
rownames(annot) <- annot$Row.names
annot$Row.names <- NULL

pdf(file=paste0(outdir,"pheatmap_JCHAIN.pdf"))
p <- pheatmap(cor(as.matrix(mat),method = "spearman"),
              show_colnames = F,show_rownames = F,annotation_col = annot,
              cutree_cols = 3,cutree_rows = 3)
dev.off()
clus <- cutree(p$tree_col,3)
table(clus)

seu@meta.data$JCHAIN_clusters <- NA
seu@meta.data$JCHAIN_clusters[rownames(seu@meta.data) %in% names(clus)[clus==1]] <- "J1"
seu@meta.data$JCHAIN_clusters[rownames(seu@meta.data) %in% names(clus)[clus==2]] <- "J2"
seu@meta.data$JCHAIN_clusters[rownames(seu@meta.data) %in% names(clus)[clus==3]] <- "J3"

pdf(file=paste0(outdir,"JCHAIN_clusters.pdf"),width = 10,height = 10)
DimPlot(object = seu, group.by = "JCHAIN_clusters")
dev.off()

c(sum(mat["MLANA",names(which(clus==1))]>0),
  sum(mat["MLANA",names(which(clus==2))]>0),
  sum(mat["MLANA",names(which(clus==3))]>0))

table(seu@meta.data$JCHAIN_clusters,seu@meta.data$is_doublet)

JCHAIN_blacklist <- names(which(clus==1))

seu@meta.data$JCHAIN_blacklist <- F
seu@meta.data[names(which(clus==cl)),]$JCHAIN_blacklist <- T
seu@meta.data[names(which(clus==cl)),]$is_doublet <- T

pdf(file=paste0(outdir,"JCHAIN_blacklist.pdf"))
ggplot(seu@meta.data[!is.na(seu@meta.data$JCHAIN_clusters),],
       aes(x=nFeature_RNA,y=nCount_RNA, color = JCHAIN_blacklist)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)
dev.off()

min_nGene_JCHAIN <- summary(annot[names(which(clus!=1)),]$nFeature_RNA)[["1st Qu."]]
min_nUMI_JCHAIN <- summary(annot[names(which(clus!=1)),]$nCount_RNA)[["1st Qu."]]
summary(annot[names(which(clus!=1)),]$percent.top_100)

# Final filtering --------------------------------------------------------------
library(miQC)
library(flexmix)

sce.list <- list()
sce.filtered.list <- list()
model.list <- list()
for (sname in unique(colData(sce)$orig.ident)){
  sce.list[[sname]] <- sce[,colData(sce)$orig.ident==sname]
  model.list[[sname]] <- mixtureModel(sce.list[[sname]])
  
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

sname <- "C8"
model.list[[sname]] <- mixtureModel(sce.list[[sname]],
                                    model_type = "one_dimensional")

pdf(file=paste0(outdir,"miQC_plotModel.",sname,".one_dimensional.pdf"),width = 10,height = 10)
plotModel(sce.list[[sname]], model.list[[sname]])
dev.off()

pdf(file=paste0(outdir,"miQC_plotFiltering.",sname,".one_dimensional.pdf"),width = 10,height = 10)
plotFiltering(sce.list[[sname]], model.list[[sname]],
              posterior_cutoff = 0.9, 
              keep_all_below_boundary = TRUE)
dev.off()

sce.filtered.list[[sname]] <- filterCells(sce.list[[sname]], 
                                          model.list[[sname]], 
                                          posterior_cutoff = 0.9, 
                                          keep_all_below_boundary = TRUE)
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

idx_erythrocyte <- seu@meta.data$HBA1 > 0 & 
  !seu@meta.data$hba1_blacklist & 
  seu@meta.data[["nFeature_RNA"]] > min_nGene_erythrocytes & 
  seu@meta.data[["nCount_RNA"]] > min_nUMI_erythrocytes

idx_plasma_cell <- !is.na(seu@meta.data$JCHAIN_clusters) & 
  !seu@meta.data$JCHAIN_blacklist & 
  seu@meta.data[["nFeature_RNA"]] > min_nGene_JCHAIN & 
  seu@meta.data[["nCount_RNA"]] > min_nUMI_JCHAIN

idx_keep <- ((seu@meta.data[["nFeature_RNA"]] > 200 & 
                seu@meta.data[["nCount_RNA"]] > 500 & 
                seu@meta.data[["percent.top_100"]] < 80) | 
               idx_erythrocyte | idx_plasma_cell) & 
  seu@meta.data$miqc_kept & 
  seu@meta.data[["percent.ribo"]] < 40 & 
  !seu@meta.data$is_doublet & 
  !seu@meta.data$hba1_blacklist & 
  !seu@meta.data$JCHAIN_blacklist

table(seu@meta.data$miqc_kept[idx_erythrocyte])
table(seu@meta.data$miqc_kept[idx_plasma_cell])

table(idx_keep)
table(idx_keep,seu@meta.data$orig.ident)

cells_keep <- rownames(seu@meta.data[idx_keep,])

seu@meta.data$cells_filtered <- ! rownames(seu@meta.data) %in% 
  cells_keep
table(seu@meta.data$cells_filtered)
table(seu@meta.data$cells_filtered)["TRUE"]/nrow(seu@meta.data)
table(seu@meta.data$cells_filtered,
      seu@meta.data$orig.ident)

table(seu@meta.data$cells_filtered[idx_erythrocyte])
table(seu@meta.data$cells_filtered[idx_plasma_cell])

seu@meta.data$erythrocyte <- idx_erythrocyte
seu@meta.data$plasma_cell <- idx_plasma_cell

pdf(file=paste0(outdir,"cells_filtered.pdf"),width = 25,height = 10)
DimPlot(object = seu, group.by = "cells_filtered", split.by = "cells_filtered")
dev.off()

pdf(file=paste0(outdir,"cells_filtered.erythrocyte.pdf"),width = 10,height = 10)
DimPlot(object = subset(seu,erythrocyte==T), group.by = "cells_filtered")
dev.off()

pdf(file=paste0(outdir,"cells_filtered.plasma_cell.pdf"),width = 10,height = 10)
DimPlot(object = subset(seu,plasma_cell==T), group.by = "cells_filtered")
dev.off()

seu.filt <- subset(seu,cells_filtered==F)
saveRDS(seu.filt, paste0(outdir, "seu.filt.rda"))

# FastMNN ----------------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
seu.list <- SplitObject(seu.filt, split.by = "orig.ident")
seu.filt <- RunFastMNN(object.list = seu.list,
                  features = 4000, auto.merge = T, d = 150, k = 15)
seu.filt <- RunUMAP(seu.filt, reduction = "mnn", dims = 1:150)
#saveRDS(seu.filt, paste0(outdir, "seu.filt.rda"))

pdf(file=paste0(outdir,"markers_post_filtering_mnn.pdf"),width = 40,height = 40)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("MLANA","CD3D","NCR1","CD4","CD8A","CD14","HBA1",
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

pdf(file=paste0(outdir,"percent_mt_ribo_top100_post_filtering_mnn.pdf"),width = 30,height = 30)
FeaturePlot(object = seu.filt, reduction = "umap", pt.size = .1,
            features = c("percent.mt","percent.ribo","percent.top_100"))
dev.off()

#seu.filt <- readRDS(paste0(outdir, "seu.filt.rda"))
seu.filt <- FindNeighbors(seu.filt,reduction = "mnn", dims=1:150)
#seu.filt.bak <- seu.filt
seu.filt <- FindClusters(seu.filt,resolution = 2)
#saveRDS(seu.filt,paste0(outdir, "seu.filt.rda"))

pdf(file=paste0(outdir,"clusters.pdf"),width = 20,height = 20)
DimPlot(seu.filt,label = T)
dev.off()

markers.all <- FindAllMarkers(seu.filt)

WriteXLS(markers.all,ExcelFileName = paste0(outdir,"markers.all.res_2.xlsx"),
         BoldHeaderRow = T,AutoFilter = T,AdjWidth = T)

#seu.filt.new <- FindNeighbors(seu.filt.new,reduction = "mnn", dims=1:150)
#seu.filt.new <- FindClusters(seu.filt.new,resolution = 2)
#DimPlot(seu.filt.new,label = T)

# Annotate clusters ------------------------------------------------------------
seu.filt@meta.data$CD3D <- NULL
seu.filt@meta.data$NCAM1 <- NULL
seu.filt@meta.data$HBA1 <- NULL
seu.filt@meta.data$PMEL <- NULL
seu.filt@meta.data$MLANA <- NULL

# macrophage markers: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8785783/
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

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4260088/
#https://www.sciencedirect.com/science/article/pii/S001448001730076X?via%3Dihub

melanoma_clusters <- c(0,1,2,3,4,6,7,13,15,19,21,29,36,37,38)
seu.filt.melanoma <- subset(seu.filt,idents = melanoma_clusters)

markers.melanoma <- FindAllMarkers(seu.filt.melanoma)
WriteXLS(markers.melanoma,ExcelFileName = paste0(outdir,"markers.melanoma.xlsx"),
         BoldHeaderRow = T,AutoFilter = T,AdjWidth = T)

table(markers.melanoma$cluster)
idx <- markers.melanoma$avg_log2FC > 0.35 & 
  markers.melanoma$pct.1 > 0.4 & 
  markers.melanoma$p_val_adj < 0.05
table(markers.melanoma[idx,]$cluster)

m <- markers.melanoma[idx,]
m <- m[order(m$pct.2,decreasing = F),]
m2 <- m[order(m$avg_log2FC,decreasing = F),]
m3 <- m[order(m$avg_log2FC,decreasing = T),]
melanoma_cluster.markers <- list()
for (clus in melanoma_clusters){
  melanoma_cluster.markers[[as.character(clus)]] <- c(
    head(m[m$cluster==clus,]$gene,n=min(c(3,length(m[m$cluster==clus,]$gene)))),
    head(m2[m2$cluster==clus,]$gene,n=min(c(3,length(m2[m2$cluster==clus,]$gene)))),
    head(m3[m3$cluster==clus,]$gene,n=min(c(3,length(m3[m3$cluster==clus,]$gene)))))
}
melanoma_cluster.markers <- unique(unlist(melanoma_cluster.markers))

pdf(file=paste0(outdir,"markers_melanoma.pdf"),width = 20,height = 10)
DotPlot(seu.filt.melanoma,
        features = melanoma_cluster.markers,
        cluster.idents = T, assay = "RNA",scale = F) + vertical_xlabels
dev.off()

#t_nk_clusters <- c(0,2,6,12,13,15,17,18,19,20,21,25,28,29,31)
t_nk_clusters <- c(22,16,28,18,34,8,9,14,25,11,27,12,10,23,20,5,17,31,33,24)
seu.filt.t_nk <- subset(seu.filt,seurat_clusters %in% t_nk_clusters)

DimPlot(seu.filt.t_nk,group.by="seurat_clusters")
#Idents(seu.filt.t_nk) <- seu.filt.t_nk@meta.data$seurat_clusters
markers.t_nk <- FindAllMarkers(seu.filt.t_nk)
WriteXLS(markers.t_nk,ExcelFileName = paste0(outdir,"markers.t_nk.xlsx"),
         BoldHeaderRow = T,AutoFilter = T,AdjWidth = T)

table(markers.t_nk$cluster)
idx <- markers.t_nk$avg_log2FC > 0.35 & 
  markers.t_nk$pct.1 > 0.65 &
  markers.t_nk$pct.2 > 0.3 & 
  markers.t_nk$p_val_adj < 0.05
table(markers.t_nk[idx,]$cluster)

m <- markers.t_nk[idx,]
m <- m[order(m$pct.2,decreasing = F),]
m2 <- m[order(m$avg_log2FC,decreasing = F),]
m3 <- m[order(m$avg_log2FC,decreasing = T),]
t_nk_cluster.markers <- list()
for (clus in t_nk_clusters){
  t_nk_cluster.markers[[as.character(clus)]] <- c(
    head(m[m$cluster==clus,]$gene,n=min(c(3,length(m[m$cluster==clus,]$gene)))),
    head(m2[m2$cluster==clus,]$gene,n=min(c(3,length(m2[m2$cluster==clus,]$gene)))),
    head(m3[m3$cluster==clus,]$gene,n=min(c(3,length(m3[m3$cluster==clus,]$gene)))))
}
t_nk_cluster.markers <- unique(unlist(t_nk_cluster.markers))

DotPlot(seu.filt.t_nk,
        features = t_nk_cluster.markers,
        cluster.idents = T, assay = "RNA",scale = F) + vertical_xlabels

marker_features <- c("PTPRC","CD3D","CD8A","CD8B","CD4","FOXP3",
                     "NCR1","NCAM1",
                     "CD69","PDCD1","TNFRSF9","LAG3","TIGIT","HAVCR2",
                     "IL7R","TCF7","KLRB1",
                     "NKG7","GZMA","CCL5","IL32",
                     "PRF1","TRBC2","CCL4","GZMK","CD27","ITM2A","TRAC",
                     "HLA-DRA","IFNG","ISG15","CD5",
                     "GZMH","GZMB","KLF2","GNLY",
                     "FCRL6","FGFBP2","FCGR3A","CX3CR1",
                     "CD38",
                     "CCL3","CXCL13","CTLA4","KIR2DL4",
                     "TNFSF14","CCR5",
                     "CD160",
                     "XCL2","XCL1","GNG4","EGR2",
                     "KLRC2","CD244","ZNF683",
                     "MKI67")

DotPlot(seu.filt.t_nk,
        features = marker_features,
        cluster.idents = F, assay = "RNA",scale = F) + vertical_xlabels

markers.t_nk <- markers.t_nk[order(markers.t_nk$avg_log2FC,decreasing = T),]
setdiff(markers.t_nk$gene[markers.t_nk$avg_log2FC > 0 & markers.t_nk$cluster==17],
        markers.t_nk$gene[markers.t_nk$avg_log2FC > 0 & markers.t_nk$cluster==31])

setdiff(markers.t_nk$gene[markers.t_nk$avg_log2FC > 0 & markers.t_nk$cluster==31],
        markers.t_nk$gene[markers.t_nk$avg_log2FC > 0 & markers.t_nk$cluster==17])


#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5438859/:

#Idents(seu.filt) <- seu.filt@meta.data$seurat_clusters
# seu.filt <- RenameIdents(seu.filt,
#                          `0` = "CD8+ T cell (Dysfynctional/exhausted)",
#                          `1` = "UM",
#                          `2` = "CD8+ T cell (IL7R+)",
#                          `3` = "UM",
#                          `4` = "UM",
#                          `5` = "UM",
#                          `6` = "CD8+ T cell (IL7R+)",
#                          `7` = "UM",
#                          `8` = "UM",
#                          `9` = "UM",
#                          `10` = "UM",
#                          `11` = "UM",
#                          `12` = "CD4+ T cell",
#                          `13` = "CD4+ T cell",
#                          `14` = "UM",
#                          `15` = "CD8+ T cell (FCGR3A+)",
#                          `16` = "UM (KIT-)",
#                          `17` = "NK cell (I)",
#                          `18` = "NK cell (I)",
#                          `19` = "CD8+ T cell (IL7R+)",
#                          `20` = "Treg (FOXP3+)",
#                          `21` = "CD8+ T cell (Dysfynctional/exhausted)",
#                          `22` = "UM (EGR2+)",
#                          `23` = "Monocyte I (ITGAX+S100A8+FCN1+)",
#                          `24` = "UM",
#                          `25` = "CD8+ T cell (IL7R+)",
#                          `26` = "UM (S100A11-)",
#                          `27` = "Monocyte II",
#                          `28` = "CD8+ T cell (CD69-)",
#                          `29` = "CD8+ T cell (Dysfynctional/exhausted, proliferative)",
#                          `30` = "Endothelial cell",
#                          `31` = "NK cell (II)",
#                          `32` = "UM (CPXM1+VGF+SPP1+)",
#                          `33` = "pDC",
#                          `34` = "Plasma cell",
#                          `35` = "B cell",
#                          `36` = "Erythrocyte")
# seu.filt <- RenameIdents(seu.filt,
#                          `0` = "UM",
#                          `1` = "UM",
#                          `2` = "UM",
#                          `3` = "UM",
#                          `4` = "UM (MXRA8+)",
#                          `5` = "",
#                          `6` = "UM",
#                          `7` = "UM (VGF+)",
#                          `8` = "",
#                          `9` = "",
#                          `10` = "",
#                          `11` = "",
#                          `12` = "",
#                          `13` = "UM (KIT-BCN2-)",
#                          `14` = "",
#                          `15` = "UM",
#                          `16` = "",
#                          `17` = "",
#                          `18` = "",
#                          `19` = "UM (SPARCL1+)",
#                          `20` = "",
#                          `21` = "UM (EGR2+IRF8+)",
#                          `22` = "",
#                          `23` = "",
#                          `24` = "",
#                          `25` = "",
#                          `26` = "Monocyte I (CD1C+ITGAM+CLEC10A+)",
#                          `27` = "",
#                          `28` = "",
#                          `29` = "UM (S100A11-)",
#                          `30` = "Monocyte II (CD1C-ITGAM-CLEC10A-)",
#                          `31` = "",
#                          `32` = "Endothelial cell",
#                          `33` = "",
#                          `34` = "",
#                          `35` = "Monocyte III (CD1C-ITGAM+CLEC10A-)",
#                          `36` = "UM (EGR2+)",
#                          `37` = "UM",
#                          `38` = "UM (VGF+CPXM1+)",
#                          `39` = "pDC",
#                          `40` = "Plasma cell",
#                          `41` = "B cell",
#                          `42` = "Erythrocyte")
Idents(seu.filt) <- seu.filt@meta.data$seurat_clusters
# seu.filt <- RenameIdents(seu.filt,
#                          `0` = "UM",
#                          `1` = "UM",
#                          `2` = "UM",
#                          `3` = "UM",
#                          `4` = "UM",
#                          `5` = "CD8 T cell (Exhausted / dysfunctional)",
#                          `6` = "UM",
#                          `7` = "UM",
#                          `8` = "CD8 T cell (Naive-like)",
#                          `9` = "CD8 T cell (Early activated)",
#                          `10` = "CD8 T cell (Early activated)",
#                          `11` = "CD4 T cell",
#                          `12` = "CD8 T cell (Exhausted / dysfunctional)",
#                          `13` = "UM (KIT-BCN2-)",
#                          `14` = "CD8 T cell (Naive-like)",
#                          `15` = "UM",
#                          `16` = "NK cell (GZMK-CX3CR1+)",
#                          `17` = "Treg (FOXP3+)",
#                          `18` = "NK cell (GZMK+CD160+)",
#                          `19` = "UM",
#                          `20` = "CD8 T cell (Exhausted / dysfunctional)",
#                          `21` = "UM",
#                          `22` = "CD8 T cell (Activated / cytotoxic, FCGR3A+)",
#                          `23` = "CD4 T cell",
#                          `24` = "CD8 T cell (Exhausted / dysfunctional, proliferative)",
#                          `25` = "CD4 T cell",
#                          `26` = "Monocyte I (CD1C+ITGAM+CLEC10A+)",
#                          `27` = "CD8 T cell (KLRB1+)",
#                          `28` = "NK cell (GZMK+CD160-)",
#                          `29` = "UM (S100A11-)",
#                          `30` = "Monocyte II (CD1C-ITGAM-CLEC10A-)",
#                          `31` = "CD8 T cell (Exhausted / dysfunctional)",
#                          `32` = "Endothelial cell",
#                          `33` = "CD8 T cell (CD69-)",
#                          `34` = "NK cell (TRBC2-)",
#                          `35` = "Monocyte III (CD1C-ITGAM+CLEC10A-)",
#                          `36` = "UM",
#                          `37` = "UM",
#                          `38` = "UM (CPXM1+VGF+)",
#                          `39` = "pDC",
#                          `40` = "Plasma cell",
#                          `41` = "B cell",
#                          `42` = "Erythrocyte")

seu.filt <- RenameIdents(seu.filt,
                         `0` = "UM",
                         `1` = "UM",
                         `2` = "UM",
                         `3` = "UM",
                         `4` = "UM",
                         `5` = "CD8 T cell",
                         `6` = "UM",
                         `7` = "UM",
                         `8` = "CD8 T cell",
                         `9` = "CD8 T cell",
                         `10` = "CD8 T cell",
                         `11` = "CD4 T cell",
                         `12` = "CD8 T cell",
                         `13` = "UM (KIT-BCN2-)",
                         `14` = "CD8 T cell",
                         `15` = "UM",
                         `16` = "NK cell (GZMK-CX3CR1+)",
                         `17` = "Treg (FOXP3+)",
                         `18` = "NK cell (GZMK+CD160+)",
                         `19` = "UM",
                         `20` = "CD8 T cell",
                         `21` = "UM",
                         `22` = "CD8 T cell (FCGR3A+)",
                         `23` = "CD4 T cell",
                         `24` = "CD8 T cell (proliferative)",
                         `25` = "CD4 T cell",
                         `26` = "Monocyte I (CD1C+ITGAM+CLEC10A+)",
                         `27` = "CD8 T cell",
                         `28` = "NK cell (GZMK+CD160-)",
                         `29` = "UM (S100A11-)",
                         `30` = "Monocyte II (CD1C-ITGAM-CLEC10A-)",
                         `31` = "CD8 T cell",
                         `32` = "Endothelial cell",
                         `33` = "CD8 T cell",
                         `34` = "NK cell (TRBC2-)",
                         `35` = "Monocyte III (CD1C-ITGAM+CLEC10A-)",
                         `36` = "UM",
                         `37` = "UM",
                         `38` = "UM (CPXM1+VGF+)",
                         `39` = "pDC",
                         `40` = "Plasma cell",
                         `41` = "B cell",
                         `42` = "Erythrocyte")

pdf(file=paste0(outdir,"clusters_annotated.pdf"),width = 19,height = 10)
DimPlot(seu.filt,label = T)
dev.off()

pdf(file=paste0(outdir,"clusters_annotated_by_patient.pdf"),width = 11,height = 10)
DimPlot(seu.filt,label = F,group.by = "orig.ident")
dev.off()

#saveRDS(seu.filt,file=paste0(outdir,"seu.filt.rda"))

#DimPlot(seu.filt,label = T,group.by = "seurat_clusters")


# Sub-clustering ---------------------------------------------------------------

cd8_t_clusters <- sort(as.numeric(as.character(unique(seu.filt@meta.data$seurat_clusters[grep("^CD8",seu.filt@active.ident)]))))
seu.filt.cd8_t <- subset(seu.filt,seurat_clusters %in% cd8_t_clusters)
seu.filt.cd8_t <- subset(seu.filt.cd8_t,CD4 == 0 & NCR1 == 0 & 
                           (CD3D > 0 | CD3G > 0 | CD8A > 0 | CD8B >0) & 
                           MLANA == 0 & CD19 == 0 & CD14 == 0)
seu.filt.cd8_t <- FindNeighbors(seu.filt.cd8_t,reduction = "mnn", k.param = 15, dims=1:150)
#seu.filt.cd8_t <- FindClusters(seu.filt.cd8_t,algorithm = 3, resolution = 1.3)
seu.filt.cd8_t <- FindClusters(seu.filt.cd8_t,algorithm = 1, resolution = 1.3)
seu.filt.cd8_t <- RunUMAP(seu.filt.cd8_t, reduction = "mnn", dims=1:150, 
                          min.dist = 0.05, n.neighbors = 15)

pdf(file=paste0(outdir,"cd8_t_clusters.pdf"),width = 11,height = 10)
DimPlot(seu.filt.cd8_t, label = T)
dev.off()

marker_features_extended <- c("PTPRC","CD3D","CD3G","CD8A","CD8B","FOXP3","IL2RA","FAS",
                     "CD69","PDCD1","CTLA4","TNFRSF9","ICOS","HLA-DRA",
                     "LAG3","TIGIT","HAVCR2",
                     "PRF1","GZMA","GZMB","GZMK","GZMH","GNLY",
                     "FCGR3A","NKG7","KLRB1","KLRC2","KLRD1","KLRG1","KIR2DL4","NCAM1",
                     "IL7R","TCF7",
                     "CCL3","CCL4","CCL5","CCR5","CCR7","CXCL8","CXCL13","CX3CR1",
                     "XCL2","XCL1",
                     "CD38","IFNG","ISG15","CD5",
                     "IL2","IL6","IL32","TRBC2","CD27","CD28","ITM2A","TRAC",
                     "KLF2","FCRL6","FGFBP2",
                     "TNFSF14","CD160",
                     "GNG4","EGR2","CD244","ZNF683",
                     "MKI67","CDKN1A","CDKN2A","ATM","H2AFX","TP53","MAPK14","GLB1","TNF",
                     "LDHA","MAPK8","MAPK1","TERT","TRGV9","TRDV2","CCNE2","CCNA2","CCNB1","CCNB2","CCND1","MDM2")

marker_features <- c("PTPRC","CD3D","CD8A","CD8B",
                              "CD69","PDCD1","CTLA4","TNFRSF9","ICOS","HLA-DRA",
                              "LAG3","TIGIT","HAVCR2",
                              "PRF1","GZMA","GZMB","GZMK","GZMH","GNLY",
                              "FCGR3A","NKG7","KLRB1","KLRC2","KLRG1","KIR2DL4","NCAM1",
                              "IL7R","TCF7","IL2RG","SELL",
                              "CCL3","CCL4","CCL5","CCR5","CCR7","CXCL8","CXCL13","CX3CR1",
                              "XCL2","XCL1","CD44",
                              "CD38","IFNG","ISG15","CD5",
                              "IL2","IL6","IL32","TRBC2",
                              "CD27","CD28","CD70","ITM2A","TRAC",
                              "KLF2","FCRL6","FGFBP2",
                              "TNFSF14","CD160",
                              "GNG4","EGR2","CD244","ZNF683",
                              "MKI67","CDKN1A","ATM","H2AFX","TNF","FOXP1","VSIR",
                     "PECAM1","ITGAE","S1PR1","ITGAL","BTLA","TRAV24","TRAJ18","TRBV11-1")


#markers.cd8_t <- FindAllMarkers(seu.filt.cd8_t)
#WriteXLS(markers.cd8_t,ExcelFileName = paste0(outdir,"markers.cd8_t.xlsx"),
#         BoldHeaderRow = T,AutoFilter = T,AdjWidth = T)

#https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.22351
DotPlot(seu.filt.cd8_t,
        features = union(marker_features,marker_features_extended),
        cluster.idents = T, assay = "RNA",scale = T) + vertical_xlabels


#https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00202-4/fulltext
#https://www.cell.com/trends/immunology/fulltext/S1471-4906(16)30123-5
#https://pubmed.ncbi.nlm.nih.gov/30595452/
#https://www.nature.com/articles/s41423-019-0344-8/tables/1

Idents(seu.filt.cd8_t) <- seu.filt.cd8_t@meta.data$seurat_clusters
markers.cd8_t <- FindAllMarkers(seu.filt.cd8_t)
WriteXLS(markers.cd8_t,ExcelFileName = paste0(outdir,"markers.cd8_t.xlsx"),
         BoldHeaderRow = T,AutoFilter = T,AdjWidth = T)


# seu.filt.cd8_t <- RenameIdents(seu.filt.cd8_t,
#                                `0` = "Exhausted",
#                                `1` = "Late activated",
#                                `2` = "Exhausted",
#                                `3` = "Aactivated",
#                                `4` = "Exhausted",
#                                `5` = "Activated / cytotoxic (FCGR3A+)",
#                                `6` = "Naive-like",
#                                `7` = "Exhausted, proliferative",
#                                `8` = "Activated",
#                                `9` = "Naive-like",
#                                `10` = "Activated, KLRB1+",
#                                `11` = "Naive-like",
#                                `12` = "Activated",
#                                `13` = "Exhausted, TCF7+",
#                                `14` = "Late activated / transitional",
#                                `15` = "CD45 low, CD69- IL7R-")

seu.filt.cd8_t <- RenameIdents(seu.filt.cd8_t,
                               `0` = "Late activated, IL7R+",
                               `1` = "Exhausted",
                               `2` = "Exhausted",
                               `3` = "Activated",
                               `4` = "Exhausted",
                               `5` = "Exhausted",
                               `6` = "Activated, FCGR3A+",
                               `7` = "Naive/memory-like",
                               `8` = "Exhausted, proliferative",
                               `9` = "Activated, CDKN1A+",
                               `10` = "Naive/memory-like",
                               `11` = "Activated, KLRB1+. NCR3+",
                               `12` = "Naive/memory-like, early activated",
                               `13` = "Activated",
                               `14` = "Exhausted, IL7R+",
                               `15` = "Late activated, IL7R-",
                               `16` = "CD45 low, CD69-, IL7R-")

DimPlot(seu.filt.cd8_t,label = T)
DimPlot(seu.filt.cd8_t,label = T,group.by="seurat_clusters")
#FeaturePlot(seu.filt.cd8_t,features = c("CCNA2","CCNE1","CCNB2","MKI67","TUBB","MR1"))

saveRDS(seu.filt.cd8_t,file=paste0(outdir,"seu.filt.cd8_t.rda"))

# Early activated and activated defined based on NR4A1+TNF 
# (see discussion of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522381/; Nur77)
# Essentially, most liver T cells express CD69, irregardless of whether recent activation or not;
# Nur77 is supposedly a marker for recent activation for these instead. 
# "That the expression of CD69 was not promoted by recent TCR cross-linking"
# "in the liver immune activation is not the primary mechanism of HLA-DR/CD38 upregulation"

seu.filt <- AddMetaData(seu.filt,metadata = seu.filt.cd8_t@active.ident,col.name = "CD8_detailed")

saveRDS(seu.filt,file=paste0(outdir,"seu.filt.rda"))

# marker_features_plot <- c("PTPRC","CD3D","CD8A","CD8B",
#                      "CD69","PDCD1","CTLA4","TNFRSF9","ICOS","HLA-DRA",
#                      "LAG3","TIGIT","HAVCR2",
#                      "PRF1","GZMA","GZMB","GZMK","GZMH","GNLY",
#                      "FCGR3A","NKG7","KLRB1","KLRC2","KLRG1","KIR2DL4","NCAM1",
#                      "IL7R","TCF7","IL2RG","SELL",
#                      "CCL3","CCL4","CCL5","CCR5","CCR7","CXCL8","CXCL13","CX3CR1",
#                      "XCL2","XCL1","CD44",
#                      "CD38","IFNG","ISG15","CD5",
#                      "IL2","IL6","IL32","TRBC2",
#                      "CD27","CD28","CD70","ITM2A","TRAC",
#                      "KLF2","FCRL6","FGFBP2",
#                      "TNFSF14","CD160",
#                      "GNG4","EGR2","CD244","ZNF683",
#                      "MKI67","CDKN1A","ATM","H2AFX","TNF","FOXP1","VSIR",
#                      "PECAM1","ITGAE","S1PR1","ITGAL","BTLA")
# 
# markers_features_old_plot <- c("CD3D","CD8A","CD8B","NKG7","GZMA","CCL5",
#                                "IL32","PRF1","TRBC2","CCL4","CD69","GZMK",
#                                "CD27","ITM2A","TRAC","LAG3","TIGIT","HLA-DRA",
#                                "IFNG","PDCD1","ISG15","CD5","GZMH","GZMB",
#                                "KLF2","GNLY","FCRL6","FGFBP2","FCGR3A",
#                                "CX3CR1","HAVCR2","CD38","TNFRSF9","CCL3",
#                                "CXCL13","CTLA4","KIR2DL4","IL7R","KLRB1",
#                                "TCF7","TNFSF14","CCR5","CD160","XCL2","XCL1",
#                                "GNG4","EGR2","KLRC2","CD244","ZNF683","MKI67")
#"RORC"
#DotPlot(seu.filt.cd8_t,
#        features = markers_features_old_plot,
#        cluster.idents = T, assay = "RNA",scale = T) + vertical_xlabels

# markers_plot <- c("CD3D","CD8A","CD8B","NKG7","GZMA","CCL5",
#                   "IL32","PRF1","TRBC2","CCL4","CD69","GZMK",
#                   "CD27","ITM2A","TRAC","LAG3","TIGIT","HLA-DRA",
#                   "IFNG","PDCD1","ISG15","CD5","GZMH","GZMB",
#                   "KLF2","GNLY","FCRL6","FGFBP2","FCGR3A",
#                   "CX3CR1","ENTPD1","TOX","EOMES","HAVCR2","CD38","TNFRSF9","CCL3",
#                   "CXCL13","CTLA4","KIR2DL4","IL7R","TCF7","KLRB1","RORC",
#                   "TNFSF14","CCR5","CD160","XCL2","XCL1",
#                   "GNG4","EGR2","KLRC2","CD244","ZNF683","MKI67","CD28","CDKN1A",
#                   "NCAM1","IL2","IL2RA","JUN","FOS","RNF128","TNF","IL10","NFATC1","NFATC2","NFATC3")

# markers_plot <- c("PTPRC","CD3D","CD8A","CD8B",
#                   "CD69","NR4A1","PDCD1","ENTPD1","TNFRSF9","HAVCR2","CTLA4","KIR2DL4","TOX","EOMES","CD38",
#                   "NKG7","GZMA","CCL5",
#                   "IL32","PRF1","TRBC2","CCL4","GZMK",
#                   "CD27","ITM2A","TRAC","LAG3","TIGIT","HLA-DRA",
#                   "IFNG","ISG15","CD5","GZMH","GZMB",
#                   "KLF2","GNLY","FCRL6","FGFBP2","FCGR3A",
#                   "CX3CR1","CCL3","CXCL13","IL7R","TCF7","SATB1","LEF1","KLRB1","RORC",
#                   "TNFSF14","CCR5","CD160","XCL2","XCL1",
#                   "GNG4","EGR2","KLRC2","CD244","ZNF683","MKI67","CD28","CDKN1A",
#                   "NCAM1","IL2","IL2RA","JUN","FOS","RNF128","TNF","IL10","NFATC1","NFATC2","NFATC3",
#                   "HIF1A","KLRG1","CCR7","BCL2","TBX21","IL17A",
#                   "CXCR6","ITGA1","RCAN2","PDZD8","CRIP1","ITGAE","PRDM1","S1PR1",
#                   "SELL","CD44","CXCR3","RUNX3","ETV5","IL15RA","ITGAL","EPAS1",
#                   "FABP1","FABP4","FABP5","AREG","NCR2","NCR3","TRDV1","LTA",
#                   "TRAJ33","TRAJ12","TRAJ20")

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
                  "CD44","RUNX3","S1PR1","SELL")

# Naive-like: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7767439/
# TRM: https://rupress.org/jem/article/219/4/e20212059/213066/Modulation-of-tissue-resident-memory-T-cells-by
# TRM: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7696520/
# Liver T cells: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522381/

pdf(file=paste0(outdir,"dotplot_cd8_markers.clusters.pdf"), width = 13, height = 5)
DotPlot(seu.filt.cd8_t,
        features =markers_plot,
        cluster.idents = T, assay = "RNA",scale = T,group.by="seurat_clusters") + vertical_xlabels
dev.off()

pdf(file=paste0(outdir,"dotplot_cd8_markers.labels.pdf"), width = 15, height = 4)
DotPlot(seu.filt.cd8_t,
        features =markers_plot,
        cluster.idents = T, assay = "RNA",scale = T) + vertical_xlabels
dev.off()

pdf(file=paste0(outdir,"featureplot_cd8_markers.pdf"), width = 15, height = 10)
FeaturePlot(seu.filt.cd8_t,
            features = c("PDCD1","HAVCR2","TNFRSF9","ICOS",
                         "HLA-DRA","IL7R","TCF7","KLF2",
                         "PRF1","GZMB","KLRB1","CCR7"),ncol = 4)
dev.off()

#seu.filt.cd8_t <- AddMetaData(seu.filt.cd8_t,
#            metadata = seu.filt@meta.data[,
#              c("first.labels.zhao_liver","labels.zhao_liver","pruned.labels.zhao_liver")])



# SingleR ----------------------------------------------------------------------
library(SingleR)
library(celldex)
library(scRNAseq)
library(BiocParallel)

#seu.filt <- readRDS(paste0(outdir,"seu.filt.rda"))

# dice <- DatabaseImmuneCellExpressionData(ensembl = F)
# 
# seu.filt.non_um <- subset(seu.filt,idents=grep("^UM",unique(as.character(seu.filt@active.ident)),value=T,invert = T))
# 
# pred.non_um <- SingleR(test = seu.filt.non_um@assays$RNA@data,
#                       ref = dice, labels = colData(dice)$label.fine)
# 
# stopifnot(identical(rownames(pred.non_um),rownames(seu.filt.non_um@meta.data)))
# 
# seu.filt.non_um <- AddMetaData(seu.filt.non_um,metadata = as.data.frame(pred.non_um[,c("first.labels",
#                                                       "labels","pruned.labels")]))
# 
# DimPlot(seu.filt.non_um,group.by="labels",label = T)
# 
# DotPlot(seu.filt.non_um,
#         features = c(marker_features,c("CD4","FOXP3","CD14")),
#         cluster.idents = T, assay = "RNA",scale = T,group.by = "labels") + vertical_xlabels

zhao_liver <- ZhaoImmuneLiverData(location=F,filter=T)
zhao_liver <- scuttle::logNormCounts(zhao_liver)
rownames(zhao_liver) <- rowData(zhao_liver)$Symbol

pred.non_um.zhao_liver <- SingleR(test = seu.filt.non_um@assays$RNA@data,
                       ref = zhao_liver, labels = colData(zhao_liver)$fine,
                       BPPARAM=MulticoreParam(workers=12))
df <- as.data.frame(pred.non_um.zhao_liver[,c("first.labels","labels","pruned.labels")])
WriteXLS(df, ExcelFileName = paste0(outdir,"pred.non_um.zhao_liver.xlsx"))
saveRDS(pred.non_um.zhao_liver,paste0(outdir,"pred.non_um.zhao_liver.rda"))

stopifnot(identical(rownames(df),rownames(seu.filt.non_um@meta.data)))
colnames(df) <- paste0(colnames(df),".zhao_liver")
 
seu.filt.non_um <- AddMetaData(seu.filt.non_um, metadata = df)

DimPlot(seu.filt.non_um,group.by="labels.zhao_liver",label = T)
DimPlot(seu.filt.non_um,group.by="pruned.labels.zhao_liver",label = T)
DimPlot(seu.filt.non_um,group.by="labels",label = T)
DimPlot(seu.filt.non_um,group.by="seurat_clusters",label = T)

#stoeck <- StoeckiusHashingData(mode='human',ensembl=F,location=F)
#kot <- KotliarovPBMCData(mode="rna",ensembl=F,location=F)
#mair <- MairPBMCData(mode = "rna", ensembl = FALSE, location = F)
bacher <- BacherTCellData(filtered = TRUE, ensembl = FALSE, location = F)
bacher <- scuttle::logNormCounts(bacher)

pred.bacher <- SingleR(test = seu.filt@assays$RNA@data,
                       ref = bacher, labels = colData(bacher)$new_cluster_names,
                       BPPARAM=MulticoreParam(workers=12))
#saveRDS(pred.bacher,paste0(outdir,"pred.bacher.rda"))
#pred.bacher <- readRDS(paste0(outdir,"pred.bacher.rda"))

df <- as.data.frame(pred.bacher[,c("first.labels","labels","pruned.labels")])
stopifnot(identical(rownames(df),rownames(seu.filt@meta.data)))
colnames(df) <- paste0(colnames(df),".bacher")

seu.filt <- AddMetaData(seu.filt, metadata = df)

DimPlot(seu.filt,group.by="labels.bacher",label = T)

pred.zhao_liver <- SingleR(test = seu.filt@assays$RNA@data,
                           ref = zhao_liver, 
                           labels = colData(zhao_liver)$fine,
                           BPPARAM=MulticoreParam(workers=12))
#saveRDS(pred.zhao_liver,paste0(outdir,"pred.zhao_liver.rda"))
pred.zhao_liver <- readRDS(paste0(outdir,"pred.zhao_liver.rda"))

df <- as.data.frame(pred.zhao_liver[,c("first.labels","labels","pruned.labels")])
stopifnot(identical(rownames(df),rownames(seu.filt@meta.data)))
colnames(df) <- paste0(colnames(df),".zhao_liver")
seu.filt <- AddMetaData(seu.filt, metadata = df)

DimPlot(seu.filt,group.by="labels.zhao_liver",label = T)

DimPlot(subset(seu.filt,
               ident = union(grep("T cell",seu.filt@active.ident,value = T),
                             grep("NK cell",seu.filt@active.ident,value = T))),
        group.by="labels.zhao_liver",label = T)

DimPlot(subset(seu.filt,
               ident = union(grep("T cell",seu.filt@active.ident,value = T),
                             grep("NK cell",seu.filt@active.ident,value = T))),
        label = T)

DotPlot(subset(seu.filt,
               ident = union(grep("T cell",seu.filt@active.ident,value = T),
                             grep("NK cell",seu.filt@active.ident,value = T))),
        features = c("CX3CR1","CXCR6","CCR6","NCR1","MR1"))

# infercnv ---------------------------------------------------------------------

library(infercnv)

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
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=paste0(outdir,"infercnv"), 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)



# Old --------------------------------------------------------------------------
lig.reannotated <- readRDS(
  "~/proj/um_ss/Investigations/seurat/results/liger_all/lig.reannotated.rda")

all.plots.reannotated <- plotByDatasetAndCluster(lig.tumor, 
                                                 axis.labels = c("UMAP 1", "UMAP 2"),
                                                 return.plots = T, 
                                                 clusters = lig.reannotated)

vasu_phenotype <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/vasu_phenotype.rda")
all.plots.vasu <- plotByDatasetAndCluster(lig.tumor, 
                                          axis.labels = c('UMAP 1', 'UMAP 2'), 
                                          return.plots = T, 
                                          clusters = vasu_phenotype, 
                                          do.shuffle = F, reorder.idents = T, 
                                          new.order = c("None", "Unidentified", 
                                                        "MART1", "CD3+41bb+"))
all.plots.reannotated[[2]] + all.plots.vasu[[2]] + all.plots.louvain[[2]]

# Investigate components
gene_loadings <- plotGeneLoadings(lig.tumor.louvain, do.spec.plot = T, 
                                  return.plots = TRUE, 
                                  factor.share.thresh = Inf)
setdiff(which(unlist(lapply(gene_loadings, length))>0), dims.remove)


i <- 12
gene_loadings[[i]]
g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]), decreasing = T),]
#head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10)

# all.plots.reannotated[[2]] + gene_loadings[[i]]


# Reannotate clusters

seu <- ligerToSeurat(lig.tumor.louvain)
Idents(seu) <- lig.tumor.louvain@clusters
seu@meta.data$clusters <- Idents(seu) 

# DimPlot(seu,label = T)
# DotPlot(seu,features = c(
#   "CD3D","CD8A","CD4","FOXP3","NCAM1","NCR1","PMEL","MLANA","TYR",
#   "HBA1","HBA2","HBB","JCHAIN","CD19","MS4A1","CD14","TNFRSF21",
#   "ITGAX","MZB1","TNFRSF9",
#   "HAVCR2","ENTPD1","PDCD1","TRGV9","TRDV2"
# ),assay = "RNA",group.by = "clusters") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

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

s <- subset(seu, clusters.annotated != "Mixed")

pdf("~/proj/um_ss/Manuscript/s.dimplot.pdf", width = 10, height = 7)
DimPlot(s, label = T)
dev.off()

s.tnk <- subset(s, clusters.annotated %in% 
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
cells_keep <- setdiff(colnames(s.tnk), cells_remove)

s.tnk <- subset(s.tnk, cells = cells_keep)

cells_keep_liger <- str_split_fixed(cells_keep, "GEX_", 2)[,2]

lig.tnk <- subsetLiger(lig.tumor.louvain, cells.use = cells_keep_liger)
lig.tnk <- runUMAP(lig.tnk, 
                   distance = "cosine", dims.use = dims.use, 
                   n_neighbors = 20, min_dist = 0.05)

s.tnk.reumap <- ligerToSeurat(lig.tnk)

stopifnot(identical(names(Idents(s.tnk.reumap)), names(Idents(s.tnk))))
Idents(s.tnk.reumap) <- Idents(s.tnk)

stopifnot(identical(rownames(s.tnk.reumap@meta.data),
                    rownames(s.tnk@meta.data)))
s.tnk.reumap@meta.data$clusters <- s.tnk@meta.data$clusters
s.tnk.reumap@meta.data$clusters.annotated <- Idents(s.tnk.reumap)
Idents(s.tnk.reumap) <- s.tnk.reumap@meta.data$clusters

cells_keep <- rownames(s.tnk.reumap@meta.data)[
  which(s.tnk.reumap@meta.data$clusters.annotated %in% 
          c("CD8 T cell (activated/exhausted)","CD8 T cell"))]

s.t <- subset(s.tnk, cells = cells_keep)

cells_keep <- str_split_fixed(cells_keep, "GEX_", 2)[,2]
lig.t <- subsetLiger(lig.tnk, cells.use = cells_keep)
lig.t <- runUMAP(lig.t, distance = "cosine", dims.use = dims.use, 
                 n_neighbors = 20, min_dist = 0.05)

s.t.reumap <- ligerToSeurat(lig.t)
stopifnot(identical(names(Idents(s.t.reumap)), names(Idents(s.t))))
Idents(s.t.reumap) <- Idents(s.t)
stopifnot(identical(rownames(s.t.reumap@meta.data),
                    rownames(s.t@meta.data)))
s.t.reumap@meta.data$clusters <- s.t@meta.data$clusters

s.t.reumap@meta.data$clusters.annotated <- Idents(s.t.reumap)
Idents(s.t.reumap) <- s.t.reumap@meta.data$clusters
s.t <- s.t.reumap

marker_features <- c("CD3D","CD8A","CD8B","NKG7","GZMA","CCL5","IL32",
                     "PRF1","TRBC2","CCL4","CD69","GZMK","CD27","ITM2A","TRAC","LAG3",
                     "TIGIT","HLA-DRA","IFNG","PDCD1","ISG15","CD5",
                     "GZMH","GZMB","KLF2","GNLY",
                     "FCRL6","FGFBP2","FCGR3A","CX3CR1",
                     "HAVCR2","CD38","TNFRSF9",
                     "CCL3","CXCL13","CTLA4","KIR2DL4",
                     "IL7R","KLRB1","TCF7","TNFSF14","CCR5",
                     "CD160",
                     "XCL2","XCL1","GNG4","EGR2",
                     "KLRC2","CD244","ZNF683",
                     "MKI67")

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5438859/:
# HLA-DR ... appears at the late stages of activation on T and NK cells; 
# thus, it is considered to be a very late activation marker.
# The earliest activation marker is CD69, which is an inducible cell surface 
# glycoprotein expressed upon activation via the TCR or the IL-2 receptor (CD25)

Idents(s.t) <- s.t@meta.data$clusters
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
                    `28` = "Exhausted / dysfunctional")

levels(s.t) <- c("Naive-like",
                 "Naive-like / early activated",
                 "Early activated (KLRC2+, CD38+)",
                 "Activated, cytotoxic (FCGR3A+)",
                 "Late activated",
                 "Late activated, proliferative",
                 "Exhausted / dysfunctional")

vasu_phenotype.t <- vasu_phenotype[intersect(names(vasu_phenotype), 
                                             str_split_fixed(rownames(s.t@meta.data), "GEX_", 2)[,2])]
stopifnot(identical(names(vasu_phenotype.t), 
                    str_split_fixed(rownames(s.t@meta.data), "GEX_", 2)[,2]))
vasu_phenotype.t <- setNames(vasu_phenotype.t, rownames(s.t@meta.data))
s.t <- AddMetaData(s.t, metadata = vasu_phenotype.t, col.name = "vasu_phenotype")
s.t@meta.data$vasu_phenotype <- factor(s.t@metsa.data$vasu_phenotype,
                                       levels=c("None", "Unidentified", 
                                                "CD3+41bb+", "MART1"))

pdf("~/proj/um_ss/Manuscript/CD8_umap.pdf", width = 9, height = 6)
DimPlot(s.t, label = T)
dev.off()

pdf("~/proj/um_ss/Manuscript/CD8_umap.vasu.pdf", width = 7, height = 6)
DimPlot(s.t,
        label = T, group.by = "vasu_phenotype",
        order = T, cols = c("gray", "gray", "blue","red"))
dev.off()

pdf("~/proj/um_ss/Manuscript/CD8_umap.vasu.v2.pdf", width = 15, height = 6)
DimPlot(subset(s.t,orig.ident %in% c("A2-GEX","C7-GEX","D4-GEX")),
        label = T, group.by = "vasu_phenotype",
        order = T, cols = c("gray", "gray", "blue","red"), 
        split.by = "orig.ident")
dev.off()

pdf("~/proj/um_ss/Manuscript/CD8_umap.vasu.v3.pdf", width = 7, height = 6)
DimPlot(subset(s.t, orig.ident %in% c("A2-GEX","C7-GEX","D4-GEX")),
        label = T, group.by = "vasu_phenotype",
        order = T, cols = c("gray", "gray", "blue","red"))
dev.off()

pdf("~/proj/um_ss/Manuscript/CD8_dotplot.pdf", width = 12, height = 3)
DotPlot(s.t, features = marker_features, cluster.idents = F, assay = "RNA") + 
  vertical_xlabels
dev.off()

pdf("~/proj/um_ss/Manuscript/CD8_tsne_features.pdf", width = 14, height = 12)
FeaturePlot(s.t,features = c("CD69",
                             "TCF7","IL7R",
                             "HLA-DRA","ICOS",
                             "PDCD1","HAVCR2","LAG3","TIGIT"),
            ncol = 3, order = F, slot = "data", keep.scale = "all")
dev.off()

saveRDS(s, file = "~/proj/um_ss/Manuscript/s.rda")
saveRDS(s.t, file = "~/proj/um_ss/Manuscript/s.t.rda")

saveRDS(lig.tumor.louvain, file = "~/proj/um_ss/Manuscript/lig.tumor.louvain.rda")
saveRDS(lig.t, file = "~/proj/um_ss/Manuscript/lig.t.rda")

markers.t <- FindAllMarkers(s.t, assay = "RNA", slot = "data", only.pos = T)

WriteXLS(markers.t,
         ExcelFileName = "~/proj/um_ss/Manuscript/Supplementary_Table_X.CD8_T_markers.xlsx",
         AdjWidth = T, row.names = F, AutoFilter = T, BoldHeaderRow = T)