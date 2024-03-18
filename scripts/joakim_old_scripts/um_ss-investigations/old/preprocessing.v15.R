library(Seurat)
library(tidyverse)
library(rliger)
#library(WriteXLS)
library(scater)
library(djvdj)
library(scDblFinder)

vertical_xlabels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                                     hjust=1))

outdir <- "/data/proj/um_ss/Investigations/seurat/results/v15/"
dir.create(outdir, showWarnings = F, recursive = T)

outdir <- "~/proj/um_ss/Investigations/seurat/results/v15/"
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

lig.tumor <- rliger::createLiger(seu)
lig.tumor <- rliger::normalize(lig.tumor)
lig.tumor <- selectGenes(lig.tumor)
lig.tumor <- scaleNotCenter(lig.tumor)

lig.tumor@cell.data[["project"]] <- ""
lig.tumor@cell.data[["project"]][
  lig.tumor@cell.data$dataset %in% 
    c("A6","A8","B6","B9C1","C3CY","D5","D9","E5")] <- "GWA-JN-605"
lig.tumor@cell.data[["project"]][
  lig.tumor@cell.data$dataset %in% 
    c("A2","B3","C5","C7","C8","D2","D4","D6")] <- "GWA-JN-388"

#saveRDS(lig.tumor,file = paste0(outdir,"/lig.tumor.rda"))

k <- 90
lig.tumor <- optimizeALS(lig.tumor, k = k, lambda = 5, verbose = T)

dims.use <- 1:k
lig.tumor <- quantile_norm(lig.tumor, dims.use = dims.use, 
                           quantiles = 100, eps = 0.01)

lig.tumor <- rliger::runUMAP(lig.tumor, dims.use = dims.use, 
                             n_neighbors = 20, min_dist = 0.05)

#saveRDS(lig.tumor,file = paste0(outdir,"/lig.tumor.rda"))

pdf(file=paste0(outdir,"suggestK_l5.cellbender.pdf"),width = 10,height = 10)
sk <- suggestK(lig.tumor, k.test = seq(10,100,10), lambda = 5, 
               return.data = T, num.cores = 1)
dev.off()

#lig.tumor <- readRDS(paste0(outdir,"/lig.tumor.rda"))

cells_keep <- setdiff(rownames(lig.tumor@cell.data),
                      setdiff(rownames(lig.tumor@cell.data),
                              unique(unlist(lapply(lig.tumor@scale.data,rownames)))))

lig.tumor <- rliger::subsetLiger(lig.tumor, cells.use = cells_keep)

# Add VDJ data -----------------------------------------------------------------
add_vdj_to_liger <- function(lig.tumor, vdj_paths){
  library(djvdj)
  fnames <- vdj_paths
  vdj <- list()
  for (nm in names(fnames)) {
    vdj[[nm]] <- djvdj::import_vdj(input = NULL, vdj_dir = fnames[[nm]],
                                   filter_chains = T, define_clonotypes = "cdr3_gene")
  }
  
  for (nm in names(fnames)) {
    #rownames(vdj[[nm]]) <- paste0(nm,"_",str_split_fixed(rownames(vdj[[nm]]),"-",2)[,1])
    rownames(vdj[[nm]]) <- paste0(nm,"_",rownames(vdj[[nm]]))
  }
  
  vdj <- do.call("rbind",vdj)
  rownames(vdj) <- str_split_fixed(rownames(vdj),"\\.",2)[,2]
  
  orig.rownames <- rownames(lig.tumor@cell.data)
  meta <- merge(lig.tumor@cell.data,vdj,by="row.names",all.x=T)
  rownames(meta) <- meta$Row.names
  meta$Row.names <- NULL
  meta <- meta[orig.rownames,]
  stopifnot(identical(rownames(meta),rownames(lig.tumor@cell.data)))
  lig.tumor@cell.data <- meta
  return(lig.tumor)
}

fnames <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/multi/*/outs/per_sample_outs/*/vdj_t/")
fnames <- fnames[basename(dirname(fnames)) %in% unique(lig.tumor@cell.data$dataset)]
names(fnames) <- basename(dirname(fnames))
fnames <- fnames[names(fnames) %in% lig.tumor@cell.data$dataset]

lig.tumor <- add_vdj_to_liger(lig.tumor, vdj_paths = fnames)

saveRDS(lig.tumor,file = paste0(outdir,"/lig.tumor.rda"))

# Doublets ---------------------------------------------------------------------
library(scDblFinder)
lig.tumor <- readRDS(file = paste0(outdir,"/lig.tumor.rda"))

# The authors do not recommend further QC filtering than this before doublet detection
idx_keep <- lig.tumor@cell.data[["nUMI"]] > 200

cells_keep <- rownames(lig.tumor@cell.data[idx_keep,])

lig.tumor.filt <- rliger::subsetLiger(lig.tumor, cells.use = cells_keep)

seu <- rliger::ligerToSeurat(lig.tumor.filt)
sce <- as.SingleCellExperiment(seu)

library(BiocParallel)
sce.standard <- scDblFinder(sce, samples = "orig.ident", 
                            BPPARAM = MulticoreParam(20))

# Investigate outcome of doublet detection -------------------------------------
mat <- counts(sce)
colnames(mat) <- str_split_fixed(colnames(mat),"_",2)[,2]

stopifnot(identical(colnames(mat),rownames(lig.tumor.filt@cell.data)))
idx <- (!is.na(lig.tumor.filt@cell.data$cdr3) | 
          colSums(mat[c("CD3D","CD4","CD8A","CD8B"),]>0) > 0) & 
  colSums(mat[c("MLANA","PMEL","TYR","HBA1","HBA2","HBB"),]>0) > 0
idx <- idx | ((!is.na(lig.tumor.filt@cell.data$cdr3) & 
                 colSums(mat[c("CD14","CD19","MS4A1","JCHAIN"),]>0) > 0))
idx <- idx | ((colSums(mat[c("MLANA","PMEL","TYR"),]>0) > 0) & 
                (colSums(mat[c("CD14","CD19","MS4A1","JCHAIN","NCR1"),]>0) > 0))
knownDoublets <- names(which(idx))
lig.tumor.filt@cell.data$knownDoublets <- rownames(lig.tumor.filt@cell.data) %in% knownDoublets

knownDoublets <- str_split_fixed(rownames(colData(sce)),"_",2)[,2] %in% knownDoublets
colData(sce) <- cbind(colData(sce),knownDoublets)

table(colData(sce.standard)[["scDblFinder.class"]])

stopifnot(identical(str_split_fixed(rownames(colData(sce.standard)),"_",2)[,2],
                    rownames(lig.tumor.filt@cell.data)))

lig.tumor.filt@cell.data$scDblFinder.standard.score <- 
  colData(sce.standard)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.standard.class <- 
  colData(sce.standard)[["scDblFinder.class"]]

pdf(file=paste0(outdir,"scDblFinder.standard.score.pdf"),width = 10,height = 10)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.standard.score", by.dataset = F)
dev.off()
pdf(file=paste0(outdir,"scDblFinder.standard.class.pdf"),width = 10,height = 10)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.standard.class", by.dataset = F,
                    discrete = T)
dev.off()
pdf(file=paste0(outdir,"knownDoublets.pdf"),width = 10,height = 10)
rliger::plotFeature(object = lig.tumor.filt, feature = "knownDoublets", 
                    by.dataset = F,discrete = T,pt.size = 0.01)
dev.off()

table(doublet_agreement=lig.tumor.filt@cell.data$scDblFinder.standard.class)["doublet"] / 
  nrow(lig.tumor.filt@cell.data)

table(lig.tumor.filt@cell.data$knownDoublets,
      lig.tumor.filt@cell.data$scDblFinder.standard.class)

lig.tumor.filt@cell.data$is_doublet <- lig.tumor.filt@cell.data$scDblFinder.standard.class=="doublet" | 
  lig.tumor.filt@cell.data$knownDoublets

lig.tumor.filt <- louvainCluster(lig.tumor.filt, resolution = 5, k = 5)

#p.louvain <- plotByDatasetAndCluster(lig.tumor.filt, return.plots = T)
#p.louvain[[2]]

doublet_cluster_fisher <- list()
for (clus in unique(lig.tumor.filt@clusters)){
  cells_clus <- names(lig.tumor.filt@clusters)[
    which(lig.tumor.filt@clusters == clus)]
  
  n_is_doublet_clus <- sum(lig.tumor.filt@cell.data[cells_clus,
                                               "is_doublet"])
  n_is_not_doublet_clus <- sum(!lig.tumor.filt@cell.data[cells_clus,
                                                    "is_doublet"])
  
  n_is_doublet_rest <- sum(lig.tumor.filt@cell.data[
    ! rownames(lig.tumor.filt@cell.data) %in% cells_clus,
    "is_doublet"])
  n_is_not_doublet_rest <- sum(!lig.tumor.filt@cell.data[
    ! rownames(lig.tumor.filt@cell.data) %in% cells_clus,
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
for (clus in unique(lig.tumor.filt@clusters)){
  idx <- lig.tumor.filt@clusters == clus
  median_nUMI <- median(lig.tumor.filt@cell.data$nUMI[idx])
  doublet_cluster_fisher[doublet_cluster_fisher$clus == clus,]$median_nUMI <- median_nUMI
  median_nGene <- median(lig.tumor.filt@cell.data$nGene[idx])
  doublet_cluster_fisher[doublet_cluster_fisher$clus == clus,]$median_nGene <- median_nGene
}

thresh <- quantile(lig.tumor.filt@cell.data$nUMI,seq(0,1,0.05))["90%"]

thresh_doublet_perc <- quantile(doublet_cluster_fisher$percent_doublets,seq(0,1,0.05))["95%"]

pdf(file=paste0(outdir,"doublet_cluster_fisher.pdf"),width = 6,height = 6)
plot(doublet_cluster_fisher$percent_doublets,doublet_cluster_fisher$median_nUMI)
abline(h=thresh,v=thresh_doublet_perc)
dev.off()

doublet_clusters <- as.numeric(doublet_cluster_fisher[which(
  doublet_cluster_fisher$q < 0.05 & 
    doublet_cluster_fisher$odds_ratio < 1 & 
    (doublet_cluster_fisher$percent_doublets >= thresh_doublet_perc) | 
    (doublet_cluster_fisher$median_nUMI > thresh)),]$clus)

cells_doublet_clusters <- names(lig.tumor.filt@clusters[
  lig.tumor.filt@clusters %in% doublet_clusters])

lig.tumor.filt@cell.data$in_doublet_cluster <- 
  rownames(lig.tumor.filt@cell.data) %in% cells_doublet_clusters

pdf(file=paste0(outdir,"in_doublet_cluster.pdf"),width = 10,height = 10)
rliger::plotFeature(object = lig.tumor.filt, 
                    feature = "in_doublet_cluster", 
                    by.dataset = F,discrete = T,pt.size = 0.01)
dev.off()
pdf(file=paste0(outdir,"is_doublet.pdf"),width = 10,height = 10)
rliger::plotFeature(object = lig.tumor.filt, 
                    feature = "is_doublet", 
                    by.dataset = F,discrete = T,pt.size = 0.01)
dev.off()

lig.tumor.filt@cell.data$is_doublet <- lig.tumor.filt@cell.data$scDblFinder.standard.class=="doublet" | 
  lig.tumor.filt@cell.data$knownDoublets | lig.tumor.filt@cell.data$in_doublet_cluster

pdf(file=paste0(outdir,"is_doublet_incl_cluster.pdf"),width = 10,height = 10)
rliger::plotFeature(object = lig.tumor.filt, 
                    feature = "is_doublet", 
                    by.dataset = F,discrete = T,pt.size = 0.01)
dev.off()

# Add some QC parameters -------------------------------------------------------
getProportionMito.human <- function (object, use.norm = FALSE) 
{
  all.genes <- Reduce(union, lapply(object@raw.data, rownames))
  mito.genes <- grep(pattern = "^MT-", x = all.genes, value = TRUE)
  data.use <- object@raw.data
  if (use.norm) {
    data.use <- object@norm.data
  }
  percent_mito <- unlist(lapply(unname(data.use), function(x) {
    colSums(x[mito.genes, ])/colSums(x)
  }), use.names = TRUE)
  return(percent_mito)
}

getProportionRibo.human <- function (object, use.norm = FALSE) 
{
  all.genes <- Reduce(union, lapply(object@raw.data, rownames))
  ribo.genes <- grep(pattern = "^RP[SL]", x = all.genes, value = TRUE)
  data.use <- object@raw.data
  if (use.norm) {
    data.use <- object@norm.data
  }
  percent_ribo <- unlist(lapply(unname(data.use), function(x) {
    ribo.genes.isect <- intersect(ribo.genes,rownames(x))
    colSums(x[ribo.genes.isect, ])/colSums(x)
  }), use.names = TRUE)
  return(percent_ribo)
}

lig.tumor.filt@cell.data[["percent_mito"]] <- getProportionMito.human(lig.tumor.filt)
lig.tumor.filt@cell.data[["percent_ribo"]] <- getProportionRibo.human(lig.tumor.filt)

#https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html
seu <- rliger::ligerToSeurat(lig.tumor.filt)
sce <- as.SingleCellExperiment(seu)

all.genes <- Reduce(union, lapply(lig.tumor.filt@raw.data, rownames))
mt.genes <- grep(pattern = "^MT-", x = all.genes, value = TRUE)

sce <- addPerCellQC(sce, subsets = list(mito = mt.genes),percent.top=100)

stopifnot(identical(str_split_fixed(rownames(colData(sce)),"_",2)[,2],
                    rownames(lig.tumor.filt@cell.data)))
lig.tumor.filt@cell.data$percent.top_100 <- colData(sce)[,"percent.top_100"]

# Plots to establish thresholds ------------------------------------------------
rliger::plotFeature(object = lig.tumor.filt,feature = "percent_mito", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt,feature = "percent_ribo", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt,feature = "nGene", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt,feature = "nUMI", by.dataset = F)

lig.tumor.filt@cell.data$percent_ribo_40 <- lig.tumor.filt@cell.data$percent_ribo > 0.4
rliger::plotFeature(object = lig.tumor.filt,feature = "percent_ribo_40", 
                    by.dataset = F,discrete = T)

lig.tumor.filt@cell.data$percent_mito_30 <- lig.tumor.filt@cell.data$percent_mito > 0.3
rliger::plotFeature(object = lig.tumor.filt,feature = "percent_mito_30", 
                    by.dataset = F,discrete = T)

hist(lig.tumor.filt@cell.data$percent_mito)
hist(lig.tumor.filt@cell.data$percent_ribo)

ggplot(lig.tumor.filt@cell.data,aes(x=percent_ribo,y=percent_mito,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels
ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent_mito,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent_mito,color=gene_val > 0)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$gene_val > 0,],
       aes(x=nGene,y=percent_mito)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent_ribo,color=gene_val > 0)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)


ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$gene_val > 0,],
       aes(x=nGene,y=percent_ribo)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$gene_val > 0,],
       aes(x=nGene,y=percent_ribo, color = gene_val)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$gene_val > 0,],
       aes(x=nGene,y=percent_mito, color = gene_val)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)


HBA1 <- getGeneValues(lig.tumor.filt@norm.data,gene = "HBA1")
HBA1_raw <- getGeneValues(lig.tumor.filt@raw.data,gene = "HBA1")
CD3D <- getGeneValues(lig.tumor.filt@norm.data,gene = "CD3D")
MLANA <- getGeneValues(lig.tumor.filt@norm.data,gene = "MLANA")
PMEL <- getGeneValues(lig.tumor.filt@norm.data,gene = "PMEL")
NCAM1 <- getGeneValues(lig.tumor.filt@norm.data,gene = "NCAM1")
stopifnot(identical(names(HBA1),rownames(lig.tumor.filt@cell.data)))
stopifnot(identical(names(CD3D),rownames(lig.tumor.filt@cell.data)))
stopifnot(identical(names(MLANA),rownames(lig.tumor.filt@cell.data)))
stopifnot(identical(names(HBA1_raw),rownames(lig.tumor.filt@cell.data)))
lig.tumor.filt@cell.data <- cbind(lig.tumor.filt@cell.data, HBA1)
lig.tumor.filt@cell.data <- cbind(lig.tumor.filt@cell.data, HBA1_raw)
lig.tumor.filt@cell.data <- cbind(lig.tumor.filt@cell.data, CD3D)
lig.tumor.filt@cell.data <- cbind(lig.tumor.filt@cell.data, MLANA)
lig.tumor.filt@cell.data <- cbind(lig.tumor.filt@cell.data, PMEL)
lig.tumor.filt@cell.data <- cbind(lig.tumor.filt@cell.data, NCAM1)

ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$HBA1 > 0 & 
                                  lig.tumor.filt@cell.data$CD3D == 0 & 
                                  lig.tumor.filt@cell.data$MLANA == 0,],
       aes(x=nGene,y=nUMI, color = is_doublet)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)


### Distinguish erythrocytes
mat <- seu@assays$RNA@data[,seu@assays$RNA@data["HBA1",] > 0]
mat <- mat[which(rowSums(mat) > 0),]
dim(mat)

library(pheatmap)
colnames(mat) <- str_split_fixed(colnames(mat),"_",2)[,2]
annot <- as.data.frame(t(mat[c("MLANA","PMEL","CD3D","NCAM1","HBA1"),]))
annot <- merge(annot, lig.tumor.filt@cell.data[,c("nGene","nUMI","is_doublet","percent.top_100")],
      all.x=T, all.y=F, by="row.names")
annot$is_doublet <- ifelse(annot$is_doublet,1,0)
rownames(annot) <- annot$Row.names
annot$Row.names <- NULL
pdf(file=paste0(outdir,"pheatmap_HBA1.pdf"))
p <- pheatmap(cor(as.matrix(mat),method = "spearman"),show_colnames = F,show_rownames = F,annotation_col = annot)
dev.off()
clus <- cutree(p$tree_col,2)
table(clus)

cl <- which.max(c(sum(mat["MLANA",names(which(clus==1))]>0),
                  sum(mat["MLANA",names(which(clus==2))]>0)))

hba1_blacklist <- names(which(clus==cl))

lig.tumor.filt@cell.data$hba1_blacklist <- F
lig.tumor.filt@cell.data[names(which(clus==cl)),]$hba1_blacklist <- T
lig.tumor.filt@cell.data[names(which(clus==cl)),]$is_doublet <- T

pdf(file=paste0(outdir,"HBA1_blacklist.pdf"))
ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$HBA1 > 0,],
       aes(x=nGene,y=nUMI, color = hba1_blacklist)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)
dev.off()

min_nGene_erythrocytes <- summary(annot[names(which(clus==2)),]$nGene)[["1st Qu."]]
min_nUMI_erythrocytes <- summary(annot[names(which(clus==2)),]$nUMI)[["1st Qu."]]
summary(annot[names(which(clus==2)),]$percent.top_100)

ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$HBA1 > 0 & !lig.tumor.filt@cell.data$hba1_blacklist,],
       aes(x=nGene,y=nUMI, color = hba1_blacklist)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500) + 
  geom_vline(xintercept = min_nGene_erythrocytes) + geom_hline(yintercept = min_nUMI_erythrocytes)

ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$HBA1 > 0 & !lig.tumor.filt@cell.data$hba1_blacklist,],
       aes(x=nGene,y=nUMI, color = is_doublet)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500) + 
  geom_vline(xintercept = min_nGene_erythrocytes) + geom_hline(yintercept = min_nUMI_erythrocytes)

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent_ribo,color=percent_mito >= 0.3)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

df_perc_mito <- cor(lig.tumor.filt@H.norm,lig.tumor.filt@cell.data$percent_mito)
df_perc_mito <- data.frame(fact=1:90, cr=df_perc_mito)
df_perc_mito <- df_perc_mito[order(df_perc_mito$cr,decreasing = T),]

df_perc_ribo <- cor(lig.tumor.filt@H.norm,lig.tumor.filt@cell.data$percent_ribo)
df_perc_ribo <- data.frame(fact=1:90, cr=df_perc_ribo)
df_perc_ribo <- df_perc_ribo[order(df_perc_ribo$cr,decreasing = T),]

df_mito_clusters <- data.frame(percent_mito = lig.tumor.filt@cell.data$percent_mito,
                               cluster=lig.tumor.filt@clusters)

df_mito_clusters_mn <- aggregate(df_mito_clusters$percent_mito,
                                 by=list(df_mito_clusters$cluster),FUN=mean)
colnames(df_mito_clusters_mn) <- c("cluster","mean")
df_mito_clusters_sd <- aggregate(df_mito_clusters$percent_mito,
                                 by=list(df_mito_clusters$cluster),FUN=sd)
colnames(df_mito_clusters_sd) <- c("cluster","sd")
df_mito_clusters_aggr <- merge(df_mito_clusters_mn,df_mito_clusters_sd,
                               by="cluster")

df_ribo_clusters <- data.frame(percent_ribo = lig.tumor.filt@cell.data$percent_ribo,
                               cluster=lig.tumor.filt@clusters)

df_ribo_clusters_mn <- aggregate(df_ribo_clusters$percent_ribo,
                                 by=list(df_ribo_clusters$cluster),FUN=mean)
colnames(df_ribo_clusters_mn) <- c("cluster","mean")
df_ribo_clusters_sd <- aggregate(df_ribo_clusters$percent_ribo,
                                 by=list(df_ribo_clusters$cluster),FUN=sd)
colnames(df_ribo_clusters_sd) <- c("cluster","sd")
df_ribo_clusters_aggr <- merge(df_ribo_clusters_mn,df_ribo_clusters_sd,
                               by="cluster")

ggplot(df_mito_clusters_aggr,aes(x=reorder(cluster,-mean),y=mean)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  theme_classic() + vertical_xlabels

ggplot(df_ribo_clusters_aggr,aes(x=reorder(cluster,-mean),y=mean)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  theme_classic() + vertical_xlabels

rliger::plotFeature(object = lig.tumor.filt, feature = "percent.top_100", by.dataset = F)

hist(lig.tumor.filt@cell.data$percent.top_100) # percentage of library size counts coming from the top 100 overall genes

ggplot(lig.tumor.filt@cell.data,aes(x=percent_mito,y=percent.top_100,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels

ggplot(lig.tumor.filt@cell.data,aes(x=percent_ribo,y=percent.top_100,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels

ggplot(lig.tumor.filt@cell.data,aes(x=percent_mito,y=percent.top_100,color=percent_ribo)) + 
  geom_point() + theme_classic() + vertical_xlabels

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent.top_100,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent.top_100,color=gene_val)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=nUMI,color=gene_val)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=nUMI,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data,aes(x=percent_mito,y=percent_ribo,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 0.3) + 
  geom_hline(yintercept = 0.4) + 
  facet_wrap(~ dataset)

# Final filtering --------------------------------------------------------------
library(miQC)
library(flexmix)
seu <- rliger::ligerToSeurat(lig.tumor)
sce <- as.SingleCellExperiment(seu)

all.genes <- Reduce(union, lapply(lig.tumor@raw.data, rownames))
mt.genes <- grep(pattern = "^MT-", x = all.genes, value = TRUE)
ribo.genes <- grep(pattern = "^RP[SL]", x = all.genes, value = TRUE)

sce <- addPerCellQC(sce, subsets = list(mito = mt.genes,
                                        ribo = ribo.genes),
                    percent.top=100)
#sce.bak <- sce
#sce <- sce[,colData(sce)$subsets_ribo_percent < 0.4]

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

model.list[[sname]] <- mixtureModel(sce.list[[sname]],
                                    model_type = "spline")

pdf(file=paste0(outdir,"miQC_plotModel.",sname,".spline.pdf"),width = 10,height = 10)
plotModel(sce.list[[sname]], model.list[[sname]])
dev.off()

pdf(file=paste0(outdir,"miQC_plotFiltering.",sname,".spline.pdf"),width = 10,height = 10)
plotFiltering(sce.list[[sname]], model.list[[sname]],
              posterior_cutoff = 0.9, 
              keep_all_below_boundary = TRUE)
dev.off()

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
  miqc.list[[sname]] <- str_split_fixed(rownames(colData(
    sce.filtered.list[[sname]])),"_",2)[,2]
}
miqc <- as.character(unlist(miqc.list))

lig.tumor@cell.data$miqc_kept <- rownames(lig.tumor@cell.data) %in% miqc

pdf(file=paste0(outdir,"miQC_kept.pdf"),width = 10,height = 10)
rliger::plotFeature(object = lig.tumor,feature = "miqc_kept", 
                    by.dataset = F,discrete = T)
dev.off()

#lig.tumor.fixed@cell.data$nUMI_ge_2000 <- lig.tumor.fixed@cell.data$nUMI >= 500
#rliger::plotFeature(object = lig.tumor.fixed,feature = "nUMI_ge_2000", 
#                    by.dataset = F,discrete = T)


lig.tumor.filt@cell.data$miqc_kept <- rownames(lig.tumor.filt@cell.data) %in% 
  rownames(lig.tumor@cell.data)[lig.tumor@cell.data$miqc_kept]

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent_mito,color=miqc_kept)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500) + 
  geom_hline(yintercept = 0.3) + 
  facet_wrap(~ dataset)
ggplot(lig.tumor.filt@cell.data,aes(x=nUMI,y=percent_mito,color=miqc_kept)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 2000) + 
  geom_hline(yintercept = 0.3) + 
  facet_wrap(~ dataset)

ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$miqc_kept & 
                                  !lig.tumor.filt@cell.data$is_doublet & 
                                  lig.tumor.filt@cell.data$percent.top_100 < 80 & 
                                  lig.tumor.filt@cell.data$percent_ribo < 40,], 
       aes(color=dataset, x=nUMI, fill= dataset)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

idx_erythrocyte <- lig.tumor.filt@cell.data$HBA1 > 0 & 
  !lig.tumor.filt@cell.data$hba1_blacklist & 
  lig.tumor.filt@cell.data[["nGene"]] > min_nGene_erythrocytes & 
  lig.tumor.filt@cell.data[["nUMI"]] > min_nUMI_erythrocytes

# idx_keep <- ((lig.tumor.filt@cell.data[["nGene"]] > 500 & 
#   lig.tumor.filt@cell.data[["nUMI"]] > 500 & 
#   lig.tumor.filt@cell.data[["percent.top_100"]] < 80) | idx_erythrocyte) & 
#   lig.tumor.filt@cell.data[["percent_mito"]] < 0.3 & 
#   lig.tumor.filt@cell.data[["percent_ribo"]] < 0.4 & 
#   !lig.tumor.filt@cell.data$is_doublet

idx_keep <- ((lig.tumor.filt@cell.data[["nGene"]] > 500 & 
                lig.tumor.filt@cell.data[["nUMI"]] > 500 & 
                lig.tumor.filt@cell.data[["percent.top_100"]] < 80) | idx_erythrocyte) & 
  lig.tumor.filt@cell.data$miqc_kept & 
  lig.tumor.filt@cell.data[["percent_ribo"]] < 0.4 & 
  !lig.tumor.filt@cell.data$is_doublet & 
  !lig.tumor.filt@cell.data$hba1_blacklist

table(lig.tumor.filt@cell.data$miqc_kept[idx_erythrocyte])

table(idx_keep)
table(idx_keep,lig.tumor.filt@cell.data$dataset)

cells_keep <- rownames(lig.tumor.filt@cell.data[idx_keep,])

lig.tumor.filt@cell.data$cells_filtered <- ! rownames(lig.tumor.filt@cell.data) %in% 
  cells_keep
table(lig.tumor.filt@cell.data$cells_filtered)
table(lig.tumor.filt@cell.data$cells_filtered)["TRUE"]/nrow(lig.tumor.filt@cell.data)
table(lig.tumor.filt@cell.data$cells_filtered,
      lig.tumor.filt@cell.data$dataset)

pdf(file=paste0(outdir,"cells_filtered.pdf"),width = 10,height = 10)
rliger::plotFeature(object = lig.tumor.filt, feature = "cells_filtered", 
                    by.dataset = F, discrete = T)
dev.off()

lig.tumor.filt <- rliger::subsetLiger(lig.tumor.filt, cells.use = cells_keep)
saveRDS(lig.tumor.filt, paste0(outdir, "lig.tumor.filt.rda"))

# Reprocess filtered dataset ---------------------------------------------------

#lig.tumor.filt <- readRDS(paste0(outdir, "lig.tumor.filt.rda"))

k <- 60
lig.tumor.filt <- optimizeALS(lig.tumor.filt, k = k, lambda = 5, verbose = T)

pdf(file=paste0(outdir,"suggestK_l5.filt.pdf"),width = 10,height = 10)
sk <- suggestK(lig.tumor.filt, k.test = seq(10,60,5), lambda = 5, 
               return.data = T, num.cores = 1)
dev.off()

#saveRDS(lig.tumor.filt, paste0(outdir, "lig.tumor.filt.rda"))

lig.tumor.filt_k50 <- optimizeALS(lig.tumor.filt, k = 50, lambda = 5, verbose = T)
#saveRDS(lig.tumor.filt_k50, paste0(outdir, "lig.tumor.filt_k50.rda"))

#lig.tumor.filt <- readRDS(paste0(outdir, "lig.tumor.filt.rda"))
table(lig.tumor.filt@cell.data$hba1_blacklist)

k <- 60
dims.use <- 1:k
lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, 
                                quantiles = 100, eps = 0.01)
lig.tumor.filt <- rliger::runUMAP(lig.tumor.filt, dims.use = dims.use, 
                                  n_neighbors = 20, min_dist = 0.05)

#saveRDS(lig.tumor.filt, paste0(outdir, "lig.tumor.filt.rda"))

# Select factors, reprocess, Annotate clusters ---------------------------------
#lig.tumor.filt <- readRDS(paste0(outdir, "lig.tumor.filt.rda"))

p <- plotByDatasetAndCluster(lig.tumor.filt, return.plots = T)
p[[1]] + p[[2]]

#pw <- fgsea::gmtPathways(gmt.file = "~/proj/um_ss/Investigations/data/MSigDB/c5.go.v2022.1.Hs.symbols.gmt")
#pw <- fgsea::gmtPathways(gmt.file = "~/proj/um_ss/Investigations/data/MSigDB/c5.go.v2022.1.Hs.entrez.gmt")
#gsea_output <- rliger::runGSEA(lig.tumor.filt, mat_v = 1:16)
gsea_output <- rliger::runGSEA(lig.tumor.filt)

gene_loadings <- plotGeneLoadings(lig.tumor.filt, do.spec.plot = T,
                                  return.plots = TRUE, factor.share.thresh = Inf)
#for (i in c(47,21,51,54,19)) {
#for (i in c(7,17,40,41,50,52,56)) {
#for (i in c(7,17,40,41,50)) {
# 9
# 10: Sex?
# 19
# 
# 33
# 
# 44
# 
# #55
# 53
# 57

20
56

#Weird outliers: 16, 24, 26

#General cell stress?: 19, 54, 3?, 5, 9?

df <- data.frame(comp=1:60,
                 cr_mito=cor(lig.tumor.filt@H.norm,lig.tumor.filt@cell.data$percent_mito),
                 cr_ribo=cor(lig.tumor.filt@H.norm,lig.tumor.filt@cell.data$percent_ribo),
                 cr_100=cor(lig.tumor.filt@H.norm,lig.tumor.filt@cell.data$percent.top_100),
                 cr_nUMI=cor(lig.tumor.filt@H.norm,lig.tumor.filt@cell.data$nUMI),
                 cr_nGene=cor(lig.tumor.filt@H.norm,lig.tumor.filt@cell.data$nGene))
rownames(df) <- df$comp
df$comp <- NULL

library(pheatmap)
pheatmap(df,border_color = NA)

df[order(df$cr_mito,decreasing = T),]
df[order(df$cr_ribo,decreasing = T),]

c(3,45,10,23)
c(21,40,25,54,9,30)
37

42
43

45

53

x <- lig.tumor.filt@H.norm[lig.tumor.filt@cell.data$HBA1_raw>0 | lig.tumor.filt@cell.data$CD3D >0,]
df2 <- data.frame(comp=1:ncol(x),cv=apply(x,2,sd)/apply(x,2,mean))
df2[order(df2$cv),]
head(df2[order(df2$cv),],n=10)

x <- lig.tumor.filt@H.norm
df2 <- data.frame(comp=1:ncol(x),cv=apply(x,2,sd)/apply(x,2,mean))
df2[order(df2$cv),]
head(df2[order(df2$cv),],n=20)

intersect(dims.use,head(df2[order(df2$cv),],n=10)$comp)


#for (i in intersect(dims.use,head(df2[order(df2$cv),],n=10)$comp)){
for (i in dims.use){
  #  i <- 10
  print(i)
  g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
  print(head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10))
  print(gene_loadings[[i]])
  i <- invisible(readline(prompt="Press [enter] to continue"))
}
#51
#55
dims.use <- setdiff(1:60,c(20)) # cell cycle
dims.use <- setdiff(dims.use,c(56)) # ribo
dims.use <- setdiff(dims.use,c(23)) # mitochondrial related, very high in only one or two samples (all cells roughly in C3CY and D5)
dims.use <- setdiff(dims.use,c(19,54)) # stress, general
dims.use <- setdiff(dims.use,c(5)) # rna processing, general
dims.use <- setdiff(dims.use,c(16,24)) # weird outliers, uninteresting categories
dims.use <- setdiff(dims.use,c(53,40,50,34,3,45,29,33,36,44,55,2)) # test
#dims.use <- setdiff(dims.use,c(56,52,33,9)) # mito + plotFactor general expr
#dims.use <- setdiff(dims.use,c(16,60)) # weird outliers
lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, 
                                quantiles = 100, eps = 0.01)
lig.tumor.filt <- rliger::runUMAP(lig.tumor.filt, dims.use = dims.use, 
                                  n_neighbors = 15, min_dist = 0.05)

p <- plotByDatasetAndCluster(lig.tumor.filt, return.plots = T)
p[[1]] + p[[2]]


library(msigdbr)
m_set <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
m_set <- m_set[grepl("MITOCHON", m_set$gs_name),]
m_set <- split(m_set$entrez_gene, f = m_set$gs_name)
gsea_output_mito <- runGSEA(lig.tumor.filt, custom_gene_sets = m_set)

for (i in 1:length(gsea_output_mito)){
  g <- gsea_output_mito[[i]][order(unlist(gsea_output_mito[[i]][,5]),decreasing = T),]
  g <- as.data.frame(g)
  g$leadingEdge <- NULL
  gsea_output_mito[[i]] <- g[which(g$padj < 0.05 & g$NES > 0),]
}
n_mito_sets <- data.frame(fact=1:60,n_mito_sets=unlist(lapply(gsea_output_mito,nrow)))
n_mito_sets$fact
n_mito_sets$n_mito_sets
head(n_mito_sets)

ggplot(n_mito_sets,aes(x=reorder(fact,-n_mito_sets),y=n_mito_sets)) +
  geom_bar(stat="identity",position="dodge") +
  theme_classic() +
  vertical_xlabels
# 
# n_mito_sets <- n_mito_sets[order(n_mito_sets$n_mito_sets,decreasing = T),]
# n_mito_sets[n_mito_sets$n_mito_sets>0,]
# 
# gsea_output_reactome <- list()
# for (i in 1:length(gsea_output)){
#   g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
#   g <- as.data.frame(g)
#   g$leadingEdge <- NULL
#   gsea_output_reactome[[i]] <- g[which(g$padj < 0.05 & g$NES > 0),]
# }
# 
# i <- 81
# gsea_output_reactome[[i]]
# gsea_output_reactome[[i]][order(unlist(gsea_output_reactome[[i]]$pval),decreasing = F),]
# gene_loadings[[i]]
# n_mito_sets[n_mito_sets$n_mito_sets>0 & n_mito_sets$fact==i,]
# 
# n_mito_sets[n_mito_sets$fact %in% sex_p[1:7,]$fact,]
# 
# rliger::plotGene(lig.tumor.filt,"XIST",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

pdf("~/proj/um_ss/Investigations/seurat/results/v15/plotFactors.pdf",width = 20)
plotFactors(lig.tumor.filt)
dev.off()

# table(lig.tumor.filt@cell.data[["project"]],lig.tumor.filt@cell.data$dataset)

get_gene_loadings <- function(liger_object,
                              liger_factor){
  W <- t(liger_object@W)
  rownames(W) <- colnames(liger_object@scale.data[[1]])
  gene_loadings <- W[order(W[, liger_factor], decreasing = TRUE),liger_factor]
  
  return(gene_loadings)
}

# gene_val <- getGeneValues(lig.tumor.filt@norm.data, "XIST")
# stopifnot(identical(names(gene_val),rownames(lig.tumor.filt@H.norm)))
# 
# cr <- cor(gene_val, lig.tumor.filt@H.norm, method = "spearman")
# cr <- as.data.frame(t(cr))
# colnames(cr) <- "correlation"
# cr$fact <- 1:nrow(cr)
# head(cr)
# 
# cr <- cr[order(abs(cr$correlation),decreasing = T),]
# 
# ggplot(cr,aes(x=reorder(fact,correlation),y=correlation)) + geom_point() + 
#   theme_classic() + vertical_xlabels
# 
# ggplot(cr,aes(x=reorder(fact,correlation),y=abs(correlation))) + geom_point() + 
#   theme_classic() + vertical_xlabels + geom_abline(intercept = 0.1,slope = 0) + 
#   xlab(NULL)

# find celltype-unspecific factors
plot_factor_vs_marker <- function(lig.tumor.filt, marker, do_plot=T){
  gene_val <- getGeneValues(lig.tumor.filt@norm.data,gene = marker)
  
  x <- lig.tumor.filt@H.norm[which(gene_val>0),]
  colnames(x) <- 1:ncol(x)
  
  y <- lig.tumor.filt@H.norm[which(gene_val==0),]
  colnames(y) <- 1:ncol(y)
  
  df <- data.frame(perc_ct=colSums(x>0)/nrow(x),
                   perc_non_ct=colSums(y>0)/nrow(y),
                   fact=1:ncol(x))
  df$ratio=df$perc_ct/df$perc_non_ct
  
  if (do_plot){
    print(ggplot(df,aes(x=reorder(fact,-ratio),y=ratio)) + geom_point() + 
            theme_classic() + vertical_xlabels + xlab(NULL) + 
            geom_bar(aes(x=reorder(fact,-ratio),y=perc_ct),stat="identity",position="dodge",fill="green") + 
            geom_bar(aes(x=reorder(fact,-ratio),y=perc_non_ct),stat="identity",position="dodge",fill="red") + 
            geom_abline(intercept = 1,slope = 0))
  }
  return(df)
}

markers <- c("MLANA","CD3D","HBA1","HBB","NCR1","CD14",
             "JCHAIN","MS4A1","CD19","CD4","CD8A","CD34","MKI67")
gene_val <- list()
for (marker in markers){
  print(marker)
  gene_val[[marker]] <- plot_factor_vs_marker(lig.tumor.filt, marker, do_plot=F)
}

#df <- do.call("cbind",lapply(gene_val,function(x){as.numeric(x$ratio)}))
df <- do.call("cbind",lapply(gene_val,function(x){as.numeric(x$perc_ct)}))
df <- as.data.frame(df)
rownames(df) <- 1:nrow(df)

pheatmap(df,border_color = F,cluster_rows = F,cluster_cols = F)
# pheatmap(df,border_color = F)
# 
# 
# pheatmap(df[order(apply(df,1,mad),decreasing = T),],
#          cluster_cols = T,cluster_rows = F,border_color = NA)
# pheatmap(df[order(apply(df,1,sd),decreasing = T),],
#          cluster_cols = T,cluster_rows = F,border_color = NA)
# pheatmap(df[order(apply(df,1,function(x){
#   v <- sort(x,decreasing = T);
#   v <- v[1]-v[2];
#   return(v)}),decreasing = T),],
#   cluster_cols = T,cluster_rows = F,border_color = NA)
# 
# apply(df,1,function(x){
#   v <- max(abs(x-median(as.numeric(x))));
#   return(v)})
# 
# pheatmap(df[order(apply(df,1,function(x){
#   v <- max(abs(x-median(as.numeric(x))));
#   return(v)}),decreasing = T),],
#   cluster_cols = T,cluster_rows = F,border_color = NA)
# 
# pheatmap(df[apply(df,1,function(x){
#   v <- max(abs(x-median(as.numeric(x))));
#   return(v)})>1,],
#   cluster_cols = T,cluster_rows = T,border_color = NA)
# 
# df2 <- df
# df2[df2<1] <- 0
# pheatmap(df2[order(rowSums(df2==0),decreasing = T),],
#          cluster_cols = T,cluster_rows = F,border_color = NA)

pheatmap(df[order(apply(df,1,max),decreasing = T),],
         cluster_rows = F,cluster_cols = T,border_color = NA)

df <- df[order(apply(df,1,max),decreasing = T),colnames(df)!="MKI67"]
#pheatmap(df[apply(df,1,max) > summary(sort(apply(df,1,max)))[["1st Qu."]],],
#         cluster_rows = F,cluster_cols = T,border_color = NA)
pheatmap(df[apply(df,1,max) > summary(sort(apply(df,1,max)))[["1st Qu."]],],
         cluster_rows = F,cluster_cols = T,border_color = NA)

dims.use <- sort(as.numeric(rownames(df[apply(df,1,max) > summary(sort(apply(df,1,max)))[["Median"]],])))


hba1 <- plot_factor_vs_marker(lig.tumor.filt = lig.tumor.filt, marker = "HBA1", do_plot = T)
hba2 <- plot_factor_vs_marker(lig.tumor.filt = lig.tumor.filt, marker = "HBA2", do_plot = T)
hbb <- plot_factor_vs_marker(lig.tumor.filt = lig.tumor.filt, marker = "HBB", do_plot = T)
mki67 <- plot_factor_vs_marker(lig.tumor.filt = lig.tumor.filt, marker = "MKI67", do_plot = T)

#dims.remove <- c(c(16,58,59,87),4,88,18,8)
dims.remove <- c(49,55)
dims.use <- setdiff(dims.use,dims.remove)
lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, 
                                quantiles = 100, eps = 0.01)
lig.tumor.filt <- rliger::runUMAP(lig.tumor.filt, dims.use = dims.use, 
                                  n_neighbors = 10, min_dist = 0.05) #distance = "cosine", 

p <- plotByDatasetAndCluster(lig.tumor.filt, return.plots = T)
p[[1]] + p[[2]]

rliger::plotGene(lig.tumor.filt,"TRBV30",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"TRGV10",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"TNFRSF9",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"HAVCR2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"PDCD1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"ENTPD1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"LAG3",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"TRGV9",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"TRDV2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"XIST",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"MKI67",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"CD3D",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CD4",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CD8A",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CD8B",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"FOXP3",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"NCAM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"NCR1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"PMEL",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"MLANA",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"TYR",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CPXM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"HBA1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"HBA2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"HBB",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"JCHAIN",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CD19",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"MS4A1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"KRT14",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"CD14",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CD68",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"ITGAM",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"ITGAX",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"ITGA4",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"FCGR3A",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"FCGR3B",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"CD69",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"ACTA2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"ALB",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"PDGFRA",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"PECAM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CD34",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotFeature(object = lig.tumor.filt, feature = "percent_ribo", by.dataset = F, discrete = T)
rliger::plotFeature(object = lig.tumor.filt, feature = "percent_mito", by.dataset = F, discrete = T)
rliger::plotFeature(object = lig.tumor.filt, feature = "percent.top_100", by.dataset = F, discrete = T)
rliger::plotFeature(object = lig.tumor.filt, feature = "nUMI", by.dataset = F, discrete = T)
rliger::plotFeature(object = lig.tumor.filt, feature = "hba1_blacklist", 
                    by.dataset = F, discrete = T,do.shuffle = F)



seu <- rliger::ligerToSeurat(lig.tumor.filt)

DotPlot(seu, 
        features = c("CD3D","CD3G","CD8A","CD4",
                     "NCAM1","NCR1",
                     "CD19","MS4A1",
                     "CD14","TNFRSF21","JCHAIN","ITGAX",
                     "ITGAM",
                     "MNDA","ELANE","CD24","FCGR3B","FCGR3A",
                     "MZB1","CD68",
                     "HBA1","HBA2","HBB",
                     "PMEL","MLANA","TYR","MKI67"),
        cluster.idents = F, assay = "RNA") + vertical_xlabels


# lig.tumor <- louvainCluster(lig.tumor, resolution = 1, dims.use = dims.use)
# #saveRDS(lig.tumor, paste0(outdir, "lig.tumor.rda"))
# 
# p.louvain  <- plotByDatasetAndCluster(lig.tumor, 
#                                       axis.labels = c("UMAP 1", "UMAP 2"),
#                                       return.plots = T)
# 
# p.louvain[[2]]

# Test with Harmony instead ----------------------------------------------------

library(harmony)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seu <- rliger::ligerToSeurat(lig.tumor.filt)
seu <- Seurat::NormalizeData(seu, verbose = T)
#seu[["percent.mt"]] <- PercentageFeatureSet(object = seu, pattern = "^MT-")
#seu[["percent.ribo"]] <- PercentageFeatureSet(object = seu, pattern = "^RP[SL]")
#seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, 
#                        set.ident = TRUE)
seu <- FindVariableFeatures(seu, selection.method = "dispersion", nfeatures = 4000)
seu <- ScaleData(seu, verbose = F, do.scale = F, do.center = T)
seu <- RunPCA(seu, pc.genes = seu@var.genes, npcs = 90, verbose = T)
seu <- RunHarmony(seu, "orig.ident", plot_convergence = TRUE, lambda = 3)#,
                  #epsilon.cluster = -Inf, 
                  #epsilon.harmony = -Inf)
seu <- Seurat::RunUMAP(seu, reduction = "harmony", dims = 1:90)

FeaturePlot(object = seu, reduction = "umap", pt.size = .1,
            features = c("MLANA","CD3D","NCR1","CD4","CD8A","CD14","HBA1",
                         "JCHAIN","MS4A1","FOXP3","CD34","MKI67"))

FeaturePlot(object = seu, reduction = "umap", pt.size = .1,
            features = c("IL7R","TCF7","LAG3","HAVCR2","TIGIT","CD68",
                         "HLA-DRA","PDCD1","GZMA","GZMH"))

# Harmony with sctransform -----------------------------------------------------

seu.list <- SplitObject(seu, split.by="orig.ident")
seu.list <- lapply(X = seu.list, 
                       FUN = SCTransform, 
                       method = "glmGamPoi", 
                       return.only.var.genes = FALSE)
var.features <- SelectIntegrationFeatures(object.list = seu.list, 
                                          nfeatures = 4000)

seu.sct <- merge(x = seu.list[[1]], y = seu.list[2:length(seu.list)], 
                 merge.data=TRUE)
VariableFeatures(seu.sct) <- var.features
seu.sct <- RunPCA(seu.sct, verbose = FALSE)
seu.sct <- RunHarmony(seu.sct, assay.use="SCT", group.by.vars = "orig.ident")
seu.sct <- RunUMAP(seu.sct, reduction = "harmony", dims = 1:30)

# FastMNN ----------------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
seu <- rliger::ligerToSeurat(lig.tumor.filt)
#seu.list <- SplitObject(seu, split.by = "orig.ident")
#seu.list <- lapply(seu.list,"SCTransform")
seu <- Seurat::NormalizeData(seu)
#seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 4000)
seu <- RunFastMNN(object.list = SplitObject(seu, split.by = "orig.ident"),
                  features = 4000, auto.merge = T, d = 150, k = 15, prop.k=0.01)
#seu <- RunFastMNN(object.list = seu.list,
#                  features = 4000, auto.merge = T)
seu <- RunUMAP(seu, reduction = "mnn", dims = 1:150)

#d = 90 n = 3000, k = 20 or 15
#d = 150 n = 4000, k = 15
#prop.k=0.05

FeaturePlot(object = seu, reduction = "umap", pt.size = .1,
            features = c("MLANA","CD3D","NCR1","CD4","CD8A","CD14","HBA1",
                         "JCHAIN","MS4A1","FOXP3","CD34","MKI67",
                         "CD1C","THBD"))

FeaturePlot(object = seu, reduction = "umap", pt.size = .1,
            features = c("IL7R","TCF7","LAG3","HAVCR2","TIGIT","CD68",
                         "HLA-DRA","PDCD1","GZMA","GZMH","TNFRSF9"))

#cos.norm=TRUE
#prop.k

# scVI attempt -----------------------------------------------------------------
library(reticulate)
use_condaenv("scvi-env", required=TRUE)
library(sceasy)
library(Seurat)
library(rliger)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

outdir <- "/data/proj/um_ss/Investigations/seurat/results/v15/"
lig.tumor.filt <- readRDS(paste0(outdir,"lig.tumor.filt.rda"))

k <- 60
dims.use <- 1:k
lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, 
                                quantiles = 100, eps = 0.01)

seu <- rliger::ligerToSeurat(lig.tumor.filt)
seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 4000)
top4000 <- head(VariableFeatures(seu), 4000)
seu <- seu[top4000]

adata <- convertFormat(seu, from="seurat", to="anndata", main_layer="counts", 
                       drop_single_values=FALSE)
scvi$model$SCVI$setup_anndata(adata)
model = scvi$model$SCVI(adata)
model$train()



# Investigate marker genes of LIGER clusters with Seurat DotPlot and reannotate ----
lig.tumor_clusters <- readRDS(
  "~/proj/um_ss/Investigations/seurat/results/liger_all/lig.tumor_clusters.rda")
lig.tumor_clusters <- setNames(as.character(lig.tumor_clusters), 
                               names(lig.tumor_clusters))

dat.cca <- AddMetaData(dat.cca,metadata = lig.tumor_clusters, 
                       col.name = "lig.tumor_clusters")
dat.cca@meta.data$lig.tumor_clusters[
  is.na(dat.cca@meta.data$lig.tumor_clusters)] <- "X"

dat.cca.lig <- subset(dat.cca, lig.tumor_clusters!="X")

pdf("~/proj/um_ss/Investigations/dat.cca.lig.dotplot.pdf", width = 6, height = 8)
DotPlot(dat.cca.lig, 
        features = c("CD3D","CD3G","CD8A","CD4",
                     "NCAM1","NCR1","CD19","MS4A1",
                     "CD14","TNFRSF21","JCHAIN","ITGAX",
                     "MZB1",
                     "HBA1","HBA2","HBB",
                     "PMEL","MLANA","TYR"),
        cluster.idents = F, assay = "RNA") + vertical_xlabels
dev.off()

dat.cca.lig@meta.data$clus <- Idents(dat.cca.lig)
Idents(dat.cca.lig) <- dat.cca.lig@meta.data$lig.tumor_clusters

dat.cca.lig <- RenameIdents(dat.cca.lig,
                            `2` = "CD8 T cell",
                            `3` = "B cell",
                            `4` = "CD8 T cell",
                            `5` = "Monocyte/DC",
                            `7` = "CD4 T cell",
                            `8` = "NK cell",
                            `11` = "CD8 T cell",
                            `12` = "Melanocytic",
                            `13` = "pDC (I)",
                            `14` = "pDC (II)",
                            `15` = "Melanocytic",
                            `16` = "Melanocytic",
                            `17` = "Melanocytic",
                            `18` = "Melanocytic",
                            `19` = "Melanocytic",
                            `20` = "Melanocytic",
                            `21` = "CD8 T cell",
                            `22` = "Melanocytic",
                            `24` = "Monocyte/DC",
                            `25` = "CD8 T cell",
                            `26` = "CD8 T cell",
                            `30` = "CD8 T cell",
                            `33` = "NK cell",
                            `34` = "CD8 T cell",
                            `37` = "Melanocytic",
                            `38` = "Melanocytic",
                            `39` = "Melanocytic",
                            `40` = "CD4 T cell",
                            `41` = "Monocyte/DC",
                            `42` = "Melanocytic",
                            `43` = "CD8 T cell",
                            `44` = "Erythrocyte",
                            `45` = "Melanocytic",
                            `46` = "Monocyte/DC",
                            `47` = "CD8 T cell",
                            `48` = "Melanocytic",
                            `49` = "pDC (III)",
                            `50` = "CD8 T cell")

#saveRDS(Idents(dat.cca.lig),
#        "~/proj/um_ss/Investigations/seurat/results/liger_all/lig.reannotated.rda")

# ----------------
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

# 16: Senescence
# 28,31,32,36,42: Translation, NMD
# 35: DNA repair
# 43: apoptosis, DNA damage
# 40: Tregs
# 12: The outlier melanoma
i <- 12
gene_loadings[[i]]
g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]), decreasing = T),]
#head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10)

# all.plots.reannotated[[2]] + gene_loadings[[i]]

# Plot specific genes
# rliger::plotGene(lig.tumor,"TRBV30",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"TRGV10",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"TNFRSF9",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"HAVCR2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"PDCD1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"ENTPD1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"LAG3",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"TRGV9",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"TRDV2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"XIST",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# 
# rliger::plotGene(lig.tumor,"CD3D",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"CD4",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"CD8A",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"CD8B",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"FOXP3",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# 
# rliger::plotGene(lig.tumor,"NCAM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"NCR1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# 
# rliger::plotGene(lig.tumor,"PMEL",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"MLANA",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"TYR",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# 
# rliger::plotGene(lig.tumor,"CPXM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"APOE",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# 
# rliger::plotGene(lig.tumor,"HBA1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"HBA2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"HBB",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# 
# rliger::plotGene(lig.tumor,"JCHAIN",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"CD19",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"MS4A1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
# rliger::plotGene(lig.tumor,"CD14",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

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