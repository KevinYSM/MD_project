library(Seurat)
library(tidyverse)
library(rliger)
#library(WriteXLS)
library(scater)
library(djvdj)
library(scDblFinder)

vertical_xlabels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                                     hjust=1))

outdir <- "/data/proj/um_ss/Investigations/seurat/results/v14/"
dir.create(outdir, showWarnings = F, recursive = T)

outdir <- "~/proj/um_ss/Investigations/seurat/results/v14/"
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

p <- plotByDatasetAndCluster(lig.tumor, return.plots = T)
p[[1]] + p[[2]]

gene_val <- getGeneValues(lig.tumor@norm.data,gene = "HBA1")
stopifnot(identical(names(gene_val),rownames(lig.tumor@cell.data)))
lig.tumor@cell.data <- cbind(lig.tumor@cell.data, gene_val)


# Add VDJ data -----------------------------------------------------------------
add_vdj_to_liger <- function(lig.tumor, vdj_paths){
  library(djvdj)
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

fnames <- Sys.glob("~/nimbus/data/proj/um_ss/Pipelines/10x/results/multi/*/outs/per_sample_outs/*/vdj_t/")
fnames <- fnames[basename(dirname(fnames)) %in% unique(lig.tumor@cell.data$dataset)]
names(fnames) <- basename(dirname(fnames))
fnames <- fnames[names(fnames) %in% lig.tumor@cell.data$dataset]

lig.tumor <- add_vdj_to_liger(lig.tumor, vdj_paths = fnames)

saveRDS(lig.tumor,file = paste0(outdir,"/lig.tumor.rda"))

# Doublets ---------------------------------------------------------------------
library(scDblFinder)

#outdir <- "~/nimbus/data/proj/um_ss/Investigations/seurat/results/v10/"
#lig.tumor <- readRDS(paste0(outdir, "lig.tumor.rda"))

# The authors do not recommend further QC filtering than this before doublet detection
idx_keep <- lig.tumor@cell.data[["nUMI"]] > 200

cells_keep <- rownames(lig.tumor@cell.data[idx_keep,])

# For some reason, one cell has been removed from scale.data, but is kept in the 
# rest of the object, causing an error with subsetLiger
setdiff(rownames(lig.tumor@cell.data),rownames(lig.tumor@H.norm))

setdiff(rownames(lig.tumor@cell.data),unique(unlist(lapply(lig.tumor@raw.data,colnames))))
setdiff(unique(unlist(lapply(lig.tumor@raw.data,colnames))),rownames(lig.tumor@cell.data))
length(unique(unlist(lapply(lig.tumor@raw.data,colnames))))
length(rownames(lig.tumor@cell.data))

setdiff(unique(unlist(lapply(lig.tumor@norm.data,colnames))),rownames(lig.tumor@cell.data))
setdiff(rownames(lig.tumor@cell.data),unique(unlist(lapply(lig.tumor@norm.data,colnames))))
length(unique(unlist(lapply(lig.tumor@norm.data,colnames))))
length(rownames(lig.tumor@cell.data))

setdiff(unique(unlist(lapply(lig.tumor@scale.data,rownames))),rownames(lig.tumor@cell.data))
setdiff(rownames(lig.tumor@cell.data),unique(unlist(lapply(lig.tumor@scale.data,rownames))))
length(unique(unlist(lapply(lig.tumor@scale.data,rownames))))
length(rownames(lig.tumor@cell.data))

intersect(cells_keep,
          setdiff(rownames(lig.tumor@cell.data),unique(unlist(lapply(lig.tumor@scale.data,rownames)))))

cells_keep <- setdiff(cells_keep,
                      setdiff(rownames(lig.tumor@cell.data),unique(unlist(lapply(lig.tumor@scale.data,rownames)))))

lig.tumor.filt <- rliger::subsetLiger(lig.tumor, cells.use = cells_keep)

seu <- rliger::ligerToSeurat(lig.tumor.filt)
sce <- as.SingleCellExperiment(seu)

library(BiocParallel)
sce.standard <- scDblFinder(sce, samples = "orig.ident", BPPARAM = MulticoreParam(4))

# Investigate outcome of doublet detection -------------------------------------
mat <- counts(sce)
colnames(mat) <- str_split_fixed(colnames(mat),"_",2)[,2]

stopifnot(identical(colnames(mat),rownames(lig.tumor.filt@cell.data)))
idx <- (!is.na(lig.tumor.filt@cell.data$cdr3) | colSums(mat[c("CD3D","CD4","CD8A","CD8B"),]>0) > 0) & 
  colSums(mat[c("MLANA","PMEL","TYR","HBA1","HBA2","HBB"),]>0) > 0
idx <- idx | ((!is.na(lig.tumor.filt@cell.data$cdr3) & 
                 colSums(mat[c("CD14","CD19","MS4A1","JCHAIN"),]>0) > 0))
idx <- idx | ((colSums(mat[c("MLANA","PMEL","TYR"),]>0) > 0) & 
                (colSums(mat[c("CD14","CD19","MS4A1","JCHAIN"),]>0) > 0))
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

rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.standard.score", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.standard.class", by.dataset = F,
                    discrete = T)
rliger::plotFeature(object = lig.tumor.filt, feature = "knownDoublets", 
                    by.dataset = F,discrete = T,pt.size = 0.01)

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
  
  tab_clus <- table(lig.tumor.filt@cell.data[cells_clus,
                                             "is_doublet"])
  tab_rest <- table(lig.tumor.filt@cell.data[
    ! rownames(lig.tumor.filt@cell.data) %in% cells_clus,
    "is_doublet"])
  
  stopifnot(sum(tab_clus) + sum(tab_rest) ==
              nrow(lig.tumor.filt@cell.data))
  
  tab <- rbind(tab_clus,tab_rest)
  ft <- fisher.test(tab)
  
  doublet_cluster_fisher[[clus]] <- data.frame(
    clus=clus, 
    tab=paste0(as.character(tab),collapse="_"), 
    p=ft$p.value, 
    odds_ratio=ft$estimate,
    percent_doublets=tab_clus["TRUE"]/sum(tab_clus))
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

plot(doublet_cluster_fisher$percent_doublets,doublet_cluster_fisher$median_nUMI)
abline(h=thresh,v=0.35)

#plot(doublet_cluster_fisher$percent_doublets,doublet_cluster_fisher$median_nGene)
#abline(v=0.35)

doublet_clusters <- as.numeric(doublet_cluster_fisher[which(
  doublet_cluster_fisher$q < 0.05 & 
    doublet_cluster_fisher$odds_ratio < 1 & 
    (doublet_cluster_fisher$percent_doublets >= 0.35) | 
    (doublet_cluster_fisher$median_nUMI > thresh)),]$clus)

cells_doublet_clusters <- names(lig.tumor.filt@clusters[
  lig.tumor.filt@clusters %in% doublet_clusters])

lig.tumor.filt@cell.data$in_doublet_cluster <- 
  rownames(lig.tumor.filt@cell.data) %in% cells_doublet_clusters

rliger::plotFeature(object = lig.tumor.filt, 
                    feature = "in_doublet_cluster", 
                    by.dataset = F,discrete = T,pt.size = 0.01)
rliger::plotFeature(object = lig.tumor.filt, 
                    feature = "is_doublet", 
                    by.dataset = F,discrete = T,pt.size = 0.01)

lig.tumor.filt@cell.data$is_doublet <- lig.tumor.filt@cell.data$scDblFinder.standard.class=="doublet" | 
  lig.tumor.filt@cell.data$knownDoublets | lig.tumor.filt@cell.data$in_doublet_cluster

rliger::plotFeature(object = lig.tumor.filt, 
                    feature = "is_doublet", 
                    by.dataset = F,discrete = T,pt.size = 0.01)

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


rownames(lig.tumor.filt@cell.data)

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
p <- pheatmap(cor(as.matrix(mat),method = "spearman"),show_colnames = F,show_rownames = F,annotation_col = annot)
clus <- cutree(p$tree_col,2)
table(clus)

cl <- which.max(c(sum(mat["MLANA",names(which(clus==1))]>0),sum(mat["MLANA",names(which(clus==2))]>0)))

hba1_blacklist <- names(which(clus==cl))

lig.tumor.filt@cell.data$hba1_blacklist <- F
lig.tumor.filt@cell.data[names(which(clus==cl)),]$hba1_blacklist <- T

ggplot(lig.tumor.filt@cell.data[lig.tumor.filt@cell.data$HBA1 > 0,],
       aes(x=nGene,y=nUMI, color = hba1_blacklist)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

min_nGene_erythrocytes <- summary(annot[names(which(clus==2)),]$nGene)[["1st Qu."]]
min_nUMI_erythrocytes <- summary(annot[names(which(clus==2)),]$nUMI)[["1st Qu."]]
summary(annot[names(which(clus==2)),]$percent.top_100)
###

### Potential neutrophils
#seu@assays$RNA@data["NAMPT",] > 0 | 
idx <- (seu@assays$RNA@data["CSF3R",] > 0 | 
  seu@assays$RNA@data["FPR1",] > 0 | 
  seu@assays$RNA@data["FCGR3B",] > 0 | 
    seu@assays$RNA@data["MNDA",] > 0) & ! seu@assays$RNA@data["CD14",] > 0
#   ! seu@assays$RNA@data["CD19",] > 0 & 
#   ! seu@assays$RNA@data["CD79A",] > 0 & 
#   ! seu@assays$RNA@data["CD79B",] > 0 & 
#   ! seu@assays$RNA@data["CD3D",] > 0 & 
#   ! seu@assays$RNA@data["CD3E",] > 0 & 
#   ! seu@assays$RNA@data["CD4",] > 0 & 
#   ! seu@assays$RNA@data["CD8A",] > 0 & 
#   ! seu@assays$RNA@data["CD8B",] > 0 & 
#   ! seu@assays$RNA@data["MLANA",] > 0 & 
#   ! seu@assays$RNA@data["PMEL",] > 0
mat <- seu@assays$RNA@data[,idx]
mat <- mat[which(rowSums(mat) > 0),]
colnames(mat) <- str_split_fixed(colnames(mat),"_",2)[,2]
dim(mat)

#CD68 <- getGeneValues(lig.tumor.filt@norm.data,gene = "CD68")
#MS4A1 <- getGeneValues(lig.tumor.filt@norm.data,gene = "MS4A1")
#stopifnot(identical(names(CD68),rownames(lig.tumor.filt@cell.data)))
#lig.tumor.filt@cell.data <- cbind(lig.tumor.filt@cell.data, CD68)
#lig.tumor.filt@cell.data <- cbind(lig.tumor.filt@cell.data, MS4A1)

ggplot(lig.tumor.filt@cell.data[colnames(mat),],
       aes(x=nGene,y=nUMI, color = MLANA)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

annot <- as.data.frame(t(mat[intersect(rownames(mat),c("CD14","CD19","CD79A","CD3D","NCAM1","MLANA","TYR",
                               "CSF3R","FPR1","FCGR3B","MNDA","NAMPT","IL3RA")),]))
annot <- merge(annot, lig.tumor.filt@cell.data[,c("nGene","nUMI","is_doublet","percent.top_100")],
               all.x=T, all.y=F, by="row.names")
#annot <- lig.tumor.filt@cell.data[,c("nGene","nUMI","is_doublet","percent.top_100")]
annot$is_doublet <- factor(ifelse(annot$is_doublet,1,0))
#annot <- annot[colnames(mat),]
rownames(annot) <- annot$Row.names
annot$Row.names <- NULL
annot <- merge(annot,lig.tumor.filt@tsne.coords,by="row.names")
rownames(annot) <- annot$Row.names
annot$Row.names <- NULL
cr <- cor(as.matrix(mat),method = "spearman")
p <- pheatmap(cr,show_colnames = F,
              show_rownames = F, annotation_col = annot, cutree_cols = 19)

clus <- cutree(p$tree_col,19)
table(clus)

blacklist_cluster <- c()
for (i in 1:19){
  print(i)
  print(table(annot[names(which(clus==i)),]$is_doublet))
  doublet_enriched <- F
  cd3d_enriched <- F
  if (any(as.numeric(as.character(annot[names(which(clus==i)),]$is_doublet))==1)){
    doublet_enriched <- table(annot[names(which(clus==i)),]$is_doublet)[1] /
      sum(table(annot[names(which(clus==i)),]$is_doublet)) > 0.1
  }
  if (any(annot[names(which(clus==i)),]$CD3D > 0)){
    cd3d_enriched <- table(annot[names(which(clus==i)),]$CD3D > 0)["TRUE"]/
      sum(table(annot[names(which(clus==i)),]$CD3D > 0)) > 0.1
  }
  blacklist_cluster[i] <- doublet_enriched | cd3d_enriched
}

names(clus)[clus %in% which(!blacklist_cluster)]


p <- pheatmap(cr,show_colnames = F,
              show_rownames = F, annotation_col = annot, cutree_cols = 19)

annot$cluster <- NA
for (i in 1:19){
  annot$cluster[rownames(df) %in% names(which(clus==i))] <- i
}
annot$cluster <- factor(df$cluster)

p <- pheatmap(cr,show_colnames = F,
              show_rownames = F, annotation_col = annot, cutree_cols = 19)

df <- annot
df$perc_CD3D <- NA
df$perc_MLANA <- NA
df$perc_CD79A <- NA
df$perc_NCAM1 <- NA
df$mean_nGene <- NA
df$mean_nUMI <- NA
df$perc_is_doublet <- NA
for (i in 1:19){
  idx <- rownames(df) %in% names(which(clus==i))
  df[idx,]$perc_CD3D <- sum(df[idx,]$CD3D>0)/sum(idx)
  df[idx,]$perc_MLANA <- sum(df[idx,]$MLANA>0)/sum(idx)
  df[idx,]$perc_CD79A <- sum(df[idx,]$CD79A>0)/sum(idx)
  df[idx,]$perc_NCAM1 <- sum(df[idx,]$NCAM1>0)/sum(idx)
  df[idx,]$mean_nGene <- median(df[idx,]$nGene)
  df[idx,]$mean_nUMI <- median(df[idx,]$nUMI)
  df[idx,]$perc_is_doublet <- sum(as.numeric(as.character(df[idx,]$is_doublet)))/sum(idx)
}
df <- unique(df[,c("cluster","perc_CD3D","perc_MLANA","perc_CD79A",
      "perc_NCAM1","mean_nGene","mean_nUMI","perc_is_doublet")])

ggplot(df,aes(x = perc_is_doublet, y=perc_CD3D, color = perc_MLANA > 0.1 | 
                perc_CD79A > 0.1 | perc_NCAM1 > 0.1)) + 
  geom_text(aes(label=cluster))

ggplot(df,aes(x = perc_is_doublet, y=perc_CD3D, color = mean_nGene > 1000)) + 
  geom_text(aes(label=cluster))


annot$whitelist <- F
annot$whitelist[annot$cluster %in% as.numeric(as.character(df$cluster[df$perc_is_doublet < 0.1]))] <- T
annot$whitelist[annot$CD3D > 0 | annot$MLANA > 0 | annot$NCAM1 > 0 | 
                  annot$CD79A > 0 | annot$CD19 > 0 | annot$is_doublet == 1] <- F
annot$whitelist <- factor(annot$whitelist)

p <- pheatmap(cr,show_colnames = F,
              show_rownames = F, annotation_col = annot, cutree_cols = 19)



#lig.tumor.filt@cell.data <- lig.tumor.filt@cell.data[,grep("neutrophil",colnames(lig.tumor.filt@cell.data),invert = T)]
lig.tumor.filt@cell.data$neutrophil_blacklist <- rownames(lig.tumor.filt@cell.data) %in% rownames(annot)[annot$whitelist!=TRUE]
lig.tumor.filt@cell.data$neutrophil_whitelist <- rownames(lig.tumor.filt@cell.data) %in% rownames(annot)[annot$whitelist==TRUE]

rliger::plotFeature(object = lig.tumor.filt,feature = "neutrophil_blacklist", 
                    by.dataset = F,discrete = T)

rliger::plotFeature(object = lig.tumor.filt,feature = "neutrophil_whitelist", 
                    by.dataset = F,discrete = T)

# Myeloblasts and other granulocytic precursors do not express CD14, but neutrophils and a small proportion of B lymphocytes may weakly express CD14. T cells, dendritic cells, and platelets are CD14-negative. 
#https://www.sciencedirect.com/book/9780128098431/atlas-of-hematopathology

# mat <- seu@assays$RNA@data[,(seu@assays$RNA@data["CD14",]>0 | 
#                              str_split_fixed(colnames(seu@assays$RNA@data),"_",2)[,2] %in% 
#                              rownames(lig.tumor.filt@cell.data)[
#                                lig.tumor.filt@cell.data$neutrophil_whitelist & 
#                                  ! lig.tumor.filt@cell.data$is_doublet
#                              ]) & 
#                              str_split_fixed(colnames(seu@assays$RNA@data),"_",2)[,2] %in% 
#                              rownames(lig.tumor.filt@cell.data)[
#                                  ! lig.tumor.filt@cell.data$is_doublet
#                              ]
#                            ]
# mat <- seu@assays$RNA@data[,(seu@assays$RNA@data["CD14",]>0 | 
#                                str_split_fixed(colnames(seu@assays$RNA@data),"_",2)[,2] %in% 
#                                rownames(lig.tumor.filt@cell.data)[
#                                  lig.tumor.filt@cell.data$neutrophil_whitelist
#                                ])
# ]

mat <- seu@assays$RNA@data[,(
                               str_split_fixed(colnames(seu@assays$RNA@data),"_",2)[,2] %in% 
                               rownames(lig.tumor.filt@cell.data)[
                                 lig.tumor.filt@cell.data$neutrophil_whitelist & 
                                   ! lig.tumor.filt@cell.data$is_doublet
                               ])
]


mat <- mat[which(rowSums(mat) > 0),]
colnames(mat) <- str_split_fixed(colnames(mat),"_",2)[,2]
dim(mat)

cr <- cor(as.matrix(mat),method = "spearman")


#In humans, neutrophils are distinguished from eosinophils and monocytes based on the expression of both CD15 and CD16/Fc gamma RIII on human neutrophils, along with the lack of expression of CD14. In addition, CD66b/CEACAM-8, CD11b/Integrin alpha M, CD33, and the cytoplasmic marker, myeloperoxidase, are other common markers that are used to identify human neutrophils.
#https://www.rndsystems.com/resources/cell-markers/immune-cells/granulocytes/neutrophil-cell-markers

#CD15: FUT4
#CD16: FCGR3A
#CD11b: ITGAM
#Myeloperoxidase: MPO

#Both PBMCs and NGs preparations contained cells that were positive for CD-68 and either neutrophil elastase (NE), or myeloperoxidase (MPO). CD-68(+)/NE(-)/MPO(-) cells were regarded as monocytes.
#https://pubmed.ncbi.nlm.nih.gov/23573303/

annot <- as.data.frame(t(mat[intersect(rownames(mat),c("CD14",
                                                       "CSF3R","FPR1","FCGR3B",
                                                       "MNDA",
                                                       "CD33","FUT4","FCGR3A","ITGAM",
                                                       "MPO","ELANE","CD24","CEACAM8",
                                                       "CD68")),]))
annot <- merge(annot, lig.tumor.filt@cell.data[,c("nGene","nUMI","percent.top_100",
                                                  "neutrophil_whitelist")],
               all.x=T, all.y=F, by="row.names")
annot$neutrophil_whitelist <- factor(ifelse(annot$neutrophil_whitelist,1,0))
rownames(annot) <- annot$Row.names
annot$Row.names <- NULL

p <- pheatmap(cr,show_colnames = F,
              show_rownames = F, annotation_col = annot, cutree_cols = 19)
clus <- cutree(p$tree_col,19)

lig.tumor.filt.mono_neutro_macro <- rliger::subsetLiger(lig.tumor.filt,
                                                        cells.use = 
                                                          rownames(lig.tumor.filt@cell.data)[
                                                            lig.tumor.filt@cell.data$neutrophil_whitelist & 
                                                              ! lig.tumor.filt@cell.data$is_doublet
                                                          ])

dims.use <- 1:90
lig.tumor.filt.mono_neutro_macro <- quantile_norm(lig.tumor.filt.mono_neutro_macro, 
                                                  dims.use = dims.use, 
                                quantiles = 100, eps = 0.01,knn_k = 3)
lig.tumor.filt.mono_neutro_macro <- rliger::runUMAP(lig.tumor.filt.mono_neutro_macro, 
                                                    dims.use = dims.use, 
                                  n_neighbors = 5, min_dist = 0.05)

p <- plotByDatasetAndCluster(lig.tumor.filt.mono_neutro_macro, return.plots = T,pt.size = 1)
p[[1]] + p[[2]]

rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"FUT4",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"CD14",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"CD68",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"ITGAM",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"HLA-DRA",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"ITGA4",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"CD24",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"CD83",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"CEACAM8",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"CD15",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"CD33",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"MNDA",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"NAMPT",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"FCGR3A",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"FCGR3B",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotGene(lig.tumor.filt.mono_neutro_macro,"CSF3R",
                 axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none",pt.size = 1)
rliger::plotFeature(lig.tumor.filt.mono_neutro_macro,"nGene",
                 axis.labels = c('UMAP 1', 'UMAP 2'),by.dataset = F,discrete = T,pt.size = 1)
rliger::plotFeature(lig.tumor.filt.mono_neutro_macro,"is_doublet",
                    axis.labels = c('UMAP 1', 'UMAP 2'),by.dataset = F,discrete = T,pt.size = 1)
rliger::plotFeature(lig.tumor.filt.mono_neutro_macro,"dataset",
                    axis.labels = c('UMAP 1', 'UMAP 2'),by.dataset = F,discrete = T,pt.size = 1)
rliger::plotFeature(lig.tumor.filt.mono_neutro_macro,"neutrophil_whitelist",
                    axis.labels = c('UMAP 1', 'UMAP 2'),by.dataset = F,discrete = T,pt.size = 1)


#ggplot(annot, aes(x=CD14,y=CSF3R)) + geom_point()
#ggplot(annot, aes(x=CD14,y=CD68)) + geom_point()

lig.tumor.filt.mono_neutro_macro <- louvainCluster(lig.tumor.filt.mono_neutro_macro)

whitelist <- setNames(annot$neutrophil_whitelist,rownames(annot))
p <- plotByDatasetAndCluster(lig.tumor.filt.mono_neutro_macro, 
                             return.plots = T,clusters = as.factor(clus),pt.size = 1)
p[[2]]
p <- plotByDatasetAndCluster(lig.tumor.filt.mono_neutro_macro, 
                             return.plots = T,clusters = as.factor(whitelist),pt.size = 1)
p[[2]]


wc <- runWilcoxon(lig.tumor.filt.mono_neutro_macro,compare.method = "clusters")
table(wc$group)
wc.sig <- wc[which(wc$padj < 0.05 & wc$logFC>1 & wc$auc > 0.5),]
wc.sig <- wc.sig[order(wc.sig$logFC,decreasing = T),]
wc.sig[wc.sig$group==0,]
wc.sig[wc.sig$group==1,]
wc.sig[wc.sig$group==2,]
wc.sig[wc.sig$group==3,]
wc.sig[wc.sig$group==4,]
wc.sig[wc.sig$group==5,]
wc.sig[wc.sig$group==6,]
wc.sig[wc.sig$group==7,]
wc.sig[wc.sig$group==8,]


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
  facet_wrap(~ dataset)

# Final filtering --------------------------------------------------------------
idx_erythrocyte <- lig.tumor.filt@cell.data$HBA1 > 0 & 
  !lig.tumor.filt@cell.data$hba1_blacklist & 
  lig.tumor.filt@cell.data[["nGene"]] > min_nGene_erythrocytes & 
  lig.tumor.filt@cell.data[["nUMI"]] > min_nUMI_erythrocytes

idx_keep <- ((lig.tumor.filt@cell.data[["nGene"]] > 500 & 
  lig.tumor.filt@cell.data[["nUMI"]] > 500 & 
  lig.tumor.filt@cell.data[["percent.top_100"]] < 80) | idx_erythrocyte) & 
  lig.tumor.filt@cell.data[["percent_mito"]] < 0.3 & 
  lig.tumor.filt@cell.data[["percent_ribo"]] < 0.4 & 
  !lig.tumor.filt@cell.data$is_doublet

table(idx_keep)

cells_keep <- rownames(lig.tumor.filt@cell.data[idx_keep,])

lig.tumor.filt@cell.data$cells_filtered <- ! rownames(lig.tumor.filt@cell.data) %in% cells_keep
table(lig.tumor.filt@cell.data$cells_filtered)
table(lig.tumor.filt@cell.data$cells_filtered)["TRUE"]/nrow(lig.tumor.filt@cell.data)
table(lig.tumor.filt@cell.data$cells_filtered,
      lig.tumor.filt@cell.data$dataset)

rliger::plotFeature(object = lig.tumor.filt, feature = "cells_filtered", by.dataset = F, discrete = T)

lig.tumor.filt <- rliger::subsetLiger(lig.tumor.filt, cells.use = cells_keep)
saveRDS(lig.tumor.filt, paste0(outdir, "lig.tumor.filt.rda"))

# Reprocess filtered dataset ---------------------------------------------------

#lig.tumor.filt <- readRDS(paste0(outdir, "lig.tumor.filt.rda"))

k <- 60
lig.tumor.filt <- optimizeALS(lig.tumor.filt, k = k, lambda = 5, verbose = T)

pdf(file=paste0(outdir,"suggestK_l5.filt.pdf"),width = 10,height = 10)
sk <- suggestK(lig.tumor.filt, k.test = seq(10,100,10), lambda = 5, 
               return.data = T, num.cores = 5)
dev.off()

#saveRDS(lig.tumor.filt, paste0(outdir, "lig.tumor.filt.rda"))


#lig.tumor.filt <- readRDS(paste0(outdir, "lig.tumor.filt.rda"))

dims.use <- 1:k
lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, 
                                quantiles = 100, eps = 0.01)
lig.tumor.filt <- rliger::runUMAP(lig.tumor.filt, dims.use = dims.use, 
                                  n_neighbors = 10, min_dist = 0.05)

#saveRDS(lig.tumor.filt, paste0(outdir, "lig.tumor.filt.rda"))

# Select factors, reprocess, Annotate clusters ---------------------------------
#lig.tumor.filt <- readRDS(paste0(outdir, "lig.tumor.filt.rda"))

p <- plotByDatasetAndCluster(lig.tumor.filt, return.plots = T)
p[[1]] + p[[2]]

gsea_output <- rliger::runGSEA(lig.tumor.filt)

gene_loadings <- plotGeneLoadings(lig.tumor.filt, do.spec.plot = T,
                                  return.plots = TRUE, factor.share.thresh = Inf)

for (i in dims.use) {
  #  i <- 10
  print(i)
  g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
  print(head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10))
  #print(gene_loadings[[i]])
  invisible(readline(prompt="Press [enter] to continue"))
}

dims.use <- setdiff(1:60,c(8,42)) # cell cycle
lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, 
                                quantiles = 100, eps = 0.01)
lig.tumor.filt <- rliger::runUMAP(lig.tumor.filt, dims.use = dims.use, 
                                  n_neighbors = 10, min_dist = 0.05)

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

pdf("~/proj/um_ss/Investigations/seurat/results/v13/plotFactors.pdf",width = 20)
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

rliger::plotGene(lig.tumor.filt,"CD14",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"CD69",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"CD68",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"ACTA2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"ALB",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"VIM",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"PDGFRA",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor.filt,"PECAM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CD34",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")


rliger::plotFeature(object = lig.tumor.filt, feature = "percent_ribo", by.dataset = F, discrete = T)
rliger::plotFeature(object = lig.tumor.filt, feature = "percent_mito", by.dataset = F, discrete = T)

length(intersect(names(lig.tumor.filt@clusters)[which(lig.tumor.filt@clusters==7)],
                 intersect(
                   names(which(getGeneValues(lig.tumor@norm.data,gene = "CPXM1")>0)),
                   names(which(getGeneValues(lig.tumor@norm.data,gene = "CD3D")>0))
                 )
))/
  length(names(lig.tumor.filt@clusters)[which(lig.tumor.filt@clusters==7)])


seu <- rliger::ligerToSeurat(lig.tumor.filt)

DotPlot(seu, 
        features = c("CD3D","CD3G","CD8A","CD4",
                     "NCAM1","NCR1","CD19","MS4A1",
                     "CD14","TNFRSF21","JCHAIN","ITGAX",
                     "MZB1",
                     "HBA1","HBA2","HBB",
                     "PMEL","MLANA","TYR"),
        cluster.idents = F, assay = "RNA") + vertical_xlabels


# lig.tumor <- louvainCluster(lig.tumor, resolution = 1, dims.use = dims.use)
# #saveRDS(lig.tumor, paste0(outdir, "lig.tumor.rda"))
# 
# p.louvain  <- plotByDatasetAndCluster(lig.tumor, 
#                                       axis.labels = c("UMAP 1", "UMAP 2"),
#                                       return.plots = T)
# 
# p.louvain[[2]]


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