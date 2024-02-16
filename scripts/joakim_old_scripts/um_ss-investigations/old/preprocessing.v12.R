library(Seurat)
library(tidyverse)
library(rliger)
#library(WriteXLS)
library(scater)
library(djvdj)
library(scDblFinder)
library(scCustomize)

vertical_xlabels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                                     hjust=1))

outdir <- "/data/proj/um_ss/Investigations/seurat/results/v13/"
#outdir <- "~/proj/um_ss/Investigations/seurat/results/v12/"
dir.create(outdir, showWarnings = F, recursive = T)

# Read cellbender data ---------------------------------------------------------
fnames <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/multi/*/outs/multi/count/cellbender_feature_bc_matrix/cellbender_feature_bc_matrix_filtered.h5")

snames <- basename(dirname(dirname(dirname(dirname(dirname(fnames))))))
idx <- ! snames %in% c("1A","1B","2A","2B","6A","6B","7A","7B") & 
  ! grepl("^SampleID",snames) & ! grepl("^PDX",snames)
fnames <- fnames[idx]
snames <- snames[idx]

names(fnames) <- snames
fnames <- as.list(fnames)

seu <- Read10X_h5_Multi_Directory(
  base_path="/data/proj/um_ss/Pipelines/10x/results/multi",
  secondary_path = "outs/multi/count/cellbender_feature_bc_matrix/",
  h5_filename = "cellbender_feature_bc_matrix_filtered.h5",
  default_10X_path = F,
  cell_bender = F,
  sample_list = names(fnames)
)

seu.bak <- seu

for (sname in names(seu)){
  colnames(seu[[sname]]) <- paste0(sname,"_",colnames(seu[[sname]]))
}

lig.tumor.cellbender <- rliger::createLiger(seu)
lig.tumor.cellbender <- rliger::normalize(lig.tumor.cellbender)
lig.tumor.cellbender <- selectGenes(lig.tumor.cellbender)
lig.tumor.cellbender <- scaleNotCenter(lig.tumor.cellbender)

#saveRDS(lig.tumor.cellbender,file = paste0(outdir,"/lig.tumor.cellbender.rda"))

lig.tumor.cellbender <- optimizeALS(lig.tumor.cellbender, k = 90, lambda = 5, verbose = T)

dims.use <- 1:90
lig.tumor.cellbender <- quantile_norm(lig.tumor.cellbender, dims.use = dims.use, 
                                      quantiles = 100, eps = 0.01)

lig.tumor.cellbender <- rliger::runUMAP(lig.tumor.cellbender, 
                                        dims.use = dims.use, 
                                        n_neighbors = 20, min_dist = 0.05)

#lig.tumor.cellbender <- readRDS(file = paste0(outdir,"/lig.tumor.cellbender.rda"))

pdf(file=paste0(outdir,"suggestK_l5.cellbender.pdf"),width = 10,height = 10)
sk <- suggestK(lig.tumor.cellbender, k.test = seq(10,100,10), lambda = 5, 
               return.data = T, num.cores = 1)
dev.off()


saveRDS(lig.tumor.cellbender,file = paste0(outdir,"/lig.tumor.cellbender.rda"))

p <- plotByDatasetAndCluster(lig.tumor.cellbender, return.plots = T)
p[[1]] + p[[2]]

lig.tumor <- lig.tumor.cellbender

gene_val <- getGeneValues(lig.tumor@norm.data,gene = "HBA1")
stopifnot(identical(names(gene_val),rownames(lig.tumor@cell.data)))
lig.tumor@cell.data <- cbind(lig.tumor@cell.data, gene_val)

# Read data and process with LIGER ---------------------------------------------

dnames <- Sys.glob("/data/proj/um_ss/Pipelines/10x/results/multi/*/outs/per_sample_outs/*/count/sample_filtered_feature_bc_matrix")
snames <- basename(dirname(dirname(dnames)))
idx <- ! snames %in% c("1A","1B","2A","2B","6A","6B","7A","7B") & ! grepl("^SampleID",snames)
dnames <- dnames[idx]
snames <- snames[idx]

matrix_list <- rliger::read10X(sample.dirs = dnames, 
                               sample.names = snames, merge = F)
matrix_list <- lapply(matrix_list,function(x){x[[1]]})
for (sname in names(matrix_list)){
  colnames(matrix_list[[sname]]) <- paste0(sname,"_",colnames(matrix_list[[sname]]))
}

lig.tumor <- createLiger(matrix_list)
lig.tumor <- normalize(lig.tumor)
lig.tumor <- selectGenes(lig.tumor)
lig.tumor <- scaleNotCenter(lig.tumor)

#vasu_phenotype <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/vasu_phenotype.rda")
pdf(file=paste0(outdir,"suggestK_l5.pdf"),width = 10,height = 10)
sk <- suggestK(lig.tumor.filt, k.test = seq(10,100,10), lambda = 5, return.data = T, num.cores = 5)
dev.off()
lig.tumor <- optimizeALS(lig.tumor, k = 90, lambda = 5, verbose = T)

lig.tumor@cell.data[["project"]] <- ""
lig.tumor@cell.data[["project"]][
  lig.tumor@cell.data$dataset %in% 
    c("A6","A8","B6","B9C1","C3CY","D5","D9","E5")] <- "GWA-JN-605"
lig.tumor@cell.data[["project"]][
  lig.tumor@cell.data$dataset %in% 
    c("A2","B3","C5","C7","C8","D2","D4","D6")] <- "GWA-JN-388"

dims.use <- 1:90
lig.tumor <- quantile_norm(lig.tumor, #ref_dataset = "C8-GEX", 
                           dims.use = dims.use, quantiles = 100, eps = 0.01)
lig.tumor <- runUMAP(lig.tumor, distance = "cosine", dims.use = dims.use, 
                     n_neighbors = 20, min_dist = 0.05)

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

#saveRDS(lig.tumor, paste0(outdir, "lig.tumor.rda"))

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

lig.tumor.cellbender.filt <- rliger::subsetLiger(lig.tumor, cells.use = cells_keep)

seu <- rliger::ligerToSeurat(lig.tumor.cellbender.filt)
sce <- as.SingleCellExperiment(seu)

library(BiocParallel)
sce.standard <- scDblFinder(sce, samples = "orig.ident", BPPARAM = MulticoreParam(4))

# Investigate outcome of doublet detection -------------------------------------
mat <- exprs(sce)
colnames(mat) <- str_split_fixed(colnames(mat),"_",2)[,2]

stopifnot(identical(colnames(mat),rownames(lig.tumor.cellbender.filt@cell.data)))
idx <- (!is.na(lig.tumor.cellbender.filt@cell.data$cdr3) | colSums(mat[c("CD3D","CD4","CD8A","CD8B"),]>0) > 0) & 
  colSums(mat[c("MLANA","PMEL","TYR","HBA1","HBA2","HBB"),]>0) > 0
idx <- idx | ((!is.na(lig.tumor.cellbender.filt@cell.data$cdr3) & 
                 colSums(mat[c("CD14","CD19","MS4A1","JCHAIN"),]>0) > 0))
idx <- idx | ((colSums(mat[c("MLANA","PMEL","TYR"),]>0) > 0) & 
                (colSums(mat[c("CD14","CD19","MS4A1","JCHAIN"),]>0) > 0))
knownDoublets <- names(which(idx))
lig.tumor.cellbender.filt@cell.data$knownDoublets <- rownames(lig.tumor.cellbender.filt@cell.data) %in% knownDoublets

knownDoublets <- str_split_fixed(rownames(colData(sce)),"_",2)[,2] %in% knownDoublets
colData(sce) <- cbind(colData(sce),knownDoublets)

table(colData(sce.standard)[["scDblFinder.class"]])

stopifnot(identical(str_split_fixed(rownames(colData(sce.standard)),"_",2)[,2],
                    rownames(lig.tumor.cellbender.filt@cell.data)))

lig.tumor.cellbender.filt@cell.data$scDblFinder.standard.score <- 
  colData(sce.standard)[["scDblFinder.score"]]
lig.tumor.cellbender.filt@cell.data$scDblFinder.standard.class <- 
  colData(sce.standard)[["scDblFinder.class"]]

rliger::plotFeature(object = lig.tumor.cellbender.filt, feature = "scDblFinder.standard.score", by.dataset = F)
rliger::plotFeature(object = lig.tumor.cellbender.filt, feature = "scDblFinder.standard.class", by.dataset = F,discrete = T)
rliger::plotFeature(object = lig.tumor.cellbender.filt, feature = "knownDoublets", 
                    by.dataset = F,discrete = T,pt.size = 0.01)

table(doublet_agreement=lig.tumor.cellbender.filt@cell.data$scDblFinder.standard.class)["doublet"] / 
  nrow(lig.tumor.cellbender.filt@cell.data)

table(lig.tumor.cellbender.filt@cell.data$knownDoublets,
      lig.tumor.cellbender.filt@cell.data$scDblFinder.standard.class)

p[[2]]

lig.tumor.cellbender.filt.louvain <- louvainCluster(lig.tumor.cellbender.filt,
                                                    resolution = 1.3)

p.louvain <- plotByDatasetAndCluster(lig.tumor.cellbender.filt.louvain, return.plots = T)
p.louvain[[2]]

lig.tumor.cellbender.filt.louvain@cell.data$is_doublet <- lig.tumor.cellbender.filt.louvain@cell.data$scDblFinder.standard.class=="doublet" | 
  lig.tumor.cellbender.filt.louvain@cell.data$knownDoublets

doublet_cluster_fisher <- list()
for (clus in unique(lig.tumor.cellbender.filt.louvain@clusters)){
  cells_clus <- names(lig.tumor.cellbender.filt.louvain@clusters)[
    which(lig.tumor.cellbender.filt.louvain@clusters == clus)]
  
  tab_clus <- table(lig.tumor.cellbender.filt.louvain@cell.data[cells_clus,
                                                                "is_doublet"])
  tab_rest <- table(lig.tumor.cellbender.filt.louvain@cell.data[
    ! rownames(lig.tumor.cellbender.filt.louvain@cell.data) %in% cells_clus,
    "is_doublet"])
  
  stopifnot(sum(tab_clus) + sum(tab_rest) ==
              nrow(lig.tumor.cellbender.filt.louvain@cell.data))
  
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

doublet_clusters <- as.numeric(doublet_cluster_fisher[which(
  doublet_cluster_fisher$q < 0.05 & 
    doublet_cluster_fisher$odds_ratio < 1 & 
    doublet_cluster_fisher$percent_doublets >= 0.3),]$clus)

cells_doublet_clusters <- names(lig.tumor.cellbender.filt.louvain@clusters[
  lig.tumor.cellbender.filt.louvain@clusters %in% doublet_clusters])

lig.tumor.cellbender.filt.louvain@cell.data$in_doublet_cluster <- 
  rownames(lig.tumor.cellbender.filt.louvain@cell.data) %in% cells_doublet_clusters

rliger::plotFeature(object = lig.tumor.cellbender.filt.louvain, 
                    feature = "in_doublet_cluster", 
                    by.dataset = F,discrete = T,pt.size = 0.01)

lig.tumor.filt <- lig.tumor.cellbender.filt
stopifnot(identical(rownames(lig.tumor.filt@cell.data),rownames(lig.tumor.cellbender.filt.louvain@cell.data)))
lig.tumor.filt@cell.data$in_doublet_cluster <- lig.tumor.cellbender.filt.louvain@cell.data$in_doublet_cluster
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

lig.tumor.filt@cell.data$percent_ribo_40 <- lig.tumor.filt@cell.data$percent_ribo > 0.4
rliger::plotFeature(object = lig.tumor.filt,feature = "percent_ribo_40", 
                    by.dataset = F,discrete = T)

hist(lig.tumor.filt@cell.data$percent_mito)
hist(lig.tumor.filt@cell.data$percent_ribo)

ggplot(lig.tumor.filt@cell.data,aes(x=percent_ribo,y=percent_mito,color=dataset)) + geom_point() + theme_classic() + vertical_xlabels
ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent_mito,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent_mito,color=gene_val > 0)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor.filt@cell.data,aes(x=nGene,y=percent_ribo,color=gene_val > 0)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

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
#idx_remove <- lig.tumor.filt@cell.data$knownDoublets | lig.tumor.filt@cell.data$scDblFinder.standard.class == "doublet"
idx_remove <- lig.tumor.filt@cell.data$is_doublet
cells_remove <- rownames(lig.tumor.filt@cell.data[idx_remove,])

idx_keep <- (lig.tumor.filt@cell.data[["nGene"]] > 500 | lig.tumor.filt@cell.data[["gene_val"]] > 0) & 
  lig.tumor.filt@cell.data[["nUMI"]] > 500 & 
  lig.tumor.filt@cell.data[["percent_mito"]] < 0.3 & 
  lig.tumor.filt@cell.data[["percent_ribo"]] < 0.4 &
  lig.tumor.filt@cell.data[["percent.top_100"]] < 80

table(idx_keep)

cells_keep <- rownames(lig.tumor.filt@cell.data[idx_keep,])
cells_keep <- setdiff(cells_keep, cells_remove)

lig.tumor.filt@cell.data$cells_filtered <- ! rownames(lig.tumor.filt@cell.data) %in% cells_keep
table(lig.tumor.filt@cell.data$cells_filtered)
table(lig.tumor.filt@cell.data$cells_filtered)["TRUE"]/nrow(lig.tumor.filt@cell.data)
table(lig.tumor.filt@cell.data$cells_filtered,
      lig.tumor.filt@cell.data$dataset)

lig.tumor.filt@cell.data$cells_qc_filtered <- ! rownames(lig.tumor.filt@cell.data) %in% rownames(lig.tumor.filt@cell.data[idx_keep,])

rliger::plotFeature(object = lig.tumor.filt, feature = "cells_filtered", by.dataset = F, discrete = T)
rliger::plotFeature(object = lig.tumor.filt, feature = "cells_qc_filtered", by.dataset = F, discrete = T)

lig.tumor.filt <- rliger::subsetLiger(lig.tumor.filt, cells.use = cells_keep)
saveRDS(lig.tumor.filt, paste0(outdir, "lig.tumor.cellranger.filt.rda"))

# Reprocess filtered dataset ---------------------------------------------------
#lig.tumor.filt <- readRDS(paste0("/data/proj/um_ss/Investigations/seurat/results/v10/", "lig.tumor.filt.rda"))
# lig.tumor.filt <- optimizeALS(lig.tumor.filt, k = 90, lambda = 5, verbose = T)
# 
# dims.use <- 1:90
# lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, 
#                                 quantiles = 100, eps = 0.01)
# lig.tumor.filt <- rliger::runUMAP(lig.tumor.filt, dims.use = dims.use, 
#                                   n_neighbors = 10, min_dist = 0.05)

pdf(file=paste0(outdir,"suggestK_l5.pdf"),width = 10,height = 10)
sk <- suggestK(lig.tumor.filt, k.test = seq(10,100,10), lambda = 5, 
               return.data = T, num.cores = 5)
dev.off()

lig.tumor.filt <- optimizeALS(lig.tumor.filt, k = 60, lambda = 5, verbose = T)

dims.use <- 1:60
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
49
(55)
(60)
for (i in dims.use) {
  #  i <- 10
  print(i)
  g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
  print(head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10))
  print(gene_loadings[[i]])
  invisible(readline(prompt="Press [enter] to continue"))
}

# library(msigdbr)
# m_set <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") 
# m_set <- m_set[grepl("MITOCHON", m_set$gs_name),]
# m_set <- split(m_set$entrez_gene, f = m_set$gs_name)
# gsea_output_mito <- runGSEA(lig.tumor.filt, custom_gene_sets = m_set)
# 
# for (i in 1:length(gsea_output_mito)){
#   g <- gsea_output_mito[[i]][order(unlist(gsea_output_mito[[i]][,5]),decreasing = T),]
#   g <- as.data.frame(g)
#   g$leadingEdge <- NULL
#   gsea_output_mito[[i]] <- g[which(g$padj < 0.05 & g$NES > 0),]
# }
# n_mito_sets <- data.frame(fact=1:90,n_mito_sets=unlist(lapply(gsea_output_mito,nrow)))
# n_mito_sets$fact
# n_mito_sets$n_mito_sets
# head(n_mito_sets)
# 
# ggplot(n_mito_sets,aes(x=reorder(fact,-n_mito_sets),y=n_mito_sets)) + 
#   geom_bar(stat="identity",position="dodge") + 
#   theme_classic() + 
#   vertical_xlabels
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

pdf("~/proj/um_ss/Investigations/seurat/results/v10/plotFactors.pdf",width = 20)
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
lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, quantiles = 100, eps = 0.01)
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