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
#outdir <- "~/proj/um_ss/Investigations/seurat/results/tils/"
#dir.create(outdir,showWarnings = F,recursive = T)

dat.blood <- extract_raw(dat.cca=dat.cca, project.name = "G18-049")
dat.tumor <- extract_raw(dat.cca=dat.cca, project.name = "GWA-JN-388")
dat.blood_tumor <- extract_raw(dat.cca=dat.cca, project.name = c("G18-049","GWA-JN-388"))
dat.blood_tumor_tils <- extract_raw(dat.cca=dat.cca, project.name = unique(dat.cca@meta.data$project.name))

# LIGER ------------------------------------------------------------------------

lig.blood <- seuratToLiger(objects = dat.blood,
                     combined.seurat = F,
                     names = "use-meta",
                     meta.var = "Sample.ID",
                     renormalize=F,
                     use.idents=T,
                     use.seurat.genes=F,
                     remove.missing=F)

lig.tumor <- seuratToLiger(objects = dat.tumor,
                           combined.seurat = F,
                           names = "use-meta",
                           meta.var = "Sample.ID",
                           renormalize=F,
                           use.idents=T,
                           use.seurat.genes=F,
                           remove.missing=F)

lig.blood_tumor <- seuratToLiger(objects = dat.blood_tumor,
                           combined.seurat = F,
                           names = "use-meta",
                           meta.var = "Sample.ID",
                           renormalize=F,
                           use.idents=T,
                           use.seurat.genes=F,
                           remove.missing=F)

lig.blood_tumor_tils <- seuratToLiger(objects = dat.blood_tumor_tils[which(names(dat.blood_tumor_tils)!="2B")],
                                 combined.seurat = F,
                                 names = "use-meta",
                                 meta.var = "Sample.ID",
                                 renormalize=F,
                                 use.idents=T,
                                 use.seurat.genes=F,
                                 remove.missing=F)

lig.blood <- normalize(lig.blood)
lig.blood <- selectGenes(lig.blood, combine = "union", do.plot = T)
lig.blood <- scaleNotCenter(lig.blood)

lig.tumor <- normalize(lig.tumor)
lig.tumor <- selectGenes(lig.tumor, combine = "union", do.plot = T)
lig.tumor <- scaleNotCenter(lig.tumor)

lig.blood_tumor <- normalize(lig.blood_tumor)
lig.blood_tumor <- selectGenes(lig.blood_tumor, combine = "union", do.plot = T)
lig.blood_tumor <- scaleNotCenter(lig.blood_tumor)

lig.blood_tumor_tils <- normalize(lig.blood_tumor_tils)
lig.blood_tumor_tils <- selectGenes(lig.blood_tumor_tils, combine = "union", do.plot = T)
lig.blood_tumor_tils <- scaleNotCenter(lig.blood_tumor_tils)

lig.blood <- optimizeALS(lig.blood, k = 50, lambda = 5, verbose = T)
lig.tumor <- optimizeALS(lig.tumor, k = 50, lambda = 5, verbose = T)
lig.blood_tumor <- optimizeALS(lig.blood_tumor, k = 100, lambda = 5, verbose = T)
lig.blood_tumor_tils <- optimizeALS(lig.blood_tumor_tils, k = 200, lambda = 5, verbose = T)

saveRDS(lig.blood,paste0(outdir,"lig.blood.rda"))
saveRDS(lig.tumor,paste0(outdir,"lig.tumor.rda"))
saveRDS(lig.blood_tumor,paste0(outdir,"lig.blood_tumor.rda"))
saveRDS(lig.blood_tumor_tils,paste0(outdir,"lig.blood_tumor_tils.rda"))

gsea_output <- rliger::runGSEA(lig.blood_tumor_tils)

# Investigate iNMF components --------------------------------------------------


celltypes <- unlist(lapply(dat.blood_tumor_tils[which(names(dat.blood_tumor_tils)!="2B")],function(x){
  setNames(x@meta.data$active.ident.updated_2,rownames(x@meta.data))
}))
names(celltypes) <- str_split_fixed(names(celltypes),"\\.",2)[,2]
identical(names(celltypes),rownames(lig.blood_tumor_tils@cell.data))

celltypes[grep("CD4 T cells",celltypes)] <- "CD4 T cells"
celltypes[grep("CD8 T cells",celltypes)] <- "CD8 T cells"
celltypes[grep("Melanocytic",celltypes)] <- "Melanocytic"
celltypes[grep("CD14 Monocytes",celltypes)] <- "CD14 Monocytes"
celltypes[grep("NK cells",celltypes)] <- "NK cells"

tils <- readRDS(file=paste0("~/proj/um_ss/Investigations/seurat/results/tils/","tils.rda"))
snames_all <- names(tils)
tils.liger <- merge(tils[[which(names(tils) %in% snames_all)[1]]],
                                      y = tils[setdiff(which(names(tils) %in% snames_all),
                                                       which(names(tils) %in% snames_all)[1])],
                                      add.cell.ids = names(tils[which(names(tils) %in% snames_all)]),
                                      project = "tils")
tils.liger <- AddMetaData(tils.liger,
            metadata = tils.liger@active.ident,
            col.name = "old.ident")

all(gsub("SampleID_[0-9]_11june18_","",rownames(tils.liger@meta.data)) %in% 
      names(lig.blood_tumor_tils@clusters))

stopifnot(identical(rownames(lig.blood_tumor_tils@cell.data)[
  grep("SampleID_",lig.blood_tumor_tils@cell.data$dataset)],
          gsub("SampleID_[0-9]_11june18_","",rownames(tils.liger@meta.data))))

rownames(lig.blood_tumor_tils@cell.data)[
  grep("SampleID_",lig.blood_tumor_tils@cell.data$dataset)]

stopifnot(identical(names(celltypes)[grep("SampleID_",lig.blood_tumor_tils@cell.data$dataset)],
          gsub("SampleID_[0-9]_11june18_","",rownames(tils.liger@meta.data))))

celltypes[grep("SampleID_",lig.blood_tumor_tils@cell.data$dataset)] <- as.character(tils.liger@meta.data$old.ident)
celltypes[celltypes=="Gamma-delta T cells (TRGV9, TRDV2)"] <- "Gamma-delta T cells (TRDV2, TRGV9)"
celltypes[grep("CD4 T cells",celltypes)] <- "T cells (CD4)"

lig.blood_tumor_tils@cell.data$celltypes <- celltypes

lig.blood_tumor_tils@cell.data$project <- ""
lig.blood_tumor_tils@cell.data$project[grep("^SampleID",lig.blood_tumor_tils@cell.data$dataset)] <- "TILs"
lig.blood_tumor_tils@cell.data$project[grep("^[0-9]",lig.blood_tumor_tils@cell.data$dataset)] <- "Blood"
lig.blood_tumor_tils@cell.data$project[grep("^[A-D][0-9]",lig.blood_tumor_tils@cell.data$dataset)] <- "Tumor"
project <- lig.blood_tumor_tils@cell.data$project
names(project) <- rownames(lig.blood_tumor_tils@cell.data)

lig.blood_tumor_tils@cell.data$pemdac_post <- "Not pemdac"
lig.blood_tumor_tils@cell.data$pemdac_post[grep("^[0-9]A",lig.blood_tumor_tils@cell.data$dataset)] <- "Pre"
lig.blood_tumor_tils@cell.data$pemdac_post[grep("^[0-9]B",lig.blood_tumor_tils@cell.data$dataset)] <- "Post"
pemdac_post <- lig.blood_tumor_tils@cell.data$pemdac_post
names(pemdac_post) <- rownames(lig.blood_tumor_tils@cell.data)

#Unsure: 198,9,11,15,23, 26,27,32,37,42,48,54,63,88,90,94,96,97,105,128,134,148,149,154,162,167,169,170,173,176,186,193,199,200,5,6,19,36,99,119,152,174
#Unsure, dataset specific: 2, 3, 4, 7, 12,16,17,18,20,21,22,29,30,31,33,34,35,38,40,41,43,
# 44,45,46,49,50,51,53,55,56,57,58,59,60,61,62,64,65,67,69,70,73,74,75,76,77,78,80,81,83,84,85,86,89,
# 91,93,95,98,100,102,103,106,107,108,109,110,112,114,118,120,121,122,124,125,126,127,130,131,135,138,
# 141,142,144,145,146,150,151,153,156,158,160,161,164,165,168,171,172,178,179,181,183,188,190,192,194,
# 195,196,197,30
dims.remove <- c(2,3,4,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,32,
                 33,34,35,36,37,38,39,40,41,42,43,44,45,46,48,49,50,51,52,53,54,55,56,57,
                 58,59,60,61,62,63,64,65,66,67,69,70,73,74,75,76,77,78,80,81,83,84,85,86,
                 88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,105,106,107,108,109,110,
                 112,113,114,118,119,120,121,122,123,124,125,126,127,128,130,131,134,135,137,
                 138,139,140,141,142,143,144,145,146,148,149,150,151,153,154,155,156,158,
                 160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,
                 178,179,181,182,183,185,186,188,190,192,194,195,196,197,193,198,199,200)
dims.use <- setdiff(1:200, dims.remove)

lig.blood_tumor_tils <- quantile_norm(lig.blood_tumor_tils, quantiles = 100, ref_dataset = "C7-GEX", 
                     eps = 0.9, min_cells = 20, knn_k = 20, max_sample = 1000,
                     do.center = T,refine.knn = T,
                     dims.use = dims.use) #knn_k = 2500,min_cells = 3000, max_sample = 10000, 
lig.blood_tumor_tils <- runUMAP(lig.blood_tumor_tils, distance = 'euclidean',
                                dims.use = dims.use) #min_dist = 0.0000000001, n_neighbors = 10)#, n_neighbors = 20, min_dist = 0.000001,

all.plots <- plotByDatasetAndCluster(lig.blood_tumor_tils, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T, clusters = celltypes)
all.plots.proj <- plotByDatasetAndCluster(lig.blood_tumor_tils, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T, clusters = project)
all.plots.pemdac_post <- plotByDatasetAndCluster(lig.blood_tumor_tils, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                          return.plots = T, clusters = pemdac_post)

all.plots[[2]] + all.plots.proj[[2]]

gene_loadings <- plotGeneLoadings(lig.blood_tumor_tils, do.spec.plot = T, return.plots = TRUE,
                                  factor.share.thresh = Inf)
setdiff(which(unlist(lapply(gene_loadings,length))>0),dims.remove)
dims.remove

i <- 2
head(unlist(gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),c(1)]),n=30)
gene_loadings[[i]]
all.plots[[2]] + gene_loadings[[i]]
all.plots.proj[[2]] + gene_loadings[[i]]
all.plots.pemdac_post[[2]] + gene_loadings[[i]]

lig.blood_tumor_tils@clusters
stopifnot(identical(rownames(lig.blood_tumor_tils@cell.data),names(lig.blood_tumor_tils@clusters)))

saveRDS(lig.blood_tumor_tils,file = paste0(outdir,"lig.blood_tumor_tils.integrated.rda"))

celltypes <- lig.blood_tumor_tils@cell.data$celltypes
names(celltypes) <- rownames(lig.blood_tumor_tils@cell.data)
celltypes <- factor(celltypes)

lig.blood_tumor_tils@cell.data$clusters.old <- lig.blood_tumor_tils@clusters
lig.blood_tumor_tils@clusters <- celltypes

m.ct <- plotClusterFactors(lig.blood_tumor_tils,return.data = T)

p <- pheatmap(m.ct,cutree_cols = 25)

lig.blood_tumor_tils@clusters <- project
m.proj <- plotClusterFactors(lig.blood_tumor_tils,return.data = T)
p.proj <- pheatmap(m.proj,cutree_cols = 20)

#proj_spec <- apply(m.proj,2,sd)/apply(m.proj,2,mean)
m.proj2 <- apply(m.proj,2,function(x){x/mean(x)})
proj_spec <- apply(m.proj2,2,function(x){max(x)-max(setdiff(x,max(x)))})
proj_spec <- proj_spec/mean(proj_spec)
#ct_spec <- apply(m.ct,2,sd)/apply(m.ct,2,mean)
m.ct2 <- apply(m.ct,2,function(x){x/mean(x)})
ct_spec <- apply(m.ct2,2,function(x){max(x)-max(setdiff(x,max(x)))})
ct_spec <- ct_spec/mean(ct_spec)

pheatmap(rbind(m.ct2,m.proj2))
whitelisted <- c(47,104,152,116,180,
                 136,147,169,81,144,100,117,19,13,132,
                 129,184,56,97,178,17,139,181,186,82,107,
                 1,72,
                 191,119,187,79,189,193,
                 61,
                 159,
                 87,24,68,157,59,111,5,98,115,110,160)

p <- pheatmap(rbind(m.ct2,m.proj2),cutree_cols = 15)
clus <- cutree(p$tree_col,15)
annot <- as.data.frame(clus)
annot$clus <- as.character(annot$clus)
annot$selected <- as.character(annot$clus %in% c(10,1,4,8,14,6,3,15,12))
pheatmap(rbind(m.ct2,m.proj2),cutree_cols = 15,annotation_col = annot)

blacklist <- c(136,147,169,81,144,100,
               17,104,152)
#dims.use <- p$tree_col$order[6:100]
#dims.use <- as.numeric(names(which(ct_spec-proj_spec > 
#                                     quantile(sort(ct_spec-proj_spec,decreasing = T),seq(0,1,0.05))
#                                   [["80%"]])))
# dims.use <- as.numeric(names(which(ct_spec > 
#                                      quantile(sort(ct_spec,decreasing = T),seq(0,1,0.05))
#                                    [["80%"]])))

#dims.use <- whitelisted

dims.use <- as.numeric(rownames(annot)[annot$selected=="TRUE"])
dims.use <- setdiff(dims.use,blacklist)

#gsea_output

lig.blood_tumor_tils <- quantile_norm(lig.blood_tumor_tils, quantiles = 100, ref_dataset = "C7-GEX", 
                                      eps = 0.9, min_cells = 20, knn_k = 20, max_sample = 1000,
                                      do.center = T,refine.knn = T,
                                      dims.use = dims.use) #knn_k = 2500,min_cells = 3000, max_sample = 10000, 
lig.blood_tumor_tils <- runUMAP(lig.blood_tumor_tils, distance = 'euclidean',
                                dims.use = dims.use) #min_dist = 0.0000000001, n_neighbors = 10)#, n_neighbors = 20, min_dist = 0.000001,

all.plots <- plotByDatasetAndCluster(lig.blood_tumor_tils, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T, clusters = celltypes)
all.plots.proj <- plotByDatasetAndCluster(lig.blood_tumor_tils, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                          return.plots = T, clusters = project)
all.plots[[2]] + all.plots.proj[[2]]

#lapply(gsea_output[dims.remove],function(x){head(unlist(x[order(unlist(x[,5]),decreasing = T),c(1)]),n=3)})

# Candidates for removal: 19, 23, 148, 26,27,32,37,42,48,54,63,88,90,94,96,
# 97,99,105,128,134,148,149,154,162,167,169,170,173,176,186, 193, 200,99

# Might add back: 169,176,193,19,36,119,80,107,178,179,61

# # Subset TILs
# lig.blood_tumor_tils.tils <- subsetLiger(lig.blood_tumor_tils, 
#                                          cells.use = rownames(lig.blood_tumor_tils@cell.data[
#                                            grep("SampleID",lig.blood_tumor_tils@cell.data$dataset),]),
#                                          remove.missing = F)
# 
# celltypes.tils <- lig.blood_tumor_tils.tils@cell.data$celltypes
# 
# all.plots <- plotByDatasetAndCluster(lig.blood_tumor_tils.tils, axis.labels = c('UMAP 1', 'UMAP 2'), 
#                                      return.plots = T, clusters = celltypes.tils)
# all.plots[[1]] + all.plots[[2]]
# 
# # Subset Blood
# cnames <- rownames(lig.blood_tumor_tils@cell.data[
#   lig.blood_tumor_tils@cell.data$dataset %in% 
#     c("1A","1B","2A","6A","6B","7A","7B"),])
# 
# lig.blood_tumor_tils.blood <- subsetLiger(lig.blood_tumor_tils, 
#                                          cells.use = cnames,
#                                          remove.missing = F)
# 
# celltypes.blood <- lig.blood_tumor_tils.blood@cell.data$celltypes
# 
# all.plots <- plotByDatasetAndCluster(lig.blood_tumor_tils.blood, axis.labels = c('UMAP 1', 'UMAP 2'), 
#                                      return.plots = T, clusters = celltypes.blood)
# all.plots[[1]] + all.plots[[2]]
# 
# # Subset Tumor
# cnames <- rownames(lig.blood_tumor_tils@cell.data[
#   lig.blood_tumor_tils@cell.data$dataset %in% 
#     c("A2-GEX","B3-GEX","C5-GEX","C7-GEX","C8-GEX","D2-GEX","D4-GEX","D6-GEX"),])
# 
# lig.blood_tumor_tils.tumor <- subsetLiger(lig.blood_tumor_tils, 
#                                           cells.use = cnames,
#                                           remove.missing = F)
# 
# celltypes.tumor <- lig.blood_tumor_tils.tumor@cell.data$celltypes
# 
# all.plots <- plotByDatasetAndCluster(lig.blood_tumor_tils.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
#                                      return.plots = T, clusters = celltypes.tumor)
# all.plots[[1]] + all.plots[[2]]

# Tumor only -------------------------------------------------------------------
outdir <- "~/proj/um_ss/Investigations/seurat/results/liger_all/"
lig.blood_tumor_tils <- readRDS(file = paste0(outdir,"lig.blood_tumor_tils.integrated.rda"))
celltypes <- lig.blood_tumor_tils@cell.data$celltypes
names(celltypes) <- rownames(lig.blood_tumor_tils@cell.data)
celltypes <- factor(celltypes)
#rm(lig.blood_tumor_tils)

lig.tumor <- readRDS(paste0(outdir,"lig.tumor.rda"))


gsea_output <- rliger::runGSEA(lig.tumor)

dims.remove <- c(6,23,29)
dims.remove <- c(dims.remove,28,31,32,35,36,42,43)

dims.use <- setdiff(1:50,dims.remove)
#lig.tumor <- quantile_norm(lig.tumor,ref_dataset = "C7-GEX",dims.use=dims.use)#,
lig.tumor <- quantile_norm(lig.tumor,ref_dataset = "C8-GEX",dims.use=dims.use,
                           quantiles = 100, eps = 0.01)
lig.tumor <- runUMAP(lig.tumor, distance = "cosine", dims.use=dims.use, 
                     n_neighbors = 20, min_dist = 0.05)

all.plots.reannotated <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                                 return.plots = T, clusters = lig.reannotated)

all.plots.reannotated[[2]]

all.plots.clusters <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T)
all.plots.clusters[[2]]
clus <- c(50,2,43,34,26,30)
cells <- names(lig.tumor@clusters)[lig.tumor@clusters %in% clus]
saveRDS(cells,paste0(outdir,"cells.rda"))

saveRDS(lig.tumor@clusters,paste0(outdir,"lig.tumor_clusters.rda"))

stopifnot(identical(names(lig.tumor@clusters),names(celltypes3)))
saveRDS(celltypes3,paste0(outdir,"celltypes3.rda"))

vasu_phenotype <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/vasu_phenotype.rda")
all.plots.vasu <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T, clusters = vasu_phenotype,do.shuffle = F,
                                     reorder.idents = T,new.order = c("None","Unidentified","MART1","CD3+41bb+"))

all.plots.reannotated[[2]] + all.plots.vasu[[2]] + all.plots.clusters[[2]]

gene_loadings <- plotGeneLoadings(lig.tumor, do.spec.plot = T, return.plots = TRUE,
                                  factor.share.thresh = Inf)
setdiff(which(unlist(lapply(gene_loadings,length))>0),dims.remove)


# 16: Senescence
# 28,31,32,36,42: Translation, NMD
# 35: DNA repair
# 43: apoptosis, DNA damage
i <- 3
g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10)
gene_loadings[[i]]

all.plots.reannotated[[2]]
all.plots.reannotated[[2]] + gene_loadings[[i]]

all.plots.clusters[[2]]





# louvain

lig.tumor.louvain <- louvainCluster(lig.tumor, resolution = 1)
all.plots.louvain  <- plotByDatasetAndCluster(lig.tumor.louvain, 
                                              axis.labels = c('UMAP 1', 'UMAP 2'), 
                                              return.plots = T)
all.plots.louvain[[2]]

#all.plots.clusters[[2]] + all.plots.louvain[[2]]
all.plots.reannotated[[2]] + all.plots.louvain[[2]]

# 40: Tregs

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

rliger::plotGene(lig.tumor,"HBA1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"HBA2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"HBB",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"JCHAIN",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD19",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"MS4A1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD14",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

lig.tumor
length(celltypes)
celltypes <- celltypes[rownames(lig.tumor@cell.data)]
celltypes <- setNames(as.character(celltypes),names(celltypes))
stopifnot(identical(rownames(lig.tumor@cell.data),names(celltypes)))
unique(celltypes)
celltypes[celltypes %in% c("Gamma-delta T cells (TRDV2, TRGV9)","Gamma-delta T cells (mixed)")] <- "CD8 T cells"

snames <- names(lig.tumor@raw.data)
for (sname in snames){
  # quantile(lig.tumor@raw.data[[1]]["TRGV9",colnames(lig.tumor@raw.data[[1]]) %in% names(celltypes)[celltypes=="CD8 T cells"]],
  #          seq(0,1,by=0.01))
  # quantile(lig.tumor@raw.data[[1]]["TRDV2",colnames(lig.tumor@raw.data[[1]]) %in% names(celltypes)[celltypes=="CD8 T cells"]],
  #          seq(0,1,by=0.01))
  # 
  # table(lig.tumor@raw.data[[1]]["TRGV9",colnames(lig.tumor@raw.data[[1]]) %in% names(celltypes)[celltypes=="CD8 T cells"]] > 0 & 
  #         lig.tumor@raw.data[[1]]["TRDV2",colnames(lig.tumor@raw.data[[1]]) %in% names(celltypes)[celltypes=="CD8 T cells"]] > 0)
  # 
  # table(lig.tumor@raw.data[[1]]["TRDV2",colnames(lig.tumor@raw.data[[1]]) %in% names(celltypes)[celltypes=="CD8 T cells"]] > 0)
  # 
  celltypes[names(which(lig.tumor@raw.data[[sname]]["TRGV9",colnames(lig.tumor@raw.data[[sname]]) %in% 
                                                      names(celltypes)[celltypes=="CD8 T cells"]] > 0 & 
    lig.tumor@raw.data[[sname]]["TRDV2",colnames(lig.tumor@raw.data[[sname]]) %in% 
                                  names(celltypes)[celltypes=="CD8 T cells"]] > 0))] <- "CD8 T cells (TRGV9+TRDV2+)"
  celltypes[names(which(lig.tumor@raw.data[[sname]]["TRGV9",colnames(lig.tumor@raw.data[[sname]]) %in% 
                                                      names(celltypes)[celltypes=="CD8 T cells"]] > 0))] <- "CD8 T cells (TRGV9+)"
}


lig.reannotated <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/lig.reannotated.rda")
lig.reannotated <- setNames(as.character(lig.reannotated),names(lig.reannotated))

all.plots.reannotated <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T, clusters = lig.reannotated)


all.plots[[2]] + all.plots.clusters[[2]] + all.plots.reannotated[[2]]

all.plots[[1]] + all.plots.reannotated[[2]]

tmp <- lig.reannotated
tmp[tmp!="Melanocytic"] <- "Other"

all.plots.tmp <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                                 return.plots = T, clusters = tmp,
                                         do.shuffle = F,reorder.idents = T,
                                         new.order = c("Other","Melanocytic"))

all.plots.tmp[[2]]

#saveRDS(lig.tumor,paste0(outdir,"lig.tumor.reannotated.rda"))
#lig.tumor.bak <- lig.tumor

suspicious <- list()
for (i in 1:length(lig.tumor@raw.data)){
  suspicious[[i]] <- names(which(lig.tumor@raw.data[[i]]["CD3D",] > 1 & 
                              (lig.tumor@raw.data[[i]]["MLANA",] > 1 & 
                                 lig.tumor@raw.data[[i]]["PMEL",] > 1 & 
                                 lig.tumor@raw.data[[i]]["TYR",] > 1)))
}
suspicious <- unlist(suspicious)
suspicious <- setNames(rep("Suspicious",length(suspicious)),suspicious)
length(suspicious)

all.plots.suspicious <- plotByDatasetAndCluster(lig.tumor, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                                 return.plots = T, clusters = suspicious)

all.plots.suspicious[[2]]
