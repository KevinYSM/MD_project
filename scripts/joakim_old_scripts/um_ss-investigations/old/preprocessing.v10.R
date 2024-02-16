library(Seurat)
library(tidyverse)
library(rliger)
library(WriteXLS)
library(scater)

vertical_xlabels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                                     hjust=1))

#outdir <- "/data/proj/um_ss/Investigations/seurat/results/v10/"
outdir <- "~/proj/um_ss/Investigations/seurat/results/v10/"
dir.create(outdir, showWarnings = F, recursive = T)

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
sk <- suggestK(lig.tumor, k.test = seq(10,100,10), lambda = 5, return.data = T, num.cores = 5)
dev.off()
lig.tumor <- optimizeALS(lig.tumor, k = 90, lambda = 5, verbose = T)

saveRDS(lig.tumor, paste0(outdir, "lig.tumor.rda"))

# Tumor only -------------------------------------------------------------------
#outdir <- "~/nimbus/data/proj/um_ss/Investigations/seurat/results/v10/"
lig.tumor <- readRDS(paste0(outdir, "lig.tumor.rda"))
gsea_output <- rliger::runGSEA(lig.tumor)

dims.use <- 1:90
lig.tumor <- quantile_norm(lig.tumor, #ref_dataset = "C8-GEX", 
                           dims.use = dims.use, quantiles = 100, eps = 0.01)
lig.tumor <- runUMAP(lig.tumor, distance = "cosine", dims.use = dims.use, 
                     n_neighbors = 20, min_dist = 0.05)

p <- plotByDatasetAndCluster(lig.tumor, return.plots = T)
p[[1]] + p[[2]]

gene_loadings <- plotGeneLoadings(lig.tumor,  do.spec.plot = T, return.plots = TRUE,
                                  factor.share.thresh = Inf)

for (i in dims.use) {
#  i <- 10
  print(i)
  g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
  print(head(unlist(g[unlist(g[,3]) < 0.05 & unlist(g[,5]) > 0,c(1)]),n=10))
  print(gene_loadings[[i]])
  invisible(readline(prompt="Press [enter] to continue"))
}

(26,30,33,34,45,48,60,70,74,77,85)
#(77)
(76)
(84)
(35,60)

(15,23)
(79)

(40)

#2,5,7

#4,6,46,86
#(45)
#3,65,(26)
#(7,15)
#69
#dims.remove <- c(16,58,59,87)

library(msigdbr)
m_set <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC") 
m_set <- m_set[grepl("MITOCHON", m_set$gs_name),]
m_set <- split(m_set$entrez_gene, f = m_set$gs_name)
gsea_output_mito <- runGSEA(lig.tumor, custom_gene_sets = m_set)

for (i in 1:length(gsea_output_mito)){
  g <- gsea_output_mito[[i]][order(unlist(gsea_output_mito[[i]][,5]),decreasing = T),]
  g <- as.data.frame(g)
  g$leadingEdge <- NULL
  gsea_output_mito[[i]] <- g[which(g$padj < 0.05 & g$NES > 0),]
}
n_mito_sets <- data.frame(fact=1:90,n_mito_sets=unlist(lapply(gsea_output_mito,nrow)))
n_mito_sets$fact
n_mito_sets$n_mito_sets
head(n_mito_sets)

ggplot(n_mito_sets,aes(x=reorder(fact,-n_mito_sets),y=n_mito_sets)) + 
  geom_bar(stat="identity",position="dodge") + 
  theme_classic() + 
  vertical_xlabels

n_mito_sets <- n_mito_sets[order(n_mito_sets$n_mito_sets,decreasing = T),]
n_mito_sets[n_mito_sets$n_mito_sets>0,]




gsea_output_reactome <- list()
for (i in 1:length(gsea_output)){
  g <- gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),]
  g <- as.data.frame(g)
  g$leadingEdge <- NULL
  gsea_output_reactome[[i]] <- g[which(g$padj < 0.05 & g$NES > 0),]
}

i <- 81
gsea_output_reactome[[i]]
gsea_output_reactome[[i]][order(unlist(gsea_output_reactome[[i]]$pval),decreasing = F),]
gene_loadings[[i]]
n_mito_sets[n_mito_sets$n_mito_sets>0 & n_mito_sets$fact==i,]

#8,9,28,31,39,51,56,62,63,75,80,83,89
#21
#25,32


#sex: 82,73,24,81
n_mito_sets[n_mito_sets$fact %in% sex_p[1:7,]$fact,]

rliger::plotGene(lig.tumor,"XIST",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

pdf("~/proj/um_ss/Investigations/seurat/results/v10/plotFactors.pdf",width = 20)
plotFactors(lig.tumor)
dev.off()

lig.tumor@cell.data[["project"]] <- ""
lig.tumor@cell.data[["project"]][
  lig.tumor@cell.data$dataset %in% 
    c("A6","A8","B6","B9C1","C3CY","D5","D9","E5")] <- "GWA-JN-605"

lig.tumor@cell.data[["project"]][
  lig.tumor@cell.data$dataset %in% 
    c("A2","B3","C5","C7","C8","D2","D4","D6")] <- "GWA-JN-388"

table(lig.tumor@cell.data[["project"]],lig.tumor@cell.data$dataset)


get_gene_loadings <- function(liger_object,
                                  liger_factor){
  W <- t(liger_object@W)
  rownames(W) <- colnames(liger_object@scale.data[[1]])
  gene_loadings <- W[order(W[, liger_factor], decreasing = TRUE),liger_factor]
  
  return(gene_loadings)
}


df <- data.frame(sname=lig.tumor@cell.data$dataset, gene=getGeneValues(lig.tumor@norm.data, "XIST"))
perc_expressed_gene_sample <- list()
for (sname in unique(df$sname)){
  perc_expressed_gene_sample[[sname]] <- sum(df$gene[df$sname==sname] > 0)/length(df$gene[df$sname==sname])
}

sex <- table(names(unlist(perc_expressed_gene_sample)),unlist(perc_expressed_gene_sample) > 0.01)
colnames(sex)[colnames(sex)=="TRUE"] <- "F"
colnames(sex)[colnames(sex)=="FALSE"] <- "M"
sex <- as.data.frame.matrix(sex)

sex$F <- ifelse(sex$F==1,"F","M")
sex$M <- NULL
colnames(sex) <- "sex"

dim(lig.tumor@H.norm)
df$barcode <- rownames(df)
df2 <- merge(df,sex,by.x="sname",by.y="row.names",all.x=T,all.y=F)
rownames(df2) <- df$barcode
df <- df2[rownames(df),]
stopifnot(identical(rownames(df),rownames(lig.tumor@H.norm)))

p_sex <- list()
for (fact in 1:90){
  p_sex[[fact]] <- wilcox.test(lig.tumor@H.norm[df$sex=="M",fact],lig.tumor@H.norm[df$sex=="F",fact])$p.value
  #boxplot(lig.tumor@H.norm[df$sex=="M",fact],lig.tumor@H.norm[df$sex=="F",fact])
}

idx <- order(unlist(p_sex),decreasing = F)
sex_p <- data.frame(fact=(1:90)[idx],p=unlist(p_sex)[idx])
sex_p$sig <- sex_p$p < 0.05

ggplot(sex_p,aes(x=reorder(fact,-p),y=-log10(p))) + geom_point() + vertical_xlabels

sex_p[1:7,]


gene_val <- getGeneValues(lig.tumor@norm.data, "XIST")
stopifnot(identical(names(gene_val),rownames(lig.tumor@H.norm)))

cr <- cor(gene_val, lig.tumor@H.norm, method = "spearman")
cr <- as.data.frame(t(cr))
colnames(cr) <- "correlation"
cr$fact <- 1:nrow(cr)
head(cr)

cr <- cr[order(abs(cr$correlation),decreasing = T),]

ggplot(cr,aes(x=reorder(fact,correlation),y=correlation)) + geom_point() + 
  theme_classic() + vertical_xlabels

ggplot(cr,aes(x=reorder(fact,correlation),y=abs(correlation))) + geom_point() + 
  theme_classic() + vertical_xlabels + geom_abline(intercept = 0.1,slope = 0) + 
  xlab(NULL)

#10,19,62,51,27,37,83

gene_loadings <- list()
for (i in 1:90){
  gene_loadings[[i]] <- get_gene_loadings(lig.tumor,i)
}

lapply(gene_loadings,function(x){which(names(x)=="XIST")})



dim(lig.tumor@V[[1]])
dim(lig.tumor@V[[2]])

dim(lig.tumor@W)

lig.tumor@W[1:10,1:10]

sort(rowSums(lig.tumor@W>0))

head(lig.tumor@norm.data[[1]])

gene_loadings <- list()
for (i in dims.use){
  gene_loadings[[i]] <- get_gene_loadings(liger_object=lig.tumor, liger_factor=i)
}

lapply(gene_loadings,function(x){which(names(x)=="HBA1")})


# find celltype-unspecific factors

markers <- c("MLANA","CD3D","HBA1","NCR1","CD14",
             "JCHAIN","MS4A1","CD19","CD4","CD8A")
gene_val <- list()
for (marker in markers){
  print(marker)
  gene_val[[marker]] <- plot_factor_vs_marker(lig.tumor, marker, do_plot=F)
}

df <- do.call("cbind",lapply(gene_val,function(x){as.numeric(x$ratio)}))
df <- as.data.frame(df)
rownames(df) <- 1:nrow(df)

pheatmap(df,border_color = F,cluster_rows = F,cluster_cols = F)
pheatmap(df,border_color = F)


pheatmap(df[order(apply(df,1,mad),decreasing = T),],
         cluster_cols = T,cluster_rows = F,border_color = NA)
pheatmap(df[order(apply(df,1,sd),decreasing = T),],
         cluster_cols = T,cluster_rows = F,border_color = NA)
pheatmap(df[order(apply(df,1,function(x){
  v <- sort(x,decreasing = T);
  v <- v[1]-v[2];
  return(v)}),decreasing = T),],
         cluster_cols = T,cluster_rows = F,border_color = NA)

apply(df,1,function(x){
  v <- max(abs(x-median(as.numeric(x))));
  return(v)})

pheatmap(df[order(apply(df,1,function(x){
  v <- max(abs(x-median(as.numeric(x))));
  return(v)}),decreasing = T),],
  cluster_cols = T,cluster_rows = F,border_color = NA)

pheatmap(df[apply(df,1,function(x){
  v <- max(abs(x-median(as.numeric(x))));
  return(v)})>1,],
  cluster_cols = T,cluster_rows = T,border_color = NA)

df2 <- df
df2[df2<1] <- 0
pheatmap(df2[order(rowSums(df2==0),decreasing = T),],
         cluster_cols = T,cluster_rows = F,border_color = NA)

pheatmap(df[order(apply(df,1,max),decreasing = T),],
         cluster_rows = F,cluster_cols = T,border_color = NA)

df <- df[order(apply(df,1,max),decreasing = T),]
pheatmap(df[apply(df,1,max) > summary(sort(apply(df,1,max)))[["1st Qu."]],],
         cluster_rows = F,cluster_cols = T,border_color = NA)

dims.use <- sort(as.numeric(rownames(df[apply(df,1,max) > summary(sort(apply(df,1,max)))[["1st Qu."]],])))

plot_factor_vs_marker <- function(lig.tumor, marker, do_plot=T){
  gene_val <- getGeneValues(lig.tumor@norm.data,gene = marker)
  
  x <- lig.tumor@H.norm[which(gene_val>0),]
  colnames(x) <- 1:ncol(x)
  
  y <- lig.tumor@H.norm[which(gene_val==0),]
  colnames(y) <- 1:ncol(y)
  
  df <- data.frame(perc_ct=colSums(x>0)/nrow(x),
             perc_non_ct=colSums(y>0)/nrow(y),
             fact=1:ncol(x))
  df$ratio=df$perc_ct/df$perc_non_ct
  
  if (do_plot){
    ggplot(df,aes(x=reorder(fact,-ratio),y=ratio)) + geom_point() + 
      theme_classic() + vertical_xlabels + xlab(NULL) + 
      geom_bar(aes(x=reorder(fact,-ratio),y=perc_ct),stat="identity",position="dodge",fill="green") + 
      geom_bar(aes(x=reorder(fact,-ratio),y=perc_non_ct),stat="identity",position="dodge",fill="red") + 
      geom_abline(intercept = 1,slope = 0)
  }
  return(df)
}

hba1 <- plot_factor_vs_marker(lig.tumor = lig.tumor, marker = "HBA1")
hba2 <- plot_factor_vs_marker(lig.tumor = lig.tumor, marker = "HBA2")
hbb <- plot_factor_vs_marker(lig.tumor = lig.tumor, marker = "HBB")


#dims.remove <- c(16,58,59,87) # cell cycle
dims.remove <- c(c(c(16,58,59,87),3,65,69),4,6,46,86) # cell cycle + few specific outlier cells
dims.remove <- sort(union(dims.remove,c(2,5,7))) # metabolic effect separating cell types that should belong together
dims.remove <- sort(union(dims.remove,c(8,9,12,28,31,39,51,56,62,63,75,80,83,89))) # plotFactors
dims.remove <- sort(union(dims.remove,c(21,25,32))) # plotFactors, less certain
dims.remove <- sort(union(dims.remove,c(19,62,83))) # sex
dims.remove <- sort(union(dims.remove,c(32,37,39,46))) # top loading lncrnas
dims.remove <- sort(union(dims.remove,c(47,52,55,57,67,17,36,50,
                                        26,30,33,34,45,48,60,70,74,77,85))) # other suspicious


#hba1 <- c(46,88,14,42,79,7,9,1,66,56,62)
#hba2 <- c(79,46,88,14,64,56,7,62,68,78)
#hbb <- c(7,46,14,79,88,56,62,74,10)

hba1_fact <- hba1$fact[hba1$ratio > 1.4]
hba2_fact <- hba2$fact[hba2$ratio > 1.4]
hbb_fact <- hbb$fact[hba2$ratio > 1.4]
erythrocyte_rescue <- sort(intersect(intersect(hba1_fact,hba2_fact),hbb_fact))
erythrocyte_rescue


#dims.use <- unique(sort(c(setdiff(1:90,dims.remove),erythrocyte_rescue)))
df_perc_mito[order(df_perc_mito$cr,decreasing = T),]

dims.remove <- c(c(16,58,59,87),4,88,18,8)
dims.use <- setdiff(1:90,dims.remove)
lig.tumor <- quantile_norm(lig.tumor, dims.use = dims.use, quantiles = 100, eps = 0.01)
lig.tumor <- rliger::runUMAP(lig.tumor, dims.use = dims.use, 
                     n_neighbors = 10, min_dist = 0.05) #distance = "cosine", 

p <- plotByDatasetAndCluster(lig.tumor, return.plots = T)
p[[1]] + p[[2]]

m <- plotClusterFactors(object = lig.tumor,use.aligned = T,return.data = T)
rliger::plotWordClouds(lig.tumor)

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
rliger::plotGene(lig.tumor,"MKI67",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

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

rliger::plotGene(lig.tumor,"HBA1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"HBA2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"HBB",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"JCHAIN",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD19",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"MS4A1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"CD14",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"CD69",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"ACTA2",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"ALB",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"VIM",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"PDGFRA",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotGene(lig.tumor,"PECAM1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor,"CD34",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

lig.tumor <- louvainCluster(lig.tumor, resolution = 1, dims.use = dims.use)
#saveRDS(lig.tumor, paste0(outdir, "lig.tumor.rda"))

p.louvain  <- plotByDatasetAndCluster(lig.tumor, 
                                      axis.labels = c("UMAP 1", "UMAP 2"),
                                      return.plots = T)

p.louvain[[2]]

# Filter cells -----------------------------------------------------------------

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


lig.tumor@cell.data[["percent_mito"]] <- getProportionMito.human(lig.tumor)
lig.tumor@cell.data[["percent_ribo"]] <- getProportionRibo.human(lig.tumor)

#rliger::plotClusterFactors(lig.tumor,use.aligned = T)
rliger::plotClusterProportions(lig.tumor)

rliger::plotFeature(object = lig.tumor,feature = "percent_mito", by.dataset = F)
rliger::plotFeature(object = lig.tumor,feature = "percent_ribo", by.dataset = F)
rliger::plotFeature(object = lig.tumor,feature = "nGene", by.dataset = F)

df_perc_mito <- cor(lig.tumor@H.norm,lig.tumor@cell.data$percent_mito)
df_perc_mito <- data.frame(fact=1:90, cr=df_perc_mito)
df_perc_mito <- df_perc_mito[order(df_perc_mito$cr,decreasing = T),]

df_perc_ribo <- cor(lig.tumor@H.norm,lig.tumor@cell.data$percent_ribo)
df_perc_ribo <- data.frame(fact=1:90, cr=df_perc_ribo)
df_perc_ribo <- df_perc_ribo[order(df_perc_ribo$cr,decreasing = T),]

hist(lig.tumor@cell.data$percent_mito)
hist(lig.tumor@cell.data$percent_ribo)

ggplot(lig.tumor@cell.data,aes(x=percent_mito,y=percent_mito,color=dataset)) + geom_point() + theme_classic() + vertical_xlabels
ggplot(lig.tumor@cell.data,aes(x=percent_ribo,y=percent_mito,color=dataset)) + geom_point() + theme_classic() + vertical_xlabels
ggplot(lig.tumor@cell.data,aes(x=nGene,y=percent_mito,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

gene_val <- getGeneValues(lig.tumor@norm.data,gene = "HBA1")
stopifnot(identical(names(gene_val),rownames(lig.tumor@cell.data)))
lig.tumor@cell.data <- cbind(lig.tumor@cell.data,gene_val)

ggplot(lig.tumor@cell.data,aes(x=nGene,y=percent_mito,color=gene_val > 0)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

#thresh_nGene <- 500 | gene_val_HBA1 > 0

df_mito_clusters <- data.frame(percent_mito = lig.tumor@cell.data$percent_mito,
           cluster=lig.tumor@clusters)

df_mito_clusters_mn <- aggregate(df_mito_clusters$percent_mito,
          by=list(df_mito_clusters$cluster),FUN=mean)
colnames(df_mito_clusters_mn) <- c("cluster","mean")
df_mito_clusters_sd <- aggregate(df_mito_clusters$percent_mito,
          by=list(df_mito_clusters$cluster),FUN=sd)
colnames(df_mito_clusters_sd) <- c("cluster","sd")
df_mito_clusters_aggr <- merge(df_mito_clusters_mn,df_mito_clusters_sd,
                               by="cluster")

ggplot(df_mito_clusters_aggr,aes(x=reorder(cluster,-mean),y=mean)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  theme_classic() + vertical_xlabels

thresh_mito <- 0.3

df_ribo_clusters <- data.frame(percent_ribo = lig.tumor@cell.data$percent_ribo,
                               cluster=lig.tumor@clusters)

df_ribo_clusters_mn <- aggregate(df_ribo_clusters$percent_ribo,
                                 by=list(df_ribo_clusters$cluster),FUN=mean)
colnames(df_ribo_clusters_mn) <- c("cluster","mean")
df_ribo_clusters_sd <- aggregate(df_ribo_clusters$percent_ribo,
                                 by=list(df_ribo_clusters$cluster),FUN=sd)
colnames(df_ribo_clusters_sd) <- c("cluster","sd")
df_ribo_clusters_aggr <- merge(df_ribo_clusters_mn,df_ribo_clusters_sd,
                               by="cluster")

ggplot(df_ribo_clusters_aggr,aes(x=reorder(cluster,-mean),y=mean)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  theme_classic() + vertical_xlabels

#thresh_ribo <- 0.3

#https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html
seu <- rliger::ligerToSeurat(lig.tumor)
sce <- as.SingleCellExperiment(seu)

all.genes <- Reduce(union, lapply(lig.tumor@raw.data, rownames))
mt.genes <- grep(pattern = "^MT-", x = all.genes, value = TRUE)

#sce <- addPerCellQCMetrics(sce, subsets = list(mito = mt.genes))
sce <- addPerCellQC(sce, subsets = list(mito = mt.genes),percent.top=100)
colData(sce)
plotColData(sce, x = "sum", y="detected", colour_by="orig.ident")
plotHighestExprs(sce, exprs_values = "counts")

sce <- logNormCounts(sce)
vars <- getVarianceExplained(sce,
  variables=c("nCount_RNA", "nFeature_RNA","sum","subsets_mito_percent",
              "percent.top_100","detected"))
head(vars)

plotExplanatoryVariables(vars)

sce <- scater::runPCA(sce, use_coldata = TRUE,
              detect_outliers = TRUE)

plotReducedDim(sce, use_dimred="PCA_coldata", colour_by = "orig.ident")


stopifnot(identical(str_split_fixed(rownames(colData(sce)),"_",2)[,2],rownames(lig.tumor@cell.data)))
lig.tumor@cell.data$percent.top_100 <- colData(sce)[,"percent.top_100"]

rliger::plotFeature(object = lig.tumor, feature = "percent.top_100", by.dataset = F)


hist(lig.tumor@cell.data$percent.top_100) # percentage of library size counts coming from the top 100 overall genes

#thresh_percent_top_100 <- 80 | HBA1 > 0

ggplot(lig.tumor@cell.data,aes(x=percent_mito,y=percent.top_100,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels

ggplot(lig.tumor@cell.data,aes(x=nGene,y=percent.top_100,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor@cell.data,aes(x=nGene,y=percent.top_100,color=gene_val)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)


ggplot(lig.tumor@cell.data,aes(x=nGene,y=nUMI,color=gene_val)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)

ggplot(lig.tumor@cell.data,aes(x=nGene,y=nUMI,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 500)


ggplot(lig.tumor@cell.data,aes(x=percent_mito,y=percent_ribo,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 0.3) + 
  facet_wrap(~ dataset)


idx_keep <- (lig.tumor@cell.data[["nGene"]] > 500 & 
    lig.tumor@cell.data[["percent_mito"]] < 0.3 & 
    lig.tumor@cell.data[["percent.top_100"]] < 80) | 
  lig.tumor@cell.data[["gene_val"]] > 0
table(idx_keep)

cells_keep <- rownames(lig.tumor@cell.data[idx_keep,])

# cells_keep_list <- list()
# for (d in unique(lig.tumor@cell.data[idx_keep,]$dataset)){
#   cells_keep_list[[d]] <- rownames(lig.tumor@cell.data[idx_keep,][lig.tumor@cell.data[idx_keep,]$dataset == d,])
# }
# lig.tumor.filt <- optimizeSubset(lig.tumor, cell.subset = cells_keep_list, 
#                                  datasets.scale = names(cells_keep_list))

# lig.tumor.filt <- optimizeSubset(lig.tumor, cell.subset = cells_keep_list, 
#                                  datasets.scale = names(cells_keep_list))

lig.tumor.filt <- rliger::subsetLiger(lig.tumor,cells.use = cells_keep)
lig.tumor.filt <- optimizeALS(lig.tumor.filt, k = 90, lambda = 5, verbose = T)
#saveRDS(lig.tumor.filt, paste0(outdir, "lig.tumor.filt.rda"))
lig.tumor.filt <- readRDS(paste0(outdir, "lig.tumor.filt.rda"))

# Reprocess filtered dataset ---------------------------------------------------

dims.use <- 1:90
lig.tumor.filt <- quantile_norm(lig.tumor.filt, dims.use = dims.use, 
                                quantiles = 100, eps = 0.01)
lig.tumor.filt <- rliger::runUMAP(lig.tumor.filt, dims.use = dims.use, 
                     n_neighbors = 10, min_dist = 0.05)

p <- plotByDatasetAndCluster(lig.tumor.filt, return.plots = T)
p[[1]] + p[[2]]


ggplot(lig.tumor.filt@cell.data,aes(x=percent_mito,y=percent_ribo,color=dataset)) + 
  geom_point() + theme_classic() + vertical_xlabels + geom_vline(xintercept = 0.3) + 
  facet_wrap(~ dataset)

rliger::plotFeature(object = lig.tumor.filt, feature = "percent.top_100", by.dataset = F)
rliger::plotGene(lig.tumor.filt,"HBA1",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

rliger::plotFeature(object = lig.tumor.filt, feature = "percent_mito", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "percent_ribo", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "nGene", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "nUMI", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "dataset", by.dataset = F,discrete = T)

rliger::plotGene(lig.tumor.filt,"CD3D",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"MLANA",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")
rliger::plotGene(lig.tumor.filt,"CD14",axis.labels = c('UMAP 1', 'UMAP 2'), plot.by = "none")

# Add VDJ data -----------------------------------------------------------------

add_vdj_to_liger <- function(lig.tumor, vdj_paths){
  library(djvdj)
  vdj <- list()
  for (nm in names(fnames)) {
    vdj[[nm]] <- djvdj::import_vdj(input = NULL, vdj_dir = fnames[[nm]],
                            filter_chains = T, define_clonotypes = "cdr3_gene")
  }
  
  for (nm in names(fnames)) {
    rownames(vdj[[nm]]) <- paste0(nm,"_",str_split_fixed(rownames(vdj[[nm]]),"-",2)[,1])
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

lig.tumor <- add_vdj_to_liger(lig.tumor, vdj_paths = fnames)

lig.tumor.filt <- add_vdj_to_liger(lig.tumor.filt, vdj_paths = fnames)

# Doublets ---------------------------------------------------------------------
library(scDblFinder)

idx_keep <- lig.tumor@cell.data[["nGene"]] > 200 | lig.tumor@cell.data[["gene_val"]] > 0

cells_keep <- rownames(lig.tumor@cell.data[idx_keep,])

lig.tumor.filt <- rliger::subsetLiger(lig.tumor,cells.use = cells_keep)

seu <- rliger::ligerToSeurat(lig.tumor.filt)
sce <- as.SingleCellExperiment(seu)

sce.standard <- scDblFinder(sce, samples = "orig.ident")
sce.clusters <- scDblFinder(sce, samples = "orig.ident", clusters="ident")
sce.clusters_t <- scDblFinder(sce, samples = "orig.ident", clusters=T)

#knownDoublets <- names(which(colSums(exprs(sce)[c("MLANA","CD3D"),]>0) > 1))

mat <- exprs(sce)
colnames(mat) <- str_split_fixed(colnames(mat),"_",2)[,2]

stopifnot(identical(colnames(mat),rownames(lig.tumor.filt@cell.data)))


idx <- (!is.na(lig.tumor.filt@cell.data$cdr3) | colSums(mat[c("CD3D","CD4","CD8A","CD8B"),]>0) > 0) & 
  colSums(mat[c("MLANA","PMEL","TYR","HBA1","HBA2","HBB"),]>0) > 0

idx <- idx | ((!is.na(lig.tumor.filt@cell.data$cdr3) & 
  colSums(mat[c("CD14","CD19","MS4A1","JCHAIN"),]>0) > 0))

idx <- idx | ((colSums(mat[c("MLANA","PMEL","TYR"),]>0) > 0) & 
                (colSums(mat[c("CD14","CD19","MS4A1","JCHAIN"),]>0) > 0))

knownDoublets <- names(which(idx))


knownDoublets <- str_split_fixed(rownames(colData(sce)),"_",2)[,2] %in% knownDoublets
colData(sce) <- cbind(colData(sce),knownDoublets)

sce.standard.knownDoublets <- scDblFinder(sce, samples = "orig.ident", 
                            knownDoublets = "knownDoublets", knownUse = "discard")
sce.clusters.knownDoublets <- scDblFinder(sce, samples = "orig.ident", clusters="ident", 
                            knownDoublets = "knownDoublets", knownUse = "discard")
sce.clusters_t.knownDoublets <- scDblFinder(sce, samples = "orig.ident", clusters=T, 
                                          knownDoublets = "knownDoublets", 
                                          knownUse = "discard")

sce.standard.knownDoublets.positive <- scDblFinder(sce, samples = "orig.ident", 
                                          knownDoublets = "knownDoublets", knownUse = "positive")
sce.clusters.knownDoublets.positive <- scDblFinder(sce, samples = "orig.ident", clusters="ident", 
                                          knownDoublets = "knownDoublets", knownUse = "positive")
sce.clusters_t.knownDoublets.positive <- scDblFinder(sce, samples = "orig.ident", clusters=T, 
                                            knownDoublets = "knownDoublets", 
                                            knownUse = "positive")


table(colData(sce.standard)[["scDblFinder.class"]])
table(colData(sce.clusters)[["scDblFinder.class"]])
table(colData(sce.clusters_t)[["scDblFinder.class"]])
table(colData(sce.standard.knownDoublets)[["scDblFinder.class"]])
table(colData(sce.clusters.knownDoublets)[["scDblFinder.class"]])
table(colData(sce.clusters_t.knownDoublets)[["scDblFinder.class"]])

stopifnot(identical(str_split_fixed(rownames(colData(sce.standard)),"_",2)[,2],
          rownames(lig.tumor.filt@cell.data)))
stopifnot(identical(str_split_fixed(rownames(colData(sce.clusters)),"_",2)[,2],
                    rownames(lig.tumor.filt@cell.data)))
stopifnot(identical(str_split_fixed(rownames(colData(sce.clusters_t)),"_",2)[,2],
                    rownames(lig.tumor.filt@cell.data)))
stopifnot(identical(str_split_fixed(rownames(colData(sce.standard.knownDoublets)),"_",2)[,2],
                    rownames(lig.tumor.filt@cell.data)))
stopifnot(identical(str_split_fixed(rownames(colData(sce.clusters.knownDoublets)),"_",2)[,2],
                    rownames(lig.tumor.filt@cell.data)))
stopifnot(identical(str_split_fixed(rownames(colData(sce.clusters_t.knownDoublets)),"_",2)[,2],
                    rownames(lig.tumor.filt@cell.data)))

lig.tumor.filt@cell.data$scDblFinder.standard.score <- 
  colData(sce.standard)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.standard.class <- 
  colData(sce.standard)[["scDblFinder.class"]]

lig.tumor.filt@cell.data$scDblFinder.clusters.score <- 
  colData(sce.clusters)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.clusters.class <- 
  colData(sce.clusters)[["scDblFinder.class"]]

lig.tumor.filt@cell.data$scDblFinder.clusters_t.score <- 
  colData(sce.clusters_t)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.clusters_t.class <- 
  colData(sce.clusters_t)[["scDblFinder.class"]]

lig.tumor.filt@cell.data$scDblFinder.standard.knownDoublets.score <- 
  colData(sce.standard.knownDoublets)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.standard.knownDoublets.class <- 
  colData(sce.standard.knownDoublets)[["scDblFinder.class"]]

lig.tumor.filt@cell.data$scDblFinder.clusters.knownDoublets.score <- 
  colData(sce.clusters.knownDoublets)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.clusters.knownDoublets.class <- 
  colData(sce.clusters.knownDoublets)[["scDblFinder.class"]]

lig.tumor.filt@cell.data$scDblFinder.clusters_t.knownDoublets.score <- 
  colData(sce.clusters_t.knownDoublets)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.clusters_t.knownDoublets.class <- 
  colData(sce.clusters_t.knownDoublets)[["scDblFinder.class"]]

lig.tumor.filt@cell.data$scDblFinder.standard.knownDoublets.positive.score <- 
  colData(sce.standard.knownDoublets.positive)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.standard.knownDoublets.positive.class <- 
  colData(sce.standard.knownDoublets.positive)[["scDblFinder.class"]]

lig.tumor.filt@cell.data$scDblFinder.clusters.knownDoublets.positive.score <- 
  colData(sce.clusters.knownDoublets.positive)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.clusters.knownDoublets.positive.class <- 
  colData(sce.clusters.knownDoublets.positive)[["scDblFinder.class"]]

lig.tumor.filt@cell.data$scDblFinder.clusters_t.knownDoublets.positive.score <- 
  colData(sce.clusters_t.knownDoublets.positive)[["scDblFinder.score"]]
lig.tumor.filt@cell.data$scDblFinder.clusters_t.knownDoublets.positive.class <- 
  colData(sce.clusters_t.knownDoublets.positive)[["scDblFinder.class"]]


rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.standard.score", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.standard.class", by.dataset = F,discrete = T)

rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.clusters.score", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.clusters.class", by.dataset = F,discrete = T)

rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.clusters_t.score", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.clusters_t.class", by.dataset = F,discrete = T)

rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.standard.knownDoublets.score", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.standard.knownDoublets.class", by.dataset = F,discrete = T)

rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.clusters.knownDoublets.score", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.clusters.knownDoublets.class", by.dataset = F,discrete = T)

rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.clusters_t.knownDoublets.score", by.dataset = F)
rliger::plotFeature(object = lig.tumor.filt, feature = "scDblFinder.clusters_t.knownDoublets.class", by.dataset = F,discrete = T)

knownDoublets <- names(which(idx))
lig.tumor.filt@cell.data$knownDoublets <- rownames(lig.tumor.filt@cell.data) %in% knownDoublets

rliger::plotFeature(object = lig.tumor.filt, feature = "knownDoublets", 
                    by.dataset = F,discrete = T,pt.size = 0.01)

ggplot(lig.tumor.filt@cell.data,aes(x=scDblFinder.standard.score,y=scDblFinder.clusters.score)) + 
  geom_point()

ggplot(lig.tumor.filt@cell.data,
       aes(x=scDblFinder.standard.score,
           y=scDblFinder.standard.knownDoublets.score,
           color=paste0(scDblFinder.standard.class,"_",scDblFinder.standard.knownDoublets.class))) + 
  geom_point()

ggplot(lig.tumor.filt@cell.data,
       aes(x=scDblFinder.clusters.score,
           y=scDblFinder.clusters.knownDoublets.score,
           color=paste0(scDblFinder.clusters.class,"_",scDblFinder.clusters.knownDoublets.class))) + 
  geom_point()

ggplot(lig.tumor.filt@cell.data,
       aes(x=scDblFinder.clusters.score,
           y=scDblFinder.clusters_t.score,
           color=paste0(scDblFinder.clusters.class,"_",scDblFinder.clusters_t.class))) + 
  geom_point()

ggplot(lig.tumor.filt@cell.data,
       aes(x=scDblFinder.standard.knownDoublets.score,
           y=scDblFinder.clusters.knownDoublets.score,
           color=paste0(scDblFinder.standard.knownDoublets.class,"_",scDblFinder.clusters.knownDoublets.class))) + 
  geom_point()

ggplot(lig.tumor.filt@cell.data,
       aes(x=scDblFinder.clusters_t.knownDoublets.score,
           y=scDblFinder.clusters.knownDoublets.score,
           color=paste0(scDblFinder.clusters_t.knownDoublets.class,"_",scDblFinder.clusters.knownDoublets.class))) + 
  geom_point()

lig.tumor.filt@cell.data$doublet_agreement <- rowSums(lig.tumor.filt@cell.data[,grepl("^scDblFinder.*class$",colnames(lig.tumor.filt@cell.data)) & 
                                                                                 !grepl("positive",colnames(lig.tumor.filt@cell.data))] == "doublet")

max(lig.tumor.filt@cell.data$doublet_agreement)
rliger::plotFeature(object = lig.tumor.filt, feature = "doublet_agreement", by.dataset = F)

table(lig.tumor.filt@cell.data$doublet_agreement==6)

table(lig.tumor.filt@cell.data$doublet_agreement==6)[["TRUE"]]/nrow(lig.tumor.filt@cell.data)
table(lig.tumor.filt@cell.data$doublet_agreement>=1)[["TRUE"]]/nrow(lig.tumor.filt@cell.data)

table(lig.tumor.filt@cell.data$doublet_agreement >= 6,
      lig.tumor.filt@cell.data$dataset)

table(doublet_agreement=lig.tumor.filt@cell.data$scDblFinder.standard.class)["doublet"] / 
  nrow(lig.tumor.filt@cell.data)
table(doublet_agreement=lig.tumor.filt@cell.data$scDblFinder.standard.knownDoublets.class)["doublet"] / 
  nrow(lig.tumor.filt@cell.data)

table(doublet_agreement=lig.tumor.filt@cell.data$scDblFinder.clusters.class)["doublet"] / 
  nrow(lig.tumor.filt@cell.data)
table(doublet_agreement=lig.tumor.filt@cell.data$scDblFinder.clusters.knownDoublets.class)["doublet"] / 
  nrow(lig.tumor.filt@cell.data)

table(doublet_agreement=lig.tumor.filt@cell.data$scDblFinder.clusters_t.class)["doublet"] / 
  nrow(lig.tumor.filt@cell.data)
table(doublet_agreement=lig.tumor.filt@cell.data$scDblFinder.clusters_t.knownDoublets.class)["doublet"] / 
  nrow(lig.tumor.filt@cell.data)


df <- lig.tumor.filt@cell.data[,c(
  "scDblFinder.standard.class",
  "scDblFinder.standard.knownDoublets.class",
  "scDblFinder.clusters.class",
  "scDblFinder.clusters.knownDoublets.class",
  "scDblFinder.clusters_t.class",
  "scDblFinder.clusters_t.knownDoublets.class"
)]

df <- df=="doublet"
df <- ifelse(df==T,1,0)
library(pheatmap)
annot <- lig.tumor.filt@cell.data[,c("knownDoublets","scDblFinder.standard.class")]
annot$scDblFinder.standard.class <- NULL
x <- df[rowSums(df)>0,]
annot$dummy <- ""
annot <- annot[rownames(annot) %in% rownames(x),]
annot$dummy <- NULL
annot$knownDoublets <- factor(annot$knownDoublets)

pheatmap(x,show_rownames = F)
pheatmap(x,show_rownames = F,annotation_row = annot)


table(lig.tumor.filt@cell.data$doublet_agreement == 9,
      lig.tumor.filt@cell.data$dataset)

table(lig.tumor.filt@cell.data$knownDoublets,
      lig.tumor.filt@cell.data$scDblFinder.standard.class)

table(lig.tumor.filt@cell.data$knownDoublets,
      lig.tumor.filt@cell.data$scDblFinder.standard.knownDoublets.class)

table(lig.tumor.filt@cell.data$knownDoublets,
      lig.tumor.filt@cell.data$scDblFinder.clusters.class)

table(lig.tumor.filt@cell.data$knownDoublets,
      lig.tumor.filt@cell.data$scDblFinder.clusters.knownDoublets.class)

table(lig.tumor.filt@cell.data$knownDoublets,
      lig.tumor.filt@cell.data$scDblFinder.clusters_t.class)

table(lig.tumor.filt@cell.data$knownDoublets,
      lig.tumor.filt@cell.data$scDblFinder.clusters_t.knownDoublets.class)


# final filtering

idx_remove <- lig.tumor.filt@cell.data$knownDoublets | lig.tumor.filt@cell.data$scDblFinder.standard.class == "doublet"



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