source("~/proj/um_ss/Investigations/tcr_analysis/bin/TCRmatch_common.R")

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/tcr/tcr_overlap_analysis/"
dir.create(outdir, recursive = T, showWarnings = F)

# Read data --------------------------------------------------------------------
outdir.old <- "~/proj/um_ss/Investigations/seurat/results/v16/"

dat.tumor <- readRDS(file=paste0(outdir.old,"seu.filt.rda"))
dat.tumor@meta.data$Sample.ID <- dat.tumor@meta.data$orig.ident
dat.tumor@meta.data$project.name <- "tumor"

seu.filt.cd8_t <- readRDS(file=paste0(outdir.old,"seu.filt.cd8_t.rda"))
seu.filt.cd8_t@meta.data$Sample.ID <- seu.filt.cd8_t@meta.data$orig.ident

dat.til <- readRDS(file = "~/proj/um_ss/Investigations/seurat/results/v16/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
dat.til <- subset(dat.til, subset = project.name %in% c("G18-023"))

tils.unintegrated <- readRDS(file = paste0(
  "~/proj/um_ss/Investigations/seurat/results/tils/", "tils.rda"))

names(tils.unintegrated) <- unlist(lapply(tils.unintegrated, function(x){
  unique(x@meta.data$UM.ID)
}))

tils.unintegrated <- lapply(tils.unintegrated, function(x){
  if (any(Idents(x) == "Gamma-delta T cells (TRGV9, TRDV2)")){
    x <- RenameIdents(x, `Gamma-delta T cells (TRGV9, TRDV2)` = "CD8 T cells")
  }
  x <- subset(x, cells = names(x@active.ident[x@active.ident != "CD8 T / NK cells"]))
})

vdj.takara <- read_vdj.takara()

tcrs.split_ab.tumor <- get_tcrs_split_ab_non_unique(dat.tumor)
tcrs.split_ab.til <- get_tcrs_split_ab_non_unique(dat.til)

# Create plots mapping biopsy to TIL and vice versa [Fig. 4d, 4g, S9c-d] -------

# TIL to biopsy

# cells_matching_biopsy <- rownames(tcrs.split_ab.til)[
#   which(tcrs.split_ab.til$CDR3b %in% tcrs.split_ab.tumor$CDR3b & 
#         tcrs.split_ab.til$TRBV %in% tcrs.split_ab.tumor$TRBV & 
#         tcrs.split_ab.til$TRBJ %in% tcrs.split_ab.tumor$TRBJ)]

tcrs.split_ab.til$TRB_clonotype <- paste0(tcrs.split_ab.til$CDR3b,":",
                                          tcrs.split_ab.til$TRBV,":",
                                          tcrs.split_ab.til$TRBV)

tcrs.split_ab.tumor$TRB_clonotype <- paste0(tcrs.split_ab.tumor$CDR3b,":",
                                          tcrs.split_ab.tumor$TRBV,":",
                                          tcrs.split_ab.tumor$TRBV)

cells_matching_biopsy <- rownames(tcrs.split_ab.til)[
  which(tcrs.split_ab.til$TRB_clonotype %in% tcrs.split_ab.tumor$TRB_clonotype)]

for (i in 1:length(tils.unintegrated)){
  cells <- rownames(tils.unintegrated[[i]]@meta.data)
  cells <- setNames(rep(F,length(cells)),cells)
  cells[names(cells) %in% cells_matching_biopsy] <- T

  tils.unintegrated[[i]] <- AddMetaData(tils.unintegrated[[i]], metadata = cells,
                                        col.name = "tcr_matching_any_biopsy")
}

p_um1 <- DimPlot(tils.unintegrated[["UM1"]],
                 group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                 cols = c("gray", "blue")) + ggtitle("UM1") + theme_void() + 
  theme(legend.position="none")
p_um2 <- DimPlot(tils.unintegrated[["UM2"]],
                 group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                 cols = c("gray", "blue")) + ggtitle("UM2") + theme_void() + 
  theme(legend.position="none")
p_um3 <- DimPlot(tils.unintegrated[["UM3"]],
                 group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                 cols = c("gray", "blue")) + ggtitle("UM3") + theme_void() + 
  theme(legend.position="none")
p_um9 <- DimPlot(tils.unintegrated[["UM9"]],
                 group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                 cols = c("gray", "blue")) + ggtitle("UM9") + theme_void() + 
  theme(legend.position="none")
p_um10 <- DimPlot(tils.unintegrated[["UM10"]],
                  group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM10") + theme_void() + 
  theme(legend.position="none")
p_um19 <- DimPlot(tils.unintegrated[["UM19"]],
                  group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM19") + theme_void() + 
  theme(legend.position="none")
p_um22 <- DimPlot(tils.unintegrated[["UM22"]],
                  group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM22") + theme_void() + 
  theme(legend.position="none")
p_um24 <- DimPlot(tils.unintegrated[["UM24"]],
                  group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM24") + theme_void() + 
  theme(legend.position="none")

# Fig. S9d
pdf(file = paste0(outdir,"til_tcr_matching_any_biopsy.pdf"),width = 36, height = 5)
p <- grid.arrange(grobs = c(list(p_um1),
                            list(p_um2),
                            list(p_um3),
                            list(p_um9),
                            list(p_um10),
                            list(p_um19),
                            list(p_um22),
                            list(p_um24)), ncol = 8, aseu.filt.cd8_table = FALSE)
dev.off()

# Biopsy to TIL, CD8
cells_matching_til <- rownames(tcrs.split_ab.tumor)[
  which(tcrs.split_ab.tumor$TRB_clonotype %in% tcrs.split_ab.til$TRB_clonotype)]

cells <- rownames(seu.filt.cd8_t@meta.data)
cells <- setNames(rep(F,length(cells)),cells)
cells[names(cells) %in% cells_matching_til] <- T

seu.filt.cd8_t <- AddMetaData(seu.filt.cd8_t, metadata = cells, 
                              col.name = "tcr_matching_any_til")

# Create plots mapping Takara to biopsy and TIL --------------------------------
vdj.takara.trb <- vdj.takara[vdj.takara$chains=="TRB",]
tcrs.split_ab.takara <- vdj.takara.trb
colnames(tcrs.split_ab.takara)[colnames(tcrs.split_ab.takara)=="cdr3"] <- "CDR3b"
colnames(tcrs.split_ab.takara)[colnames(tcrs.split_ab.takara)=="v_gene"] <- "TRBV"
colnames(tcrs.split_ab.takara)[colnames(tcrs.split_ab.takara)=="j_gene"] <- "TRBJ"
colnames(tcrs.split_ab.takara)[colnames(tcrs.split_ab.takara)=="d_gene"] <- "TRBD"
colnames(tcrs.split_ab.takara)[colnames(tcrs.split_ab.takara)=="c_gene"] <- "TRBC"
tcrs.split_ab.takara$Sample <- NULL
tcrs.split_ab.takara$chains <- NULL

tcrs.split_ab.takara$TRB_clonotype <- paste0(tcrs.split_ab.takara$CDR3b,":",
                                             tcrs.split_ab.takara$TRBV,":",
                                             tcrs.split_ab.takara$TRBJ)
tcrs.split_ab.tumor$TRB_clonotype <- paste0(tcrs.split_ab.tumor$CDR3b,":",
                                             tcrs.split_ab.tumor$TRBV,":",
                                             tcrs.split_ab.tumor$TRBJ)
tcrs.split_ab.tumor$cell <- rownames(tcrs.split_ab.tumor)

tcrs.split_ab.tumor.takara.merged <- merge(
  tcrs.split_ab.tumor, 
  unique(tcrs.split_ab.takara[,c("TRB_clonotype","vasu_phenotype")]), 
  by = "TRB_clonotype")

meta <- seu.filt.cd8_t@meta.data
meta <- merge(meta, tcrs.split_ab.tumor.takara.merged, by.x = "row.names", 
      by.y = "cell", all.x = T, all.y = F)
rownames(meta) <- meta$Row.names
meta <- meta[rownames(seu.filt.cd8_t@meta.data),]
seu.filt.cd8_t@meta.data <- meta
seu.filt.cd8_t@meta.data$vasu_phenotype[
  is.na(seu.filt.cd8_t@meta.data$vasu_phenotype)] <- "None"
seu.filt.cd8_t@meta.data$vasu_phenotype <- factor(
  seu.filt.cd8_t@meta.data$vasu_phenotype, 
  levels=unique(seu.filt.cd8_t@meta.data$vasu_phenotype))

tcrs.split_ab.til$TRB_clonotype <- paste0(tcrs.split_ab.til$CDR3b,":",
                                          tcrs.split_ab.til$TRBV,":",
                                          tcrs.split_ab.til$TRBJ)
tcrs.split_ab.til$cell <- rownames(tcrs.split_ab.til)
tcrs.split_ab.til.takara.merged <- merge(
  tcrs.split_ab.til, 
  unique(tcrs.split_ab.takara[,c("TRB_clonotype","vasu_phenotype")]), 
  by = "TRB_clonotype")

for (i in 1:length(tils.unintegrated)){
  meta <- tils.unintegrated[[i]]@meta.data
  meta <- merge(meta, tcrs.split_ab.til.takara.merged, by.x = "row.names", 
                by.y = "cell", all.x = T, all.y = F)
  rownames(meta) <- meta$Row.names
  meta <- meta[rownames(tils.unintegrated[[i]]@meta.data),]
  tils.unintegrated[[i]]@meta.data <- meta
  tils.unintegrated[[i]]@meta.data$vasu_phenotype[
    is.na(tils.unintegrated[[i]]@meta.data$vasu_phenotype)] <- "None"
  tils.unintegrated[[i]]@meta.data$vasu_phenotype <- factor(
    tils.unintegrated[[i]]@meta.data$vasu_phenotype, 
    levels = unique(tils.unintegrated[[i]]@meta.data$vasu_phenotype))
}

p_um1 <- DimPlot(tils.unintegrated[["UM1"]],
                 group.by = "vasu_phenotype", order = T, pt.size = 1) + 
  ggtitle("UM1") + theme_void()
p_um2 <- DimPlot(tils.unintegrated[["UM2"]],
                 group.by = "vasu_phenotype", order = T, pt.size = 1) + 
  ggtitle("UM2") + theme_void()
p_um3 <- DimPlot(tils.unintegrated[["UM3"]],
                 group.by = "vasu_phenotype", order = T, pt.size = 1) + 
  ggtitle("UM3") + theme_void()
p_um9 <- DimPlot(tils.unintegrated[["UM9"]],
                 group.by = "vasu_phenotype", order = T, pt.size = 1) + 
  ggtitle("UM9") + theme_void()
p_um10 <- DimPlot(tils.unintegrated[["UM10"]],
                  group.by = "vasu_phenotype", order = T, pt.size = 1) + 
  ggtitle("UM10") + theme_void()
p_um19 <- DimPlot(tils.unintegrated[["UM19"]],
                  group.by = "vasu_phenotype", order = T, pt.size = 1) + 
  ggtitle("UM19") + theme_void()
p_um22 <- DimPlot(tils.unintegrated[["UM22"]],
                  group.by = "vasu_phenotype", order = T, pt.size = 1) + 
  ggtitle("UM22") + theme_void()
p_um24 <- DimPlot(tils.unintegrated[["UM24"]],
                  group.by = "vasu_phenotype", order = T, pt.size = 1) + 
  ggtitle("UM24") + theme_void()

# Fig. S6g
pdf(file = paste0(outdir,"til_vasu_phenotype.pdf"),width = 35, height = 5)
p <- grid.arrange(grobs = c(list(p_um1),
                            list(p_um2),
                            list(p_um3),
                            list(p_um9),
                            list(p_um10),
                            list(p_um19),
                            list(p_um22),
                            list(p_um24)), ncol = 8, aseu.filt.cd8_table = FALSE)
dev.off()

pdf(file = paste0(outdir,"til_vasu_phenotype.UM1.UM9.pdf"), width = 10, height = 5)
p <- grid.arrange(grobs = c(list(p_um1),
                            list(p_um9)), ncol = 2, aseu.filt.cd8_table = FALSE)
dev.off()

# Create plots and statistical test for clusters vs Takara-biopsy vs TIL -------
get_proportions_til_vs_bulk_tcr <- function(seu.filt.cd8_t, var = "vasu_phenotype"){
  proportion_til <- list()
  proportion_vasu <- list()
  for (ident in unique(seu.filt.cd8_t@active.ident)){
    proportion_til[[ident]] <- sum(seu.filt.cd8_t@active.ident == ident &
                                     seu.filt.cd8_t@meta.data$tcr_matching_any_til) /
      sum(seu.filt.cd8_t@meta.data$tcr_matching_any_til)
    
    proportion_vasu[[ident]] <- sum(
      seu.filt.cd8_t@active.ident == ident & seu.filt.cd8_t@meta.data[[var]] != "None") /
      sum(seu.filt.cd8_t@meta.data[[var]] != "None")
  }
  
  df <- as.data.frame(cbind(unlist(proportion_til), unlist(proportion_vasu)))
  colnames(df) <- c("proportion_of_tils", "proportion_of_vasu")
  df$cluster <- rownames(df)
  df <- reshape2::melt(df)
  
  res <- test_cluster_overrepresantation_bulk_tcr_vs_til(seu.filt.cd8_t)
  
  df <- merge(df, res, by.x = "cluster", by.y = "row.names")
  clusters <- unique(df$cluster)
  df$is_max <- F
  for (cluster in clusters){
    df$is_max[df$cluster == cluster][
      which.max(df$value[df$cluster == cluster])] <- T
  }
  
  return(df)
}

test_cluster_overrepresantation_bulk_tcr_vs_til <- function(seu.filt.cd8_t){
  bt <- list()
  for (ident in unique(seu.filt.cd8_t@active.ident)){
    expected_probability <- sum(seu.filt.cd8_t@active.ident == ident & 
                                  seu.filt.cd8_t@meta.data$tcr_matching_any_til) /
      sum(seu.filt.cd8_t@meta.data$tcr_matching_any_til)
    n_rolls <- sum(seu.filt.cd8_t@meta.data$vasu_phenotype != "None")
    n_successes <- sum(seu.filt.cd8_t@meta.data$vasu_phenotype != "None" &
                         seu.filt.cd8_t@active.ident == ident)
    bt[[ident]] <- binom.test(n_successes, n_rolls, expected_probability,
                              alternative = "two.sided")
  }
  
  p_null <- unlist(lapply(bt, function(x){ as.numeric(x$null.value) }))
  p_estimate <- unlist(lapply(bt, function(x){ as.numeric(x$estimate) }))
  p <- unlist(lapply(bt, function(x){ x$p.value }))
  q <- p.adjust(p, method = "bonferroni") # probably not independent tests (the same dice rolls for all)
  res <- data.frame(p_null, p_estimate, p, q, vasu_greater = p_estimate > p_null,
                    sig = q < 0.05)
  return(res)
}

plot_cluster_overrepresentation <- function(df, 
  title = "Proportion of either that is present in a given cluster"){
  g <- ggplot(df, aes(x = cluster, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") + theme_classic() +
    vertical_xlabels + xlab(NULL) + ylab("Proportion") +
    geom_text(aes(label = ifelse(q < 0.05, "*", "")),
              data = subset(df, is_max), nudge_y = 0.02) +
    ggtitle(title)
  return(g)
}

get_proportions_bulk_tcr_vs_cluster_sample <- function(df, um_id){
  x <- table(df$vasu_phenotype[df$UM.ID == um_id],df$cluster[df$UM.ID == um_id])
  x <- t(t(x)/colSums(x))
  x <- as.data.frame.matrix(x)
  x$phenotype <- rownames(x)
  x <- reshape2::melt(x)
  x$UM.ID <- um_id
  return(x)
}

get_proportions_clusters_vs_bulk_tcr_sample <- function(
    seu.filt.cd8_t, um_id, lvls = sort(unique(levels(seu.filt.cd8_t@active.ident)))){
  
  tab <- table(seu.filt.cd8_t@meta.data$UM.ID,seu.filt.cd8_t@active.ident)
  tab <- as.data.frame.matrix(tab)
  tab.perc <- tab[um_id,]/sum(tab[um_id,])
  
  tab.vasu <- table(seu.filt.cd8_t@meta.data$vasu_phenotype[
    seu.filt.cd8_t@meta.data$UM.ID == um_id],
    seu.filt.cd8_t@active.ident[seu.filt.cd8_t@meta.data$UM.ID == um_id])
  
  tab.vasu.perc <- colSums(tab.vasu[c("MART1","CD3+41bb+"),]) / 
    sum(colSums(tab.vasu[c("MART1","CD3+41bb+"),]))
  
  tab.combined <- rbind(all_sample = tab.perc, 
                        all_vasu_matches_sample = tab.vasu.perc)
  tab.combined$categ <- rownames(tab.combined)
  tab.combined <- reshape2::melt(tab.combined)
  colnames(tab.combined) <- c("Category","Cluster","Proportion")
  
  missing_levels <- setdiff(lvls,as.character(tab.combined$Cluster))
  if (length(missing_levels)>0){
    df_missing_all_sample <- data.frame(Category = "all_sample", 
                                        Cluster = missing_levels, Proportion = 0)
    df_missing_all_vasu_matches_sample <- data.frame(Category = "all_vasu_matches_sample", 
                                                     Cluster = missing_levels, Proportion = 0)
    
    tab.combined <- rbind(tab.combined, 
                          df_missing_all_sample, 
                          df_missing_all_vasu_matches_sample)
  }
  
  tab.combined$Cluster <- as.character(tab.combined$Cluster)
  tab.combined <- tab.combined[order(tab.combined$Cluster),]
  tab.combined$Cluster <- factor(tab.combined$Cluster, 
                                 levels = lvls)
  
  return(tab.combined)
}

res <- test_cluster_overrepresantation_bulk_tcr_vs_til(seu.filt.cd8_t)

library(WriteXLS)
WriteXLS(res,
         ExcelFileName = paste0(outdir,"proportions_til_vs_vasu_matches_per_cluster.xlsx"),
         AdjWidth = T, AutoFilter = T, BoldHeaderRow = T, row.names = T)

# Plot
proportion_til <- list()
proportion_vasu <- list()
for (ident in unique(seu.filt.cd8_t@active.ident)){
  proportion_til[[ident]] <- sum(seu.filt.cd8_t@active.ident == ident &
                                   seu.filt.cd8_t@meta.data$tcr_matching_any_til) /
    sum(seu.filt.cd8_t@meta.data$tcr_matching_any_til)
  
  proportion_vasu[[ident]] <- sum(
    seu.filt.cd8_t@active.ident == ident & seu.filt.cd8_t@meta.data$vasu_phenotype != "None") /
    sum(seu.filt.cd8_t@meta.data$vasu_phenotype != "None")
}

df <- as.data.frame(cbind(unlist(proportion_til), unlist(proportion_vasu)))
colnames(df) <- c("proportion_of_tils", "proportion_of_vasu")
df$cluster <- rownames(df)
df <- reshape2::melt(df)

df <- merge(df, res, by.x = "cluster", by.y = "row.names")
clusters <- unique(df$cluster)
df$is_max <- F
for (cluster in clusters){
  df$is_max[df$cluster == cluster][
    which.max(df$value[df$cluster == cluster])] <- T
}

# Fig. 4d-f
vertical_xlabels <- theme(axis.text.x = element_text(
  angle = 90, vjust = 0.5, hjust=1))

p1 <- DimPlot(seu.filt.cd8_t, group.by = "tcr_matching_any_til", order = T, 
              cols = c("gray","blue"), pt.size = 1) + theme_void()
p2 <- DimPlot(seu.filt.cd8_t, group.by = "vasu_phenotype", order = T,
          cols = c("gray","blue","red","green"), pt.size = 1) + theme_void()
  
p3 <- plot_cluster_overrepresentation(df)

p1 + p2 + p3

ggsave(filename = paste0(outdir, "proportions_til_vs_vasu_matches_per_cluster.pdf"),
       width = 18, height = 7)

seu.filt.cd8_t <- add_anatomic_location_meta(
  seu.filt.cd8_t, pth = "~/proj/um_ss/Investigations/samples_10x.xlsx")
stopifnot(identical(seu.filt.cd8_t@meta.data$UM.ID.y, 
                    seu.filt.cd8_t@meta.data$UM.ID.x))
seu.filt.cd8_t@meta.data$UM.ID <- seu.filt.cd8_t@meta.data$UM.ID.x
seu.filt.cd8_t@meta.data$UM.ID.x <- NULL
seu.filt.cd8_t@meta.data$UM.ID.y <- NULL

df <- data.frame(vasu_phenotype = seu.filt.cd8_t@meta.data$vasu_phenotype,
                 cluster = seu.filt.cd8_t@active.ident,
                 Original_biopsy_tissue_site = 
                   seu.filt.cd8_t@meta.data$Original_biopsy_tissue_site,
                 UM.ID = seu.filt.cd8_t@meta.data$UM.ID)

x.UM1 <- get_proportions_bulk_tcr_vs_cluster_sample(df, "UM1")
x.UM9 <- get_proportions_bulk_tcr_vs_cluster_sample(df, "UM9")
x.UM46 <- get_proportions_bulk_tcr_vs_cluster_sample(df, "UM46")
x <- rbind(x.UM1,x.UM9,x.UM46)
x$UM.ID <- factor(x$UM.ID, levels = c("UM1","UM9","UM46"))
x$variable <- as.character(x$variable)
x$variable <- factor(x$variable,levels=sort(unique(x$variable)))

ggplot(x, aes(x = variable, y = value, fill = phenotype)) + 
  geom_bar(stat = "identity") + theme_classic() + facet_wrap(~ UM.ID, ncol = 1) + 
  vertical_xlabels + xlab(NULL)
ggsave(file = paste0(outdir, "bulk_tcr_matches_vs_cluster_per_sample.pdf"),width = 3.5,height = 5)

seu.filt.cd8_t@meta.data$vasu_phenotype_41bb <- ifelse(
  as.character(seu.filt.cd8_t@meta.data$vasu_phenotype) == "CD3+41bb+", "CD3+41bb+", "None")

seu.filt.cd8_t@meta.data$vasu_phenotype_MART1 <- ifelse(
  as.character(seu.filt.cd8_t@meta.data$vasu_phenotype) == "MART1", "MART1","None")

df.41bb <- get_proportions_til_vs_bulk_tcr(seu.filt.cd8_t, var = "vasu_phenotype_41bb")
df.MART1 <- get_proportions_til_vs_bulk_tcr(seu.filt.cd8_t, var = "vasu_phenotype_MART1")

df.41bb$variable <- as.character(df.41bb$variable)
df.41bb$variable[df.41bb$variable=="proportion_of_vasu"] <- "proportion_of_41bb"

df.MART1$variable <- as.character(df.MART1$variable)
df.MART1$variable[df.MART1$variable=="proportion_of_vasu"] <- "proportion_of_MART1"

df <- rbind(df.41bb[df.41bb$variable!="proportion_of_tils",],
      df.MART1[df.MART1$variable!="proportion_of_tils",])
df$q <- 1
plot_cluster_overrepresentation(df, title = "41bb_vs_MART1")
ggsave(file = paste0(outdir,"41bb_vs_MART1_matches_clusters.pdf"),width = 4,height = 3.5)

tab.UM1 <- get_proportions_clusters_vs_bulk_tcr_sample(
  seu.filt.cd8_t, um_id = "UM1", 
  lvls = sort(unique(levels(seu.filt.cd8_t@active.ident))))

tab.UM9 <- get_proportions_clusters_vs_bulk_tcr_sample(
  seu.filt.cd8_t, um_id = "UM9", 
  lvls = sort(unique(levels(seu.filt.cd8_t@active.ident))))

tab.UM46 <- get_proportions_clusters_vs_bulk_tcr_sample(
  seu.filt.cd8_t, um_id = "UM46", 
  lvls = sort(unique(levels(seu.filt.cd8_t@active.ident))))

tab.UM1$UM.ID <- "UM1"
tab.UM9$UM.ID <- "UM9"
tab.UM46$UM.ID <- "UM46"

tab.all.vasu <- rbind(tab.UM1[tab.UM1$Category != "all_sample",],
                 tab.UM9[tab.UM9$Category != "all_sample",],
                 tab.UM46[tab.UM46$Category != "all_sample",])
tab.all.vasu$UM.ID <- factor(tab.all.vasu$UM.ID, levels = c("UM1", "UM9", "UM46"))

tab.all.vasu$Cluster <- as.character(tab.all.vasu$Cluster)
tab.all.vasu$Cluster <- factor(tab.all.vasu$Cluster,levels=sort(unique(tab.all.vasu$Cluster),decreasing = T))

tab.all <- rbind(tab.UM1[tab.UM1$Category == "all_sample",],
                      tab.UM9[tab.UM9$Category == "all_sample",],
                      tab.UM46[tab.UM46$Category == "all_sample",])
tab.all$UM.ID <- factor(tab.all$UM.ID, levels = c("UM1", "UM9", "UM46"))

tab.all$Cluster <- as.character(tab.all$Cluster)
tab.all$Cluster <- factor(tab.all$Cluster,levels=sort(unique(tab.all$Cluster),decreasing = T))

stopifnot(identical(dim(tab.all.vasu),dim(tab.all)))
stopifnot(identical(tab.all.vasu$Cluster,tab.all$Cluster))
stopifnot(identical(tab.all.vasu$UM.ID,tab.all$UM.ID))

tab.diff <- tab.all.vasu
tab.diff$Proportion <- tab.all.vasu$Proportion - tab.all$Proportion

p1 <- ggplot(tab.all.vasu, aes(x = UM.ID, y = Cluster, fill = Proportion)) + geom_tile() + 
  theme_minimal() + xlab(NULL) + ylab(NULL) + ggtitle("CD3+4-1BB+ or MART1 match") + 
  scale_fill_gradient(low = "white", high = "red") + vertical_xlabels

p2 <- ggplot(tab.all, aes(x = UM.ID, y = Cluster, fill = Proportion)) + geom_tile() + 
  theme_minimal() + xlab(NULL) + ylab(NULL) + ggtitle("All biopsy cells") + 
  scale_fill_gradient(low = "white", high = "red") + vertical_xlabels

p3 <- ggplot(tab.diff, aes(x = UM.ID, y = Cluster, fill = Proportion)) + geom_tile() + 
  theme_minimal() + xlab(NULL) + ylab(NULL) + scale_fill_distiller(palette = "RdGy") + 
  ggtitle("Difference") + vertical_xlabels

p1 + p2 + p3
ggsave(filename = paste0(outdir,"tcr_bulk_matches_vs_all_biopsy_cells.pdf"),
       width = 11, height = 2.5)

# Fig. S9c
seu.filt.cd8_t@meta.data$UM.ID <- factor(seu.filt.cd8_t@meta.data$UM.ID, 
                              levels = unique(paste0("UM",sort(
                                as.numeric(gsub(pattern = "UM", replacement = "", 
                                                seu.filt.cd8_t@meta.data$UM.ID))))))

DimPlot(subset(seu.filt.cd8_t, UM.ID %in% c("UM1","UM2","UM3","UM9","UM19","UM24")), 
        group.by = "tcr_matching_any_til", order = T, split.by = "UM.ID",
        cols = c("gray","blue"), pt.size = 1) + theme_void()
ggsave(filename = paste0(outdir,"tcr_biopsy_matcting_til.pdf"),
       width = 26, height = 5)

DimPlot(seu.filt.cd8_t, 
        group.by = "tcr_matching_any_til", order = T, split.by = "UM.ID",
        cols = c("gray","blue"), pt.size = 1) + theme_void()
ggsave(filename = paste0(outdir,"tcr_biopsy_matcting_til.all.pdf"),
       width = 70, height = 5,limitsize = F)

# Create plots on overlaps between TIL and biopsy ------------------------------
tcrs.split_ab.tumor$UM.ID <- str_split_fixed(tcrs.split_ab.tumor$`subject:condition`,":",2)[,1]
tcrs.split_ab.til$UM.ID <- str_split_fixed(tcrs.split_ab.til$`subject:condition`,":",2)[,1]

stopifnot(identical(unique(as.character(dat.tumor@meta.data$UM.ID)), unique(tcrs.split_ab.tumor$UM.ID)))
stopifnot(identical(names(tils.unintegrated), unique(tcrs.split_ab.til$UM.ID)))

tcrs.split_ab.tumor$UM.ID <- paste0(tcrs.split_ab.tumor$UM.ID,"_Biopsy")
tcrs.split_ab.til$UM.ID <- paste0(tcrs.split_ab.til$UM.ID,"_TILs")

unique_um_id_tumor <- unique(tcrs.split_ab.tumor$UM.ID)
unique_um_id_til <- unique(tcrs.split_ab.til$UM.ID)

unique_um_id_tumor <- unique_um_id_tumor[
  order(as.numeric(gsub("_Biopsy","",gsub("UM","",unique_um_id_tumor))),
        decreasing = F)]
unique_um_id_til <- unique_um_id_til[
  order(as.numeric(gsub("_TILs","",gsub("UM","",unique_um_id_til))),
        decreasing = F)]

unique_um_id <- c(unique_um_id_tumor, unique_um_id_til)

tcrs.split_ab_tumor_til <- rbind(tcrs.split_ab.tumor,tcrs.split_ab.til)

overlap.tumor_til <- matrix(NA, nrow = length(unique_um_id), 
                            ncol = length(unique_um_id), 
                            dimnames = list(unique_um_id, unique_um_id))
for (um_id_1 in unique_um_id){
  for (um_id_2 in unique_um_id){
    overlap.tumor_til[um_id_1,um_id_2] <- length(intersect(
      tcrs.split_ab_tumor_til$TRB_clonotype[tcrs.split_ab_tumor_til$UM.ID==um_id_1],
      tcrs.split_ab_tumor_til$TRB_clonotype[tcrs.split_ab_tumor_til$UM.ID==um_id_2]))
  }
}

overlap.tumor_til[lower.tri(overlap.tumor_til)] <- NA
diag(overlap.tumor_til) <- NA

annot <- data.frame(sname = unique_um_id, 
                    type = str_split_fixed(unique_um_id,"_",2)[,2])
rownames(annot) <- annot$sname
annot$sname <- NULL

# Fig. 4g
pdf(file=paste0(outdir,"tcr_overlap.pdf"), width = 6, height = 5)
pheatmap(overlap.tumor_til, cluster_rows = F, cluster_cols = F, 
         display_numbers = overlap.tumor_til, border_color = NA,
         show_rownames = T, annotation_row = annot, annotation_col = annot,
         show_colnames = T, na_col = "white", colorRampPalette(c("white", "red"))(100))
dev.off()

# Create plots on overlaps between TIL and biopsy with respect to takara -------
TRB_clonotypes_takara <- tcrs.split_ab.takara$TRB_clonotype

overlap.tumor_til.takara <- matrix(NA, nrow = length(unique_um_id), 
                            ncol = length(unique_um_id), 
                            dimnames = list(unique_um_id, unique_um_id))
for (um_id_1 in unique_um_id){
  for (um_id_2 in unique_um_id){
    overlap.tumor_til.takara[um_id_1, um_id_2] <- length(
      intersect(
        intersect(
          tcrs.split_ab_tumor_til$TRB_clonotype[tcrs.split_ab_tumor_til$UM.ID == um_id_1],
          tcrs.split_ab_tumor_til$TRB_clonotype[tcrs.split_ab_tumor_til$UM.ID == um_id_2]),
        TRB_clonotypes_takara))
  }
}

overlap.tumor_til.takara[lower.tri(overlap.tumor_til.takara)] <- NA

annot <- data.frame(sname = unique_um_id, 
                    type = str_split_fixed(unique_um_id,"_",2)[,2])
rownames(annot) <- annot$sname
annot$sname <- NULL

# Fig. S9a
pdf(file=paste0(outdir,"tcr_overlap.takara.pdf"), width = 6, height = 5)
pheatmap(overlap.tumor_til.takara, cluster_rows = F, cluster_cols = F, 
         display_numbers = overlap.tumor_til.takara, border_color = NA,
         show_rownames = T, annotation_row = annot, annotation_col = annot,
         show_colnames = T, na_col = "white", colorRampPalette(c("white", "red"))(100))
dev.off()

unique_um_ids <- unique(tcrs.split_ab.takara$UM.ID)
unique_um_ids <- unique_um_ids[
  order(as.numeric(gsub("_TILs","",gsub("UM","",unique_um_ids))),
        decreasing = F)]

overlap.takara <- matrix(NA, nrow=length(unique_um_ids),
                         ncol = length(unique_um_ids), 
                         dimnames = list(unique_um_ids,unique_um_ids))
for (um_id_1 in unique_um_ids){
  for (um_id_2 in unique_um_ids){
    overlap.takara[um_id_1, um_id_2] <- length(intersect(
      tcrs.split_ab.takara$TRB_clonotype[tcrs.split_ab.takara$UM.ID == um_id_1],
      tcrs.split_ab.takara$TRB_clonotype[tcrs.split_ab.takara$UM.ID == um_id_2]))
  }
}

overlap.takara[lower.tri(overlap.takara)] <- NA

# Fig. S9b
pdf(file=paste0(outdir,"tcr_overlap.takara_only.pdf"), width = 3, height = 2.4)
pheatmap(overlap.takara, cluster_rows = F, cluster_cols = F, 
         display_numbers = overlap.takara, border_color = NA,
         show_rownames = T,
         show_colnames = T, na_col = "white", colorRampPalette(c("white", "red"))(100))
dev.off()