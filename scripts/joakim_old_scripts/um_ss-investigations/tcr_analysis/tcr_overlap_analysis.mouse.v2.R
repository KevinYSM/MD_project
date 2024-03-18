library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(WriteXLS)
library(scater)
library(djvdj)
library(scDblFinder)
library(Matrix)
library(foreach)
library(doParallel)

source("~/proj/um_ss/Investigations/tcr_analysis/bin/TCRmatch_common.R")

vertical_xlabels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                     hjust=1))

#outdir <- "/data/proj/um_ss/Investigations/pdx_human/results/v1/"

outdir.old <- "~/proj/um_ss/Investigations/pdx_human/results/v1/"
outdir <- "~/proj/um_ss/Investigations/pdx_human/results/tcr/tcr_overlap_analysis/"

dir.create(outdir,recursive=T,showWarnings=F)

# Load data --------------------------------------------------------------------
seu.filt <- readRDS(file=paste0(outdir.old,"seu.filt.rda"))
seu.filt.cd8_t <- readRDS(file=paste0(outdir.old,"seu.filt.cd8_t.rda"))

seu.filt.biopsy <- readRDS(file=paste0("~/proj/um_ss/Investigations/seurat/results/v16/","seu.filt.rda"))
seu.filt.cd8_t.biopsy <- readRDS(file=paste0("~/proj/um_ss/Investigations/seurat/results/v16/","seu.filt.cd8_t.rda"))

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

# ------------------------------------------------------------------------------
tcrs.split_ab.pdx <- get_tcrs_split_ab_non_unique(seu.filt.cd8_t)
tcrs.split_ab.biopsy <- get_tcrs_split_ab_non_unique(seu.filt.cd8_t.biopsy)
tcrs.split_ab.til <- get_tcrs_split_ab_non_unique(dat.til)

# No matches between biopsy and PDX on basis of any CDR3b sequence
common <- intersect(unique(unlist(strsplit(seu.filt.cd8_t.biopsy@meta.data$cdr3,";"))),
                    unique(unlist(strsplit(seu.filt.cd8_t@meta.data$cdr3,";"))))
common <- common[!is.na(common)]
common[1];seu.filt.cd8_t.biopsy@meta.data[grep(common[1],seu.filt.cd8_t.biopsy@meta.data$cdr3),c("chains","cdr3")]
common[2];seu.filt.cd8_t.biopsy@meta.data[grep(common[2],seu.filt.cd8_t.biopsy@meta.data$cdr3),c("chains","cdr3")]
common[3];seu.filt.cd8_t.biopsy@meta.data[grep(common[3],seu.filt.cd8_t.biopsy@meta.data$cdr3),c("chains","cdr3")]

common[1];seu.filt.cd8_t@meta.data[grep(common[1],seu.filt.cd8_t@meta.data$cdr3),c("chains","cdr3")]
common[2];seu.filt.cd8_t@meta.data[grep(common[2],seu.filt.cd8_t@meta.data$cdr3),c("chains","cdr3")]
common[3];seu.filt.cd8_t@meta.data[grep(common[3],seu.filt.cd8_t@meta.data$cdr3),c("chains","cdr3")]

# meta.pdx <- seu.filt.cd8_t@meta.data
# meta.pdx.pair_ab <- meta.pdx[which(meta.pdx$chains=="TRA;TRB"),]
# meta.pdx.pair_ab$cdr3a <- str_split_fixed(meta.pdx.pair_ab$cdr3,";",2)[,1]
# meta.pdx.pair_ab$cdr3b <- str_split_fixed(meta.pdx.pair_ab$cdr3,";",2)[,2]
# meta.biopsy <- seu.filt.cd8_t.biopsy@meta.data
# 
# beta <- list()
# for (cdr in common){
#   beta[[cdr]] <- unique(meta.pdx.pair_ab[meta.pdx.pair_ab$cdr3a==cdr,]$cdr3b)
# }
# tcrs.split_ab.biopsy[tcrs.split_ab.biopsy$CDR3b %in% as.character(unlist(beta)),]

tcrs.split_ab.biopsy$TRB_clonotype <- paste0(tcrs.split_ab.biopsy$CDR3b,":",
                                             tcrs.split_ab.biopsy$TRBV,":",
                                             tcrs.split_ab.biopsy$TRBJ)
tcrs.split_ab.pdx$TRB_clonotype <- paste0(tcrs.split_ab.pdx$CDR3b,":",
                                             tcrs.split_ab.pdx$TRBV,":",
                                             tcrs.split_ab.pdx$TRBJ)

cells_matching_biopsy <- rownames(tcrs.split_ab.biopsy)[
  which(tcrs.split_ab.biopsy$TRB_clonotype %in% tcrs.split_ab.pdx$TRB_clonotype)]

common <- intersect(unique(unlist(strsplit(dat.til@meta.data$cdr3,";"))),
                    unique(unlist(strsplit(seu.filt.cd8_t@meta.data$cdr3,";"))))
common <- common[!is.na(common)]

common[1];dat.til@meta.data[grep(common[1],dat.til@meta.data$cdr3),c("chains","cdr3")]
common[2];dat.til@meta.data[grep(common[2],dat.til@meta.data$cdr3),c("chains","cdr3")]
common[3];dat.til@meta.data[grep(common[3],dat.til@meta.data$cdr3),c("chains","cdr3")]

tcrs.split_ab.til$TRB_clonotype <- paste0(tcrs.split_ab.til$CDR3b,":",
                                          tcrs.split_ab.til$TRBV,":",
                                          tcrs.split_ab.til$TRBJ)

cells_matching_til <- rownames(tcrs.split_ab.til)[
  which(tcrs.split_ab.til$TRB_clonotype %in% tcrs.split_ab.pdx$TRB_clonotype)]

# Only UM22 TILs match any mouse TCR
for (i in 1:length(tils.unintegrated)){
  cells <- rownames(tils.unintegrated[[i]]@meta.data)
  cells <- setNames(rep(F,length(cells)),cells)
  cells[names(cells) %in% cells_matching_til] <- T
  
  tils.unintegrated[[i]] <- AddMetaData(tils.unintegrated[[i]], metadata = cells,
                                        col.name = "tcr_matching_any_til")
}

p_um22 <- DimPlot(tils.unintegrated[["UM22"]],
                  group.by = "tcr_matching_any_til", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM22") + theme_void() + 
  theme(legend.position="none")

pdf(file = paste0(outdir,"til_tcr_matching_any_pdx.pdf"),width = 5, height = 5)
p_um22
dev.off()

# FeaturePlot(tils.unintegrated[["UM22"]],features = c("CD4","CD8A"),split.by = "tcr_matching_any_til")
# 
# x <- subset(tils.unintegrated[["UM22"]], (CD8A > 0 | CD8B > 0) & CD4 == 0)
# dim(x)
# FeaturePlot(x,features = c("CD4","CD8A"),split.by = "tcr_matching_any_til")

# res <- FindMarkers(x, ident.1 = T, ident.2 = F, group.by = "tcr_matching_any_til")
# res <- res[order(res$avg_log2FC,decreasing = T),]
# res <- res[res$p_val_adj < 0.05,]

# Map Takara to PDX ------------------------------------------------------------
vdj.takara <- read_vdj.takara()
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

tcrs.split_ab.pdx.takara.merged <- merge(
  tcrs.split_ab.pdx, 
  unique(tcrs.split_ab.takara[,c("TRB_clonotype","vasu_phenotype")]), 
  by = "TRB_clonotype")

# No overlap with Takara
intersect(tcrs.split_ab.takara$TRB_clonotype, tcrs.split_ab.pdx$TRB_clonotype)
intersect(tcrs.split_ab.takara$CDR3b, tcrs.split_ab.pdx$CDR3b)

# Two CDRa overlap
common <- intersect(vdj.takara[vdj.takara$chains=="TRA",]$cdr3,
                    unique(unlist(strsplit(seu.filt.cd8_t@meta.data$cdr3,";"))))

tmp1 <- seu.filt.cd8_t@meta.data[grep(common[1],seu.filt.cd8_t@meta.data$cdr3),c("chains","cdr3")]
tmp2 <- seu.filt.cd8_t@meta.data[grep(common[2],seu.filt.cd8_t@meta.data$cdr3),c("chains","cdr3")]

paired_cdr3b <- list()
for (cdr3a in common){
  tmp <- seu.filt.cd8_t@meta.data[grep(cdr3a,seu.filt.cd8_t@meta.data$cdr3),c("chains","cdr3")]
  idx <- which(unlist(lapply(strsplit(tmp$chains,";"),function(x){length(x)==2 && x[2]=="TRB"})))
  paired_cdr3b[[cdr3a]] <- unique(str_split_fixed(tmp[idx,]$cdr3,";",2)[,2])
}

# The corresponding CDR3b chain was not found in the Takara data
tcrs.split_ab.takara[tcrs.split_ab.takara$CDR3b %in% unlist(paired_cdr3b),]

# -------------------

# Overlap of TCRs between PDX groups 
tcrs.split_ab.pdx <- merge(tcrs.split_ab.pdx,seu.filt.cd8_t@meta.data, by="row.names",all.x=T)
unique_groups <- sort(unique(tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection))
res <- matrix(NA,
              nrow=length(unique_groups),
              ncol=length(unique_groups),
              dimnames = list(unique_groups, unique_groups))
for (group_1 in unique_groups){
  for (group_2 in unique_groups){
    res[group_1,group_2] <- length(intersect(
      tcrs.split_ab.pdx$TRB_clonotype[tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection == group_1],
      tcrs.split_ab.pdx$TRB_clonotype[tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection == group_2]
    ))
  }
}

pdf(file = paste0(outdir, "TCR_overlap_PDX_groups.pdf"), width =2.6, height = 2.22)
pheatmap(res, cluster_rows = F, cluster_cols = F, 
         display_numbers = res, border_color = NA,
         show_rownames = T,
         show_colnames = T, na_col = "white", 
         colorRampPalette(c("white", "red"))(100))
dev.off()

# Percentage of unique clonotypes per group ------------------------------------
# 
# tcrb.nonunique <- list()
# for (group in unique_groups){
#   tcrb.nonunique[[group]] <- unique(tcrs.split_ab.pdx$TRB_clonotype[
#     tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection == group])
# }
# 
# tcrs.split_ab.til$TRB_clonotype <- paste0(tcrs.split_ab.til$CDR3b,":",
#                                              tcrs.split_ab.til$TRBV,":",
#                                              tcrs.split_ab.til$TRBJ)
# 
# cells_til_nonunique_clonotypes <- list()
# for (group in unique_groups){
#   cells_til_nonunique_clonotypes[[group]] <- rownames(tcrs.split_ab.til)[
#     tcrs.split_ab.til$TRB_clonotype %in% tcrb.nonunique[[group]]]
# }
# 
# 
# unique_cell_id <- unique(unlist(cells_til_nonunique_clonotypes))
# 
# matches <- list()
# for (cell_id in unique_cell_id){
#   matches[[cell_id]] <- names(which(unlist(lapply(cells_til_nonunique_clonotypes,
#                                                   function(x){cell_id %in% x}))))
# }
# matches <- lapply(matches,function(x){
#   if (length(x)>1){
#     x <- paste0(x,collapse = ", ")
#   } 
#   return(x)
# })
# matches <- unlist(matches)
# 
# 
# sname <- "UM22"
# cells <- rownames(tils.unintegrated[[sname]]@meta.data)
# cells <- setNames(rep(NA,length(cells)),cells)
# 
# cells[names(matches)] <- matches
# 
# 
# tils.unintegrated[[sname]] <- AddMetaData(tils.unintegrated[[sname]], metadata = cells,
#                                       col.name = "PDX_group_nonunique_TCRb")
# 
# tils.unintegrated[[sname]]@meta.data$PDX_group_nonunique_TCRb_Liver_MART1 <- NA
# tils.unintegrated[[sname]]@meta.data$PDX_group_nonunique_TCRb_Liver_MART1[
#   grep("Liver_MART1",tils.unintegrated[[sname]]@meta.data$PDX_group_nonunique_TCRb)] <- "Liver_MART1"
# tils.unintegrated[[sname]]@meta.data$PDX_group_nonunique_TCRb_Liver_MART1[
#   ! grepl("Liver_MART1",tils.unintegrated[[sname]]@meta.data$PDX_group_nonunique_TCRb) & 
#     !is.na(tils.unintegrated[[sname]]@meta.data$PDX_group_nonunique_TCRb)] <- "Other"
# 
# DimPlot(tils.unintegrated[[sname]], group.by = "PDX_group_nonunique_TCRb",order = T)
# DimPlot(tils.unintegrated[[sname]], group.by = "PDX_group_nonunique_TCRb_Liver_MART1",order = T)
# 
# 
# 
# UM22.PDX_group_nonunique_TCRb_Liver_MART1 <- subset(tils.unintegrated[[sname]],
#                                                     cells=rownames(tils.unintegrated[[sname]]@meta.data)[
#                                                       !is.na(tils.unintegrated[[sname]]@meta.data$PDX_group_nonunique_TCRb_Liver_MART1)])
# DimPlot(UM22.PDX_group_nonunique_TCRb_Liver_MART1,group.by = "PDX_group_nonunique_TCRb_Liver_MART1",order = T)
# 
# res <- FindMarkers(UM22.PDX_group_nonunique_TCRb_Liver_MART1,
#                    ident.1 = "Liver_MART1",
#                    ident.2 = "Other",
#                    group.by = "PDX_group_nonunique_TCRb_Liver_MART1")
#res <- res[order(res$avg_log2FC,decreasing = T),]
#res <- res[res$p_val_adj < 0.05,]


# res.specific <- matrix(NA,
#               nrow=length(unique(tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection)),
#               ncol=1,
#               dimnames = list(unique(tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection),
#                               "Others"))
# for (group_1 in unique(tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection)){
#   #for (group_2 in unique(tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection)){
#     res.specific[group_1,] <- length(setdiff(
#       tcrs.split_ab.pdx$TRB_clonotype[tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection==group_1],
#       tcrs.split_ab.pdx$TRB_clonotype[tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection!=group_1]
#     ))/length(unique(tcrs.split_ab.pdx$TRB_clonotype))
#   #}
# }
# 
# # Unique clonotypes in each group
# unique_clonotypes <- list()
# for (group_1 in unique(tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection)){
#   unique_clonotypes[[group_1]] <- setdiff(
#     tcrs.split_ab.pdx$TRB_clonotype[tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection==group_1],
#     tcrs.split_ab.pdx$TRB_clonotype[tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection!=group_1]
#   )
# }
# 
# tcrs.split_ab.til$TRB_clonotype <- paste0(tcrs.split_ab.til$CDR3b,":",
#                                              tcrs.split_ab.til$TRBV,":",
#                                              tcrs.split_ab.til$TRBJ)
# 
# cells_til_unique_clonotypes <- list()
# for (group in unique(tcrs.split_ab.pdx$Mouse_biopsy_site_and_TIL_selection)){
#   cells_til_unique_clonotypes[[group]] <- rownames(tcrs.split_ab.til)[
#     tcrs.split_ab.til$TRB_clonotype %in% unique_clonotypes[[group]]]  
# }
# 
# sname <- "UM22"
# cells <- rownames(tils.unintegrated[[sname]]@meta.data)
# cells <- setNames(rep(NA,length(cells)),cells)
# for (group in names(cells_til_unique_clonotypes)){
#   cells[names(cells) %in% cells_til_unique_clonotypes[[group]]] <- group
# }
# 
# tils.unintegrated[[sname]] <- AddMetaData(tils.unintegrated[[sname]], metadata = cells,
#                                       col.name = "PDX_group_unique_TCRb")
# 
# DimPlot(tils.unintegrated[[sname]], group.by = "PDX_group_unique_TCRb",order = T)
# 
# UM22.PDX_group_unique_TCRb <- subset(tils.unintegrated[[sname]],cells=names(cells)[!is.na(cells)])
# 
# FeaturePlot(UM22.PDX_group_unique_TCRb, features = c("CD4","CD8A", "NCAM1")) # None express CD4 or NCAM1
# 
# res <- FindMarkers(UM22.PDX_group_unique_TCRb, 
#                    ident.1 = "Liver_MART1", 
#                    ident.2 = "Liver_Bulk", 
#                    group.by = "PDX_group_unique_TCRb", 
#                    min.cells.group = 2)
# res <- res[order(res$avg_log2FC,decreasing = T),]
# res <- res[res$p_val_adj < 0.05,]
# 
# FeaturePlot(UM22.PDX_group_unique_TCRb, 
#             features = c("DUSP4","KRT86", "PTPN22","LAYN"),
#             ncol = 3,order = T,split.by = "PDX_group_unique_TCRb")
# 
# # rownames(tils.unintegrated[[sname]]@meta.data)[
# #   unlist(lapply(tils.unintegrated[[sname]]@meta.data$cdr3,function(x){any(unlist(strsplit(x,";")) %in% 
# #     str_split_fixed(unique_clonotypes$Liver_MART1,":",3)[,1])}))]
# 
