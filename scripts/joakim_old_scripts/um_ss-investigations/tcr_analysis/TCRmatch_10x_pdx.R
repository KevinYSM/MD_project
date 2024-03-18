source("~/proj/um_ss/Investigations/tcr_analysis/bin/TCRmatch_common.R")

outdir <- "~/proj/um_ss/Investigations/pdx_human/results/tcr/TCRmatch/TCRmatch_10x_pdx/"
dir.create(outdir, showWarnings = F, recursive = T)

# Functions for plotting -------------------------------------------------------
add_metadata <- function(tcrmatch_output, seu.filt){
  tcrs.split_ab <- get_tcrs_split_ab_non_unique(seu.filt)
  
  tcrs.split_ab$tcrmatch_input <- ""
  tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$CDR3b
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "^C", replacement = "", 
                                       tcrs.split_ab$tcrmatch_input)
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "[FW]$", replacement = "", 
                                       tcrs.split_ab$tcrmatch_input)
  tcrs.split_ab$cell <- rownames(tcrs.split_ab)
  
  tcrs.split_ab.merged <- merge(tcrs.split_ab, tcrmatch_output, 
                                by.x="tcrmatch_input", 
                                by.y="input_sequence")#, all.x=T)
  tcrs.split_ab.merged$epitope[is.na(tcrs.split_ab.merged$epitope)] <- ""
  tcrs.split_ab.merged$antigen[is.na(tcrs.split_ab.merged$antigen)] <- ""
  tcrs.split_ab.merged$organism[is.na(tcrs.split_ab.merged$organism)] <- ""
  
  df <- data.frame(celltype = seu.filt.cd8_t@active.ident, 
                   TIL_selection = seu.filt.cd8_t@meta.data$TIL_selection, 
                   Mouse_biopsy_site_and_TIL_selection = seu.filt.cd8_t@meta.data$Mouse_biopsy_site_and_TIL_selection)
  #colnames(df) <- "celltype"
  tcrs.split_ab.merged <- merge(tcrs.split_ab.merged,df,by.x="cell",by.y="row.names")
  
  tcrs.split_ab.merged$celltype <- as.character(tcrs.split_ab.merged$celltype)
  tcrs.split_ab.merged$antigen[tcrs.split_ab.merged$antigen==""] <- "Unknown protein"
  tcrs.split_ab.merged$organism[tcrs.split_ab.merged$organism==""] <- "Unknown species"
  
  # Remove [organism] from antigen names
  tcrs.split_ab.merged$antigen <- unlist(lapply(tcrs.split_ab.merged$antigen,function(x){
    if (any(grepl("[",x,fixed=T))){
      x <- unlist(strsplit(x,";"))
      x <- paste0(unlist(lapply(x,subfun)),collapse=";")
    }
    return(x)
  }))
  
  tcrs.split_ab.merged$subject <- str_split_fixed(tcrs.split_ab.merged$`subject:condition`,":",2)[,1]
  tcrs.split_ab.merged$condition <- str_split_fixed(tcrs.split_ab.merged$`subject:condition`,":",2)[,2]
  
  tcrs.split_ab.merged <- rename_organisms(tcrs.split_ab.merged)
  tcrs.split_ab.merged <- rename_antigens(tcrs.split_ab.merged)
  
  tcrs.split_ab.merged <- tcrs.split_ab.merged[tcrs.split_ab.merged$antigen!="Unknown protein",]
  
  nms_dup <- names(which(table(unique(tcrs.split_ab.merged[,c("antigen","organism")])$antigen) > 1))
  if (length(nms_dup)>0){
    idx <- tcrs.split_ab.merged$antigen %in% nms_dup
    tcrs.split_ab.merged[idx,]$antigen <- paste0(
      tcrs.split_ab.merged[idx,]$antigen," (", tcrs.split_ab.merged[idx,]$organism,")")
  }
  
  return(tcrs.split_ab.merged)
}

add_unknown_tcrs <- function(tcrs.split_ab.merged, tcrs.split_ab, seu.filt.cd8_t){
  unknown_cells <- setdiff(tcrs.split_ab$cell, tcrs.split_ab.merged$cell)
  tcrs.split_ab <- tcrs.split_ab[tcrs.split_ab$cell %in% unknown_cells,]
  tcrs.split_ab$antigen <- "Unknown protein"
  tcrs.split_ab$organism <- "Unknown organism"
  rows_add <- tcrs.split_ab
  rows_add$subject <- str_split_fixed(rows_add$`subject:condition`,":",2)[,1]
  rows_add$condition <- str_split_fixed(rows_add$`subject:condition`,":",2)[,2]
  
  df <- data.frame(
    celltype = seu.filt.cd8_t@active.ident, 
    TIL_selection = seu.filt.cd8_t@meta.data$TIL_selection, 
    Mouse_biopsy_site_and_TIL_selection = seu.filt.cd8_t@meta.data$Mouse_biopsy_site_and_TIL_selection)
  #colnames(df) <- "celltype"
  rows_add <- merge(rows_add,df,by.x="cell",by.y="row.names",all.x=T,all.y = F)
  
  rows_add$celltype <- as.character(rows_add$celltype)
  
  cnames_missing <- setdiff(colnames(tcrs.split_ab.merged),colnames(rows_add))
  for (cname in cnames_missing){
    rows_add[[cname]] <- NA
  }
  rows_add <- rows_add[,colnames(tcrs.split_ab.merged)]
  tcrs.split_ab.merged <- rbind(tcrs.split_ab.merged,rows_add)
  
  return(tcrs.split_ab.merged)
}


# Write input files for TCRmatch -----------------------------------------------

# PDX
seu.filt.cd8_t <- readRDS(paste0("~/proj/um_ss/Investigations/pdx_human/results/v1/", "seu.filt.cd8_t.rda"))

# See, for trimming: https://www.frontiersin.org/articles/10.3389/fimmu.2021.640725/full
tcrs.split_ab <- get_tcrs_split_ab_non_unique(seu.filt.cd8_t)
tcrs.split_ab$tcrmatch_input <- ""
tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$CDR3b
tcrs.split_ab$tcrmatch_input <- gsub(pattern="^C", replacement="", 
                                     tcrs.split_ab$tcrmatch_input)
tcrs.split_ab$tcrmatch_input <- gsub(pattern="[FW]$", replacement="",
                                     tcrs.split_ab$tcrmatch_input)
tcrs.split_ab$cell <- rownames(tcrs.split_ab)

write.table(unique(tcrs.split_ab$tcrmatch_input), file = paste0(outdir,"seu.filt.cd8_t.tcrmatch_input.txt"),
            sep = "\t", col.names = F, quote = F, row.names = F)

# Read TCRmatch output ---------------------------------------------------------
tcrmatch_output <- read_tcrmatch(pth = paste0(outdir,"seu.filt.cd8_t.tcrmatch_output_2.txt"))

tcrs.split_ab.merged <- add_metadata(tcrmatch_output, seu.filt.cd8_t)
tcrs.split_ab.merged <- add_unknown_tcrs(tcrs.split_ab.merged, tcrs.split_ab, 
                                         seu.filt.cd8_t)

df <- tcrs.split_ab.merged

# Make Sankey plots ------------------------------------------------------------

df.sankey <- format_for_ggsankey(
  tcrs.split_ab.merged = df, 
  cols_include = c("cell", "subject", "condition", "antigen", "organism", "human_antigen"),
  cols_plot = c("subject", "condition", "organism", "human_antigen"))

df.sankey <- reorder_sankey_factors(df.sankey, df)

pdf(file = paste0(outdir,"sankey.pdf"), width = 10, height = 6.2)
plot_sankey(df.sankey, flow.alpha = .8)
dev.off()

# Compare TCR representation in different PDX groups ---------------------------
df$TCRb <- paste0(df$CDR3b,":",df$TRBV,":",df$TRBJ)
df$antigen_organism <- paste0(df$antigen,":",df$organism)
tab <- table(df$Mouse_biopsy_site_and_TIL_selection, df$TCRb)

dat.til <- readRDS(file = "~/proj/um_ss/Investigations/seurat/results/v16/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
dat.til <- subset(dat.til, subset = project.name %in% c("G18-023"))

tcrs.split_ab.til <- get_tcrs_split_ab_non_unique(dat.til)
tcrs.split_ab.til$TCRb <- paste0(tcrs.split_ab.til$CDR3b,":",tcrs.split_ab.til$TRBV,":",tcrs.split_ab.til$TRBJ)

#cells_til <- rownames(tcrs.split_ab.til[tcrs.split_ab.til$TCRb %in% df$TCRb,])
cells_pdx <- df[df$TCRb %in% tcrs.split_ab.til$TCRb,]$cell

df.pdx_and_til <- df[df$cell %in% cells_pdx,]

annot <- unique(data.frame(antigen_organism = df$antigen_organism, TCRb = df$TCRb))
annot$clonotype_in_til <- ifelse(annot$TCRb %in% tcrs.split_ab.til$TCRb,"Yes","No")
rownames(annot) <- annot$TCRb
annot$TCRb <- NULL

pdf(file = paste0(outdir, "TCR_vs_PDX_group.pdf"), width = 16, height = 5)
pheatmap(tab, border_color = NA, display_numbers = tab, 
         number_color = "#000000", annotation_col = annot, cluster_rows = T, 
         cluster_cols = T, na_col = "white", 
         colorRampPalette(c("white", "red"))(100))
dev.off()

tab <- table(df$celltype, df$TCRb)
pdf(file = paste0(outdir, "TCR_vs_cluster_group.pdf"), width = 16, height = 5)
pheatmap(tab, border_color = NA, display_numbers = tab, 
         number_color = "#000000", annotation_col = annot, cluster_rows = T, 
         cluster_cols = T, na_col = "white", 
         colorRampPalette(c("white", "red"))(100))
dev.off()

tab.cluster_vs_PDX_group <- table(
  df$celltype, df$Mouse_biopsy_site_and_TIL_selection)
pdf(file = paste0(outdir, "cluster_vs_PDX_group.pdf"), width = 3.3, height = 2.22)
pheatmap(tab.cluster_vs_PDX_group, border_color = NA, display_numbers = tab.cluster_vs_PDX_group, 
         number_color = "#000000", cluster_rows = T, 
         cluster_cols = F, na_col = "white", treeheight_row = 3,
         colorRampPalette(c("white", "red"))(100))
dev.off()

tab.cluster_vs_PDX_group.pdx_and_til <- table(
  df.pdx_and_til$celltype, df.pdx_and_til$Mouse_biopsy_site_and_TIL_selection)
pdf(file = paste0(outdir, "cluster_vs_PDX_group.pdx_and_til.pdf"), width = 3.3, height = 2.22)
pheatmap(tab.cluster_vs_PDX_group.pdx_and_til, border_color = NA, 
         display_numbers = tab.cluster_vs_PDX_group.pdx_and_til, 
         number_color = "#000000", cluster_rows = T, 
         cluster_cols = F, na_col = "white", treeheight_row = 3,
         colorRampPalette(c("white", "red"))(100))
dev.off()

stopifnot(all(tab.cluster_vs_PDX_group >= tab.cluster_vs_PDX_group.pdx_and_til))

tab.diff <- tab.cluster_vs_PDX_group - tab.cluster_vs_PDX_group.pdx_and_til

pdf(file = paste0(outdir, "cluster_vs_PDX_group.diff.pdf"), width = 3.3, height = 2.22)
pheatmap(tab.diff, border_color = NA, 
         display_numbers = tab.diff, 
         number_color = "#000000", cluster_rows = T, 
         cluster_cols = F, na_col = "white", treeheight_row = 3,
         colorRampPalette(c("white", "red"))(100))
dev.off()
