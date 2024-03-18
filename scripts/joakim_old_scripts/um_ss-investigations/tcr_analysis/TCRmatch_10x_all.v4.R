source("~/proj/um_ss/Investigations/tcr_analysis/bin/TCRmatch_common.R")

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/tcr/TCRmatch/TCRmatch_10x_all/"
dir.create(outdir, showWarnings = F, recursive = T)

# Functions for plotting -------------------------------------------------------
add_metadata <- function(tcrmatch_output, seu.filt, seu.filt.cd8_t){
  tcrs.split_ab <- get_tcrs_split_ab_non_unique(seu.filt)
  
  tcrs.split_ab$tcrmatch_input <- ""
  tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$CDR3b
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "^C", replacement = "", tcrs.split_ab$tcrmatch_input)
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "[FW]$", replacement = "", tcrs.split_ab$tcrmatch_input)
  
  tcrs.split_ab$cell <- rownames(tcrs.split_ab)
  
  tcrs.split_ab.merged <- merge(tcrs.split_ab, tcrmatch_output, 
                                by.x="tcrmatch_input", 
                                by.y="input_sequence")#, all.x=T)
  tcrs.split_ab.merged$epitope[is.na(tcrs.split_ab.merged$epitope)] <- ""
  tcrs.split_ab.merged$antigen[is.na(tcrs.split_ab.merged$antigen)] <- ""
  tcrs.split_ab.merged$organism[is.na(tcrs.split_ab.merged$organism)] <- ""
  
  df <- data.frame(seu.filt.cd8_t@active.ident)
  colnames(df) <- "celltype"
  tcrs.split_ab.merged <- merge(tcrs.split_ab.merged, df, by.x = "cell", 
                                by.y = "row.names")
  
  seu.filt.cd8_t <- add_anatomic_location_meta(
    seu.filt.cd8_t, pth = "~/proj/um_ss/Investigations/samples_10x.xlsx")
  
  tcrs.split_ab.merged <- merge(tcrs.split_ab.merged,
                                data.frame(cell = rownames(seu.filt.cd8_t@meta.data),
                                           tissue_site = seu.filt.cd8_t@meta.data$Original_biopsy_tissue_site), 
                                by = "cell")
  
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

add_metadata_til <- function(tcrmatch_output, dat.til){
  tcrs.split_ab <- get_tcrs_split_ab_non_unique(dat.til, condition = "TIL")
  
  tcrs.split_ab$tcrmatch_input <- ""
  tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$CDR3b
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "^C", replacement = "", tcrs.split_ab$tcrmatch_input)
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "[FW]$", replacement = "", tcrs.split_ab$tcrmatch_input)
  
  tcrs.split_ab$cell <- rownames(tcrs.split_ab)
  
  tcrs.split_ab.merged <- merge(tcrs.split_ab, tcrmatch_output, 
                                by.x="tcrmatch_input", 
                                by.y="input_sequence")#, all.x=T)
  
  tcrs.split_ab.merged$epitope[is.na(tcrs.split_ab.merged$epitope)] <- ""
  tcrs.split_ab.merged$antigen[is.na(tcrs.split_ab.merged$antigen)] <- ""
  tcrs.split_ab.merged$organism[is.na(tcrs.split_ab.merged$organism)] <- ""
  
  #tcrs.split_ab.merged$celltype <- as.character(tcrs.split_ab.merged$celltype)
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

# plot_heatmap <- function(df, subset_organism = NULL, invert = F,
#                          group_top = c("Homo sapiens", "Unknown species"),
#                          norm = F){
#   lst.biopsy <- df %>% separate_longer_delim(c(antigen, organism), delim = ";")
#   
#   if (!is.null(subset_organism)){
#     if (invert){
#       lst.biopsy <- lst.biopsy[! lst.biopsy$organism %in% subset_organism,]
#     } else {
#       lst.biopsy <- lst.biopsy[lst.biopsy$organism %in% subset_organism,]  
#     }
#   }
#   
#   tmp <- unique(lst.biopsy[,c("antigen","organism")])
#   names_dup <- names(which(table(tmp$antigen) > 1))
#   if (length(names_dup) > 0){
#     idx <- lst.biopsy$antigen %in% names_dup
#     lst.biopsy[idx,]$antigen <- paste0(lst.biopsy[idx,]$antigen," (",lst.biopsy[idx,]$organism,")")
#   }
#   annot <- unique(lst.biopsy[,c("antigen","organism")])
#   rownames(annot) <- annot$antigen
#   
#   if (norm){
#     # For each CDR3b that maps to multiple antigens, divide by the number it maps to
#     unique_samples <- unique(lst.biopsy$`subject:condition`)
#     unique_antigens <- unique(lst.biopsy$antigen)
#     x <- matrix(NA, nrow = length(unique_antigens), ncol = length(unique_samples),
#                 dimnames = list(unique_antigens, unique_samples))
#     for (i in 1:length(unique_antigens)){
#       for (j in 1:length(unique_samples)){
#         idx <- lst.biopsy$antigen == unique_antigens[i] & 
#           lst.biopsy$`subject:condition` == unique_samples[j]
#         tmp <- lst.biopsy[which(idx),]
#         if (nrow(tmp)>=1){
#           unique_CDR3b <- unique(tmp$CDR3b)
#           res <- list()
#           for (CDR3b in unique_CDR3b){
#             res[[CDR3b]] <- length(unique(lst.biopsy$antigen[lst.biopsy$CDR3b==CDR3b]))
#           }
#           res <- as.data.frame(unlist(res))
#           colnames(res) <- "n_different_antigens"
#           tmp <- merge(tmp,res,by.x="CDR3b",by.y="row.names")
#           x[i,j] <- sum(1/tmp$n_different_antigens)
#         } else {
#           x[i,j] <- nrow(tmp)
#         }
#       }
#     }
#     number_format <- "%.2f"
#   } else {
#     x <- table(lst.biopsy$antigen, lst.biopsy$`subject:condition`)
#     number_format <- "%.0f"
#   }
#   
#   annot_hum <- data.frame()
#   annot_nonhum <- data.frame()
#   if (any(annot$organism %in% group_top)){
#     annot_hum <- annot[annot$organism %in% group_top,]
#     annot_hum$organism <- factor(annot_hum$organism, level = group_top)
#     annot_hum <- annot_hum[order(annot_hum$organism),]
#   }
#   if (any(! annot$organism %in% group_top)){
#     annot_nonhum <- annot[! annot$organism %in% group_top,]
#     annot_nonhum <- annot_nonhum[order(annot_nonhum$organism),]
#   }
#   if (nrow(annot_hum) > 0 && nrow(annot_nonhum) > 0){
#     annot <- rbind(annot_hum,annot_nonhum)
#   } else if (nrow(annot_hum) > 0){
#     annot <- annot_hum
#   } else if (nrow(annot_nonhum) > 0){
#     annot <- annot_nonhum
#   }
#   
#   annot$antigen <- NULL
#   nms <- rownames(annot)[order(annot$organism)]
#   snames <- colnames(x)[order(as.numeric(gsub("[aA-zZ:]","",colnames(x))))]
#   x <- x[nms,snames]
#   x[is.na(x)] <- NA
#   x_non_na <- x
#   
#   p <- pheatmap(x, border_color = NA, display_numbers = round(x_non_na,1), 
#                 number_format = number_format, number_color = "#000000",
#                 annotation_row = annot, cluster_rows = F, cluster_cols = F,
#                 na_col = "white", colorRampPalette(c("white", "red"))(100))
#   return(p)
# }

# Write input files for TCRmatch -----------------------------------------------

# Biopsies
seu.filt <- readRDS(paste0("~/proj/um_ss/Investigations/seurat/results/v16/", "seu.filt.rda"))

# See, for trimming: https://www.frontiersin.org/articles/10.3389/fimmu.2021.640725/full
tcrs.split_ab <- get_tcrs_split_ab_non_unique(seu.filt)
tcrs.split_ab$tcrmatch_input <- ""
tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$CDR3b
tcrs.split_ab$tcrmatch_input <- gsub(pattern = "^C", replacement = "", 
                                     tcrs.split_ab$tcrmatch_input)
tcrs.split_ab$tcrmatch_input <- gsub(pattern = "[FW]$", replacement = "", 
                                     tcrs.split_ab$tcrmatch_input)
tcrs.split_ab$cell <- rownames(tcrs.split_ab)

write.table(unique(tcrs.split_ab$tcrmatch_input), 
            file = paste0(outdir,"seu.filt.tcrmatch_input.txt"),
            sep = "\t", col.names = F, quote = F, row.names = F)

# TILs
dat.til <- readRDS(file = "~/proj/um_ss/Investigations/seurat/results/v16/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
dat.til <- subset(dat.til, subset = project.name %in% c("G18-023"))

tcrs.split_ab.til <- get_tcrs_split_ab_non_unique(dat.til, "TIL")
tcrs.split_ab.til$tcrmatch_input <- ""
tcrs.split_ab.til$tcrmatch_input <- tcrs.split_ab.til$CDR3b
tcrs.split_ab.til$tcrmatch_input <- gsub(pattern = "^C", replacement = "", 
                                         tcrs.split_ab.til$tcrmatch_input)
tcrs.split_ab.til$tcrmatch_input <- gsub(pattern = "[FW]$", replacement = "", 
                                         tcrs.split_ab.til$tcrmatch_input)
tcrs.split_ab.til$cell <- rownames(tcrs.split_ab.til)

write.table(unique(tcrs.split_ab.til$tcrmatch_input), 
            file = paste0(outdir,"til.tcrmatch_input.txt"),
            sep = "\t", col.names = F, quote = F, row.names = F)

# Read TCRmatch output ---------------------------------------------------------
# outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/10x_all_TCRmatch_v4/"
# seu.filt <- readRDS(paste0("~/proj/um_ss/Investigations/seurat/results/v16/","seu.filt.rda"))
# dat.til <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/v16/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
# dat.til <- subset(dat.til, subset = project.name %in% c("G18-023"))

tcrmatch_output <- read_tcrmatch(pth = paste0(outdir,"seu.filt.tcrmatch_output_2.txt"))
tcrmatch_output.til <- read_tcrmatch(pth = paste0(outdir,"til.tcrmatch_output_2.txt"))

# TIL and Biopsy combined
seu.filt.cd8_t <- readRDS(paste0("~/proj/um_ss/Investigations/seurat/results/v16/","seu.filt.cd8_t.rda"))
tcrs.split_ab.merged <- add_metadata(tcrmatch_output, seu.filt, seu.filt.cd8_t)

tcrs.split_ab.til.merged <- add_metadata_til(tcrmatch_output.til, dat.til)

tcrs.split_ab.til.merged <- merge(tcrs.split_ab.til.merged,
                                  unique(tcrs.split_ab.merged[,c("subject","tissue_site")]),
                                  by="subject", all.x = T, all.y = F)

tcrs.split_ab.til.merged$tissue_site[tcrs.split_ab.til.merged$subject == "UM10"] <- "Liver"
tcrs.split_ab.til.merged$tissue_site[tcrs.split_ab.til.merged$subject == "UM22"] <- "Subcutaneous"

tcrs.split_ab.merged <- add_unknown_tcrs(tcrs.split_ab.merged, tcrs.split_ab)
tcrs.split_ab.til.merged <- add_unknown_tcrs(tcrs.split_ab.til.merged, tcrs.split_ab.til)

cnames <- intersect(colnames(tcrs.split_ab.til.merged), colnames(tcrs.split_ab.merged))
df <- rbind(tcrs.split_ab.til.merged[,cnames], tcrs.split_ab.merged[,cnames])

df.subcutaneous <- df[df$tissue_site == "Subcutaneous",]
df.liver <- df[df$tissue_site == "Liver",]

# Make Sankey plots ------------------------------------------------------------
df.sankey.biopsy_til.subcutaneous <- format_for_ggsankey(
  tcrs.split_ab.merged = df.subcutaneous,
  cols_include = c("cell", "subject", "condition", "antigen", "organism", "human_antigen"),
  cols_plot = c("subject", "condition", "organism", "human_antigen"))

df.sankey.biopsy_til.liver <- format_for_ggsankey(
  tcrs.split_ab.merged = df.liver, 
  cols_include = c("cell", "subject", "condition", "antigen", "organism", "human_antigen"),
  cols_plot = c("subject", "condition", "organism", "human_antigen"))

df.sankey.biopsy_til.subcutaneous <- reorder_sankey_factors(
  df.sankey.biopsy_til.subcutaneous, df)
df.sankey.biopsy_til.liver <- reorder_sankey_factors(
  df.sankey.biopsy_til.liver, df)

pdf(file = paste0(outdir,"sankey.biopsy_til.subcutaneous.pdf"), width = 10, height = 6.2)
plot_sankey(df.sankey.biopsy_til.subcutaneous, flow.alpha = .8)
dev.off()

pdf(file = paste0(outdir,"sankey.biopsy_til.liver.pdf"), width = 10, height = 6.2)
plot_sankey(df.sankey.biopsy_til.liver, flow.alpha = .8)
dev.off()