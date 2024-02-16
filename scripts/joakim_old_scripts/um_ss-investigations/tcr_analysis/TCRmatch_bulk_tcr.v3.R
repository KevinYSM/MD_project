source("~/proj/um_ss/Investigations/tcr_analysis/bin/TCRmatch_common.R")

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/tcr/TCRmatch/TCRmatch_bulk_tcr/"
dir.create(outdir, showWarnings = F, recursive = T)

# Functions for TCRmatch I/O----------------------------------------------------
convert_to_tcrs_split_format <- function(vdj.takara){
  tcrs.split_ab <- vdj.takara[vdj.takara$chains=="TRB",]
  tcrs.split_ab$tcrmatch_input <- ""
  tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$cdr3
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "^C", replacement = "", tcrs.split_ab$tcrmatch_input)
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "[FW]$", replacement = "", tcrs.split_ab$tcrmatch_input)
  colnames(tcrs.split_ab)[colnames(tcrs.split_ab) == "cdr3"] <- "CDR3b"
  colnames(tcrs.split_ab)[colnames(tcrs.split_ab) == "v_gene"] <- "TRBV"
  colnames(tcrs.split_ab)[colnames(tcrs.split_ab) == "j_gene"] <- "TRBJ"
  tcrs.split_ab[["subject:condition"]] <- paste0(tcrs.split_ab$UM.ID,":",tcrs.split_ab$vasu_phenotype)
  tcrs.split_ab$clone_id <- paste0(tcrs.split_ab$CDR3b,";",
                                   tcrs.split_ab$TRBV,";",
                                   tcrs.split_ab$d_gene,";",
                                   tcrs.split_ab$TRBJ,";",
                                   tcrs.split_ab$c_gene,";",
                                   tcrs.split_ab$`subject:condition`)
  tcrs.split_ab <- tcrs.split_ab[,c("CDR3b","TRBV","TRBJ","subject:condition","clone_id","tcrmatch_input")]
  tcrs.split_ab <- unique(tcrs.split_ab)
  rownames(tcrs.split_ab) <- tcrs.split_ab$clone_id
  tcrs.split_ab$cell <- rownames(tcrs.split_ab)
  
  return(tcrs.split_ab)
}

add_metadata <- function(tcrmatch_output, tcrs.split_ab){
  tcrs.split_ab.merged <- merge(tcrs.split_ab, tcrmatch_output, 
                                by.x="tcrmatch_input", 
                                by.y="input_sequence")#, all.x=T)
  tcrs.split_ab.merged$epitope[is.na(tcrs.split_ab.merged$epitope)] <- ""
  tcrs.split_ab.merged$antigen[is.na(tcrs.split_ab.merged$antigen)] <- ""
  tcrs.split_ab.merged$organism[is.na(tcrs.split_ab.merged$organism)] <- ""
  
  tcrs.split_ab.merged$antigen[tcrs.split_ab.merged$antigen==""] <- "Unknown protein"
  tcrs.split_ab.merged$organism[tcrs.split_ab.merged$organism==""] <- "Unknown species"
  
  tcrs.split_ab.merged$antigen <- unlist(lapply(tcrs.split_ab.merged$antigen,function(x){
    if (any(grepl("[",x,fixed=T))){
      x <- unlist(strsplit(x,";"))
      x <- paste0(unlist(lapply(x,subfun)),collapse=";")
    }
    return(x)
  }))
  
  tcrs.split_ab.merged$subject <- str_split_fixed(
    tcrs.split_ab.merged$`subject:condition`,":",2)[,1]
  tcrs.split_ab.merged$condition <- str_split_fixed(
    tcrs.split_ab.merged$`subject:condition`,":",2)[,2]
  
  tcrs.split_ab.merged <- rename_organisms(tcrs.split_ab.merged)
  tcrs.split_ab.merged <- rename_antigens(tcrs.split_ab.merged)
  
  tcrs.split_ab.merged <- tcrs.split_ab.merged[
    tcrs.split_ab.merged$antigen!="Unknown protein",]
  
  nms_dup <- names(which(table(unique(tcrs.split_ab.merged[
    ,c("antigen","organism")])$antigen) > 1))
  if (length(nms_dup)>0){
    idx <- tcrs.split_ab.merged$antigen %in% nms_dup
    tcrs.split_ab.merged[idx,]$antigen <- paste0(
      tcrs.split_ab.merged[idx,]$antigen," (", tcrs.split_ab.merged[idx,]$organism,")")
  }
  
  return(tcrs.split_ab.merged)
}

# Write input to TCRmatch ------------------------------------------------------
vdj.takara <- read_vdj.takara()
tcrs.split_ab <- convert_to_tcrs_split_format(vdj.takara)

write.table(unique(tcrs.split_ab$tcrmatch_input), 
            file = paste0(outdir,"vdj.takara.tcrmatch_input.txt"),
            sep = "\t", col.names = F, quote = F, row.names = F)

# Read TCRmatch output ---------------------------------------------------------
tcrmatch_output <- read_tcrmatch(paste0(outdir, "vdj.takara.tcrmatch_output_2.txt"))

tcrs.split_ab.merged <- add_metadata(tcrmatch_output, tcrs.split_ab)

tcrs.split_ab.merged$tissue_site <- ""
tcrs.split_ab.merged$tissue_site[tcrs.split_ab.merged$subject == "UM1"] <- "Liver"
tcrs.split_ab.merged$tissue_site[tcrs.split_ab.merged$subject == "UM9"] <- "Liver"
tcrs.split_ab.merged$tissue_site[tcrs.split_ab.merged$subject == "UM46"] <- "Subcutaneous"

tcrs.split_ab.merged <- add_unknown_tcrs(tcrs.split_ab.merged, tcrs.split_ab)

tcrs.split_ab.merged$condition <- paste0(tcrs.split_ab.merged$condition, "_condition")

df <- tcrs.split_ab.merged

df.subcutaneous <- df[df$tissue_site == "Subcutaneous",]
df.liver <- df[df$tissue_site == "Liver",]

# Make Saneky plots ------------------------------------------------------------
df.sankey.biopsy_til.subcutaneous <- format_for_ggsankey(
  tcrs.split_ab.merged = df.subcutaneous, 
  cols_include = c("cell", "subject", "antigen", "organism", "human_antigen"),
  cols_plot = c("subject", "organism", "human_antigen"))

df.sankey.biopsy_til.liver <- format_for_ggsankey(
  tcrs.split_ab.merged = df.liver, 
  cols_include = c("cell", "subject", "antigen", "organism", "human_antigen"),
  cols_plot = c("subject", "organism", "human_antigen"))

df.sankey.biopsy_til.subcutaneous <- reorder_sankey_factors(
  df.sankey.biopsy_til.subcutaneous, df.subcutaneous, condition = "MART1_condition")
df.sankey.biopsy_til.liver <- reorder_sankey_factors(
  df.sankey.biopsy_til.liver, df.liver, condition = "CD3+41bb+_condition")

pdf(file = paste0(outdir, "sankey.takara.mart1_subcutaneous.pdf"), width = 8, height = 5.2)
plot_sankey(df.sankey.biopsy_til.subcutaneous, flow.alpha = .8)
dev.off()

pdf(file = paste0(outdir, "sankey.takara.41bb_liver.pdf"), width = 8, height = 4)
plot_sankey(df.sankey.biopsy_til.liver, flow.alpha = .8)
dev.off()

# Create plots on overlap between bulk TCR-seq, biopsy and TIL -----------------
# read_vdj.takara_all <- function(){
#   fnames <- Sys.glob("~/proj/um_ss/Pipelines/takara_smarter/preprocessing/results_both/report/results_both_*_mig_cdr3_report.xlsx")
#   vdj.takara <- list()
#   for (fname in fnames){
#     sname <- str_split_fixed(basename(fname),"_",4)[,3]
#     tmp_A <- as.data.frame(read_excel(fname,sheet = paste0(sname,"_TRA_clone")),
#                            stringsAsFactors=F)
#     tmp_B <- as.data.frame(read_excel(fname,sheet = paste0(sname,"_TRB_clone")),
#                            stringsAsFactors=F)
#     
#     tmp_A$Sample <- sname
#     tmp_B$Sample <- sname
#     
#     tmp_A$Chain <- "TRA"
#     tmp_B$Chain <- "TRB"
#     
#     tmp_A$`Read Count` <- as.numeric(tmp_A$`Read Count`)
#     tmp_B$`Read Count` <- as.numeric(tmp_B$`Read Count`)
#     
#     tmp_A <- tmp_A[order(tmp_A$`Read Count`,decreasing = T),]
#     tmp_B <- tmp_B[order(tmp_B$`Read Count`,decreasing = T),]
#     
#     tmp_A$reads.rank <- 1:nrow(tmp_A)
#     tmp_B$reads.rank <- 1:nrow(tmp_B)
#     
#     vdj.takara[[basename(fname)]] <- unique(rbind(tmp_A[,c("Sample","Chain",
#                                                            "CDR3 Amino Acid Sequence",
#                                                            "V segment","D segment",
#                                                            "J segment","C segment",
#                                                            "Read Count","reads.rank")],
#                                                   tmp_B[,c("Sample","Chain",
#                                                            "CDR3 Amino Acid Sequence",
#                                                            "V segment","D segment",
#                                                            "J segment","C segment",
#                                                            "Read Count","reads.rank")]))
#   }
#   vdj.takara <- do.call("rbind",vdj.takara)
#   rownames(vdj.takara) <- NULL
#   colnames(vdj.takara)[colnames(vdj.takara)=="Read Count"] <- "reads"
#   colnames(vdj.takara)[colnames(vdj.takara)=="Chain"] <- "chains"
#   colnames(vdj.takara)[colnames(vdj.takara)=="CDR3 Amino Acid Sequence"] <- "cdr3"
#   colnames(vdj.takara)[colnames(vdj.takara)=="V segment"] <- "v_gene"
#   colnames(vdj.takara)[colnames(vdj.takara)=="D segment"] <- "d_gene"
#   colnames(vdj.takara)[colnames(vdj.takara)=="J segment"] <- "j_gene"
#   colnames(vdj.takara)[colnames(vdj.takara)=="C segment"] <- "c_gene"
#   
#   vdj.takara$v_gene[is.na(vdj.takara$v_gene)] <- "None"
#   vdj.takara$d_gene[is.na(vdj.takara$d_gene)] <- "None"
#   vdj.takara$j_gene[is.na(vdj.takara$j_gene)] <- "None"
#   vdj.takara$c_gene[is.na(vdj.takara$c_gene)] <- "None"
#   
#   vdj.takara$UM.ID <- ""
#   vdj.takara$UM.ID[vdj.takara$Sample=="A"] <- "UM1"
#   vdj.takara$UM.ID[vdj.takara$Sample=="B"] <- "UM9"
#   vdj.takara$UM.ID[vdj.takara$Sample=="C"] <- "UM46"
#   vdj.takara$UM.ID[vdj.takara$Sample=="D"] <- "UM22"
#   
#   vdj.takara$ID_2 <- ""
#   vdj.takara$ID_2[vdj.takara$Sample=="A"] <- "UM170208"
#   vdj.takara$ID_2[vdj.takara$Sample=="B"] <- "UM160411"
#   vdj.takara$ID_2[vdj.takara$Sample=="C"] <- "4418-7"
#   vdj.takara$ID_2[vdj.takara$Sample=="D"] <- "UM170322"
#   
#   vdj.takara$vasu_phenotype <- ""
#   vdj.takara$vasu_phenotype[vdj.takara$Sample=="A"] <- "CD3+41bb+"
#   vdj.takara$vasu_phenotype[vdj.takara$Sample=="B"] <- "CD3+41bb+"
#   vdj.takara$vasu_phenotype[vdj.takara$Sample=="C"] <- "MART1"
#   vdj.takara$vasu_phenotype[vdj.takara$Sample=="D"] <- "MART1"
#   
#   vdj.takara$clonotype <- apply(
#     vdj.takara[,c("chains","cdr3","v_gene","d_gene","j_gene","c_gene")], 1,
#     paste0, collapse = ":")
#   
#   vdj.takara <- vdj.takara[vdj.takara$UM.ID!="",]
#   vdj.takara <- vdj.takara[!grepl(pattern="*",vdj.takara$cdr3,fixed=T),]
#   vdj.takara <- vdj.takara[!grepl(pattern="_",vdj.takara$cdr3),]
#   
#   return(vdj.takara)
# }

#vdj.takara <- read_vdj.takara_all()
vdj.takara <- read_vdj.takara()
tcrs.split_ab <- convert_to_tcrs_split_format(vdj.takara)

dat.til <- readRDS(file = "~/proj/um_ss/Investigations/seurat/results/v16/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
dat.til <- subset(dat.til, subset = project.name %in% c("G18-023"))

tcrs.split_ab.til <- get_tcrs_split_ab_non_unique(dat.til,"TIL")
tcrs.split_ab.til$tcrmatch_input <- ""
tcrs.split_ab.til$tcrmatch_input <- tcrs.split_ab.til$CDR3b
tcrs.split_ab.til$tcrmatch_input <- gsub(pattern = "^C", replacement = "", 
                                         tcrs.split_ab.til$tcrmatch_input)
tcrs.split_ab.til$tcrmatch_input <- gsub(pattern = "[FW]$", replacement = "", 
                                         tcrs.split_ab.til$tcrmatch_input)
tcrs.split_ab.til$cell <- rownames(tcrs.split_ab.til)
tcrs.split_ab.til$subject <- str_split_fixed(tcrs.split_ab.til$`subject:condition`,":",2)[,1]
tcrs.split_ab.til$TRB_clonotype <- paste0(tcrs.split_ab.til$CDR3b,":",
                                          tcrs.split_ab.til$TRBV,":",
                                          tcrs.split_ab.til$TRBJ)

tcrs.split_ab$subject <- str_split_fixed(tcrs.split_ab$`subject:condition`,":",2)[,1]

tcrs.split_ab$TRB_clonotype <- paste0(tcrs.split_ab$CDR3b,":",
                                          tcrs.split_ab$TRBV,":",
                                          tcrs.split_ab$TRBJ)

seu.filt <- readRDS(paste0("~/proj/um_ss/Investigations/seurat/results/v16/", "seu.filt.rda"))

tcrs.split_ab.biopsy <- get_tcrs_split_ab_non_unique(seu.filt)
tcrs.split_ab.biopsy$tcrmatch_input <- ""
tcrs.split_ab.biopsy$tcrmatch_input <- tcrs.split_ab.biopsy$CDR3b
tcrs.split_ab.biopsy$tcrmatch_input <- gsub(pattern="^C", replacement="", 
                                            tcrs.split_ab.biopsy$tcrmatch_input)
tcrs.split_ab.biopsy$tcrmatch_input <- gsub(pattern="[FW]$", replacement="", 
                                            tcrs.split_ab.biopsy$tcrmatch_input)
tcrs.split_ab.biopsy$cell <- rownames(tcrs.split_ab.biopsy)
tcrs.split_ab.biopsy$subject <- str_split_fixed(tcrs.split_ab.biopsy$`subject:condition`,":",2)[,1]
tcrs.split_ab.biopsy$TRB_clonotype <- paste0(tcrs.split_ab.biopsy$CDR3b,":",
                                             tcrs.split_ab.biopsy$TRBV,":",
                                             tcrs.split_ab.biopsy$TRBJ)

tcrs.split_ab.biopsy$subject <- paste0(tcrs.split_ab.biopsy$subject,"_Biopsy")
tcrs.split_ab.til$subject <- paste0(tcrs.split_ab.til$subject,"_TILs")

tcrs.split_ab.biopsy_til <- rbind(tcrs.split_ab.biopsy, tcrs.split_ab.til)

m <-  matrix(NA,
             nrow=length(unique(tcrs.split_ab.biopsy_til$subject)),
             ncol=length(unique(tcrs.split_ab$subject)),
             dimnames=list(unique(tcrs.split_ab.biopsy_til$subject),
                           unique(tcrs.split_ab$subject)))

m <- m[order(as.numeric(gsub("[aA-zZ_]","",rownames(m))),decreasing = F),]
m <- rbind(m[grep("_Biopsy",rownames(m)),],
      m[grep("_TILs",rownames(m)),])

for (um_1 in unique(tcrs.split_ab.biopsy_til$subject)){
  for (um_2 in unique(tcrs.split_ab$subject)){
    m[um_1,um_2] <- length(intersect(
      tcrs.split_ab.biopsy_til$TRB_clonotype[tcrs.split_ab.biopsy_til$subject == um_1],
      tcrs.split_ab$TRB_clonotype[tcrs.split_ab$subject == um_2]
    ))
  }
}

# pdf(file = paste0(outdir, "TCR_overlap_Bulk_TCR_biopsy_TILs.pdf"), width = 2.5, height = 3.1)
# pheatmap(m, cluster_rows = F, cluster_cols = F,
#          display_numbers = m, border_color = NA,
#          show_rownames = T,
#          show_colnames = T, na_col = "white",
#          colorRampPalette(c("white", "red"))(100))
# dev.off()

samples_10x <- read_excel(path = "~/proj/um_ss/Investigations/samples_10x.xlsx")
samples_10x <- data.frame(samples_10x)
samples_10x <- unique(samples_10x[!is.na(samples_10x$Original_biopsy_tissue_site),])

tcrs.split_ab.biopsy_til$UM.ID <- str_split_fixed(tcrs.split_ab.biopsy_til$`subject:condition`,":",2)[,1]

tcrs.split_ab.biopsy_til <- merge(tcrs.split_ab.biopsy_til, samples_10x, by="UM.ID",all.x=T,all.y=F)

annot <- unique(data.frame(
  subject = tcrs.split_ab.biopsy_til$subject,
  Original_biopsy_tissue_site = tcrs.split_ab.biopsy_til$Original_biopsy_tissue_site))

annot[annot$subject=="UM10_TILs",]$Original_biopsy_tissue_site <- "Liver"

rownames(annot) <- annot$subject
annot$subject <- NULL

pdf(file = paste0(outdir, "TCR_overlap_Bulk_TCR_biopsy_TILs.no_UM22.pdf"), width = 4.5, height = 4.5)
pheatmap(m, cluster_rows = F, cluster_cols = F, 
         display_numbers = m, border_color = NA,
         show_rownames = T,
         show_colnames = T, na_col = "white",
         colorRampPalette(c("white", "red"))(100),
         annotation_row = annot)
dev.off()
