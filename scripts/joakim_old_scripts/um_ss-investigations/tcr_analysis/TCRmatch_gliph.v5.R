source("~/proj/um_ss/Investigations/tcr_analysis/bin/TCRmatch_common.R")

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/tcr/TCRmatch/TCRmatch_gliph/"
dir.create(outdir, showWarnings = F, recursive = T)

# Functions for TCRmatch I/O----------------------------------------------------
read_gliph <- function(){
  gliph <- read.table(paste0("~/proj/um_ss/Investigations/seurat/results/v16/",
                             "seu.filt.cd8_t_hla_ref_v2.csv"), sep=",", header = T)
  gliph <- gliph[!is.na(gliph$index),]
  gliph$original_gliph_row <- rownames(gliph)
  return(gliph)
}

convert_to_tcrs_split_format <- function(gliph){
  tcrs.split_ab <- gliph
  tcrs.split_ab$tcrmatch_input <- ""
  tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$TcRb
  tcrs.split_ab$tcrmatch_input <- gsub(pattern="^C", replacement="", 
                                       tcrs.split_ab$tcrmatch_input)
  tcrs.split_ab$tcrmatch_input <- gsub(pattern="[FW]$", replacement="", 
                                       tcrs.split_ab$tcrmatch_input)
  colnames(tcrs.split_ab)[colnames(tcrs.split_ab) == "TcRb"] <- "CDR3b"
  colnames(tcrs.split_ab)[colnames(tcrs.split_ab) == "V"] <- "TRBV"
  colnames(tcrs.split_ab)[colnames(tcrs.split_ab) == "J"] <- "TRBJ"
  colnames(tcrs.split_ab)[colnames(tcrs.split_ab) == "Sample"] <- "subject:condition"
  
  tcrs.split_ab$cluster <- paste0("C", tcrs.split_ab$index)
  
  tcrs.split_ab$clone_id <- paste0(tcrs.split_ab$CDR3b,";",
                                   tcrs.split_ab$TRBV,";",
                                   tcrs.split_ab$TRBJ,";",
                                   tcrs.split_ab$`subject:condition`)
  tcrs.split_ab <- tcrs.split_ab[,c("CDR3b","TRBV","TRBJ","subject:condition",
                                    "clone_id","tcrmatch_input",
                                    "original_gliph_row","cluster")]
  rownames(tcrs.split_ab) <- tcrs.split_ab$original_gliph_row
  
  tcrs.split_ab$cell <- rownames(tcrs.split_ab)
  
  return(tcrs.split_ab)
}

add_metadata <- function(tcrmatch_output, tcrs.split_ab, gliph, seu.filt.cd8_t){
  gliph.merged <- merge(tcrs.split_ab, gliph, by = "original_gliph_row", all=T)
  stopifnot(nrow(gliph.merged) == nrow(tcrs.split_ab))
  tcrs.split_ab <- gliph.merged
  
  tcrs.split_ab.merged <- merge(tcrs.split_ab, tcrmatch_output, 
                                by.x="tcrmatch_input", 
                                by.y="input_sequence")#, all.x=T)
  tcrs.split_ab.merged$epitope[is.na(tcrs.split_ab.merged$epitope)] <- ""
  tcrs.split_ab.merged$antigen[is.na(tcrs.split_ab.merged$antigen)] <- ""
  tcrs.split_ab.merged$organism[is.na(tcrs.split_ab.merged$organism)] <- ""
  
  tcrs.split_ab.merged$antigen[tcrs.split_ab.merged$antigen==""] <- "Unknown protein"
  tcrs.split_ab.merged$organism[tcrs.split_ab.merged$organism==""] <- "Unknown species"
  
  tcrs.split_ab.merged$antigen <- unlist(lapply(tcrs.split_ab.merged$antigen, function(x){
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
  
  seu.filt.cd8_t <- add_anatomic_location_meta(
    seu.filt.cd8_t, pth = "~/proj/um_ss/Investigations/samples_10x.xlsx")
  
  tcrs.split_ab.merged <- merge(
    tcrs.split_ab.merged, 
    unique(data.frame(UM.ID = seu.filt.cd8_t@meta.data$UM.ID.x,
                      tissue_site = seu.filt.cd8_t@meta.data$Original_biopsy_tissue_site)), 
    by.x = "subject", by.y = "UM.ID",all.x=T)
  
  tcrs.split_ab.merged <- rename_organisms(tcrs.split_ab.merged)
  tcrs.split_ab.merged <- rename_antigens(tcrs.split_ab.merged)
  
  tcrs.split_ab.merged <- tcrs.split_ab.merged[
    tcrs.split_ab.merged$antigen!="Unknown protein",]
  
  nms_dup <- names(which(table(unique(
    tcrs.split_ab.merged[,c("antigen","organism")])$antigen) > 1))
  if (length(nms_dup)>0){
    idx <- tcrs.split_ab.merged$antigen %in% nms_dup
    tcrs.split_ab.merged[idx,]$antigen <- paste0(
      tcrs.split_ab.merged[idx,]$antigen," (", tcrs.split_ab.merged[idx,]$organism,")")
  }
  
  return(tcrs.split_ab.merged)
}

write_gliph_tcrmatch_results_xlsx <- function(gliph,tcrmatch_output, outdir){
  gliph$tcrmatch_input <- ""
  gliph$tcrmatch_input <- gliph$TcRb
  gliph$tcrmatch_input <- gsub(pattern="^C", replacement="", gliph$tcrmatch_input)
  gliph$tcrmatch_input <- gsub(pattern="[FW]$", replacement="", gliph$tcrmatch_input)
  gliph.merged <- merge(gliph,tcrmatch_output, by.x = "tcrmatch_input",
                        by.y = "input_sequence", all.x = T)
  gliph.merged <- gliph.merged[order(gliph.merged$index, decreasing = F),]
  gliph.merged$epitope[is.na(gliph.merged$epitope)] <- ""
  gliph.merged$antigen[is.na(gliph.merged$antigen)] <- ""
  gliph.merged$organism[is.na(gliph.merged$organism)] <- ""
  
  WriteXLS(gliph.merged, ExcelFileName = paste0(outdir,"gliph.merged.xlsx"),
           BoldHeaderRow = T, AutoFilter = T, AdjWidth = T, row.names = F)
}

get_tcrs_split_ab_non_unique_gliph <- function(seu.filt.cd8_t){
  tcrs <- seu.filt.cd8_t@meta.data[which(seu.filt.cd8_t@meta.data$n_chains==2),
                                   c("cdr3","v_gene","j_gene","chains","UM.ID")]
  tcrs$`subject:condition` <- paste0(tcrs$UM.ID,":","biopsy")
  tcrs$UM.ID.x <- NULL
  
  stopifnot(identical(unlist(lapply(strsplit(tcrs$chain,";"),function(x){length(x)})),
                      unlist(lapply(strsplit(tcrs$cdr3,";"),function(x){length(x)}))))
  stopifnot(identical(unlist(lapply(strsplit(tcrs$chain,";"),function(x){length(x)})),
                      unlist(lapply(strsplit(tcrs$v_gene,";"),function(x){length(x)}))))
  stopifnot(identical(unlist(lapply(strsplit(tcrs$chain,";"),function(x){length(x)})),
                      unlist(lapply(strsplit(tcrs$j_gene,";"),function(x){length(x)}))))
  
  tcrs.split_ab <- list()
  for (i in 1:nrow(tcrs)){
    chains <- unlist(strsplit(tcrs[i,]$chains,";"))
    cdr3 <- unlist(strsplit(tcrs[i,]$cdr3,";"))
    v_gene <- unlist(strsplit(tcrs[i,]$v_gene,";"))
    j_gene <- unlist(strsplit(tcrs[i,]$j_gene,";"))
    cdr3a <- cdr3[chains=="TRA"]
    v_gene_a <- v_gene[chains=="TRA"]
    j_gene_a <- j_gene[chains=="TRA"]
    cdr3b <- cdr3[chains=="TRB"]
    v_gene_b <- v_gene[chains=="TRB"]
    j_gene_b <- j_gene[chains=="TRB"]
    cdr3a <- ifelse(length(cdr3a)>1,paste0(cdr3a,collapse=";"),cdr3a)
    v_gene_a <- ifelse(length(v_gene_a)>1,paste0(v_gene_a,collapse=";"),v_gene_a)
    j_gene_a <- ifelse(length(j_gene_a)>1,paste0(j_gene_a,collapse=";"),j_gene_a)
    cdr3b <- ifelse(length(cdr3b)>1,paste0(cdr3b,collapse=";"),cdr3b)
    v_gene_b <- ifelse(length(v_gene_b)>1,paste0(v_gene_b,collapse=";"),v_gene_b)
    j_gene_b <- ifelse(length(j_gene_b)>1,paste0(j_gene_b,collapse=";"),j_gene_b)
    
    tcrs.split_ab[[i]] <- data.frame(cdr3a,v_gene_a,j_gene_a,cdr3b,v_gene_b,j_gene_b,
                                     tcrs[i,]$`subject:condition`)
  }
  
  tcrs.split_ab <- do.call("rbind",tcrs.split_ab)
  
  tcrs.split_ab <- tcrs.split_ab[,c("cdr3b","v_gene_b","j_gene_b","cdr3a",
                                    "tcrs.i.....subject.condition.")]
  
  colnames(tcrs.split_ab) <- c("CDR3b", "TRBV", "TRBJ", "CDR3a", "subject:condition")
  
  tcrs.split_ab$clone_id <- apply(tcrs.split_ab, 1, function(x){
    paste0(x, collapse=";")
  })
  
  rownames(tcrs.split_ab) <- rownames(tcrs)
  return(tcrs.split_ab)
}

# Write input files for TCRmatch -----------------------------------------------
gliph <- read_gliph()
tcrs.split_ab <- convert_to_tcrs_split_format(gliph)

write.table(unique(tcrs.split_ab$tcrmatch_input),
            file = paste0(outdir,"gliph.cd8_t_hla_ref_v2.tcrmatch_input.txt"),
            sep = "\t",col.names = F,quote = F,row.names = F)

# Read TCRmatch output ---------------------------------------------------------
tcrmatch_output <- read_tcrmatch(paste0(outdir, "gliph.cd8_t_hla_ref_v2.tcrmatch_output_2.txt"))

write_gliph_tcrmatch_results_xlsx(gliph, tcrmatch_output, outdir)

seu.filt.cd8_t <- readRDS(paste0("~/proj/um_ss/Investigations/seurat/results/v16/",
                                 "seu.filt.cd8_t.rda"))

tcrs.split_ab.merged <- add_metadata(tcrmatch_output, tcrs.split_ab, gliph, seu.filt.cd8_t)
tcrs.split_ab.merged <- add_unknown_tcrs(tcrs.split_ab.merged, tcrs.split_ab)

df <- tcrs.split_ab.merged

high_confidence <- gliph$vb_score < 0.05 & 
  gliph$index %in% names(which(table(gliph$index)>2)) & 
  gliph$index %in% names(which(rowSums(table(gliph$index, gliph$Sample)>0)>2))
gliph.high_confidence <- gliph[high_confidence,]
nodes_high_confidence <- paste0("C",unique(gliph.high_confidence$index))

df.high_confidence <- df[df$cluster %in% nodes_high_confidence,]

df.subcutaneous <- df[df$tissue_site == "Subcutaneous",]
df.liver <- df[df$tissue_site == "Liver",]

df.high_confidence.subcutaneous <- df.high_confidence[
  df.high_confidence$tissue_site == "Subcutaneous",]
df.high_confidence.liver <- df.high_confidence[
  df.high_confidence$tissue_site == "Liver",]

# Make Sankey plots ------------------------------------------------------------
df.sankey.high_confidence <- format_for_ggsankey(
  tcrs.split_ab.merged = df.high_confidence,
  cols_include = c("cluster", "subject", "antigen", "organism", "human_antigen"),
  cols_plot = c("cluster", "organism", "human_antigen"))

df.sankey.all <- format_for_ggsankey(
  tcrs.split_ab.merged = df,
  cols_include = c("cell", "subject", "antigen", "organism", "human_antigen"),
  cols_plot = c("subject", "organism", "human_antigen"))

df.sankey.all <- reorder_sankey_factors(df.sankey.all, df)
df.sankey.high_confidence <- reorder_sankey_factors_cluster(
  df.sankey.high_confidence, df)

pdf(file = paste0(outdir, "sankey.high_confidence.pdf"), width = 6, height = 3)
plot_sankey(df.sankey.high_confidence, flow.alpha = .8)
dev.off()

pdf(file = paste0(outdir, "sankey.all.pdf"), width = 10, height = 8.5)
plot_sankey(df.sankey.all, flow.alpha = .8)
dev.off()

# Map gliph clusters and antigen matches to tSNE -------------------------------
gliph.high_confidence$clonotype_id <- paste0(gliph.high_confidence$TcRa,":",
                                             gliph.high_confidence$TcRb,":",
                                             gliph.high_confidence$V,":",
                                             gliph.high_confidence$J)

tcrs.split_ab <- get_tcrs_split_ab_non_unique_gliph(seu.filt.cd8_t)
tcrs.split_ab$clonotype_id <- paste0(tcrs.split_ab$CDR3a,":",tcrs.split_ab$CDR3b,
                                     ":",tcrs.split_ab$TRBV,":",tcrs.split_ab$TRBJ)
gliph.high_confidence$cluster <- paste0("C",gliph.high_confidence$index)

# Plot by cluster
unique_clusters <- unique(gliph.high_confidence$cluster)
clusters.cells <- list()
for (clus in unique_clusters){
  clusters.cells[[clus]] <- rownames(tcrs.split_ab)[
    tcrs.split_ab$clonotype_id %in% gliph.high_confidence$clonotype_id[
      gliph.high_confidence$cluster == clus]]
}

# Plot by clusters with information about species matches
links_cluster_species <- df[df$cluster %in% gliph.high_confidence$cluster,
                            c("cluster","organism")]
colnames(links_cluster_species) <- c("Cluster","Species")
x <- table(links_cluster_species$Cluster, links_cluster_species$Species)
x <- as.data.frame.matrix(x)
species_nrs <- list()
for (clus in unique_clusters){
  nms <- colnames(x[clus,])[which(x[clus,]>0)]
  nms <- sort(unique(unlist(strsplit(nms, ";"))))
  nms <- paste0(nms, collapse = ";")
  species_nrs[[clus]] <- sort(nms)
}

p <- list()
for (clus in unique_clusters){
  p[[clus]] <- DimPlot(seu.filt.cd8_t, cells.highlight = clusters.cells[[clus]]) + 
    theme(legend.position="none") + ggtitle(paste0(clus,":",species_nrs[[clus]])) + 
    theme(plot.title = element_text(size = 8))
}
pdf(file = paste0(outdir, "tsne_gliph_by_cluster.species.pdf"), width = 20, height = 10)
grid.arrange(grobs = p, ncol = 4)
dev.off()

# List co-enriched alleles -----------------------------------------------------
unique_clusters <- paste0("C", unique(gliph$index))
co_enriched_hla <- list()
for (clus in unique_clusters){
  tmp <- gliph[paste0("C",gliph$index) == clus,]
  x.a <- unlist(strsplit(tmp$HLA.A,"/"))
  x.a <- unique(x.a[grep(pattern="!",x.a)])
  x.b <- unlist(strsplit(tmp$HLA.B,"/"))
  x.b <- unique(x.b[grep(pattern="!",x.b)])
  x.c <- unlist(strsplit(tmp$HLA.C,"/"))
  x.c <- unique(x.c[grep(pattern="!",x.c)])
  
  x.a <- ifelse(length(x.a)>0,x.a,"")
  x.b <- ifelse(length(x.b)>0,x.b,"")
  x.c <- ifelse(length(x.c)>0,x.c,"")
  
  co_enriched_hla[[clus]] <- data.frame(HLA.A = x.a, HLA.B = x.b, HLA.C = x.c)
}

co_enriched_hla <- do.call("rbind", co_enriched_hla)
co_enriched_hla$is_high_confidence <- rownames(co_enriched_hla) %in% 
  gliph.high_confidence$cluster

WriteXLS(co_enriched_hla, ExcelFileName = paste0(outdir, "co_enriched_hla.xlsx"),
         AutoFilter = T, BoldHeaderRow = T, row.names = T, col.names = T,
         AdjWidth = T)

# Plot all MART1 matching CDR3 -------------------------------------------------
tcrs.split_ab <- get_tcrs_split_ab_non_unique_gliph(seu.filt.cd8_t)

lst <- df %>% separate_longer_delim(c(antigen, organism), delim = ";")
lst$clone_id <-  paste0(lst$TcRb,";",
                        lst$V,";",
                        lst$J,";",
                        lst$TcRa,";",
                        lst$`subject:condition`)

cells <- rownames(tcrs.split_ab)[tcrs.split_ab$clone_id %in% 
                                   lst[lst$antigen == "MART1",]$clone_id]

# all_cells <- rownames(seu.filt.cd8_t@meta.data)
# 
# mart1_matches_gliph <- setNames(rep("None",length(all_cells)),all_cells)
# mart1_matches_gliph[cells] <- "MART1_high_conf_gliph_cluster"

# Individual cells with MART1 matches in gliph clusters
pdf(file = paste0(outdir,"tsne_gliph_all_mart1_matches_cells.pdf"), 
    width = 10, height = 10)
DimPlot( seu.filt.cd8_t, cells.highlight = cells) + 
  theme(legend.position="none") + ggtitle("MART1 matches") + 
  theme(plot.title = element_text(size = 8))
dev.off()

gliph$clonotype_id <- paste0(gliph$TcRb,";",
                             gliph$V,";",
                             gliph$J,";",
                             gliph$TcRa,";",
                             gliph$Sample)
gliph$cluster <- paste0("C",gliph$index)

unique_clusters <- unique(gliph$cluster)
clusters.cells <- list()
for (clus in unique_clusters){
  clusters.cells[[clus]] <- rownames(tcrs.split_ab)[
    tcrs.split_ab$clone_id %in% gliph$clonotype_id[
      gliph$cluster==clus]]
}

clusters.cells <- clusters.cells[which(names(clusters.cells) %in% 
                                         lst[lst$antigen == "MART1",]$cluster)]
clusters.cells <- unique(unlist(clusters.cells))

# mart1_matches_gliph[setdiff(clusters.cells,cells)] <- "MART1_low_conf_gliph_cluster"
# seu.filt.cd8_t <- AddMetaData(seu.filt.cd8_t, metadata = mart1_matches_gliph, col.name = "mart1_matches_gliph")

# All cells in clusters with any member having MART1 match
pdf(file = paste0(outdir,"tsne_gliph_all_mart1_matches_clusters.pdf"), 
    width = 10, height = 10)
DimPlot( seu.filt.cd8_t, cells.highlight = clusters.cells) + 
  theme(legend.position="none") + ggtitle("MART1 matches") + 
  theme(plot.title = element_text(size = 8))
dev.off()

# seu.filt.cd8_t@meta.data$mart1_matches_gliph <- factor(
#   seu.filt.cd8_t@meta.data$mart1_matches_gliph, 
#   levels = c("MART1_high_conf_gliph_cluster","MART1_low_conf_gliph_cluster","None"))

pdf(file = paste0(outdir,"tsne_gliph_all_mart1_matches_clusters_v2.pdf"), 
    width = 13, height = 10)
DimPlot( seu.filt.cd8_t, group.by = "mart1_matches_gliph", 
         cells.highlight = list(MART1_high_conf_gliph_cluster = cells,
                                MART1_low_conf_gliph_cluster = setdiff(clusters.cells,cells)),
         cols.highlight =  c("blue","red"),pt.size = 1) + 
  ggtitle("MART1 matches") + 
  theme(plot.title = element_text(size = 8))
dev.off()