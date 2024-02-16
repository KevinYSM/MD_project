library(readxl)
library(WriteXLS)
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(reshape2)
library(ggsankey)
library(dplyr)
library(Seurat)
library(gridExtra)
library(pheatmap)

# Functions for TCRmatch I/O and annotation ------------------------------------
read_tcrmatch <- function(pth){
  tcrmatch_output <- read.table(pth, sep = "\t", header = T, fill = T)
  stopifnot(all(is.na(tcrmatch_output[,8])))
  tcrmatch_output[,8] <- NULL
  tcrmatch_output$original_tcrmatch_row <- rownames(tcrmatch_output)
  return(tcrmatch_output)
}

rename_organisms <- function(tcrs.split_ab.merged){
  organism_rename_map <- get_organism_rename_map(tcrs.split_ab.merged)
  for (i in 1:nrow(tcrs.split_ab.merged)){
    y <- unlist(strsplit(tcrs.split_ab.merged[i,]$organism,";"))
    for (j in 1:length(y)){
      if (y[j] != "Unknown species"){
        y[j] <- organism_rename_map$new[organism_rename_map$old == y[j]]
      }
    }
    tcrs.split_ab.merged[i,]$organism <- paste0(y,collapse=";")
  }
  return(tcrs.split_ab.merged)
}

rename_antigens <- function(tcrs.split_ab.merged){
  antigen_rename_map <- get_antigen_rename_map(tcrs.split_ab.merged)
  for (i in 1:nrow(tcrs.split_ab.merged)){
    y <- unlist(strsplit(tcrs.split_ab.merged[i,]$antigen,";"))
    for (j in 1:length(y)){
      if (y[j] != "Unknown species"){
        y[j] <- antigen_rename_map$new[antigen_rename_map$old == y[j]]
      }
    }
    tcrs.split_ab.merged[i,]$antigen <- paste0(y,collapse=";")
  }
  return(tcrs.split_ab.merged)
}

get_organism_rename_map <- function(tcrmatch_output){
  y <- unique(unlist(strsplit(tcrmatch_output$organism,";")))
  #x <- sort(str_split_fixed(x," \\(",2)[,1])
  x <- str_split_fixed(y," \\(",2)[,1]
  x[grep("Human herpesvirus 5",x)] <- "CMV"
  x[grep("Human herpesvirus 4",x)] <- "EBV"
  x[grep("Human immunodeficiency virus 1",x)] <- "HIV-1"
  x[x=="Human T-cell leukemia virus type I"] <- "HTLV-1"
  x[x=="Influenza A virus"] <- "Influenza A"
  x[x=="Influenza A virus H3N2"] <- "Influenza A"
  x[x=="Yellow fever virus 17D"] <- "YFV"
  x[x=="SARS-CoV1"] <- "Human coronavirus"
  x[x=="SARS-CoV2"] <- "Human coronavirus"
  x[x=="SARS coronavirus BJ01"] <- "Human coronavirus"
  x[x=="SARS coronavirus Urbani"] <- "Human coronavirus"
  x[x=="Mycobacterium tuberculosis"] <- "Tuberculosis"
  x[x=="Mycobacterium tuberculosis H37Rv"] <- "Tuberculosis"
  x[x=="Human coronavirus 229E"] <- "Human coronavirus"
  x[x=="Human coronavirus NL63"] <- "Human coronavirus"
  x[x=="Human coronavirus HKU1"] <- "Human coronavirus"
  x[x=="Lymphocytic choriomeningitis mammarenavirus"] <- "LCMV"
  x[x=="Dengue virus 2 Thailand/16681/84"] <- "Dengue virus 2"
  x[x=="Dengue virus type 2"] <- "Dengue virus 2"
  x[x=="Dengue virus type 1"] <- "Dengue virus 1"
  x[x=="dengue virus type I"] <- "Dengue virus 1"
  x[x=="Dengue virus type 3"] <- "Dengue virus 3"
  x[x=="Hepatitis C virus"] <- "Hepatitis C"
  x[x=="hepatitis C virus genotype 1a"] <- "Hepatitis C"
  x[x=="Hepatitis B virus"] <- "Hepatitis B"
  x[x=="Human rotavirus strain WA"] <- "Human rotavirus"
  x[x=="Plasmodium falciparum NF54"] <- "Plasmodium falciparum"
  x[x=="Cytomegalovirus"] <- "CMV"
  x[x=="Epstein Barr virus"] <- "EBV"
  x[x=="Human T-lymphotropic virus 1"] <- "HTLV-1"
  x[x=="HomoSapiens"] <- "Homo sapiens"
  x[x=="Influenza"] <- "Influenza A"
  x[x=="Human immunodeficiency virus"] <- "HIV-1"
  x[x=="Murid herpesvirus 1 deltaMS94.5"] <- "Murid betaherpesvirus 1"
  
  orgamism_rename_map <- data.frame(new=x,old=y)
  
  return(orgamism_rename_map)
}

get_antigen_rename_map <- function(tcrmatch_output){
  y <- unique(unlist(strsplit(tcrmatch_output$antigen,";")))
  #x <- str_split_fixed(y," \\(",2)[,1]
  x <- y
  
  x[x=="Melanoma antigen recognized by T-cells 1"] <- "MART1"
  x[x=="Melanoma-associated antigen 1"] <- "MAGEA1"
  x[x=="G1/S-specific cyclin-D1"] <- "CCND1"
  x[x=="Coagulation factor VIII"] <- "F8"
  x[x=="phosphorylase b kinase regulatory subunit alpha_ liver isoform"] <- "PHKA2"
  x[x=="Transcriptional enhancer factor TEF-1"] <- "TEAD1"
  x[x=="Tribbles homolog 2"] <- "TRIB2"
  x[x=="islet-specific glucose-6-phosphatase-related protein isoform 1"] <- "G6PC2"
  x[x=="Myelin basic protein"] <- "MBP"
  x[x=="orexin precursor"] <- "HCRT"
  x[x=="Palmitoyltransferase ZDHHC7"] <- "ZDHHC7"
  x[x=="Insulin precursor"] <- "INS"
  x[x=="zinc transporter 8 isoform a"] <- "SLC30A8"
  x[x=="HAUS augmin-like complex subunit 3"] <- "HAUS3"
  x[x=="histocompatibility (minor) HA-1"] <- "ARHGAP45"
  x[x=="Histone acetyltransferase KAT6A"] <- "KAT6A"
  x[x=="Sterol-4-alpha-carboxylate 3-dehydrogenase_ decarboxylating"] <- "NSDHL"
  x[x=="HLA class I histocompatibility antigen_ A-2 alpha chain precursor"] <- "HLA-A"
  x[x=="Melanocyte protein PMEL"] <- "PMEL"
  x[x=="Sarcospan"] <- "SSPN"
  x[x=="Sterol regulatory element-binding protein 1"] <- "SREBF1"
  
  antigen_rename_map <- data.frame(new=x,old=y)
  
  return(antigen_rename_map)
}

subfun <- function(x){
  if (any(grepl("[",x,fixed=T))){
    # Check that all antigens with [ somewhere in the name all end with ]
    stopifnot(substr(x, nchar(x), nchar(x)) == "]")
    x <- str_split_fixed(x," *\\[",2)[,1]
  }
  return(x)
}

add_anatomic_location_meta <- function(seu.filt.cd8_t, pth){
  samples_10x <- read_excel(path=pth)
  samples_10x <- data.frame(samples_10x)
  
  seu.filt.cd8_t@meta.data$Sample.ID <- seu.filt.cd8_t@meta.data$orig.ident
  meta <- seu.filt.cd8_t@meta.data
  meta$Rownames <- rownames(meta)
  
  meta <- merge(meta, samples_10x, by="Sample.ID",all.x=T,all.y=F)
  rownames(meta) <- meta$Rownames
  meta$Rownames <- NULL
  meta <- meta[rownames(seu.filt.cd8_t@meta.data),]
  stopifnot(identical(rownames(meta),rownames(seu.filt.cd8_t@meta.data)))
  seu.filt.cd8_t@meta.data <- meta
  
  return(seu.filt.cd8_t)
}

add_unknown_tcrs <- function(tcrs.split_ab.merged, tcrs.split_ab){
  unknown_cells <- setdiff(tcrs.split_ab$cell, tcrs.split_ab.merged$cell)
  tcrs.split_ab <- tcrs.split_ab[tcrs.split_ab$cell %in% unknown_cells,]
  tcrs.split_ab$antigen <- "Unknown protein"
  tcrs.split_ab$organism <- "Unknown organism"
  rows_add <- tcrs.split_ab
  rows_add$subject <- str_split_fixed(rows_add$`subject:condition`,":",2)[,1]
  rows_add$condition <- str_split_fixed(rows_add$`subject:condition`,":",2)[,2]
  
  rows_add <- merge(rows_add, unique(tcrs.split_ab.merged[,c("subject","tissue_site")]), 
                    by = "subject", all.x=T, all.y=F)
  
  cnames_missing <- setdiff(colnames(tcrs.split_ab.merged),colnames(rows_add))
  for (cname in cnames_missing){
    rows_add[[cname]] <- NA
  }
  rows_add <- rows_add[,colnames(tcrs.split_ab.merged)]
  tcrs.split_ab.merged <- rbind(tcrs.split_ab.merged,rows_add)
  
  return(tcrs.split_ab.merged)
}

get_tcrs_split_ab_non_unique <- function(seu.filt.cd8_t, condition = "biopsy"){
  tcrs <- seu.filt.cd8_t@meta.data[which(unlist(lapply(
    seu.filt.cd8_t@meta.data$chains,function(x){length(
      which(unlist(strsplit(x,";"))=="TRB"))})) == 1),
    c("cdr3","v_gene","j_gene","chains","UM.ID")]
  tcrs$`subject:condition` <- paste0(tcrs$UM.ID,":",condition)
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
    cdr3b <- cdr3[chains=="TRB"]
    v_gene_b <- v_gene[chains=="TRB"]
    j_gene_b <- j_gene[chains=="TRB"]
    cdr3b <- ifelse(length(cdr3b)>1,paste0(cdr3b,collapse=";"),cdr3b)
    v_gene_b <- ifelse(length(v_gene_b)>1,paste0(v_gene_b,collapse=";"),v_gene_b)
    j_gene_b <- ifelse(length(j_gene_b)>1,paste0(j_gene_b,collapse=";"),j_gene_b)
    
    tcrs.split_ab[[i]] <- data.frame(cdr3b,v_gene_b,j_gene_b,
                                     tcrs[i,]$`subject:condition`)
  }
  
  tcrs.split_ab <- do.call("rbind",tcrs.split_ab)
  
  tcrs.split_ab <- tcrs.split_ab[,c("cdr3b","v_gene_b","j_gene_b",
                                    "tcrs.i.....subject.condition.")]
  
  colnames(tcrs.split_ab) <- c("CDR3b","TRBV","TRBJ","subject:condition")
  
  tcrs.split_ab$clone_id <- apply(tcrs.split_ab, 1, function(x){
    paste0(x, collapse=";")
  })
  
  rownames(tcrs.split_ab) <- rownames(tcrs)
  return(tcrs.split_ab)
}

# Read VDJ - Takara ------------------------------------------------------------
read_vdj.takara <- function(){
  fnames <- Sys.glob("~/proj/um_ss/Pipelines/takara_smarter/preprocessing/results_both/report/results_both_*_mig_cdr3_report.xlsx")
  vdj.takara <- list()
  for (fname in fnames){
    sname <- str_split_fixed(basename(fname),"_",4)[,3]
    tmp_A <- as.data.frame(read_excel(fname,sheet = paste0(sname,"_TRA_clone")),
                           stringsAsFactors=F)
    tmp_B <- as.data.frame(read_excel(fname,sheet = paste0(sname,"_TRB_clone")),
                           stringsAsFactors=F)
    
    tmp_A$Sample <- sname
    tmp_B$Sample <- sname
    
    tmp_A$Chain <- "TRA"
    tmp_B$Chain <- "TRB"
    
    tmp_A$`Read Count` <- as.numeric(tmp_A$`Read Count`)
    tmp_B$`Read Count` <- as.numeric(tmp_B$`Read Count`)
    
    tmp_A <- tmp_A[order(tmp_A$`Read Count`,decreasing = T),]
    tmp_B <- tmp_B[order(tmp_B$`Read Count`,decreasing = T),]
    
    tmp_A$reads.rank <- 1:nrow(tmp_A)
    tmp_B$reads.rank <- 1:nrow(tmp_B)
    
    vdj.takara[[basename(fname)]] <- unique(rbind(tmp_A[,c("Sample","Chain",
                                                           "CDR3 Amino Acid Sequence",
                                                           "V segment","D segment",
                                                           "J segment","C segment",
                                                           "Read Count","reads.rank")],
                                                  tmp_B[,c("Sample","Chain",
                                                           "CDR3 Amino Acid Sequence",
                                                           "V segment","D segment",
                                                           "J segment","C segment",
                                                           "Read Count","reads.rank")]))
  }
  vdj.takara <- do.call("rbind",vdj.takara)
  rownames(vdj.takara) <- NULL
  colnames(vdj.takara)[colnames(vdj.takara)=="Read Count"] <- "reads"
  colnames(vdj.takara)[colnames(vdj.takara)=="Chain"] <- "chains"
  colnames(vdj.takara)[colnames(vdj.takara)=="CDR3 Amino Acid Sequence"] <- "cdr3"
  colnames(vdj.takara)[colnames(vdj.takara)=="V segment"] <- "v_gene"
  colnames(vdj.takara)[colnames(vdj.takara)=="D segment"] <- "d_gene"
  colnames(vdj.takara)[colnames(vdj.takara)=="J segment"] <- "j_gene"
  colnames(vdj.takara)[colnames(vdj.takara)=="C segment"] <- "c_gene"
  
  vdj.takara$v_gene[is.na(vdj.takara$v_gene)] <- "None"
  vdj.takara$d_gene[is.na(vdj.takara$d_gene)] <- "None"
  vdj.takara$j_gene[is.na(vdj.takara$j_gene)] <- "None"
  vdj.takara$c_gene[is.na(vdj.takara$c_gene)] <- "None"
  
  vdj.takara$UM.ID <- ""
  vdj.takara$UM.ID[vdj.takara$Sample=="A"] <- "UM1"
  vdj.takara$UM.ID[vdj.takara$Sample=="B"] <- "UM9"
  vdj.takara$UM.ID[vdj.takara$Sample=="C"] <- "UM46"
  vdj.takara$UM.ID[vdj.takara$Sample=="D"] <- "UM22"
  
  vdj.takara$ID_2 <- ""
  vdj.takara$ID_2[vdj.takara$Sample=="A"] <- "UM170208"
  vdj.takara$ID_2[vdj.takara$Sample=="B"] <- "UM160411"
  vdj.takara$ID_2[vdj.takara$Sample=="C"] <- "4418-7"
  vdj.takara$ID_2[vdj.takara$Sample=="D"] <- "UM170322"
  
  vdj.takara$vasu_phenotype <- ""
  vdj.takara$vasu_phenotype[vdj.takara$Sample=="A"] <- "CD3+41bb+"
  vdj.takara$vasu_phenotype[vdj.takara$Sample=="B"] <- "CD3+41bb+"
  vdj.takara$vasu_phenotype[vdj.takara$Sample=="C"] <- "MART1"
  vdj.takara$vasu_phenotype[vdj.takara$Sample=="D"] <- "MART1"
  
  vdj.takara$clonotype <- apply(
    vdj.takara[,c("chains","cdr3","v_gene","d_gene","j_gene","c_gene")], 1,
    paste0, collapse = ":")
  
  vdj.takara <- vdj.takara[vdj.takara$UM.ID!="",]
  vdj.takara <- vdj.takara[!grepl(pattern="*",vdj.takara$cdr3,fixed=T),]
  vdj.takara <- vdj.takara[!grepl(pattern="_",vdj.takara$cdr3),]
  
  # Remove UM22
  vdj.takara <- vdj.takara[vdj.takara$UM.ID!="UM22",]
  
  return(vdj.takara)
}

# Functions for plotting -------------------------------------------------------
format_for_ggsankey <- function(
    tcrs.split_ab.merged, 
    cols_include = c("cell", "celltype", "tissue_site", "antigen", "organism"),
    cols_plot = c("celltype", "tissue_site", "organism"), 
    subset_organism = NULL){
  
  lst.biopsy <- tcrs.split_ab.merged %>% separate_longer_delim(
    c(antigen, organism), delim = ";")
  
  orgamism_rename_map <- get_organism_rename_map(lst.biopsy)
  
  lst.biopsy$human_antigen <- lst.biopsy$antigen
  lst.biopsy$human_antigen[lst.biopsy$organism != "Homo sapiens"] <- NA
  
  for (i in 1:nrow(lst.biopsy)){
    if (lst.biopsy[i,]$organism != "Unknown species"){
      lst.biopsy[i,]$organism <- orgamism_rename_map$new[
        orgamism_rename_map$old==lst.biopsy[i,]$organism]
    }
  }
  
  lst <- lst.biopsy[,cols_include]
  if (!is.null(subset_organism)){
    lst <- lst[lst$organism %in% subset_organism,]
  }
  lst <- unique(lst)
  df <- lst %>% make_long(all_of(cols_plot))
  
  return(df)
}

plot_sankey <- function(df.sankey, space = NULL, size = 3, flow.alpha = .6){
  g <- ggplot(df.sankey, aes(x = x, next_x = next_x, node = node, 
                             next_node = next_node, 
                             fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = flow.alpha,
                node.color = NA, space = space, na.rm = T) +
    geom_sankey_label(size = size, color = "black", fill = "white", 
                      space = space, na.rm = T) +
    scale_fill_viridis_d() +
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5))
  return(g)
}

reorder_sankey_factors <- function(df.sankey, df, condition = c("TIL","biopsy")){
  df.sankey$node <- as.character(df.sankey$node)
  subject <- 
    unique(df.sankey$node[df.sankey$node %in% df$subject])[order(
      as.numeric(
        gsub(pattern="UM",replacement = "",unique(df.sankey$node[
          df.sankey$node %in% df$subject]))),
      decreasing = T)]
  antigen <- rev(unique(df.sankey$node[df.sankey$node %in% df$antigen]))
  organism <- c(setdiff(rev(unique(df.sankey$node[df.sankey$node %in% df$organism])),
                        c("Unknown species","Homo sapiens")),
                c("Unknown species","Homo sapiens"))
  fact_order <- c(subject, condition, organism, antigen)
  df.sankey$node <- factor(df.sankey$node, levels = fact_order)
  df.sankey$next_node <- factor(df.sankey$next_node, levels = fact_order)
  
  return(df.sankey)
}

reorder_sankey_factors_cluster <- function(df.sankey, df, condition = c("TIL","biopsy")){
  df.sankey$node <- as.character(df.sankey$node)
  cluster <- 
    unique(df.sankey$node[df.sankey$node %in% df$cluster])[order(
      as.numeric(
        gsub(pattern="^C",replacement = "",unique(
          df.sankey$node[df.sankey$node %in% df$cluster]))),
      decreasing = T)]
  antigen <- rev(unique(df.sankey$node[df.sankey$node %in% df$antigen]))
  organism <- c(setdiff(rev(unique(df.sankey$node[df.sankey$node %in% df$organism])),
                        c("Unknown species","Homo sapiens")),
                c("Unknown species","Homo sapiens"))
  fact_order <- c(cluster, condition, organism, antigen)
  df.sankey$node <- factor(df.sankey$node, levels = fact_order)
  df.sankey$next_node <- factor(df.sankey$next_node, levels = fact_order)
  
  return(df.sankey)
}
