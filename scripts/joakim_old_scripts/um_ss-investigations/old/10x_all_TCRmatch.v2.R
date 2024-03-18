library(readxl)
library(WriteXLS)
library(stringr)
library(pheatmap)
library(Seurat)

get_tcrs_split_ab_non_unique <- function(seu.filt.cd8_t, condition = "biopsy"){
  tcrs <- seu.filt.cd8_t@meta.data[which(seu.filt.cd8_t@meta.data$n_chains == 2),
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
  
  colnames(tcrs.split_ab) <- c("CDR3b","TRBV","TRBJ","CDR3a","subject:condition")
  
  head(tcrs.split_ab)
  
  tcrs.split_ab$clone_id <- apply(tcrs.split_ab, 1, function(x){
    paste0(x, collapse=";")
  })
  
  rownames(tcrs.split_ab) <- rownames(tcrs)
  return(tcrs.split_ab)
}

# Write input files for TCRmatch -----------------------------------------------

# Biopsies

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/"

seu.filt <- readRDS(paste0(outdir,"seu.filt.rda"))

tcrs.split_ab <- get_tcrs_split_ab_non_unique(seu.filt)

tcrs.split_ab$tcrmatch_input <- ""
tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$CDR3b
tcrs.split_ab$tcrmatch_input <- gsub(pattern="^C",replacement="",tcrs.split_ab$tcrmatch_input)
tcrs.split_ab$tcrmatch_input <- gsub(pattern="F$",replacement="",tcrs.split_ab$tcrmatch_input)

write.table(unique(tcrs.split_ab$tcrmatch_input), file = paste0(outdir,"seu.filt.tcrmatch_input.txt"),
            sep = "\t", col.names = F, quote = F, row.names = F)

# TILs

dat.til <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/v16/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
dat.til <- subset(dat.til, subset = project.name %in% c("G18-023"))

tcrs.split_ab.til <- get_tcrs_split_ab_non_unique(dat.til,"TIL")
head(tcrs.split_ab.til)

tcrs.split_ab.til$tcrmatch_input <- ""
tcrs.split_ab.til$tcrmatch_input <- tcrs.split_ab.til$CDR3b
tcrs.split_ab.til$tcrmatch_input <- gsub(pattern="^C",replacement="",tcrs.split_ab.til$tcrmatch_input)
tcrs.split_ab.til$tcrmatch_input <- gsub(pattern="F$",replacement="",tcrs.split_ab.til$tcrmatch_input)

write.table(unique(tcrs.split_ab.til$tcrmatch_input), file = paste0(outdir,"til.tcrmatch_input.txt"),
            sep = "\t", col.names = F, quote = F, row.names = F)

# Read TCRmatch output ---------------------------------------------------------
read_tcrmatch <- function(pth){
  tcrmatch_output <- read.table(pth,sep = "\t",header = T,fill = T)
  stopifnot(all(is.na(tcrmatch_output[,8])))
  tcrmatch_output[,8] <- NULL
  tcrmatch_output$original_tcrmatch_row <- rownames(tcrmatch_output)
  return(tcrmatch_output)
}

read_iedb <- function(pth){
  iedb <- data.table::fread(pth,
                            sep = "\t", header = F, data.table = F, drop = 7)
  colnames(iedb) <- c("trimmed_seq", "original_seq", "receptor_group", "epitopes",
                      "source_organisms", "source_antigens")
  iedb$trimmed_seq <- NULL
  iedb$original_seq <- NULL
  iedb$epitopes <- NULL
  iedb <- unique(iedb)

  return(iedb)
}

replace_internal_commas <- function(iedb, pth_iedb_external){
  iedb$source_antigens.original <- iedb$source_antigens
  # Replace from IEDB main file: identify internal commas as those where 
  # only one source organism, but commas present in source antigen
  replace_from_iedb_internal <- function(iedb){
    iedb_only_antigen_comma <- unique(iedb[!grepl(",",iedb$source_organisms) & 
                                             grepl(",",iedb$source_antigens),])
    iedb_only_antigen_comma$source_antigens_internal_commas_removed <- gsub(
      ",","_",iedb_only_antigen_comma$source_antigens)
    
    for (ant in unique(iedb_only_antigen_comma$source_antigens)){
      iedb$source_antigens <- gsub(
        pattern=ant,
        replacement = unique(
          iedb_only_antigen_comma$source_antigens_internal_commas_removed[
            iedb_only_antigen_comma$source_antigens==ant]),
        iedb$source_antigens,fixed = T)
    }
    return(iedb)
  }
  
  replace_from_iedb_external <- function(
    iedb_both_comma, pth_iedb_external){
    
    tmp <- iedb_both_comma
    idx <- which(unlist(lapply(tmp$source_organisms,function(x){
      length(unlist(strsplit(x,",")))})) != 
      unlist(lapply(tmp$source_antigens,function(x){
        length(unlist(strsplit(x,",")))})))
    
    iedb_tcell_full <- read.table(
      pth_iedb_external,
      skip = 1, header = T, comment.char = "", sep=",")
    
    cols_relevant <- c("Description",	"Antigen Name","Parent Protein",
                       "Organism Name","Parent Species","Name",
                       "Immunogen Description","Immunogen Source Molecule Name",
                       "Immunogen protein parent Name","Immunogen Organism Name",
                       "Immunogen Organism Species","Immunogen Description",
                       "Immunogen Source Molecule Name","Protein Parent Name",
                       "Immunogen Organism Name","Immunogen Organism Species",
                       "Antigen Description","Antigen Source Molecule Name",
                       "Protein Parent Name","Antigen Organism Name",
                       "Organism Species Name")
    cols_relevant <- gsub(" ",".",cols_relevant)
    cols_relevant <- c(unique(cols_relevant), paste0(
      names(which(table(cols_relevant)>1)),".1"))
    
    iedb_tcell_full.relevant <- iedb_tcell_full[,cols_relevant]
    iedb_tcell_full.relevant <- unique(iedb_tcell_full.relevant)
    
    antigen_cols <- c("Antigen.Name","Parent.Protein","Immunogen.Description",
                      "Antigen.Source.Molecule.Name","Protein.Parent.Name.1")
    
    antigens <- unique(c(iedb_tcell_full.relevant[,antigen_cols[1]],
             iedb_tcell_full.relevant[,antigen_cols[2]],
             iedb_tcell_full.relevant[,antigen_cols[3]],
             iedb_tcell_full.relevant[,antigen_cols[4]],
             iedb_tcell_full.relevant[,antigen_cols[5]]))
    
    nms <- unique(antigens[grep(",",antigens)])
    
    iedb_internal_commas_2 <- data.frame(source_antigens=nms, 
               source_antigens_internal_commas_removed = gsub(",", "_", nms))
    
    x <- tmp[idx,]
    for (ant in unique(iedb_internal_commas_2$source_antigens)){
      x$source_antigens <- gsub(pattern=ant,
                                   replacement = unique(
                                     iedb_internal_commas_2$source_antigens_internal_commas_removed[
                                       iedb_internal_commas_2$source_antigens==ant]),
                                   x$source_antigens, fixed = T)
    }
    tmp[idx,] <- x
    
    stopifnot(all(unlist(lapply(tmp$source_organisms,function(x){length(unlist(strsplit(x,",")))})) == 
      unlist(lapply(tmp$source_antigens,function(x){length(unlist(strsplit(x,",")))}))))
    
    return(tmp)
  }
  
  iedb <- replace_from_iedb_internal(iedb)
  
  iedb_not_both_comma <- iedb[!(grepl(",",iedb$source_organisms) & grepl(",",iedb$source_antigens)),]
  iedb_both_comma <- iedb[grepl(",",iedb$source_organisms) & grepl(",",iedb$source_antigens),]
  
  iedb_both_comma <- replace_from_iedb_external(iedb_both_comma, pth_iedb_external)
  
  iedb <- rbind(iedb_not_both_comma, iedb_both_comma)
  iedb$source_antigens <- gsub(pattern=",",replacement=";",iedb$source_antigens)
  
  return(iedb)
}

f.check.string.equality <- function(s1, s2) {
  resEqualCheck = apply(do.call(rbind, strsplit(c(s1, s2), "")), 2, function(x) {identical(x[1],x[2])})
  return(resEqualCheck)
}

clean_tcrmatch_output <- function(tcrmatch_output, iedb){
  tmp <- tcrmatch_output
  tmp$antigen.original <- tmp$antigen
  tmp$antigen <- ""
  for (i in 1:nrow(tmp)){
    receptor_groups <- as.character(tmp$receptor_group[i])
    receptor_groups <- as.numeric(unique(unlist(strsplit(receptor_groups,","))))
    for (j in 1:length(receptor_groups)){
      if (tmp$antigen[i] == ""){
        tmp$antigen[i] <- iedb$source_antigens[
          iedb$receptor_group == receptor_groups[j]]
      } else {
        tmp$antigen[i] <- paste0(tmp$antigen[i],"|",iedb$source_antigens[
          iedb$receptor_group == receptor_groups[j]])
      }
    }
  }
  
  # Validate assumptions
  tmp2 <- tmp[tmp[,"antigen"]!=tmp[,"antigen.original"],]
  res <- list()
  for (i in 1:nrow(tmp2)){
    res[[i]] <- f.check.string.equality(tmp2[i,"antigen"],tmp2[i,"antigen.original"])
    res[[i]] <- unlist(strsplit(tmp2[i,"antigen"],""))[!res[[i]]]
  }
  stopifnot(all(unique(unlist(res)) %in% c("_",";","|")))
  
  tmp$organism <- gsub(",",";",tmp$organism)
  
  return(tmp)
}

add_anatimic_location_meta <- function(seu.filt.cd8_t, pth){
  library(readxl)
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
  
  orgamism_rename_map <- data.frame(new=x,old=y)
}

iedb <- read_iedb(pth = "~/nimbus/data/proj/um_ss/Pipelines/TCRMatch/data/IEDB_data.tsv")
iedb <- replace_internal_commas(iedb, pth_iedb_external = "~/proj/um_ss/Investigations/data/TCR_databases/IEDB/tcell_full_v3.csv")

tcrmatch_output <- read_tcrmatch(pth = "~/nimbus/data/proj/um_ss/Investigations/seurat/results/v16/seu.filt.tcrmatch_output_2.txt")
tcrmatch_output <- clean_tcrmatch_output(tcrmatch_output, iedb)

tcrmatch_output.til <- read_tcrmatch(pth = "~/proj/um_ss/Investigations/seurat/results/v16/til.tcrmatch_output_2.txt")
tcrmatch_output.til <- clean_tcrmatch_output(tcrmatch_output.til, iedb)

# ggsankey ---------------------------------------------------------------------
subfun <- function(x){
  if (any(grepl("[",x,fixed=T))){
    # Check that all antigens with [ somewhere in the name all end with ]
    stopifnot(substr(x, nchar(x), nchar(x)) == "]")
    x <- str_split_fixed(x," *\\[",2)[,1]
  }
  return(x)
}

rename_organisms <- function(tcrs.split_ab.merged){
  orgamism_rename_map <- get_organism_rename_map(tcrs.split_ab.merged)
  for (i in 1:nrow(tcrs.split_ab.merged)){
    y <- unlist(strsplit(tcrs.split_ab.merged[i,]$organism,";"))
    for (j in 1:length(y)){
      if (y[j] != "Unknown species"){
        y[j] <- orgamism_rename_map$new[orgamism_rename_map$old == y[j]]
      }
    }
    tcrs.split_ab.merged[i,]$organism <- paste0(y,collapse=";")
  }
  return(tcrs.split_ab.merged)
}

add_metadata <- function(tcrmatch_output, seu.filt, seu.filt.cd8_t){
  tcrs.split_ab <- get_tcrs_split_ab_non_unique(seu.filt)
  
  tcrs.split_ab$tcrmatch_input <- ""
  tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$CDR3b
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "^C", replacement = "", tcrs.split_ab$tcrmatch_input)
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "F$", replacement = "", tcrs.split_ab$tcrmatch_input)
  
  tcrs.split_ab$cell <- rownames(tcrs.split_ab)
  
  tcrs.split_ab.merged <- merge(tcrs.split_ab, tcrmatch_output, 
                                by.x="tcrmatch_input", 
                                by.y="input_sequence")#, all.x=T)
  tcrs.split_ab.merged$epitope[is.na(tcrs.split_ab.merged$epitope)] <- ""
  tcrs.split_ab.merged$antigen[is.na(tcrs.split_ab.merged$antigen)] <- ""
  tcrs.split_ab.merged$organism[is.na(tcrs.split_ab.merged$organism)] <- ""
  
  df <- data.frame(seu.filt.cd8_t@active.ident)
  colnames(df) <- "celltype"
  tcrs.split_ab.merged <- merge(tcrs.split_ab.merged,df,by.x="cell",by.y="row.names")
  
  seu.filt.cd8_t <- add_anatimic_location_meta(
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
  tcrs.split_ab$tcrmatch_input <- gsub(pattern = "F$", replacement = "", tcrs.split_ab$tcrmatch_input)
  
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
  
  nms_dup <- names(which(table(unique(tcrs.split_ab.merged[,c("antigen","organism")])$antigen) > 1))
  if (length(nms_dup)>0){
    idx <- tcrs.split_ab.merged$antigen %in% nms_dup
    tcrs.split_ab.merged[idx,]$antigen <- paste0(
      tcrs.split_ab.merged[idx,]$antigen," (", tcrs.split_ab.merged[idx,]$organism,")")
  }
  
  return(tcrs.split_ab.merged)
}

format_for_ggsankey <- function(tcrs.split_ab.merged, 
  cols_include = c("cell", "celltype", "tissue_site", "antigen", "organism"),
  cols_plot = c("celltype", "tissue_site", "organism"), subset_organism = NULL){
  
  lst.biopsy <- tcrs.split_ab.merged[,cols_include]
  
  lst.biopsy <- lst.biopsy %>% separate_longer_delim(c(antigen, organism), delim = ";")
  #lst.biopsy <- unique(lst.biopsy)
  
  orgamism_rename_map <- get_organism_rename_map(lst.biopsy)
  
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
  #lst <- unique(lst)
  #lst <- lst[lst$Species %in% names(which(table(lst$Species)>10)),]
  lst <- unique(lst)
  df <- lst %>% make_long(all_of(cols_plot))
  
  #celltype <- names(sort(table(lst$celltype),decreasing = F))
  #tissue_site <- names(sort(table(lst$tissue_site),decreasing = F))
  #antigen <- unique(sort(lst$antigen,decreasing = T))
  #organism <- unique(sort(lst$Species,decreasing = T))
  
  #fact_order <- c(celltype, tissue_site, antigen, organism)
  #fact_order <- c(celltype, tissue_site, organism)
  #df$node <- factor(df$node, levels = fact_order)
  #df$next_node <- factor(df$next_node, levels = fact_order)
  
  return(df)
}

plot_sankey <- function(df.sankey){
  g <- ggplot(df.sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
                      fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = .6,
                node.color = NA) +
    geom_sankey_label(size = 3, color = "black", fill = "white") +
    scale_fill_viridis_d() +
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5))
  return(g)
}

plot_heatmap <- function(df, subset_organism = NULL, invert = F,
                         group_top = c("Homo sapiens", "Unknown species"),
                         norm = F){
  lst.biopsy <- df %>% separate_longer_delim(c(antigen, organism), delim = ";")
  
  if (!is.null(subset_organism)){
    if (invert){
      lst.biopsy <- lst.biopsy[! lst.biopsy$organism %in% subset_organism,]
    } else {
      lst.biopsy <- lst.biopsy[lst.biopsy$organism %in% subset_organism,]  
    }
  }
  
  tmp <- unique(lst.biopsy[,c("antigen","organism")])
  names_dup <- names(which(table(tmp$antigen) > 1))
  if (length(names_dup) > 0){
    idx <- lst.biopsy$antigen %in% names_dup
    lst.biopsy[idx,]$antigen <- paste0(lst.biopsy[idx,]$antigen," (",lst.biopsy[idx,]$organism,")")
  }
  annot <- unique(lst.biopsy[,c("antigen","organism")])
  rownames(annot) <- annot$antigen
  
  if (norm){
    # For each CDR3b that maps to multiple antigens, divide by the number it maps to
    unique_samples <- unique(lst.biopsy$`subject:condition`)
    unique_antigens <- unique(lst.biopsy$antigen)
    x <- matrix(NA, nrow = length(unique_antigens), ncol = length(unique_samples),
                    dimnames = list(unique_antigens, unique_samples))
    for (i in 1:length(unique_antigens)){
      for (j in 1:length(unique_samples)){
        idx <- lst.biopsy$antigen == unique_antigens[i] & 
          lst.biopsy$`subject:condition` == unique_samples[j]
        tmp <- lst.biopsy[idx,]
        if (nrow(tmp)>=1){
          unique_CDR3b <- unique(tmp$CDR3b)
          res <- list()
          for (CDR3b in unique_CDR3b){
            res[[CDR3b]] <- length(unique(lst.biopsy$antigen[lst.biopsy$CDR3b==CDR3b]))
          }
          res <- as.data.frame(unlist(res))
          colnames(res) <- "n_different_antigens"
          tmp <- merge(tmp,res,by.x="CDR3b",by.y="row.names")
          x[i,j] <- sum(1/tmp$n_different_antigens)
        } else {
          x[i,j] <- nrow(tmp)
        }
      }
    }
    number_format <- "%.2f"
  } else {
    x <- table(lst.biopsy$antigen, lst.biopsy$`subject:condition`)
    number_format <- "%.0f"
  }
  
  annot_hum <- data.frame()
  annot_nonhum <- data.frame()
  if (any(annot$organism %in% group_top)){
    annot_hum <- annot[annot$organism %in% group_top,]
    annot_hum$organism <- factor(annot_hum$organism, level = group_top)
    annot_hum <- annot_hum[order(annot_hum$organism),]
  }
  if (any(! annot$organism %in% group_top)){
    annot_nonhum <- annot[! annot$organism %in% group_top,]
    annot_nonhum <- annot_nonhum[order(annot_nonhum$organism),]
  }
  if (nrow(annot_hum) > 0 && nrow(annot_nonhum) > 0){
    annot <- rbind(annot_hum,annot_nonhum)
  } else if (nrow(annot_hum) > 0){
    annot <- annot_hum
  } else if (nrow(annot_nonhum) > 0){
    annot <- annot_nonhum
  }
  
  annot$antigen <- NULL
  nms <- rownames(annot)[order(annot$organism)]
  snames <- colnames(x)[order(as.numeric(gsub("[aA-zZ:]","",colnames(x))))]
  x <- x[nms,snames]
  x[is.na(x)] <- NA
  x_non_na <- x
  
  p <- pheatmap(x, border_color = NA, display_numbers = round(x_non_na,2), 
                number_format = number_format, number_color = "#000000",
                annotation_row = annot, cluster_rows = F, cluster_cols = F,
                na_col = "white", colorRampPalette(c("white", "red"))(100))
  return(p)
}

# Biopsy

seu.filt.cd8_t <- readRDS(paste0(outdir,"seu.filt.cd8_t.rda"))

tcrs.split_ab.merged <- add_metadata(tcrmatch_output, seu.filt, seu.filt.cd8_t)

df.sankey <- format_for_ggsankey(tcrs.split_ab.merged = tcrs.split_ab.merged, 
                          cols_include = c("cell", "celltype", "tissue_site",
                                           "subject:condition", "antigen", "organism"), 
                          cols_plot = c("subject:condition", "celltype", "antigen"),
                          subset_organism = c("Homo sapiens","Unknown species"))
plot_sankey(df.sankey)

df.sankey <- format_for_ggsankey(tcrs.split_ab.merged = tcrs.split_ab.merged, 
                                 cols_include = c("cell", "celltype", "tissue_site", 
                                                  "subject:condition",  "antigen", "organism"), 
                                 cols_plot = c("subject:condition", "celltype", "organism"))
plot_sankey(df.sankey)

# TIL

tcrs.split_ab.til.merged <- add_metadata_til(tcrmatch_output.til, dat.til)

df.sankey <- format_for_ggsankey(tcrs.split_ab.merged = tcrs.split_ab.til.merged, 
                                 cols_include = c("cell","subject:condition", "antigen", "organism"),
                                 cols_plot = c("subject:condition", "antigen", "organism"))
plot_sankey(df.sankey)

# TIL and Biopsy combined

df.sankey.biopsy <- format_for_ggsankey(tcrs.split_ab.merged = tcrs.split_ab.merged, 
                                 cols_include = c("cell","subject:condition",  "antigen", "organism"), 
                                 cols_plot = c("subject:condition", "antigen", "organism"))

df.sankey.til <- format_for_ggsankey(tcrs.split_ab.merged = tcrs.split_ab.til.merged, 
                                 cols_include = c("cell","subject:condition", "antigen", "organism"),
                                 cols_plot = c("subject:condition", "antigen", "organism"))

tcrs.split_ab.til.merged <- merge(tcrs.split_ab.til.merged,
      unique(tcrs.split_ab.merged[,c("subject","tissue_site")]),
      by="subject",all.x = T,all.y = F)

tcrs.split_ab.til.merged$tissue_site[tcrs.split_ab.til.merged$subject=="UM10"] <- "Liver"
tcrs.split_ab.til.merged$tissue_site[tcrs.split_ab.til.merged$subject=="UM22"] <- "Subcutaneous"

cnames <- intersect(colnames(tcrs.split_ab.til.merged),colnames(tcrs.split_ab.merged))
df <- rbind(tcrs.split_ab.til.merged[,cnames],tcrs.split_ab.merged[,cnames])

df$human_antigen <- df$antigen
df$human_antigen[df$organism != "Homo sapiens"] <- NA

# df.sankey.biopsy_til <- format_for_ggsankey(tcrs.split_ab.merged = df,
#                           cols_include = c("cell","subject:condition", 
#                                            "antigen", "organism", "human_antigen"),
#                           cols_plot = c("subject:condition", "organism","human_antigen"))

df.sankey.biopsy_til <- format_for_ggsankey(tcrs.split_ab.merged = df,
                                            cols_include = c("cell","subject","condition", "tissue_site",
                                                             "antigen", "organism", "human_antigen"),
                                            cols_plot = c("subject","tissue_site", "condition",
                                                          "organism","human_antigen"))

df.sankey <- df.sankey.biopsy_til

# df.sankey$col <- "black"
#  
# idx <- which(str_split_fixed(df.sankey$node,":",2)[,1] %in% df$subject & 
#                 df.sankey$x=="subject" & df.sankey$next_x=="organism")
# df.sankey$col[idx] <- df.sankey$node[idx]
# 
# idx <- which(df.sankey$next_node %in% df$subject & 
#                df.sankey$x=="condition" & df.sankey$next_x=="subject")
# df.sankey$col[idx] <- df.sankey$next_node[idx]

# subject_condition <- 
#   unique(df.sankey$node[df.sankey$node %in% df$`subject:condition`])[order(
#     as.numeric(str_split_fixed(
#     gsub(pattern="UM",replacement = "",unique(df.sankey$node[df.sankey$node %in% df$`subject:condition`])),
#     ":",2)[,1]),decreasing = T)]
subject <- 
  unique(df.sankey$node[df.sankey$node %in% df$subject])[order(
    as.numeric(
      gsub(pattern="UM",replacement = "",unique(df.sankey$node[df.sankey$node %in% df$subject]))),
    decreasing = T)]
condition <- c("TIL","biopsy")
tissue_site <- rev(unique(df.sankey$node[df.sankey$node %in% df$tissue_site]))
antigen <- rev(unique(df.sankey$node[df.sankey$node %in% df$antigen]))
organism <- c(setdiff(rev(unique(df.sankey$node[df.sankey$node %in% df$organism])),
                      c("Unknown species","Homo sapiens")),
              c("Unknown species","Homo sapiens"))
fact_order <- c(subject, tissue_site, condition, organism, antigen)
df.sankey$node <- factor(df.sankey$node, levels = fact_order)
df.sankey$next_node <- factor(df.sankey$next_node, levels = fact_order)

ggplot(df.sankey, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
                           fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = NA) +
  geom_sankey_label(size = 3, color = "black", fill = "white") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))

# Make heatmap -----------------------------------------------------------------
plot_heatmap(df)
plot_heatmap(df, norm = T)
plot_heatmap(df, subset_organism = c("Homo sapiens"))
plot_heatmap(df, subset_organism = c("Homo sapiens"), norm = T)
plot_heatmap(df, subset_organism = c("Homo sapiens", "Unknown species"))
plot_heatmap(df, subset_organism = c("Homo sapiens", "Unknown species"), norm=T)
plot_heatmap(df, subset_organism = c("Homo sapiens", "Unknown species"), 
             invert = T, group_top = c())
plot_heatmap(df, subset_organism = c("Homo sapiens", "Unknown species"), 
             invert = T, group_top = c(), norm = T)

# Attempt to represent as graph ------------------------------------------------
library(igraph)

make_g <- function(df){
  unique_CDR3b <- unique(df$CDR3b)
  unique_antigen <- unique(df$antigen)
  adjmat <- matrix(0, nrow = length(unique_CDR3b), ncol = length(unique_antigen), 
                   dimnames = list(unique_CDR3b, unique_antigen))
  for (CDR3b in unique_CDR3b){
    for (antigen in unique_antigen){
      adjmat[CDR3b,antigen] <- any(df$CDR3b == CDR3b & df$antigen == antigen)
    }
  }
  g <- graph_from_incidence_matrix(adjmat, directed = T, mode = "out")
  return(g)
}

make_g2 <- function(df){
  unique_antigen <- unique(df$antigen)
  adjmat2 <- matrix(0, nrow = length(unique_antigen), ncol = length(unique_antigen), 
                   dimnames = list(unique_antigen,unique_antigen))
  options(warn=2)
  for (antigen_1 in unique_antigen){
    for (antigen_2 in unique_antigen){
      adjmat2[antigen_1,antigen_2] <- unique(
        df$organism[df$antigen == antigen_1]) == unique(df$organism[df$antigen == antigen_2])
    }
  }
 
  g2 <- graph_from_adjacency_matrix(adjmatrix = adjmat2, mode = "undirected", diag = F)
  return(g2)
}

make_g3 <- function(df){
  unique_CDR3b <- unique(df$CDR3b)
  adjmat3 <- matrix(0, nrow = length(unique_CDR3b), ncol = length(unique_CDR3b),
                    dimnames = list(unique_CDR3b,unique_CDR3b))
  options(warn=2)
  for (CDR3b_1 in unique_CDR3b){
    for (CDR3b_2 in unique_CDR3b){
      adjmat3[CDR3b_1,CDR3b_2] <- min(length(intersect(
        df$subject[df$CDR3b == CDR3b_1],
        df$subject[df$CDR3b == CDR3b_2])),1)
    }
  }
  
  g3 <- graph_from_adjacency_matrix(adjmatrix = adjmat3, mode = "undirected", diag = F)
  
  return(g3)
}

merge_to_bipartite <- function(g,g2, col_1 = "yellow", col_2 = "red", species_g, species_g2){
  g.df <- as_data_frame(g, "vertices")
  g.df <- g.df[,c("type","name")]
  g.df$type <- NULL
  g2.df <- as_data_frame(g2, "vertices")
  
  attrs <- rbind(g.df, g2.df) %>% unique()
  el <- rbind(as_data_frame(g), as_data_frame(g2))
  
  new_g <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)
  V(new_g)$frame.color <- NA
  V(new_g)$color[! V(new_g)$name %in% V(g2)$name] <- col_1
  V(new_g)$color[V(new_g)$name %in% V(g2)$name] <- col_2
  
  V(new_g)$species <- ifelse(V(new_g)$color == col_1, species_g, species_g2)
  
  return(new_g)
}

get_color_edges_bipartite <- function(new_g, non_organism_col = "red"){
  edge.start <- ends(new_g, es=E(new_g), names=F)[,1]
  edge.end <- ends(new_g, es=E(new_g), names=F)[,2]
  edge.col_1 <- V(new_g)$color[edge.start]
  edge.col_2 <- V(new_g)$color[edge.end]
  edge.col <- ifelse(edge.col_1 == edge.col_2, edge.col_1, non_organism_col)
  return(edge.col)
}

resize_degree_incoming <- function(new_g, species_g, species_g2){
  el <- as_data_frame(new_g)
  
  nm.species_g <- V(new_g)$name[V(new_g)$species==species_g]
  nm.species_g2 <- V(new_g)$name[V(new_g)$species==species_g2]
  
  stopifnot(nrow(el[el$to %in% nm.species_g & el$from %in% nm.species_g2,])==0)
  
  deg <- el[el$from %in% nm.species_g & el$to %in% nm.species_g2,]
  deg$degree <- NA
  for (antigen in unique(deg$to)){
    deg$degree[deg$to==antigen] <- length(unique(deg$from[deg$to==antigen]))
  }
  
  V(new_g)$size <- 1
  stopifnot(nrow(unique(deg[,c("to","degree")])) == length(unique(deg$to)))
  for (antigen in unique(deg$to)){
    V(new_g)$size[V(new_g)$name==antigen] <- unique(deg$degree[deg$to==antigen])
  }
  
  V(new_g)$size <- log10(V(new_g)$size*10)
  
  return(new_g)
}

color_by_organism <- function(new_g, df, non_organism_col = "red"){
  library(RColorBrewer)
  
  V(new_g)$type <- V(new_g)$color == non_organism_col
  
  tmp_col <- rep(non_organism_col,length(V(new_g)$color))
  for (antigen in intersect(unique(V(new_g)$name),colnames(adjmat))){
    tmp_col[V(new_g)$name==antigen] <- unique(df[df$antigen==antigen,]$organism)
  }
  
  n <- length(unique(tmp_col))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  col_vector <- setNames(col_vector,setdiff(unique(tmp_col),non_organism_col))
  for (nm in names(col_vector)){
    tmp_col[tmp_col == nm] <- col_vector[nm]
  }
  V(new_g)$color <- tmp_col
  
  return(new_g)
}

g <- make_g(df)
g2 <- make_g2(df)
new_g <- merge_to_bipartite(g,g2,species_g="CDR3b", species_g2="antigen")
new_g <- resize_degree_incoming(new_g, species_g="CDR3b", species_g2="antigen")
new_g <- color_by_organism(new_g, df, non_organism_col = "gray")
edge.col <- get_color_edges_bipartite(new_g, non_organism_col = "gray")

plot(new_g, vertex.label=NA, edge.arrow.size=.4, vertex.size = V(new_g)$size*2,
     layout = layout_with_fr(new_g), edge.color=edge.col)



new_g.df <- as_data_frame(new_g, "vertices")
new_g.df <- new_g.df[,c("type","name")]
new_g.df$type <- NULL
g3.df <- as_data_frame(g3, "vertices")

attrs <- rbind(new_g.df, g3.df) %>% unique()
el <- rbind(as_data_frame(new_g), as_data_frame(g3))

g_all <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)

stopifnot(identical(V(g_all)$name,V(new_g)$name))
V(g_all)$species <- V(new_g)$species
V(g_all)$color <- V(new_g)$color

g_all <- resize_degree_incoming(g_all, species_g="CDR3b", species_g2="antigen")

V(g_all)$frame.color <- V(g_all)$color

edge.col <- get_color_edges_bipartite(g_all, non_organism_col = "white")

V(g_all)$color[V(g_all)$color=="gray"] <- "black"
V(g_all)$frame.color <- V(g_all)$color

plot(g_all, vertex.label=NA, edge.arrow.size=.4, vertex.size = V(new_g)$size*2,
     layout = layout_with_fr(g_all, niter = 10^4, grid = "nogrid"), 
     edge.color=edge.col, edge.curved=1)

library("ForceAtlas2")

l <- layout.forceatlas2(g_all, directed=F, iterations = 3000, 
                        linlog = FALSE, pos = NULL, nohubs = FALSE, 
                        k = 4000, gravity=1, ks=1, ksmax=50, delta = 1,  
                        center=NULL, tolerance = 0.1, dim = 2,
                        plotstep=10, plotlabels=F)

edge.start <- ends(g_all, es=E(g_all), names=F)[,1]
edge.end <- ends(g_all, es=E(g_all), names=F)[,2]
edge.col_1 <- V(g_all)$species[edge.start]
edge.col_2 <- V(g_all)$species[edge.end]
edge.col <- c()
for (i in 1:length(edge.col_1)){
  edge.col[i] <- ifelse(edge.col_1[i] == "CDR3b" & edge.col_2[i] == "CDR3b", "gray","")
  edge.col[i] <- ifelse(edge.col_1[i] == "CDR3b" & edge.col_2[i] == "antigen", V(g_all)$color[edge.end][i],edge.col[i])
  edge.col[i] <- ifelse(edge.col_1[i] == "antigen" & edge.col_2[i] == "CDR3b", V(g_all)$color[edge.start][i],edge.col[i])
  edge.col[i] <- ifelse(edge.col_1[i] == "antigen" & edge.col_2[i] == "antigen", V(g_all)$color[edge.end][i],edge.col[i])
}

plot(g_all, vertex.label=NA, edge.arrow.size=.4, vertex.size = V(new_g)$size*2,
     layout = l, edge.color=edge.col, edge.curved=1)

