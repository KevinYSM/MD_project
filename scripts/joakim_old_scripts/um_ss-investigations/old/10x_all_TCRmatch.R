library(readxl)
library(WriteXLS)
library(stringr)
library(pheatmap)
library(Seurat)

get_tcrs_split_ab_non_unique <- function(seu.filt.cd8_t){
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
  
  colnames(tcrs.split_ab) <- c("CDR3b","TRBV","TRBJ","CDR3a","subject:condition")
  
  head(tcrs.split_ab)
  
  tcrs.split_ab$clone_id <- apply(tcrs.split_ab, 1, function(x){
    paste0(x, collapse=";")
  })
  
  rownames(tcrs.split_ab) <- rownames(tcrs)
  return(tcrs.split_ab)
}

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/"

seu.filt <- readRDS(paste0(outdir,"seu.filt.rda"))

tcrs.split_ab <- get_tcrs_split_ab_non_unique(seu.filt)

tcrs.split_ab$tcrmatch_input <- ""
tcrs.split_ab$tcrmatch_input <- tcrs.split_ab$CDR3b
tcrs.split_ab$tcrmatch_input <- gsub(pattern="^C",replacement="",tcrs.split_ab$tcrmatch_input)
tcrs.split_ab$tcrmatch_input <- gsub(pattern="F$",replacement="",tcrs.split_ab$tcrmatch_input)

write.table(unique(tcrs.split_ab$tcrmatch_input), file = paste0(outdir,"seu.filt.tcrmatch_input.txt"),
            sep = "\t", col.names = F, quote = F, row.names = F)

# Read TCRmatch output ---------------------------------------------------------
read_tcrmatch <- function(pth){
  tcrmatch_output <- read.table(pth,sep = "\t",header = T,fill = T)
  stopifnot(all(is.na(tcrmatch_output[,8])))
  tcrmatch_output[,8] <- NULL
  tcrmatch_output$original_tcrmatch_row <- rownames(tcrmatch_output)
  return(tcrmatch_output)
}

pth <- "~/nimbus/data/proj/um_ss/Investigations/seurat/results/v16/seu.filt.tcrmatch_output_2.txt"

tcrmatch_output <- read_tcrmatch(pth)



tcrs.split_ab$cell <- rownames(tcrs.split_ab)

tcrs.split_ab.merged <- merge(tcrs.split_ab, tcrmatch_output, 
                               by.x="tcrmatch_input", 
                               by.y="input_sequence")#, all.x=T)
tcrs.split_ab.merged$epitopes[is.na(tcrs.split_ab.merged$epitopes)] <- ""
tcrs.split_ab.merged$antigen[is.na(tcrs.split_ab.merged$antigen)] <- ""
tcrs.split_ab.merged$source_organism[is.na(tcrs.split_ab.merged$source_organism)] <- ""

for (i in 1:nrow(tcrs.split_ab.merged)){
  tcrs.split_ab.merged$antigen[i] <- ifelse(
    is.na(tcrs.split_ab.merged$antigen[i]),
    NA,
    paste0(unique(unlist(strsplit(tcrs.split_ab.merged$antigen[i],","))),collapse = ";"))
  tcrs.split_ab.merged$source_organism[i] <- ifelse(
    is.na(tcrs.split_ab.merged$source_organism[i]),
    NA,
    paste0(unique(unlist(strsplit(tcrs.split_ab.merged$source_organism[i],","))),collapse = ";"))
  tcrs.split_ab.merged$epitopes[i] <- ifelse(
    is.na(tcrs.split_ab.merged$epitopes[i]),
    NA,
    paste0(unique(unlist(strsplit(tcrs.split_ab.merged$epitopes[i],","))),collapse = ";"))
}

# tab <- table(tcrs.split_ab.merged$`subject:condition`,tcrs.split_ab.merged$antigen)
# pheatmap(tab,border_color = NA,fontsize_col = 3)
# 
# tab <- table(tcrs.split_ab.merged$`subject:condition`,tcrs.split_ab.merged$source_organism)
# pheatmap(tab,border_color = NA,fontsize_col = 3)
# 
# x <- tcrs.split_ab.merged[grep("Homo",tcrs.split_ab.merged$source_organism),]
# x$antigen <- gsub(pattern="\\[.+\\]",replacement = "",x$antigen, ignore.case = T,perl=T)
# tab <- table(x$`subject:condition`,x$antigen)
# pheatmap(tab,border_color = NA,fontsize_col = 5)

tcrs.split_ab.merged$source_organism <- gsub(pattern="Yellow fever virus 17D (Yellow fever virus (STRAIN 17D))",
     replacement = "YFV",
     tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Human immunodeficiency virus 1 (human immunodeficiency virus 1 HIV-1)",
                                             replacement = "HIV-1",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Influenza A virus (A/Memphis/4/1973(H3N2)) (Influenza A virus (A/Memphis/4/73(H3N2)))",
                                             replacement = "Influenza A",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Influenza A (A/Puerto Rico/8/1934(H1N1)) (Influenza A (A/PR/8/1934(H1N1)))",
                                             replacement = "Influenza A",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Influenza A virus H3N2 (A/Resvir-9 (H3N2))",
                                             replacement = "Influenza A",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Influenza A virus",
                                             replacement = "Influenza A",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Murid betaherpesvirus 1 (Murine cytomegalovirus)",
                                             replacement = "Murine CMV",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Triticum aestivum (Canadian hard winter wheat)",
                                             replacement = "Wheat",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Human herpesvirus 4 strain B95-8 (Epstein-Barr virus (strain B95-8))",
                                             replacement = "EBV",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Human herpesvirus 4 (Epstein Barr virus)",
                                             replacement = "EBV",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Human T-cell leukemia virus type I (Human adult T-cell leukemia virus)",
                                             replacement = "HTLV-1",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Human rotavirus strain WA (Human rotavirus serotype 1 / strain WA)",
                                             replacement = "Rotavirus",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Hepatitis B virus (hepatitis B virus (HBV))",
                                             replacement = "HBV",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Dengue virus 3 (Dengue virus serotype 3)",
                                             replacement = "Dengue virus 3",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Dengue virus 1 (dengue type 1 D1 virus)",
                                             replacement = "Dengue virus 1",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Dengue virus 2 Thailand/16681/84 (Dengue virus type 2 (strain 16681))",
                                             replacement = "Dengue virus 2",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Dengue virus 2 (dengue 2 virus DEN-2)",
                                             replacement = "Dengue virus 2",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="SARS coronavirus Urbani (SARS-CoV (Urbani strain))",
                                             replacement = "SARS-CoV (Urbani strain)",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Mycobacterium tuberculosis H37Rv (Mycobacterium tuberculosis str. H37Rv)",
                                             replacement = "Tuberculosis",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Mycobacterium tuberculosis",
                                             replacement = "Tuberculosis",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Murid herpesvirus 1 deltaMS94.5",
                                             replacement = "Murid herpesvirus 1",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Lymphocytic choriomeningitis mammarenavirus (Lymphocytic choriomeningitis virus)",
                                             replacement = "LCMV",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Human herpesvirus 5 strain AD169 (Human cytomegalovirus (strain AD169))",
                                             replacement = "CMV",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Human herpesvirus 5 strain Towne (Human cytomegalovirus (strain Towne))",
                                             replacement = "CMV",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Human herpesvirus 5 (Human cytomegalovirus)",
                                             replacement = "CMV",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Gallus gallus (chicken)",
                                             replacement = "Chicken",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Mus musculus (mouse)",
                                             replacement = "Mouse",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Homo sapiens (human)",
                                             replacement = "Human",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Hepatitis C virus subtype 1a (Hepatitis C virus type 1a)",
                                             replacement = "Hepatitis C",
                                             tcrs.split_ab.merged$source_organism,fixed = T)
tcrs.split_ab.merged$source_organism <- gsub(pattern="Hepatitis C virus",
                                             replacement = "Hepatitis C",
                                             tcrs.split_ab.merged$source_organism,fixed = T)

tcrs.split_ab.merged.best <- list()
for (cell in unique(tcrs.split_ab.merged$cell)){
  rw <- tcrs.split_ab.merged[tcrs.split_ab.merged$cell==cell,]
  rw <- rw[rw$score==max(rw$score),]
  rw <- rw[,c("epitopes","antigen","source_organism","score")]
  for (cname in colnames(rw)){
    rw[,cname] <- paste0(unique(rw[,cname]),collapse = "|")
  }
  rw <- unique(rw)
  rw$cell <- cell
  tcrs.split_ab.merged.best[[cell]] <- rw
}

tcrs.split_ab.merged.best <- do.call("rbind",tcrs.split_ab.merged.best)

meta <- seu.filt@meta.data
meta.merged <- merge(meta,
                     tcrs.split_ab.merged.best,
                     by.x="row.names",
                     by.y="cell",all.x = T)
rownames(meta.merged) <- meta.merged$Row.names
meta.merged <- meta.merged[rownames(meta),]

seu.filt@meta.data <- meta.merged

seu.filt@meta.data$source_organism[which(seu.filt@meta.data$source_organism=="")] <- NA

DimPlot(subset(seu.filt, 
               cells = rownames(seu.filt@meta.data[!is.na(seu.filt@meta.data$source_organism),])), 
        group.by = "source_organism")

seu.filt@meta.data$source_organism2 <- seu.filt@meta.data$source_organism

#seu.filt@meta.data$source_organism2[seu.filt@meta.data$source_organism2!="Human"] <- "Non_human_or_ambiguous"
seu.filt@meta.data$source_organism2[grep("Human",seu.filt@meta.data$source_organism2)] <- "Human_match"
seu.filt@meta.data$source_organism2[!grepl("Human",seu.filt@meta.data$source_organism2) & 
                                    !is.na(seu.filt@meta.data$source_organism2)] <- "Non_human_match"

DimPlot(subset(seu.filt, 
               cells = rownames(seu.filt@meta.data[!is.na(seu.filt@meta.data$source_organism2),])), 
        group.by = "source_organism2",order = T)

sample_info <- as.data.frame(read_excel("~/proj/um_ss/Investigations/data/SCANDIUM WGS 2022.xlsx"),
                             stringsAsFactors=F)
sample_info <- sample_info[!is.na(sample_info$BOR),]

sample_info$`SCC-PDX kod`
sample_info$subject...2
sample_info$subject...9

seu.filt@meta.data$UM.ID
seu.filt@meta.data

all_um_sample_info <- as.data.frame(read_excel("~/proj/common_data/um.sample_master_sheet.v5.xlsx"),
                                    stringsAsFactors=F)

all_um_sample_info <- all_um_sample_info[,c("scandium.PATIENT", "other.PATIENT", "UM ID")]

all_um_sample_info[all_um_sample_info$`UM ID` %in% seu.filt@meta.data$UM.ID,]

which(! (seu.filt@meta.data$UM.ID %in% all_um_sample_info$`UM ID`))

sample_info$`SCC-PDX kod`

organisms <- unique(unlist(strsplit(unique(seu.filt@meta.data$source_organism),"[;|]")))
organisms <- organisms[!is.na(organisms)]

clusters <- unique(as.character(seu.filt.cd8_t@active.ident))

get_organism_cluster_cor <- function(seu.filt, seu.filt.cd8_t, organism, cluster){
  #organims <- "CMV"
  #cluster <- "Activated, FCGR3A+"
  tab.cmv <- table(grepl(pattern = organism, seu.filt@meta.data$source_organism, fixed = T), 
                   seu.filt@meta.data$UM.ID)
  tab.cmv <- as.data.frame.matrix(tab.cmv)
  tab.cmv <- t(tab.cmv)/colSums(tab.cmv)
  colnames(tab.cmv) <- c("organism_false_perc","organism_true_perc")
  
  #seu.filt.cd8_t <- readRDS(paste0(outdir,"seu.filt.cd8_t.rda"))
  tab.nk_like <- table(seu.filt.cd8_t@active.ident == cluster, seu.filt.cd8_t@meta.data$UM.ID)
  tab.nk_like <- as.data.frame.matrix(tab.nk_like)
  tab.nk_like <- t(tab.nk_like)/colSums(tab.nk_like)
  colnames(tab.nk_like) <- c("cluster_false_perc","cluster_true_perc")
  
  tab <- merge(tab.cmv,tab.nk_like,by="row.names")
  colnames(tab)[colnames(tab)=="Row.names"] <- "UM.ID"
  
  return(tab)
}

tab.list <- list()
for (organism in organisms){
  for (cluster in clusters){
    tab.list[[paste0(organism,":",cluster)]] <- get_organism_cluster_cor(
      seu.filt, seu.filt.cd8_t, organism, cluster)
  }
}

library(ggplot2)
ggplot(tab,aes(x=cmv_true_perc,y=nk_link_true_perc)) + geom_point()



tab.cor <- lapply(tab.list, function(x){
  cor.test(x$organism_true_perc,x$cluster_true_perc)
})

tab.cor.flat <- lapply(tab.cor,function(x){
  data.frame(cor=x$estimate,p=x$p.value)
})
tab.cor.flat <- do.call("rbind",tab.cor.flat)
tab.cor.flat$organism <- str_split_fixed(rownames(tab.cor.flat),":",2)[,1]
tab.cor.flat$cluster <- str_split_fixed(rownames(tab.cor.flat),":",2)[,2]
rownames(tab.cor.flat) <- NULL

tab.cor.mat <- (tab.cor.flat)

tab.cor.mat <- reshape(tab.cor.flat[,c("cor","organism","cluster")], 
        idvar = "organism", timevar = "cluster", direction = "wide")
rownames(tab.cor.mat) <- tab.cor.mat$organism
tab.cor.mat$organism <- NULL
colnames(tab.cor.mat) <- str_split_fixed(colnames(tab.cor.mat),"\\.",2)[,2]

pheatmap(tab.cor.mat, border_color = NA)


meta <- seu.filt.cd8_t@meta.data
meta.merged <- merge(meta,
                     tcrs.split_ab.merged.best,
                     by.x="row.names",
                     by.y="cell",all.x = T)
rownames(meta.merged) <- meta.merged$Row.names
meta.merged <- meta.merged[rownames(meta),]

seu.filt.cd8_t@meta.data <- meta.merged

get_organism_cluster_perc <- function(seu.filt, seu.filt.cd8_t, organism){
  #organims <- "CMV"
  #cluster <- "Activated, FCGR3A+"
  tab.cmv <- table(grepl(pattern = organism, seu.filt.cd8_t@meta.data$source_organism, fixed = T), 
                   seu.filt.cd8_t@active.ident)
  tab.cmv <- as.data.frame.matrix(tab.cmv)
  stopifnot(identical(names(tab.cmv),names(table(seu.filt.cd8_t@active.ident))))
  
  tab.cmv <- t(tab.cmv)/colSums(tab.cmv)
  colnames(tab.cmv) <- c("organism_false_perc","organism_true_perc")
  
  tab.cmv
  
  return(tab.cmv)
}

organisms <- unique(unlist(strsplit(unique(
  seu.filt.cd8_t@meta.data$source_organism),"[;|]")))
organisms <- organisms[!is.na(organisms)]

clusters <- unique(as.character(seu.filt.cd8_t@active.ident))

tab_organism_cluste_perc.list <- list()
for (organism in organisms){
    tab_organism_cluste_perc.list[[organism]] <- get_organism_cluster_perc(
      seu.filt, seu.filt.cd8_t, organism)
}

tab_organism_cluste_perc.flat <- do.call('rbind', lapply(
  tab_organism_cluste_perc.list, function(x){x[,2]}))

p <- pheatmap(tab_organism_cluste_perc.flat,border_color = NA)
tab_organism_cluste_perc.flat.na_zero <- tab_organism_cluste_perc.flat
tab_organism_cluste_perc.flat.na_zero[tab_organism_cluste_perc.flat.na_zero==0] <- NA
pdf(file = paste0(outdir,"tcrmatch_all_cd8_t_phenotypes_perc_matches.pdf"), width = 4, height = 4.5)
pheatmap(tab_organism_cluste_perc.flat.na_zero[p$tree_row$order,p$tree_col$order],
         border_color = NA,cluster_rows = F,cluster_cols = F,na_col = NA)
dev.off()

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

# ggsankey ---------------------------------------------------------------------
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

seu.filt.cd8_t <- readRDS(paste0(outdir,"seu.filt.cd8_t.rda"))

df <- data.frame(seu.filt.cd8_t@active.ident)
colnames(df) <- "celltype"
tcrs.split_ab.merged <- merge(tcrs.split_ab.merged,df,by.x="cell",by.y="row.names")

seu.filt.cd8_t <- add_anatimic_location_meta(
  seu.filt.cd8_t, pth = "~/proj/um_ss/Investigations/samples_10x.xlsx")

orgamism_rename_map <- get_organism_rename_map(tcrmatch_output)

tcrs.split_ab.merged <- merge(tcrs.split_ab.merged,
      data.frame(cell = rownames(seu.filt.cd8_t@meta.data),
                 tissue_site = seu.filt.cd8_t@meta.data$Original_biopsy_tissue_site), 
      by = "cell")

lst.biopsy <- tcrs.split_ab.merged[,c("cell", "celltype", "tissue_site", "antigen", "organism")]
#lst.biopsy <- unique(lst.biopsy)

colnames(lst.biopsy) <- c("cell", "celltype", "tissue_site", "Antigen", "Species")
lst.biopsy$Antigen[lst.biopsy$Antigen==""] <- "Unknown protein"
lst.biopsy$Species[lst.biopsy$Species==""] <- "Unknown species"
lst.biopsy$celltype <- as.character(lst.biopsy$celltype)

lst.biopsy <- lst.biopsy %>% separate_longer_delim(c(Antigen, Species), delim = ";")
#lst.biopsy <- unique(lst.biopsy)

for (i in 1:nrow(lst.biopsy)){
  if (lst.biopsy[i,]$Species != "Unknown species"){
    lst.biopsy[i,]$Species <- orgamism_rename_map$new[
      orgamism_rename_map$old==lst.biopsy[i,]$Species]
  }
}

#lst <- unique(lst.biopsy[,c("celltype", "tissue_site", "Antigen", "Species")])
lst <- lst.biopsy[,c("cell", "celltype", "tissue_site", "Antigen", "Species")]
#lst <- unique(lst)
lst <- lst[lst$Species %in% names(which(table(lst$Species)>10)),]
lst <- unique(lst)
df <- lst %>% make_long(celltype, tissue_site, Species)

celltype <- names(sort(table(lst$celltype),decreasing = F))
tissue_site <- names(sort(table(lst$tissue_site),decreasing = F))
Antigen <- unique(sort(lst$Antigen,decreasing = T))
Species <- unique(sort(lst$Species,decreasing = T))
fact_order <- c(celltype, tissue_site, Antigen, Species)
df$node <- factor(df$node, levels = fact_order)
df$next_node <- factor(df$next_node, levels = fact_order)

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
               fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = NA) +
  geom_sankey_label(size = 3, color = "black", fill = "white") +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))

