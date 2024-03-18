library(readxl)
library(WriteXLS)
library(stringr)
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)
library(reshape2)
library(ggsankey)
library(dplyr)

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/"

# Write input to TCRmatch ------------------------------------------------------
gliph <- read.table(paste0(outdir,"seu.filt.cd8_t_hla_ref_v2.csv"),sep=",",header = T)
gliph <- gliph[!is.na(gliph$index),]
stopifnot(length(grep(pattern="^C",gliph$TcRb,invert = T))==0)
stopifnot(length(grep(pattern="F$",gliph$TcRb,invert = T))==0)
gliph$tcrmatch_input <- gliph$TcRb
gliph$tcrmatch_input <- gsub(pattern="^C",replacement="",gliph$tcrmatch_input)
gliph$tcrmatch_input <- gsub(pattern="F$",replacement="",gliph$tcrmatch_input)
write.table(gliph$tcrmatch_input,file = paste0(outdir,"gliph.cd8_t_hla_ref_v2.tcrmatch_input.txt"),
            sep = "\t",col.names = F,quote = F,row.names = F)

# Read output from TCRmatch and merge with gliph table -------------------------

pth <- "~/nimbus/data/proj/um_ss/Investigations/seurat/results/v16/gliph.cd8_t_hla_ref_v2.tcrmatch_output_2.txt"

read_tcrmatch <- function(pth){
  tcrmatch_output <- read.table(pth,sep = "\t",header = T)
  tcrmatch_output$original_tcrmatch_row <- rownames(tcrmatch_output)
  
  tcrmatch_output$source_organism <- gsub(pattern=", $",replacement = "",tcrmatch_output$source_organism)
  tcrmatch_output$antigen <- gsub(pattern=", $",replacement = "",tcrmatch_output$antigen)
  tcrmatch_output$source_organism[tcrmatch_output$source_organism==" "] <- ""
  tcrmatch_output$antigen[tcrmatch_output$antigen==" "] <- ""
  tcrmatch_output$source_organism <- gsub(pattern="^,",replacement = "",tcrmatch_output$source_organism)
  tcrmatch_output$antigen <- gsub(pattern="^,",replacement = "",tcrmatch_output$antigen)
  tcrmatch_output$source_organism <- gsub(pattern="^ ,",replacement = "",tcrmatch_output$source_organism)
  tcrmatch_output$antigen <- gsub(pattern="^ ,",replacement = "",tcrmatch_output$antigen)
  
  return(tcrmatch_output)
}

tcrmatch_output <- read_tcrmatch(pth)

gliph$original_gliph_row <- rownames(gliph)

gliph.merged <- merge(gliph,tcrmatch_output,by.x="tcrmatch_input",by.y="input_sequence",all.x=T)
gliph.merged <- gliph.merged[order(gliph.merged$index,decreasing = F),]
gliph.merged$epitopes[is.na(gliph.merged$epitopes)] <- ""
gliph.merged$antigen[is.na(gliph.merged$antigen)] <- ""
gliph.merged$source_organism[is.na(gliph.merged$source_organism)] <- ""

WriteXLS(gliph.merged, ExcelFileName = paste0(outdir,"gliph.merged.xlsx"),
         BoldHeaderRow = T, AutoFilter = T, AdjWidth = T,row.names = F)

nms <- unique(gliph.merged$original_gliph_row)

gliph.merged$epitope_antigen <- ""

# Simplify species
sort(unique(unlist(strsplit(unique(gliph.merged$source_organism),","))))

gliph.merged$source_organism <- gsub(
  pattern = "Influenza A virus (A/Puerto Rico/8/1934(H1N1)) (Influenza A virus (A/PR/8/1934(H1N1)))",
  replacement = "Influenza A virus", gliph.merged$source_organism,fixed = T)

gliph.merged$source_organism <- gsub(
  pattern = "Human herpesvirus 5 strain AD169 (Human cytomegalovirus (strain AD169))",
  replacement = "Human herpesvirus 5 (Human cytomegalovirus)", gliph.merged$source_organism,fixed = T)

gliph.merged$source_organism <- gsub(
  pattern = "Human herpesvirus 5 strain Towne (Human cytomegalovirus (strain Towne))",
  replacement = "Human herpesvirus 5 (Human cytomegalovirus)", gliph.merged$source_organism,fixed = T)

sort(unique(unlist(strsplit(unique(gliph.merged$source_organism),","))))

out <- list()
for (nm in nms){
  tmp <- gliph.merged[gliph.merged$original_gliph_row==nm,]
  for (i in 1:nrow(tmp)){
    tmp_row <- paste0(unlist(strsplit(tmp$epitopes[i],",")),":",
                      gsub(pattern="Nuclear antigen EBNA-3",replacement="Epstein-Barr nuclear antigen 3",
                           as.character(sapply(
                        gsub(pattern="^ ",replacement="",unlist(strsplit(tmp$antigen[i],","))),
                        function(x){
                          paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
                        }))),":",
                      unlist(strsplit(tmp$source_organism[i],","))
                      )
    tmp_row <- unique(tmp_row)
    tmp_row <- paste0(tmp_row,collapse = ";")
    tmp$epitope_antigen[i] <- tmp_row
  }
  out[[nm]] <- tmp
}

get_lst <- function(out){
  lst <- list()
  for (nm in names(out)){
    tmp <- unique(out[[nm]]$epitope_antigen)
    tmp <- unique(unlist(strsplit(tmp,";")))
    cell <- nm
    cluster <- unique(out[[nm]]$index)
    peptides <- str_split_fixed(tmp,":",3)[,1]
    antigens <- str_split_fixed(tmp,":",3)[,2]
    species <- str_split_fixed(tmp,":",3)[,3]
    rw <- data.frame(peptides=peptides,antigens=antigens,species=species,
                     cluster=cluster,cell=cell)
    lst[[nm]] <- rw
  }
  
  lst <- do.call("rbind",lst)
  rownames(lst) <- NULL
  
  lst$cluster <- paste0("C",lst$cluster)
  
  lst$peptides[lst$peptides==""] <- "Unknown epitope"
  lst$antigens[lst$antigens==""] <- "Unknown protein"
  lst$species[lst$species==""] <- "Unknown species"
  
  return(lst)
}

lst <- get_lst(out)

high_confidence <- gliph$vb_score < 0.05 & 
  gliph$index %in% names(which(table(gliph$index)>2)) & 
  gliph$index %in% names(which(rowSums(table(gliph$index, gliph$Sample)>0)>2))

gliph.high_confidence <- gliph[high_confidence,]

nodes_high_confidence <- paste0("C",unique(gliph.high_confidence$index))

lst <- lst[lst$cluster %in% nodes_high_confidence,]

lst$antigens <- gsub(pattern=" \\[.+\\]", replacement="", lst$antigens)

# ggsankey ---------------------------------------------------------------------
colnames(lst)[colnames(lst)=="cluster"] <- "Cluster"
colnames(lst)[colnames(lst)=="cell"] <- "CDR3b"
colnames(lst)[colnames(lst)=="antigens"] <- "Antigen"
colnames(lst)[colnames(lst)=="peptides"] <- "Epitope"
colnames(lst)[colnames(lst)=="species"] <- "Species"
lst$CDR3b <- paste0("CDR3_",lst$CDR3b)
df <- lst %>% make_long(Cluster, CDR3b, Epitope, Antigen, Species)
Clusters <- rev(unique(df$node[df$node %in% lst$Cluster]))
CDR3bs <- rev(unique(df$node[df$node %in% lst$CDR3b]))
Epitope <- rev(unique(df$node[df$node %in% lst$Epitope]))
Antigen <- rev(unique(df$node[df$node %in% lst$Antigen]))
Species <- rev(unique(df$node[df$node %in% lst$Species]))
fact_order <- c(Clusters, CDR3bs, Epitope, Antigen, Species)
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
ggsave(filename = paste0(outdir,"gliph.highconfidence.sankeyNetwork.pdf"),
       width = 14, height = 9.5)

# Map gliph clusters and antigen matches to tSNE -------------------------------
library(Seurat)
library(gridExtra)

s.t <- readRDS(file=paste0(outdir,"seu.filt.cd8_t.rda"))

gliph.high_confidence$clonotype_id <- paste0(gliph.high_confidence$TcRa,":",gliph.high_confidence$TcRb,":",
       gliph.high_confidence$V,":",gliph.high_confidence$J)

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

tcrs.split_ab <- get_tcrs_split_ab_non_unique(s.t)

head(tcrs.split_ab)

tcrs.split_ab$clonotype_id <- paste0(tcrs.split_ab$CDR3a,":",tcrs.split_ab$CDR3b,
                                     ":",tcrs.split_ab$TRBV,":",tcrs.split_ab$TRBJ)

gliph.high_confidence$cluster <- paste0("C",gliph.high_confidence$index)

# Plot by cluster
unique_clusters <- unique(gliph.high_confidence$cluster)
clusters.cells <- list()
for (clus in unique_clusters){
  clusters.cells[[clus]] <- rownames(tcrs.split_ab)[
    tcrs.split_ab$clonotype_id %in% gliph.high_confidence$clonotype_id[
      gliph.high_confidence$cluster==clus]]
}

p <- list()
for (clus in unique_clusters){
  p[[clus]] <- DimPlot(s.t,cells.highlight = clusters.cells[[clus]]) + 
    theme(legend.position="none") + ggtitle(clus)
}
pdf(file = paste0(outdir,"tsne_gliph_by_cluster.pdf"),width = 20,height = 10)
grid.arrange(grobs = p, ncol=4)
dev.off()

# Plot by clusters with information about species matches
links_cluster_species <- lst[,c("Cluster","Species")]
x <- table(links_cluster_species$Cluster,links_cluster_species$Species)
x <- as.data.frame.matrix(x)
species_nrs <- list()
for (clus in unique_clusters){
  nms <- colnames(x[clus,])[which(x[clus,]>0)]
  nrs <- as.numeric(x[clus,][which(x[clus,]>0)])
  species_nrs[[clus]] <- paste0(sort(paste0(nms,":",nrs)),collapse=";")
}

p <- list()
for (clus in unique_clusters){
  p[[clus]] <- DimPlot(s.t,cells.highlight = clusters.cells[[clus]]) + 
    theme(legend.position="none") + ggtitle(paste0(clus,":",species_nrs[[clus]])) + 
    theme(plot.title = element_text(size = 8))
}
pdf(file = paste0(outdir,"tsne_gliph_by_cluster.species.pdf"),width = 20,height = 10)
grid.arrange(grobs = p, ncol=4)
dev.off()

# Plot clusters with info on antigen matches
links_cluster_antigen <- lst[,c("Cluster","Antigen")]
x <- table(links_cluster_antigen$Cluster,links_cluster_antigen$Antigen)
x <- as.data.frame.matrix(x)
antigen_nrs <- list()
for (clus in unique_clusters){
  nms <- colnames(x[clus,])[which(x[clus,]>0)]
  nms <- gsub(pattern=" \\[.+\\]",replacement="",nms)
  nrs <- as.numeric(x[clus,][which(x[clus,]>0)])
  antigen_nrs[[clus]] <- paste0(sort(paste0(nms,":",nrs)),collapse=";")
}

p <- list()
for (clus in unique_clusters){
  p[[clus]] <- DimPlot(s.t,cells.highlight = clusters.cells[[clus]]) + 
    theme(legend.position="none") + ggtitle(paste0(clus,":",antigen_nrs[[clus]])) + 
    theme(plot.title = element_text(size = 8))
}
pdf(file = paste0(outdir,"tsne_gliph_by_cluster.antigen.pdf"),width = 20,height = 10)
grid.arrange(grobs = p, ncol=4)
dev.off()

# Plot by individual antigens instead of clusters
tmp <- lst[lst$Cluster %in% unique(gliph.high_confidence$cluster),c("CDR3b","Antigen")]
tmp$Cluster <- ""
for (i in 1:nrow(tmp)){
  tmp$Cluster <- gliph.high_confidence$cluster[
    gliph.high_confidence$original_gliph_row==gsub(pattern="CDR3_",replacement="",tmp$CDR3b[i])]
}
tmp <- unique(tmp)
tmp$Antigen <- gsub(pattern=" \\[.+\\]", replacement="", tmp$Antigen)

unique_antigens <- unique(tmp$Antigen)
antigens.cells <- list()
for (antigen in unique_antigens){
  antigens.cells[[antigen]] <- rownames(tcrs.split_ab)[
    tcrs.split_ab$clonotype_id %in% gliph.high_confidence$clonotype_id[
      gliph.high_confidence$original_gliph_row %in% tmp$cell[tmp$Antigen==antigen]
    ]]
}

p <- list()
for (antigen in unique_antigens){
  p[[antigen]] <- DimPlot(s.t, cells.highlight = antigens.cells[[antigen]]) + 
    theme(legend.position="none") + ggtitle(antigen) + 
    theme(plot.title = element_text(size = 8))
}
pdf(file = paste0(outdir,"tsne_gliph_by_antigen.pdf"),width = 20,height = 10)
grid.arrange(grobs = p, ncol=5)
dev.off()

links_cluster_species <- lst[,c("Cluster","Species")]

links_cluster_CDR3b <- lst[,c("Cluster","CDR3b")]
table(links_cluster_CDR3b$Cluster,links_cluster_CDR3b$CDR3b)

stopifnot(max(colSums(table(links_cluster_CDR3b$Cluster,links_cluster_CDR3b$CDR3b)>0))==1)

gliph$clonotype_sample <- paste0(gliph$TcRb,":", gliph$V,":", 
                                 gliph$J,":", gliph$TcRa,":", 
                                 gliph$Sample)

links_cluster_CDR3b <- gliph[,c("index", "clonotype_sample")]

nms <- names(which(colSums(table(links_cluster_CDR3b$index, 
                                 links_cluster_CDR3b$clonotype_sample)>0)>1))

# List co-enriched alleles -----------------------------------------------------

unique_clusters <- paste0("C",unique(gliph$index))
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
  
  co_enriched_hla[[clus]] <- data.frame(HLA.A=x.a, HLA.B=x.b, HLA.C=x.c)
}

co_enriched_hla <- do.call("rbind",co_enriched_hla)
head(co_enriched_hla)

co_enriched_hla$is_high_confidence <- rownames(co_enriched_hla) %in% gliph.high_confidence$cluster

WriteXLS(co_enriched_hla, ExcelFileName = paste0(outdir,"co_enriched_hla.xlsx"),
         AutoFilter = T, BoldHeaderRow = T, row.names = T, col.names = T,
         AdjWidth = T)

# Plot all MART1 matching CDR3 -------------------------------------------------
lst <- get_lst(out)
lst$antigens <- gsub(pattern=" \\[.+\\]", replacement="", lst$antigens)

lst.mart1 <- lst[grep("Melanoma",lst$antigens),]

colnames(lst.mart1)[colnames(lst.mart1)=="cluster"] <- "Cluster"
colnames(lst.mart1)[colnames(lst.mart1)=="cell"] <- "CDR3b"
colnames(lst.mart1)[colnames(lst.mart1)=="antigens"] <- "Antigen"
colnames(lst.mart1)[colnames(lst.mart1)=="peptides"] <- "Epitope"
colnames(lst.mart1)[colnames(lst.mart1)=="species"] <- "Species"

gliph$clonotype_id <- paste0(gliph$TcRa, ":", gliph$TcRb, ":",
                             gliph$V, ":" ,gliph$J)
gliph$cluster <- paste0("C",gliph$index)

unique_clusters <- unique(gliph$cluster)
clusters.cells <- list()
for (clus in unique_clusters){
  clusters.cells[[clus]] <- rownames(tcrs.split_ab)[
    tcrs.split_ab$clonotype_id %in% gliph$clonotype_id[
      gliph$cluster==clus]]
}

clusters.cells <- clusters.cells[which(names(clusters.cells) %in% lst.mart1$Cluster)]
clusters.cells <- unique(unlist(clusters.cells))

pdf(file = paste0(outdir,"tsne_gliph_all_mart1_matches.pdf"), width = 10, height = 10)
DimPlot(s.t, cells.highlight = clusters.cells) + 
  theme(legend.position="none") + ggtitle("MART1 matches") + 
  theme(plot.title = element_text(size = 8))
dev.off()

table(s.t@meta.data[clusters.cells,]$UM.ID)

lst.mart1 <- lst[grep("Melanoma",lst$antigens),]
clusters.cells <- clusters.cells[which(names(clusters.cells) %in% lst.mart1$Cluster)]
clusters.cells <- unique(unlist(clusters.cells))

