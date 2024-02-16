library(readxl)
library(WriteXLS)
library(stringr)

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
tcrmatch_output <- read.table(pth,sep = "\t",header = T)
head(tcrmatch_output)
tcrmatch_output$original_tcrmatch_row <- rownames(tcrmatch_output)

tcrmatch_output$source_organism <- gsub(pattern=", $",replacement = "",tcrmatch_output$source_organism)
tcrmatch_output$antigen <- gsub(pattern=", $",replacement = "",tcrmatch_output$antigen)
tcrmatch_output$source_organism[tcrmatch_output$source_organism==" "] <- ""
tcrmatch_output$antigen[tcrmatch_output$antigen==" "] <- ""
tcrmatch_output$source_organism <- gsub(pattern="^,",replacement = "",tcrmatch_output$source_organism)
tcrmatch_output$antigen <- gsub(pattern="^,",replacement = "",tcrmatch_output$antigen)
tcrmatch_output$source_organism <- gsub(pattern="^ ,",replacement = "",tcrmatch_output$source_organism)
tcrmatch_output$antigen <- gsub(pattern="^ ,",replacement = "",tcrmatch_output$antigen)

unique(tcrmatch_output$antigen)

gliph$original_gliph_row <- rownames(gliph)

gliph.merged <- merge(gliph,tcrmatch_output,by.x="tcrmatch_input",by.y="input_sequence",all.x=T)
dim(gliph)
dim(gliph.merged)
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

library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)
library(reshape2)

lst <- list()
for (nm in names(out)){
  tmp <- unique(out[[nm]]$epitope_antigen)
  tmp <- unique(unlist(strsplit(tmp,";")))
  cell <- nm
  cluster <- unique(out[[nm]]$index)
  peptides <- str_split_fixed(tmp,":",3)[,1]
  antigens <- str_split_fixed(tmp,":",3)[,2]
  species <- str_split_fixed(tmp,":",3)[,3]
  rw <- data.frame(peptides=peptides,antigens=antigens,species=species,cluster=cluster,cell=cell)
  lst[[nm]] <- rw
}

lst <- do.call("rbind",lst)
rownames(lst) <- NULL

lst$cluster <- paste0("C",lst$cluster)

lst$peptides[lst$peptides==""] <- "Unknown epitope"
lst$antigens[lst$antigens==""] <- "Unknown protein"
lst$species[lst$species==""] <- "Unknown species"

high_confidence <- gliph$vb_score < 0.05 & 
  gliph$index %in% names(which(table(gliph$index)>2)) & 
  gliph$index %in% names(which(rowSums(table(gliph$index, gliph$Sample)>0)>2))

gliph.high_confidence <- gliph[high_confidence,]

nodes_high_confidence <- paste0("C",unique(gliph.high_confidence$index))

lst <- lst[lst$cluster %in% nodes_high_confidence,]

lst$antigens <- gsub(pattern=" \\[.+\\]", replacement="", lst$antigens)

# Cluster to CDR3 (cell) links
links_cluster_cdr3 <- unique(lst[,c("cluster","cell")])
links_cluster_cdr3$cell <- paste0("CDR3_",links_cluster_cdr3$cell)
colnames(links_cluster_cdr3) <- c("source","target")
links_cluster_cdr3$value <- 1

# CDR3 to peptide links
links_cdr3_peptide <- unique(lst[,c("cell","peptides")])
links_cdr3_peptide$cell <- paste0("CDR3_",links_cdr3_peptide$cell)
colnames(links_cdr3_peptide) <- c("source","target")
links_cdr3_peptide$value <- 1

# Peptide to antigen (protein) links
links_peptide_antigen <- unique(lst[,c("peptides","antigens")])
colnames(links_peptide_antigen) <- c("source","target")
tab <- table(links_cdr3_peptide$source,links_cdr3_peptide$target)
tab <- as.data.frame.matrix(tab)
tab <- colSums(tab)
links_peptide_antigen$value <- NA
for (peptide in names(tab)){
  links_peptide_antigen$value[links_peptide_antigen$source==peptide] <- tab[peptide]
}

# Antigen to species links
links_antigen_species <- unique(lst[,c("antigens","species")])
colnames(links_antigen_species) <- c("source","target")
links_antigen_species$value <- NA
tab <- aggregate(links_peptide_antigen$value, by=list(antigen=links_peptide_antigen$target), FUN=sum)
colnames(tab)[colnames(tab)=="x"] <- "value"
links_antigen_species$value <- NA
for (antigen in tab$antigen){
  links_antigen_species$value[links_antigen_species$source==antigen] <- tab$value[tab$antigen==antigen]
}

links <- rbind(links_cluster_cdr3, 
               links_cdr3_peptide,
               links_peptide_antigen,
               links_antigen_species)

nodes <- data.frame(name=unique(c(links$source,links$target)))
head(nodes)

for (i in 1:nrow(links)){
  links$source[i] <- as.numeric(rownames(nodes)[nodes$name==links$source[i]])
  links$target[i] <- as.numeric(rownames(nodes)[nodes$name==links$target[i]])
}
links$source <- as.numeric(links$source)
links$target <- as.numeric(links$target)

links$source <- links$source - 1
links$target <- links$target - 1

#gliph.highconfidence.sankeyNetwork.pdf
g <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name", 
              nodePadding = 5)

library(htmlwidgets)
saveWidget(g, 
           file="/Users/00105606/proj/um_ss/Investigations/seurat/results/v16/gliph.highconfidence.sankeyNetwork.html",
           selfcontained = T)
devtools::install_github("davidsjoberg/ggsankey")# Open in web browser and save to PDF


library(webshot)
webshot("file:////Users/00105606/proj/um_ss/Investigations/seurat/results/v16/gliph.highconfidence.sankeyNetwork.html",
        "gliph.highconfidence.sankeyNetwork.R.pdf")

# Alternative with ggsankey
library(ggsankey)
library(ggplot2)
library(dplyr)

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
ggsave(filename = "gliph.highconfidence.sankeyNetwork.pdf")

# Compare HLA representation in each cluster, relative to background in participating samples

as.data.frame.matrix(table(gliph$HLA.A,gliph$index))
as.numeric(table(gliph$index))

(as.data.frame.matrix(table(gliph$HLA.A,gliph$index))/as.numeric(table(gliph$index)))[1,]

unique_hla_a <- unique(gsub(pattern="!",replacement="",unlist(strsplit(gliph$HLA.A,"/"))))
unique_clusters <- paste0("C",unique(gliph$index))

m_cluster_hla_a <- matrix(0,nrow=length(unique_clusters),
                          ncol=length(unique_hla_a),
                          dimnames=list(unique_clusters,unique_hla_a))

for (clus in unique_clusters){
  tmp <- gliph[paste0("C",gliph$index) == clus,]
  for (j in 1:length(tmp$HLA.A)){
    hlas <- gsub(pattern="!",replacement = "",unlist(strsplit(tmp$HLA.A[j],"/")))
    hlas <- unique(hlas) # Count by binary presence only (don't double count homozygous)
    for (i in 1:length(hlas)){
      m_cluster_hla_a[clus,hlas[i]] <- m_cluster_hla_a[clus,hlas[i]] + 1
    }
  }
}

m_cluster_hla_a <- as.data.frame(m_cluster_hla_a)
m_cluster_hla_a$n_members <- NA

for (clus in unique_clusters){
  m_cluster_hla_a$n_members[rownames(m_cluster_hla_a)==clus] <- sum(paste0("C",gliph$index)==clus)
}
stopifnot(all(apply(m_cluster_hla_a[,1:(ncol(m_cluster_hla_a)-1)],1,max) <= m_cluster_hla_a$n_members))

library(ggplot2)
library(reshape2)

m_cluster_hla_a$cluster <- rownames(m_cluster_hla_a)
m_cluster_hla_a.long <- melt(m_cluster_hla_a,id.vars = c("cluster","n_members"))
m_cluster_hla_a.long$percentage <- 100*m_cluster_hla_a.long$value/m_cluster_hla_a.long$n_members

ggplot(m_cluster_hla_a.long[m_cluster_hla_a.long$cluster %in% nodes_high_confidence,],
       aes(x=reorder(cluster,-n_members),y=percentage,fill=variable)) + 
  geom_bar(stat="identity",position="dodge") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab(NULL)

pheatmap(m_cluster_hla_a[,1:(ncol(m_cluster_hla_a)-2)]/m_cluster_hla_a$n_members)

# Map gliph clusters and antigen matches to tSNE -------------------------------
library(Seurat)
library(gridExtra)

s.t <- readRDS(file=paste0(outdir,"seu.filt.cd8_t.rda"))

gliph.high_confidence


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
links_cluster_species <- lst[,c("cluster","species")]
x <- table(links_cluster_species$cluster,links_cluster_species$species)
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
links_cluster_antigen <- lst[,c("cluster","antigens")]
x <- table(links_cluster_antigen$cluster,links_cluster_antigen$antigens)
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
tmp <- lst[lst$cluster %in% unique(gliph.high_confidence$cluster),c("cell","antigens")]
for (i in 1:nrow(tmp)){
  tmp$cluster <- gliph.high_confidence$cluster[gliph.high_confidence$original_gliph_row==tmp$cell[i]]
}
tmp <- unique(tmp)
tmp$antigens <- gsub(pattern=" \\[.+\\]", replacement="", tmp$antigens)

unique_antigens <- unique(tmp$antigens)
antigens.cells <- list()
for (antigen in unique_antigens){
  antigens.cells[[antigen]] <- rownames(tcrs.split_ab)[
    tcrs.split_ab$clonotype_id %in% gliph.high_confidence$clonotype_id[
      gliph.high_confidence$original_gliph_row %in% tmp$cell[tmp$antigens==antigen]
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

