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


#dups <- names(which(table(gliph.merged$original_gliph_row)>1))
nms <- unique(gliph.merged$original_gliph_row)

gliph.merged$epitope_antigen <- ""
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
out2 <- lapply(out,function(x){
  unq <- unique(x$epitope_antigen);
  tmp <- c()
  for (i in 1:length(unq)){
    tmp <- c(tmp,unique(unlist(strsplit(unq[i],";"))))
  }
  tmp <- unique(tmp)
  tmp <- table(str_split_fixed(tmp,":",2)[,2])
  tmp
  })
unique_antigens <- sort(unique(unlist(lapply(out2,names))))

m <- matrix(nrow=length(out2),ncol=length(unique_antigens),dimnames = list(names(out2),unique_antigens))
for (i in 1:length(out2)){
  for (j in 1:length(unique_antigens)){
    m[i,j] <- sum(out2[[i]][unique_antigens[j]])
  }
}

library(pheatmap)
pheatmap(m,cluster_rows = F,cluster_cols = F,show_colnames = F)

stopifnot(nrow(m)==nrow(gliph))

m <- as.data.frame(m)
nms <- unique(gliph.merged$original_gliph_row)

stopifnot(identical(nms,rownames(m)))
m$gliph_index <- gliph$index

unique_gliph_index <- unique(m$gliph_index)
n_matches_gliph_cluster <- matrix(nrow=length(unique_gliph_index),
                                  ncol=ncol(m)-1,
                                  dimnames = list(unique_gliph_index,
                                                  colnames(m)[1:(length(colnames(m))-1)]))
for (i in 1:length(unique_gliph_index)){
  gi <- unique_gliph_index[i]
  for (j in 1:(ncol(m)-1)){
    n_matches_gliph_cluster[i,j] <- sum(m[m$gliph_index == gi,j],na.rm = T)
  }
}

pheatmap(n_matches_gliph_cluster,cluster_rows = F,cluster_cols = F,show_colnames = F)

stopifnot(identical(unique(as.character(gliph$index)),rownames(n_matches_gliph_cluster)))

gliph_index_n <- setNames(as.data.frame(table(gliph$index),stringsAsFactors = F),c("gliph_index","n"))

dim(gliph_index_n)
dim(n_matches_gliph_cluster)

stopifnot(identical(as.character(gliph_index_n$gliph_index),rownames(n_matches_gliph_cluster)))

perc_matches_gliph_cluster <- n_matches_gliph_cluster/gliph_index_n$n

apply(perc_matches_gliph_cluster[,2:ncol(perc_matches_gliph_cluster)],1,function(x){
  colnames(perc_matches_gliph_cluster)[2:ncol(perc_matches_gliph_cluster)][which.max(x)]})

apply(perc_matches_gliph_cluster,1,function(x){
  x <- x[x>0]
  paste0(paste0(names(x),":",
         as.numeric(x)),collapse=";")
  
})

library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)
library(reshape2)


high_confidence <- gliph$vb_score < 0.05 & 
  gliph$index %in% names(which(table(gliph$index)>2)) & 
  gliph$index %in% names(which(rowSums(table(gliph$index, gliph$Sample)>0)>2))

table(high_confidence)
gliph.high_confidence <- gliph[high_confidence,]

nodes_high_confidence <- paste0("C",unique(gliph.high_confidence$index))

nodes_keep <- paste0("C",gliph_index_n[gliph_index_n$n > 2,]$gliph_index)

nodes_keep <- intersect(nodes_high_confidence,nodes_keep)

target_nodes <- data.frame(name=colnames(perc_matches_gliph_cluster))
source_nodes <- data.frame(name=paste0("C",rownames(perc_matches_gliph_cluster)))
nodes <- rbind(target_nodes,source_nodes)

links <- melt(perc_matches_gliph_cluster)
colnames(links) <- c("source","target","value")
links$source <- paste0("C",links$source)
links$target <- as.character(links$target)

links <- links[links$source %in% nodes_keep,]
links <- links[links$value > 0,]
links <- links[links$target!=":",]

nodes <- nodes[nodes$name %in% links$source | nodes$name %in% links$target,]
nodes <- data.frame(name=nodes)


target_nodes_species <- data.frame(name=str_split_fixed(colnames(perc_matches_gliph_cluster),":",2)[,2])
source_nodes <- data.frame(name=paste0("C",rownames(perc_matches_gliph_cluster)))
nodes_species <- rbind(target_nodes_species,source_nodes)

#stopifnot(identical(as.character(unique(links$source)),rownames(perc_matches_gliph_cluster)))
for (i in 1:nrow(links)){
  links$source[i] <- as.numeric(rownames(nodes)[nodes$name==links$source[i]])
  links$target[i] <- as.numeric(rownames(nodes)[nodes$name==links$target[i]])
}

links <- links[links$value > 0,]
links$value <- links$value*100
links$source <- as.numeric(links$source)
links$target <- as.numeric(links$target)

links$source <- links$source - 1
links$target <- links$target - 1

sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name")

# n_matches_gliph_cluster_species <- as.data.frame(n_matches_gliph_cluster)
# species <- str_split_fixed(colnames(n_matches_gliph_cluster_species),":",2)[,2]
# unique_species <- unique(species)
# 
# m_species <- matrix(nrow=nrow(n_matches_gliph_cluster_species),
#                     ncol=length(unique_species),
#                     dimnames=list(rownames(n_matches_gliph_cluster_species),
#                                   unique_species))
# nms <- unique(gliph.merged$index)
# n_matches_gliph_cluster_species$gliph_index <- nms
# 
# stopifnot(identical(str_split_fixed(colnames(n_matches_gliph_cluster_species),":",2)[,2][
#   1:(ncol(n_matches_gliph_cluster_species)-1)],species))
# for (i in 1:length(nms)){
#   gi <- nms[i]
#   for (j in 1:length(unique_species)){
#     sp <- unique_species[j]
#     m_species[i,j] <- sum(n_matches_gliph_cluster_species[
#       n_matches_gliph_cluster_species$gliph_index == gi,species==sp],na.rm = T)
#   }
# }
# 
# perc_matches_gliph_cluster_species <- m_species/gliph_index_n$n

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

#sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
#              Target = "target", Value = "value", NodeID = "name")

high_confidence <- gliph$vb_score < 0.05 & 
  gliph$index %in% names(which(table(gliph$index)>2)) & 
  gliph$index %in% names(which(rowSums(table(gliph$index, gliph$Sample)>0)>2))

gliph.high_confidence <- gliph[high_confidence,]

nodes_high_confidence <- paste0("C",unique(gliph.high_confidence$index))

lst <- lst[lst$cluster %in% nodes_high_confidence,]

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

#links$value <- 1

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

#links <- links[links$value > 0,]
#links$value <- links$value*100

#nodes_keep <- paste0("C",gliph_index_n[gliph_index_n$n > 2,]$gliph_index)
#nodes_keep <- intersect(nodes_high_confidence,nodes_keep)
#links$source

sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name")
