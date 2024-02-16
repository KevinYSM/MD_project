library(Seurat)
library(tidyverse)
library(gridExtra)
library(readxl)
library(foreach)
library(doParallel)
library(pheatmap)

# Read data --------------------------------------------------------------------
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

  return(vdj.takara)
}

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/"
dat.tumor <- readRDS(file=paste0(outdir,"seu.filt.rda"))
dat.tumor@meta.data$Sample.ID <- dat.tumor@meta.data$orig.ident
dat.tumor@meta.data$project.name <- "tumor"

s.t <- readRDS(file=paste0(outdir,"seu.filt.cd8_t.rda"))
s.t@meta.data$Sample.ID <- s.t@meta.data$orig.ident

dat.til <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
dat.til <- subset(dat.til, subset = project.name %in% c("G18-023"))

tils.unintegrated <- readRDS(file = paste0(
  "~/proj/um_ss/Investigations/seurat/results/tils/", "tils.rda"))

names(tils.unintegrated) <- unlist(lapply(tils.unintegrated, function(x){
  unique(x@meta.data$UM.ID)
}))

tils.unintegrated <- lapply(tils.unintegrated, function(x){
  if (any(Idents(x) == "Gamma-delta T cells (TRGV9, TRDV2)")){
    x <- RenameIdents(x, `Gamma-delta T cells (TRGV9, TRDV2)` = "CD8 T cells")
  }
  x <- subset(x, cells = names(x@active.ident[x@active.ident != "CD8 T / NK cells"]))
})

vdj.takara <- read_vdj.takara()

# Map TIL TCRs to biopsies -----------------------------------------------------
# Question: do all yTIL TCRs end up in the exhausted biopsy clusters?

# 1: Determing overlapping TCRs between samples
#dat.tumor@meta.data$cdr3
#dat.til@meta.data$cdr3

extract_tcr <- function(dat.cca){
  cnames <- c("Sample.ID","project.name","UM.ID",
              "v_gene","d_gene","j_gene","c_gene","chains","cdr3",
              "reads","umis","productive")

  vdj.10x <- dat.cca@meta.data[,cnames]
  vdj.10x$cell_id <- rownames(vdj.10x)

  clean_fun <- function(x){
    res <- unlist(strsplit(x,";"))
    res <- res[!is.na(res)]
    res <- setdiff(res,"None")
    if (length(res)==0){
      res <- "None"
    }
    return(res)
  }

  vdj.10x_cdr3 <- setNames(sapply(vdj.10x$cdr3, clean_fun), rownames(vdj.10x))
  vdj.10x_v_gene <- setNames(sapply(vdj.10x$v_gene, clean_fun), rownames(vdj.10x))
  vdj.10x_d_gene <- setNames(sapply(vdj.10x$d_gene, clean_fun), rownames(vdj.10x))
  vdj.10x_j_gene <- setNames(sapply(vdj.10x$j_gene, clean_fun), rownames(vdj.10x))
  vdj.10x_c_gene <- setNames(sapply(vdj.10x$c_gene, clean_fun), rownames(vdj.10x))
  vdj.10x_chains <- setNames(sapply(vdj.10x$chains, clean_fun), rownames(vdj.10x))

  return(list(
    vdj.10x_cdr3 = vdj.10x_cdr3,
    vdj.10x_v_gene = vdj.10x_v_gene,
    vdj.10x_d_gene = vdj.10x_d_gene,
    vdj.10x_j_gene = vdj.10x_j_gene,
    vdj.10x_c_gene = vdj.10x_c_gene,
    vdj.10x_chains = vdj.10x_chains
  ))
}

# Matching criteria: there must be at least one complete subset of v, d, j, c and cdr3
# between two cells in each dataset

get_matches <- function(dat.tumor.vdj, dat.til.vdj, sname){
  til_matches <- list()
  vdj.10x_cdr3.match <- list()
  vdj.10x_v_gene.match <- list()
  vdj.10x_d_gene.match <- list()
  vdj.10x_j_gene.match <- list()
  vdj.10x_c_gene.match <- list()
  for (i in 1:length(dat.tumor.vdj[[sname]]$vdj.10x_cdr3)){
    if (all(dat.tumor.vdj[[sname]]$vdj.10x_cdr3[[i]] != "None")){
      vdj.10x_cdr3.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_cdr3, function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_cdr3[[i]])
      }))
      vdj.10x_v_gene.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_v_gene, function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_v_gene[[i]])
      }))
      vdj.10x_d_gene.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_d_gene, function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_d_gene[[i]])
      }))
      vdj.10x_j_gene.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_j_gene, function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_j_gene[[i]])
      }))
      vdj.10x_c_gene.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_c_gene, function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_c_gene[[i]])
      }))

      #til_matches[[i]] <- which(vdj.10x_cdr3.match[[i]])
      til_matches[[i]] <- which(vdj.10x_cdr3.match[[i]] &
                                  vdj.10x_v_gene.match[[i]] &
                                  vdj.10x_d_gene.match[[i]] &
                                  vdj.10x_j_gene.match[[i]] &
                                  vdj.10x_c_gene.match[[i]])
    } else {
      til_matches[[i]] <- integer(0)
    }
  }
  return(list(
    til_matches = til_matches,
    vdj.10x_cdr3.match = vdj.10x_cdr3.match,
    vdj.10x_v_gene.match = vdj.10x_v_gene.match,
    vdj.10x_d_gene.match = vdj.10x_d_gene.match,
    vdj.10x_j_gene.match = vdj.10x_j_gene.match,
    vdj.10x_c_gene.match = vdj.10x_c_gene.match
  ))
}

dat.tumor.vdj <- lapply(SplitObject(dat.tumor, split.by = "UM.ID"), extract_tcr)
dat.til.vdj <- lapply(SplitObject(dat.til, split.by = "UM.ID"), extract_tcr)

snames <- intersect(names(dat.tumor.vdj),
                    names(dat.til.vdj))
res <- lapply(snames,function(sname){
  get_matches(dat.tumor.vdj, dat.til.vdj, sname)
  })
names(res) <- snames

dat.tumor@meta.data$tcr_matches_til <- F
for (sname in snames){
  dat.tumor@meta.data[dat.tumor@meta.data$UM.ID==sname,]$tcr_matches_til <-
    unlist(lapply(res[[sname]]$til_matches,length)) > 0
}

tcr_matches_til <- setNames(dat.tumor@meta.data$tcr_matches_til,
                                   rownames(dat.tumor@meta.data))

dat.tumor <- AddMetaData(dat.tumor, metadata = tcr_matches_til,
                         col.name = "tcr_matches_til")

s.t <- AddMetaData(s.t, metadata = tcr_matches_til, col.name = "tcr_matches_til")

s.t.common <- subset(s.t, UM.ID %in% intersect(unique(dat.tumor@meta.data$UM.ID),
                                              names(tils.unintegrated)))
s.t.common@meta.data$UM.ID <- factor(
  s.t.common@meta.data$UM.ID,
  levels= unique(s.t.common@meta.data$UM.ID)[order(
    as.numeric(gsub("UM","",unique(s.t.common@meta.data$UM.ID))),
    decreasing = F)])

# Fig. S9c
pdf(file = paste0(outdir,"tcr_matches_til.pdf"),width = 24,height = 5)
DimPlot(s.t.common, group.by = "tcr_matches_til", order = T, split.by = "UM.ID")
dev.off()

rm(s.t.common)
gc()

# Alternative attempt to match TCRs --------------------------------------------
get_tcr_matrix <- function(s.t){
  if (class(s.t)=="Seurat"){
    meta <- s.t@meta.data
  } else if (class(meta.tumor_til)=="data.frame"){
    meta <- s.t
  } else {
    stop("Error")
  }
  cells <- rownames(meta)

  cdr3 <- lapply(meta$cdr3, function(x){ unlist(strsplit(x, ";")) })
  v_gene <- lapply(meta$v_gene, function(x){ unlist(strsplit(x, ";")) })
  j_gene <- lapply(meta$j_gene, function(x){ unlist(strsplit(x, ";")) })
  d_gene <- lapply(meta$d_gene, function(x){ unlist(strsplit(x, ";")) })
  c_gene <- lapply(meta$c_gene, function(x){ unlist(strsplit(x, ";")) })
  names(cdr3) <- cells
  names(v_gene) <- cells
  names(j_gene) <- cells
  names(d_gene) <- cells
  names(c_gene) <- cells

  unique_cdr3 <- setdiff(unique(unlist(cdr3)),NA)
  unique_v_gene <- setdiff(unique(unlist(v_gene)),NA)
  unique_j_gene <- setdiff(unique(unlist(j_gene)),NA)
  unique_d_gene <- setdiff(unique(unlist(d_gene)),NA)
  unique_c_gene <- setdiff(unique(unlist(c_gene)),NA)

  cols <- c(unique_cdr3, unique_v_gene, unique_j_gene, unique_d_gene, unique_c_gene)

  tcr <- matrix(NA, nrow = length(cells), ncol = length(cols),
                dimnames = list(cells, cols))

  cores=14
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  tcr_list <- foreach(i=1:length(cells)) %dopar% {
    cell <- cells[i]
    tcr_cell <- tcr[cell,]
    tcr_cell[names(tcr_cell) %in% cdr3[[cell]]] <- T
    tcr_cell[names(tcr_cell) %in% v_gene[[cell]]] <- T
    tcr_cell[names(tcr_cell) %in% j_gene[[cell]]] <- T
    tcr_cell[names(tcr_cell) %in% d_gene[[cell]]] <- T
    tcr_cell[names(tcr_cell) %in% c_gene[[cell]]] <- T
    tcr_cell[is.na(tcr_cell)] <- F
    tcr_cell
  }
  names(tcr_list) <- cells
  tcr <- do.call("rbind",tcr_list)

  stopCluster(cl)

  return(tcr)
}

# Match Vasu sequences to TCR matrix, with different cutoffs -------------------
get_matching_cells_rank <- function(tcr.matches,vdj.takara,rank){
  is_cdr3_match <- colnames(tcr.matches) %in% vdj.takara$cdr3
  is_v_gene_match <- colnames(tcr.matches) %in% vdj.takara$v_gene
  is_d_gene_match <- colnames(tcr.matches) %in% vdj.takara$d_gene
  is_j_gene_match <- colnames(tcr.matches) %in% vdj.takara$j_gene
  is_c_gene_match <- colnames(tcr.matches) %in% vdj.takara$c_gene

  vdj.takara$pass <- vdj.takara$reads.rank <= rank

  is_cdr3_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$cdr3[vdj.takara$pass]
  is_v_gene_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$v_gene[vdj.takara$pass]
  is_d_gene_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$d_gene[vdj.takara$pass]
  is_j_gene_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$j_gene[vdj.takara$pass]
  is_c_gene_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$c_gene[vdj.takara$pass]

  tcr.matches.cdr3 <- tcr.matches[,is_cdr3_match & is_cdr3_match_cutoff]
  tcr.matches.v_gene <- tcr.matches[,is_v_gene_match & is_v_gene_match_cutoff]
  tcr.matches.j_gene <- tcr.matches[,is_j_gene_match & is_j_gene_match_cutoff]

  if (!any(class(tcr.matches.cdr3)=="matrix")){
    cells.matching.cdr3 <- names(which(tcr.matches.cdr3 == T))
  } else {
    cells.matching.cdr3 <- names(which(rowSums(
      tcr.matches[,is_cdr3_match & is_cdr3_match_cutoff])>0))
  }
  if (!any(class(tcr.matches.v_gene)=="matrix")){
    cells.matching.v_gene <- names(which(tcr.matches.v_gene == T))
  } else {
    cells.matching.v_gene <- names(which(rowSums(
      tcr.matches[,is_v_gene_match & is_v_gene_match_cutoff])>0))
  }
  if (!any(class(tcr.matches.j_gene)=="matrix")){
    cells.matching.j_gene <- names(which(tcr.matches.j_gene == T))
  } else {
    cells.matching.j_gene <- names(which(rowSums(
      tcr.matches[,is_j_gene_match & is_j_gene_match_cutoff])>0))
  }

  cells_matching <- intersect(intersect(
    cells.matching.cdr3,
    cells.matching.v_gene),
    cells.matching.j_gene)

  return(cells_matching)
}

get_cells_matching <- function(vdj.takara,sname,chain,tcr.matches){
  um_id <- unique(s.t@meta.data$UM.ID[s.t@meta.data$Sample.ID == sname])
  ranks <- lapply(split(vdj.takara, paste0(vdj.takara$UM.ID, ":",
                                           vdj.takara$chains)),
                  function(x){max(x$reads.rank)})
  max_rank <- ranks[[paste0(um_id, ":", chain)]]

  cells_matching <- list()
  for (rank in 1:max_rank){
    cells_matching[[rank]] <- get_matching_cells_rank(
      tcr.matches[grep(sname,rownames(tcr.matches)),],
      vdj.takara[vdj.takara$UM.ID == um_id & vdj.takara$chains == chain,],
      rank)
  }

  return(cells_matching)
}

get_cells_matching_per_cluster <- function(cells_matching){
  n_matches_rank.per_cluster <- do.call("rbind",lapply(cells_matching, function(x){
    tab <- table(s.t@active.ident[x])
    clusters <- as.character(unique(s.t@active.ident))
    missing <- setdiff(clusters,names(tab))
    tab <- as.data.frame.matrix(t(as.matrix(tab)))
    if (length(missing) > 0){
      tab <- cbind(tab,t(rep(0,length(missing))))
      colnames(tab)[(ncol(tab)-length(missing)+1):ncol(tab)] <- missing
    }
    tab <- tab[,clusters]
    tab
  }))
  n_matches_rank.per_cluster$rank <- as.numeric(rownames(n_matches_rank.per_cluster))
  n_matches_rank.per_cluster <- reshape2::melt(n_matches_rank.per_cluster,
                                               id.vars = "rank")

  return(n_matches_rank.per_cluster)
}

get_cells_per_cluster <- function(um_id){
  tab <- table(s.t@active.ident[s.t@meta.data$UM.ID==um_id])
  tab <- as.data.frame.matrix(t(as.matrix(tab)))
  return(tab)
}

# Find all matches
tcr.matches.all <- get_tcr_matrix(dat.tumor)

cells <- rownames(s.t@meta.data)
tcr.matches <- tcr.matches.all[cells,]

cells_matching <- list()
for (um_id in unique(vdj.takara$UM.ID)){
  if (um_id %in% s.t@meta.data$UM.ID){
    snames <- unique(s.t@meta.data$Sample.ID[s.t@meta.data$UM.ID==um_id])
    for (sname in snames){
      for (chain in c("TRA","TRB")){
        cells_matching[[paste0(um_id, "_", sname, ":", chain)]] <- get_cells_matching(
          vdj.takara, sname, chain, tcr.matches)
      }
    }
  }
}

ranks <- lapply(split(
  vdj.takara, paste0(vdj.takara$UM.ID, ":", vdj.takara$chains)),
  function(x){max(x$reads.rank)})
max_rank <- max(unlist(ranks))

n_matches_rank.per_cluster <- list()
for (um_id_chain in names(cells_matching)){
  n_matches_rank.per_cluster[[um_id_chain]] <- get_cells_matching_per_cluster(
    cells_matching[[um_id_chain]])
  n_matches_rank.per_cluster[[um_id_chain]]$um_id <- str_split_fixed(
    um_id_chain, ":", 2)[,1]
  n_matches_rank.per_cluster[[um_id_chain]]$chain <- str_split_fixed(
    um_id_chain, ":", 2)[,2]
}
n_matches_rank.per_cluster <- do.call("rbind", n_matches_rank.per_cluster)

cells_per_cluster <- list()
for (um_id in unique(n_matches_rank.per_cluster$um_id)){
  cells_per_cluster[[um_id]] <- get_cells_per_cluster(
    str_split_fixed(um_id, "_", 2)[,1])
}
cells_per_cluster <- do.call("rbind", cells_per_cluster)
cells_per_cluster$um_id <- rownames(cells_per_cluster)
cells_per_cluster <- reshape2::melt(cells_per_cluster)

ggplot(n_matches_rank.per_cluster[n_matches_rank.per_cluster$chain == "TRA",],
       aes(x = rank, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") + theme_classic() +
  facet_wrap(~ um_id, scales = "free") +
  theme(strip.background = element_blank()) +
  xlab(NULL) + ylab("Number of cells matching TCR")
ggsave(filename = paste0(outdir, "vasu_tcr.n_matches_rank.per_cluster.tra.pdf"),
       width = 15, height = 3)

ggplot(n_matches_rank.per_cluster[n_matches_rank.per_cluster$chain == "TRB",],
       aes(x = rank, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") + theme_classic() +
  facet_wrap(~ um_id, scales = "free") +
  theme(strip.background = element_blank()) +
  xlab(NULL) + ylab("Number of cells matching TCR")
ggsave(filename = paste0(outdir, "vasu_tcr.n_matches_rank.per_cluster.trb.pdf"),
       width = 15, height = 3)

ggplot(cells_per_cluster, aes(x = um_id, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") + theme_classic() +
  facet_wrap(~ um_id, scales = "free") +
  theme(strip.background = element_blank()) +
  xlab(NULL) + ylab("Number of cells")
ggsave(filename = paste0(outdir, "cells_per_cluster.pdf"),
       width = 5, height = 3)

# Get matches to TCRs that are the same between UM22 and UM46
vdj.takara.um22 <- vdj.takara[vdj.takara$UM.ID == "UM22",]
vdj.takara.um46 <- vdj.takara[vdj.takara$UM.ID == "UM46",]
vdj.takara.um22$cdr3_v_j <- apply(vdj.takara.um22[
  ,c("cdr3","v_gene","j_gene")],1,paste0,collapse=":")
vdj.takara.um46$cdr3_v_j <- apply(vdj.takara.um46[
  ,c("cdr3","v_gene","j_gene")],1,paste0,collapse=":")

shared_cdr3_v_j <- intersect(vdj.takara.um22$cdr3_v_j, vdj.takara.um46$cdr3_v_j)
vdj.takara.um22.shared <- vdj.takara.um22[
  vdj.takara.um22$cdr3_v_j %in% shared_cdr3_v_j,]
vdj.takara.um46.shared <- vdj.takara.um46[
  vdj.takara.um46$cdr3_v_j %in% shared_cdr3_v_j,]

vdj.takara.um46.shared$reads.rank[
  vdj.takara.um46.shared$chains=="TRA"] <- 1:length(
    vdj.takara.um46.shared[vdj.takara.um46.shared$chains=="TRA",]$reads.rank)
vdj.takara.um46.shared$reads.rank[
  vdj.takara.um46.shared$chains=="TRB"] <- 1:length(
    vdj.takara.um46.shared[vdj.takara.um46.shared$chains=="TRB",]$reads.rank)

sname <- unique(s.t@meta.data$Sample.ID[s.t@meta.data$UM.ID=="UM46"])
cells_matching.um46.shared <- list()
for (chain in c("TRA","TRB")){
  cells_matching.um46.shared[[chain]] <- get_cells_matching(
    vdj.takara.um46.shared, sname, chain, tcr.matches)
}

n_matches_rank.per_cluster.um46.shared <- list()
for (chain in names(cells_matching.um46.shared)){
  n_matches_rank.per_cluster.um46.shared[[chain]] <- get_cells_matching_per_cluster(
    cells_matching.um46.shared[[chain]])
  n_matches_rank.per_cluster.um46.shared[[chain]]$um_id <- "UM46"
  n_matches_rank.per_cluster.um46.shared[[chain]]$chain <- chain
}
n_matches_rank.per_cluster.um46.shared <- do.call(
  "rbind", n_matches_rank.per_cluster.um46.shared)

ggplot(n_matches_rank.per_cluster.um46.shared[
  n_matches_rank.per_cluster.um46.shared$chain == "TRA",],
  aes(x = rank, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") + theme_classic() +
  facet_wrap(~ um_id, scales = "free") +
  theme(strip.background = element_blank()) +
  xlab(NULL) + ylab("Number of cells matching TCR")
ggsave(filename = paste0(
  outdir, "vasu_tcr.n_matches_rank.per_cluster.um46.shared.tra.pdf"),
       width = 6, height = 3)

ggplot(n_matches_rank.per_cluster.um46.shared[
  n_matches_rank.per_cluster.um46.shared$chain == "TRB",],
  aes(x = rank, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") + theme_classic() +
  facet_wrap(~ um_id, scales = "free") +
  theme(strip.background = element_blank()) +
  xlab(NULL) + ylab("Number of cells matching TCR")
ggsave(filename = paste0(
  outdir, "vasu_tcr.n_matches_rank.per_cluster.um46.shared.trb.pdf"),
       width = 6, height = 3)

# Map Vasu sequences to biopsies -----------------------------------------------
get_cells_matching_2 <- function(tcr.matches, vdj.takara, chain){
  ranks <- lapply(split(
    vdj.takara, paste0(vdj.takara$UM.ID, ":", vdj.takara$chains)),
    function(x){max(x$reads.rank)})

  cells_matching_all <- list()
  for (sname in unique(str_split_fixed(rownames(tcr.matches),"_",2)[,1])){
    tcr.matches.sname <- tcr.matches[
      str_split_fixed(rownames(tcr.matches),"_",2)[,1] == sname,]

    cells_matching_sname <- list()
    for (um_id in unique(vdj.takara$UM.ID)){
      max_rank <- ranks[[paste0(um_id,":",chain)]]

      #cells_matching <- list()
      #for (rank in 1:max_rank){
      rank <- max_rank
  #      cells_matching[[rank]] <- get_matching_cells_rank(
      cells_matching <- get_matching_cells_rank(
          tcr.matches.sname,
          vdj.takara[vdj.takara$UM.ID==um_id & vdj.takara$chains==chain,],
          rank)
      #}
      cells_matching_sname[[um_id]] <- cells_matching
    }
    cells_matching_all[[sname]] <- cells_matching_sname
  }
  return(cells_matching_all)
}

annotate_matching_tcrs <- function(s.t, cells_matching_tra, cells_matching_trb){
  s.t@meta.data$vasu_match_tra_UM1 <- F
  s.t@meta.data$vasu_match_tra_UM9 <- F
  s.t@meta.data$vasu_match_tra_UM22 <- F
  s.t@meta.data$vasu_match_tra_UM46 <- F
  s.t@meta.data$vasu_match_tra_CD3_41bb <- F
  s.t@meta.data$vasu_match_tra_MART1 <- F

  s.t@meta.data$vasu_match_trb_UM1 <- F
  s.t@meta.data$vasu_match_trb_UM9 <- F
  s.t@meta.data$vasu_match_trb_UM22 <- F
  s.t@meta.data$vasu_match_trb_UM46 <- F
  s.t@meta.data$vasu_match_trb_CD3_41bb <- F
  s.t@meta.data$vasu_match_trb_MART1 <- F

  for (sname in unique(s.t@meta.data$orig.ident)){
    s.t@meta.data$vasu_match_tra_UM1 <- s.t@meta.data$vasu_match_tra_UM1 |
      rownames(s.t@meta.data) %in% cells_matching_tra[[sname]][["UM1"]]
    s.t@meta.data$vasu_match_tra_UM9 <- s.t@meta.data$vasu_match_tra_UM9 |
      rownames(s.t@meta.data) %in% cells_matching_tra[[sname]][["UM9"]]
    s.t@meta.data$vasu_match_tra_UM22 <- s.t@meta.data$vasu_match_tra_UM22 |
      rownames(s.t@meta.data) %in% cells_matching_tra[[sname]][["UM22"]]
    s.t@meta.data$vasu_match_tra_UM46 <- s.t@meta.data$vasu_match_tra_UM46 |
      rownames(s.t@meta.data) %in% cells_matching_tra[[sname]][["UM46"]]
    s.t@meta.data$vasu_match_tra_CD3_41bb <- s.t@meta.data$vasu_match_tra_CD3_41bb |
      rownames(s.t@meta.data) %in% cells_matching_tra[[sname]][["UM1"]] |
      rownames(s.t@meta.data) %in% cells_matching_tra[[sname]][["UM9"]]
    s.t@meta.data$vasu_match_tra_MART1 <- s.t@meta.data$vasu_match_tra_MART1 |
      rownames(s.t@meta.data) %in% cells_matching_tra[[sname]][["UM22"]] |
      rownames(s.t@meta.data) %in% cells_matching_tra[[sname]][["UM46"]]

    s.t@meta.data$vasu_match_trb_UM1 <- s.t@meta.data$vasu_match_trb_UM1 |
      rownames(s.t@meta.data) %in% cells_matching_trb[[sname]][["UM1"]]
    s.t@meta.data$vasu_match_trb_UM9 <- s.t@meta.data$vasu_match_trb_UM9 |
      rownames(s.t@meta.data) %in% cells_matching_trb[[sname]][["UM9"]]
    s.t@meta.data$vasu_match_trb_UM22 <- s.t@meta.data$vasu_match_trb_UM22 |
      rownames(s.t@meta.data) %in% cells_matching_trb[[sname]][["UM22"]]
    s.t@meta.data$vasu_match_trb_UM46 <- s.t@meta.data$vasu_match_trb_UM46 |
      rownames(s.t@meta.data) %in% cells_matching_trb[[sname]][["UM46"]]
    s.t@meta.data$vasu_match_trb_CD3_41bb <- s.t@meta.data$vasu_match_trb_CD3_41bb |
      rownames(s.t@meta.data) %in% cells_matching_trb[[sname]][["UM1"]] |
      rownames(s.t@meta.data) %in% cells_matching_trb[[sname]][["UM9"]]
    s.t@meta.data$vasu_match_trb_MART1 <- s.t@meta.data$vasu_match_trb_MART1 |
      rownames(s.t@meta.data) %in% cells_matching_trb[[sname]][["UM22"]] |
      rownames(s.t@meta.data) %in% cells_matching_trb[[sname]][["UM46"]]
  }

  s.t@meta.data$vasu_match_CD3_41bb <- s.t@meta.data$vasu_match_tra_CD3_41bb | s.t@meta.data$vasu_match_trb_CD3_41bb
  s.t@meta.data$vasu_match_MART1 <- s.t@meta.data$vasu_match_tra_MART1 | s.t@meta.data$vasu_match_trb_MART1
  s.t@meta.data$vasu_match_UM1 <- s.t@meta.data$vasu_match_tra_UM1 | s.t@meta.data$vasu_match_trb_UM1
  s.t@meta.data$vasu_match_UM9 <- s.t@meta.data$vasu_match_tra_UM9 | s.t@meta.data$vasu_match_trb_UM9
  s.t@meta.data$vasu_match_UM22 <- s.t@meta.data$vasu_match_tra_UM22 | s.t@meta.data$vasu_match_trb_UM22
  s.t@meta.data$vasu_match_UM46 <- s.t@meta.data$vasu_match_tra_UM46 | s.t@meta.data$vasu_match_trb_UM46
  s.t@meta.data$vasu_match_CD3_41bb_and_MART1 <- s.t@meta.data$vasu_match_CD3_41bb & s.t@meta.data$vasu_match_MART1

  s.t@meta.data$vasu_phenotype <- paste0(
    ifelse(s.t@meta.data$vasu_match_CD3_41bb,"CD3_41bb",""),
    ifelse(s.t@meta.data$vasu_match_CD3_41bb_and_MART1,":",""),
    ifelse(s.t@meta.data$vasu_match_MART1,"MART1",""))
  s.t@meta.data$vasu_phenotype[s.t@meta.data$vasu_phenotype==""] <- "None"
  s.t@meta.data$vasu_phenotype <- factor(
    s.t@meta.data$vasu_phenotype,
    levels=c("None","CD3_41bb","MART1","CD3_41bb:MART1"))

  s.t@meta.data$vasu_phenotype_um <- paste0(as.character(
    s.t@meta.data$vasu_phenotype),":",
         ifelse(s.t@meta.data$vasu_match_UM1,"UM1:",""),
         ifelse(s.t@meta.data$vasu_match_UM9,"UM9:",""),
         ifelse(s.t@meta.data$vasu_match_UM22,"UM22:",""),
         ifelse(s.t@meta.data$vasu_match_UM46,"UM46",""))

  s.t@meta.data$vasu_phenotype_um <- gsub(":$","",s.t@meta.data$vasu_phenotype_um)
  s.t@meta.data$vasu_phenotype_um <- factor(
    s.t@meta.data$vasu_phenotype_um,
    levels=c("None",setdiff(sort(unique(s.t@meta.data$vasu_phenotype_um)),"None")))

  return(s.t)
}

cells_matching_tra <- get_cells_matching_2(tcr.matches, vdj.takara, chain="TRA")
cells_matching_trb <- get_cells_matching_2(tcr.matches, vdj.takara, chain="TRB")

s.t <- annotate_matching_tcrs(s.t, cells_matching_tra, cells_matching_trb)

pdf(file=paste0(outdir,"biopsy.vasu_phenotype_um.pdf"), width = 35, height = 10)
DimPlot(s.t, group.by="vasu_phenotype_um", order = T, split.by = "UM.ID", ncol=7,
        pt.size = 1)
dev.off()

s.t.common <- subset(s.t, UM.ID %in% intersect(
  unique(dat.tumor@meta.data$UM.ID), names(tils.unintegrated)))
s.t.common@meta.data$UM.ID <- factor(
  s.t.common@meta.data$UM.ID,
  levels= unique(s.t.common@meta.data$UM.ID)[order(
    as.numeric(gsub("UM","",unique(s.t.common@meta.data$UM.ID))),
    decreasing = F)])

pdf(file=paste0(outdir,"biopsy.common_samples.vasu_phenotype_um.pdf"), width = 25, height = 5)
DimPlot(s.t.common, group.by="vasu_phenotype_um", order = T, split.by = "UM.ID",
        ncol=7, pt.size = 1)
dev.off()

rm(s.t.common)
gc()

tcr.matches.til.unintegrated <- lapply(tils.unintegrated, get_tcr_matrix)

for (sname in names(tcr.matches.til.unintegrated)){
  cells_matching_tra_til <- get_cells_matching_2(
    tcr.matches.til.unintegrated[[sname]], vdj.takara, chain="TRA")
  cells_matching_trb_til <- get_cells_matching_2(
    tcr.matches.til.unintegrated[[sname]], vdj.takara, chain="TRB")
  tils.unintegrated[[sname]] <- annotate_matching_tcrs(
    tils.unintegrated[[sname]], cells_matching_tra_til, cells_matching_trb_til)
}

lvls <- c("None",setdiff(sort(unique(unlist(lapply(tils.unintegrated,
  function(x){
    unique(as.character(x@meta.data$vasu_phenotype_um))
  })))),"None"))

tils.unintegrated <- lapply(tils.unintegrated, function(x){
  x@meta.data$vasu_phenotype_um <- factor(
    x@meta.data$vasu_phenotype_um, levels=lvls);
  return(x)})

p_um1 <- DimPlot(tils.unintegrated[["UM1"]], split.by = "UM.ID",
              group.by = "vasu_phenotype_um",
              order = T,pt.size = 1) +
  theme_void() + ggtitle(NULL)
p_um2 <- DimPlot(tils.unintegrated[["UM2"]], split.by = "UM.ID",
              group.by = "vasu_phenotype_um",
              order = T,pt.size = 1) +
  theme_void() + ggtitle(NULL)
p_um3 <- DimPlot(tils.unintegrated[["UM3"]],
              group.by = "vasu_phenotype_um", split.by = "UM.ID",
              order = T,pt.size = 1) +
  theme_void() + ggtitle(NULL)
p_um9 <- DimPlot(tils.unintegrated[["UM9"]], split.by = "UM.ID",
              group.by = "vasu_phenotype_um",
              order = T,pt.size = 1) +
  theme_void() + ggtitle(NULL)
p_um19 <- DimPlot(tils.unintegrated[["UM19"]], split.by = "UM.ID",
              group.by = "vasu_phenotype_um",
              order = T,pt.size = 1) +
  theme_void() + ggtitle(NULL)
p_um22 <- DimPlot(tils.unintegrated[["UM22"]], split.by = "UM.ID",
                  group.by = "vasu_phenotype_um",
                  order = T,pt.size = 1) +
  theme_void() + ggtitle(NULL)
p_um24 <- DimPlot(tils.unintegrated[["UM24"]], split.by = "UM.ID",
              group.by = "vasu_phenotype_um",
              order = T,pt.size = 1) +
  theme_void() + ggtitle(NULL)
p_um10 <- DimPlot(tils.unintegrated[["UM10"]], split.by = "UM.ID",
                  group.by = "vasu_phenotype_um",
                  order = T,pt.size = 1) +
  theme_void() + ggtitle(NULL)

pdf(file=paste0(outdir,"til.vasu_phenotype_um.pdf"), width = 35, height = 10)
p.til <- grid.arrange(grobs = c(list(p_um1),
                                list(p_um2),
                                list(p_um3),
                                list(p_um9),
                                list(p_um19),
                                list(p_um22),
                                list(p_um24),
                                list(p_um10)), ncol = 8, as.table = FALSE)
dev.off()

pdf(file=paste0(outdir,"til.common.vasu_phenotype_um.pdf"), width = 35, height = 5)
p.til <- grid.arrange(grobs = c(list(p_um1),
                                list(p_um2),
                                list(p_um3),
                                list(p_um9),
                                list(p_um19),
                                list(p_um24)), ncol = 6, as.table = FALSE)
dev.off()

# Make TCR matrix for biopsies and TILs combined -------------------------------
dat.cca.til.with_cdr3 <- subset(dat.til, cells = rownames(dat.til@meta.data[
                                  !is.na(dat.til@meta.data$cdr3),]))

dat.tumor@meta.data$Sample.ID <- dat.tumor@meta.data$orig.ident
dat.tumor.with_cdr3 <- subset(dat.tumor, cells = rownames(dat.tumor@meta.data[
                                !is.na(dat.tumor@meta.data$cdr3),]))

meta.tumor <- dat.tumor.with_cdr3@meta.data[,c("Sample.ID", "cdr3", "v_gene",
                                               "j_gene", "d_gene", "c_gene")]
meta.til <- dat.cca.til.with_cdr3@meta.data[,c("Sample.ID", "cdr3", "v_gene",
                                               "j_gene", "d_gene", "c_gene")]

meta.tumor_til <- rbind(meta.tumor, meta.til)

tcr.matches <- get_tcr_matrix(meta.tumor_til)

tcr <- tcr.matches[rownames(tcr.matches) %in% rownames(s.t@meta.data),]

#saveRDS(tcr.matches, file = paste0(outdir,"til_tumor.tcr_matric.rda"))
#saveRDS(tcr, file = paste0(outdir,"s.t.tcr_matric.rda"))

# Make UMAP plots for TILs that show TCRs matching biopsies --------------------
cells <- rownames(meta.tumor_til)
cdr3 <- lapply(meta.tumor_til$cdr3, function(x){unlist(strsplit(x,";"))})
v_gene <- lapply(meta.tumor_til$v_gene, function(x){unlist(strsplit(x,";"))})
j_gene <- lapply(meta.tumor_til$j_gene, function(x){unlist(strsplit(x,";"))})
d_gene <- lapply(meta.tumor_til$d_gene, function(x){unlist(strsplit(x,";"))})
c_gene <- lapply(meta.tumor_til$c_gene, function(x){unlist(strsplit(x,";"))})
names(cdr3) <- cells
names(v_gene) <- cells
names(j_gene) <- cells
names(d_gene) <- cells
names(c_gene) <- cells

unique_cdr3 <- unique(unlist(cdr3))
unique_v_gene <- unique(unlist(v_gene))
unique_j_gene <- unique(unlist(j_gene))
unique_d_gene <- unique(unlist(d_gene))
unique_c_gene <- unique(unlist(c_gene))

is_cdr3 <- colnames(tcr.matches) %in% unique_cdr3
is_v_gene <- colnames(tcr.matches) %in% unique_v_gene
is_d_gene <- colnames(tcr.matches) %in% unique_d_gene
is_j_gene <- colnames(tcr.matches) %in% unique_j_gene
is_c_gene <- colnames(tcr.matches) %in% unique_c_gene

cells_til <- rownames(meta.tumor_til)[grep("Sample", meta.tumor_til$Sample.ID)]
cells_tumor <- rownames(meta.tumor_til)[grep("Sample", meta.tumor_til$Sample.ID,
                                             invert = T)]

idx_til <- rownames(tcr.matches) %in% cells_til
idx_tumor <- rownames(tcr.matches) %in% cells_tumor

tcr.matches.til <- tcr.matches[idx_til,]
tcr.matches.tumor <- tcr.matches[idx_tumor,]

cores=14
cl <- makeCluster(cores)
registerDoParallel(cl)
til_matches_tumor <- foreach (i = 1:nrow(tcr.matches.tumor)) %dopar% {
  tcr.matches.til.cdr3 <- tcr.matches.til[, is_cdr3] & tcr.matches.tumor[i, is_cdr3]
  tcr.matches.til.v <- tcr.matches.til[, is_v_gene] & tcr.matches.tumor[i, is_v_gene]
  tcr.matches.til.j <- tcr.matches.til[, is_j_gene] & tcr.matches.tumor[i, is_j_gene]

  idx_til_matches_tumor_i <- which(rowSums(tcr.matches.til.cdr3) > 0 &
                                     rowSums(tcr.matches.til.v) > 0 &
                                     rowSums(tcr.matches.til.j) > 0)
  rownames(tcr.matches.til)[idx_til_matches_tumor_i]
}
stopCluster(cl)

#saveRDS(til_matches_tumor, file = paste0(outdir,"til_matches_tumor.rda"))

all_tils_matching_any_tumor <- sort(unique(unlist(til_matches_tumor)))

for (i in 1:length(tils.unintegrated)){
  cells <- rownames(tils.unintegrated[[i]]@meta.data)
  cells <- setNames(rep(F,length(cells)),cells)
  cells[names(cells) %in% all_tils_matching_any_tumor] <- T

  tils.unintegrated[[i]] <- AddMetaData(tils.unintegrated[[i]], metadata = cells,
                                        col.name = "tcr_matching_any_biopsy")
}

p_um1 <- DimPlot(tils.unintegrated[["UM1"]],
                 group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                 cols = c("gray", "blue")) + ggtitle("UM1")
p_um2 <- DimPlot(tils.unintegrated[["UM2"]],
                 group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                 cols = c("gray", "blue")) + ggtitle("UM2")
p_um3 <- DimPlot(tils.unintegrated[["UM3"]],
                 group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                 cols = c("gray", "blue")) + ggtitle("UM3")
p_um9 <- DimPlot(tils.unintegrated[["UM9"]],
                 group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                 cols = c("gray", "blue")) + ggtitle("UM9")
p_um10 <- DimPlot(tils.unintegrated[["UM10"]],
                  group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM10")
p_um19 <- DimPlot(tils.unintegrated[["UM19"]],
                  group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM19")
p_um22 <- DimPlot(tils.unintegrated[["UM22"]],
                  group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM22")
p_um24 <- DimPlot(tils.unintegrated[["UM24"]],
                  group.by = "tcr_matching_any_biopsy", order = T, pt.size = 1,
                  cols = c("gray", "blue")) + ggtitle("UM24")

library(gridExtra)
pdf(file = paste0(outdir,"til_tcr_matching_any_biopsy.pdf"),width = 55, height = 6)
p <- grid.arrange(grobs = c(list(p_um1),
                            list(p_um2),
                            list(p_um3),
                            list(p_um9),
                            list(p_um19),
                            list(p_um10),
                            list(p_um24),
                            list(p_um22)), ncol = 8, as.table = FALSE)
dev.off()

# Pair-wise overlap in TCRs among samples --------------------------------------
get_matches_pair.v3 <- function(dat.tumor_til.vdj, sname_1, sname_2,
                                match_on = c("cdr3","v","d","j","c")){
  idx_1 <- unlist(lapply(dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3, function(x){
    all(x!="None")
  }))

  idx_2 <- unlist(lapply(dat.tumor_til.vdj[[sname_2]]$vdj.10x_cdr3, function(x){
    all(x!="None")
  }))

  tmp_1.10x_cdr3 <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3[which(idx_1)]
  tmp_2.10x_cdr3 <- dat.tumor_til.vdj[[sname_2]]$vdj.10x_cdr3[which(idx_2)]

  tmp_1.10x_v_gene <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_v_gene[which(idx_1)]
  tmp_2.10x_v_gene <- dat.tumor_til.vdj[[sname_2]]$vdj.10x_v_gene[which(idx_2)]

  tmp_1.10x_d_gene <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_d_gene[which(idx_1)]
  tmp_2.10x_d_gene <- dat.tumor_til.vdj[[sname_2]]$vdj.10x_d_gene[which(idx_2)]

  tmp_1.10x_j_gene <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_j_gene[which(idx_1)]
  tmp_2.10x_j_gene <- dat.tumor_til.vdj[[sname_2]]$vdj.10x_j_gene[which(idx_2)]

  tmp_1.10x_c_gene <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_c_gene[which(idx_1)]
  tmp_2.10x_c_gene <- dat.tumor_til.vdj[[sname_2]]$vdj.10x_c_gene[which(idx_2)]

  which_sname_2.10x_cdr3_in_sname_1.10x_cdr3 <- unlist(
    lapply(tmp_1.10x_cdr3,function(x){
    which(unlist(lapply(tmp_2.10x_cdr3,function(y){any(y %in% x)})))
  }))
  which_sname_2.10x_v_gene_in_sname_1.10x_v_gene <- unlist(
    lapply(tmp_1.10x_v_gene,function(x){
    which(unlist(lapply(tmp_2.10x_v_gene,function(y){any(y %in% x)})))
  }))
  which_sname_2.10x_d_gene_in_sname_1.10x_d_gene <- unlist(
    lapply(tmp_1.10x_d_gene,function(x){
    which(unlist(lapply(tmp_2.10x_d_gene,function(y){any(y %in% x)})))
  }))
  which_sname_2.10x_j_gene_in_sname_1.10x_j_gene <- unlist(
    lapply(tmp_1.10x_j_gene,function(x){
    which(unlist(lapply(tmp_2.10x_j_gene,function(y){any(y %in% x)})))
  }))
  which_sname_2.10x_c_gene_in_sname_1.10x_c_gene <- unlist(
    lapply(tmp_1.10x_c_gene,function(x){
    which(unlist(lapply(tmp_2.10x_c_gene,function(y){any(y %in% x)})))
  }))

  all_pairs <- c(names(which_sname_2.10x_cdr3_in_sname_1.10x_cdr3),
                 names(which_sname_2.10x_v_gene_in_sname_1.10x_v_gene),
                 names(which_sname_2.10x_d_gene_in_sname_1.10x_d_gene),
                 names(which_sname_2.10x_j_gene_in_sname_1.10x_j_gene),
                 names(which_sname_2.10x_c_gene_in_sname_1.10x_c_gene))

  unique_pairs <- unique(all_pairs)

  df <- data.frame(pair = unique_pairs,
                   in_10x_cdr3 = unique_pairs %in%
                     names(which_sname_2.10x_cdr3_in_sname_1.10x_cdr3),
                   in_10x_v_gene = unique_pairs %in%
                     names(which_sname_2.10x_v_gene_in_sname_1.10x_v_gene),
                   in_10x_d_gene = unique_pairs %in%
                     names(which_sname_2.10x_d_gene_in_sname_1.10x_d_gene),
                   in_10x_j_gene = unique_pairs %in%
                     names(which_sname_2.10x_j_gene_in_sname_1.10x_j_gene),
                   in_10x_c_gene = unique_pairs %in%
                     names(which_sname_2.10x_c_gene_in_sname_1.10x_c_gene),
                   shared_cdr3 = "",
                   shared_v_gene = "",
                   shared_d_gene = "",
                   shared_j_gene = "",
                   shared_c_gene = "")

  if (! any(c("cdr3","v","d","j","c") %in% match_on)){
    stop("At least on of cdr3, v, d, j, c must be in match_on")
  }
  if (! any(c("cdr3") %in% match_on)){
    stop("cdr3 must be in match_on")
  }

  if ("cdr3" %in% match_on){
    df$match <- df$in_10x_cdr3
  }
  if ("v" %in% match_on){
    df$match <- df$match & df$in_10x_v_gene
  }
  if ("d" %in% match_on){
    df$match <- df$match & df$in_10x_d_gene
  }
  if ("j" %in% match_on){
    df$match <- df$match & df$in_10x_j_gene
  }
  if ("c" %in% match_on){
    df$match <- df$match & df$in_10x_c_gene
  }

  #df$match <- df$in_10x_cdr3 & df$in_10x_v_gene & df$in_10x_d_gene & df$in_10x_j_gene & df$in_10x_c_gene
  df$cell_1 <- str_split_fixed(df$pair,"\\.",2)[,1]
  df$cell_2 <- str_split_fixed(df$pair,"\\.",2)[,2]

  df.match <- df[df$match,]

  if (nrow(df.match)>0){
    for (i in 1:nrow(df.match)){
      cell_1 <- df.match$cell_1[i]
      cell_2 <- df.match$cell_2[i]
      df.match$shared_cdr3[i] <- paste0(
        sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3[[cell_1]],
                       dat.tumor_til.vdj[[sname_2]]$vdj.10x_cdr3[[cell_2]])),
        collapse = ";")
      df.match$shared_v_gene[i] <- paste0(
        sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_v_gene3[[cell_1]],
                       dat.tumor_til.vdj[[sname_2]]$vdj.10x_v_gene[[cell_2]])),
        collapse = ";")
      df.match$shared_d_gene[i] <- paste0(
        sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_d_gene[[cell_1]],
                       dat.tumor_til.vdj[[sname_2]]$vdj.10x_d_gene[[cell_2]])),
        collapse = ";")
      df.match$shared_j_gene[i] <- paste0(
        sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_j_gene[[cell_1]],
                       dat.tumor_til.vdj[[sname_2]]$vdj.10x_j_gene[[cell_2]])),
        collapse = ";")
      df.match$shared_c_gene[i] <- paste0(
        sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_c_gene[[cell_1]],
                       dat.tumor_til.vdj[[sname_2]]$vdj.10x_c_gene[[cell_2]])),
        collapse = ";")
    }
  }
  return(df.match)
}

dat.tumor.split <- SplitObject(dat.tumor, split.by = "UM.ID")
tils.unintegrated.split <- tils.unintegrated
names(dat.tumor.split) <- paste0(names(dat.tumor.split),"_biopsy")
names(tils.unintegrated.split) <- paste0(names(tils.unintegrated.split),"_til")
dat.tumor_til.split <- c(dat.tumor.split, tils.unintegrated.split)
rm(dat.tumor.split)
rm(tils.unintegrated.split)
gc()

dat.tumor_til.vdj <- lapply(dat.tumor_til.split, extract_tcr)
#saveRDS(dat.tumor_til.vdj,paste0(outdir,"dat.tumor_til.vdj.rda"))

snames <- names(dat.tumor_til.vdj)
snames_grid <- as.data.frame(expand.grid(snames,snames),stringsAsFactors=F)
snames_grid[,1] <- as.character(snames_grid[,1])
snames_grid[,2] <- as.character(snames_grid[,2])
snames_grid <- snames_grid[snames_grid[,1] != snames_grid[,2],]
for (i in 1:nrow(snames_grid)){
  snames_grid[i,] <- sort(as.character(snames_grid[i,]))
}
snames_grid <- unique(snames_grid)
dim(snames_grid)
rownames(snames_grid) <- NULL

gc()
gc()
gc()
gc()
gc()

library(foreach)
library(doParallel)
cores=8
cl <- makeCluster(cores)
registerDoParallel(cl)

matches_pair <- list()
matches_pair <- foreach(i=1:nrow(snames_grid),.packages = c("tidyverse")) %dopar% {
  print(snames_grid[i,])
  sname_1 <- snames_grid[i,1]
  sname_2 <- snames_grid[i,2]
  matches_pair[[paste0(sname_1,":",sname_2)]] <- get_matches_pair.v3(
    dat.tumor_til.vdj, sname_1 = sname_1, sname_2 = sname_2)
}

stopCluster(cl)

names(matches_pair) <- apply(snames_grid,1,paste0,collapse=":")
#saveRDS(matches_pair,file = paste0(outdir,"matches_pair.rda"))

matches_pair.clones <- lapply(matches_pair,function(x){
  unique(x[c("shared_cdr3","shared_v_gene","shared_d_gene","shared_j_gene",
             "shared_c_gene")])
})

matches_pair.nr <- as.data.frame(unlist(lapply(matches_pair.clones, nrow)))
colnames(matches_pair.nr) <- "nr"
matches_pair.nr$sample_1 <- str_split_fixed(rownames(matches_pair.nr),":",2)[,1]
matches_pair.nr$sample_2 <- str_split_fixed(rownames(matches_pair.nr),":",2)[,2]
rownames(matches_pair.nr) <- NULL

unique_samples <- sort(unique(c(str_split_fixed(names(matches_pair),":",2)[,1],
                                str_split_fixed(names(matches_pair),":",2)[,2])))
unique_samples.til <- unique_samples[grep("_til",unique_samples)]
unique_samples.til <- unique_samples.til[order(as.numeric(str_split_fixed(
  gsub("UM","",unique_samples.til),"_til",2)[,1]),decreasing = F)]

unique_samples.biopsy <- unique_samples[grep("_biopsy",unique_samples)]
unique_samples.biopsy <- unique_samples.biopsy[order(as.numeric(str_split_fixed(
  gsub("UM","",unique_samples.biopsy),"_biopsy",2)[,1]),decreasing = F)]

unique_samples <- c(unique_samples.til, unique_samples.biopsy)

matches_pair.nr.m <- matrix(NA,nrow = length(unique_samples),
                            ncol=length(unique_samples),
                            dimnames = list(unique_samples,unique_samples))
for (sname_1 in unique_samples){
  for (sname_2 in unique_samples){
    if (any(matches_pair.nr$sample_1==sname_1 & matches_pair.nr$sample_2==sname_2)){
      matches_pair.nr.m[sname_1,sname_2] <- matches_pair.nr$nr[
        matches_pair.nr$sample_1==sname_1 & matches_pair.nr$sample_2==sname_2]
    } else if (any(matches_pair.nr$sample_1==sname_2 & matches_pair.nr$sample_2==sname_1)){
      matches_pair.nr.m[sname_1,sname_2] <- matches_pair.nr$nr[
        matches_pair.nr$sample_1==sname_2 & matches_pair.nr$sample_2==sname_1]
    }
  }
}

library(pheatmap)

annot <- data.frame(sname = unique_samples, 
                    type = str_split_fixed(unique_samples,"_",2)[,2])
rownames(annot) <- annot$sname
annot$sname <- NULL

pdf(file=paste0(outdir,"tcr_overlap.pdf"), width = 6, height = 5)
pheatmap(matches_pair.nr.m, cluster_rows = F, cluster_cols = F,
         display_numbers = matches_pair.nr.m, border_color = NA,
         show_rownames = T, annotation_row = annot,annotation_col = annot,
         show_colnames = T, na_col = "white", colorRampPalette(c("white", "red"))(100))
dev.off()

# Matches only on CDR3 and only for subset of Vasu TCRs occuring twice or more ----
library(foreach)
library(doParallel)
cores=8
cl <- makeCluster(cores)
registerDoParallel(cl)

matches_pair <- list()
matches_pair <- foreach(i=1:nrow(snames_grid),.packages = c("tidyverse")) %dopar% {
  print(snames_grid[i,])
  sname_1 <- snames_grid[i,1]
  sname_2 <- snames_grid[i,2]
  matches_pair[[paste0(sname_1,":",sname_2)]] <- get_matches_pair.v3(
    dat.tumor_til.vdj, sname_1 = sname_1, sname_2 = sname_2, match_on = "cdr3")
}

stopCluster(cl)

names(matches_pair) <- apply(snames_grid,1,paste0,collapse=":")
saveRDS(matches_pair, file = paste0(outdir, "matches_pair.cdr3_only.rda"))

#names(matches_pair) <- apply(snames_grid, 1, paste0, collapse=":")
#matches_pair.bak <- matches_pair

# Subset to Vasu TCRs (occurring at least twice)
vasu_tcrs_at_least_twice.cdr3 <- names(which(table(vdj.takara$cdr3) >= 2))

matches_pair.vasu <- lapply(matches_pair,function(x){
  is_vasu <- unlist(lapply(strsplit(x$shared_cdr3,";"),function(y){
    any(y %in% vasu_tcrs_at_least_twice.cdr3)
  }))
  if (!is.null(is_vasu)){
    out <- x[is_vasu,]
  } else {
    out <- x[0,]
  }
  return(out)
})

matches_pair.clones <- lapply(matches_pair.vasu,function(x){
  unique(x[c("shared_cdr3","shared_v_gene","shared_d_gene","shared_j_gene",
             "shared_c_gene")])
})

matches_pair.nr <- as.data.frame(unlist(lapply(matches_pair.clones,nrow)))
colnames(matches_pair.nr) <- "nr"
matches_pair.nr$sample_1 <- str_split_fixed(rownames(matches_pair.nr),":",2)[,1]
matches_pair.nr$sample_2 <- str_split_fixed(rownames(matches_pair.nr),":",2)[,2]
rownames(matches_pair.nr) <- NULL

unique_samples <- sort(unique(c(str_split_fixed(names(matches_pair),":",2)[,1],
                                str_split_fixed(names(matches_pair),":",2)[,2])))
unique_samples.til <- unique_samples[grep("_til",unique_samples)]
unique_samples.til <- unique_samples.til[order(as.numeric(str_split_fixed(
  gsub("UM","",unique_samples.til),"_til",2)[,1]),decreasing = F)]

unique_samples.biopsy <- unique_samples[grep("_biopsy",unique_samples)]
unique_samples.biopsy <- unique_samples.biopsy[order(as.numeric(str_split_fixed(
  gsub("UM","",unique_samples.biopsy),"_biopsy",2)[,1]),decreasing = F)]

unique_samples <- c(unique_samples.til, unique_samples.biopsy)

matches_pair.nr.m <- matrix(NA,nrow = length(unique_samples),
                            ncol=length(unique_samples),
                            dimnames = list(unique_samples,unique_samples))
for (sname_1 in unique_samples){
  for (sname_2 in unique_samples){
    if (any(matches_pair.nr$sample_1==sname_1 & matches_pair.nr$sample_2==sname_2)){
      matches_pair.nr.m[sname_1,sname_2] <- matches_pair.nr$nr[
        matches_pair.nr$sample_1==sname_1 & matches_pair.nr$sample_2==sname_2]
    } else if (any(matches_pair.nr$sample_1==sname_2 & matches_pair.nr$sample_2==sname_1)){
      matches_pair.nr.m[sname_1,sname_2] <- matches_pair.nr$nr[
        matches_pair.nr$sample_1==sname_2 & matches_pair.nr$sample_2==sname_1]
    }
  }
}

annot <- data.frame(sname=unique_samples,type=str_split_fixed(unique_samples,"_",2)[,2])
rownames(annot) <- annot$sname
annot$sname <- NULL

pdf(file=paste0(outdir,"tcr_overlap.cdr3.vasu.pdf"),width = 6,height = 5)
pheatmap(matches_pair.nr.m, cluster_rows = F, cluster_cols = F,
         display_numbers = matches_pair.nr.m, border_color = NA,
         annotation_col = annot, annotation_row = annot, show_rownames = T,
         show_colnames = T, na_col = "white", colorRampPalette(c("white", "red"))(100))
dev.off()

# Do a similar overlap plot for only the Vasu data -----------------------------
snames <- c("UM1","UM9","UM22","UM46")
overlap_tcr_vasu_tra <- matrix(NA, nrow=length(snames), ncol=length(snames),
                               dimnames=list(snames,snames))
for (sname_1 in snames){
  for (sname_2 in snames){
    overlap_tcr_vasu_tra[sname_1,sname_2] <- length(
      intersect(unique(vdj.takara$cdr3[vdj.takara$UM.ID==sname_1 & vdj.takara$chains=="TRA"]),
                unique(vdj.takara$cdr3[vdj.takara$UM.ID==sname_2 & vdj.takara$chains=="TRA"])))
  }
}

overlap_tcr_vasu_trb <- matrix(NA, nrow=length(snames), ncol=length(snames),
                               dimnames=list(snames,snames))
for (sname_1 in snames){
  for (sname_2 in snames){
    overlap_tcr_vasu_trb[sname_1,sname_2] <- length(
      intersect(unique(vdj.takara$cdr3[vdj.takara$UM.ID==sname_1 & vdj.takara$chains=="TRB"]),
                unique(vdj.takara$cdr3[vdj.takara$UM.ID==sname_2 & vdj.takara$chains=="TRB"])))
  }
}

pdf(file = paste0(outdir,"overlap_tcr_vasu_tra.pdf"), width = 3, height = 2.5)
pheatmap(overlap_tcr_vasu_tra,cluster_rows = F,cluster_cols = F,
         display_numbers = overlap_tcr_vasu_tra,border_color = NA,
         show_rownames = T, show_colnames = T,
         na_col = "white",colorRampPalette(c("white", "red"))(100))
dev.off()

pdf(file = paste0(outdir,"overlap_tcr_vasu_trb.pdf"), width = 3, height = 2.5)
pheatmap(overlap_tcr_vasu_trb,cluster_rows = F,cluster_cols = F,
         display_numbers = overlap_tcr_vasu_trb,border_color = NA,
         show_rownames = T, show_colnames = T,
         na_col = "white",colorRampPalette(c("white", "red"))(100))
dev.off()

# Overrepresentation test, Vasu data -------------------------------------------

# Can be considered a 7-sided dice (7 clusters), where each side has a different probability.
# This probability is estimated by checking how often TILs end up (i.e. what fraction of all TILs
# are present) in a given cluster. Everything else equal, rolling the "dice" a number of times
# equal to the number of Vasu cells in the dataset should give a similar proportion out of
# the Vasu cells that end up in a given cluster as the proportion out of the TILs that end up there
# (the expected probability). The expected biology of the TILs and the base cluster sizes would be
# baked into this expected probability. Since the base cluster sizes are the same, the only meaningful
# difference would be a biological difference between the Vasu cells and the TILs.

bt <- list()
for (ident in unique(s.t@active.ident)){
  expected_probability <- sum(s.t@active.ident == ident & s.t@meta.data$tcr_matches_til) /
    sum(s.t@meta.data$tcr_matches_til)
  #n_rolls <- sum(s.t@meta.data$vasu_phenotype %in% c("MART1","CD3+41bb+"))
  #n_successes <- sum(s.t@meta.data$vasu_phenotype %in% c("MART1","CD3+41bb+") &
  #                     s.t@active.ident == ident)
  n_rolls <- sum(s.t@meta.data$vasu_phenotype != "None")
  n_successes <- sum(s.t@meta.data$vasu_phenotype != "None" &
                       s.t@active.ident == ident)
  bt[[ident]] <- binom.test(n_successes, n_rolls, expected_probability,
                            alternative = "two.sided")
}

p_null <- unlist(lapply(bt, function(x){ as.numeric(x$null.value) }))
p_estimate <- unlist(lapply(bt, function(x){ as.numeric(x$estimate) }))
p <- unlist(lapply(bt, function(x){ x$p.value }))
q <- p.adjust(p, method = "bonferroni") # probably not independent tests (the same dice rolls for all)

res <- data.frame(p_null, p_estimate, p, q, vasu_greater = p_estimate > p_null,
                  sig = q < 0.05)

library(WriteXLS)
WriteXLS(res,
         ExcelFileName = paste0(outdir,"proportions_til_vs_vasu_matches_per_cluster.xlsx"),
         AdjWidth = T, AutoFilter = T, BoldHeaderRow = T, row.names = T)

# Although, perhaps one would either need to subset to the set of samples included in vasus experiment
# or do it on a sample by sample basis. Or maybe not, since the overall question is whether TILs
# tend to come from exhauster clusters as a general principle. Might be worth further consideration.

proportion_til <- list()
proportion_vasu <- list()
for (ident in unique(s.t@active.ident)){
  proportion_til[[ident]] <- sum(s.t@active.ident == ident &
                                   s.t@meta.data$tcr_matches_til) /
    sum(s.t@meta.data$tcr_matches_til)
  # proportion_vasu[[ident]] <- sum(
  #   s.t@active.ident == ident & s.t@meta.data$vasu_phenotype %in%
  #     c("MART1","CD3+41bb+")) / sum(s.t@meta.data$vasu_phenotype %in%
  #                                     c("MART1","CD3+41bb+"))
  proportion_vasu[[ident]] <- sum(
    s.t@active.ident == ident & s.t@meta.data$vasu_phenotype != "None") /
    sum(s.t@meta.data$vasu_phenotype != "None")
}

df <- as.data.frame(cbind(unlist(proportion_til), unlist(proportion_vasu)))
colnames(df) <- c("proportion_of_tils", "proportion_of_vasu")
df$cluster <- rownames(df)
df <- reshape2::melt(df)

df <- merge(df, res, by.x = "cluster", by.y = "row.names")
clusters <- unique(df$cluster)
df$is_max <- F
for (cluster in clusters){
  df$is_max[df$cluster == cluster][
    which.max(df$value[df$cluster == cluster])] <- T
}

vertical_xlabels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                     hjust=1))

DimPlot(s.t, group.by = "tcr_matches_til", order = T) +
  DimPlot(s.t, group.by = "vasu_phenotype", order = T,
          cols = c("gray","blue","red","green")) +
  ggplot(df, aes(x = cluster, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  vertical_xlabels + xlab(NULL) + ylab("Proportion") +
  geom_text(aes(label = ifelse(q < 0.05, "*", "")),
            data = subset(df, is_max), nudge_y = 0.02) +
  ggtitle("Proportion of either that is present in a given cluster")
ggsave(filename = paste0(outdir,"proportions_til_vs_vasu_matches_per_cluster.pdf"),
       width = 20, height = 7)

DimPlot(s.t, group.by = "tcr_matches_til", order = T, split.by = "UM.ID", pt.size = 1,ncol=7)
ggsave(filename = paste0(outdir,"til_matches_biopsy.per_sample.pdf"),
       width = 23, height = 7)

# Highlight those of Vasu's TCRs that match relevant things in TCR databases -----
get_cells_matching_mart1_vasu_and_db <- function(dat.tumor_til.vdj, db_matches.vasu, sname_1){

  um_id <- str_split_fixed(sname_1,"_",2)[,1]

  idx_1 <- unlist(lapply(dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3, function(x){
    all(x != "None")
  }))

  matches_cdr3 <- lapply(dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3[which(idx_1)], function(x){
    if (any(x %in% db_matches.vasu[[um_id]]$takara.cdr3)){
      antigens <- list()
      for (cdr3 in intersect(x,db_matches.vasu[[um_id]]$takara.cdr3)){
        antigens[[cdr3]] <- db_matches.vasu[[um_id]][
          db_matches.vasu[[um_id]]$takara.cdr3 == cdr3,]
        antigens[[cdr3]]$tenX.cdr3 <- cdr3
      }
      antigens <- do.call("rbind", antigens)
      rownames(antigens) <- NULL
    } else {
      antigens <- ""
    }
    return(antigens)
  })

  matches_cdr3 <- matches_cdr3[which(!unlist(lapply(matches_cdr3, function(x){
    all(x == "")
  })))]
  if (length(matches_cdr3) > 0){
    for (cell in names(matches_cdr3)){
      matches_cdr3[[cell]]$cell <- paste0(sname_1, "_", cell)
    }
    matches_cdr3 <- do.call("rbind", matches_cdr3)
    rownames(matches_cdr3) <- NULL
    matches_cdr3$db.Antigen[is.na(matches_cdr3$db.Antigen)] <- "Unknown"
    matches_cdr3 <- matches_cdr3[matches_cdr3$db.Antigen != "Unknown",]

    matches_cdr3.small <- unique(
      matches_cdr3[,c("tenX.cdr3", "db.Antigen",
                      "db.Antigen_species_or_pathology", "cell")])
    matches_cdr3.small <- matches_cdr3.small[
      matches_cdr3.small$db.Antigen == "MART1",]

    cells_matching_mart1_db <- matches_cdr3.small$cell
  } else {
    cells_matching_mart1_db <- c()
  }
  return(cells_matching_mart1_db)
}

sheets <- excel_sheets("~/proj/um_ss/Investigations/db_matches.vasu.xlsx")
db_matches.vasu <- list()
for (sheet in sheets){
  db_matches.vasu[[sheet]] <- as.data.frame(
    read_excel(path = "~/proj/um_ss/Investigations/db_matches.vasu.xlsx",
               sheet = sheet), stringsAsFactors = F)
}

gc()

cells_matching_mart1_vasu_and_db <- list()
for (sname in names(dat.tumor_til.vdj)){
  cells_matching_mart1_vasu_and_db[[sname]] <- get_cells_matching_mart1_vasu_and_db(
    dat.tumor_til.vdj,db_matches.vasu,sname)
}

cells_matching_mart1_vasu_and_db <- unlist(cells_matching_mart1_vasu_and_db)
names(cells_matching_mart1_vasu_and_db) <- NULL

cells_matching_mart1_vasu_and_db <- setNames(
  rep("MART1", length(cells_matching_mart1_vasu_and_db)),
  cells_matching_mart1_vasu_and_db)

cells_matching_mart1_vasu_and_db.biopsy <- cells_matching_mart1_vasu_and_db[
  grep("biopsy",names(cells_matching_mart1_vasu_and_db))]
names(cells_matching_mart1_vasu_and_db.biopsy) <- str_split_fixed(
  names(cells_matching_mart1_vasu_and_db.biopsy),"_biopsy_",2)[,2]

s.t <- AddMetaData(s.t, metadata = cells_matching_mart1_vasu_and_db.biopsy,
                   col.name = "matches_mart1_vasu_and_db")

tmp <- setNames(as.character(s.t@meta.data$vasu_phenotype), rownames(s.t@meta.data))
tmp[which(s.t@meta.data$vasu_phenotype == "MART1" &
            s.t@meta.data$matches_mart1_vasu_and_db == "MART1")] <- "MART1_db"
tmp[which(s.t@meta.data$vasu_phenotype == "CD3_41bb" &
            s.t@meta.data$matches_mart1_vasu_and_db == "MART1")] <- "CD3_41bb:MART1_db"
tmp[which(s.t@meta.data$vasu_phenotype == "CD3_41bb:MART1" &
            s.t@meta.data$matches_mart1_vasu_and_db == "MART1")] <- "CD3_41bb:MART1:MART1_db"
tmp <- factor(tmp,levels = c("None","CD3_41bb","MART1","MART1_db","CD3_41bb:MART1"))

s.t <- AddMetaData(s.t, metadata = tmp, col.name = "vasu_phenotype_with_mart1_db")

pdf(file = paste0(outdir,"dimplot.vasu_phenotype_with_mart1_db.pdf"), width = 7,
    height = 5)
DimPlot(s.t, group.by = "vasu_phenotype_with_mart1_db", order = T,
        cols = c("gray", "blue", "red", "green", "purple"))
dev.off()

s.t@meta.data$UM.ID <- factor(s.t@meta.data$UM.ID,
                              levels = unique(paste0("UM",sort(as.numeric(
                                gsub("UM","",s.t@meta.data$UM.ID))))))

DimPlot(s.t, group.by = "vasu_phenotype_with_mart1_db", order = T,
        cols = c("gray", "blue", "red", "green", "purple"), split.by = "UM.ID", pt.size = 1,ncol=7)
ggsave(filename = paste0(outdir,"dimplot.vasu_phenotype_with_mart1_db.per_sample.pdf"),
       width = 23, height = 7)

tils.unintegrated <- lapply(tils.unintegrated, function(x){
  if (any(names(cells_matching_mart1_vasu_and_db) %in% rownames(x@meta.data))){
    x <- AddMetaData(x, metadata = cells_matching_mart1_vasu_and_db,
                     col.name = "matches_mart1_vasu_and_db")
  } else {
    x@meta.data$cells_matching_mart1_vasu_and_db <- NA
  }
  return(x)
})

for (sname in names(tils.unintegrated)){
  x <- tils.unintegrated[[sname]]
  tmp <- setNames(as.character(x@meta.data$vasu_phenotype), rownames(x@meta.data))
  tmp[which(x@meta.data$vasu_phenotype == "MART1" &
              x@meta.data$matches_mart1_vasu_and_db == "MART1")] <- "MART1_db"
  tmp[which(x@meta.data$vasu_phenotype == "CD3_41bb" &
              x@meta.data$matches_mart1_vasu_and_db == "MART1")] <- "CD3_41bb:MART1_db"
  tmp[which(x@meta.data$vasu_phenotype == "CD3_41bb:MART1" &
              x@meta.data$matches_mart1_vasu_and_db == "MART1")] <- "CD3_41bb:MART1:MART1_db"
  tmp <- factor(tmp,levels = c("None","CD3_41bb","MART1","MART1_db","CD3_41bb:MART1"))


  tils.unintegrated[[sname]] <- AddMetaData(
    x, metadata = tmp, col.name = "vasu_phenotype_with_mart1_db")
}

p_um1 <- DimPlot(tils.unintegrated[["UM1"]],
                 group.by = "vasu_phenotype_with_mart1_db", order = T, pt.size = 1,
                 cols =c("gray", "blue", "red", "green", "purple")) + ggtitle("UM1")
p_um2 <- DimPlot(tils.unintegrated[["UM2"]],
                 group.by = "vasu_phenotype_with_mart1_db", order = T, pt.size = 1,
                 cols = c("gray", "blue", "red", "green", "purple")) + ggtitle("UM2")
p_um3 <- DimPlot(tils.unintegrated[["UM3"]],
                 group.by = "vasu_phenotype_with_mart1_db", order = T, pt.size = 1,
                 cols = c("gray", "blue", "red", "green", "purple")) + ggtitle("UM3")
p_um9 <- DimPlot(tils.unintegrated[["UM9"]],
                 group.by = "vasu_phenotype_with_mart1_db", order = T, pt.size = 1,
                 cols = c("gray", "blue", "red", "green", "purple")) + ggtitle("UM9")
p_um10 <- DimPlot(tils.unintegrated[["UM10"]],
                  group.by = "vasu_phenotype_with_mart1_db", order = T, pt.size = 1,
                  cols = c("gray", "blue", "red", "green", "purple")) + ggtitle("UM10")
p_um19 <- DimPlot(tils.unintegrated[["UM19"]],
                  group.by = "vasu_phenotype_with_mart1_db", order = T, pt.size = 1,
                  cols = c("gray", "blue", "red", "green", "purple")) + ggtitle("UM19")
p_um22 <- DimPlot(tils.unintegrated[["UM22"]],
                  group.by = "vasu_phenotype_with_mart1_db", order = T, pt.size = 1,
                  cols = c("gray", "blue", "red", "green", "purple")) + ggtitle("UM22")
p_um24 <- DimPlot(tils.unintegrated[["UM24"]],
                  group.by = "vasu_phenotype_with_mart1_db", order = T, pt.size = 1,
                  cols = c("gray", "blue", "red", "green", "purple")) + ggtitle("UM24")

library(gridExtra)
pdf(file = paste0(outdir,"til_vasu_phenotype_with_mart1_db.pdf"),width = 55, height = 6)
p <- grid.arrange(grobs = c(list(p_um1),
                            list(p_um2),
                            list(p_um3),
                            list(p_um9),
                            list(p_um19),
                            list(p_um10),
                            list(p_um24),
                            list(p_um22)), ncol = 8, as.table = FALSE)
dev.off()

#saveRDS(dat.tumor,file=paste0(outdir,"dat.tumor.rda"))
#saveRDS(tils.unintegrated,file=paste0(outdir,"tils.unintegrated.rda"))
