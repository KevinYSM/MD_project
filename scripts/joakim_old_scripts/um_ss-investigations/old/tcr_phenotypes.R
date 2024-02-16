library(rliger)
library(Seurat)
library(tidyverse)
library(readxl)

outdir <- "~/proj/um_ss/Investigations/seurat/results/liger_all/"
lig <- readRDS(paste0(outdir,"lig.blood_tumor_tils.integrated.rda"))
dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")

# Read Vasu clonotypes ---------------------------------------------------------
fnames <- Sys.glob("~/proj/um_ss/Pipelines/takara_smarter/preprocessing/results_both/report/results_both_*_mig_cdr3_report.xlsx")
vdj.takara <- list()
for (fname in fnames){
  sname <- str_split_fixed(basename(fname),"_",4)[,3]
  tmp_A <- as.data.frame(read_excel(fname,sheet = paste0(sname,"_TRA_clone")),stringsAsFactors=F)
  tmp_B <- as.data.frame(read_excel(fname,sheet = paste0(sname,"_TRB_clone")),stringsAsFactors=F)

  tmp_A$Sample <- sname
  tmp_B$Sample <- sname
  
  tmp_A$Chain <- "TRA"
  tmp_B$Chain <- "TRB"
  
  vdj.takara[[basename(fname)]] <- unique(rbind(tmp_A[,c("Sample","Chain",
                                                     "CDR3 Amino Acid Sequence",
                                                     "V segment","D segment",
                                                     "J segment","C segment")],
                                            tmp_B[,c("Sample","Chain",
                                                     "CDR3 Amino Acid Sequence",
                                                     "V segment","D segment",
                                                     "J segment","C segment")]))
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
  vdj.takara[,c("chains","cdr3","v_gene","d_gene","j_gene","c_gene")],1,paste0,
  collapse = ":")

# Find clonotypes and samples shared with Vasu's data --------------------------

cnames <- c("Sample.ID","project.name","UM.ID",
            "v_gene","d_gene","j_gene","c_gene","chains","cdr3",
            "reads","umis","productive","djvdj.clone_freq","djvdj.clone_pct")

vdj.10x <- dat.cca@meta.data[,cnames]
vdj.10x$cell_id <- rownames(vdj.10x)

head(vdj.10x)
head(vdj.takara)

gsub("None","",vdj.10x$d_gene)

clean_fun <- function(x){
  res <- unlist(strsplit(x,";"))
  res <- res[!is.na(res)]
  res <- setdiff(res,"None")
  if (length(res)==0){
    res <- "None"
  }
  return(res)
}

vdj.10x_cdr3 <- setNames(sapply(vdj.10x$cdr3,clean_fun),rownames(vdj.10x))
vdj.10x_v_gene <- setNames(sapply(vdj.10x$v_gene,clean_fun),rownames(vdj.10x))
vdj.10x_d_gene <- setNames(sapply(vdj.10x$d_gene,clean_fun),rownames(vdj.10x))
vdj.10x_j_gene <- setNames(sapply(vdj.10x$j_gene,clean_fun),rownames(vdj.10x))
vdj.10x_c_gene <- setNames(sapply(vdj.10x$c_gene,clean_fun),rownames(vdj.10x))
vdj.10x_chains <- setNames(sapply(vdj.10x$chains,clean_fun),rownames(vdj.10x))

vdj.10x_v_gene.match_idx <- sapply(vdj.10x_v_gene,function(x){which(vdj.takara$v_gene %in% x)})
vdj.10x_d_gene.match_idx <- sapply(vdj.10x_d_gene,function(x){which(vdj.takara$d_gene %in% x)})
vdj.10x_j_gene.match_idx <- sapply(vdj.10x_j_gene,function(x){which(vdj.takara$j_gene %in% x)})
vdj.10x_c_gene.match_idx <- sapply(vdj.10x_c_gene,function(x){which(vdj.takara$c_gene %in% x)})
vdj.10x_cdr3.match_idx <- sapply(vdj.10x_cdr3,function(x){which(vdj.takara$cdr3 %in% x)})

vdj.10x.cdr3.match <- sapply(vdj.10x_cdr3.match_idx,length) > 0
vdj.10x.v_gene.match <- rep(F,length(vdj.10x.cdr3.match))
vdj.10x.d_gene.match <- rep(F,length(vdj.10x.cdr3.match))
vdj.10x.j_gene.match <- rep(F,length(vdj.10x.cdr3.match))
vdj.10x.c_gene.match <- rep(F,length(vdj.10x.cdr3.match))
widx <- which(vdj.10x.cdr3.match)
res <- lapply(widx,function(i){
  v_gene.match <- intersect(vdj.10x_cdr3.match_idx[[i]],vdj.10x_v_gene.match_idx[[i]])
  d_gene.match <- intersect(vdj.10x_cdr3.match_idx[[i]],vdj.10x_d_gene.match_idx[[i]])
  j_gene.match <- intersect(vdj.10x_cdr3.match_idx[[i]],vdj.10x_j_gene.match_idx[[i]])
  c_gene.match <- intersect(vdj.10x_cdr3.match_idx[[i]],vdj.10x_c_gene.match_idx[[i]])
  
  all_match <- as.numeric(names(which(table(c(v_gene.match,d_gene.match,j_gene.match,c_gene.match))==4)))
  
  return(
    list(vdj.10x.v_gene.match=any(v_gene.match),
         vdj.10x.d_gene.match=any(d_gene.match),
         vdj.10x.j_gene.match=any(j_gene.match),
         vdj.10x.c_gene.match=any(c_gene.match),
         v_gene.match=paste0(unlist(v_gene.match),collapse = ";"),
         d_gene.match=paste0(unlist(d_gene.match),collapse = ";"),
         j_gene.match=paste0(unlist(j_gene.match),collapse = ";"),
         c_gene.match=paste0(unlist(c_gene.match),collapse = ";"),
         all_match=paste0(all_match,collapse=";"))
  )
})
res <- as.data.frame(do.call("rbind",res))
vdj.10x.v_gene.match[widx] <- unlist(res$vdj.10x.v_gene.match)
vdj.10x.d_gene.match[widx] <- unlist(res$vdj.10x.d_gene.match)
vdj.10x.j_gene.match[widx] <- unlist(res$vdj.10x.j_gene.match)
vdj.10x.c_gene.match[widx] <- unlist(res$vdj.10x.c_gene.match)

vdj.10x.match <- data.frame(vdj.10x.cdr3.match,
           vdj.10x.v_gene.match,
           vdj.10x.d_gene.match,
           vdj.10x.j_gene.match,
           vdj.10x.c_gene.match)

stopifnot(identical(
  rownames(vdj.10x),
  rownames(vdj.10x.match)
))

vdj.10x$vdj.10x.cdr3.match <- vdj.10x.cdr3.match
vdj.10x$vdj.10x.v_gene.match <- vdj.10x.v_gene.match
vdj.10x$vdj.10x.d_gene.match <- vdj.10x.d_gene.match
vdj.10x$vdj.10x.j_gene.match <- vdj.10x.j_gene.match
vdj.10x$vdj.10x.c_gene.match <- vdj.10x.c_gene.match

vdj.10x$vdj.10x.v_gene.match <- vdj.10x.v_gene.match
vdj.10x$vdj.10x.d_gene.match <- vdj.10x.d_gene.match
vdj.10x$vdj.10x.j_gene.match <- vdj.10x.j_gene.match
vdj.10x$vdj.10x.c_gene.match <- vdj.10x.c_gene.match

vdj.10x$v_gene.match <- ""
vdj.10x$d_gene.match <- ""
vdj.10x$j_gene.match <- ""
vdj.10x$c_gene.match <- ""
vdj.10x$all_match <- ""

vdj.10x$v_gene.match[widx] <- unlist(res$v_gene.match)
vdj.10x$d_gene.match[widx] <- unlist(res$d_gene.match)
vdj.10x$j_gene.match[widx] <- unlist(res$j_gene.match)
vdj.10x$c_gene.match[widx] <- unlist(res$c_gene.match)
vdj.10x$all_match[widx] <- unlist(res$all_match)

#idx_all_match <- matrixStats::rowProds(as.matrix(sapply(vdj.10x.match,as.numeric)))>0
idx_all_match <- vdj.10x$all_match != ""
vdj.10x.match[idx_all_match,]

table(vdj.10x[rownames(vdj.10x.match[idx_all_match,]),]$project.name)
table(vdj.10x[rownames(vdj.10x.match[idx_all_match,]),]$UM.ID)

#vdj.10x$all_match <- rownames(vdj.10x) %in% rownames(vdj.10x.match[idx_all_match,])

head(vdj.10x[vdj.10x$all_match != "",])

vdj.10x$takara_clonotype <- ""
widx <- which(vdj.10x$all_match != "")
res <- sapply(widx,function(i){
  takara_widx <- as.numeric(unlist(strsplit(vdj.10x[i,]$all_match,";")))
  res <- unique(apply(vdj.takara[takara_widx,c("chains","cdr3","v_gene","d_gene","j_gene","c_gene")],
               1,paste0,collapse = ":"))
  res <- paste0(res,collapse=";")
  return(res)
})
vdj.10x$takara_clonotype[widx] <- res

head(vdj.10x[widx,])


table(vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype)
table(vdj.takara$UM.ID[vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype])
table(vdj.takara$vasu_phenotype[vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype])
table(vdj.takara$vasu_phenotype[vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype & vdj.takara$Sample!="Undetermined"],
      vdj.takara$UM.ID[vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype & vdj.takara$Sample!="Undetermined"])

table(vdj.takara[vdj.takara$Sample!="Undetermined",]$vasu_phenotype,
      vdj.takara[vdj.takara$Sample!="Undetermined",]$UM.ID)

stopifnot(identical(rownames(vdj.10x),rownames(dat.cca@meta.data)))
dat.cca <- AddMetaData(dat.cca,
                       metadata = vdj.10x[,
                                          setdiff(colnames(vdj.10x),
                                                  colnames(dat.cca@meta.data))])

intersect(dat.cca@meta.data$takara_clonotype,vdj.takara$clonotype)

x <- merge(dat.cca@meta.data,vdj.takara,by.x="takara_clonotype",by.y="clonotype",all.x=T,all.y=F)
x[,c("takara_clonotype","vasu_phenotype","cell_id")]

vasu_phenotype <- setNames(x$vasu_phenotype,x$cell_id)
vasu_phenotype <- vasu_phenotype[rownames(dat.cca@meta.data)]
identical(names(vasu_phenotype),rownames(dat.cca@meta.data))
vasu_phenotype[is.na(vasu_phenotype)] <- "None"
vasu_phenotype[vasu_phenotype==""] <- "Unidentified"

dat.cca <- AddMetaData(dat.cca,metadata = vasu_phenotype,col.name = "vasu_phenotype")

# Within each dataset and sample, locate the Vasu clonotypes -------------------

idx <- vdj.10x$all_match != ""
vdj.10x[idx,]
vdj.10x[!idx,]

dim(vdj.10x[idx,])
dim(vdj.10x[!idx & !vdj.10x.cdr3.match, ])

dim(vdj.10x[idx & vdj.10x$project.name=="G18-023",])
dim(vdj.10x[idx & vdj.10x$project.name=="G18-049",])
dim(vdj.10x[idx & vdj.10x$project.name=="GWA-JN-388",])

table(vdj.10x[idx & vdj.10x$project.name=="G18-023",]$UM.ID)
table(vdj.10x[idx & vdj.10x$project.name=="G18-049",]$UM.ID)
table(vdj.10x[idx & vdj.10x$project.name=="G18-049",]$Sample.ID)
table(vdj.10x[idx & vdj.10x$project.name=="GWA-JN-388",]$UM.ID)
table(vdj.10x[idx & vdj.10x$project.name=="GWA-JN-388",]$Sample.ID)

tmp <- reshape2::melt(table(vdj.10x$UM.ID,vdj.10x$Sample.ID))
tmp <- tmp[tmp$value>0,]
rownames(tmp) <- NULL
tmp <- tmp[,c(1,2)]
sample_map <- tmp

sample_map$Var1 <- as.character(sample_map$Var1)
sample_map$Var2 <- as.character(sample_map$Var2)

sample_map$in_takara <- F
sample_map[sample_map$Var1 %in% setdiff(unique(vdj.takara$UM.ID),""),]$in_takara <- T
sample_map[sample_map$in_takara,]

sample_map$vasu_phenotype <- ""
sample_map$vasu_phenotype[sample_map$Var1 %in% c("UM1","UM9")] <- "CD3+41bb+"
sample_map$vasu_phenotype[sample_map$Var1 %in% c("UM22","UM46")] <- "MART1"

# Testing (cell type agnostic) -------------------------------------------------

snames <- unique(dat.cca@meta.data$Sample.ID[dat.cca@meta.data$takara_clonotype!=""])
res <- list()
res.downsample_100 <- list()
for (sname in snames){
  print(sname)
  
  idx_clones_takara <- dat.cca@meta.data$Sample.ID == sname & 
    dat.cca@meta.data$takara_clonotype != "" & dat.cca@meta.data$vasu_phenotype == "CD3+41bb+"
  idx_clones_other <- dat.cca@meta.data$Sample.ID == sname & 
    dat.cca@meta.data$takara_clonotype == "" & !is.na(dat.cca@meta.data$cdr3) & 
    dat.cca@meta.data$vasu_phenotype != "MART1" & dat.cca@meta.data$vasu_phenotype != "CD3+41bb+"
  #dat.cca@meta.data[idx_clones_takara,]$active.ident.updated_2
  
  if (sum(idx_clones_takara)>0 & sum(idx_clones_other)>0){
    grp <- rep("grp_exclude",length(idx_clones_takara))
    grp[idx_clones_takara] <- "grp_takara"
    grp[idx_clones_other] <- "grp_other"
    print(table(grp))
    names(grp) <- rownames(dat.cca@meta.data)
    
    stopifnot(identical(names(Idents(dat.cca)),names(grp)))
    Idents(dat.cca) <- grp
    
    res[[sname]] <- FindMarkers(object = dat.cca,
                                ident.1 = "grp_takara",
                                ident.2 = "grp_other",
                                test.use = "wilcox", min.cells.group = 1, 
                                assay = "RNA", slot = "data", verbose = T)
    res.downsample_100[[sname]] <- FindMarkers(object = dat.cca, 
                                               ident.1 = "grp_takara",
                                               ident.2 = "grp_other",
                                               test.use = "wilcox", 
                                               min.cells.group = 1,
                                               assay = "RNA", slot = "data", 
                                               max.cells.per.ident = 100,
                                               verbose = T)
  }
}

nm <- names(res[which(names(res) %in% sample_map[sample_map$in_takara,]$Var2)])

genes <- list()
for (n in nm){
  # genes[[n]] <- intersect(
  #   rownames(res[[n]][
  #     res[[n]]$p_val_adj < 0.05 & res[[n]]$avg_log2FC > 0,]),
  #   rownames(res.downsample_100[[n]][
  #     res.downsample_100[[n]]$p_val_adj < 0.05 & res.downsample_100[[n]]$avg_log2FC > 0,]))
  genes[[n]] <- rownames(res[[n]][res[[n]]$p_val_adj < 0.05 & res[[n]]$avg_log2FC > 10,])
}

sort(names(table(unlist(genes))[table(unlist(genes))>=2]))

# Testing (matched cell cycle phase) -------------------------------------------
DefaultAssay(dat.cca) <- "RNA"
dat.cca <- CellCycleScoring(object = dat.cca,
                            s.features = cc.genes$s.genes, 
                            g2m.features = cc.genes$g2m.genes, set.ident = FALSE)

snames <- unique(dat.cca@meta.data$Sample.ID[dat.cca@meta.data$takara_clonotype!=""])
res <- list()
res.downsample_100 <- list()
for (sname in snames){
  print(sname)
  idx_clones_takara <- dat.cca@meta.data$Sample.ID == sname & 
    dat.cca@meta.data$takara_clonotype != ""
  idx_clones_other <- dat.cca@meta.data$Sample.ID == sname & 
    dat.cca@meta.data$takara_clonotype == "" & !is.na(dat.cca@meta.data$cdr3)
  #dat.cca@meta.data[idx_clones_takara,]$active.ident.updated_2
  
  phases_takara <- unique(dat.cca@meta.data$Phase[idx_clones_takara])
  phases_other <- unique(dat.cca@meta.data$Phase[idx_clones_other])
  
  phases_compare <- intersect(phases_takara,phases_other)
  if (length(phases_compare)>0){
    
    for (phase in phases_compare){
      grp <- rep("grp_exclude",length(idx_clones_takara))
      grp[idx_clones_takara & dat.cca@meta.data$Phase == phase] <- "grp_takara"
      grp[idx_clones_other & dat.cca@meta.data$Phase == phase] <- "grp_other"
      print(table(grp))
      
      names(grp) <- rownames(dat.cca@meta.data)
      
      stopifnot(identical(names(Idents(dat.cca)),names(grp)))
      Idents(dat.cca) <- grp
      
      res[[paste0(sname,":",phase)]] <- FindMarkers(object = dat.cca,
                                  ident.1 = "grp_takara",
                                  ident.2 = "grp_other",
                                  test.use = "wilcox", min.cells.group = 1, 
                                  assay = "RNA", slot = "data", verbose = T)
      res.downsample_100[[paste0(sname,":",phase)]] <- FindMarkers(object = dat.cca, 
                                                 ident.1 = "grp_takara",
                                                 ident.2 = "grp_other",
                                                 test.use = "wilcox", 
                                                 min.cells.group = 1,
                                                 assay = "RNA", slot = "data", 
                                                 max.cells.per.ident = 100,
                                                 verbose = T)
    }
  } else {
    warning("No matching cell cycle phases")
  }
}


nm <- names(res[which(str_split_fixed(names(res),":",2)[,1] %in% sample_map[sample_map$in_takara,]$Var2)])
snames <- str_split_fixed(nm,":",2)[,1]
phases <- str_split_fixed(nm,":",2)[,2]

genes <- list()
for (sname in snames){
  phs <- phases[snames==sname]
  genes_phases <- list()
  for (phase in phs){
    n <- paste0(sname,":",phase)
    genes_phases[[n]] <- sort(rownames(res[[n]][res[[n]]$p_val_adj < 0.05 & res[[n]]$avg_log2FC > 0,]))
  }
  if (length(phs)>1){
    genes[[sname]] <- Reduce("intersect",genes_phases)
  } else {
    genes[[sname]] <- genes_phases[[1]]
  }
}

genes

sort(names(table(unlist(genes))[table(unlist(genes))>=2]))


# Testing (matched cell cycle phase, restrict to CD8+) -------------------------------------------
DefaultAssay(dat.cca) <- "RNA"
dat.cca <- CellCycleScoring(object = dat.cca,
                            s.features = cc.genes$s.genes, 
                            g2m.features = cc.genes$g2m.genes, set.ident = FALSE)

snames <- unique(dat.cca@meta.data$Sample.ID[dat.cca@meta.data$takara_clonotype!=""])
res <- list()
res.downsample_100 <- list()
for (sname in snames){
  print(sname)
  
  idx_cd8 <- (dat.cca@assays$RNA@counts["CD8A",] > 0 | dat.cca@assays$RNA@counts["CD8B",] > 0) &
    dat.cca@assays$RNA@counts["CD4",] == 0
  
  idx_clones_takara <- dat.cca@meta.data$Sample.ID == sname & 
    dat.cca@meta.data$takara_clonotype != "" & idx_cd8
  idx_clones_other <- dat.cca@meta.data$Sample.ID == sname & 
    dat.cca@meta.data$takara_clonotype == "" & !is.na(dat.cca@meta.data$cdr3) & idx_cd8
  #dat.cca@meta.data[idx_clones_takara,]$active.ident.updated_2
  
  if (sum(idx_clones_takara)>0 & sum(idx_clones_other)>0){
    
    phases_takara <- unique(dat.cca@meta.data$Phase[idx_clones_takara])
    phases_other <- unique(dat.cca@meta.data$Phase[idx_clones_other])
    
    phases_compare <- intersect(phases_takara,phases_other)
    if (length(phases_compare)>0){
      
      for (phase in phases_compare){
        grp <- rep("grp_exclude",length(idx_clones_takara))
        grp[idx_clones_takara & dat.cca@meta.data$Phase == phase] <- "grp_takara"
        grp[idx_clones_other & dat.cca@meta.data$Phase == phase] <- "grp_other"
        print(table(grp))
        
        names(grp) <- rownames(dat.cca@meta.data)
        
        stopifnot(identical(names(Idents(dat.cca)),names(grp)))
        Idents(dat.cca) <- grp
        
        res[[paste0(sname,":",phase)]] <- FindMarkers(object = dat.cca,
                                                      ident.1 = "grp_takara",
                                                      ident.2 = "grp_other",
                                                      test.use = "wilcox", min.cells.group = 1, 
                                                      assay = "RNA", slot = "data", verbose = T)
        res.downsample_100[[paste0(sname,":",phase)]] <- FindMarkers(object = dat.cca, 
                                                                     ident.1 = "grp_takara",
                                                                     ident.2 = "grp_other",
                                                                     test.use = "wilcox", 
                                                                     min.cells.group = 1,
                                                                     assay = "RNA", slot = "data", 
                                                                     max.cells.per.ident = 100,
                                                                     verbose = T)
      }
    } else {
      warning("No matching cell cycle phases")
    }
  }
}


nm <- names(res[which(str_split_fixed(names(res),":",2)[,1] %in% sample_map[sample_map$in_takara,]$Var2)])
snames <- str_split_fixed(nm,":",2)[,1]
phases <- str_split_fixed(nm,":",2)[,2]

genes <- list()
for (sname in snames){
  phs <- phases[snames==sname]
  genes_phases <- list()
  for (phase in phs){
    n <- paste0(sname,":",phase)
    genes_phases[[n]] <- sort(rownames(res[[n]][res[[n]]$p_val_adj < 0.05 & res[[n]]$avg_log2FC > 0,]))
  }
  if (length(phs)>1){
    #genes[[sname]] <- Reduce("intersect",genes_phases)
    genes[[sname]] <- Reduce("union",genes_phases)
  } else {
    genes[[sname]] <- genes_phases[[1]]
  }
}

genes

sort(names(table(unlist(genes))[table(unlist(genes))>=3]))

# Per cluster ------------------------------------------------------------------
tils <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/tils/tils.rda")

sample_map[sample_map$in_takara,]$Var2

tils

# Exploratory ------------------------------------------------------------------
sample_map[sample_map$in_takara,]

dat.cca.split <- SplitObject(dat.cca, split.by = "Sample.ID")

get_data <- function(x,sname){
  x <- dat.cca.split[[sname]]
  x <- FindVariableFeatures(x)
  x <- ScaleData(x)
  x <- RunPCA(x,npcs = 30)
  x <- RunUMAP(x,dims = 1:30)
  x <- FindNeighbors(x,dims=1:30)
  x <- FindClusters(x)
  return(x)
}

# x <- dat.cca.split[["C7-GEX"]]
# x <- FindVariableFeatures(x)
# x <- ScaleData(x)
# x <- RunPCA(x,npcs = 30)
# x <- RunUMAP(x,dims = 1:30)
# x <- FindNeighbors(x,dims=1:30)
# x <- FindClusters(x)
# DimPlot(x)

#x <- dat.cca[["SampleID_2_11june18"]]
#cnames <- setdiff(colnames(dat.cca@meta.data),colnames(x@meta.data))
#x <- AddMetaData(x,metadata = dat.cca@meta.data[colnames(x),cnames])
##unique(vdj.takara[vdj.takara$vasu_phenotype=="",]$Sample)

#DimPlot(x,group.by = "vasu_phenotype",split.by = "active.ident.updated_2")
#DimPlot(x,group.by = "active.ident.updated_2",split.by = "vasu_phenotype")
DimPlot(x,group.by = "seurat_clusters",split.by = "vasu_phenotype")
DimPlot(x,group.by = "takara_clonotype",split.by = "vasu_phenotype")

unique(x@meta.data$vasu_phenotype)

run_tests <- function(x, sname){
  x <- get_data(x = dat.cca.split, sname = sname)
  
  res <- list()
  res.groups <- list()
  for (pheno in setdiff(unique(x@meta.data$vasu_phenotype),c("None","Unidentified"))){
    
    idx_cd8 <- (x@assays$RNA@counts["CD8A",] > 0 | x@assays$RNA@counts["CD8B",] > 0) & 
      x@assays$RNA@counts["CD4",] == 0
    
    idx_clones_takara <- x@meta.data$takara_clonotype != "" & x@meta.data$vasu_phenotype == "CD3+41bb+" & idx_cd8
    idx_clones_other <- x@meta.data$takara_clonotype == "" & !is.na(x@meta.data$cdr3) & x@meta.data$vasu_phenotype == "None" & idx_cd8
    
    if (sum(idx_clones_takara)>0 & sum(idx_clones_other)>0){
    
      phases_takara <- unique(x@meta.data$Phase[idx_clones_takara])
      phases_other <- unique(x@meta.data$Phase[idx_clones_other])
      phases_shared <- intersect(phases_takara,phases_other)
      
      if (length(phases_shared)>0){
        for (phase in phases_shared){
          
          clusters_takara <- unique(as.character(sort(x@meta.data$seurat_clusters[idx_clones_takara & x@meta.data$Phase == phase])))
          clusters_other <- unique(as.character(sort(x@meta.data$seurat_clusters[idx_clones_other & x@meta.data$Phase == phase])))
          clusters_shared <- intersect(clusters_takara,clusters_other)
          
          if (length(clusters_shared)>0){
            for (cluster in clusters_shared){
              grp <- rep("grp_exclude",length(idx_clones_takara))
              grp[idx_clones_takara & x@meta.data$Phase == phase & 
                    as.character(x@meta.data$seurat_clusters) == cluster] <- "grp_takara"
              grp[idx_clones_other & x@meta.data$Phase == phase & 
                    as.character(x@meta.data$seurat_clusters) == cluster] <- "grp_other"
              names(grp) <- rownames(x@meta.data)
              stopifnot(identical(names(Idents(x)),names(grp)))
              Idents(x) <- grp
              
              id <- paste0(pheno,":",phase,":",cluster)
              
              res.groups[[id]] <- table(Idents(x))
              print(res.groups[[id]])
              
              res[[id]] <- FindMarkers(object = x,
                                 ident.1 = "grp_takara",
                                 ident.2 = "grp_other",
                                 test.use = "wilcox", min.cells.group = 1, 
                                 assay = "RNA", slot = "data", verbose = T)
            }
          }
        }
      }
    }
  }
  return(list(res = res,
              res.groups = res.groups))
}

snames <- sample_map[sample_map$in_takara,]$Var2

res_all <- lapply(snames,function(sname){
  res <- run_tests(x = dat.cca.split, sname = sname)  
})
names(res_all) <- snames


filter_res <- function(res){
  #res <- res[res$p_val_adj<0.05 & res$avg_log2FC > 0 & res$pct.1 > 0,]
  res <- res[res$p_val<0.05 & res$avg_log2FC > 0 & res$pct.1 > 0,]
  res$gene <- rownames(res)
  return(res)
}

lapply(res,filter_res)

lapply(res_all$SampleID_7_11june18$res,filter_res)

do.call("rbind",lapply(res_all$SampleID_7_11june18$res,filter_res))

res_all_filtered <- sapply(names(res_all), function(sname){
  tmp <- do.call("rbind",lapply(res_all[[sname]]$res,filter_res))
  tmp$phenotype <- str_split_fixed(rownames(tmp),":",3)[,1]
  tmp$phase <- str_split_fixed(rownames(tmp),":",3)[,2]
  tmp$cluster <- str_split_fixed(str_split_fixed(rownames(tmp),":",3)[,3],"\\.",2)[,1]
  rownames(tmp) <- NULL
  tmp$sample <- sname
  return(tmp)
})

res_all_filtered <- res_all_filtered[which(unlist(lapply(res_all_filtered,class))=="data.frame")]

lapply(res_all_filtered,function(x){sort(table(x$gene),decreasing = T)})

sort(table(unlist(lapply(res_all_filtered,function(x){x$gene}))),decreasing = T)

outdir <- "~/proj/um_ss/Investigations/seurat/results/"
saveRDS(list(res_all = res_all,
     res_all_filtered = res_all_filtered),
     file = paste0(outdir,"tcr_clonotypes_diff_expr.phase_cluster_cd8.rda"))

# Using cell types instead of clusters -----------------------------------------

run_tests2 <- function(x, sname){
  x <- dat.cca.split[[sname]]
  
  res <- list()
  res.groups <- list()
  for (pheno in setdiff(unique(x@meta.data$vasu_phenotype),c("None","Unidentified"))){
    
    idx_cd8 <- (x@assays$RNA@counts["CD8A",] > 0 | x@assays$RNA@counts["CD8B",] > 0) & 
      x@assays$RNA@counts["CD4",] == 0
    
    idx_clones_takara <- x@meta.data$takara_clonotype != "" & x@meta.data$vasu_phenotype == "CD3+41bb+" & idx_cd8
    idx_clones_other <- x@meta.data$takara_clonotype == "" & !is.na(x@meta.data$cdr3) & x@meta.data$vasu_phenotype == "None" & idx_cd8
    
    if (sum(idx_clones_takara)>0 & sum(idx_clones_other)>0){
      
      phases_takara <- unique(x@meta.data$Phase[idx_clones_takara])
      phases_other <- unique(x@meta.data$Phase[idx_clones_other])
      phases_shared <- intersect(phases_takara,phases_other)
      
      if (length(phases_shared)>0){
        for (phase in phases_shared){
          
          clusters_takara <- unique(as.character(sort(x@meta.data$active.ident.updated_2[idx_clones_takara & x@meta.data$Phase == phase])))
          clusters_other <- unique(as.character(sort(x@meta.data$active.ident.updated_2[idx_clones_other & x@meta.data$Phase == phase])))
          clusters_shared <- intersect(clusters_takara,clusters_other)
          
          if (length(clusters_shared)>0){
            for (cluster in clusters_shared){
              grp <- rep("grp_exclude",length(idx_clones_takara))
              grp[idx_clones_takara & x@meta.data$Phase == phase & 
                    as.character(x@meta.data$active.ident.updated_2) == cluster] <- "grp_takara"
              grp[idx_clones_other & x@meta.data$Phase == phase & 
                    as.character(x@meta.data$active.ident.updated_2) == cluster] <- "grp_other"
              names(grp) <- rownames(x@meta.data)
              stopifnot(identical(names(Idents(x)),names(grp)))
              Idents(x) <- grp
              
              id <- paste0(pheno,":",phase,":",cluster)
              
              res.groups[[id]] <- table(Idents(x))
              print(res.groups[[id]])
              
              res[[id]] <- FindMarkers(object = x,
                                       ident.1 = "grp_takara",
                                       ident.2 = "grp_other",
                                       test.use = "wilcox", min.cells.group = 1, 
                                       assay = "RNA", slot = "data", verbose = T)
            }
          }
        }
      }
    }
  }
  return(list(res = res,
              res.groups = res.groups))
}

snames <- sample_map[sample_map$in_takara,]$Var2

res_all <- lapply(snames,function(sname){
  res <- run_tests2(x = dat.cca.split, sname = sname)  
})
names(res_all) <- snames

sapply(names(res_all), function(sname){
  lapply(res_all[[sname]]$res,filter_res)
})

res_all_filtered <- sapply(names(res_all), function(sname){
  if (length(res_all[[sname]]$res)>0){
    tmp <- do.call("rbind",lapply(res_all[[sname]]$res,filter_res))
    tmp$phenotype <- str_split_fixed(rownames(tmp),":",3)[,1]
    tmp$phase <- str_split_fixed(rownames(tmp),":",3)[,2]
    tmp$cluster <- str_split_fixed(str_split_fixed(rownames(tmp),":",3)[,3],"\\.",2)[,1]
    rownames(tmp) <- NULL
    if(nrow(tmp)>0){
      tmp$sample <- sname
    } else {
      tmp$sample <- character(0)
    }
  } else {
    tmp <- NULL
  }
  return(tmp)
})

res_all_filtered <- res_all_filtered[which(unlist(lapply(res_all_filtered,class))=="data.frame")]

lapply(res_all_filtered,function(x){sort(table(x$gene),decreasing = T)})

sort(table(unlist(lapply(res_all_filtered,function(x){x$gene}))),decreasing = T)

outdir <- "~/proj/um_ss/Investigations/seurat/results/"
saveRDS(list(res_all = res_all,
             res_all_filtered = res_all_filtered),
        file = paste0(outdir,"tcr_clonotypes_diff_expr.phase_cell_types_cd8.rda"))

# Interactive exploration ------------------------------------------------------

sapply(snames,function(x){
  table(dat.cca.split[[x]]@meta.data$vasu_phenotype)
})

sapply(names(dat.cca.split),function(x){
  table(dat.cca.split[[x]]@meta.data$vasu_phenotype)
})

table(dat.cca@meta.data$vasu_phenotype,dat.cca@meta.data$Sample.ID)

table(dat.cca.split[["A2-GEX"]]@meta.data$vasu_phenotype)

snames

#sname <- "SampleID_9_11june18"
sname <- "A2-GEX"
x <- dat.cca.split[[sname]]

DimPlot(x,group.by = "active.ident.updated_2")
DimPlot(x,group.by = "vasu_phenotype")

table(x@meta.data$vasu_phenotype,x@meta.data$active.ident.updated_2)
table(x@meta.data$Phase,x@meta.data$active.ident.updated_2)
table(x@meta.data$Phase,x@meta.data$vasu_phenotype)


dat <- log2(
  x@assays$RNA@data[
    ,
    x@meta.data$active.ident.updated_2=="CD8 T cells (I)"
  ]+1)

pheatmap(cor(as.matrix(dat)),show_rownames = F,show_colnames = F,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])

dat <- t(scale(t(log2(
  x@assays$RNA@data[
    genes,
    x@meta.data$active.ident.updated_2=="CD8 T cells (I)"
  ]+1)),center = T,scale = T))

dat <- dat[apply(dat,1,function(x){all(is.finite(x))}),]
dim(dat)

pheatmap(dat[,
             c(intersect(colnames(dat),rownames(x@meta.data)[x@meta.data$vasu_phenotype=="MART1"]),
               intersect(colnames(dat),rownames(x@meta.data)[x@meta.data$vasu_phenotype=="None"]))],
         show_colnames = F,cluster_cols = F,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])

pheatmap(dat[,
             c(intersect(colnames(dat),rownames(x@meta.data)[x@meta.data$vasu_phenotype=="MART1"]),
               intersect(colnames(dat),rownames(x@meta.data)[x@meta.data$vasu_phenotype=="None"]))],
         show_colnames = F,cluster_cols = T,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])



dat <- x@assays$RNA@scale.data[,grepl("^CD8 T cells",x@meta.data$active.ident.updated_2)]
#dat <- dat[rowSums(dat>0)>=100,]
dim(dat)

cv <- apply(dat,1,sd)/rowMeans(dat)
max(cv,na.rm = T)
hist(cv)
idx <- which(cv >= quantile(cv,na.rm = T,probs = seq(0,1,0.01))[["99%"]])
#idx <- which(cv >= min(sort(cv,decreasing = T)[1:500]))
length(idx)

#dat <- t(scale(t(log2(dat[idx,]+1)),center = T,scale = T))
#dim(dat)

dat <- dat[idx,]

pheatmap(dat[,
             c(intersect(colnames(dat),rownames(x@meta.data)[x@meta.data$vasu_phenotype=="MART1"]),
               intersect(colnames(dat),rownames(x@meta.data)[x@meta.data$vasu_phenotype=="None"]))],
         show_colnames = F,cluster_cols = T,show_rownames = F,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])

pheatmap(cor(dat),
         show_colnames = F,show_rownames = F,cluster_cols = T,cluster_rows = T,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])

FeaturePlot(x,features = "TNFRSF9")

DotPlot(subset(x, subset = active.ident.updated_2 == "CD8 T cells (VI)"),
        features=c("TNFRSF9","CCL5","GTSF1","MSC","GPAT3","MOSPD3","NCR3","REG4",
                   "CD8A","CD8B","CD4"), 
        group.by = "vasu_phenotype",assay = "RNA",scale = F)

genes <- names(head(sort(table(unlist(lapply(
  res_all_filtered,function(x){x$gene}))),decreasing = T),n=100))

DotPlot(x,features=genes,group.by = "vasu_phenotype",assay = "RNA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DoHeatmap(object = x,
          features = genes, 
          group.by = "vasu_phenotype",
          assay = "RNA",slot = "data")

library(pheatmap)
dim(x@assays$RNA@data[genes,])
pheatmap(x@assays$RNA@data[genes,],show_colnames = F,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])

pheatmap(t(scale(t(x@assays$RNA@data[genes,]),center = T,scale = T)),show_colnames = F,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])


pheatmap(t(scale(t(log2(x@assays$RNA@data[genes,x@meta.data$Phase=="G1"]+1)),center = T,scale = T)),
         show_colnames = F,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])

table(x@meta.data[x@meta.data$Phase=="G1",]$vasu_phenotype,
      x@meta.data[x@meta.data$Phase=="G1",]$active.ident.updated_2)

dat <- t(scale(t(log2(
  x@assays$RNA@data[
    genes,
    x@meta.data$Phase=="G2M" & x@meta.data$active.ident.updated_2=="CD8 T cells (VI)"
  ]+1)),center = T,scale = T))

dat <- t(scale(t(log2(
  x@assays$RNA@data[
    genes,
    x@meta.data$active.ident.updated_2=="CD8 T cells (VI)"
  ]+1)),center = T,scale = T))

dat <- dat[apply(dat,1,function(x){all(is.finite(x))}),]

pheatmap(dat[,
             c(intersect(colnames(dat),rownames(x@meta.data)[x@meta.data$vasu_phenotype=="CD3+41bb+"]),
               intersect(colnames(dat),rownames(x@meta.data)[x@meta.data$vasu_phenotype=="None"]))],
         show_colnames = F,cluster_cols = F,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])


dat <- 
  x@assays$RNA@data[
    genes,
    x@meta.data$Phase=="G2M" & x@meta.data$active.ident.updated_2=="CD8 T cells (VI)"]

dat <- cor(as.matrix(dat),method = "spearman")

pheatmap(dat,
         show_colnames = F,
         show_rownames = F,
         annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])


# Correlation tests ------------------------------------------------------------
cl <- x@meta.data$cdr3[which(x@meta.data$vasu_phenotype=="MART1")]
cl <- cl[!is.na(cl)]
any(cl=="")
any(cl=="NA")
cl <- unique(unlist(strsplit(cl,";")))

cl2 <- x@meta.data$cdr3[which(x@meta.data$vasu_phenotype!="MART1")]
nms <- rownames(x@meta.data)[which(x@meta.data$vasu_phenotype!="MART1")]
cl2 <- setNames(cl2,nms)
cl2 <- cl2[which(!is.na(cl2))]
cl2 <- strsplit(cl2,";")
nms_shared <- names(which(unlist(lapply(cl2,function(x){any(x %in% cl)}))==T))

# stopifnot(length(intersect(rownames(x@meta.data[x@meta.data$vasu_phenotype=="MART1",]),
#           rownames(x@meta.data[x@meta.data$vasu_phenotype!="MART1",])))==0)
# 
# stopifnot(length(intersect(rownames(x@meta.data[x@meta.data$vasu_phenotype=="MART1",]),
#   rownames(x@meta.data[(! rownames(x@meta.data) %in% nms_shared) & x@meta.data$vasu_phenotype!="MART1",])))==0)

idx_mart1 <- x@meta.data$vasu_phenotype=="MART1"
idx_not_mart1 <-  (! rownames(x@meta.data) %in% nms_shared) & (x@meta.data$vasu_phenotype!="MART1") & 
  !is.na(x@meta.data$cdr3) & x@meta.data$takara_clonotype==""

# intersect(rownames(x@meta.data[idx_mart1,]),
#           rownames(x@meta.data[idx_not_mart1,]))
# 
# stopifnot(length(intersect(x@meta.data[idx_mart1,]$cdr3,
#           x@meta.data[idx_not_mart1,]$cdr3))==0)
# 
# stopifnot(length(intersect(unlist(strsplit(x@meta.data[idx_mart1,]$cdr3,";")),
#           unlist(strsplit(x@meta.data[idx_not_mart1,]$cdr3,";"))))==0)

idx_cells <- grepl("^CD8 T cells",x@meta.data$active.ident.updated_2) & 
  (idx_mart1 | idx_not_mart1) & x@meta.data$vasu_phenotype != "Unidentified"
# table(idx_cells)
# 
meta <- x@meta.data[idx_cells,]
# 
# intersect(meta$cdr3[meta$vasu_phenotype=="MART1"],
#           meta$cdr3[meta$vasu_phenotype=="None"])
# table(meta$vasu_phenotype)

dat <- x@assays$RNA@scale.data[grep("\\.",rownames(x@assays$RNA@scale.data),invert = T),idx_cells]
dat <- dat[grep("LINC",rownames(dat),invert = T),]
dim(dat)
idx <- rowSums(dat==0) < ncol(dat)
dat <- dat[idx,]
idx <- rowMeans(dat) >= quantile(rowMeans(dat),seq(0,1,0.05))[["10%"]]
dat <- dat[idx,]
cv <- apply(dat,1,sd)/apply(dat,1,mean)
idx <- cv >= quantile(cv,seq(0,1,0.01))[["10%"]]
dat <- dat[idx,]
dim(dat)

genes_keep <- rownames(dat)



cr.s <- abs(cor(t(dat),meta$S.Score))
cr.g2m <- abs(cor(t(dat),meta$G2M.Score))
# cr.a <- cor(t(dat),as.numeric(meta$vasu_phenotype %in% c("MART1")))
# 
# plot(apply(cbind(cr.s,cr.g2m),1,max),cr.a)

genes_keep <- names(which(apply(cbind(cr.s,cr.g2m),1,max) < 0.2))

# genes <- cr.a[cr.a>0.1 & apply(cbind(cr.s,cr.g2m),1,max) < 0.2,]
# genes
# 
# pheatmap(dat[head(rownames(cr.a)[order(cr.a,decreasing = T)],n=100),],
#          show_rownames = F,show_colnames = F,
#          annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")])
# 
# 
# FeaturePlot(object = x, 
#             features = names(genes)[1], 
#             split.by = "vasu_phenotype", cells = colnames(x)[idx_cells])
# 
# sort(table(as.character(str_split_fixed(meta[meta$vasu_phenotype=="MART1",]$takara_clonotype,":",6)[,3:6])),
#      decreasing = T)[names(genes)[1]]
# 
# meta[meta$vasu_phenotype=="MART1",]$takara_clonotype[grep(names(genes)[1],meta[meta$vasu_phenotype=="MART1",]$takara_clonotype)]
# unique(meta[meta$vasu_phenotype=="None",]$takara_clonotype)


y <- x
DefaultAssay(y) <- "RNA"
y <- FindVariableFeatures(y)
y <- RunPCA(y,npcs = 30)
y <- RunUMAP(y,dims = 1:30)
y <- FindNeighbors(y,dims=1:30)
y <- FindClusters(y)

# FeaturePlot(object = y, 
#             features = names(genes)[21], 
#             split.by = "vasu_phenotype", cells = colnames(x)[idx_cells])
# 
# FeaturePlot(object = y, 
#             features = names(genes)[1], 
#             split.by = "vasu_phenotype")
# 
# DimPlot(y, group.by = "active.ident.updated_2", label = T)
# DimPlot(y, group.by = "vasu_phenotype", label = T)

#dat.cca.split[[sname]]
# 
# d <- x@assays$RNA@scale.data[grep("^TR[AB]",names(genes),invert = T,value = T),idx_cells]
# d <- d[,which(apply(d,2,sd)!=0)]
# 
# dim(d)
# 
# cr <- cor(t(as.matrix(d)))
# dim(cr)
# 
# pheatmap(cr)
# 
# cr <- cor(as.matrix(d))
# 
# pheatmap(cr,annotation_col = x@meta.data[,c("vasu_phenotype","Phase","active.ident.updated_2")],
#          show_rownames = F,show_colnames = F)
# 
# DotPlot(object = y,features = grep("^TR[AB]",names(genes),invert = T,value = T), 
#         assay = "RNA", group.by = "vasu_phenotype") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

z <- subset(y, cells = which(idx_cells))

# DotPlot(object = z,features = grep("^TR[AB]",names(genes),invert = T,value = T), 
#         assay = "RNA",split.by = "vasu_phenotype") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# tab <- table(Idents(z),z@meta.data$vasu_phenotype)
# tab <- data.frame(MART1=tab[,1],None=tab[,2])
# 
# tab$pct.mart1 <- tab[,1]/rowSums(tab)
# tab$pct.none <- tab[,2]/rowSums(tab)
# tab <- tab[order(tab$pct.mart1,decreasing = T),]
# 
# tab$pct_of_mart1 <- tab$MART1/sum(tab$MART1)
# 
# tab$ratio <- tab$pct.mart1/tab$pct_of_mart1
# 
# tab[order(tab$pct_of_mart1,decreasing = T),]
# tab[order(tab$pct.mart1,decreasing = T),]
# 
# DimPlot(z,group.by = "active.ident.updated_2")
# DimPlot(z,label = T)
# 
# 
# FeaturePlot(object = z, features = c("CD8A","NCAM1","CD14","CD4","TNFRSF9"))
# 
# FeatureScatter(object = z, feature1 = "CD8A", feature2 = "TNFRSF9", 
#                group.by = "vasu_phenotype",shuffle = T,slot = "scale.data")


extract_features <- function(x,expr_cols=NULL,meta_cols=NULL){
  ex <- as.data.frame(t(x@assays$RNA@data[expr_cols,]))
  stopifnot(identical(rownames(ex),rownames(x@meta.data)))
  dat <- cbind(ex,x@meta.data[,meta_cols])
  return(dat)
}

dat <- extract_features(z,expr_cols = unique(c(genes_keep,"CD8A","CD8B")),
                        meta_cols = c("percent_mito","percent_ribo",
                                      "S.Score","G2M.Score","vasu_phenotype"))

dat <- dat[(! z@meta.data$seurat_clusters %in% c(11,14,15,16,17)) & 
             (z@assays$RNA@counts["CD8A",] > 0) & (z@assays$RNA@counts["CD8B",] > 0) & 
             (z@assays$RNA@counts["CD4",] == 0) & 
             z@meta.data$percent_mito <= summary(dat$percent_mito)[["3rd Qu."]] & 
             z@meta.data$percent_mito >= summary(dat$percent_mito)[["1st Qu."]],]
# table(dat$vasu_phenotype)
# 
# table((z@assays$RNA@counts["CD8A",] > 0) & (z@assays$RNA@counts["CD4",] == 0))
# table((z@assays$RNA@counts["CD8B",] > 0) & (z@assays$RNA@counts["CD4",] == 0))
# table((z@assays$RNA@counts["CD8A",] > 0) & (z@assays$RNA@counts["CD8B",] > 0) & (z@assays$RNA@counts["CD4",] == 0),dat$vasu_phenotype)
# 
# summary(lm(log2(TNFRSF9+1) ~ log2(CD8A+1) + percent_mito + S.Score + G2M.Score, data = dat))
# mod <- lm(as.numeric(vasu_phenotype=="MART1") ~ log2(TNFRSF9+1) + log2(CD8A+1) + log2(CD8B+1) + percent_mito + S.Score + G2M.Score, data = dat)
# summary(mod)
# 
# DimPlot(z,label = T,cells = rownames(dat))
# FeaturePlot(object = z, features = c("CD8A","NCAM1","TNFRSF9"), cells = rownames(dat), 
#             split.by = "vasu_phenotype")
# tab <- table(TNFRSF9=dat$TNFRSF9>0,MART1=dat$vasu_phenotype=="MART1")
# tab
# fisher.test(tab)

genes_keep <- setdiff(genes_keep,c("CD8A","CD8B"))

res <- list()
for (gene in genes_keep){
  tmp <- dat[,unique(c(setdiff(colnames(dat),setdiff(genes_keep,gene)),"CD8A","CD8B"))]
  colnames(tmp)[colnames(tmp)==gene] <- "gene"
  mod <- lm(as.numeric(vasu_phenotype=="MART1") ~ 
              log2(gene+1) + log2(CD8A+1) + log2(CD8B+1) + 
              percent_mito + percent_ribo + 
              S.Score + G2M.Score, data = tmp)
  mod <- summary(mod)
  res[[gene]] <- mod$coefficients
}

res.filtered <- res[which(unlist(lapply(res,function(x){
  as.data.frame(x)["log2(gene + 1)","Pr(>|t|)"] < 0.05 & 
    as.data.frame(x)["percent_mito","Pr(>|t|)"] >= 0.05 & 
    as.data.frame(x)["percent_ribo","Pr(>|t|)"] >= 0.05 & 
    as.data.frame(x)["S.Score","Pr(>|t|)"] >= 0.05 & 
    as.data.frame(x)["G2M.Score","Pr(>|t|)"] >= 0.05 & 
    as.data.frame(x)["log2(gene + 1)","Estimate"] > 0 
})))]
length(res.filtered)

res.filtered
names(res.filtered)
#FeaturePlot(object = z, features = names(res.filtered)[1:4], cells = rownames(dat),split.by = "vasu_phenotype")
#FeaturePlot(object = z, features = names(res.filtered)[5:length(names(res.filtered))], cells = rownames(dat),split.by = "vasu_phenotype")

# Implement the above test for all samples -------------------------------------
filter_cells <- function(x,pheno){
  cl <- x@meta.data$cdr3[which(x@meta.data$vasu_phenotype==pheno)]
  cl <- cl[!is.na(cl)]
  any(cl=="")
  any(cl=="NA")
  cl <- unique(unlist(strsplit(cl,";")))
  
  cl2 <- x@meta.data$cdr3[which(x@meta.data$vasu_phenotype!=pheno)]
  nms <- rownames(x@meta.data)[which(x@meta.data$vasu_phenotype!=pheno)]
  cl2 <- setNames(cl2,nms)
  cl2 <- cl2[which(!is.na(cl2))]
  cl2 <- strsplit(cl2,";")
  nms_shared <- names(which(unlist(lapply(cl2,function(x){any(x %in% cl)}))==T))
  
  idx_mart1 <- x@meta.data$vasu_phenotype==pheno
  idx_not_mart1 <-  (! rownames(x@meta.data) %in% nms_shared) & (x@meta.data$vasu_phenotype!=pheno) & 
    !is.na(x@meta.data$cdr3) & x@meta.data$takara_clonotype==""
  
   #idx_cells <- grepl("^CD8 T cells",x@meta.data$active.ident.updated_2) & 
  #   (idx_mart1 | idx_not_mart1) & x@meta.data$vasu_phenotype != "Unidentified"
  
  idx_cells <- (idx_mart1 | idx_not_mart1) & x@meta.data$vasu_phenotype != "Unidentified"
  
  idx_cells <- idx_cells &
    ((x@assays$RNA@counts["CD8A",] > 0) | (x@assays$RNA@counts["CD8B",] > 0)) & 
    (x@assays$RNA@counts["CD4",] == 0) & (x@assays$RNA@counts["CD14",] == 0) & 
    (x@assays$RNA@counts["CD19",] == 0) & (x@assays$RNA@counts["MLANA",] == 0) & 
    (x@assays$RNA@counts["PMEL",] == 0) & (x@assays$RNA@counts["TYR",] == 0)
  
  idx_cells <- idx_cells & 
    x@meta.data$percent_mito <= summary(x@meta.data[idx_cells,]$percent_mito)[["3rd Qu."]] & 
    x@meta.data$percent_mito >= summary(x@meta.data[idx_cells,]$percent_mito)[["1st Qu."]]
  
  return(idx_cells)
}

reprocess_sample <- function(x){
  y <- x
  DefaultAssay(y) <- "RNA"
  y <- FindVariableFeatures(y)
  y <- RunPCA(y,npcs = 30)
  y <- RunUMAP(y,dims = 1:30)
  y <- FindNeighbors(y,dims=1:30)
  y <- FindClusters(y)
}

extract_features <- function(x,expr_cols=NULL,meta_cols=NULL){
  ex <- as.data.frame(t(x@assays$RNA@data[expr_cols,]))
  stopifnot(identical(rownames(ex),rownames(x@meta.data)))
  dat <- cbind(ex,x@meta.data[,meta_cols])
  return(dat)
}

run_tests3 <- function(dat, genes_keep, pheno){
  genes_keep <- setdiff(genes_keep,c("CD8A","CD8B"))
  res <- list()
  for (gene in genes_keep){
    tmp <- dat[,unique(c(setdiff(colnames(dat),setdiff(genes_keep,gene)),"CD8A","CD8B"))]
    colnames(tmp)[colnames(tmp)==gene] <- "gene"
    mod <- lm(as.numeric(vasu_phenotype == pheno) ~ 
                log2(gene+1) + log2(CD8A+1) + log2(CD8B+1) + 
                percent_mito + percent_ribo + 
                S.Score + G2M.Score + nCount_RNA, data = tmp)
    mod <- summary(mod)
    res[[gene]] <- mod$coefficients
  }
  return(res)
}

filter_res <- function(res){
  res.filtered <- res[which(unlist(lapply(res,function(x){
    as.data.frame(x)["log2(gene + 1)","Pr(>|t|)"] < 0.05 & 
      as.data.frame(x)["percent_mito","Pr(>|t|)"] >= 0.05 & 
      as.data.frame(x)["percent_ribo","Pr(>|t|)"] >= 0.05 & 
      as.data.frame(x)["nCount_RNA","Pr(>|t|)"] >= 0.05 & 
      as.data.frame(x)["S.Score","Pr(>|t|)"] >= 0.05 & 
      as.data.frame(x)["G2M.Score","Pr(>|t|)"] >= 0.05 & 
      as.data.frame(x)["log2(gene + 1)","Estimate"] > 0 
  })))]
  return(res.filtered)
}

filter_genes <- function(x,idx_cells){
  dat <- x@assays$RNA@scale.data[grep("\\.",rownames(x@assays$RNA@scale.data),invert = T),idx_cells]
  dat <- dat[grep("LINC",rownames(dat),invert = T),]
  idx <- rowSums(dat==0) < ncol(dat)
  dat <- dat[idx,]
  idx <- rowMeans(dat) >= quantile(rowMeans(dat),seq(0,1,0.05))[["10%"]]
  dat <- dat[idx,]
  cv <- apply(dat,1,sd)/apply(dat,1,mean)
  idx <- cv >= quantile(cv,seq(0,1,0.01))[["10%"]]
  dat <- dat[idx,]
  return(dat)
}

prepare_data <- function(dat.cca.split.rescaled, sname, pheno){
  x <- dat.cca.split.rescaled[[sname]];
  y <- reprocess_sample(x)
  
  idx_cells <- filter_cells(x, pheno = pheno)
  
  z <- subset(y, cells = which(idx_cells))
  
  dat <- filter_genes(x,idx_cells)
  
  cr.s <- abs(cor(t(dat),x@meta.data$S.Score[idx_cells]))
  cr.g2m <- abs(cor(t(dat),x@meta.data$G2M.Score[idx_cells]))
  
  genes_keep <- names(which(apply(cbind(cr.s,cr.g2m),1,max) < 0.2))
  
  dat <- extract_features(z,expr_cols = unique(c(genes_keep,"CD8A","CD8B")),
                          meta_cols = c("percent_mito","percent_ribo",
                                        "S.Score","G2M.Score","vasu_phenotype","nCount_RNA"))
  
  out <- list(x=x,
              y=y,
              z=z,
              dat=dat,
              genes_keep=genes_keep,
              idx_cells=idx_cells)
  
  return(out)
}

make_cluster_plots <- function(out,sname,outdir){
  z <- out[[sname]]$z
  y <- out[[sname]]$y
  dat <- out[[sname]]$dat
  
  dir.create(paste0(outdir,"plots_sample_clusters/"),showWarnings = F,recursive = T)
  
  pdf(file=paste0(outdir,"plots_sample_clusters/","dimplot.",sname,".pdf"),width=10,height=10)
  print(DimPlot(y,label = T))
  dev.off()
  
  pdf(file=paste0(outdir,"plots_sample_clusters/","featureplot.",sname,".pdf"),width=10*2,height=10*2)
  print(FeaturePlot(y, features = c("CD8A","CD8B","CD4","CD3D","CD3G","CD14", "CD19","NCAM1","PMEL"),order = T))
  dev.off()
}

# dat.cca.split.rescaled <- lapply(dat.cca.split, function(x){
#   x <- ScaleData(x, verbose = T, vars.to.regress = c("S.Score", "G2M.Score"))
#   return(x)
# })

#saveRDS(dat.cca.split.rescaled,file = paste0(outdir,"dat.cca.split.rescaled.rda"))
#dat.cca.split.rescaled <- readRDS(file = paste0(outdir,"dat.cca.split.rescaled.rda"))

snames.mart1 <- sample_map[sample_map$in_takara & sample_map$vasu_phenotype == "MART1",]$Var2
snames.41bb <- sample_map[sample_map$in_takara & sample_map$vasu_phenotype == "CD3+41bb+",]$Var2

out.mart1 <- list()
for (sname in snames.mart1){
  out.mart1[[sname]] <- prepare_data(dat.cca.split.rescaled = dat.cca.split.rescaled, 
                      sname = sname, 
                      pheno = "MART1")
}

out.41bb <- list()
for (sname in snames.41bb){
  out.41bb[[sname]] <- prepare_data(dat.cca.split.rescaled = dat.cca.split.rescaled, 
                      sname = sname, 
                      pheno = "CD3+41bb+")
}

#saveRDS(out.mart1,file = paste0(outdir,"out.mart1.rda"))
#saveRDS(out.41bb,file = paste0(outdir,"out.41bb.rda"))

for (sname in names(out.mart1)){
  make_cluster_plots(out = out.mart1,
                     sname = sname,
                     outdir = outdir)
}

for (sname in names(out.41bb)){
  make_cluster_plots(out = out.41bb,
                     sname = sname,
                     outdir = outdir)
}

new.idents <- list()

new.idents[["A2-GEX"]] <- c("0" = "CD8 T cells",
                "1" = "CD8 T cells", 
                "2" = "CD4 T cells",
                "3" = "CD8 T cells",
                "4" = "NK cells",
                "5" = "Melanocytic",
                "6" = "NK cells",
                "7" = "CD8 T cells",
                "8" = "CD8 T cells",
                "9" = "CD8 T cells",
                "10" = "Melanocytic",
                "11" = "CD4 T cells",
                "12" = "CD8 T cells",
                "13" = "CD8 T cells",
                "14" = "14_Mono_DC",
                "15" = "15_Mono_DC",
                "16" = "B cells",
                "17" = "17_Mono_DC")

new.idents[["SampleID_2_11june18"]] <- c("0" = "NK cells",
                                         "1" = "NK cells",
                                         "2" = "CD8 T cells",
                                         "3" = "CD8 T cells",
                                         "4" = "CD8 T cells",
                                         "5" = "NK cells",
                                         "6" = "CD8 T cells",
                                         "7" = "CD8 T cells",
                                         "8" = "CD4 T cells",
                                         "9" = "CD8 T cells",
                                         "10" = "CD8 T cells")

new.idents[["C7-GEX"]] <- c("0" = "Melanocytic",
                            "1" = "Melanocytic",
                            "2" = "CD8 T cells",
                            "3" = "CD4 T cells",
                            "4" = "CD8 T cells",
                            "5" = "5_Mono_DC_B",
                            "6" = "CD4 T cells",
                            "7" = "Melanocytic / CD8 T cells",
                            "8" = "NK cells",
                            "9" = "CD4 T cells",
                            "10" = "Melanocytic")

new.idents[["D4-GEX"]] <- c("0" = "Melanocytic",
                            "1" = "Melanocytic",
                            "2" = "CD8 / CD4 / NK cells",
                            "3" ="Melanocytic",
                            "4" ="Melanocytic")

new.idents[["SampleID_7_11june18"]] <- c("0" = "CD8 T cells",
                                         "1" = "CD8 T cells",
                                         "2" = "CD8 T cells",
                                         "3" = "CD8 T cells",
                                         "4" = "CD8 T cells",
                                         "5" = "CD4 T cells",
                                         "6" = "CD8 T cells",
                                         "7" = "CD4 T cells",
                                         "8" = "CD8 T cells",
                                         "9" = "CD8 T cells",
                                         "10" = "CD4 T cells",
                                         "11" = "NK cells")

new.idents[["SampleID_9_11june18"]] <- c("0" = "CD8 T cells",
                                         "1" = "CD8 T cells",
                                         "2" = "CD8 T cells",
                                         "3" = "CD8 T cells",
                                         "4" = "CD8 T cells",
                                         "5" = "CD8 T cells",
                                         "6" = "CD8 T cells",
                                         "7" = "CD8 T cells",
                                         "8" = "CD4 T cells",
                                         "9" = "CD4 T cells",
                                         "10" = "CD8 T cells",
                                         "11" = "NK cells")

for (sname in names(out.mart1)){
  out.mart1[[sname]]$y.renamed <- RenameIdents(out.mart1[[sname]]$y,new.idents[[sname]])
  pdf(file=paste0(outdir,"plots_sample_clusters/","dimplot.renamed.",sname,".pdf"),width=10,height=10)
  print(DimPlot(out.mart1[[sname]]$y.renamed))
  dev.off()
}

for (sname in names(out.41bb)){
  out.41bb[[sname]]$y.renamed <- RenameIdents(out.41bb[[sname]]$y,new.idents[[sname]])
  pdf(file=paste0(outdir,"plots_sample_clusters/","dimplot.renamed.",sname,".pdf"),width=10,height=10)
  print(DimPlot(out.41bb[[sname]]$y.renamed))
  dev.off()
}

#saveRDS(out.mart1,file = paste0(outdir,"out.mart1.rda"))
#saveRDS(out.41bb,file = paste0(outdir,"out.41bb.rda"))

# Manual step: determine clusters to remove from analysis

#DimPlot(z,label = T)
#FeaturePlot(z, features = c("CD8A","CD8B","CD3D","CD3G"),order = T)
# mrk <- FindAllMarkers(y,only.pos = T)
# mrk.bak <- mrk
# mrk <- mrk[mrk$p_val_adj<0.05,]
# mrk <- mrk[order(mrk$avg_log2FC,decreasing = T),]
# head(mrk[mrk$cluster==15,])

# head(mrk[mrk$cluster==17,])



# Back to automated analysis
#lapply(out.mart1,function(x){levels(Idents(x$y.renamed))})
#lapply(out.41bb,function(x){levels(Idents(x$y.renamed))})
clusters_retain <- c("CD8 T cells","CD8 / CD4 / NK cells")
for (sname in names(out.mart1)){
  dat <- out.mart1[[sname]]$dat
  y <- out.mart1[[sname]]$y.renamed
  genes_keep <- out.mart1[[sname]]$genes_keep
  
  dat <- dat[intersect(rownames(dat),names(Idents(y))[Idents(y) %in% clusters_retain]),]
  
  res <- run_tests3(dat,genes_keep, pheno = "MART1")
  res.filtered <- filter_res(res)
  
  out.mart1[[sname]]$res <- res
  out.mart1[[sname]]$res.filtered <- res.filtered
  out.mart1[[sname]]$res.filtered.genes <- names(res.filtered)
}

get_pathways <- function(expr, pw_path="~/proj/pemdac/Investigations/data/gene_sets/c2.cp.reactome.v7.1.symbols.gmt"){
  library(SeqGSEA)
  
  gene.set <- list()
  gene.set[["perturb"]] <- loadGenesets(pw_path,
                                        rownames(expr), geneID.type="gene.symbol",
                                        genesetsize.min = 0, genesetsize.max = 10000)
  pathways <- list()
  for (pw in names(gene.set)){
    for (i in 1:length(gene.set[[pw]]@GS)){
      gene_set_name <- gene.set[[pw]]@GSNames[[i]]
      gene_set_genes <- gene.set[[pw]]@geneList[gene.set[[pw]]@GS[[i]]]
      pathways[[gene_set_name]] <- gene_set_genes
    }
  }
  return(pathways)
}

for (sname in names(out.41bb)){
  dat <- out.41bb[[sname]]$dat
  y <- out.41bb[[sname]]$y.renamed
  genes_keep <- out.41bb[[sname]]$genes_keep
  
  #dat <- dat[intersect(rownames(dat),names(Idents(y))[Idents(y) %in% clusters_retain]),]
  if (sname == "SampleID_9_11june18"){
    clusters_retain <- c(3,4,5,6,7,10)
    dat <- dat[intersect(rownames(dat),colnames(y)[y@meta.data$seurat_clusters %in% clusters_retain]),]
    genes_remove <- intersect(genes_keep,names(which(apply(dat,2,function(x){length(unique(x))})==1)))
    genes_remove <- setdiff(genes_remove,c("CD8A","CD8B"))
    dat <- dat[,! colnames(dat) %in% genes_remove]
    genes_keep <- setdiff(genes_keep,genes_remove)
  }
  stopifnot(identical(rownames(dat),colnames(subset(y,cells=rownames(dat)))))
  cr.s <- abs(cor(as.matrix(dat[,genes_keep]),subset(y,cells=rownames(dat))@meta.data$S.Score))
  cr.g2m <- abs(cor(as.matrix(dat[,genes_keep]),subset(y,cells=rownames(dat))@meta.data$G2M.Score))
  
  tmp <- names(which(apply(cbind(cr.s,cr.g2m),1,max,na.rm=T) < 0.2))
  genes_remove <- setdiff(genes_keep,tmp)
  genes_remove <- setdiff(genes_remove,c("CD8A","CD8B"))
  dat <- dat[,! colnames(dat) %in% genes_remove]
  genes_keep <- setdiff(genes_keep,genes_remove)
  
  stopifnot(identical(sort(setdiff(genes_keep,genes_remove)),sort(c(tmp,"CD8A","CD8B"))))
  all(genes_keep %in% colnames(dat))
  
  res <- run_tests3(dat,genes_keep, pheno = "CD3+41bb+")
  res.filtered <- filter_res(res)
  names(res.filtered)
  
  RidgePlot(subset(y,cells=intersect(rownames(dat),colnames(y)[y@meta.data$seurat_clusters %in% clusters_retain])),features = c("nCount_RNA"),group.by = "vasu_phenotype")
  FeaturePlot(subset(y,cells=intersect(rownames(dat),colnames(y)[y@meta.data$seurat_clusters %in% clusters_retain])),features = c("nCount_RNA"),split.by = "vasu_phenotype")
  
  FeaturePlot(subset(y,cells=intersect(rownames(dat),colnames(y)[y@meta.data$seurat_clusters %in% clusters_retain])),features = c("TNFRSF9"),split.by = "vasu_phenotype")
  
  sub <- subset(y,cells=intersect(rownames(dat),colnames(y)[y@meta.data$seurat_clusters %in% clusters_retain]))
  
  DimPlot(y, group.by = "cdr3",split.by = "vasu_phenotype",label = F) + NoLegend()
  DimPlot(y, group.by = "Phase",split.by = "vasu_phenotype",label = F)
  
  RidgePlot(subset(y,cells=colnames(y)[!is.na(y@meta.data$cdr3)]),
            features = c("nCount_RNA"),group.by = "vasu_phenotype")
  
  expr <- y@assays$RNA@counts[genes_keep,rownames(dat)]
  expr <- as.matrix(expr)
  
  library(GSVA)
  
  pathways <- list()
  pathways[["reactome"]] <- get_pathways(expr, pw_path="~/proj/pemdac/Investigations/data/gene_sets/c2.cp.reactome.v7.1.symbols.gmt")
  pathways[["c1_all"]] <- get_pathways(expr, pw_path="~/proj/pemdac/Investigations/data/gene_sets/c1.all.v7.1.symbols.gmt")
  pathways[["go_bp"]] <- get_pathways(expr, pw_path="~/proj/pemdac/Investigations/data/gene_sets/c5.bp.v7.1.symbols.gmt")
  pathways[["c8.all.v7.5.1"]] <- get_pathways(expr, pw_path="~/proj/common_data/MSigDB/c8.all.v7.5.1.symbols.gmt")
  pathways[["c7.all.v7.5.1"]] <- get_pathways(expr, pw_path="~/proj/common_data/MSigDB/c7.all.v7.5.1.symbols.gmt")
  
  ssgsea.res.mat <- list()
  ssgsea.res.mat[["reactome"]] <- gsva(expr = as.matrix(expr), method="ssgsea", 
                         kcdf="Poisson", gset.idx.list = pathways[["reactome"]])
  ssgsea.res.mat[["c1_all"]] <- gsva(expr = as.matrix(expr), method="ssgsea", 
                                       kcdf="Poisson", gset.idx.list = pathways[["c1_all"]])
  ssgsea.res.mat[["go_bp"]] <- gsva(expr = as.matrix(expr), method="ssgsea", 
                                       kcdf="Poisson", gset.idx.list = pathways[["go_bp"]])
  ssgsea.res.mat[["c8.all.v7.5.1"]] <- gsva(expr = as.matrix(expr), method="ssgsea", 
                                    kcdf="Poisson", gset.idx.list = pathways[["c8.all.v7.5.1"]])
  ssgsea.res.mat[["c7.all.v7.5.1"]] <- gsva(expr = as.matrix(expr), method="ssgsea", 
                                            kcdf="Poisson", gset.idx.list = pathways[["c7.all.v7.5.1"]])
  
  
  cv <- apply(ssgsea.res.mat,1,sd)/apply(ssgsea.res.mat,1,mean)
  pw.keep <- names(sort(cv,decreasing = T))[1:25]
  
  pheatmap(ssgsea.res.mat[pw.keep,],show_rownames = F,show_colnames = F,
           annotation_col = y@meta.data[colnames(expr),c("S.Score","G2M.Score","vasu_phenotype")])
  
  pheatmap(ssgsea.res.mat[pw.keep,],show_rownames = T,show_colnames = F,
           annotation_col = y@meta.data[colnames(expr),c("S.Score","G2M.Score","vasu_phenotype")])
  
  mean_pw <- cbind(rowMeans(ssgsea.res.mat[["c1_all"]][,rownames(dat)[dat$vasu_phenotype=="CD3+41bb+"]]),
                   rowMeans(ssgsea.res.mat[["c1_all"]][,rownames(dat)[dat$vasu_phenotype=="None"]]))
  
  pheatmap(mean_pw[order(mean_pw[,1]-mean_pw[,2],decreasing = T),][1:25,],
        show_rownames = T,show_colnames = F,cluster_rows = F,cluster_cols = F)
  
  
  out.41bb[[sname]]$res <- res
  out.41bb[[sname]]$res.filtered <- res.filtered
  out.41bb[[sname]]$res.filtered.genes <- names(res.filtered)
}

lapply(out.mart1,function(x){sort(x$res.filtered.genes)})
lapply(out.41bb,function(x){sort(x$res.filtered.genes)})

head(sort(table(unlist(lapply(out.41bb,function(x){x$res.filtered.genes}))),decreasing = T))

# Instead of clusters, compare clonotypes?
sname <- "SampleID_9_11june18"
dat <- out.41bb[[sname]]$dat
y <- out.41bb[[sname]]$y.renamed
genes_keep <- out.41bb[[sname]]$genes_keep

DimPlot(y, group.by = "cdr3",split.by = "vasu_phenotype",label = F) + NoLegend()
DimPlot(y, group.by = "Phase",split.by = "vasu_phenotype",label = F)

idx_cells <- filter_cells(y, "CD3+41bb+")
idx_cells <- idx_cells & !is.na(y@meta.data$cdr3)
idx_cells <- idx_cells & y@meta.data$vasu_phenotype %in% c("CD3+41bb+","None")

DimPlot(subset(y,cells=names(idx_cells)[idx_cells]))

DimPlot(subset(y,cells=names(idx_cells)[idx_cells]),group.by = "takara_clonotype") + NoLegend()
DimPlot(subset(y,cells=names(idx_cells)[idx_cells]),group.by = "cdr3") + NoLegend()
head(y@meta.data)

DimPlot(subset(y,cells=names(idx_cells)[idx_cells]),group.by = "takara_clonotype",split.by = "vasu_phenotype") + NoLegend()
DimPlot(subset(y,cells=names(idx_cells)[idx_cells]),group.by = "cdr3",split.by = "vasu_phenotype") + NoLegend()

cdrs <- table(subset(y,cells=names(idx_cells)[idx_cells])@meta.data$cdr3)
sort(cdrs,decreasing = T)
cdrs[cdrs>5]




s <- subset(y,cells=names(idx_cells)[idx_cells])
unique(s@meta.data$takara_clonotype[s@meta.data$vasu_phenotype=="None"])

tst <- FindMarkers(object = s, test.use = "LR" , only.pos = T, group.by = "vasu_phenotype",
                   ident.1 = "CD3+41bb+",ident.2 = "None",
                   latent.vars = c("S.Score","G2M.Score","percent_mito","percent_ribo","nCount_RNA"))


s2 <- subset(s, cells=colnames(s)[s@meta.data$cdr3 %in% names(cdrs[cdrs>5])])
unique(s2$cdr3)

table(s2$cdr3,s2$vasu_phenotype)
Idents(s2) <- s2@meta.data$cdr3

tst2 <- FindAllMarkers(s2, test.use = "LR", only.pos = T, 
                       latent.vars = c("S.Score","G2M.Score","percent_mito","percent_ribo","nCount_RNA"))

table(subset(y,cells=names(idx_cells)[idx_cells])@meta.data$takara_clonotype)

cls <- setdiff(names(table(subset(s,cells=names(idx_cells)[idx_cells])@meta.data$takara_clonotype)[
  table(subset(s,cells=names(idx_cells)[idx_cells])@meta.data$takara_clonotype)>2]),"")
tst3 <- list()
for (cl in cls){
  tst3[[cl]] <- FindMarkers(object = s, test.use = "LR" , only.pos = T, 
                            group.by = "takara_clonotype",
                     ident.1 = cl, ident.2 = "",
                     latent.vars = c("S.Score","G2M.Score","percent_mito","percent_ribo","nCount_RNA"))
}


tst4 <- FindAllMarkers(s2, test.use = "LR", only.pos = T, 
                       latent.vars = c("S.Score","G2M.Score",
                                       "percent_mito","percent_ribo","nCount_RNA",
                                       "seurat_cluster"))

tst5 <- FindMarkers(object = s, test.use = "LR" , only.pos = T, group.by = "vasu_phenotype",
                   ident.1 = "CD3+41bb+",ident.2 = "None",
                   latent.vars = c("S.Score","G2M.Score","percent_mito",
                                   "percent_ribo","nCount_RNA",
                                   "seurat_cluster"))

tst6 <- FindMarkers(object = s, test.use = "LR" , only.pos = T, group.by = "vasu_phenotype",
                    ident.1 = "CD3+41bb+",ident.2 = "None",
                    latent.vars = c("S.Score","G2M.Score","percent_mito",
                                    "percent_ribo","nCount_RNA",
                                    "seurat_cluster","cdr3"))

tst7 <- FindMarkers(object = s, test.use = "negbinom" , only.pos = T, group.by = "vasu_phenotype",
                    ident.1 = "CD3+41bb+",ident.2 = "None", assay = "RNA", slot = "data",
                    latent.vars = c("nCount_RNA"))

head(tst7)

gene <- "KIAA1671"

FeaturePlot(s, features = gene, order = T) + DimPlot(s, group.by = "vasu_phenotype",order = T) + 
  FeaturePlot(s, features = "nCount_RNA",order = T) + FeaturePlot(s, features = "G2M.Score",order = T)

summary(lm(as.numeric(s@meta.data$vasu_phenotype == "CD3+41bb+") ~ 
     log2(s@assays$RNA@data[gene,]+1) + log2(s@meta.data$nCount_RNA+1) + 
       s@meta.data$S.Score + s@meta.data$G2M.Score + s@meta.data$percent_mito + 
       s@meta.data$percent_ribo + s@meta.data$seurat_clusters))

r.41bb <- lm(s@meta.data$nCount_RNA[s@meta.data$vasu_phenotype == "CD3+41bb+"] ~ 
             s@assays$RNA@data[gene,s@meta.data$vasu_phenotype == "CD3+41bb+"])

r.none <- lm(s@meta.data$nCount_RNA[s@meta.data$vasu_phenotype == "None"] ~ 
             s@assays$RNA@data[gene,s@meta.data$vasu_phenotype == "None"])

FeatureScatter(s,feature1 = gene,feature2 = "nCount_RNA",
               group.by = "vasu_phenotype",shuffle = T,slot = "data",
               plot.cor = T) + 
  geom_abline(slope = coef(r.none)[[2]], 
              intercept = coef(r.none)[["(Intercept)"]]) + 
  geom_abline(slope = coef(r.41bb)[[2]], 
              intercept = coef(r.41bb)[["(Intercept)"]])

  

test_gene_wx <- function(gene){
  p.vasu_phenotype <- wilcox.test(
    s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "CD3+41bb+"],
    s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "None"])$p.value
  
  p.g1s <- wilcox.test(s@assays$RNA@data[gene,][s@meta.data$Phase == "G1"],
              s@assays$RNA@data[gene,][s@meta.data$Phase == "S"])$p.value
  
  p.g1g2m <- wilcox.test(s@assays$RNA@data[gene,][s@meta.data$Phase == "G1"],
              s@assays$RNA@data[gene,][s@meta.data$Phase == "G2M"])$p.value
  
  p.sg2m <- wilcox.test(s@assays$RNA@data[gene,][s@meta.data$Phase == "S"],
              s@assays$RNA@data[gene,][s@meta.data$Phase == "G2M"])$p.value
  
  p.rna <- cor.test(s@assays$RNA@data[gene,],s@meta.data$nCount_RNA)$p.value
  p.mito <- cor.test(s@assays$RNA@data[gene,],s@meta.data$percent_mito)$p.value
  p.ribo <- cor.test(s@assays$RNA@data[gene,],s@meta.data$percent_ribo)$p.value
  
  s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "CD3+41bb+"]
  s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "None"]
  pct.41bb <- sum(s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "CD3+41bb+"] > 0)/sum(s@meta.data$vasu_phenotype == "CD3+41bb+")
  pct.none <- sum(s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "None"] > 0)/sum(s@meta.data$vasu_phenotype == "None")
  
  res <- data.frame(p.vasu_phenotype = p.vasu_phenotype,
                    p.g1s = p.g1s,
                    p.g1g2m = p.g1g2m,
                    p.sg2m = p.sg2m, 
                    p.rna = p.rna, 
                    p.mito = p.mito, 
                    p.ribo = p.ribo,
                    pct.41bb = pct.41bb,
                    pct.none = pct.none)
  return(res)
}

test_gene_wx(gene)

rs <- lapply(genes_keep,test_gene_wx)
rs <- do.call("rbind",rs)
rownames(rs) <- genes_keep
head(rs)
rs$q.vasu_phenotype <- p.adjust(rs$p.vasu_phenotype,method = "BH")

rs.filtered <- rs[which(rs$p.vasu_phenotype < 0.05 & rs$p.rna >= 0.05 & rs$p.mito >= 0.05 & 
                rs$p.ribo >= 0.05 & rs$p.g1s >= 0.05 & rs$p.g1g2m >= 0.05 & 
                rs$p.sg2m >= 0.05 & rs$q.vasu_phenotype < 0.05 & 
                rs$pct.41bb > rs$pct.none),]
dim(rs.filtered)
rs.filtered

#wilcox.test(s@meta.data$nCount_RNA[s@meta.data$vasu_phenotype == "CD3+41bb+"],
#            s@meta.data$nCount_RNA[s@meta.data$vasu_phenotype == "None"])$p.value

#wilcox.test(s@meta.data$percent_mito[s@meta.data$vasu_phenotype == "CD3+41bb+"],
#            s@meta.data$percent_mito[s@meta.data$vasu_phenotype == "None"])$p.value

#wilcox.test(s@meta.data$percent_ribo[s@meta.data$vasu_phenotype == "CD3+41bb+"],
#            s@meta.data$percent_ribo[s@meta.data$vasu_phenotype == "None"])$p.value






########################
head(s@meta.data[,c("v_gene","d_gene","j_gene","c_gene","takara_clonotype")])

lst <- apply(s@meta.data[,c("v_gene","j_gene","c_gene")],1,function(x){unlist(strsplit(x,";"))})
d <- str_split_fixed(s@meta.data$takara_clonotype,":",6)[,c(3,5,6)]
d.lst <- list()
for (i in 1:length(lst)){
  dd <- c()
  for (j in 1:nrow(d)){
    dd <- any(dd,all(d[j,] %in% lst[[i]]))
  }
  d.lst[[i]] <- dd
}

length(d.lst)

table(s@meta.data$takara_clonotype=="" & unlist(d.lst))
table(s@meta.data$takara_clonotype!="" & unlist(d.lst))
s@meta.data[which(s@meta.data$takara_clonotype=="" & unlist(d.lst)),]
