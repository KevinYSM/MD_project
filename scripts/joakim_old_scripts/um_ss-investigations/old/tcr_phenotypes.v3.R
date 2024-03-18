library(rliger)
library(Seurat)
library(tidyverse)
library(readxl)

outdir <- "~/proj/um_ss/Investigations/seurat/results/liger_all/"
lig <- readRDS(paste0(outdir,"lig.blood_tumor_tils.integrated.rda"))
dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
DefaultAssay(dat.cca) <- "RNA"
dat.cca <- CellCycleScoring(object = dat.cca,
                            s.features = cc.genes$s.genes, 
                            g2m.features = cc.genes$g2m.genes, set.ident = FALSE)

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
                                                     "J segment","C segment","Read Count")],
                                            tmp_B[,c("Sample","Chain",
                                                     "CDR3 Amino Acid Sequence",
                                                     "V segment","D segment",
                                                     "J segment","C segment","Read Count")]))
}
vdj.takara <- do.call("rbind",vdj.takara)
rownames(vdj.takara) <- NULL
colnames(vdj.takara)[colnames(vdj.takara)=="Read Count"] <- "vasu_read_count"
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

# Plot clonotype read count bar plot (pre merging with 10x) --------------------
vdj.takara$vasu_read_count <- as.numeric(vdj.takara$vasu_read_count)

vdj.takara$id <- apply(vdj.takara[,setdiff(colnames(vdj.takara),"vasu_read_count")],1,paste0,collapse="#")

vdj.takara.aggr <- aggregate(vdj.takara$vasu_read_count,
                             by=list(vdj.takara$id),sum)

vdj.takara.aggr <- cbind(as.data.frame(str_split_fixed(vdj.takara.aggr$Group.1,"#",(ncol(vdj.takara)-2))),vdj.takara.aggr$x)

colnames(vdj.takara.aggr) <- c(setdiff(colnames(vdj.takara),
                                       c("vasu_read_count","id")),"vasu_read_count")

get_clonotype_cutoff <- function(vdj.takara.aggr.sample, pct){
  qt <- quantile(sort(vdj.takara.aggr.sample$vasu_read_count),seq(0,1,by=0.05))
  print(qt)
  return(qt[[paste0(pct,"%")]])
}

thresh_A <- get_clonotype_cutoff(vdj.takara.aggr[vdj.takara.aggr$Sample=="A",],95)
thresh_B <- get_clonotype_cutoff(vdj.takara.aggr[vdj.takara.aggr$Sample=="B",],95)
thresh_C <- get_clonotype_cutoff(vdj.takara.aggr[vdj.takara.aggr$Sample=="C",],95)
thresh_D <- get_clonotype_cutoff(vdj.takara.aggr[vdj.takara.aggr$Sample=="D",],95)

plot_clonotypes <- function(vdj.takara.aggr,sname,thresh){
  dat <- vdj.takara.aggr[vdj.takara.aggr$Sample==sname,]
  dat <- dat[dat$vasu_read_count >= thresh,]
  g <- ggplot(dat,
         aes(x=reorder(clonotype,-vasu_read_count),y=vasu_read_count)) + theme_classic() +
    geom_bar(stat="identity",position="dodge") +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + xlab(NULL)
  return(g)
}

plot_clonotypes(vdj.takara.aggr,"A",thresh_A)
plot_clonotypes(vdj.takara.aggr,"B",thresh_B)
plot_clonotypes(vdj.takara.aggr,"C",thresh_C)
plot_clonotypes(vdj.takara.aggr,"D",thresh_D)

#plot_clonotypes(vdj.takara.aggr,"A",0) + ylim(0,200)

vdj.takara <- vdj.takara.aggr

vdj.takara$pass_threshold <- F
vdj.takara[vdj.takara$Sample == "A" & vdj.takara$vasu_read_count > thresh_A,]$pass_threshold <- T
vdj.takara[vdj.takara$Sample == "B" & vdj.takara$vasu_read_count > thresh_B,]$pass_threshold <- T
vdj.takara[vdj.takara$Sample == "C" & vdj.takara$vasu_read_count > thresh_C,]$pass_threshold <- T
vdj.takara[vdj.takara$Sample == "D" & vdj.takara$vasu_read_count > thresh_D,]$pass_threshold <- T

# Find clonotypes and samples shared with Vasu's data --------------------------

cnames <- c("Sample.ID","project.name","UM.ID",
            "v_gene","d_gene","j_gene","c_gene","chains","cdr3",
            "reads","umis","productive","djvdj.clone_freq","djvdj.clone_pct")

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

#table(vdj.10x[rownames(vdj.10x.match[idx_all_match,]),]$project.name)
#table(vdj.10x[rownames(vdj.10x.match[idx_all_match,]),]$UM.ID)

#vdj.10x$all_match <- rownames(vdj.10x) %in% rownames(vdj.10x.match[idx_all_match,])

#head(vdj.10x[vdj.10x$all_match != "",])

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

#head(vdj.10x[widx,])


# table(vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype)
# table(vdj.takara$UM.ID[vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype])
# table(vdj.takara$vasu_phenotype[vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype])
# table(vdj.takara$vasu_phenotype[vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype & vdj.takara$Sample!="Undetermined"],
#       vdj.takara$UM.ID[vdj.takara$clonotype %in% vdj.10x[widx,]$takara_clonotype & vdj.takara$Sample!="Undetermined"])
# 
# table(vdj.takara[vdj.takara$Sample!="Undetermined",]$vasu_phenotype,
#       vdj.takara[vdj.takara$Sample!="Undetermined",]$UM.ID)

stopifnot(identical(rownames(vdj.10x),rownames(dat.cca@meta.data)))
dat.cca <- AddMetaData(dat.cca,
                       metadata = vdj.10x[,
                                          setdiff(colnames(vdj.10x),
                                                  colnames(dat.cca@meta.data))])

intersect(dat.cca@meta.data$takara_clonotype,vdj.takara$clonotype)

x <- merge(dat.cca@meta.data,vdj.takara[,c("clonotype","vasu_phenotype",
                                           "vasu_read_count","pass_threshold")],
           by.x="takara_clonotype",by.y="clonotype",all.x=T,all.y=F)

vasu_phenotype <- setNames(x$vasu_phenotype,x$cell_id)
vasu_phenotype <- vasu_phenotype[rownames(dat.cca@meta.data)]
identical(names(vasu_phenotype),rownames(dat.cca@meta.data))
vasu_phenotype[is.na(vasu_phenotype)] <- "None"
vasu_phenotype[vasu_phenotype==""] <- "Unidentified"

vasu_read_count <- setNames(x$vasu_read_count,x$cell_id)
vasu_read_count <- vasu_read_count[rownames(dat.cca@meta.data)]
identical(names(vasu_read_count),rownames(dat.cca@meta.data))
vasu_read_count[is.na(vasu_read_count)] <- 0

pass_threshold <- setNames(x$pass_threshold,x$cell_id)
pass_threshold <- pass_threshold[rownames(dat.cca@meta.data)]
identical(names(pass_threshold),rownames(dat.cca@meta.data))
pass_threshold[is.na(pass_threshold)] <- F

dat.cca <- AddMetaData(dat.cca,metadata = vasu_phenotype,col.name = "vasu_phenotype")
dat.cca <- AddMetaData(dat.cca,metadata = vasu_read_count,col.name = "vasu_read_count")
dat.cca <- AddMetaData(dat.cca,metadata = pass_threshold,col.name = "pass_threshold")

#saveRDS(dat.cca,paste0(outdir,"dat.cca.modified.rda"))
#dat.cca <- readRDS(paste0(outdir,"dat.cca.modified.rda"))

# Within each dataset and sample, locate the Vasu clonotypes -------------------

idx <- vdj.10x$all_match != ""
vdj.10x[idx,]
vdj.10x[!idx,]

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

#saveRDS(sample_map,paste0(outdir,"sample_map.vasu.rda"))
#sample_map <- readRDS(paste0(outdir,"sample_map.vasu.rda"))

# # Investigate the frequency of each clonotype in Vasu's data (per patient and experiment) -----
# 
# df <- data.frame(clonotype = dat.cca@meta.data$takara_clonotype,
#            experiment = dat.cca@meta.data$vasu_phenotype,
#            sample = dat.cca@meta.data$Sample.ID)
# 
# df <- df[df$sample %in% sample_map[sample_map$vasu_phenotype!="",]$Var2,]
# df <- df[df$clonotype!="",]
# df <- df[df$clonotype!="Unidentified",]
# 
# head(df)
# 
# table(df$clonotype,df$sample)
# table(df$experiment,df$sample)

# Redefine passing clonotypes --------------------------------------------------

# Clonotypes that:
# 1) Exist in both vasu, biopsy and TIL culture
# 2) n:th quantile of these, per sample

stopifnot(max(unlist(lapply(strsplit(dat.cca@meta.data$takara_clonotype,";"),length)))==3)

tmp <- as.data.frame(do.call("rbind",lapply(strsplit(dat.cca@meta.data$takara_clonotype,";"),function(x){
  tra <- grep("^TRA",x,value = T)
  trb <- grep("^TRB",x,value = T)
  
  tra_1 <- ifelse(is.na(tra[1]),"",tra[1])
  tra_2 <- ""
  if (length(tra)==2){
    tra_2 <- tra[2]
  }
  trb_1 <- ifelse(is.na(trb[1]),"",trb[1])
  trb_2 <- ""
  if (length(trb)==2){
    trb_2 <- trb[2]
  }
  
  return(list(
    tra_1 = tra_1,
    tra_2 = tra_2,
    trb_1 = trb_1,
    trb_2 = trb_2
  ))
})))

dat.cca@meta.data$takara_clonotype_TRA_1 <- unlist(tmp$tra_1)
dat.cca@meta.data$takara_clonotype_TRA_2 <- unlist(tmp$tra_2)
dat.cca@meta.data$takara_clonotype_TRB_1 <- unlist(tmp$trb_1)
dat.cca@meta.data$takara_clonotype_TRB_2 <- unlist(tmp$trb_2)

dat.cca@meta.data$takara_clonotype.common_biopsy_til <- F

vdj.takara$in_til <- F
vdj.takara$in_biopsy <- F
vdj.takara$in_both <- F
vdj.takara$represented <- ""

vdj.takara$percent <- NA

snames <- setdiff(unique(vdj.takara$UM.ID),"")
for (sname in snames){
  sname.til <- unique(dat.cca@meta.data$Sample.ID[
    dat.cca@meta.data$UM.ID==sname & dat.cca@meta.data$project.name=="G18-023"])
  sname.biopsy <- unique(dat.cca@meta.data$Sample.ID[
    dat.cca@meta.data$UM.ID==sname & dat.cca@meta.data$project.name=="GWA-JN-388"])
  
  tra_1.clonotypes.til <- dat.cca@meta.data$takara_clonotype_TRA_1[dat.cca@meta.data$Sample.ID==sname.til]
  tra_2.clonotypes.til <- dat.cca@meta.data$takara_clonotype_TRA_2[dat.cca@meta.data$Sample.ID==sname.til]
  trb_1.clonotypes.til <- dat.cca@meta.data$takara_clonotype_TRB_1[dat.cca@meta.data$Sample.ID==sname.til]
  trb_2.clonotypes.til <- dat.cca@meta.data$takara_clonotype_TRB_2[dat.cca@meta.data$Sample.ID==sname.til]
  
  tra_1.clonotypes.biopsy <- dat.cca@meta.data$takara_clonotype_TRA_1[dat.cca@meta.data$Sample.ID==sname.biopsy]
  tra_2.clonotypes.biopsy <- dat.cca@meta.data$takara_clonotype_TRA_2[dat.cca@meta.data$Sample.ID==sname.biopsy]
  trb_1.clonotypes.biopsy <- dat.cca@meta.data$takara_clonotype_TRB_1[dat.cca@meta.data$Sample.ID==sname.biopsy]
  trb_2.clonotypes.biopsy <- dat.cca@meta.data$takara_clonotype_TRB_2[dat.cca@meta.data$Sample.ID==sname.biopsy]
  
  clonotypes.common <- setdiff(intersect(
    c(tra_1.clonotypes.til,
          tra_2.clonotypes.til,
          trb_1.clonotypes.til,
          trb_2.clonotypes.til),
    c(tra_1.clonotypes.biopsy,
          tra_2.clonotypes.biopsy,
          trb_1.clonotypes.biopsy,
          trb_2.clonotypes.biopsy)
    ),"")
  
  vdj.takara.clonotypes.common <- vdj.takara[which(vdj.takara$clonotype %in% 
                                                     clonotypes.common & 
                                                     vdj.takara$UM.ID == sname),]
  
  dat.cca@meta.data$takara_clonotype.common_biopsy_til[dat.cca@meta.data$UM.ID == sname & 
    (dat.cca@meta.data$takara_clonotype_TRA_1 %in% vdj.takara.clonotypes.common$clonotype |
    dat.cca@meta.data$takara_clonotype_TRA_2 %in% vdj.takara.clonotypes.common$clonotype | 
    dat.cca@meta.data$takara_clonotype_TRB_1 %in% vdj.takara.clonotypes.common$clonotype | 
    dat.cca@meta.data$takara_clonotype_TRB_2 %in% vdj.takara.clonotypes.common$clonotype)] <- T
  
  vdj.takara.sname.til <- vdj.takara[vdj.takara$clonotype %in% 
               c(dat.cca@meta.data$takara_clonotype_TRA_1[dat.cca@meta.data$Sample.ID == sname.til],
                 dat.cca@meta.data$takara_clonotype_TRA_2[dat.cca@meta.data$Sample.ID == sname.til],
                 dat.cca@meta.data$takara_clonotype_TRB_1[dat.cca@meta.data$Sample.ID == sname.til],
                 dat.cca@meta.data$takara_clonotype_TRA_2[dat.cca@meta.data$Sample.ID == sname.til]) & 
               vdj.takara$UM.ID == sname,]
  
  vdj.takara.sname.biopsy <- vdj.takara[vdj.takara$clonotype %in% 
                c(dat.cca@meta.data$takara_clonotype_TRA_1[dat.cca@meta.data$Sample.ID == sname.biopsy],
                  dat.cca@meta.data$takara_clonotype_TRA_2[dat.cca@meta.data$Sample.ID == sname.biopsy],
                  dat.cca@meta.data$takara_clonotype_TRB_1[dat.cca@meta.data$Sample.ID == sname.biopsy],
                  dat.cca@meta.data$takara_clonotype_TRA_2[dat.cca@meta.data$Sample.ID == sname.biopsy]) & 
                  vdj.takara$UM.ID == sname,]
  
  vdj.takara.sname.til <- vdj.takara.sname.til[order(vdj.takara.sname.til$vasu_read_count,decreasing = T),]
  vdj.takara.sname.biopsy <- vdj.takara.sname.biopsy[order(vdj.takara.sname.biopsy$vasu_read_count,decreasing = T),]
  
  vdj.takara.sname <- vdj.takara[vdj.takara$UM.ID==sname,]
  vdj.takara.sname$in_til <- vdj.takara.sname$clonotype %in% vdj.takara.sname.til$clonotype
  vdj.takara.sname$in_biopsy <- vdj.takara.sname$clonotype %in% vdj.takara.sname.biopsy$clonotype
  vdj.takara.sname$in_both <- vdj.takara.sname$in_til & vdj.takara.sname$in_biopsy
  
  vdj.takara.sname$represented <- ifelse(vdj.takara.sname$in_both,"both","")
  vdj.takara.sname$represented[vdj.takara.sname$represented==""] <- ifelse(vdj.takara.sname$in_til[vdj.takara.sname$represented==""],"til","")
  vdj.takara.sname$represented[vdj.takara.sname$represented==""] <- ifelse(vdj.takara.sname$in_biopsy[vdj.takara.sname$represented==""],"biopsy","")
  vdj.takara.sname$represented[vdj.takara.sname$represented==""] <- "none"
  
  stopifnot(identical(colnames(vdj.takara),colnames(vdj.takara.sname)))
  
  vdj.takara[vdj.takara$UM.ID==sname,] <- vdj.takara.sname
  
  vdj.takara[vdj.takara$UM.ID==sname,]$percent <- vdj.takara[vdj.takara$UM.ID==sname,]$vasu_read_count / 
    sum(vdj.takara[vdj.takara$UM.ID==sname,]$vasu_read_count)
}




ggplot(vdj.takara.sname, 
       aes(x = reorder(clonotype,-vasu_read_count), y = vasu_read_count, fill = represented)) + 
  geom_bar(stat="identity") + theme(axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank())



ggplot(vdj.takara[vdj.takara$UM.ID!="",], 
       aes(x = reorder(clonotype,-vasu_read_count), y = vasu_read_count, fill = represented)) + 
  geom_bar(stat="identity") + facet_wrap(~ UM.ID, scales = "free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + geom_hline(yintercept=10)

ggplot(vdj.takara[vdj.takara$UM.ID!="",], 
       aes(x = reorder(clonotype,-percent), y = percent, fill = represented)) + 
  geom_bar(stat="identity") + facet_wrap(~ UM.ID, scales = "free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + geom_hline(yintercept=0.1)

vdj.takara$shared <- vdj.takara$clonotype %in% shared
ggplot(vdj.takara[vdj.takara$UM.ID %in% c("UM46","UM22"),], 
       aes(x = reorder(clonotype,-percent), y = percent, fill = shared)) + 
  geom_bar(stat="identity") + facet_wrap(~ UM.ID, scales = "free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + geom_hline(yintercept=0.1)



vdj.takara$pass_threshold <- F
vdj.takara$pass_threshold <- vdj.takara$percent > 0.01

#vdj.takara[vdj.takara$UM.ID!="" & vdj.takara$represented!="none",][
#  vdj.takara[vdj.takara$UM.ID!="" & vdj.takara$represented!="none",]$pass_threshold,]

dat.cca@meta.data$pass_threshold <- NULL

x <- merge(dat.cca@meta.data,vdj.takara[,c("clonotype","pass_threshold")],
           by.x="takara_clonotype",by.y="clonotype",all.x=T,all.y=F)

pass_threshold <- setNames(x$pass_threshold,x$cell_id)
pass_threshold <- pass_threshold[rownames(dat.cca@meta.data)]
identical(names(pass_threshold),rownames(dat.cca@meta.data))
pass_threshold[is.na(pass_threshold)] <- F

dat.cca <- AddMetaData(dat.cca, metadata = pass_threshold, col.name = "pass_threshold")

#saveRDS(dat.cca,paste0(outdir,"dat.cca.modified.rda"))
#dat.cca <- readRDS(paste0(outdir,"dat.cca.modified.rda"))

# Which clonotypes are 41bb+ in biopsy data? -----------------------------------

table(dat.cca@assays$RNA@counts["TNFRSF9",] > 0 & dat.cca@meta.data$takara_clonotype != "" ,
      dat.cca@meta.data$Sample.ID)

dat.cca@meta.data$is_41bb_pos_biopsy <- dat.cca@assays$RNA@counts["TNFRSF9",] > 0 & 
  dat.cca@meta.data$takara_clonotype != "" & dat.cca@meta.data$project.name == "GWA-JN-388"

#dat.cca@meta.data$is_41bb_pos_biopsy <- dat.cca@meta.data$takara_clonotype %in% 
#  dat.cca@meta.data$takara_clonotype[dat.cca@meta.data$is_41bb_pos_biopsy]

dat.cca@meta.data$is_41bb_pos_biopsy <- 
  (dat.cca@meta.data$takara_clonotype_TRA_1 != "" & 
  dat.cca@meta.data$takara_clonotype_TRA_1 %in% c(
    dat.cca@meta.data$takara_clonotype_TRA_1[dat.cca@meta.data$is_41bb_pos_biopsy],
    dat.cca@meta.data$takara_clonotype_TRA_2[dat.cca@meta.data$is_41bb_pos_biopsy])) | 
  (dat.cca@meta.data$takara_clonotype_TRA_2 != "" & 
  dat.cca@meta.data$takara_clonotype_TRA_2 %in% c(
    dat.cca@meta.data$takara_clonotype_TRA_1[dat.cca@meta.data$is_41bb_pos_biopsy],
    dat.cca@meta.data$takara_clonotype_TRA_2[dat.cca@meta.data$is_41bb_pos_biopsy])) | 
  (dat.cca@meta.data$takara_clonotype_TRB_1 != "" & 
  dat.cca@meta.data$takara_clonotype_TRB_1 %in% c(
    dat.cca@meta.data$takara_clonotype_TRB_1[dat.cca@meta.data$is_41bb_pos_biopsy],
    dat.cca@meta.data$takara_clonotype_TRB_2[dat.cca@meta.data$is_41bb_pos_biopsy])) | 
  (dat.cca@meta.data$takara_clonotype_TRB_2 != "" & 
  dat.cca@meta.data$takara_clonotype_TRB_2 %in% c(
    dat.cca@meta.data$takara_clonotype_TRB_1[dat.cca@meta.data$is_41bb_pos_biopsy],
    dat.cca@meta.data$takara_clonotype_TRB_2[dat.cca@meta.data$is_41bb_pos_biopsy]))

dat.cca@meta.data$pass_threshold_2 <- dat.cca@meta.data$pass_threshold | dat.cca@meta.data$is_41bb_pos_biopsy

table(dat.cca@meta.data$pass_threshold_2,dat.cca@meta.data$Sample.ID)
table(dat.cca@meta.data$pass_threshold,dat.cca@meta.data$Sample.ID)

length(unique(dat.cca@meta.data$takara_clonotype[dat.cca@meta.data$pass_threshold_2]))
length(unique(dat.cca@meta.data$takara_clonotype[dat.cca@meta.data$pass_threshold]))

vdj.takara$pass_threshold_2 <- F
vdj.takara$pass_threshold_2[vdj.takara$clonotype %in% c(
  dat.cca@meta.data[dat.cca@meta.data$pass_threshold_2,]$takara_clonotype_TRA_1,
  dat.cca@meta.data[dat.cca@meta.data$pass_threshold_2,]$takara_clonotype_TRA_2,
  dat.cca@meta.data[dat.cca@meta.data$pass_threshold_2,]$takara_clonotype_TRB_1,
  dat.cca@meta.data[dat.cca@meta.data$pass_threshold_2,]$takara_clonotype_TRB_2)] <- T

# Implement test for all samples -----------------------------------------------
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
  # idx_not_mart1 <- (! rownames(x@meta.data) %in% nms_shared) & (x@meta.data$vasu_phenotype!=pheno) & 
  #   !is.na(x@meta.data$cdr3) & x@meta.data$takara_clonotype==""
  idx_not_mart1 <- (! rownames(x@meta.data) %in% nms_shared) & (x@meta.data$vasu_phenotype!=pheno) & 
    (!is.na(x@meta.data$cdr3) | ((x@assays$RNA@counts["CD3D",] > 0) | (x@assays$RNA@counts["CD3G",] > 0))) & 
    x@meta.data$takara_clonotype==""
  
   #idx_cells <- grepl("^CD8 T cells",x@meta.data$active.ident.updated_2) & 
  #   (idx_mart1 | idx_not_mart1) & x@meta.data$vasu_phenotype != "Unidentified"
  
  idx_cells <- (idx_mart1 | idx_not_mart1) & x@meta.data$vasu_phenotype != "Unidentified"
  
  idx_cells <- idx_cells &
    ((x@assays$RNA@counts["CD8A",] > 0) | (x@assays$RNA@counts["CD8B",] > 0)) &
    (x@assays$RNA@counts["CD4",] == 0) & (x@assays$RNA@counts["CD14",] == 0) & 
    (x@assays$RNA@counts["CD19",] == 0) & (x@assays$RNA@counts["MLANA",] == 0) & 
    (x@assays$RNA@counts["PMEL",] == 0) & (x@assays$RNA@counts["TYR",] == 0)
  
  # ((x@assays$RNA@counts["CD3D",] > 0) | (x@assays$RNA@counts["CD3G",] > 0)) &
  
  idx_cells <- idx_cells & (x@meta.data$pass_threshold_2 | x@meta.data$vasu_phenotype == "None")
  
  idx_cells <- idx_cells & (! x@meta.data$takara_clonotype %in% 
                              vdj.takara$clonotype[!vdj.takara$pass_threshold_2 & 
                                                     vdj.takara$UM.ID == unique(x@meta.data$UM.ID)])
  
  #idx_cells <- idx_cells & 
  #  x@meta.data$percent_mito <= summary(x@meta.data[idx_cells,]$percent_mito)[["3rd Qu."]] & 
  #  x@meta.data$percent_mito >= summary(x@meta.data[idx_cells,]$percent_mito)[["1st Qu."]]
  
  
  cl <- dat.cca@meta.data$cdr3[dat.cca@assays$RNA@counts["TNFRSF9",] > 0 & dat.cca@meta.data$project.name=="GWA-JN-388"]
  cl <- cl[!is.na(cl)]
  any(cl=="")
  any(cl=="NA")
  cl <- unique(unlist(strsplit(cl,";")))
  
  cl2 <- dat.cca@meta.data$cdr3[which(dat.cca@meta.data$vasu_phenotype=="None")]
  nms <- rownames(dat.cca@meta.data)[which(dat.cca@meta.data$vasu_phenotype=="None")]
  cl2 <- setNames(cl2,nms)
  cl2 <- cl2[which(!is.na(cl2))]
  cl2 <- strsplit(cl2,";")
  nms_shared <- names(which(unlist(lapply(cl2,function(x){any(x %in% cl)}))==T))
  
  idx_cells <- idx_cells & ! colnames(x) %in% nms_shared
  
  cl <- vdj.takara$cdr3
  cl2 <- dat.cca@meta.data$cdr3[which(dat.cca@meta.data$vasu_phenotype=="None")]
  nms <- rownames(dat.cca@meta.data)[which(dat.cca@meta.data$vasu_phenotype=="None")]
  cl2 <- setNames(cl2,nms)
  cl2 <- cl2[which(!is.na(cl2))]
  cl2 <- strsplit(cl2,";")
  nms_shared <- names(which(unlist(lapply(cl2,function(x){any(x %in% cl)}))==T))
  
  idx_cells <- idx_cells & ! colnames(x) %in% nms_shared
  
  return(idx_cells)
}

reprocess_sample <- function(x){
  y <- x
  DefaultAssay(y) <- "RNA"
  y <- FindVariableFeatures(y)
  y <- ScaleData(y,vars.to.regress = c("S.Score","G2M.Score","nCount_RNA"))
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

filter_genes <- function(x,idx_cells){
  dat <- x@assays$RNA@data[grep("\\.",rownames(x@assays$RNA@data),invert = T),idx_cells]
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

prepare_data <- function(dat.cca.split, sname, pheno){
  x <- dat.cca.split[[sname]];
  y <- reprocess_sample(x)
  
  out <- list(x=x,
              y=y)
  
  return(out)
}

test_gene_wx <- function(s,gene,pheno){
  if (sum(s@meta.data$vasu_phenotype == pheno) > 2 && 
      sum(s@meta.data$vasu_phenotype == "None")){
    
    pct.41bb <- sum(s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == pheno] > 0) / 
      sum(s@meta.data$vasu_phenotype == pheno)
    pct.none <- sum(s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "None"] > 0) / 
      sum(s@meta.data$vasu_phenotype == "None")
    
    if (pct.41bb > 0.2 | pct.none > 0.2){
    
      p.vasu_phenotype <- wilcox.test(
        s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == pheno],
        s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "None"],exact = F)$p.value
      
      # p.g1s <- wilcox.test(s@assays$RNA@data[gene,][s@meta.data$Phase == "G1"],
      #             s@assays$RNA@data[gene,][s@meta.data$Phase == "S"],exact = F)$p.value
      # 
      # p.g1g2m <- wilcox.test(s@assays$RNA@data[gene,][s@meta.data$Phase == "G1"],
      #             s@assays$RNA@data[gene,][s@meta.data$Phase == "G2M"],exact = F)$p.value
      # 
      # p.sg2m <- wilcox.test(s@assays$RNA@data[gene,][s@meta.data$Phase == "S"],
      #             s@assays$RNA@data[gene,][s@meta.data$Phase == "G2M"],exact = F)$p.value
      # 
      # p.rna <- cor.test(s@assays$RNA@data[gene,],s@meta.data$nCount_RNA)$p.value
      # p.mito <- cor.test(s@assays$RNA@data[gene,],s@meta.data$percent_mito)$p.value
      # p.ribo <- cor.test(s@assays$RNA@data[gene,],s@meta.data$percent_ribo)$p.value
      
      #s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "CD3+41bb+"]
      #s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "None"]
      pct.41bb <- sum(s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == pheno] > 0) / 
        sum(s@meta.data$vasu_phenotype == pheno)
      pct.none <- sum(s@assays$RNA@data[gene,][s@meta.data$vasu_phenotype == "None"] > 0) / 
        sum(s@meta.data$vasu_phenotype == "None")
      
      mean.nonzero.41bb <- mean(s@assays$RNA@data[gene,][
        s@meta.data$vasu_phenotype == pheno & s@assays$RNA@data[gene,]>0])
      mean.nonzero.none <- mean(s@assays$RNA@data[gene,][
        s@meta.data$vasu_phenotype == "None" & s@assays$RNA@data[gene,]>0])
      
      res <- data.frame(p.vasu_phenotype = p.vasu_phenotype,
                        pct.41bb = pct.41bb,
                        pct.none = pct.none,
                        mean.nonzero.41bb = mean.nonzero.41bb,
                        mean.nonzero.none = mean.nonzero.none)
      } else {
        res <- data.frame(p.vasu_phenotype = NA,
                          pct.41bb = NA,
                          pct.none = NA,
                          mean.nonzero.41bb = NA,
                          mean.nonzero.none = NA)
      }
    } else {
      res <- data.frame(p.vasu_phenotype = NA,
                        pct.41bb = NA,
                        pct.none = NA,
                        mean.nonzero.41bb = NA,
                        mean.nonzero.none = NA)
    }
  return(res)
}

# Prepare sample data
get_subset <- function(sname,out,pheno){
  y <- out[[sname]]$y
  
  idx_cells <- filter_cells(y, pheno)
  idx_cells <- idx_cells & (!is.na(y@meta.data$cdr3) | 
                              ((y@assays$RNA@counts["CD3D",] > 0) | (y@assays$RNA@counts["CD3G",] > 0)))
  idx_cells <- idx_cells & y@meta.data$vasu_phenotype %in% c(pheno,"None")
  
  x <- dat.cca.split[[sname]]
  dat <- filter_genes(x,idx_cells)
  genes_keep <- rownames(dat)
  
  s <- subset(y,cells=names(idx_cells)[idx_cells])
  return(list(s = s,
              genes_keep = genes_keep))
}

# Test (Wilcoxon)
run_test_gene_wx <- function(s,genes_keep,pheno){
  rs <- lapply(genes_keep,function(x){test_gene_wx(s,x,pheno=pheno)})
  rs <- do.call("rbind",rs)
  rownames(rs) <- genes_keep
  rs$q.vasu_phenotype <- p.adjust(rs$p.vasu_phenotype,method = "BH")
  return(rs)
}

dat.cca.split <- SplitObject(dat.cca, split.by = "Sample.ID")
saveRDS(dat.cca.split,file = paste0(outdir,"dat.cca.split.rda"))
#dat.cca.split <- readRDS(file = paste0(outdir,"dat.cca.split.rda"))

#snames.mart1 <- sample_map[sample_map$in_takara & sample_map$vasu_phenotype == "MART1",]$Var2
snames.41bb <- sample_map[sample_map$in_takara & sample_map$vasu_phenotype == "CD3+41bb+",]$Var2

# out.mart1 <- list()
# for (sname in snames.mart1){
#   out.mart1[[sname]] <- prepare_data(dat.cca.split = dat.cca.split, 
#                                      sname = sname, 
#                                      pheno = "MART1")
# }

out.41bb <- list()
for (sname in snames.41bb){
  out.41bb[[sname]] <- prepare_data(dat.cca.split = dat.cca.split, 
                                    sname = sname, 
                                    pheno = "CD3+41bb+")
}

#s.mart1 <- lapply(snames.mart1,function(x){get_subset(x,out.mart1,pheno="MART1")})
s.41bb <- lapply(grep("^Sample",snames.41bb,value=T),function(x){get_subset(x,out.41bb,pheno="CD3+41bb+")})
#names(s.mart1) <- snames.mart1
names(s.41bb) <- grep("^Sample",snames.41bb,value=T)

#tst.wx.mart1 <- lapply(s.mart1,function(x){run_test_gene_wx(x$s,x$genes_keep,pheno="MART1")})
tst.wx.41bb <- lapply(s.41bb,function(x){run_test_gene_wx(x$s,x$genes_keep,pheno="CD3+41bb+")})
#names(tst.wx.mart1) <- snames.mart1
names(tst.wx.41bb) <- grep("^Sample",snames.41bb,value=T)

filter_fun <- function(rs){
  # rs.filtered <- rs[which(rs$p.vasu_phenotype < 0.05 & rs$p.rna >= 0.05 & rs$p.mito >= 0.05 & 
  #                           rs$p.ribo >= 0.05 & rs$p.g1s >= 0.05 & rs$p.g1g2m >= 0.05 & 
  #                           rs$p.sg2m >= 0.05 & rs$q.vasu_phenotype < 0.05 & 
  #                           rs$pct.41bb > rs$pct.none),]
  #rs.filtered <- rs[which(rs$p.vasu_phenotype < 0.05 & rs$q.vasu_phenotype < 0.05 & 
  #                          rs$pct.41bb > rs$pct.none & 
  #                          abs(rs$pct.41bb - rs$pct.none) > 0.3),]
  rs.filtered <- rs[which(rs$p.vasu_phenotype < 0.05 &
                            rs$pct.41bb > rs$pct.none & rs$pct.41bb > 0.1 & 
                            rs$mean.nonzero.41bb > rs$mean.nonzero.none),]
  #rs.filtered <- rs[which(rs$q.vasu_phenotype < 0.05 &
  #                          rs$pct.41bb > rs$pct.none & rs$pct.41bb > 0.1),]
  rs.filtered <- rs[which(rs$q.vasu_phenotype < 0.05 &
                            rs$pct.41bb > rs$pct.none & 
                            rs$mean.nonzero.41bb > rs$mean.nonzero.none),]
  return(rs.filtered)
}

lapply(tst.wx.41bb,filter_fun)

write.table(sort(names(which(table(unlist(
  lapply(lapply(tst.wx.41bb,filter_fun),rownames)))>=1))),
  quote = F,row.names = F,col.names = F)

write.table(sort(names(which(table(unlist(
  lapply(lapply(tst.wx.41bb[grep("Sample",names(tst.wx.41bb))],filter_fun),rownames)))>=2))),
  quote = F,row.names = F,col.names = F)




genes.sig <- sort(names(which(table(unlist(
  lapply(lapply(tst.wx.41bb[grep("Sample",names(tst.wx.41bb))],filter_fun),rownames)))>=1)))[1:10]

genes.sig <- rownames(lapply(tst.wx.41bb,filter_fun)[["SampleID_9_11june18"]])

FeaturePlot(s.41bb[["SampleID_9_11june18"]]$s,
            features = c(genes.sig[1:10],
              "nCount_RNA","TRGV10","ICOS","HAVCR2","TNFRSF9"),order = T) + 
  DimPlot(s.41bb[["SampleID_9_11june18"]]$s,group.by = "vasu_phenotype") + 
  DimPlot(s.41bb[["SampleID_9_11june18"]]$s,split.by = "vasu_phenotype",group.by = "vasu_phenotype") + 
  DimPlot(s.41bb[["SampleID_9_11june18"]]$s,group.by = "Phase")


subs <- subset(s.41bb[["SampleID_9_11june18"]]$s,vasu_phenotype=="CD3+41bb+")
FeaturePlot(subs,
            features = c(genes.sig[1:10],
                         "nCount_RNA","TRGV10","ICOS","HAVCR2","TNFRSF9"),order = T) + 
  DimPlot(subs,group.by = "vasu_phenotype") + 
  DimPlot(subs,group.by = "Phase")

subs <- subset(s.41bb[["SampleID_9_11june18"]]$s,vasu_phenotype=="None")
subs_none <- subset(subs,pass_threshold_2==F)
FeaturePlot(subs_none,
            features = c(genes.sig[1:10],
                         "nCount_RNA","TRGV10","ICOS","HAVCR2","TNFRSF9"),order = T) + 
  DimPlot(subs_none,group.by = "vasu_phenotype") + 
  DimPlot(subs_none,group.by = "Phase")


subs <- subset(s.41bb[["SampleID_9_11june18"]]$s,pass_threshold==T)
subs_pass <- subset(subs,pass_threshold==T)
FeaturePlot(subs,
            features = c(genes.sig[1:10],
                         "nCount_RNA","TRGV10","ICOS","HAVCR2","TNFRSF9"),order = T) + 
  DimPlot(subs,group.by = "vasu_phenotype") + 
  DimPlot(subs,group.by = "Phase")


subs_keep <- subset(s.41bb[["SampleID_9_11june18"]]$s,cells = c(colnames(subs_none),colnames(subs_pass)))
subs_keep <- subset(subs_keep,nCount_RNA > 1000)
subs_keep <- subset(subs_keep,nFeature_RNA > 1000)
subs_keep <- reprocess_sample(subs_keep)

DimPlot(subs_keep)

FeaturePlot(subs_keep,
            features = c(genes.sig[1:10],
                         "nCount_RNA","TRGV10","ICOS","HAVCR2","TNFRSF9","CRIP1"),order = T) + 
  DimPlot(subs_keep,group.by = "vasu_phenotype") + 
  DimPlot(subs_keep,group.by = "vasu_phenotype", split.by = "vasu_phenotype") + 
  DimPlot(subs_keep,group.by = "Phase")

FeaturePlot(subs_keep,
            features = c("TRBV30","TRAV41","TRGV10","TNFRSF9","ICOS"),split.by = "vasu_phenotype",order = T)

FeaturePlot(subset(subs_keep,is_41bb_pos_biopsy==F),
            features = c("TRBV30","TRAV41","TRGV10","TNFRSF9","ICOS"),split.by = "vasu_phenotype",order = T)





mark <- FindMarkers(subs_keep,group.by = "vasu_phenotype",ident.1 = "CD3+41bb+",ident.2 = "None")
mark <- mark[order(mark$avg_log2FC,decreasing = T),]
head(mark[which(mark$p_val_adj < 0.1),],n=100)

mark7 <- FindMarkers(s.41bb[["SampleID_7_11june18"]]$s,
                    group.by = "vasu_phenotype",ident.1 = "CD3+41bb+",ident.2 = "None",test.use = "roc")
mark7 <- mark7[order(mark7$avg_log2FC,decreasing = T),]
head(mark7[which(mark7$p_val_adj < 0.1),],n=100)

mark7 <- mark7[order(mark7$avg_diff,decreasing = T),]
head(mark7[mark7$myAUC > 0.5,])

mark9 <- FindMarkers(s.41bb[["SampleID_9_11june18"]]$s,
                     group.by = "vasu_phenotype",ident.1 = "CD3+41bb+",ident.2 = "None",test.use = "roc")
mark9 <- mark9[order(mark9$avg_log2FC,decreasing = T),]
head(mark9[which(mark9$p_val_adj < 0.1),],n=100)

mark9 <- mark9[order(mark9$myAUC,decreasing = T),]
head(mark9[which(mark9$myAUC > 0.5),],n=50)


intersect(rownames(mark7[mark7$myAUC > 0.5 & mark7$avg_diff > 0,]),
          rownames(mark9[mark9$myAUC > 0.5 & mark9$avg_diff > 0,]))

intersect(rownames(mark7[mark7$myAUC > 0.5,]),
          rownames(mark9[mark9$myAUC > 0.5,]))


FeaturePlot(s.41bb[["SampleID_9_11june18"]]$s,
            features = names(tab.cca)[idx][1:5],
            split.by = "vasu_phenotype",order = T)



table(dat.cca@assays$RNA@counts["CD4",]>0,dat.cca@meta.data$vasu_phenotype)
table(dat.cca@assays$RNA@counts["CD8A",]>0,dat.cca@meta.data$vasu_phenotype)
table(dat.cca@assays$RNA@counts["CD8B",]>0,dat.cca@meta.data$vasu_phenotype)
table(dat.cca@assays$RNA@counts["TNFRSF9",]>0,dat.cca@meta.data$vasu_phenotype)
table(dat.cca@assays$RNA@counts["TRGV10",]>0,dat.cca@meta.data$vasu_phenotype)

genes_test <- names(which(rowSums(
  dat.cca@assays$RNA@counts[,dat.cca@meta.data$vasu_phenotype %in% c("CD3+41bb+","None")] > 0) > 100))
length(genes_test)


idx <- subs_keep@meta.data$Sample.ID %in% c("SampleID_7_11june18","SampleID_9_11june18") & 
  subs_keep@meta.data$vasu_phenotype %in% c("CD3+41bb+","None")

tab <- list()
for (gene in genes_test){
  tab[[gene]] <- table(subs_keep@assays$RNA@counts[gene,idx]>0, subs_keep@meta.data$vasu_phenotype[idx])
}


idx <- dat.cca@meta.data$Sample.ID %in% c("SampleID_7_11june18","SampleID_9_11june18") & 
  dat.cca@meta.data$vasu_phenotype %in% c("CD3+41bb+","None")

tab.cca <- list()
for (gene in genes_test){
  tab.cca[[gene]] <- table(dat.cca@assays$RNA@counts[gene,idx]>0,dat.cca@meta.data$vasu_phenotype[idx])
}

idx <- which(unlist(lapply(tab.cca,function(x){
  x["TRUE","CD3+41bb+"]/x["FALSE","CD3+41bb+"] > 
    5*x["TRUE","None"]/x["FALSE","None"] & 
    x["TRUE","MART1"]/x["FALSE","MART1"] > 
    5*x["TRUE","None"]/x["FALSE","None"]
})))


names(tab.cca)[idx]











FeaturePlot(s.41bb[["SampleID_7_11june18"]]$s,
            features = c(sort(names(which(table(unlist(
              lapply(lapply(tst.wx.41bb[grep("Sample",names(tst.wx.41bb))],filter_fun),rownames)))>=2)))[1:10],
              "nCount_RNA","TRGV10","ICOS","HAVCR2","TNFRSF9"),order = T) + 
  DimPlot(s.41bb[["SampleID_7_11june18"]]$s,group.by = "vasu_phenotype") + 
  DimPlot(s.41bb[["SampleID_7_11june18"]]$s,split.by = "vasu_phenotype",group.by = "vasu_phenotype") + 
  DimPlot(s.41bb[["SampleID_7_11june18"]]$s,group.by = "Phase")




FeaturePlot(s.mart1[["A2-GEX"]]$s,
            features = sort(names(which(table(unlist(
              lapply(lapply(tst.wx.mart1,filter_fun),rownames)))>=1)))[1:10],order = T) + 
  DimPlot(s.mart1[["A2-GEX"]]$s,group.by = "vasu_phenotype")

# Exlude cells that are in the same cluster as the target (41bb) clonotypes?
# Subcluster to better understand sources of gene expression?

# From liger -----------------------

#cells <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/cells.rda")
#table(dat.cca@meta.data[cells,]$vasu_phenotype)/table(dat.cca@meta.data$vasu_phenotype)

# Calculate vasu_phenotype percentage per liger cluster
lig.tumor_clusters <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/lig.tumor_clusters.rda")
lig.tumor_clusters <- setNames(as.character(lig.tumor_clusters),names(lig.tumor_clusters))

#celltypes3 <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/celltypes3.rda")
#celltypes3 <- setNames(as.character(celltypes3),names(celltypes3))

pct <- list()
clus <- unique(lig.tumor_clusters)
for (cl in clus){
  cells <- names(lig.tumor_clusters[lig.tumor_clusters==cl])
  tmp <- table(dat.cca@meta.data[cells,]$vasu_phenotype)
  tmp2 <- table(dat.cca@meta.data$vasu_phenotype)
  missing <- setdiff(names(tmp2),names(tmp))
  if (length(missing)>0){
    tmp <- c(tmp,rep(0,length(missing)))
    names(tmp)[names(tmp)==""] <- missing
  }
  tmp <- tmp[order(names(tmp))]
  tmp2 <- tmp2[order(names(tmp2))]
  pct[[cl]] <- print(tmp/tmp2)
}
pct <- do.call("rbind",pct)

library(pheatmap)
#pheatmap(pct)
#pheatmap(scale(pct,center = T,scale = T))
#pheatmap(t(scale(t(pct),center = T,scale = T)))
#pheatmap(table(lig.tumor_clusters,celltypes3))

tcell_clusters <- as.character(c(21,25,47,50,2,4,11,30,26,34))

#pheatmap(pct[tcell_clusters,])
#pheatmap(scale(pct[tcell_clusters,]))
#pheatmap(t(scale(t(pct[tcell_clusters,]))))
#pheatmap(cor(scale(pct[tcell_clusters,]),t(scale(t(pct[tcell_clusters,])))))
#pheatmap(cor(t(scale(pct[tcell_clusters,])),t(t(scale(t(pct[tcell_clusters,]))))))

vasu_phenotype <- setNames(dat.cca@meta.data$vasu_phenotype,rownames(dat.cca@meta.data))
saveRDS(vasu_phenotype,"~/proj/um_ss/Investigations/seurat/results/liger_all/vasu_phenotype.rda")

dat.cca <- AddMetaData(dat.cca,metadata = lig.tumor_clusters,col.name = "lig.tumor_clusters")
table(is.na(dat.cca@meta.data$lig.tumor_clusters))
dim(dat.cca@meta.data)
dat.cca@meta.data$lig.tumor_clusters[is.na(dat.cca@meta.data$lig.tumor_clusters)] <- "X"
dat.cca.lig <- subset(dat.cca,lig.tumor_clusters!="X")
table(dat.cca.lig@meta.data$lig.tumor_clusters)

dat.cca.lig.tcells <- subset(dat.cca.lig, lig.tumor_clusters %in% tcell_clusters)
table(dat.cca.lig.tcells@meta.data$lig.tumor_clusters)
dat.cca.lig.tcells <- subset(dat.cca.lig.tcells,cells=names(which(dat.cca.lig.tcells@assays$RNA@counts["CD4",] == 0)))
table(dat.cca.lig.tcells@meta.data$lig.tumor_clusters)

dat.cca.lig.tcells.re <- reprocess_sample(dat.cca.lig.tcells)
DimPlot(dat.cca.lig.tcells.re) + DimPlot(dat.cca.lig.tcells.re,group.by = "lig.tumor_clusters",label = T)

FeaturePlot(dat.cca.lig.tcells.re,
            features = c("CD3D","CD8A","CD8B","NCAM1","TRGV9","TNFRSF9",
                         "ICOS","HAVCR2","FOXP3","CTLA4","TRGV10","PDCD1",
                         "LAG3","nCount_RNA"),order = T) + 
  DimPlot(dat.cca.lig.tcells.re,group.by = "vasu_phenotype") + 
  DimPlot(dat.cca.lig.tcells.re,group.by = "Phase")

mark <- FindMarkers(dat.cca.lig.tcells.re, 
                    group.by = "lig.tumor_clusters",ident.1 = c("2","4","50"),
                    ident.2 = c("21","25","47"))

mark.biopsy <- FindMarkers(dat.cca.lig.tcells.re, 
                    group.by = "lig.tumor_clusters",ident.1 = c("2","4","11","26","34","50"),
                    ident.2 = c("21","25","30","47"))

head(mark.biopsy[mark.biopsy$p_val_adj < 0.05 & mark.biopsy$pct.1 > mark.biopsy$pct.2 & mark.biopsy$avg_log2FC > 0,], n = 100)
head(mark.biopsy[mark.biopsy$p_val_adj < 0.05 & mark.biopsy$pct.2 > mark.biopsy$pct.1 & mark.biopsy$avg_log2FC < 0,], n = 100)

# Same, but for TILs -----------------------------------------------------------

dat.cca.til <- subset(dat.cca, Sample.ID %in% c("SampleID_7_11june18","SampleID_9_11june18"))
dat.cca.til <- reprocess_sample(dat.cca.til)
DimPlot(dat.cca.til) + DimPlot(dat.cca.til, group.by = "Sample.ID")

FeaturePlot(dat.cca.til,
            features = c("CD3D","CD8A","CD8B","NCAM1","TRGV9","TNFRSF9",
                         "ICOS","HAVCR2","FOXP3","CTLA4","TRGV10","PDCD1",
                         "LAG3","nCount_RNA"),order = T) + 
  DimPlot(dat.cca.til,group.by = "vasu_phenotype") + 
  DimPlot(dat.cca.til,group.by = "Phase")

FeaturePlot(dat.cca.til,
            features = c("CD3D","CD3D","CD8A","NCAM1","CD4"),order = T) + DimPlot(dat.cca.til, label = T)

dat.cca.til@meta.data$clus <- Idents(dat.cca.til)

dat.cca.til <- RenameIdents(dat.cca.til,
             `0` = "CD8 T cells",
             `1` = "CD4 T cells",
             `2` = "CD8 T cells",
             `3` = "CD8 T cells",
             `4` = "CD8 T cells",
             `5` = "CD8 T cells",
             `6` = "CD8 T cells",
             `7` = "CD8 T cells",
             `8` = "CD8 T cells",
             `9` = "CD8 T cells",
             `10` = "CD8 T cells",
             `11` = "CD8 T cells",
             `12` = "CD8 T cells",
             `13` = "CD4 T cells",
             `14` = "CD8 T cells",
             `15` = "NK cells")

DimPlot(dat.cca.til, group.by = "clus",label = T) + 
  DimPlot(dat.cca.til) + FeaturePlot(dat.cca.til, features = "CD8A") + 
  DimPlot(dat.cca.til, group.by = "vasu_phenotype")

tcell_clusters <- c(0,2,3,4,5,6,7,8,9,10,11,12,14)

calculate_percentage <- function(x,clusters){
  pct <- list()
  clus <- unique(clusters)
  for (cl in clus){
    cells <- names(clusters[clusters==cl])
    tmp <- table(x@meta.data[cells,]$vasu_phenotype)
    tmp2 <- table(x@meta.data$vasu_phenotype)
    missing <- setdiff(names(tmp2),names(tmp))
    if (length(missing)>0){
      tmp <- c(tmp,rep(0,length(missing)))
      names(tmp)[names(tmp)==""] <- missing
    }
    tmp <- tmp[order(names(tmp))]
    tmp2 <- tmp2[order(names(tmp2))]
    pct[[cl]] <- tmp/tmp2
  }
  pct <- do.call("rbind",pct)
  return(pct)
}

pct.til <- calculate_percentage(dat.cca.til,
                                setNames(as.character(dat.cca.til@meta.data$clus),rownames(dat.cca.til@meta.data)))

pheatmap(pct.til)
pheatmap(scale(pct.til,center = T,scale = T))
pheatmap(t(scale(t(pct.til),center = T,scale = T)))

pheatmap(pct.til[as.character(tcell_clusters),])
pheatmap(scale(pct.til[as.character(tcell_clusters),],center = T,scale = T))
pheatmap(t(scale(t(pct.til[as.character(tcell_clusters),]),center = T,scale = T)))



clusters_exclude <- c("2","5","6")
setdiff(setdiff(as.character(tcell_clusters),c("3","8","9","11")),clusters_exclude)


mark.til <- FindMarkers(dat.cca.til, 
                    group.by = "clus",ident.1 = c("3","8","9","11"),
                    ident.2 =setdiff(setdiff(as.character(tcell_clusters),c("3","8","9","11")),clusters_exclude))

head(mark.til[mark.til$p_val_adj < 0.05 & mark.til$pct.1 > mark.til$pct.2 & mark.til$avg_log2FC > 0,], n = 100)
head(mark.til[mark.til$p_val_adj < 0.05 & mark.til$pct.2 > mark.til$pct.1 & mark.til$avg_log2FC < 0,], n = 100)

mark.til.7 <- FindMarkers(subset(dat.cca.til,Sample.ID == "SampleID_7_11june18"), 
                        group.by = "clus",ident.1 = c("3","8","9","11"),
                        ident.2 = setdiff(setdiff(as.character(tcell_clusters),c("3","8","9","11")),clusters_exclude))

mark.til.9 <- FindMarkers(subset(dat.cca.til,Sample.ID == "SampleID_9_11june18"), 
                          group.by = "clus",ident.1 = c("3","8","9","11"),
                          ident.2 = setdiff(setdiff(as.character(tcell_clusters),c("3","8","9","11")),clusters_exclude))

head(mark.til.7[mark.til.7$p_val_adj < 0.05 & mark.til.7$pct.1 > mark.til.7$pct.2 & mark.til.7$avg_log2FC > 0,], n = 100)
head(mark.til.9[mark.til.9$p_val_adj < 0.05 & mark.til.9$pct.1 > mark.til.9$pct.2 & mark.til.9$avg_log2FC > 0,], n = 100)

genes <- intersect(rownames(mark.til.7[mark.til.7$p_val_adj < 0.05 & mark.til.7$pct.1 > mark.til.7$pct.2 & mark.til.7$avg_log2FC > 0,]),
          rownames(mark.til.9[mark.til.9$p_val_adj < 0.05 & mark.til.9$pct.1 > mark.til.9$pct.2 & mark.til.9$avg_log2FC > 0,]))


FeaturePlot(dat.cca.til, features = genes,order = T) + 
  DimPlot(dat.cca.til) + 
  DimPlot(dat.cca.til, group.by = "vasu_phenotype")


genes <- rownames(mark.til[mark.til$p_val_adj < 0.05 & mark.til$pct.1 > mark.til$pct.2 & mark.til$avg_log2FC > 0,])[11:20]
FeaturePlot(dat.cca.til, features = genes,order = T) + DimPlot(dat.cca.til, group.by = "clus") + 
  DimPlot(dat.cca.til) + 
  DimPlot(dat.cca.til, group.by = "vasu_phenotype") + DimPlot(dat.cca.til, group.by = "vasu_phenotype",
                                                              split.by = "vasu_phenotype")

DimPlot(dat.cca.til, group.by = "vasu_phenotype",
        split.by = "vasu_phenotype")
