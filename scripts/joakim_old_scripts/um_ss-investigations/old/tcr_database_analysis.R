vdjdb <- read.table("~/proj/um_ss/Investigations/data/TCR_databases/VDJdb/SearchTable-2022-08-30 10 35 20.945.tsv",
                    sep="\t",header=T,fill=T)

mcpas_tcr <- read.table("~/proj/um_ss/Investigations/data/TCR_databases/McPAS-TCR/McPAS-TCR.csv",sep=",",header=T)

tcr3d_cancer <- read.table("~/proj/um_ss/Investigations/data/TCR_databases/TCR3d/cancer_targeting_tcrs.csv",sep=",",header=T)

tcr3d_virus <- read.table("~/proj/um_ss/Investigations/data/TCR_databases/TCR3d/virus_targeting-tcrs.csv",sep=",",header=T)

vdjdb.simplified <- vdjdb[,c("Gene","CDR3","V","J","Epitope.gene","Species","Epitope.species")]
vdjdb.simplified$CDR3.alpha <- ""
vdjdb.simplified$CDR3.beta <- ""
vdjdb.simplified$CDR3.alpha[vdjdb.simplified$Gene=="TRA"] <- vdjdb.simplified$CDR3[vdjdb.simplified$Gene=="TRA"]
vdjdb.simplified$CDR3.beta[vdjdb.simplified$Gene=="TRB"] <- vdjdb.simplified$CDR3[vdjdb.simplified$Gene=="TRB"]
vdjdb.simplified$CDR3 <- NULL
vdjdb.simplified$TRAV <- ""
vdjdb.simplified$TRAJ <- ""
vdjdb.simplified$TRBV <- ""
vdjdb.simplified$TRBJ <- ""
vdjdb.simplified$TRBD <- ""
vdjdb.simplified$TRAV[vdjdb.simplified$Gene=="TRA"] <- vdjdb.simplified$V[vdjdb.simplified$Gene=="TRA"]
vdjdb.simplified$TRAJ[vdjdb.simplified$Gene=="TRA"] <- vdjdb.simplified$J[vdjdb.simplified$Gene=="TRA"]
vdjdb.simplified$TRBV[vdjdb.simplified$Gene=="TRB"] <- vdjdb.simplified$V[vdjdb.simplified$Gene=="TRB"]
vdjdb.simplified$TRBJ[vdjdb.simplified$Gene=="TRB"] <- vdjdb.simplified$J[vdjdb.simplified$Gene=="TRB"]
vdjdb.simplified$V <- NULL
vdjdb.simplified$J <- NULL
vdjdb.simplified$Gene <- NULL
vdjdb.simplified <- vdjdb.simplified[,c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD",
                    "Epitope.gene","Epitope.species","Species")]
colnames(vdjdb.simplified) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD",
                                "Antigen","Antigen_species_or_pathology","Species")


mcpas_tcr.simplified <- mcpas_tcr[,c("CDR3.alpha.aa","CDR3.beta.aa","TRAV","TRAJ","TRBV","TRBD",
                                     "TRBJ","Antigen.protein","Species","Pathology")]
mcpas_tcr.simplified <- mcpas_tcr.simplified[,c("CDR3.alpha.aa","CDR3.beta.aa",
  "TRAV","TRAJ","TRBV","TRBJ","TRBD","Antigen.protein","Pathology","Species")]
colnames(mcpas_tcr.simplified) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD",
                                "Antigen","Antigen_species_or_pathology","Species")

tcr3d_cancer.simplified <- tcr3d_cancer[,c("CDR3.alpha.","CDR3.beta.","TRAV","TRBV","Antigen","Cancer.BR.Type")]
tcr3d_cancer.simplified$TRAJ <- ""
tcr3d_cancer.simplified$TRBJ <- ""
tcr3d_cancer.simplified$TRBD <- ""
tcr3d_cancer.simplified$Species <- ""
tcr3d_cancer.simplified <- tcr3d_cancer.simplified[,
  c("CDR3.alpha.","CDR3.beta.","TRAV","TRAJ","TRBV","TRBJ","TRBD","Antigen","Cancer.BR.Type","Species")]
colnames(tcr3d_cancer.simplified) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD",
                                       "Antigen","Antigen_species_or_pathology","Species")

tcr3d_virus.simplified <- tcr3d_virus[,c("CDR3.alpha.","CDR3.beta.","TRAV","TRBV","Antigen","Virus.BR.Type")]
tcr3d_virus.simplified$TRAJ <- ""
tcr3d_virus.simplified$TRBJ <- ""
tcr3d_virus.simplified$TRBD <- ""
tcr3d_virus.simplified$Species <- ""
tcr3d_virus.simplified <- tcr3d_virus.simplified[,c("CDR3.alpha.","CDR3.beta.","TRAV","TRAJ","TRBV","TRBJ","TRBD","Antigen","Virus.BR.Type","Species")]
colnames(tcr3d_virus.simplified) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD",
                                       "Antigen","Antigen_species_or_pathology","Species")

vdjdb.simplified$database <- "vdjdb"
mcpas_tcr.simplified$database <- "mcpas_tcr"
tcr3d_cancer.simplified$database <- "tcr3d_cancer"
tcr3d_virus.simplified$database <- "tcr3d_virus"

head(vdjdb.simplified)
head(mcpas_tcr.simplified)
head(tcr3d_cancer.simplified)
head(tcr3d_virus.simplified)

db <- rbind(vdjdb.simplified,
      mcpas_tcr.simplified,
      tcr3d_cancer.simplified,
      tcr3d_virus.simplified)
db[db==""] <- NA

db$TRBV <- gsub("^ ","",db$TRBV)
db$TRAV <- gsub(",$","",db$TRAV)

split_entries <- function(db,coln,delim=",",idx_split=NULL){
  if (is.null(idx_split)){
    idx <- grepl(delim,db[[coln]])
  } else {
    idx <- idx_split
  }
  if (length(which(idx))>0){
    db.tmp_1 <- db[!idx,]
    db.tmp_2 <- db[idx,]
    tmp.list <- list()
    for (i in 1:nrow(db.tmp_2)){
      entries <- unlist(strsplit(db.tmp_2[i,][[coln]],delim))
      entries <- gsub(" ","",entries)
      tmp <- db.tmp_2[i,]
      for (j in 1:(length(entries)-1)){
        tmp <- rbind(tmp,db.tmp_2[i,])
      }
      tmp[[coln]] <- entries
      tmp.list[[i]] <- tmp
    }
    db.new <- rbind(db.tmp_1,do.call("rbind",tmp.list))
  } else {
    db.new <- db
  }
  return(db.new)
}

dim(db)
db <- split_entries(db,coln = "TRAV")
db <- split_entries(db,coln = "TRAJ")
db <- split_entries(db,coln = "TRBV")
db <- split_entries(db,coln = "TRBJ")
db <- split_entries(db,coln = "TRBD")
dim(db)

library(stringr)
db$TRAV <- str_split_fixed(db$TRAV,"\\*",2)[,1]
db$TRAJ <- str_split_fixed(db$TRAJ,"\\*",2)[,1]
db$TRBV <- str_split_fixed(db$TRBV,"\\*",2)[,1]
db$TRBJ <- str_split_fixed(db$TRBJ,"\\*",2)[,1]
db$TRBD <- str_split_fixed(db$TRBD,"\\*",2)[,1]

db$TRAV <- str_split_fixed(db$TRAV,":",2)[,1]
db$TRAJ <- str_split_fixed(db$TRAJ,":",2)[,1]
db$TRBV <- str_split_fixed(db$TRBV,":",2)[,1]
db$TRBJ <- str_split_fixed(db$TRBJ,":",2)[,1]
db$TRBD <- str_split_fixed(db$TRBD,":",2)[,1]

db$TRAV <- str_split_fixed(db$TRAV,"\xa0",2)[,1]
db$TRAJ <- str_split_fixed(db$TRAJ,"\xa0",2)[,1]
db$TRBV <- str_split_fixed(db$TRBV,"\xa0",2)[,1]
db$TRBJ <- str_split_fixed(db$TRBJ,"\xa0",2)[,1]
db$TRBD <- str_split_fixed(db$TRBD,"\xa0",2)[,1]

db$TRAV <- str_split_fixed(db$TRAV," ",2)[,1]
db$TRAJ <- str_split_fixed(db$TRAJ," ",2)[,1]
db$TRBV <- str_split_fixed(db$TRBV," ",2)[,1]
db$TRBJ <- str_split_fixed(db$TRBJ," ",2)[,1]
db$TRBD <- str_split_fixed(db$TRBD," ",2)[,1]

db$TRAV <- str_split_fixed(db$TRAV,"\\.",2)[,1]
db$TRAJ <- str_split_fixed(db$TRAJ,"\\.",2)[,1]
db$TRBV <- str_split_fixed(db$TRBV,"\\.",2)[,1]
db$TRBJ <- str_split_fixed(db$TRBJ,"\\.",2)[,1]
db$TRBD <- str_split_fixed(db$TRBD,"\\.",2)[,1]

db$TRAV <- gsub("/$","",db$TRAV)
db$TRAJ <- gsub("/$","",db$TRAJ)
db$TRBV <- gsub("/$","",db$TRBV)
db$TRBJ <- gsub("/$","",db$TRBJ)
db$TRBD <- gsub("/$","",db$TRBD)

db$TRBD[which(db$TRBD=="unknown")] <- NA
db$TRBD[which(db$TRBD=="na")] <- NA
db[db==""] <- NA

db <- db[which(db$TRAJ!="CATSESSGQTYEQYF"),] # malformed entry in original csv

unique(db$TRAV)
unique(db$TRAJ)
unique(db$TRBV)

gsub("^TRBV","",unique(db$TRBV))

db$TRBV <- gsub("^TRBV","",db$TRBV)
idx_split <- grepl("/",db$TRBV) & !grepl("DV",db$TRBV)

tst <- split_entries(db,coln="TRBV",delim="/",idx_split=idx_split)
paste0("TRBV",tst$TRBV[!is.na(tst$TRBV)])

table(rowSums(table(db$CDR3.alpha,db$database) > 0))

head(db)

dat.tumor_til.vdj <- readRDS(file="~/proj/um_ss/Manuscript/dat.tumor_til.vdj.rda")

sname_1 <- 1

idx_1 <- unlist(lapply(dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3, function(x){
  all(x!="None")
}))

tmp_1.10x_cdr3 <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3[which(idx_1)]
tmp_1.10x_v_gene <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_v_gene[which(idx_1)]
tmp_1.10x_d_gene <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_d_gene[which(idx_1)]
tmp_1.10x_j_gene <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_j_gene[which(idx_1)]
tmp_1.10x_c_gene <- dat.tumor_til.vdj[[sname_1]]$vdj.10x_c_gene[which(idx_1)]

which_db_tra_cdr3_in_sname_1.10x_cdr3 <- lapply(tmp_1.10x_cdr3,function(x){
  which(db$CDR3.alpha %in% x)
})

which_db_trb_cdr3_in_sname_1.10x_cdr3 <- lapply(tmp_1.10x_cdr3,function(x){
  which(db$CDR3.beta %in% x)
})

which_db_tra_v_gene_in_sname_1.10x_v_gene <- lapply(tmp_1.10x_v_gene,function(x){
  which(db$TRAV %in% x)
})

which_db_trb_v_gene_in_sname_1.10x_v_gene <- lapply(tmp_1.10x_v_gene,function(x){
  which(db$TRBV %in% x)
})

which_db_trb_d_gene_in_sname_1.10x_d_gene <- lapply(tmp_1.10x_d_gene,function(x){
  which(db$TRBD %in% x)
})

which_db_tra_j_gene_in_sname_1.10x_j_gene <- lapply(tmp_1.10x_j_gene,function(x){
  which(db$TRAJ %in% x)
})

which_db_trb_j_gene_in_sname_1.10x_j_gene <- lapply(tmp_1.10x_j_gene,function(x){
  which(db$TRBJ %in% x)
})

#cells that match db, the matching ones among their own
# CDR3a
# CDR3b
# TRAV
# TRBV
# TRAJ
# TRBJ
# TRBD


cells_matching_tra_cdr3 <- names(which(lapply(which_db_tra_cdr3_in_sname_1.10x_cdr3,length) > 0))
cells_matching_tra_cdr3.seq <- list()
for (cell in cells_matching_tra_cdr3){
  idx <- which_db_tra_cdr3_in_sname_1.10x_cdr3[[cell]]
  cells_matching_tra_cdr3.seq[[cell]] <- unique(db$CDR3.alpha[idx])
}

cells_matching_trb_cdr3 <- names(which(lapply(which_db_trb_cdr3_in_sname_1.10x_cdr3,length) > 0))
cells_matching_trb_cdr3.seq <- list()
for (cell in cells_matching_trb_cdr3){
  idx <- which_db_trb_cdr3_in_sname_1.10x_cdr3[[cell]]
  cells_matching_trb_cdr3.seq[[cell]] <- unique(db$CDR3.beta[idx])
}

cells_matching_tra_v_gene <- names(which(lapply(which_db_tra_v_gene_in_sname_1.10x_v_gene,length) > 0))
cells_matching_tra_v_gene.seq <- list()
for (cell in cells_matching_tra_v_gene){
  idx <- which_db_tra_v_gene_in_sname_1.10x_v_gene[[cell]]
  cells_matching_tra_v_gene.seq[[cell]] <- unique(db$TRAV[idx])
}

cells_matching_trb_v_gene <- names(which(lapply(which_db_trb_v_gene_in_sname_1.10x_v_gene,length) > 0))
cells_matching_trb_v_gene.seq <- list()
for (cell in cells_matching_trb_v_gene){
  idx <- which_db_trb_v_gene_in_sname_1.10x_v_gene[[cell]]
  cells_matching_trb_v_gene.seq[[cell]] <- unique(db$TRBV[idx])
}

cells_matching_trb_d_gene <- names(which(lapply(which_db_trb_d_gene_in_sname_1.10x_d_gene,length) > 0))
cells_matching_trb_d_gene.seq <- list()
for (cell in cells_matching_trb_d_gene){
  idx <- which_db_trb_d_gene_in_sname_1.10x_d_gene[[cell]]
  cells_matching_trb_d_gene.seq[[cell]] <- unique(db$TRBD[idx])
}

cells_matching_tra_j_gene <- names(which(lapply(which_db_tra_j_gene_in_sname_1.10x_j_gene,length) > 0))
cells_matching_tra_j_gene.seq <- list()
for (cell in cells_matching_tra_j_gene){
  idx <- which_db_tra_j_gene_in_sname_1.10x_j_gene[[cell]]
  cells_matching_tra_j_gene.seq[[cell]] <- unique(db$TRAJ[idx])
}


cells_matching_trb_j_gene <- names(which(lapply(which_db_trb_j_gene_in_sname_1.10x_j_gene,length) > 0))
cells_matching_trb_j_gene.seq <- list()
for (cell in cells_matching_trb_j_gene){
  idx <- which_db_trb_j_gene_in_sname_1.10x_j_gene[[cell]]
  cells_matching_trb_j_gene.seq[[cell]] <- unique(db$TRBJ[idx])
}

tra_cdr3 <- data.frame(cell=names(cells_matching_tra_cdr3.seq),tra_cdr3=unlist(lapply(cells_matching_tra_cdr3.seq,paste0,collapse=";")))
trb_cdr3 <- data.frame(cell=names(cells_matching_trb_cdr3.seq),trb_cdr3=unlist(lapply(cells_matching_trb_cdr3.seq,paste0,collapse=";")))
tra_v_gene <- data.frame(cell=names(cells_matching_tra_v_gene.seq),tra_v_gene=unlist(lapply(cells_matching_tra_v_gene.seq,paste0,collapse=";")))
trb_v_gene <- data.frame(cell=names(cells_matching_trb_v_gene.seq),trb_v_gene=unlist(lapply(cells_matching_trb_v_gene.seq,paste0,collapse=";")))
trb_d_gene <- data.frame(cell=names(cells_matching_trb_d_gene.seq),trb_d_gene=unlist(lapply(cells_matching_trb_d_gene.seq,paste0,collapse=";")))
tra_j_gene <- data.frame(cell=names(cells_matching_tra_j_gene.seq),tra_j_gene=unlist(lapply(cells_matching_tra_j_gene.seq,paste0,collapse=";")))
trb_j_gene <- data.frame(cell=names(cells_matching_trb_j_gene.seq),trb_j_gene=unlist(lapply(cells_matching_trb_j_gene.seq,paste0,collapse=";")))

df.merge <- data.frame(cell=unique(c(cells_matching_tra_cdr3,
                         cells_matching_trb_cdr3,
                         cells_matching_tra_v_gene,
                         cells_matching_trb_v_gene,
                         cells_matching_trb_d_gene,
                         cells_matching_tra_j_gene,
                         cells_matching_trb_j_gene)))

if (nrow(tra_cdr3)>0)
  df.merge <- merge(df.merge,tra_cdr3,all=T,by="cell")
if (nrow(trb_cdr3)>0)
  df.merge <- merge(df.merge,trb_cdr3,all=T,by="cell")
if (nrow(tra_v_gene)>0)
  df.merge <- merge(df.merge,tra_v_gene,all=T,by="cell")
if (nrow(trb_v_gene)>0)
  df.merge <- merge(df.merge,trb_v_gene,all=T,by="cell")
if (nrow(trb_d_gene)>0)
  df.merge <- merge(df.merge,trb_d_gene,all=T,by="cell")
if (nrow(tra_j_gene)>0)
  df.merge <- merge(df.merge,tra_j_gene,all=T,by="cell")
if (nrow(trb_j_gene)>0)
  df.merge <- merge(df.merge,trb_j_gene,all=T,by="cell")

rownames(df.merge) <- df.merge$cell
df.merge$cell <- NULL
missing <- setdiff(c("tra_cdr3","trb_cdr3","tra_v_gene","trb_v_gene","tra_j_gene","trb_j_gene","trb_d_gene"),colnames(df.merge))
df.merge[[missing]] <- NA
df.merge <- df.merge[,c("tra_cdr3","trb_cdr3","tra_v_gene","tra_j_gene","trb_v_gene","trb_j_gene","trb_d_gene")]
colnames(df.merge) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD")
head(df.merge)
dim(df.merge)

df.merge.candidates <- df.merge[!is.na(df.merge$CDR3.alpha) | !is.na(df.merge$CDR3.beta),]

for (i in 1:nrow(df.merge.candidates)){
  CDR3.alpha <- df.merge.candidates$CDR3.alpha[i]
  CDR3.beta <- df.merge.candidates$CDR3.beta[i]
  TRAV <- df.merge.candidates$TRAV[i]
  TRAJ <- df.merge.candidates$TRAJ[i]
  TRBV <- df.merge.candidates$TRBV[i]
  TRBJ <- df.merge.candidates$TRBJ[i]
  TRBD <- df.merge.candidates$TRBD[i]
  db[which(
    (db$CDR3.alpha==CDR3.alpha & db$TRAV==TRAV & db$TRAJ==TRAJ) | 
      (db$CDR3.beta==CDR3.beta & db$TRBV==TRBV & db$TRBJ==TRBJ)),]
  
  db[which(
    db$CDR3.alpha==CDR3.alpha & db$TRAV==TRAV & db$TRAJ==TRAJ),]
  
  db[which(db$CDR3.alpha==CDR3.alpha),]
  
  
}

db[which(db$CDR3.alpha=="CAESRGNTGKLIF"),]


dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")

"1_AAAGATGGTCCCTACT_1"



which_db_tra_cdr3_in_sname_1.10x_cdr3[names(which(unlist(lapply(which_db_tra_cdr3_in_sname_1.10x_cdr3,length)>0)))]
db[unlist(which_db_tra_cdr3_in_sname_1.10x_cdr3[names(which(unlist(lapply(which_db_tra_cdr3_in_sname_1.10x_cdr3,length)>0)))]),]

which_db_trb_cdr3_in_sname_1.10x_cdr3[names(which(unlist(lapply(which_db_trb_cdr3_in_sname_1.10x_cdr3,length)>0)))]
db[unlist(which_db_trb_cdr3_in_sname_1.10x_cdr3[names(which(unlist(lapply(which_db_trb_cdr3_in_sname_1.10x_cdr3,length)>0)))]),]
