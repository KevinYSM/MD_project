library(pheatmap)
library(stringr)

read_databases <- function(vdjdb_path = "~/proj/um_ss/Investigations/data/TCR_databases/VDJdb/SearchTable-2022-08-30 10 35 20.945.tsv",
                           mcpas_path = "~/proj/um_ss/Investigations/data/TCR_databases/McPAS-TCR/McPAS-TCR.csv",
                           tcr3d_cancer_path = "~/proj/um_ss/Investigations/data/TCR_databases/TCR3d/cancer_targeting_tcrs.csv",
                           tcr3d_virus_path = "~/proj/um_ss/Investigations/data/TCR_databases/TCR3d/virus_targeting-tcrs.csv"){

  vdjdb <- read.table(vdjdb_path, sep = "\t", header = T, fill = T)
  mcpas_tcr <- read.table(mcpas_path, sep = ",", header = T)
  tcr3d_cancer <- read.table(tcr3d_cancer_path, sep = ",", header = T)
  tcr3d_virus <- read.table(tcr3d_virus_path, sep = ",", header = T)

  vdjdb.simplified <- vdjdb[,c("Gene","CDR3","V","J","Epitope.gene","Species",
                               "Epitope.species","complex.id")]
  vdjdb.simplified$CDR3.alpha <- ""
  vdjdb.simplified$CDR3.beta <- ""
  vdjdb.simplified$CDR3.alpha[vdjdb.simplified$Gene=="TRA"] <- vdjdb.simplified$CDR3[
    vdjdb.simplified$Gene=="TRA"]
  vdjdb.simplified$CDR3.beta[vdjdb.simplified$Gene=="TRB"] <- vdjdb.simplified$CDR3[
    vdjdb.simplified$Gene=="TRB"]
  vdjdb.simplified$CDR3 <- NULL
  vdjdb.simplified$TRAV <- ""
  vdjdb.simplified$TRAJ <- ""
  vdjdb.simplified$TRBV <- ""
  vdjdb.simplified$TRBJ <- ""
  vdjdb.simplified$TRBD <- ""
  vdjdb.simplified$TRAV[vdjdb.simplified$Gene == "TRA"] <- vdjdb.simplified$V[
    vdjdb.simplified$Gene == "TRA"]
  vdjdb.simplified$TRAJ[vdjdb.simplified$Gene == "TRA"] <- vdjdb.simplified$J[
    vdjdb.simplified$Gene == "TRA"]
  vdjdb.simplified$TRBV[vdjdb.simplified$Gene == "TRB"] <- vdjdb.simplified$V[
    vdjdb.simplified$Gene == "TRB"]
  vdjdb.simplified$TRBJ[vdjdb.simplified$Gene == "TRB"] <- vdjdb.simplified$J[
    vdjdb.simplified$Gene == "TRB"]
  vdjdb.simplified$V <- NULL
  vdjdb.simplified$J <- NULL
  vdjdb.simplified$Gene <- NULL
  vdjdb.simplified <- vdjdb.simplified[,c("CDR3.alpha","CDR3.beta","TRAV","TRAJ",
                                          "TRBV","TRBJ","TRBD","Epitope.gene",
                                          "Epitope.species","Species","complex.id")]
  vdjdb.simplified <- unique(vdjdb.simplified)
  dim(vdjdb.simplified)

  max(table(vdjdb.simplified$complex.id[vdjdb.simplified$complex.id!=0]))

  vdjdb.simplified_0 <- vdjdb.simplified[vdjdb.simplified$complex.id==0,]
  head(vdjdb.simplified_0)

  vdjdb.simplified_paired <- vdjdb.simplified[vdjdb.simplified$complex.id!=0,]
  head(vdjdb.simplified_paired)

  unique_complex_id <- unique(vdjdb.simplified_paired$complex.id)
  paired <- list()
  for (id in unique_complex_id){
    tmp <- vdjdb.simplified_paired[vdjdb.simplified_paired$complex.id == id,]
    tmp_a <- tmp[tmp$CDR3.alpha!="",]
    tmp_b <- tmp[tmp$CDR3.beta!="",]
    stopifnot(nrow(tmp_a)==1 && nrow(tmp_b)==1)
    tmp_a$CDR3.beta <- tmp_b$CDR3.beta
    tmp_a$TRBV <- tmp_b$TRBV
    tmp_a$TRBJ <- tmp_b$TRBJ
    tmp_a$TRBD <- tmp_b$TRBD
    paired[[id]] <- tmp_a
  }
  paired <- do.call("rbind",paired)
  rownames(paired) <- NULL
  head(paired)
  vdjdb.simplified <- rbind(vdjdb.simplified_0,paired)
  vdjdb.simplified$complex.id <- NULL

  colnames(vdjdb.simplified) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD",
                                  "Antigen","Antigen_species_or_pathology","Species")

  mcpas_tcr.simplified <- mcpas_tcr[,c("CDR3.alpha.aa","CDR3.beta.aa","TRAV","TRAJ","TRBV","TRBD",
                                       "TRBJ","Antigen.protein","Species","Pathology")]
  mcpas_tcr.simplified <- mcpas_tcr.simplified[,c("CDR3.alpha.aa","CDR3.beta.aa",
    "TRAV","TRAJ","TRBV","TRBJ","TRBD","Antigen.protein","Pathology","Species")]
  colnames(mcpas_tcr.simplified) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD",
                                  "Antigen","Antigen_species_or_pathology","Species")

  tcr3d_cancer.simplified <- tcr3d_cancer[,c("CDR3.alpha.","CDR3.beta.","TRAV",
                                             "TRBV","Antigen","Cancer.BR.Type")]
  tcr3d_cancer.simplified$TRAJ <- ""
  tcr3d_cancer.simplified$TRBJ <- ""
  tcr3d_cancer.simplified$TRBD <- ""
  tcr3d_cancer.simplified$Species <- ""
  tcr3d_cancer.simplified <- tcr3d_cancer.simplified[,
    c("CDR3.alpha.","CDR3.beta.","TRAV","TRAJ","TRBV","TRBJ","TRBD","Antigen","Cancer.BR.Type","Species")]
  colnames(tcr3d_cancer.simplified) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD",
                                         "Antigen","Antigen_species_or_pathology","Species")

  tcr3d_virus.simplified <- tcr3d_virus[,c("CDR3.alpha.","CDR3.beta.","TRAV",
                                           "TRBV","Antigen","Virus.BR.Type")]
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

  db <- db[which(db$TRAJ!="CATSESSGQTYEQYF"),]

  db$TRBV <- gsub("^TRBV","",db$TRBV)
  idx_split <- grepl("/",db$TRBV) & !grepl("DV",db$TRBV)
  db <- split_entries(db,coln="TRBV",delim="/",idx_split=idx_split)
  db$TRBV[!is.na(db$TRBV)] <- paste0("TRBV",db$TRBV[!is.na(db$TRBV)])

  db <- split_entries(db,coln = "TRAV")
  db <- split_entries(db,coln = "TRAJ")
  db <- split_entries(db,coln = "TRBV")
  db <- split_entries(db,coln = "TRBJ")
  db <- split_entries(db,coln = "TRBD")

  unique(db$TRAV)
  unique(db$TRAJ)
  unique(db$TRBV)
  unique(db$TRBJ)
  unique(db$TRBD)

  db$TRAV <- str_split_fixed(db$TRAV,":",2)[,1]
  db$TRAV <- gsub("^TRAV","",db$TRAV)
  idx_split <- grepl("/",db$TRAV) & !grepl("DV",db$TRAV)
  db <- split_entries(db,coln="TRAV",delim="/",idx_split=idx_split)
  db$TRAV[!is.na(db$TRAV)] <- paste0("TRAV",db$TRAV[!is.na(db$TRAV)])
  unique(db$TRAV)

  db$TRAV <- gsub("/$","",db$TRAV)
  db$TRAJ <- gsub("/$","",db$TRAJ)
  db$TRBV <- gsub("/$","",db$TRBV)
  db$TRBJ <- gsub("/$","",db$TRBJ)
  db$TRBD <- gsub("/$","",db$TRBD)

  db$TRAJ <- gsub("^TRAJ","",db$TRAJ)
  idx_split <- grepl("/",db$TRAJ) & !grepl("DV",db$TRAJ)
  db$TRAJ[idx_split]
  db <- split_entries(db,coln="TRAJ",delim="/",idx_split=idx_split)
  db$TRAJ[!is.na(db$TRAJ)] <- paste0("TRAJ",db$TRAJ[!is.na(db$TRAJ)])
  unique(db$TRAJ)
  unique(db$TRBJ)
  unique(db$TRBD)

  db$TRAV <- str_split_fixed(db$TRAV,"\\*",2)[,1]
  db$TRAJ <- str_split_fixed(db$TRAJ,"\\*",2)[,1]
  db$TRBV <- str_split_fixed(db$TRBV,"\\*",2)[,1]
  db$TRBJ <- str_split_fixed(db$TRBJ,"\\*",2)[,1]
  db$TRBD <- str_split_fixed(db$TRBD,"\\*",2)[,1]

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

  db$TRBD[which(db$TRBD=="unknown")] <- NA
  db$TRBD[which(db$TRBD=="na")] <- NA
  db[db==""] <- NA

  # malformed entry in original csv

  unique(db$TRAV)
  unique(db$TRAJ)
  unique(db$TRBV)
  unique(db$TRBJ)
  unique(db$TRBD)
  db <- db[!grepl("Donor",db$TRBJ),] # Malformed original entry
  db[db==""] <- NA

  db$Antigen[db$Antigen %in% c("MelanA/MART1","MLANA","Melan-A/MART-1")] <- "MART1"
  db$Antigen[db$Antigen %in% c("M1","Matrix protein (M1)","M")] <- "M1"
  db$Antigen[db$Antigen %in% c("Gag polyprotein RQ13","Gag")] <- "Gag"
  db$Antigen[db$Antigen %in% c("BZLF-1","BZLF1")] <- "BZLF1"
  db$Antigen[db$Antigen %in% c("DQ2.5-glia-?2")] <- "DQ2.5-glia-alpha2"
  db$Antigen[db$Antigen %in% c("DQ8.5-GLIA-GAMMA1")] <- "DQ8.5-glia-gamma1"
  db$Antigen[db$Antigen %in% c("TT")] <- "Tetanus toxin"
  db$Antigen[db$Antigen %in% c("EphA2")] <- "EPHA2"
  db$Antigen[db$Antigen %in% c("gp100")] <- "PMEL"
  db$Antigen[db$Antigen %in% c("p53")] <- "TP53"
  db$Antigen[db$Antigen %in% c("BMLF-1")] <- "BMLF1"
  db$Antigen[db$Antigen %in% c("EBNA-6")] <- "EBNA6"

  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology=="CMV"] <- "Cytomegalovirus (CMV)"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology=="InfluenzaA"] <- "Influenza"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology=="Melanoma"] <- "HomoSapiens"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology=="Human"] <- "HomoSapiens"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("Celiac disease","TriticumAestivum")] <- "Celiac disease"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("Yellow fever virus","YFV")] <- "Yellow fever virus"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("Epstein Barr virus (EBV)" ,"EBV")] <- "Epstein Barr virus (EBV)"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("Human immunodeficiency virus (HIV)","HIV-1")] <- "Human immunodeficiency virus (HIV)"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("HomoSapiens","Tumor associated antigen (TAA)")] <- "HomoSapiens"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("Alzheimer's disease") & db$Antigen=="BZLF1"] <- "Epstein Barr virus (EBV)"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("HCV")] <- "Hepatitis C virus (HCV)"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("HPV")] <- "Human papillomavirus (HPV)"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("HTLV-1")] <- "Human T-lymphotropic virus 1 (HTLV-1)"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("StreptomycesKanamyceticus")] <- "Streptomyces kanamyceticus"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("M.tuberculosis")] <- "Mycobacterium tuberculosis"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("Hepatitis E virus infection (cHEV)")] <- "Hepatitis E virus (HEV)"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("synthetic") & db$Antigen=="M1-G4E"] <- "Influenza"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("M.Tuberculosis","M. tuberculosis")] <- "Mycobacterium tuberculosis"
  db$Antigen_species_or_pathology[db$Antigen_species_or_pathology %in% c("COVID-19")] <- "SARS-CoV-2"
  db$Species[db$Species=="HomoSapiens"] <- "Human"

  db <- db[!(is.na(db$CDR3.alpha) & !is.na(db$CDR3.beta)),]

  db <- unique(db)
  rownames(db) <- NULL
  head(db)

  #table(rowSums(table(db$CDR3.alpha,db$database) > 0))

  return(db)
}

get_db_matches <- function(dat.tumor_til.vdj,db,sname){
  idx_1 <- unlist(lapply(dat.tumor_til.vdj[[sname]]$vdj.10x_cdr3, function(x){
    all(x!="None")
  }))

  tmp_1.10x_cdr3 <- dat.tumor_til.vdj[[sname]]$vdj.10x_cdr3[which(idx_1)]
  tmp_1.10x_v_gene <- dat.tumor_til.vdj[[sname]]$vdj.10x_v_gene[which(idx_1)]
  tmp_1.10x_d_gene <- dat.tumor_til.vdj[[sname]]$vdj.10x_d_gene[which(idx_1)]
  tmp_1.10x_j_gene <- dat.tumor_til.vdj[[sname]]$vdj.10x_j_gene[which(idx_1)]
  tmp_1.10x_c_gene <- dat.tumor_til.vdj[[sname]]$vdj.10x_c_gene[which(idx_1)]

  which_db_tra_cdr3_in_sname.10x_cdr3 <- lapply(tmp_1.10x_cdr3,function(x){
    which(db$CDR3.alpha %in% x)
  })

  which_db_trb_cdr3_in_sname.10x_cdr3 <- lapply(tmp_1.10x_cdr3,function(x){
    which(db$CDR3.beta %in% x)
  })

  which_db_tra_v_gene_in_sname.10x_v_gene <- lapply(tmp_1.10x_v_gene,function(x){
    which(db$TRAV %in% x)
  })

  which_db_trb_v_gene_in_sname.10x_v_gene <- lapply(tmp_1.10x_v_gene,function(x){
    which(db$TRBV %in% x)
  })

  which_db_trb_d_gene_in_sname.10x_d_gene <- lapply(tmp_1.10x_d_gene,function(x){
    which(db$TRBD %in% x)
  })

  which_db_tra_j_gene_in_sname.10x_j_gene <- lapply(tmp_1.10x_j_gene,function(x){
    which(db$TRAJ %in% x)
  })

  which_db_trb_j_gene_in_sname.10x_j_gene <- lapply(tmp_1.10x_j_gene,function(x){
    which(db$TRBJ %in% x)
  })

  cells_matching_tra_cdr3 <- names(which(lapply(which_db_tra_cdr3_in_sname.10x_cdr3,length) > 0))
  cells_matching_tra_cdr3.seq <- list()
  for (cell in cells_matching_tra_cdr3){
    idx <- which_db_tra_cdr3_in_sname.10x_cdr3[[cell]]
    cells_matching_tra_cdr3.seq[[cell]] <- unique(db$CDR3.alpha[idx])
  }

  cells_matching_trb_cdr3 <- names(which(lapply(which_db_trb_cdr3_in_sname.10x_cdr3,length) > 0))
  cells_matching_trb_cdr3.seq <- list()
  for (cell in cells_matching_trb_cdr3){
    idx <- which_db_trb_cdr3_in_sname.10x_cdr3[[cell]]
    cells_matching_trb_cdr3.seq[[cell]] <- unique(db$CDR3.beta[idx])
  }

  cells_matching_tra_v_gene <- names(which(lapply(which_db_tra_v_gene_in_sname.10x_v_gene,length) > 0))
  cells_matching_tra_v_gene.seq <- list()
  for (cell in cells_matching_tra_v_gene){
    idx <- which_db_tra_v_gene_in_sname.10x_v_gene[[cell]]
    cells_matching_tra_v_gene.seq[[cell]] <- unique(db$TRAV[idx])
  }

  cells_matching_trb_v_gene <- names(which(lapply(which_db_trb_v_gene_in_sname.10x_v_gene,length) > 0))
  cells_matching_trb_v_gene.seq <- list()
  for (cell in cells_matching_trb_v_gene){
    idx <- which_db_trb_v_gene_in_sname.10x_v_gene[[cell]]
    cells_matching_trb_v_gene.seq[[cell]] <- unique(db$TRBV[idx])
  }

  cells_matching_trb_d_gene <- names(which(lapply(which_db_trb_d_gene_in_sname.10x_d_gene,length) > 0))
  cells_matching_trb_d_gene.seq <- list()
  for (cell in cells_matching_trb_d_gene){
    idx <- which_db_trb_d_gene_in_sname.10x_d_gene[[cell]]
    cells_matching_trb_d_gene.seq[[cell]] <- unique(db$TRBD[idx])
  }

  cells_matching_tra_j_gene <- names(which(lapply(which_db_tra_j_gene_in_sname.10x_j_gene,length) > 0))
  cells_matching_tra_j_gene.seq <- list()
  for (cell in cells_matching_tra_j_gene){
    idx <- which_db_tra_j_gene_in_sname.10x_j_gene[[cell]]
    cells_matching_tra_j_gene.seq[[cell]] <- unique(db$TRAJ[idx])
  }

  cells_matching_trb_j_gene <- names(which(lapply(which_db_trb_j_gene_in_sname.10x_j_gene,length) > 0))
  cells_matching_trb_j_gene.seq <- list()
  for (cell in cells_matching_trb_j_gene){
    idx <- which_db_trb_j_gene_in_sname.10x_j_gene[[cell]]
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
  for (m in missing){
    df.merge[[m]] <- NA
  }
  df.merge <- df.merge[,c("tra_cdr3","trb_cdr3","tra_v_gene","tra_j_gene","trb_v_gene","trb_j_gene","trb_d_gene")]
  colnames(df.merge) <- c("CDR3.alpha","CDR3.beta","TRAV","TRAJ","TRBV","TRBJ","TRBD")
  head(df.merge)
  dim(df.merge)

  df.merge.candidates <- df.merge[!is.na(df.merge$CDR3.alpha) | !is.na(df.merge$CDR3.beta),]

  df.merge.candidates.db_match <- list()
  for (i in 1:nrow(df.merge.candidates)){
    CDR3.alpha <- df.merge.candidates$CDR3.alpha[i]
    CDR3.beta <- df.merge.candidates$CDR3.beta[i]
    TRAV <- df.merge.candidates$TRAV[i]
    TRAJ <- df.merge.candidates$TRAJ[i]
    TRBV <- df.merge.candidates$TRBV[i]
    TRBJ <- df.merge.candidates$TRBJ[i]
    TRBD <- df.merge.candidates$TRBD[i]
     # df.merge.candidates.db_match[[i]] <- db[which(
     #   (db$CDR3.alpha==CDR3.alpha & db$TRAV==TRAV & db$TRAJ==TRAJ) |
     #     (db$CDR3.beta==CDR3.beta & db$TRBV==TRBV & db$TRBJ==TRBJ)),]
    df.merge.candidates.db_match[[i]] <- db[which(
      db$CDR3.alpha==CDR3.alpha | db$CDR3.beta==CDR3.beta &
        !(is.na(db$CDR3.alpha==CDR3.alpha) & is.na(db$CDR3.beta==CDR3.beta))),]
    #df.merge.candidates.db_match[[i]] <- db[which(
    #  db$CDR3.alpha==CDR3.alpha & db$CDR3.beta==CDR3.beta),]
    #df.merge.candidates.db_match[[i]] <- db[which(db$CDR3.beta==CDR3.beta),]
  }
  names(df.merge.candidates.db_match) <- rownames(df.merge.candidates)

  df.merge.candidates.db_match <- df.merge.candidates.db_match[which(unlist(lapply(df.merge.candidates.db_match,nrow)) > 0)]
  for (cell in names(df.merge.candidates.db_match)){
    df.merge.candidates.db_match[[cell]]$cell <- cell
  }
  df.merge.candidates.db_match <- do.call("rbind",df.merge.candidates.db_match)
  rownames(df.merge.candidates.db_match) <- NULL

  return(df.merge.candidates.db_match)
}

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/"
dat.tumor_til.vdj <- readRDS(file=paste0(outdir,"dat.tumor_til.vdj.rda"))

db <- read_databases()

db <- db[!(is.na(db$CDR3.alpha) & is.na(db$CDR3.beta)),]
db <- db[db$Species!="Mouse",]

dat.tumor <- readRDS(file=paste0(outdir,"seu.filt.rda"))

trav <- grep(pattern="TRAV",rownames(dat.tumor@assays$RNA@counts),value = T)
traj <- grep(pattern="TRAJ",rownames(dat.tumor@assays$RNA@counts),value = T)
trac <- grep(pattern="TRAC",rownames(dat.tumor@assays$RNA@counts),value = T)
trbv <- grep(pattern="TRBV",rownames(dat.tumor@assays$RNA@counts),value = T)
trbj <- grep(pattern="TRBJ",rownames(dat.tumor@assays$RNA@counts),value = T)
trbd <- grep(pattern="TRBD",rownames(dat.tumor@assays$RNA@counts),value = T)
trbc <- grep(pattern="TRBC",rownames(dat.tumor@assays$RNA@counts),value = T)

db$TRAV[db$TRAV=="TRAV14/DV4"] <- "TRAV14DV4"
db$TRAV[db$TRAV=="TRAV23/DV6"] <- "TRAV23DV6"
db$TRAV[db$TRAV=="TRAV38-2/DV8"] <- "TRAV38-2DV8"
db$TRAV[db$TRAV=="TRAV29/DV5"] <- "TRAV29DV5"
db$TRAV[db$TRAV=="TRAV36/DV7"] <- "TRAV36DV7"

sort(setdiff(db$TRAV,trav))
sort(trav)

db[db$TRAV=="TRAV6-01",]

db.matches <- list()
for (sname in names(dat.tumor_til.vdj)){
  print(sname)
  db.matches[[sname]] <- get_db_matches(dat.tumor_til.vdj,db,sname)
}

for (sname in names(db.matches)){
  db.matches[[sname]]$Sample.ID <- sname
  db.matches[[sname]]$UM.ID <- str_split_fixed(sname,"_",2)[,1]
  db.matches[[sname]]$project.name <- str_split_fixed(sname,"_",2)[,2]
  db.matches[[sname]]$project.name[db.matches[[sname]]$project.name=="til"] <- "TIL"
  db.matches[[sname]]$project.name[db.matches[[sname]]$project.name=="biopsy"] <- "Biopsy"
}

db.matches <- do.call("rbind",db.matches)
rownames(db.matches) <- NULL
head(db.matches)

table(db.matches$Species)

tmp <- db.matches[db.matches$Species!="Mouse",c("Antigen","UM.ID","project.name","cell")]

dim(tmp)
dim(unique(tmp))
tmp <- unique(tmp)
m <- table(tmp$Antigen,paste0(tmp$UM.ID,":",tmp$project.name))

annot_row <- data.frame(Antigen=db.matches$Antigen,origin=db.matches$Antigen_species_or_pathology,Species=db.matches$Species)
annot_row$origin[annot_row$origin=="HomoSapiens"] <- "Homo sapiens"
annot_row <- annot_row[annot_row$Species!="Mouse",]
annot_row$Species <- NULL
annot_row <- unique(annot_row)
annot_row <- annot_row[!is.na(annot_row$Antigen),]
annot_row[annot_row$Antigen %in% names(which(table(annot_row$Antigen)>1)),]
annot_row <- unique(annot_row)
annot_row[annot_row$Antigen %in% names(which(table(annot_row$Antigen)>1)),]
annot_row <- unique(annot_row)
annot_row[annot_row$Antigen=="NY-ESO-1" & annot_row$origin=="Tumor",]$origin <- "Homo sapiens"
#annot_row[annot_row$Antigen=="NS3" & annot_row$origin %in% c("Hepatitis C virus (HCV)","DENV1","DENV3/4"),]$origin <- "Hepatitis C virus (HCV) / DENV1/3/4"
annot_row[annot_row$Antigen=="NS3" & annot_row$origin %in% c("Hepatitis C virus (HCV)","DENV1","DENV3/4"),]$origin <- "Hepatitis C virus (HCV)"
annot_row <- unique(annot_row)
rownames(annot_row) <- annot_row$Antigen
annot_row$Antigen <- NULL
annot_row$antigen <- rownames(annot_row)
annot_row <- annot_row[order(annot_row$antigen,decreasing = T),]
annot_row <- annot_row[order(annot_row$origin,decreasing = T),]
annot_row$antigen <- NULL
annot_row$origin[annot_row$origin=="SARS-CoV-2"] <- "Coronaviruses / SARS-CoV-2"
annot_row

rn <- rownames(annot_row)[order(annot_row$origin)]
m <- m[rn,]
m <- m[,order(as.numeric(str_split_fixed(gsub("UM","",colnames(m)),":",2)[,1]))]

pdf(file=paste0(outdir,"tcr_database_matches.pdf"),width = 10.5,height = 10)
pheatmap(m,border_color = NA,na_col = "white",colorRampPalette(c("white", "red"))(100),
         annotation_row = annot_row,display_numbers = m,cluster_rows = F,
         cluster_cols = F)
dev.off()

m2 <- list()
for (orig in unique(annot_row$origin)){
  idx <- rownames(m) %in% rownames(annot_row)[annot_row$origin==orig]
  if (sum(idx)==1){
    m2[[orig]] <- m[idx,]
  } else {
    m2[[orig]] <- colSums(m[idx,])
  }
}
m2 <- do.call("rbind",m2)
m2 <- m2[order(-rowSums(m2)),]

pdf(file=paste0(outdir,"tcr_database_matches_per_origin.pdf"),width = 8.5,height = 4)
pheatmap(m2,border_color = NA,na_col = "white",colorRampPalette(c("white", "red"))(100),
         display_numbers = m2,cluster_rows = F,
         cluster_cols = F)
dev.off()

m3 <- m[rownames(m) %in% rownames(annot_row)[annot_row$origin=="Homo sapiens"],]
m3 <- m3[order(-rowSums(m3)),]

pdf(file=paste0(outdir,"tcr_database_matches_homo_sapiens_only.pdf"),width = 4.5,height = 4)
pheatmap(m3,border_color = NA,na_col = "white",colorRampPalette(c("white", "red"))(100),
         display_numbers = m3, cluster_rows = F,
         cluster_cols = F)
dev.off()

# Sanity checks -------------
# dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
#

db.matches[which(db.matches$Antigen=="KRAS"),]

dat.cca@meta.data["2_CCATTCGAGTGGTAGC_3",c("cdr3","v_gene","d_gene","j_gene","c_gene")]

# Compare Vasu clonotypes against TCR databases --------------------------------
library(readxl)
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
  vdj.takara[,c("chains","cdr3","v_gene","d_gene","j_gene","c_gene")],1,paste0,
  collapse = ":")

vdj.takara <- vdj.takara[vdj.takara$UM.ID!="",]

head(vdj.takara)

vdj.takara <- split(vdj.takara,vdj.takara$UM.ID)

vdj.takara <- lapply(vdj.takara,function(x){
  rownames(x) <- paste0(1:nrow(x),"_",unique(x$UM.ID))
  return(x)
})

get_db_matches.vasu <- function(vdj.takara,db,sname){
  vdj <- vdj.takara[[sname]]
  db.matches <- list()
  for (i in 1:nrow(vdj)){
    if (vdj[i,]$chains=="TRA"){
      tra <- vdj[i,]
      db.matches[[rownames(vdj)[i]]] <- db[which(db$CDR3.alpha == tra$cdr3),]
    } else if (vdj[i,]$chains=="TRB"){
      trb <- vdj[i,]
      db.matches[[rownames(vdj)[i]]] <- db[which(db$CDR3.beta == trb$cdr3),]
    } else {
      stop("Unexpected chain")
    }
  }

  db.matches <- do.call("rbind",db.matches)
  db.matches$id <- str_split_fixed(rownames(db.matches),"\\.",2)[,1]
  rownames(db.matches) <- NULL
  return(db.matches)
}

db_matches.vasu <- list()
db_matches.vasu <- lapply(names(vdj.takara),function(sname){
  get_db_matches.vasu(vdj.takara,db,sname)
})
names(db_matches.vasu) <- names(vdj.takara)

for (sname in names(db_matches.vasu)){
  out <- list()
  for (i in 1:nrow(db_matches.vasu[[sname]])){
    id <- db_matches.vasu[[sname]][i,]$id
    takara <- vdj.takara[[sname]][id,]
    db_matches <- db_matches.vasu[[sname]][i,]
    colnames(takara) <- paste0("takara.",colnames(takara))
    colnames(db_matches) <- paste0("db.",colnames(db_matches))
    out[[i]] <- cbind(db_matches,takara)
  }
  db_matches.vasu[[sname]] <- do.call("rbind",out)
}

#lapply(db_matches.vasu,function(x){table(x$db.Antigen)})
#lapply(db_matches.vasu,function(x){table(x$db.Antigen_species_or_pathology)})
#lapply(db_matches.vasu,function(x){table(x$db.Species)})
#lapply(db_matches.vasu,function(x){x[x$db.Antigen_species_or_pathology=="HomoSapiens",]})

unique_samples <- names(db_matches.vasu)
overlaps_cdr3 <- list()
for (sname_1 in unique_samples){
  for (sname_2 in unique_samples){
    if (sname_1!=sname_2){
      overlaps_cdr3[[paste0(sname_1,":",sname_2)]] <- intersect(db_matches.vasu[[sname_1]]$takara.cdr3,
                                                              db_matches.vasu[[sname_2]]$takara.cdr3)
    }
  }
}

nms <- sapply(names(overlaps_cdr3),function(pair){
  paste0(sort(str_split_fixed(pair,":",2)),collapse=":")
})

stopifnot(all(unique(nms) %in% names(overlaps_cdr3)))
overlaps_cdr3 <- overlaps_cdr3[unique(nms)]

for (sname in names(db_matches.vasu)){
  db_matches.vasu[[sname]]$takara.shared_cdr3 <- ""
  for (cdr3 in db_matches.vasu[[sname]]$takara.cdr3){
    matching_pairs <- which(unlist(lapply(overlaps_cdr3, function(x){cdr3 %in% x})))
    if (length(matching_pairs)>0){
      matching_pairs <- paste0(names(matching_pairs),collapse=";")
      db_matches.vasu[[sname]]$takara.shared_cdr3[db_matches.vasu[[sname]]$takara.cdr3 == cdr3] <- matching_pairs
    }
  }
}

for (sname in names(db_matches.vasu)){
  db_matches.vasu[[sname]]$db.CDR3.alpha.match <- db_matches.vasu[[sname]]$db.CDR3.alpha==db_matches.vasu[[sname]]$takara.cdr3
  db_matches.vasu[[sname]]$db.CDR3.beta.match <- db_matches.vasu[[sname]]$db.CDR3.beta==db_matches.vasu[[sname]]$takara.cdr3
  db_matches.vasu[[sname]]$db.TRAV.match <- db_matches.vasu[[sname]]$db.TRAV==db_matches.vasu[[sname]]$takara.v_gene
  db_matches.vasu[[sname]]$db.TRAJ.match <- db_matches.vasu[[sname]]$db.TRAJ==db_matches.vasu[[sname]]$takara.j_gene
  db_matches.vasu[[sname]]$db.TRBV.match <- db_matches.vasu[[sname]]$db.TRBV==db_matches.vasu[[sname]]$takara.v_gene
  db_matches.vasu[[sname]]$db.TRBJ.match <- db_matches.vasu[[sname]]$db.TRBJ==db_matches.vasu[[sname]]$takara.j_gene

  db_matches.vasu[[sname]]$db.CDR3.alpha_and_CDR3.alpha_in_takara <-
    db_matches.vasu[[sname]]$db.CDR3.alpha %in% db_matches.vasu[[sname]]$takara.cdr3 &
    db_matches.vasu[[sname]]$db.CDR3.beta %in% db_matches.vasu[[sname]]$takara.cdr3
}

library(WriteXLS)
WriteXLS(db_matches.vasu,ExcelFileName = paste0(outdir,"db_matches.vasu.xlsx"),
         AdjWidth = T,AutoFilter = T,BoldHeaderRow = T,row.names = F)
