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

read_iedb <- function(pth){
  iedb <- data.table::fread(pth,
                            sep = "\t", header = F, data.table = F, drop = 7)
  colnames(iedb) <- c("trimmed_seq", "original_seq", "receptor_group", "epitopes",
                      "source_organisms", "source_antigens")
  #iedb$trimmed_seq <- NULL
  #iedb$original_seq <- NULL
  #iedb$epitopes <- NULL
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

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/"

iedb <- read_iedb(pth = "~/nimbus/data/proj/um_ss/Pipelines/TCRMatch/data/IEDB_data.tsv")
iedb <- replace_internal_commas(iedb, pth_iedb_external = "~/proj/um_ss/Investigations/data/TCR_databases/IEDB/tcell_full_v3.csv")

db <- read_databases()
db <- db[db$Species!="Mouse",]
stopifnot(length(grep(";",db$Antigen))==0)
stopifnot(length(grep(",",db$source_organisms))==0)
db.iedb_format <- db[,c("CDR3.beta","Antigen","Antigen_species_or_pathology")]
db.iedb_format <- db.iedb_format[!is.na(db.iedb_format$CDR3.beta),]
db.iedb_format <- unique(db.iedb_format)

idx_nonconforming_CDR3b <- substr(db.iedb_format$CDR3.beta,1,1)!="C" | 
  (substr(db.iedb_format$CDR3.beta,nchar(db.iedb_format$CDR3.beta),nchar(db.iedb_format$CDR3.beta))!="F" & 
     substr(db.iedb_format$CDR3.beta,nchar(db.iedb_format$CDR3.beta),nchar(db.iedb_format$CDR3.beta))!="W")
db.iedb_format <- db.iedb_format[!idx_nonconforming_CDR3b,]

db.iedb_format$trimmed_seq <- gsub(pattern="^C",replacement="",db.iedb_format$CDR3.beta)
db.iedb_format$trimmed_seq <- gsub(pattern="[FW]$",replacement="",db.iedb_format$trimmed_seq)
db.iedb_format$original_seq <- db.iedb_format$trimmed_seq
db.iedb_format$CDR3.beta <- NULL
colnames(db.iedb_format)[colnames(db.iedb_format)=="Antigen"] <- "source_antigens"
colnames(db.iedb_format)[colnames(db.iedb_format)=="Antigen_species_or_pathology"] <- "source_organisms"

db.iedb_format <- db.iedb_format[!is.na(db.iedb_format$source_antigens),]
db.iedb_format <- db.iedb_format[db.iedb_format$source_antigens!="synthetic",]
stopifnot(!any(is.na(db.iedb_format$source_organisms)))

#db.iedb_format$receptor_group <- (max(iedb$receptor_group)+1):(max(iedb$receptor_group)+nrow(db.iedb_format))
#stopifnot(length(intersect(db.iedb_format$receptor_group,iedb$receptor_group))==0)

# https://www.frontiersin.org/articles/10.3389/fimmu.2021.640725/full
# Following trimming, the dataset consisted of 24,973 receptor groups, each defined by a unique CDR3β sequence

unique_original_seq <- unique(db.iedb_format$original_seq)
length(unique_original_seq)
length(db.iedb_format$original_seq)

receptor_groups <- list()
for (i in 1:length(unique_original_seq)){
  seq <- unique_original_seq[i]
  tmp <- unique(db.iedb_format[
    db.iedb_format$original_seq == seq,])
  source_antigens <- paste0(tmp$source_antigens,collapse = ";")
  source_organisms <- paste0(tmp$source_organisms,collapse = ";")
  receptor_groups[[i]] <- data.frame(trimmed_seq = seq, original_seq = seq, 
                                     receptor_group = (max(iedb$receptor_group)+1)+i,
                                     source_antigens = source_antigens, 
                                     source_organisms = source_organisms)
}
db.iedb_format <- do.call("rbind",receptor_groups)
db.iedb_format <- db.iedb_format[!(db.iedb_format$trimmed_seq %in% iedb$trimmed_seq),]
db.iedb_format$epitopes <- "Unknown"

head(db.iedb_format)
head(iedb)

iedb$source_organisms <- unlist(lapply(iedb$source_organisms,function(x){
  paste0(unlist(strsplit(x,",")),collapse = ";")
}))

cnames <- c("trimmed_seq","original_seq","receptor_group","epitopes",
            "source_organisms","source_antigens")

amino_acids <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

malformed_cdr3 <- lapply(db.iedb_format$trimmed_seq,function(x){
  any(! unlist(strsplit(x,"")) %in% amino_acids)
})

db.iedb_format <- db.iedb_format[! unlist(malformed_cdr3),]

iedb.expanded <- rbind(iedb[,cnames],db.iedb_format[,cnames])

write.table(iedb.expanded,
            file = "~/nimbus/data/proj/um_ss/Pipelines/TCRMatch/data/IEDB_data.expanded.tsv",
            sep = "\t", quote = F, col.names = T, row.names = F)

# tmp -------------------

# db <- read.table("~/nimbus/data/proj/um_ss/Pipelines/TCRMatch/data/IEDB_data.expanded.tsv", sep="\t",header=T)
# unknown_epitopes <- unique(db[db$source_antigens=="",]$epitopes)
# 
# db.known_antigen <- db[db$source_antigens!="",]
# 
# matches <- list()
# for (epitope in unknown_epitopes){
#   matches[[epitope]] <- db.known_antigen[db.known_antigen$epitopes == epitope, 
#                                          c("source_organisms","source_antigens")]
# }
# 
# do.call("rbind",matches) # None

# what does the tcrmatch score mean?