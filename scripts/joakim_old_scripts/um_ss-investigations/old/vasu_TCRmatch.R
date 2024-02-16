library(readxl)
library(WriteXLS)
library(stringr)

outdir <- "~/proj/um_ss/Investigations/seurat/results/v16/"

# Write input to TCRmatch ------------------------------------------------------
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

vdj.takara <- read_vdj.takara()
vdj.takara.trb <- vdj.takara[vdj.takara$chains=="TRB",]
vdj.takara.trb$tcrmatch_input <- ""
vdj.takara.trb$tcrmatch_input <- vdj.takara.trb$cdr3
vdj.takara.trb$tcrmatch_input <- gsub(pattern="^C",replacement="",vdj.takara.trb$tcrmatch_input)
vdj.takara.trb$tcrmatch_input <- gsub(pattern="F$",replacement="",vdj.takara.trb$tcrmatch_input)
write.table(unique(vdj.takara.trb$tcrmatch_input), file = paste0(outdir,"vdj.takara.tcrmatch_input.txt"),
            sep = "\t", col.names = F, quote = F, row.names = F)

# Read TCRmatch output ---------------------------------------------------------
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

pth <- "~/nimbus/data/proj/um_ss/Investigations/seurat/results/v16/vdj.takara.tcrmatch_output_2.txt"

tcrmatch_output <- read_tcrmatch(pth)

vdj.takara.trb.merged <- merge(vdj.takara.trb, tcrmatch_output, 
                               by.x="tcrmatch_input", 
                               by.y="input_sequence", all.x=T)
vdj.takara.trb.merged$epitopes[is.na(vdj.takara.trb.merged$epitopes)] <- ""
vdj.takara.trb.merged$antigen[is.na(vdj.takara.trb.merged$antigen)] <- ""
vdj.takara.trb.merged$source_organism[is.na(vdj.takara.trb.merged$source_organism)] <- ""

for (i in 1:nrow(vdj.takara.trb.merged)){
  vdj.takara.trb.merged$antigen[i] <- ifelse(
    is.na(vdj.takara.trb.merged$antigen[i]),
    NA,
    paste0(unique(unlist(strsplit(vdj.takara.trb.merged$antigen[i],","))),collapse = ";"))
  vdj.takara.trb.merged$source_organism[i] <- ifelse(
    is.na(vdj.takara.trb.merged$source_organism[i]),
    NA,
    paste0(unique(unlist(strsplit(vdj.takara.trb.merged$source_organism[i],","))),collapse = ";"))
  vdj.takara.trb.merged$epitopes[i] <- ifelse(
    is.na(vdj.takara.trb.merged$epitopes[i]),
    NA,
    paste0(unique(unlist(strsplit(vdj.takara.trb.merged$epitopes[i],","))),collapse = ";"))
}

library(pheatmap)
tab <- table(vdj.takara.trb.merged$vasu_phenotype,vdj.takara.trb.merged$antigen)
pheatmap(tab, border_color = NA)

tab <- table(vdj.takara.trb.merged$vasu_phenotype,vdj.takara.trb.merged$source_organism)
pheatmap(tab, border_color = NA)
