dat.cca <- readRDS("~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")

# Compare with TCR data from Takara --------------------------------------------
library(readxl)
library(stringr)
takara <- list()
fnames <- Sys.glob("~/proj/um_ss/Pipelines/takara_smarter/preprocessing/results_both/report/results_both_*_mig_cdr3_report.xlsx")
for (fname in fnames){
  sname <- str_split_fixed(basename(fname),"_",4)[,3]
  tmp_A <- as.data.frame(read_excel(fname,sheet = paste0(sname,"_TRA_clone")),stringsAsFactors=F)
  tmp_B <- as.data.frame(read_excel(fname,sheet = paste0(sname,"_TRB_clone")),stringsAsFactors=F)
  takara[[basename(fname)]] <- unique(rbind(tmp_A[,c("CDR3 Amino Acid Sequence","V segment","D segment","J segment")],
                                            tmp_B[,c("CDR3 Amino Acid Sequence","V segment","D segment","J segment")]))
}

takara <- do.call("rbind",takara)
rownames(takara) <- NULL

dat.cca@meta.data[,c("cdr3","v_gene","d_gene","j_gene")]

#intersect(dat.cca@meta.data$cdr3,takara$`CDR3 Amino Acid Sequence`)

tcr_ss <- data.frame(cell_id=rownames(dat.cca@meta.data),cdr3=dat.cca@meta.data$cdr3)
tcr_ss <- tcr_ss[!is.na(tcr_ss$cdr3),]

library(tidyr)
library(dplyr)

tcr_ss <- tcr_ss |> mutate(cdr3 = strsplit(as.character(cdr3), ";")) |> unnest(cdr3)
tcr_ss <- as.data.frame(tcr_ss)
tcr_ss <- unique(tcr_ss)

tcr_ss$in_takara <- tcr_ss$cdr3 %in% takara$`CDR3 Amino Acid Sequence`
table(tcr_ss$in_takara)
length(unique(tcr_ss$cell_id[tcr_ss$in_takara]))

in_takara <- unique(tcr_ss$cell_id[tcr_ss$in_takara])
in_takara <- ifelse(rownames(dat.cca@meta.data) %in% in_takara,"Yes","No")
names(in_takara) <- rownames(dat.cca@meta.data)
dat.cca <- AddMetaData(dat.cca, metadata = in_takara, col.name = "TCR_found")
dat.cca@meta.data$TCR_found



pdf(file = paste0(outdir,"dat.",nm,".TCR_found.pdf"),width = 12,height = 10)
DimPlot(dat.cca, label = T, reduction = "tsne",group.by = "TCR_found")
dev.off()

pdf(file = paste0(outdir,"dat.",nm,".TCR_found.project.pdf"),width = 33,height = 10)
DimPlot(dat.cca, label = T, reduction = "tsne",group.by = "TCR_found",split.by = "project.name")
dev.off()

dat.cca <- AddMetaData(dat.cca,
                       metadata = setNames(paste0(dat.cca@meta.data$Sample.ID,":",dat.cca@meta.data$UM.ID),rownames(dat.cca@meta.data)),
                       col.name = "Sample.ID.UM.ID")

pdf(file = paste0(outdir,"dat.",nm,".TCR_found.SampleID.pdf"),width = 220,height = 10)
DimPlot(dat.cca, label = T, reduction = "tsne",group.by = "TCR_found",split.by = "Sample.ID.UM.ID")
dev.off()

x <- as.data.frame(table(dat.cca@meta.data$TCR_found,dat.cca@meta.data$Sample.ID.UM.ID))
colnames(x) <- c("TCR_found","Sample","Frequency")
x <- x[x$TCR_found=="Yes",]

x <- x[order(x$Frequency,decreasing = T),]

x$Most_frequenct_CDR3 <- ""
for (nm in unique(dat.cca@meta.data$Sample.ID.UM.ID)){
  cdr3 <- sort(table(dat.cca@meta.data$cdr3[dat.cca@meta.data$Sample.ID.UM.ID == nm  & dat.cca@meta.data$TCR_found == "Yes"]),decreasing = T)[1]
  if (!is.null(cdr3)){
    cdr3 <- paste0(names(cdr3),": ",as.numeric(cdr3))
    x$Most_frequenct_CDR3[x$Sample == nm] <- cdr3
  }
}

