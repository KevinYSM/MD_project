# Vasu's sequences should be compared/visualised in all yTILs and biopsies ------------
dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
s <- readRDS(file = "~/proj/um_ss/Manuscript/s.rda")
s.t <- readRDS(file = "~/proj/um_ss/Manuscript/s.t.rda")
vasu_phenotype <- readRDS(
  "~/proj/um_ss/Investigations/seurat/results/liger_all/vasu_phenotype.rda")

dat.cca.vasu_phenotype <- vasu_phenotype[intersect(names(vasu_phenotype),
                                                   rownames(dat.cca@meta.data))]
dat.cca <- AddMetaData(dat.cca, metadata = dat.cca.vasu_phenotype, col.name = "vasu_phenotype")
table(dat.cca@meta.data$vasu_phenotype, dat.cca@meta.data$UM.ID)

dat.cca@meta.data$vasu_phenotype <- factor(dat.cca@meta.data$vasu_phenotype,
                                           levels = c("None", "Unidentified", 
                                                      "CD3+41bb+", "MART1"))

UM.ID <- setNames(dat.cca@meta.data$UM.ID, paste0(
  dat.cca@meta.data$Sample.ID, "_", rownames(dat.cca@meta.data)))
s.t <- AddMetaData(s.t, metadata = UM.ID, col.name = "UM.ID")

s.t$UM.ID <- factor(s.t$UM.ID, levels = paste0(
  "UM", sort(as.numeric(gsub("UM", "", unique(s.t@meta.data$UM.ID))))))
levels(s.t$UM.ID)

dat.tumor <- subset(dat.cca, project.name == "GWA-JN-388")
dat.til <- subset(dat.cca, project.name == "G18-023")


p1 <- DimPlot(subset(s.t,UM.ID %in% c("UM1","UM2","UM3","UM9","UM19","UM24")), 
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "green", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p2 <- DimPlot(subset(dat.til,UM.ID %in% c("UM1","UM2","UM3","UM9","UM19","UM24")), 
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "green", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p1 + p2



p1 <- DimPlot(s.t, 
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "green", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p2 <- DimPlot(dat.til,
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "green", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p1 + p2


p1 <- DimPlot(s.t,
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "green", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p2 <- DimPlot(subset(dat.til,UM.ID %in% c("UM1","UM2","UM3","UM9","UM19","UM24")), 
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "green", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p1 + p2

library(rliger)
lig <- readRDS(paste0("~/proj/um_ss/Investigations/seurat/results/tils/","tils.liger.rda"))

tils <- ligerToSeurat(lig)

tils@meta.data$cell_id <- str_split_fixed(rownames(tils@meta.data),"11june18_",2)[,2]

meta <- dat.til@meta.data
rownames(meta) <- paste0(meta$Sample.ID,"_",rownames(meta))

tils <- AddMetaData(tils,metadata = meta)
rm(meta)
gc()

DimPlot(tils, group.by = "vasu_phenotype", cols = c("gray", "gray", "blue","red"))
FeaturePlot(tils, features = c("CD8A","CD4","NCAM1"))

p1 <- DimPlot(s.t,
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "green", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p2 <- DimPlot(subset(tils,UM.ID %in% c("UM1","UM2","UM3","UM9","UM19","UM24")), 
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "green", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p1 + p2

tils.unintegrated <- readRDS(file=paste0("~/proj/um_ss/Investigations/seurat/results/tils/","tils.rda"))

tils.unintegrated <- lapply(tils.unintegrated,function(x){
  x <- AddMetaData(x, metadata = dat.cca.vasu_phenotype, col.name = "vasu_phenotype")  
})

names(tils.unintegrated) <- unlist(lapply(tils.unintegrated,function(x){
  unique(x@meta.data$UM.ID)  
}))

tils.unintegrated <- lapply(tils.unintegrated,function(x){
  x@meta.data$vasu_phenotype <- factor(x@meta.data$vasu_phenotype,
                                       levels=c("None","Unidentified","CD3+41bb+","MART1"))
  return(x)
})

tils.unintegrated <- lapply(tils.unintegrated,function(x){
  if (any(Idents(x) == "Gamma-delta T cells (TRGV9, TRDV2)")){
    x <- RenameIdents(x,`Gamma-delta T cells (TRGV9, TRDV2)` = "CD8 T cells")
  }
  x <- subset(x, cells = names(x@active.ident[x@active.ident != "CD8 T / NK cells"]))
})


p1 <- DimPlot(tils.unintegrated[["UM1"]],
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "gray", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p2 <- DimPlot(tils.unintegrated[["UM2"]],
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "red"),order = T) + theme_void() + ggtitle(NULL)
p3 <- DimPlot(tils.unintegrated[["UM3"]],
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "gray", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p4 <- DimPlot(tils.unintegrated[["UM9"]],
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "gray", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p5 <- DimPlot(tils.unintegrated[["UM19"]],
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray", "gray", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p_um22 <- DimPlot(tils.unintegrated[["UM22"]],
                  split.by = "UM.ID", group.by = "vasu_phenotype",
                  cols = c("gray", "gray", "blue","red"),order = T) + theme_void() + ggtitle(NULL)
p6 <- DimPlot(tils.unintegrated[["UM24"]],
              split.by = "UM.ID", group.by = "vasu_phenotype",
              cols = c("gray","red"),order = T) + theme_void() + ggtitle(NULL)
p_um10 <- DimPlot(tils.unintegrated[["UM10"]],
                  split.by = "UM.ID", group.by = "vasu_phenotype",
                  cols = c("gray", "gray", "blue","red"),order = T) + theme_void() + ggtitle(NULL)




library(gridExtra)
p.til <- grid.arrange(grobs = c(list(p1), 
                                list(p2),
                                list(p3),
                                list(p4),
                                list(p5),
                                list(p_um22),
                                list(p6),
                                list(p_um10)), ncol = 8, as.table = FALSE)


p.biopsy <- DimPlot(s.t,
                    split.by = "UM.ID", group.by = "vasu_phenotype",
                    cols = c("gray", "gray", "blue","red"),order = T) + theme_void() + ggtitle(NULL)

p.biopsy + p.til

p1 <- DimPlot(tils.unintegrated[["UM1"]]) + theme_void() + ggtitle(NULL)
p2 <- DimPlot(tils.unintegrated[["UM2"]]) + theme_void() + ggtitle(NULL)
p3 <- DimPlot(tils.unintegrated[["UM3"]]) + theme_void() + ggtitle(NULL)
p4 <- DimPlot(tils.unintegrated[["UM9"]]) + theme_void() + ggtitle(NULL)
p5 <- DimPlot(tils.unintegrated[["UM19"]]) + theme_void() + ggtitle(NULL)
p6 <- DimPlot(tils.unintegrated[["UM24"]]) + theme_void() + ggtitle(NULL)

p.til.celltype <- grid.arrange(grobs = c(list(p1), 
                                         list(p2),
                                         list(p3),
                                         list(p4),
                                         list(p5),
                                         list(p6)), ncol = 6, as.table = FALSE)

p.biopsy + p.til + p.til.celltype

pdf(file = "~/proj/um_ss/Manuscript/biopsy_tils_vasu.pdf",width = 30,height = 7)
p.biopsy + p.til
dev.off()

# yTIL TCRs should be mapped to biopsies (visualised) --------------------------
# Question: do all yTIL TCRs end up in the exhausted biopsy clusters?

# 1: Determing overlapping TCRs between samples
#dat.tumor@meta.data$cdr3
#dat.til@meta.data$cdr3

extract_tcr <- function(dat.cca){
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
  
  return(list(
    vdj.10x_cdr3 = vdj.10x_cdr3,
    vdj.10x_v_gene = vdj.10x_v_gene,
    vdj.10x_d_gene = vdj.10x_d_gene,
    vdj.10x_j_gene = vdj.10x_j_gene,
    vdj.10x_c_gene = vdj.10x_c_gene,
    vdj.10x_chains = vdj.10x_chains
  ))
}

dat.tumor.split <- SplitObject(dat.tumor, split.by = "UM.ID")
dat.til.split <- SplitObject(dat.til, split.by = "UM.ID")

dat.tumor.vdj <- lapply(dat.tumor.split,extract_tcr)
dat.til.vdj <- lapply(dat.til.split,extract_tcr)

# Matching criteria: there must be at least one complete subset of v, d, j, c and cdr3 
# between two cells in each dataset

get_matches <- function(dat.tumor.vdj,dat.til.vdj,sname){
  til_matches <- list()
  vdj.10x_cdr3.match <- list()
  vdj.10x_v_gene.match <- list()
  vdj.10x_d_gene.match <- list()
  vdj.10x_j_gene.match <- list()
  vdj.10x_c_gene.match <- list()
  for (i in 1:length(dat.tumor.vdj[[sname]]$vdj.10x_cdr3)){
    if (all(dat.tumor.vdj[[sname]]$vdj.10x_cdr3[[i]]!="None")){
      vdj.10x_cdr3.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_cdr3,function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_cdr3[[i]])
      }))
      vdj.10x_v_gene.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_v_gene,function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_v_gene[[i]])
      }))
      vdj.10x_d_gene.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_d_gene,function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_d_gene[[i]])
      }))
      vdj.10x_j_gene.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_j_gene,function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_j_gene[[i]])
      }))
      vdj.10x_c_gene.match[[i]] <- unlist(lapply(dat.til.vdj[[sname]]$vdj.10x_c_gene,function(x){
        any(x %in% dat.tumor.vdj[[sname]]$vdj.10x_c_gene[[i]])
      }))
      
      #til_matches[[i]] <- which(vdj.10x_cdr3.match[[i]])
      til_matches[[i]] <- which(vdj.10x_cdr3.match[[i]] & vdj.10x_v_gene.match[[i]] &
                                  vdj.10x_d_gene.match[[i]] & vdj.10x_j_gene.match[[i]] &
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

snames <- intersect(names(dat.tumor.vdj),
                    names(dat.til.vdj))
res <- lapply(snames,
              function(sname){get_matches(dat.tumor.vdj,dat.til.vdj,sname)})
names(res) <- snames

# dat.til@meta.data$tcr_matches_til <- F
# for (sname in snames){
#   matches <- res[[sname]]$til_matches[which(unlist(lapply(res[[sname]]$til_matches,length)) > 0)]
#   if (length(matches)>0){
#     for (i in 1:length(matches)){
#       idx <- matches[[i]]
#       dat.til@meta.data[dat.til@meta.data$UM.ID==sname,][idx,]$tcr_matches_til <- T
#     }
#   }
# }
# 
# table(dat.til@meta.data$tcr_matches_til,dat.til@meta.data$UM.ID)

dat.tumor@meta.data$tcr_matches_til <- F
for (sname in snames){
  dat.tumor@meta.data[dat.tumor@meta.data$UM.ID==sname,]$tcr_matches_til <- 
    unlist(lapply(res[[sname]]$til_matches,length)) > 0
}
table(dat.tumor@meta.data$tcr_matches_til,
      dat.tumor@meta.data$UM.ID)

# sname <- "UM1"
# i_match <- 3
# widx <- which(unlist(lapply(res[[sname]]$til_matches,function(x){length(x) > 0})))
# dat.tumor.vdj[[sname]]$vdj.10x_cdr3[widx[i_match]]
# dat.til.vdj[[sname]]$vdj.10x_cdr3[as.numeric(res[[sname]]$til_matches[widx][[i_match]])]
# 
# dat.tumor.vdj[[sname]]$vdj.10x_v_gene[widx[i_match]]
# dat.til.vdj[[sname]]$vdj.10x_v_gene[as.numeric(res[[sname]]$til_matches[widx][[i_match]])]
# 
# dat.tumor.vdj[[sname]]$vdj.10x_d_gene[widx[i_match]]
# dat.til.vdj[[sname]]$vdj.10x_d_gene[as.numeric(res[[sname]]$til_matches[widx][[i_match]])]
# 
# dat.tumor.vdj[[sname]]$vdj.10x_j_gene[widx[i_match]]
# dat.til.vdj[[sname]]$vdj.10x_j_gene[as.numeric(res[[sname]]$til_matches[widx][[i_match]])]
# 
# dat.tumor.vdj[[sname]]$vdj.10x_c_gene[widx[i_match]]
# dat.til.vdj[[sname]]$vdj.10x_c_gene[as.numeric(res[[sname]]$til_matches[widx][[i_match]])]
# 
# dat.tumor.vdj[[sname]]$vdj.10x_chains[widx[i_match]]
# dat.til.vdj[[sname]]$vdj.10x_chains[as.numeric(res[[sname]]$til_matches[widx][[i_match]])]
# 
# for (sname in names(dat.tumor.vdj)){
#   print(sname)
#   print(length(which(unlist(lapply(lapply(dat.tumor.vdj[[sname]]$vdj.10x_cdr3, function(x){
#     x %in% intersect(setdiff(unlist(dat.tumor.vdj[[sname]]$vdj.10x_cdr3),"None"),
#                      setdiff(unlist(dat.til.vdj[[sname]]$vdj.10x_cdr3),"None"))
#   }),function(x){length(which(x))})) > 0)))
# }

tcr_matches_til <- setNames(dat.tumor@meta.data$tcr_matches_til,
                            paste0(dat.tumor@meta.data$Sample.ID,"_",
                                   rownames(dat.tumor@meta.data)))

s.t <- AddMetaData(s.t, metadata = tcr_matches_til,col.name = "tcr_matches_til")

proportion_til <- list()
proportion_vasu <- list()
for (ident in unique(s.t@active.ident)){
  proportion_til[[ident]] <- sum(s.t@active.ident == ident & 
                                   s.t@meta.data$tcr_matches_til) / 
    sum(s.t@meta.data$tcr_matches_til)
  proportion_vasu[[ident]] <- sum(s.t@active.ident == ident & 
                                    s.t@meta.data$vasu_phenotype %in% c("MART1","CD3+41bb+")) / 
    sum(s.t@meta.data$vasu_phenotype %in% c("MART1","CD3+41bb+"))
}

df <- as.data.frame(cbind(unlist(proportion_til),unlist(proportion_vasu)))
colnames(df) <- c("proportion_of_tils","proportion_of_vasu")
df$cluster <- rownames(df)
df <- reshape2::melt(df)

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
  n_rolls <- sum(s.t@meta.data$vasu_phenotype %in% c("MART1","CD3+41bb+"))
  n_successes <- sum(s.t@meta.data$vasu_phenotype %in% c("MART1","CD3+41bb+") & s.t@active.ident == ident)
  bt[[ident]] <- binom.test(n_successes, n_rolls, expected_probability, alternative = "two.sided")
}

p_null <- unlist(lapply(bt,function(x){as.numeric(x$null.value)}))
p_estimate <- unlist(lapply(bt,function(x){as.numeric(x$estimate)}))
p <- unlist(lapply(bt,function(x){x$p.value}))
q <- p.adjust(p,method = "bonferroni") # probably not independent tests (the same dice rolls for all)

res <- data.frame(p_null,p_estimate,p,q,vasu_greater=p_estimate > p_null, sig = q < 0.05)
res
WriteXLS(res,ExcelFileName = "~/proj/um_ss/Manuscript/proportions_til_vs_vasu_matches_per_cluster.xlsx",
         AdjWidth = T,AutoFilter = T,BoldHeaderRow = T,row.names = T)

# Although, perhaps one would either need to subset to the set of samples included in vasus experiment
# or do it on a sample by sample basis. Or maybe not, since the overall question is whether TILs
# tend to come from exhauster clusters as a general principle. Might be worth further consideration.

df <- merge(df,res,by.x="cluster",by.y="row.names")

clusters <- unique(df$cluster)
df$is_max <- F
for (cluster in clusters){
  df$is_max[df$cluster == cluster][which.max(df$value[df$cluster == cluster])] <- T
}

DimPlot(s.t,group.by = "tcr_matches_til", order = T) + 
  DimPlot(s.t,group.by = "vasu_phenotype", order = T, cols = c("gray","gray","red","green")) + 
  ggplot(df,aes(x=cluster,y=value,fill=variable)) + 
  geom_bar(stat="identity",position = "dodge") + theme_classic() + 
  vertical_xlabels + xlab(NULL) + ylab("Proportion") + 
  geom_text(aes(label=ifelse(q < 0.05,"*","NS")), 
            data=subset(df,is_max),nudge_y = 0.02) + 
  ggtitle("Proportion of either that is present in a given cluster")
ggsave(filename = paste0("~/proj/um_ss/Manuscript/proportions_til_vs_vasu_matches_per_cluster.pdf"),
       width = 20,height = 7)


DimPlot(s.t,group.by = "tcr_matches_til", order = T,split.by = "UM.ID") +
  DimPlot(s.t,group.by = "vasu_phenotype", order = T, 
          cols = c("gray","gray","red","green"), split.by = "UM.ID")

DimPlot(s.t,group.by = "tcr_matches_til", order = T,split.by = "UM.ID")
ggsave(filename = "~/proj/um_ss/Manuscript/til_matches_biopsy.per_sample.pdf",
       width = 40,height = 6)

# Highlight those of Vasu's TCRs that match relevant things in TCR databases ---
library(readxl)

sheets <- excel_sheets("~/proj/um_ss/Investigations/db_matches.vasu.xlsx")
db_matches.vasu <- list()
for (sheet in sheets){
  db_matches.vasu[[sheet]] <- as.data.frame(
    read_excel(path = "~/proj/um_ss/Investigations/db_matches.vasu.xlsx",
               sheet = sheet),stringsAsFactors=F)
}

db_matches.vasu[["UM46"]]$takara.cdr3

head(s.t@meta.data)

DimPlot(s.t, group.by = "vasu_phenotype")


meta <- dat.tumor@meta.data
rownames(meta) <- paste0(meta$Sample.ID,"_",rownames(meta))

s.t <- AddMetaData(s.t,metadata = meta)
rm(meta)
gc()


get_cells_matching_mart1_vasu_and_db <- function(dat.tumor_til.vdj,db_matches.vasu,sname_1){
  
  um_id <- unique(dat.cca@meta.data$UM.ID[dat.cca@meta.data$Sample.ID==sname_1])
  
  idx_1 <- unlist(lapply(dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3, function(x){
    all(x!="None")
  }))
  
  matches_cdr3 <- lapply(dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3[which(idx_1)],function(x){
    if (any(x %in% db_matches.vasu[[um_id]]$takara.cdr3)){
      antigens <- list()
      for (cdr3 in intersect(x,db_matches.vasu[[um_id]]$takara.cdr3)){
        antigens[[cdr3]] <- db_matches.vasu[[um_id]][db_matches.vasu[[um_id]]$takara.cdr3 == cdr3,]
        antigens[[cdr3]]$tenX.cdr3 <- cdr3
      }
      antigens <- do.call("rbind",antigens)
      rownames(antigens) <- NULL
    } else {
      antigens <- ""
    }
    return(antigens)
  })
  
  matches_cdr3 <- matches_cdr3[which(!unlist(lapply(matches_cdr3,function(x){all(x=="")})))]
  if (length(matches_cdr3)>0){
    for (cell in names(matches_cdr3)){
      matches_cdr3[[cell]]$cell <- paste0(sname_1,"_",cell)
    }
    matches_cdr3 <- do.call("rbind",matches_cdr3)
    rownames(matches_cdr3) <- NULL
    matches_cdr3$db.Antigen[is.na(matches_cdr3$db.Antigen)] <- "Unknown"
    matches_cdr3 <- matches_cdr3[matches_cdr3$db.Antigen!="Unknown",]
    
    matches_cdr3.small <- unique(matches_cdr3[,c("tenX.cdr3","db.Antigen",
                                                 "db.Antigen_species_or_pathology","cell")])
    matches_cdr3.small <- matches_cdr3.small[matches_cdr3.small$db.Antigen == 
                                               "MART1",]
    
    cells_matching_mart1_db <- matches_cdr3.small$cell
    #cells_matching_mart1_db <- paste0(sname_1,"_",cells_matching_mart1_db)
  } else {
    cells_matching_mart1_db <- c()
  }
  return(cells_matching_mart1_db)
}

cells_matching_mart1_vasu_and_db <- list()
for (sname in names(dat.tumor_til.vdj)){
  cells_matching_mart1_vasu_and_db[[sname]] <- get_cells_matching_mart1_vasu_and_db(
    dat.tumor_til.vdj,db_matches.vasu,sname)
}

cells_matching_mart1_vasu_and_db <- unlist(cells_matching_mart1_vasu_and_db)
names(cells_matching_mart1_vasu_and_db) <- NULL

cells_matching_mart1_vasu_and_db <- setNames(rep("MART1",
                                                 length(cells_matching_mart1_vasu_and_db)),cells_matching_mart1_vasu_and_db)

s.t <- AddMetaData(s.t, metadata = cells_matching_mart1_vasu_and_db, 
                   col.name = "matches_mart1_vasu_and_db")

tmp <- setNames(as.character(s.t@meta.data$vasu_phenotype),rownames(s.t@meta.data))
tmp[which(s.t@meta.data$vasu_phenotype == "MART1" & 
            s.t@meta.data$matches_mart1_vasu_and_db == "MART1")] <- "MART1_db"
tmp <- factor(tmp,levels=c("None","Unidentified","CD3+41bb+","MART1","MART1_db"))

s.t <- AddMetaData(s.t,metadata = tmp, col.name = "vasu_phenotype_with_mart1_db")

pdf(file = "~/proj/um_ss/Manuscript/dimplot.vasu_phenotype_with_mart1_db.pdf",
    width = 6,height = 5)
DimPlot(s.t, group.by = "vasu_phenotype_with_mart1_db",order = T, 
        cols = c("gray","gray","red","green","purple"))
dev.off()

names(cells_matching_mart1_vasu_and_db) <- gsub("SampleID_[0-9]_11june18_","",
                                                names(cells_matching_mart1_vasu_and_db))

tils.unintegrated <- lapply(tils.unintegrated, function(x){
  if (any(names(cells_matching_mart1_vasu_and_db) %in% rownames(x@meta.data))){
    x <- AddMetaData(x,metadata = cells_matching_mart1_vasu_and_db, 
                     col.name = "matches_mart1_vasu_and_db")
  } else {
    x@meta.data$cells_matching_mart1_vasu_and_db <- NA
  }
  return(x)
})

for (sname in names(tils.unintegrated)){
  x <- tils.unintegrated[[sname]]
  tmp <- setNames(as.character(x@meta.data$vasu_phenotype),rownames(x@meta.data))
  tmp[which(x@meta.data$vasu_phenotype == "MART1" & 
              x@meta.data$matches_mart1_vasu_and_db == "MART1")] <- "MART1_db"
  tmp <- factor(tmp,levels=c("None","Unidentified","CD3+41bb+","MART1","MART1_db"))
  tils.unintegrated[[sname]] <- AddMetaData(x,metadata = tmp, 
                                            col.name = "vasu_phenotype_with_mart1_db")
}

DimPlot(tils.unintegrated[["UM22"]], group.by = "vasu_phenotype_with_mart1_db",
        order = T, cols = c("gray","gray","red","green","purple"))

# Alternative attempt to match TCRs --------------------------------------------

get_tcr_matrix <- function(s.t){
  cells <- rownames(s.t@meta.data)
  
  cdr3 <- lapply(s.t@meta.data$cdr3,function(x){unlist(strsplit(x,";"))})
  v_gene <- lapply(s.t@meta.data$v_gene,function(x){unlist(strsplit(x,";"))})
  j_gene <- lapply(s.t@meta.data$j_gene,function(x){unlist(strsplit(x,";"))})
  d_gene <- lapply(s.t@meta.data$d_gene,function(x){unlist(strsplit(x,";"))})
  c_gene <- lapply(s.t@meta.data$c_gene,function(x){unlist(strsplit(x,";"))})
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
  
  cols <- c(unique_cdr3, unique_v_gene, unique_j_gene, unique_d_gene, unique_c_gene)
  
  tcr <- matrix(NA, nrow = length(cells), ncol = length(cols), 
                dimnames = list(cells, cols))
  #  dim(tcr)
  
  library(foreach)
  library(doParallel)
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
  #tcr[is.na(tcr)] <- F
  tcr <- do.call("rbind",tcr_list)
  
  stopCluster(cl)
  
  return(tcr)
}

meta <- dat.cca@meta.data
rownames(meta) <- paste0(meta$Sample.ID, "_", rownames(meta))
s.t <- AddMetaData(s.t, metadata = meta)

tcr <- get_tcr_matrix(s.t)

saveRDS(tcr, file = "~/proj/um_ss/Investigations/s.t.tcr_matric.rda")

# Make TCR matrix for biopsies and TILs combined -------------------------------

dat.cca.til_tumor <- subset(dat.cca, subset = project.name %in% c("G18-023","GWA-JN-388"))
dat.cca.til_tumor.with_cdr3 <- subset(dat.cca.til_tumor,
                                      cells=rownames(dat.cca.til_tumor@meta.data[
                                        !is.na(dat.cca.til_tumor@meta.data$cdr3),]))

tcr_til_tumors <- get_tcr_matrix(dat.cca.til_tumor.with_cdr3)

saveRDS(tcr_til_tumors, file = "~/proj/um_ss/Investigations/til_tumor.tcr_matric.rda")

# Make UMAP plots for TILs that show TCRs matching biopsies --------------------

# Find all matches
tcr.matches <- tcr_til_tumors

cells <- rownames(dat.cca.til_tumor.with_cdr3@meta.data)
cdr3 <- lapply(dat.cca.til_tumor.with_cdr3@meta.data$cdr3,function(x){unlist(strsplit(x,";"))})
v_gene <- lapply(dat.cca.til_tumor.with_cdr3@meta.data$v_gene,function(x){unlist(strsplit(x,";"))})
j_gene <- lapply(dat.cca.til_tumor.with_cdr3@meta.data$j_gene,function(x){unlist(strsplit(x,";"))})
d_gene <- lapply(dat.cca.til_tumor.with_cdr3@meta.data$d_gene,function(x){unlist(strsplit(x,";"))})
c_gene <- lapply(dat.cca.til_tumor.with_cdr3@meta.data$c_gene,function(x){unlist(strsplit(x,";"))})
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

#colnames(tcr.matches[,is_v_gene])

cells_til <- rownames(dat.cca.til_tumor.with_cdr3@meta.data)[
  dat.cca.til_tumor.with_cdr3@meta.data$project.name == "G18-023"]

cells_tumor <- rownames(dat.cca.til_tumor.with_cdr3@meta.data)[
  dat.cca.til_tumor.with_cdr3@meta.data$project.name == "GWA-JN-388"]

idx_til <- rownames(tcr.matches) %in% cells_til
idx_tumor <- rownames(tcr.matches) %in% cells_tumor

tcr.matches.til <- tcr.matches[idx_til,]
tcr.matches.tumor <- tcr.matches[idx_tumor,]

library(foreach)
library(doParallel)
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
  #til_matches_tumor[[rownames(tcr.matches.tumor)[i]]] <- 
  rownames(tcr.matches.til)[idx_til_matches_tumor_i]
}
stopCluster(cl)

saveRDS(til_matches_tumor,file = "~/proj/um_ss/Investigations/til_matches_tumor.rda")

til_matches_tumor.bak <- til_matches_tumor

all_tils_matching_any_tumor <- sort(unique(unlist(til_matches_tumor)))

for (i in 1:length(tils.unintegrated)){
  cells <- rownames(tils.unintegrated[[i]]@meta.data)
  cells <- setNames(rep(F,length(cells)),cells)
  cells[names(cells) %in% all_tils_matching_any_tumor] <- T
  
  tils.unintegrated[[i]] <- AddMetaData(tils.unintegrated[[i]], metadata = cells, 
                                        col.name = "tcr_matching_any_biopsy")
}

p_um1 <- DimPlot(tils.unintegrated[["UM1"]],
                 group.by = "tcr_matching_any_biopsy", order = T, 
                 cols = c("gray", "blue")) + ggtitle("UM1")
p_um2 <- DimPlot(tils.unintegrated[["UM2"]], 
                 group.by = "tcr_matching_any_biopsy", order = T, 
                 cols = c("gray", "blue")) + ggtitle("UM2")
p_um3 <- DimPlot(tils.unintegrated[["UM3"]], 
                 group.by = "tcr_matching_any_biopsy", order = T, 
                 cols = c("gray", "blue")) + ggtitle("UM3")
p_um9 <- DimPlot(tils.unintegrated[["UM9"]], 
                 group.by = "tcr_matching_any_biopsy", order = T,
                 cols = c("gray", "blue")) + ggtitle("UM9")
p_um10 <- DimPlot(tils.unintegrated[["UM10"]], 
                  group.by = "tcr_matching_any_biopsy", order = T,
                  cols = c("gray", "blue")) + ggtitle("UM10")
p_um19 <- DimPlot(tils.unintegrated[["UM19"]], 
                  group.by = "tcr_matching_any_biopsy", order = T,
                  cols = c("gray", "blue")) + ggtitle("UM19")
p_um22 <- DimPlot(tils.unintegrated[["UM22"]], 
                  group.by = "tcr_matching_any_biopsy", order = T,
                  cols = c("gray", "blue")) + ggtitle("UM22")
p_um24 <- DimPlot(tils.unintegrated[["UM24"]], 
                  group.by = "tcr_matching_any_biopsy", order = T,
                  cols = c("gray", "blue")) + ggtitle("UM24")

library(gridExtra)
pdf(file = "~/proj/um_ss/Manuscript/til_tcr_matching_any_biopsy.pdf",width = 55, height = 6)
p <- grid.arrange(grobs = c(list(p_um1), 
                            list(p_um2),
                            list(p_um3),
                            list(p_um9),
                            list(p_um19),
                            list(p_um10),
                            list(p_um24),
                            list(p_um22)), ncol = 8, as.table = FALSE)
dev.off()



# Match Vasu sequences to TCR matrix, with different cutoffs -------------------
read_vdj.takara <- function(){
  library(readxl)
  library(stringr)
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
  
  return(vdj.takara)
}

vdj.takara <- read_vdj.takara()

# Define cutoff



# Find all matches
tcr.matches <- tcr

cells <- rownames(s.t@meta.data)
cdr3 <- lapply(s.t@meta.data$cdr3,function(x){unlist(strsplit(x,";"))})
v_gene <- lapply(s.t@meta.data$v_gene,function(x){unlist(strsplit(x,";"))})
j_gene <- lapply(s.t@meta.data$j_gene,function(x){unlist(strsplit(x,";"))})
d_gene <- lapply(s.t@meta.data$d_gene,function(x){unlist(strsplit(x,";"))})
c_gene <- lapply(s.t@meta.data$c_gene,function(x){unlist(strsplit(x,";"))})
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

is_cdr3_match <- colnames(tcr.matches) %in% vdj.takara$cdr3
is_v_gene_match <- colnames(tcr.matches) %in% vdj.takara$v_gene
is_d_gene_match <- colnames(tcr.matches) %in% vdj.takara$d_gene
is_j_gene_match <- colnames(tcr.matches) %in% vdj.takara$j_gene
is_c_gene_match <- colnames(tcr.matches) %in% vdj.takara$c_gene

is_cdr3 <- colnames(tcr.matches) %in% unique_cdr3
is_v_gene <- colnames(tcr.matches) %in% unique_v_gene
is_d_gene <- colnames(tcr.matches) %in% unique_d_gene
is_j_gene <- colnames(tcr.matches) %in% unique_j_gene
is_c_gene <- colnames(tcr.matches) %in% unique_c_gene

is_cdr3_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$cdr3[vdj.takara$pass]
is_v_gene_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$v_gene[vdj.takara$pass]
is_d_gene_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$d_gene[vdj.takara$pass]
is_j_gene_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$j_gene[vdj.takara$pass]
is_c_gene_match_cutoff <- colnames(tcr.matches) %in% vdj.takara$c_gene[vdj.takara$pass]

#rowSums(tcr.matches.cutoff.cdr3)

get_matching_cells_rank <- function(tcr.matches,vdj.takara,rank){
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
    cells.matching.cdr3 <- names(which(rowSums(tcr.matches[,is_cdr3_match & is_cdr3_match_cutoff])>0))
  }
  if (!any(class(tcr.matches.v_gene)=="matrix")){
    cells.matching.v_gene <- names(which(tcr.matches.v_gene == T))
  } else {
    cells.matching.v_gene <- names(which(rowSums(tcr.matches[,is_v_gene_match & is_v_gene_match_cutoff])>0))
  }
  if (!any(class(tcr.matches.j_gene)=="matrix")){
    cells.matching.j_gene <- names(which(tcr.matches.j_gene == T))
  } else {
    cells.matching.j_gene <- names(which(rowSums(tcr.matches[,is_j_gene_match & is_j_gene_match_cutoff])>0))
  }
  
  cells_matching <- intersect(intersect(
    cells.matching.cdr3,
    cells.matching.v_gene),
    cells.matching.j_gene)
  
  return(cells_matching)
}

ranks <- lapply(split(vdj.takara,paste0(vdj.takara$UM.ID,":",vdj.takara$chains)),function(x){max(x$reads.rank)})
max_rank <- max(unlist(ranks))

get_cells_matching <- function(vdj.takara,sname,chain){
  #sname <- "C7-GEX"
  um_id <- unique(s.t@meta.data$UM.ID[s.t@meta.data$Sample.ID==sname])
  #chain <- "TRA"
  ranks <- lapply(split(vdj.takara,paste0(vdj.takara$UM.ID,":",vdj.takara$chains)),function(x){max(x$reads.rank)})
  max_rank <- ranks[[paste0(um_id,":",chain)]]
  
  cells_matching <- list()
  for (rank in 1:max_rank){
    cells_matching[[rank]] <- get_matching_cells_rank(
      tcr.matches[grep(sname,rownames(tcr.matches)),],
      vdj.takara[vdj.takara$UM.ID==um_id & vdj.takara$chains==chain,],
      rank)
  }
  
  return(cells_matching)
}

cells_matching <- list()
for (um_id in unique(vdj.takara$UM.ID)){
  if (um_id %in% s.t@meta.data$UM.ID){
    sname <- unique(s.t@meta.data$Sample.ID[s.t@meta.data$UM.ID==um_id])
    for (chain in c("TRA","TRB")){
      cells_matching[[paste0(um_id,":",chain)]] <- get_cells_matching(vdj.takara,sname,chain)
    }
  }
}

#n_matches_rank <- data.frame(rank=1:max_rank,n_matches=unlist(lapply(cells_matching,length)))

get_cells_matching_per_cluster <- function(cells_matching){
  n_matches_rank.per_cluster <- do.call("rbind",lapply(cells_matching,function(x){
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
  n_matches_rank.per_cluster <- reshape2::melt(n_matches_rank.per_cluster,id.vars="rank")
  
  return(n_matches_rank.per_cluster)
}

n_matches_rank.per_cluster <- list()
for (um_id_chain in names(cells_matching)){
  n_matches_rank.per_cluster[[um_id_chain]] <- get_cells_matching_per_cluster(cells_matching[[um_id_chain]])
  n_matches_rank.per_cluster[[um_id_chain]]$um_id <- str_split_fixed(um_id_chain,":",2)[,1]
  n_matches_rank.per_cluster[[um_id_chain]]$chain <- str_split_fixed(um_id_chain,":",2)[,2]
}
n_matches_rank.per_cluster <- do.call("rbind",n_matches_rank.per_cluster)


get_cells_per_cluster <- function(um_id){
  tab <- table(s.t@active.ident[s.t@meta.data$UM.ID==um_id])
  tab <- as.data.frame.matrix(t(as.matrix(tab)))
  return(tab)
}

cells_per_cluster <- list()
for (um_id in unique(n_matches_rank.per_cluster$um_id)){
  cells_per_cluster[[um_id]] <- get_cells_per_cluster(um_id)
}
cells_per_cluster <- do.call("rbind",cells_per_cluster)
cells_per_cluster$um_id <- rownames(cells_per_cluster)
cells_per_cluster <- reshape2::melt(cells_per_cluster)

# n_matches_rank.per_sample <- do.call("rbind",lapply(cells_matching,function(x){
#   samples <- str_split_fixed(x,"_",2)[,1]
#   tab <- table(samples)
#   
#   missing <- setdiff(unique(s.t@meta.data$Sample.ID),names(tab))
#   tab <- as.data.frame.matrix(t(as.matrix(tab)))
#   if (length(missing) > 0){
#     tab <- cbind(tab,t(rep(0,length(missing))))
#     colnames(tab)[(ncol(tab)-length(missing)+1):ncol(tab)] <- missing
#   }
#   tab <- tab[,unique(s.t@meta.data$Sample.ID)]
#   tab
# }))
# n_matches_rank.per_sample$rank <- as.numeric(rownames(n_matches_rank.per_sample))
# n_matches_rank.per_sample <- reshape2::melt(n_matches_rank.per_sample,id.vars="rank")

# ggplot(n_matches_rank,aes(x = rank, y = n_matches)) + 
#   geom_bar(stat = "identity", position = "dodge")

n_matches_rank.per_cluster$variable <- factor(as.character(n_matches_rank.per_cluster$variable),
                                              levels = c("Exhausted / dysfunctional",
                                                         "Late activated, proliferative",
                                                         "Late activated",
                                                         "Activated, cytotoxic (FCGR3A+)",
                                                         "Early activated (KLRC2+, CD38+)",
                                                         "Naive-like / early activated",
                                                         "Naive-like"))

ggplot(n_matches_rank.per_cluster[n_matches_rank.per_cluster$chain=="TRA",],
       aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "stack") + theme_classic() + 
  facet_wrap(~ um_id, scales = "free") + theme(strip.background = element_blank()) + 
  xlab(NULL) + ylab("Number of cells matching TCR")
ggsave(filename = "~/proj/um_ss/Manuscript/vasu_tcr.n_matches_rank.per_cluster.tra.pdf",
       width = 15, height = 3)

ggplot(n_matches_rank.per_cluster[n_matches_rank.per_cluster$chain=="TRB",],
       aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "stack") + theme_classic() + 
  facet_wrap(~ um_id, scales = "free") + theme(strip.background = element_blank()) + 
  xlab(NULL) + ylab("Number of cells matching TCR")
ggsave(filename = "~/proj/um_ss/Manuscript/vasu_tcr.n_matches_rank.per_cluster.trb.pdf",
       width = 15, height = 3)

cells_per_cluster$variable <- factor(as.character(cells_per_cluster$variable),
                                     levels = c("Exhausted / dysfunctional",
                                                "Late activated, proliferative",
                                                "Late activated",
                                                "Activated, cytotoxic (FCGR3A+)",
                                                "Early activated (KLRC2+, CD38+)",
                                                "Naive-like / early activated",
                                                "Naive-like"))

ggplot(cells_per_cluster,
       aes(x = um_id, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "stack") + theme_classic() + 
  facet_wrap(~ um_id, scales = "free") + theme(strip.background = element_blank()) + 
  xlab(NULL) + ylab("Number of cells")
ggsave(filename = "~/proj/um_ss/Manuscript/cells_per_cluster.pdf",
       width = 5, height = 3) 

# ggplot(n_matches_rank.per_sample,aes(x = rank, y = value, fill = variable)) + 
#   geom_bar(stat = "identity", position = "stack")

# Get matches to TCRs that are the same between UM22 and UM46

vdj.takara.um22 <- vdj.takara[vdj.takara$UM.ID=="UM22",]#c("cdr3","v_gene","j_gene")]
vdj.takara.um46 <- vdj.takara[vdj.takara$UM.ID=="UM46",]#c("cdr3","v_gene","j_gene")]
vdj.takara.um22$cdr3_v_j <- apply(vdj.takara.um22[,c("cdr3","v_gene","j_gene")],1,paste0,collapse=":")
vdj.takara.um46$cdr3_v_j <- apply(vdj.takara.um46[,c("cdr3","v_gene","j_gene")],1,paste0,collapse=":")
shared_cdr3_v_j <- intersect(vdj.takara.um22$cdr3_v_j,
                             vdj.takara.um46$cdr3_v_j)
vdj.takara.um22.shared <- vdj.takara.um22[vdj.takara.um22$cdr3_v_j %in% shared_cdr3_v_j,]
vdj.takara.um46.shared <- vdj.takara.um46[vdj.takara.um46$cdr3_v_j %in% shared_cdr3_v_j,]

vdj.takara.um46.shared$reads.rank[vdj.takara.um46.shared$chains=="TRA"] <- 1:length(vdj.takara.um46.shared[vdj.takara.um46.shared$chains=="TRA",]$reads.rank)
vdj.takara.um46.shared$reads.rank[vdj.takara.um46.shared$chains=="TRB"] <- 1:length(vdj.takara.um46.shared[vdj.takara.um46.shared$chains=="TRB",]$reads.rank)

sname <- unique(s.t@meta.data$Sample.ID[s.t@meta.data$UM.ID=="UM46"])
cells_matching.um46.shared <- list()
for (chain in c("TRA","TRB")){
  cells_matching.um46.shared[[chain]] <- get_cells_matching(vdj.takara.um46.shared,sname,chain)
}

n_matches_rank.per_cluster.um46.shared <- list()
for (chain in names(cells_matching.um46.shared)){
  n_matches_rank.per_cluster.um46.shared[[chain]] <- get_cells_matching_per_cluster(cells_matching.um46.shared[[chain]])
  n_matches_rank.per_cluster.um46.shared[[chain]]$um_id <- "UM46"
  n_matches_rank.per_cluster.um46.shared[[chain]]$chain <- chain
}
n_matches_rank.per_cluster.um46.shared <- do.call("rbind",n_matches_rank.per_cluster.um46.shared)

n_matches_rank.per_cluster.um46.shared$variable <- factor(as.character(n_matches_rank.per_cluster.um46.shared$variable),
                                                          levels = c("Exhausted / dysfunctional",
                                                                     "Late activated, proliferative",
                                                                     "Late activated",
                                                                     "Activated, cytotoxic (FCGR3A+)",
                                                                     "Early activated (KLRC2+, CD38+)",
                                                                     "Naive-like / early activated",
                                                                     "Naive-like"))

ggplot(n_matches_rank.per_cluster.um46.shared[
  n_matches_rank.per_cluster.um46.shared$chain=="TRA",],
  aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "stack") + theme_classic() + 
  facet_wrap(~ um_id, scales = "free") + theme(strip.background = element_blank()) + 
  xlab(NULL) + ylab("Number of cells matching TCR")
ggsave(filename = "~/proj/um_ss/Manuscript/vasu_tcr.n_matches_rank.per_cluster.um46.shared.tra.pdf",
       width = 6, height = 3)

ggplot(n_matches_rank.per_cluster.um46.shared[
  n_matches_rank.per_cluster.um46.shared$chain=="TRB",],
  aes(x = rank, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "stack") + theme_classic() + 
  facet_wrap(~ um_id, scales = "free") + theme(strip.background = element_blank()) + 
  xlab(NULL) + ylab("Number of cells matching TCR")
ggsave(filename = "~/proj/um_ss/Manuscript/vasu_tcr.n_matches_rank.per_cluster.um46.shared.trb.pdf",
       width = 6, height = 3)

# Sample and cluster proportions -----------------------------------------------

make_pies2 <- function(tab){
  library(ggplot2)
  
  tab <- as.data.frame.matrix(tab)
  
  for (i in 1:nrow(tab)){
    tab[i,] <- tab[i,]/sum(tab[i,])
  }
  tab <- as.data.frame(t(tab))
  tab$cluster <- rownames(tab)
  tab.melt <- reshape2::melt(tab)
  
  print(ggplot(tab.melt, aes(x = 1, y = value, fill = cluster)) +
          facet_grid(cols = vars(variable)) + 
          coord_polar(theta = "y") +
          geom_col(position = position_stack(reverse = TRUE), 
                   size = 3, show.legend = T) + 
          theme_void())
  
  return(tab)
}

plot_cluster_proportion_correlation <- function(tab, correction="bonferroni"){
  
  p <- matrix(NA,nrow=ncol(tab),ncol=ncol(tab))
  for (i in 1:ncol(tab)){
    for (j in 1:ncol(tab)){
      p[i,j] <- cor.test(tab[,i],tab[,j],method="spearman",exact = F)$p.value
    }
  }
  p[lower.tri(p)] <- NA
  diag(p) <- NA
  p.uncorrected <- p
  p <- matrix(p.adjust(p,method=correction),ncol(tab),ncol(tab))
  
  for (i in 1:ncol(p)){
    for (j in 1:ncol(p)){
      p[j,i] <- p[i,j]
    }
  }
  
  p.sig <- ifelse(p < 0.05,"*","")
  diag(p.sig) <- ""
  
  cr <- cor(tab,method = "spearman")
  pl <- pheatmap(cr,border_color = NA,display_numbers = p.sig,
                 cluster_rows = T,cluster_cols = T)
  return(list(pl=pl,
              cr=cr,
              p.sig=p.sig,
              q=p,
              p=p.uncorrected))
}

pdf(file = "~/proj/um_ss/Manuscript/cluster_sample_contribution.pdf",width = 15,height = 7)
DimPlot(s, group.by = "UM.ID") + DimPlot(s.t,group.by = "UM.ID")
dev.off()

tab <- table(s.t@meta.data$UM.ID,s.t@active.ident)
tab <- as.data.frame.matrix(tab)
for (i in 1:nrow(tab)){
  tab[i,] <- tab[i,]/sum(tab[i,])
}
pdf(file = "~/proj/um_ss/Manuscript/cluster_proportion_correlation.t.pdf",width = 5.5,height = 5)
cluster_proportion_correlation.t <- plot_cluster_proportion_correlation(tab,correction = "BH")
dev.off()

tab <- table(s@meta.data$UM.ID,s@active.ident)
tab <- tab[,setdiff(colnames(tab),c("Melanoma","Melanoma (CPXM1+)"))]
tab <- as.data.frame.matrix(tab)
for (i in 1:nrow(tab)){
  tab[i,] <- tab[i,]/sum(tab[i,])
}
pdf(file = "~/proj/um_ss/Manuscript/cluster_proportion_correlation.all_no_melanoma.pdf",
    width = 5.5,height = 5)
cluster_proportion_correlation.all_no_melanoma <- 
  plot_cluster_proportion_correlation(tab,correction = "BH")
dev.off()

tab <- table(s.t@meta.data$UM.ID,s.t@active.ident)
make_pies2(tab)
ggsave(filename = "~/proj/um_ss/Manuscript/cluster_sample_contribution.t.pie2.pdf",
       width = 15,height = 3)

tab <- table(s@meta.data$UM.ID,s@active.ident)
make_pies2(tab)
ggsave(filename = "~/proj/um_ss/Manuscript/cluster_sample_contribution.all.pie2.pdf",
       width = 15,height = 3)

# Pair-wise overlap in TCRs among samples --------------------------------------
get_matches_pair.v3 <- function(dat.tumor_til.vdj,sname_1,sname_2,match_on = c("cdr3","v","d","j","c")){
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
  
  
  which_sname_2.10x_cdr3_in_sname_1.10x_cdr3 <- unlist(lapply(tmp_1.10x_cdr3,function(x){
    which(unlist(lapply(tmp_2.10x_cdr3,function(y){any(y %in% x)})))
  })) 
  which_sname_2.10x_v_gene_in_sname_1.10x_v_gene <- unlist(lapply(tmp_1.10x_v_gene,function(x){
    which(unlist(lapply(tmp_2.10x_v_gene,function(y){any(y %in% x)})))
  })) 
  which_sname_2.10x_d_gene_in_sname_1.10x_d_gene <- unlist(lapply(tmp_1.10x_d_gene,function(x){
    which(unlist(lapply(tmp_2.10x_d_gene,function(y){any(y %in% x)})))
  })) 
  which_sname_2.10x_j_gene_in_sname_1.10x_j_gene <- unlist(lapply(tmp_1.10x_j_gene,function(x){
    which(unlist(lapply(tmp_2.10x_j_gene,function(y){any(y %in% x)})))
  }))
  which_sname_2.10x_c_gene_in_sname_1.10x_c_gene <- unlist(lapply(tmp_1.10x_c_gene,function(x){
    which(unlist(lapply(tmp_2.10x_c_gene,function(y){any(y %in% x)})))
  }))
  
  all_pairs <- c(names(which_sname_2.10x_cdr3_in_sname_1.10x_cdr3),
                 names(which_sname_2.10x_v_gene_in_sname_1.10x_v_gene),
                 names(which_sname_2.10x_d_gene_in_sname_1.10x_d_gene),
                 names(which_sname_2.10x_j_gene_in_sname_1.10x_j_gene),
                 names(which_sname_2.10x_c_gene_in_sname_1.10x_c_gene))
  
  unique_pairs <- unique(all_pairs)
  
  df <- data.frame(pair = unique_pairs,
                   in_10x_cdr3 = unique_pairs %in% names(which_sname_2.10x_cdr3_in_sname_1.10x_cdr3),
                   in_10x_v_gene = unique_pairs %in% names(which_sname_2.10x_v_gene_in_sname_1.10x_v_gene),
                   in_10x_d_gene = unique_pairs %in% names(which_sname_2.10x_d_gene_in_sname_1.10x_d_gene),
                   in_10x_j_gene = unique_pairs %in% names(which_sname_2.10x_j_gene_in_sname_1.10x_j_gene),
                   in_10x_c_gene = unique_pairs %in% names(which_sname_2.10x_c_gene_in_sname_1.10x_c_gene),
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
    df$match <- df$match & df$in_v_gene
  }
  if ("d" %in% match_on){
    df$match <- df$match & df$in_d_gene
  }
  if ("j" %in% match_on){
    df$match <- df$match & df$in_j_gene
  }
  if ("c" %in% match_on){
    df$match <- df$match & df$in_c_gene
  }
  
  #df$match <- df$in_10x_cdr3 & df$in_10x_v_gene & df$in_10x_d_gene & df$in_10x_j_gene & df$in_10x_c_gene
  df$cell_1 <- str_split_fixed(df$pair,"\\.",2)[,1]
  df$cell_2 <- str_split_fixed(df$pair,"\\.",2)[,2]
  
  df.match <- df[df$match,]
  
  if (nrow(df.match)>0){
    for (i in 1:nrow(df.match)){
      cell_1 <- df.match$cell_1[i]
      cell_2 <- df.match$cell_2[i]
      df.match$shared_cdr3[i] <- paste0(sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_cdr3[[cell_1]],
                                                       dat.tumor_til.vdj[[sname_2]]$vdj.10x_cdr3[[cell_2]])),collapse = ";")
      df.match$shared_v_gene[i] <- paste0(sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_v_gene3[[cell_1]],
                                                         dat.tumor_til.vdj[[sname_2]]$vdj.10x_v_gene[[cell_2]])),collapse = ";")
      df.match$shared_d_gene[i] <- paste0(sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_d_gene[[cell_1]],
                                                         dat.tumor_til.vdj[[sname_2]]$vdj.10x_d_gene[[cell_2]])),collapse = ";")
      df.match$shared_j_gene[i] <- paste0(sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_j_gene[[cell_1]],
                                                         dat.tumor_til.vdj[[sname_2]]$vdj.10x_j_gene[[cell_2]])),collapse = ";")
      df.match$shared_c_gene[i] <- paste0(sort(intersect(dat.tumor_til.vdj[[sname_1]]$vdj.10x_c_gene[[cell_1]],
                                                         dat.tumor_til.vdj[[sname_2]]$vdj.10x_c_gene[[cell_2]])),collapse = ";")
    }
  }
  return(df.match)
}

# dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")
# 
# s <- readRDS(file = "~/proj/um_ss/Manuscript/s.rda")
# s.t <- readRDS(file = "~/proj/um_ss/Manuscript/s.t.rda")
# vasu_phenotype <- readRDS("~/proj/um_ss/Investigations/seurat/results/liger_all/vasu_phenotype.rda")
# 
# dat.cca.vasu_phenotype <- vasu_phenotype[intersect(names(vasu_phenotype),
#                                                    rownames(dat.cca@meta.data))]
# dat.cca <- AddMetaData(dat.cca, metadata = dat.cca.vasu_phenotype, col.name = "vasu_phenotype")
# table(dat.cca@meta.data$vasu_phenotype, dat.cca@meta.data$UM.ID)
# 
# dat.cca@meta.data$vasu_phenotype <- factor(dat.cca@meta.data$vasu_phenotype,
#                                            levels=c("None","Unidentified","CD3+41bb+","MART1"))
# 
# UM.ID <- setNames(dat.cca@meta.data$UM.ID,paste0(dat.cca@meta.data$Sample.ID,"_",rownames(dat.cca@meta.data)))
# s.t <- AddMetaData(s.t,metadata = UM.ID,col.name = "UM.ID")
# 
# s.t$UM.ID <- factor(s.t$UM.ID,levels=paste0("UM",sort(as.numeric(gsub("UM","",unique(s.t@meta.data$UM.ID))))))
# levels(s.t$UM.ID)
# 
# UM.ID <- setNames(dat.cca@meta.data$UM.ID,paste0(dat.cca@meta.data$Sample.ID,"_",rownames(dat.cca@meta.data)))
# s <- AddMetaData(s,metadata = UM.ID,col.name = "UM.ID")
# s$UM.ID <- factor(s$UM.ID,levels=paste0("UM",sort(as.numeric(gsub("UM","",unique(s@meta.data$UM.ID))))))
# levels(s$UM.ID)

#dat.tumor_til <- subset(dat.cca,project.name %in% c("G18-023","GWA-JN-388"))
#dat.tumor_til.split <- SplitObject(dat.tumor_til, split.by = "Sample.ID")
#dat.tumor_til.vdj <- lapply(dat.tumor_til.split,extract_tcr)
dat.tumor_til.vdj <- readRDS(file="~/proj/um_ss/Manuscript/dat.tumor_til.vdj.rda")

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
#saveRDS(matches_pair,file = "~/proj/um_ss/Investigations/matches_pair.rda")
#matches_pair <- readRDS(file = "~/proj/um_ss/Investigations/matches_pair.rda")

matches_pair.clones <- lapply(matches_pair,function(x){
  unique(x[c("shared_cdr3","shared_v_gene","shared_d_gene","shared_j_gene","shared_c_gene")])  
})

matches_pair.nr <- as.data.frame(unlist(lapply(matches_pair.clones,nrow)))
colnames(matches_pair.nr) <- "nr"
matches_pair.nr$sample_1 <- str_split_fixed(rownames(matches_pair.nr),":",2)[,1]
matches_pair.nr$sample_2 <- str_split_fixed(rownames(matches_pair.nr),":",2)[,2]
rownames(matches_pair.nr) <- NULL

unique_samples <- c("SampleID_7_11june18","SampleID_8_11june18","SampleID_6_11june18",
                    "SampleID_9_11june18","SampleID_1_11june18","SampleID_4_11june18",
                    "SampleID_2_11june18","SampleID_3_11june18",
                    "D4-GEX","D6-GEX","B3-GEX","C7-GEX","D2-GEX","C8-GEX",
                    "C5-GEX","A2-GEX")
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

dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")

annot <- dat.cca@meta.data[,c("Sample.ID","UM.ID","project.name")]
annot <- unique(annot)
rownames(annot) <- annot$Sample.ID
annot$Sample.ID <- NULL
annot <- annot[annot$project.name!="G18-049",]
annot$UM.ID <- factor(annot$UM.ID,
                      levels=paste0("UM",sort(as.numeric(gsub("UM","",unique(annot$UM.ID))))))

rownames(annot)[annot$project.name=="G18-023"][
  order(as.numeric(gsub("UM","",unique(as.character(annot$UM.ID)))))]

annot$UM.ID[annot$project.name=="G18-023"][
  order(as.numeric(gsub("UM","",unique(as.character(annot$UM.ID)))))]

pdf(file="~/proj/um_ss/Manuscript/tcr_overlap.pdf", width = 7, height = 6)
pheatmap(matches_pair.nr.m, cluster_rows = F, cluster_cols = F,
         display_numbers = matches_pair.nr.m, border_color = NA,
         annotation_col = annot, annotation_row = annot, show_rownames = F,
         show_colnames = F, na_col = "white", colorRampPalette(c("white", "red"))(100))
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
#saveRDS(matches_pair,file = "~/proj/um_ss/Investigations/matches_pair.rda")
#matches_pair <- readRDS(file = "~/proj/um_ss/Investigations/matches_pair.rda")

matches_pair.bak <- matches_pair

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
  unique(x[c("shared_cdr3","shared_v_gene","shared_d_gene","shared_j_gene","shared_c_gene")])  
})

matches_pair.nr <- as.data.frame(unlist(lapply(matches_pair.clones,nrow)))
colnames(matches_pair.nr) <- "nr"
matches_pair.nr$sample_1 <- str_split_fixed(rownames(matches_pair.nr),":",2)[,1]
matches_pair.nr$sample_2 <- str_split_fixed(rownames(matches_pair.nr),":",2)[,2]
rownames(matches_pair.nr) <- NULL

unique_samples <- c("SampleID_7_11june18","SampleID_8_11june18","SampleID_6_11june18",
                    "SampleID_9_11june18","SampleID_1_11june18","SampleID_4_11june18",
                    "SampleID_2_11june18","SampleID_3_11june18",
                    "D4-GEX","D6-GEX","B3-GEX","C7-GEX","D2-GEX","C8-GEX",
                    "C5-GEX","A2-GEX")
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

dat.cca <- readRDS(file="~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")

annot <- dat.cca@meta.data[,c("Sample.ID","UM.ID","project.name")]
annot <- unique(annot)
rownames(annot) <- annot$Sample.ID
annot$Sample.ID <- NULL
annot <- annot[annot$project.name!="G18-049",]
annot$UM.ID <- factor(annot$UM.ID,
                      levels=paste0("UM",sort(as.numeric(gsub("UM","",unique(annot$UM.ID))))))

rownames(annot)[annot$project.name=="G18-023"][
  order(as.numeric(gsub("UM","",unique(as.character(annot$UM.ID)))))]

annot$UM.ID[annot$project.name=="G18-023"][
  order(as.numeric(gsub("UM","",unique(as.character(annot$UM.ID)))))]

pdf(file="~/proj/um_ss/Manuscript/tcr_overlap.cdr3.vasu.pdf",width = 7,height = 6)
pheatmap(matches_pair.nr.m,cluster_rows = F,cluster_cols = F,
         display_numbers = matches_pair.nr.m,border_color = NA,
         annotation_col = annot,annotation_row = annot,show_rownames = F,
         show_colnames = F,na_col = "white",colorRampPalette(c("white", "red"))(100))
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

pdf(file = "~/proj/um_ss/Manuscript/overlap_tcr_vasu_tra.pdf", width = 3, height = 2.5)
pheatmap(overlap_tcr_vasu_tra,cluster_rows = F,cluster_cols = F,
         display_numbers = overlap_tcr_vasu_tra,border_color = NA,
         show_rownames = T, show_colnames = T, 
         na_col = "white",colorRampPalette(c("white", "red"))(100))
dev.off()

pdf(file = "~/proj/um_ss/Manuscript/overlap_tcr_vasu_trb.pdf", width = 3, height = 2.5)
pheatmap(overlap_tcr_vasu_trb,cluster_rows = F,cluster_cols = F,
         display_numbers = overlap_tcr_vasu_trb,border_color = NA,
         show_rownames = T, show_colnames = T, 
         na_col = "white",colorRampPalette(c("white", "red"))(100))
dev.off()