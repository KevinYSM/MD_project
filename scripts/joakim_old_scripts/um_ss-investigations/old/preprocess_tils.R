library(Seurat)
library(future)
library(tidyverse)
library(djvdj)
library(plyranges)
library(rtracklayer)
library(DoubletFinder)
library(pheatmap)

#plan("multiprocess", workers = 63)
#plan("multiprocess", workers = 15)
#options(future.globals.maxSize = 8000 * 1024^2)
options(future.globals.maxSize = 80000 * 1024^2)

outdir <- "~/proj/um_ss/Investigations/seurat/results/tils/"
dir.create(outdir)

dat.cca <- readRDS("~/proj/um_ss/Investigations/seurat/results/dat.integrated.doubletfiltered.reprocecessed.filtered.reprocessed.rda")

tils <- list()
for (sname in unique(dat.cca@meta.data$Sample.ID[dat.cca@meta.data$project.name == "G18-023"])){
  dat <- subset(dat.cca, Sample.ID == sname)
  meta <- dat.cca@meta.data[dat.cca$Sample.ID == sname,]
  dat <- GetAssayData(dat, slot = "counts", assay = "RNA")
  dat <- CreateSeuratObject(counts = dat, meta.data = meta)
  dat <- NormalizeData(dat,verbose = T)
  dat <- FindVariableFeatures(dat, selection.method = "vst", 
                              nfeatures = 4000)
  dat <- CellCycleScoring(dat, s.features = cc.genes$s.genes, 
                          g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  dat <- ScaleData(dat, verbose = T, vars.to.regress = c("S.Score", "G2M.Score"))
  #dat <- ScaleData(dat, verbose = T)
  dat <- RunPCA(dat, pc.genes = VariableFeatures(dat), 
                npcs = 30, verbose = T)
  dat <- RunUMAP(dat, dims = 1:30)
  dat <- FindNeighbors(dat, dims = 1:30)
  dat <- FindClusters(dat, resolution = 0.5)
  tils[[sname]] <- dat
}
tils.bak <- tils

sname <- names(tils)[3]
dat <- tils[[sname]]

# For modularity.fxn=2, algorithm = 4
resolutions <- list()
resolutions[["SampleID_1_11june18"]] <- 0.1
resolutions[["SampleID_2_11june18"]] <- 0.4
resolutions[["SampleID_3_11june18"]] <- 0.8
resolutions[["SampleID_4_11june18"]] <- 6.6
resolutions[["SampleID_6_11june18"]] <- 5.2
resolutions[["SampleID_7_11june18"]] <- 0.6
resolutions[["SampleID_8_11june18"]] <- 0.1
resolutions[["SampleID_9_11june18"]] <- 0.4

dat <- FindClusters(dat, resolution = resolutions[[sname]], modularity.fxn=2, algorithm = 4)
pdf(file=paste0(outdir,sname,".umap.pdf"),width = 10,height = 10)
DimPlot(dat,
        reduction = "umap", 
        label = T)
dev.off()

pdf(file=paste0(outdir,sname,".features.pdf"),width = 30,height = 20)
FeaturePlot(dat, features = c("CD3D","CD3G","CD3E","CD8A","CD8B","CD4","NCAM1",
                              "FCGR3A","KLRC1","NCR1"),order = T)
dev.off()

# FeaturePlot(dat, features = union(grep("^TRGV",rownames(dat@assays$RNA@data),value = T),
#                                   grep("^TRDV",rownames(dat@assays$RNA@data),value = T)),
#             min.cutoff = 1,slot = "counts",order = T)

pdf(file=paste0(outdir,sname,".features.TRGV9.TRDV2.pdf"),width = 10,height = 5)
FeaturePlot(dat, features = c("TRGV9","TRDV2"),
            min.cutoff = 1,slot = "counts",order = T)
dev.off()


# markers <- FindAllMarkers(dat,assay = "RNA",slot = "data",
#                           latent.vars = c("S.Score", "G2M.Score"),
#                           test.use = "LR", only.pos=T)

#markers <- FindAllMarkers(dat,assay = "RNA",slot = "data", only.pos=T)

# 1
# dat <- RenameIdents(dat,
#              `1` = "CD4 T cells",
#              `2` = "NK cells",
#              `3` = "CD8 T cells")

#2
# dat <- RenameIdents(dat,
#                     `1` = "NK cells",
#                     `2` = "CD8 T cells",
#                     `3` = "CD8 T cells",
#                     `4` = "CD4 T cells",
#                     `5` = "CD8 T cells",
#                     `6` = "CD8 T cells",
#                     `7` = "CD8 T cells",
#                     `8` = "CD8 T / NK cells")

#3
# dat <- RenameIdents(dat,
#                     `1` = "CD4 T cells",
#                     `2` = "CD4 T cells",
#                     `3` = "CD4 T cells",
#                     `4` = "Gamma-delta T cells (TRGV9, TRDV2)",
#                     `5` = "CD4 T cells",
#                     `6` = "CD4 T cells",
#                     `7` = "CD4 T cells",
#                     `8` = "CD8 T cells",
#                     `9` = "CD4 T cells",
#                     `10` = "CD4 T cells",
#                     `11` = "CD4 T cells",
#                     `12` = "CD8 T cells",
#                     `13` = "NK cells",
#                     `14` = "CD4 T cells")

#4
# dat <- RenameIdents(dat,
#                     `1` = "CD4 T cells",
#                     `2` = "CD4 T cells",
#                     `3` = "CD4 T cells",
#                     `4` = "CD4 T cells",
#                     `5` = "CD4 T cells",
#                     `6` = "CD4 T cells",
#                     `7` = "CD4 T cells",
#                     `8` = "CD4 T cells",
#                     `9` = "CD4 T cells",
#                     `10` = "CD4 T cells",
#                     `11` = "CD4 T cells",
#                     `12` = "CD4 T cells",
#                     `13` = "CD4 T cells",
#                     `14` = "CD8 T cells",
#                     `15` = "CD4 T cells",
#                     `16` = "CD4 T cells",
#                     `17` = "CD4 T cells",
#                     `18` = "CD4 T cells",
#                     `19` = "CD4 T cells",
#                     `20` = "CD4 T cells",
#                     `21` = "CD4 T cells",
#                     `22` = "CD4 T cells",
#                     `23` = "CD4 T cells",
#                     `24` = "CD4 T cells",
#                     `25` = "CD4 T cells",
#                     `26` = "CD4 T cells",
#                     `27` = "CD4 T cells",
#                     `28` = "CD4 T cells",
#                     `29` = "CD4 T cells",
#                     `30` = "CD4 T cells",
#                     `31` = "CD4 T cells",
#                     `32` = "CD4 T cells",
#                     `33` = "CD4 T cells",
#                     `34` = "CD4 T cells",
#                     `35` = "CD4 T cells",
#                     `36` = "CD4 T cells",
#                     `37` = "CD4 T cells",
#                     `38` = "CD4 T cells",
#                     `39` = "CD8 T cells",
#                     `40` = "CD4 T cells",
#                     `41` = "CD4 T cells",
#                     `42` = "CD4 T cells",
#                     `43` = "CD4 T cells",
#                     `44` = "CD4 T cells",
#                     `45` = "CD4 T cells",
#                     `46` = "CD4 T cells",
#                     `47` = "CD4 T cells",
#                     `48` = "CD4 T cells",
#                     `49` = "CD4 T cells",
#                     `50` = "NK cells",
#                     `51` = "CD4 T cells")

#5
# dat <- RenameIdents(dat,
#                     `1` = "CD8 T cells",
#                     `2` = "CD4 T cells",
#                     `3` = "CD4 T cells",
#                     `4` = "CD4 T cells",
#                     `5` = "CD4 T cells",
#                     `6` = "CD8 T cells",
#                     `7` = "CD4 T cells",
#                     `8` = "CD4 T cells",
#                     `9` = "CD4 T cells",
#                     `10` = "CD8 T cells",
#                     `11` = "CD4 T cells",
#                     `12` = "CD8 T cells",
#                     `13` = "CD8 T cells",
#                     `14` = "CD8 T cells",
#                     `15` = "CD4 T cells",
#                     `16` = "CD4 T cells",
#                     `17` = "CD4 T cells",
#                     `18` = "CD8 T cells",
#                     `19` = "CD4 T cells",
#                     `20` = "CD8 T cells",
#                     `21` = "CD8 T cells",
#                     `22` = "CD8 T cells",
#                     `23` = "CD8 T cells",
#                     `24` = "CD4 T cells",
#                     `25` = "CD4 T cells",
#                     `26` = "CD8 T cells",
#                     `27` = "CD4 T cells",
#                     `28` = "CD4 T cells",
#                     `29` = "CD4 T cells",
#                     `30` = "CD8 T cells",
#                     `31` = "CD4 T cells",
#                     `32` = "CD4 T cells",
#                     `33` = "CD4 T cells",
#                     `34` = "CD8 T cells",
#                     `35` = "CD4 T cells",
#                     `36` = "CD4 T cells",
#                     `37` = "CD8 T cells",
#                     `38` = "CD4 T cells",
#                     `39` = "CD8 T cells",
#                     `40` = "CD4 T cells",
#                     `41` = "NK cells",
#                     `42` = "CD8 T / NK cells")

#6
# dat <- RenameIdents(dat,
#                     `1` = "CD8 T cells",
#                     `2` = "CD8 T cells",
#                     `3` = "CD8 T cells",
#                     `4` = "CD4 T cells",
#                     `5` = "CD8 T cells",
#                     `6` = "CD8 T cells",
#                     `7` = "CD4 T cells",
#                     `8` = "CD4 T cells",
#                     `9` = "CD8 T cells",
#                     `10` = "CD4 T cells",
#                     `11` = "CD8 T cells",
#                     `12` = "CD8 T cells",
#                     `13` = "CD8 T cells",
#                     `14` = "NK cells")

#7
# dat <- RenameIdents(dat,
#                     `1` = "CD8 T cells",
#                     `2` = "CD4 T cells",
#                     `3` = "CD8 T cells",
#                     `4` = "NK cells")

#8
# dat <- RenameIdents(dat,
#                     `1` = "CD8 T cells",
#                     `2` = "CD8 T cells",
#                     `3` = "CD8 T cells",
#                     `4` = "CD8 T cells",
#                     `5` = "CD4 T cells",
#                     `6` = "CD4 T cells",
#                     `7` = "CD8 T cells",
#                     `8` = "NK cells",
#                     `9` = "CD8 T cells")

tils[[sname]] <- dat

saveRDS(tils,file=paste0(outdir,"tils.rda"))

for (sname in names(tils)){
  pdf(file = paste0(outdir,sname,".umap.annotated.pdf"), width = 10, height = 10)
  print(DimPlot(tils[[sname]],
          reduction = "umap", 
          label = T))
  dev.off()
}

# Integrate --------------------------------------------------------------------
library(harmony)

snames_m <- names(which(unlist(lapply(tils,function(x){unique(x@meta.data$Sex)})) == "M"))
snames_f <- names(which(unlist(lapply(tils,function(x){unique(x@meta.data$Sex)})) == "F"))

tils.merged <- merge(tils[[which(names(tils) %in% snames_m)[1]]], 
                     y = tils[setdiff(which(names(tils) %in% snames_m),which(names(tils) %in% snames_m)[1])], 
                     add.cell.ids = names(tils[which(names(tils) %in% snames_m)]), project = "tils")
tils.merged <- NormalizeData(tils.merged, verbose = T)
tils.merged <- CellCycleScoring(tils.merged, s.features = cc.genes$s.genes, 
                        g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
tils.merged <- ScaleData(tils.merged, verbose = T, vars.to.regress = c("S.Score", "G2M.Score"))
#tils.merged <- ScaleData(tils.merged, verbose = T)
tils.merged <- FindVariableFeatures(tils.merged, selection.method = "vst", nfeatures = 3000)
tils.merged <- RunPCA(tils.merged, pc.genes = VariableFeatures(tils.merged), npcs = 30, verbose = T)
#unique(tils.merged@meta.data$Sample.ID)

DimPlot(tils.merged,
        reduction = "pca", 
        label = T, dims = c(2,3))

DimPlot(tils.merged,
        reduction = "pca", 
        label = T, group.by = "old.ident")

DimPlot(tils.merged,
        reduction = "pca", group.by = "Sample.ID",
        label = T)

FeaturePlot(tils.merged,features = "HAVCR2",
        reduction = "pca", dims = c(1,2),
        label = T)

FeaturePlot(tils.merged,features = "CD8A",
            reduction = "pca", dims = c(1,2),
            label = T)

FeaturePlot(tils.merged,features = "CD4",
            reduction = "pca", dims = c(1,2),
            label = T)

#tils.merged.bak <- tils.merged

Idents(tils.merged) <- tils.merged$old.ident


tils_harmony <- RunHarmony(tils.merged, group.by.vars = "Sample.ID", 
                           plot_convergence = TRUE)

harmony_embeddings <- Embeddings(tils_harmony, 'harmony')

DimPlot(tils_harmony,
        reduction = "harmony", 
        label = T, dims = c(3,4))

# library(psych)
# color = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],
#                length(levels(Idents(tils_harmony))))
# 
# pairs.panels(harmony_embeddings[,1:30],
#              gap=0,
#              bg = color[Idents(tils_harmony)],
#              pch=21,
#              smooth = F,density = F,ellipses = F,hist.col = NA,rug = F)

tils_harmony <- RunUMAP(tils_harmony, reduction = "harmony", dims = 1:30, verbose = T)
tils_harmony <- RunTSNE(tils_harmony, reduction = "harmony", dims = 1:30, verbose = T)
tils_harmony <- FindNeighbors(tils_harmony, reduction = "harmony", dims = 1:30)
#tils_harmony <- FindClusters(tils_harmony, resolution = 1)

DimPlot(tils_harmony,
        reduction = "umap", 
        label = T)

DimPlot(tils_harmony,
        reduction = "tsne", 
        label = T)

# Try to integrate samples that are either CD4 or CD8 dominated separately -----

snames_cd8_m <- c("SampleID_9_11june18")
snames_cd4_m <- c("SampleID_1_11june18","SampleID_3_11june18")
snames_balanced_m <- c("SampleID_7_11june18","SampleID_8_11june18")

tils.merged.balanced <- merge(tils[[which(names(tils) %in% snames_balanced_m)[1]]], 
                     y = tils[setdiff(which(names(tils) %in% snames_balanced_m),which(names(tils) %in% snames_balanced_m)[1])], 
                     add.cell.ids = names(tils[which(names(tils) %in% snames_balanced_m)]), project = "tils")
tils.merged.balanced <- NormalizeData(tils.merged.balanced, verbose = T)
tils.merged.balanced <- CellCycleScoring(tils.merged.balanced, s.features = cc.genes$s.genes, 
                                g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
tils.merged.balanced <- ScaleData(tils.merged.balanced, verbose = T, vars.to.regress = c("S.Score", "G2M.Score"))
#tils.merged.balanced <- ScaleData(tils.merged.balanced, verbose = T)
tils.merged.balanced <- FindVariableFeatures(tils.merged.balanced, selection.method = "vst", nfeatures = 3000)
tils.merged.balanced <- RunPCA(tils.merged.balanced, pc.genes = VariableFeatures(tils.merged.balanced), npcs = 30, verbose = T)

tils.merged.balanced_harmony <- RunHarmony(tils.merged.balanced, group.by.vars = "Sample.ID", 
                           plot_convergence = TRUE, max.iter.harmony = 100)
tils.merged.balanced_harmony <- RunUMAP(tils.merged.balanced_harmony, reduction = "harmony", 
                                        dims = 1:30, verbose = T, min.dist=0.001, negative.sample.rate = 100L)
tils.merged.balanced_harmony <- RunTSNE(tils.merged.balanced_harmony, reduction = "harmony", 
                                        dims = 1:30, verbose = T)
tils.merged.balanced_harmony <- FindNeighbors(tils.merged.balanced_harmony, reduction = "harmony", dims = 1:30)

Idents(tils.merged.balanced_harmony) <- tils.merged.balanced_harmony$old.ident

DimPlot(tils.merged.balanced_harmony,
        reduction = "umap", 
        label = T)
DimPlot(tils.merged.balanced_harmony,
        reduction = "tsne", 
        label = T)

FeaturePlot(tils.merged.balanced_harmony, reduction = "umap", features = c("CD4","CD8A"), order = T)
FeaturePlot(tils.merged.balanced_harmony, reduction = "umap", features = c("HAVCR2","CCNB2"), order = T)
FeaturePlot(tils.merged.balanced_harmony, reduction = "umap", features = c("IL2RA","IL2RB"), order = T)

DotPlot(tils.merged.balanced_harmony, 
        features = c("CD3D","CD4","CD8A","CD8B","NCAM1","NCR1","FCGR3A"),
        cols = c("blue", "red"), dot.scale = 8,assay = "RNA",scale = F) +
  RotatedAxis()

idx <- tils.merged.balanced_harmony@assays$RNA@counts["CD4",] > 0 & (tils.merged.balanced_harmony@assays$RNA@counts["CD8A",] > 0 | tils.merged.balanced_harmony@assays$RNA@counts["CD8A",] > 0)

table(tils.merged.balanced_harmony@active.ident,idx)

Idents(tils.merged.balanced_harmony,cells = which(idx)) <- "T cells"

idx <- tils.merged.balanced_harmony@assays$RNA@counts["CD4",] > 0 & tils.merged.balanced_harmony@active.ident == "CD8 T cells"
Idents(tils.merged.balanced_harmony,cells = which(idx)) <- "CD4 T cells"

idx <- (tils.merged.balanced_harmony@assays$RNA@counts["CD8A",] > 0 | 
          tils.merged.balanced_harmony@assays$RNA@counts["CD8B",] > 0) & 
  tils.merged.balanced_harmony@active.ident == "CD4 T cells"

Idents(tils.merged.balanced_harmony,cells = which(idx)) <- "CD8 T cells"

# Attempt with LIGER -----------------------------------------------------------
library(SeuratWrappers)
library(rliger)

tils.merged.balanced.liger <- merge(tils[[which(names(tils) %in% snames_balanced_m)[1]]], 
                              y = tils[setdiff(which(names(tils) %in% snames_balanced_m),which(names(tils) %in% snames_balanced_m)[1])], 
                              add.cell.ids = names(tils[which(names(tils) %in% snames_balanced_m)]), project = "tils")
tils.merged.balanced.liger <- NormalizeData(tils.merged.balanced.liger, verbose = T)
tils.merged.balanced.liger <- FindVariableFeatures(tils.merged.balanced.liger, verbose = T)
tils.merged.balanced.liger <- ScaleData(tils.merged.balanced.liger, verbose = T, split.by = "Sample.ID", do.center = F)
tils.merged.balanced.liger <- RunOptimizeALS(tils.merged.balanced.liger, k = 20, lambda = 5, split.by = "Sample.ID")
tils.merged.balanced.liger <- RunQuantileNorm(tils.merged.balanced.liger, split.by = "Sample.ID")
tils.merged.balanced.liger <- FindNeighbors(tils.merged.balanced.liger, reduction = "iNMF", dims = 1:20)
# tils.merged.balanced.liger <- RunUMAP(tils.merged.balanced.liger, 
#                                       dims = 1:ncol(tils.merged.balanced.liger[["iNMF"]]), 
#                                       reduction = "iNMF",
                                  #    min.dist=0.001, negative.sample.rate = 100L)
tils.merged.balanced.liger <- RunUMAP(tils.merged.balanced.liger, 
                                      dims = 1:ncol(tils.merged.balanced.liger[["iNMF"]]), 
                                      reduction = "iNMF")
stopifnot(identical(rownames(tils.merged.balanced.liger@meta.data),rownames(tils.merged.balanced_harmony@meta.data)))
Idents(tils.merged.balanced.liger) <- tils.merged.balanced_harmony@active.ident

DimPlot(tils.merged.balanced.liger, group.by = c("Sample.ID", "ident"), 
        ncol = 2, reduction = "umap")


tils.merged.balanced.liger.cc <- merge(tils[[which(names(tils) %in% snames_balanced_m)[1]]], 
                                    y = tils[setdiff(which(names(tils) %in% snames_balanced_m),which(names(tils) %in% snames_balanced_m)[1])], 
                                    add.cell.ids = names(tils[which(names(tils) %in% snames_balanced_m)]), project = "tils")
tils.merged.balanced.liger.cc <- NormalizeData(tils.merged.balanced.liger.cc, verbose = T)
tils.merged.balanced.liger.cc <- FindVariableFeatures(tils.merged.balanced.liger.cc, verbose = T)
tils.merged.balanced.liger.cc <- CellCycleScoring(tils.merged.balanced.liger.cc, s.features = cc.genes$s.genes, 
                                         g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
tils.merged.balanced.liger.cc <- ScaleData(tils.merged.balanced.liger.cc, verbose = T, 
                                           split.by = "Sample.ID", do.center = F,
                                           vars.to.regress = c("S.Score", "G2M.Score"))
tils.merged.balanced.liger.cc <- RunOptimizeALS(tils.merged.balanced.liger.cc, k = 20, lambda = 5, split.by = "Sample.ID")
tils.merged.balanced.liger.cc <- RunQuantileNorm(tils.merged.balanced.liger.cc, split.by = "Sample.ID")
tils.merged.balanced.liger.cc <- FindNeighbors(tils.merged.balanced.liger.cc, reduction = "iNMF", dims = 1:20)
# tils.merged.balanced.liger.cc <- RunUMAP(tils.merged.balanced.liger.cc, 
#                                       dims = 1:ncol(tils.merged.balanced.liger.cc[["iNMF"]]), 
#                                       reduction = "iNMF",
#                                       min.dist=0.001, negative.sample.rate = 100L)
tils.merged.balanced.liger.cc <- RunUMAP(tils.merged.balanced.liger.cc, 
                                         dims = 1:ncol(tils.merged.balanced.liger.cc[["iNMF"]]), 
                                         reduction = "iNMF")

stopifnot(identical(rownames(tils.merged.balanced.liger.cc@meta.data),rownames(tils.merged.balanced_harmony@meta.data)))
Idents(tils.merged.balanced.liger.cc) <- tils.merged.balanced_harmony@active.ident

DimPlot(tils.merged.balanced.liger.cc, group.by = c("Sample.ID", "ident"), 
        ncol = 2, reduction = "umap")

# Try LIGER on all male --------------------------------------------------------

tils.merged.balanced.liger.m <- merge(tils[[which(names(tils) %in% snames_m)[1]]], 
                                    y = tils[setdiff(which(names(tils) %in% snames_m),which(names(tils) %in% snames_m)[1])], 
                                    add.cell.ids = names(tils[which(names(tils) %in% snames_m)]), project = "tils")
tils.merged.balanced.liger.m <- NormalizeData(tils.merged.balanced.liger.m, verbose = T)
tils.merged.balanced.liger.m <- FindVariableFeatures(tils.merged.balanced.liger.m, verbose = T)
# tils.merged.balanced.liger.m <- CellCycleScoring(tils.merged.balanced.liger.m, s.features = cc.genes$s.genes,
#                                                   g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
# tils.merged.balanced.liger.m <- ScaleData(tils.merged.balanced.liger.m, verbose = T,
#                                            split.by = "Sample.ID", do.center = F,
#                                            vars.to.regress = c("S.Score", "G2M.Score"))
tils.merged.balanced.liger.m <- ScaleData(tils.merged.balanced.liger.m, verbose = T, split.by = "Sample.ID", do.center = F)
tils.merged.balanced.liger.m <- RunOptimizeALS(tils.merged.balanced.liger.m, k = 20, lambda = 5, split.by = "Sample.ID")
tils.merged.balanced.liger.m <- RunQuantileNorm(tils.merged.balanced.liger.m, split.by = "Sample.ID")
tils.merged.balanced.liger.m <- FindNeighbors(tils.merged.balanced.liger.m, reduction = "iNMF", dims = 1:20)
# tils.merged.balanced.liger.m <- RunUMAP(tils.merged.balanced.liger.m, 
#                                       dims = 1:ncol(tils.merged.balanced.liger.m[["iNMF"]]), 
#                                       reduction = "iNMF",
#    min.dist=0.001, negative.sample.rate = 100L)
tils.merged.balanced.liger.m <- RunUMAP(tils.merged.balanced.liger.m, 
                                      dims = 1:ncol(tils.merged.balanced.liger.m[["iNMF"]]), 
                                      reduction = "iNMF")

tils.merged.balanced.liger.m <- AddMetaData(tils.merged.balanced.liger.m, metadata = tils.merged$old.ident,col.name = "celltype")
Idents(tils.merged.balanced.liger.m) <- tils.merged.balanced.liger.m@meta.data$celltype

DimPlot(tils.merged.balanced.liger.m, group.by = c("Sample.ID", "ident"), 
        ncol = 2, reduction = "umap")

# Try LIGER on all--------------------------------------------------------------
tils <- readRDS(file=paste0(outdir,"tils.rda"))

snames_all <- names(tils)

tils.allerged.balanced.liger.all <- merge(tils[[which(names(tils) %in% snames_all)[1]]], 
                                      y = tils[setdiff(which(names(tils) %in% snames_all),
                                                       which(names(tils) %in% snames_all)[1])], 
                                      add.cell.ids = names(tils[which(names(tils) %in% snames_all)]), 
                                      project = "tils")
tils.allerged.balanced.liger.all <- AddMetaData(tils.allerged.balanced.liger.all,
            metadata = tils.allerged.balanced.liger.all@active.ident,
            col.name = "old.ident")
tils.allerged.balanced.liger.all <- NormalizeData(tils.allerged.balanced.liger.all, 
                                                  normalization.method = "RC", verbose = T)
#tils.allerged.balanced.liger.all <- FindVariableFeatures(tils.allerged.balanced.liger.all, 
#                                                         verbose = T,nfeatures = 10000)

variable_features <- list()
for (s in names(tils)){
  tmp <- subset(tils.allerged.balanced.liger.all,Sample.ID == s)
  tmp <- FindVariableFeatures(tmp, verbose = T,nfeatures = 10000)
  variable_features[[s]] <- VariableFeatures(tmp)
}

vf_intersect <- Reduce("intersect",variable_features)
vf_t <- sort(table(unlist(variable_features)),decreasing = T)

#x <- cor(t(as.matrix(tils.allerged.balanced.liger.all@assays$RNA@data[vf_intersect,])))
x <- cor(t(as.matrix(tils.allerged.balanced.liger.all@assays$RNA@data[names(vf_t[vf_t > 4]),])))
p <- pheatmap(x, show_rownames = F, show_colnames = F, cutree_rows = 50, silent = T)
gene_clusters <- cutree(p$tree_row,50)
table(gene_clusters)
gene_clusters.subsample <- list()
for (i in unique(gene_clusters)){
  gene_clusters.subsample[[i]] <- names(sample(gene_clusters[gene_clusters==i],min(table(gene_clusters))))
}
gene_clusters.subsample <- unlist(gene_clusters.subsample)

y <- tils.allerged.balanced.liger.all@assays$RNA@data[names(which(gene_clusters == 2)),]
#y <- tils.allerged.balanced.liger.all@assays$RNA@data[gene_clusters.subsample,]
y <- aggregate(t(y),list(as.character(tils.allerged.balanced.liger.all@active.ident)),mean)
rownames(y) <- y$Group.1
y$Group.1 <- NULL
pheatmap(log2(y+1),show_colnames = F,border_color = NA,scale = "column")
#pheatmap(log2(y+1),show_colnames = F,border_color = NA)

#y <- tils.allerged.balanced.liger.all@assays$RNA@data[gene_clusters.subsample,]
y <- tils.allerged.balanced.liger.all@assays$RNA@data[names(which(gene_clusters == 1)),]
y <- aggregate(t(y),list(as.character(tils.allerged.balanced.liger.all@meta.data$Phase)),mean)
rownames(y) <- y$Group.1
y$Group.1 <- NULL
pheatmap(log2(y+1),show_colnames = F,border_color = NA,scale = "column")

cr <- list()
#for (i in unique(gene_clusters[gene_clusters %in% names(which(table(gene_clusters) > 1))])){
clusters_rm <- c()
#for (i in unique(gene_clusters[! gene_clusters %in% c(clusters_rm,2) & gene_clusters %in% names(which(table(gene_clusters) > 1))])){
for (i in unique(gene_clusters[ gene_clusters %in% names(which(table(gene_clusters) > 1))])){
  y <- tils.allerged.balanced.liger.all@assays$RNA@data[names(which(gene_clusters == i)),]
  cr[[as.character(i)]] <- data.frame(G2M_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$G2M.Score),
             S_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$S.Score),
             Sex_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sex=="M"),
             CD4_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@active.ident=="CD4 T cells"),
             CD8_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@active.ident=="CD8 T cells"),
             NK_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@active.ident=="NK cells"),
             CD8_NK_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@active.ident=="CD8 T / NK cells"),
             GDT_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@active.ident=="Gamma-delta T cells (TRGV9, TRDV2)"),
             S1_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sample.ID=="SampleID_1_11june18"),
             S2_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sample.ID=="SampleID_2_11june18"),
             S3_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sample.ID=="SampleID_3_11june18"),
             S4_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sample.ID=="SampleID_4_11june18"),
             S6_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sample.ID=="SampleID_6_11june18"),
             S7_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sample.ID=="SampleID_7_11june18"),
             S8_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sample.ID=="SampleID_8_11june18"),
             S9_cor=cor(colMeans(y),tils.allerged.balanced.liger.all@meta.data$Sample.ID=="SampleID_9_11june18"))
}
cr <- do.call("rbind",cr)
pheatmap(cr,cluster_rows = F,cluster_cols = F,border_color = NA,display_numbers = T)

clusters_rm <- c(1,3,4,6,7,8,9,12,14,15,16,17,20,21,22,23,25,26,28,29,36,37,40,41,42,43,50)
gene_clusters.subsample <- list()
for (i in setdiff(unique(gene_clusters),clusters_rm)){
  gene_clusters.subsample[[i]] <- names(sample(gene_clusters[gene_clusters==i],
                                               min(table(gene_clusters[! gene_clusters %in% clusters_rm]))))
}
gene_clusters.subsample <- unlist(gene_clusters.subsample)

sz <- mean(table(gene_clusters[! gene_clusters %in% c(clusters_rm,2)]))
sample(gene_clusters[gene_clusters==2],sz,replace = F)

gene_clusters.selected <- union(names(gene_clusters[! gene_clusters %in% c(clusters_rm,2)]),
                                names(sample(gene_clusters[gene_clusters==2],sz,replace = F)))

# tils.allerged.balanced.liger.all <- CellCycleScoring(tils.allerged.balanced.liger.all, s.features = cc.genes$s.genes,
#                                                   g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
# tils.allerged.balanced.liger.all <- ScaleData(tils.allerged.balanced.liger.all, verbose = T,
#                                            split.by = "Sample.ID", do.center = F,
#                                            vars.to.regress = c("S.Score", "G2M.Score"))
#VariableFeatures(tils.allerged.balanced.liger.all) <- names(vf_t[vf_t > 4])
#VariableFeatures(tils.allerged.balanced.liger.all) <- names(gene_clusters[! gene_clusters %in% clusters_rm])
VariableFeatures(tils.allerged.balanced.liger.all) <- gene_clusters.selected
#y <- aggregate(t(x),list(as.character(tils.allerged.balanced.liger.all@active.ident)),mean)
# rownames(y) <- y$Group.1
# y$Group.1 <- NULL
# z <- y[,colMeans(y) >= summary(colMeans(y))[["Median"]]]
# z <- z[rownames(z) %in% c("CD8 T cells","CD4 T cells"),]
# cd4_cd8_genes <- names(which(apply(z,2,sd)/apply(z,2,mean) > 0.8))
# VariableFeatures(tils.allerged.balanced.liger.all) <- cd4_cd8_genes

tils.allerged.balanced.liger.all <- ScaleData(tils.allerged.balanced.liger.all, verbose = T, 
                                              split.by = "Sample.ID", do.center = F,
                                              features = VariableFeatures(tils.allerged.balanced.liger.all))
tils.allerged.balanced.liger.all <- RunOptimizeALS(tils.allerged.balanced.liger.all, 
                                                   k = 20, lambda = 1, split.by = "Sample.ID",
                                                   nrep = 1)
#dims_use <- c(1,14,37,13,40,4,16,32,27,2,7,22)
#dims_use <- setdiff(1:20,c(11,10,17,9,7,18,8,16))
#dims_use <- setdiff(1:20,c(12,8,14,6,15,4,7,20,13,19,1,9,10))
#dims_use <- setdiff(1:40,c(3,19,11,28,7,18,26,1,23,25,5,8,22,4,10,12,20,31))
dims_use <- setdiff(1:20,components_exclude)
tils.allerged.balanced.liger.all <- RunQuantileNorm(tils.allerged.balanced.liger.all, 
                                                    split.by = "Sample.ID",
                                                    eps = 0.01, 
                                                    ref_dataset = "SampleID_6_11june18",
                                                    do.center = T,
                                                    dims.use = dims_use,
                                                    knn_k = 2500,
                                                    quantiles = 100,
                                                    min_cells = 3000, max_sample = 10000)
# knn_k: "Number of nearest neighbors for within-dataset knn graph"

tils.allerged.balanced.liger.all <- FindNeighbors(tils.allerged.balanced.liger.all, 
                                                  reduction = "iNMF", 
                                                  dims = dims_use,k.param = 200)
tils.allerged.balanced.liger.all <- RunUMAP(tils.allerged.balanced.liger.all, 
                                        dims = dims_use, 
                                        reduction = "iNMF",n.neighbors = 200)

Idents(tils.allerged.balanced.liger.all) <- tils.allerged.balanced.liger.all@meta.data$old.ident

DimPlot(tils.allerged.balanced.liger.all, group.by = c("Sample.ID", "ident"),
        ncol = 2, reduction = "umap")

FeaturePlot(tils.allerged.balanced.liger.all,features=c("CD8A","CD4"),order = T,slot = "counts")



lig <- seuratToLiger(objects = tils,
                     combined.seurat = F,
                     names = "use-meta",
                     meta.var = "Sample.ID",
                     renormalize=F,
                     use.idents=T,
                     use.seurat.genes=F,
                     remove.missing=F)
lig <- normalize(lig)
lig <- selectGenes(lig, combine = "union", do.plot = T)
lig <- scaleNotCenter(lig)
lig <- optimizeALS(lig, k = 50, lambda = 5,verbose = T)

dims.use.vetted <- c(45, 50, 25, 14, 37, 35, 29)
dims.remove <- c(5,6,8,10,11,12,15,17,20,24,27,28,31,32,34,38,39,40,41,42,48,21,22,49,19,2,4)
dims.remove.vetted <- setdiff(c(45, 50, 3, 47, 25, 33, 16, 14, 37, 1, 7, 43, 9, 35, 13, 29, 23),dims.use.vetted)
dims.not_yet_vetted <- setdiff(1:50,union(dims.use.vetted,dims.remove.vetted))
dims.use <- setdiff(union(dims.use.vetted,setdiff(dims.not_yet_vetted,dims.remove)),c(18,29,36,30,35)) #44
dims.use <- c(dims.use,19,49,34)

lig <- quantile_norm(lig, quantiles = 100, ref_dataset = "SampleID_6_11june18", 
                     eps = 0.01, min_cells = 200, knn_k = 2500, max_sample = 2000,
                     do.center = T,refine.knn = T,
                     dims.use = dims.use) #knn_k = 2500,min_cells = 3000, max_sample = 10000, 
lig <- runUMAP(lig, distance = 'euclidean',dims.use = dims.use, 
               min_dist = 0.0000000001, n_neighbors = 10)#, n_neighbors = 20, min_dist = 0.000001,

stopifnot(identical(names(lig@clusters),
                    gsub("SampleID_[0-9]_11june18_","",rownames(tils.allerged.balanced.liger.all@meta.data))))
celltypes <- setNames(tils.allerged.balanced.liger.all@meta.data$old.ident,names(lig@clusters))

all.plots <- plotByDatasetAndCluster(lig, axis.labels = c('UMAP 1', 'UMAP 2'), 
                                     return.plots = T, clusters = celltypes)
all.plots[[2]]

all.plots[[1]] + all.plots[[2]]

gene_loadings <- plotGeneLoadings(lig, do.spec.plot = T, return.plots = TRUE)

dims.use

i <- 34
all.plots[[2]] + gene_loadings[[i]]
all.plots[[1]] + gene_loadings[[i]]
gene_loadings[[i]]
head(unlist(gsea_output[[i]][order(unlist(gsea_output[[i]][,5]),decreasing = T),c(1)]),n=10)

# 34: Possible Treg signature, although also related to apoptosis
# 49: TCF7, IL7R
# 19: LAG3
# 35: SampleID_7-specific junk

plotGenes(lig,genes=c("CD8A","CD8B","CD4"),scale.by = "none", plot.by = "none",
          keep.scale = T, set.dr.lims = T)


gout <- lapply(gsea_output,function(x){as.data.frame(x[,c(1,5)])})
for (i in 1:length(gout)){gout[[i]]$nr <- i}
gout <- do.call("rbind",gout)
gout <- spread(gout, key = nr, value = NES)
rownames(gout) <- gout$pathway
gout$pathway <- NULL
tmp <- sapply(gout,as.numeric)
rownames(tmp) <- rownames(gout)
gout <- tmp
rm(tmp)

gout[is.na(gout)] <- 0

pdf(paste0(outdir,"liger_gsea_heatmap.pdf"),width = 11,height = 21)
pheatmap(gout,fontsize_row = 5)
dev.off()
pdf(paste0(outdir,"liger_gsea_heatmap.cor_1.pdf"),width = 10,height = 10)
pheatmap(cor(gout),border_color = NA)
dev.off()
pdf(paste0(outdir,"liger_gsea_heatmap.cor_2.pdf"),width = 25,height = 21)
p <- pheatmap(cor(t(gout)),fontsize_row = 5,show_colnames = F,
              cutree_rows = 10,cutree_cols = 10)
dev.off()

clus <- cutree(p$tree_row,10)
clus[clus==1] <- "Stress_translation"
clus[clus==2] <- "Unknown_1 (Notch 1 signaling)"
clus[clus==3] <- "Unknown_2 (Interleukin signaling)"
clus[clus==4] <- "Cell cycle"
clus[clus==5] <- "MAPK_RAS_RAF_signaling"
clus[clus==6] <- "Unknown_3 (PD-1 signaling)"
clus[clus==7] <- "Unknown_4 (Interferon signaling)"
clus[clus==8] <- "Toll-like receptors, senescence, DNA repair"
clus[clus==9] <- "Unknown_5 (Chromatin organization, glycosylation)"
clus[clus==10] <- "Unknown_6 (G-protein signaling)"

cluster_means <- list()
for (i in unique(clus)){
  cluster_means[[i]] <- colMeans(gout[names(clus)[clus==i],])
}

cluster_means <- do.call("rbind",cluster_means)
pheatmap(cluster_means,border_color = NA)
pdf(paste0(outdir,"liger_gsea_heatmap.module_means.pdf"),width = 11,height = 2.5)
pheatmap(t(scale(t(cluster_means),center = T,scale = T)),
         border_color = NA,cutree_cols = 10)
dev.off()
