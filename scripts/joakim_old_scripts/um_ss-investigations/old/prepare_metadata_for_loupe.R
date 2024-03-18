# For Loupe --------------------------------------------------------------------
dat.cca.tmp.sct$active.ident.updated.simplified2 <- dat.cca.tmp.sct$active.ident.updated.simplified
dat.cca.tmp.sct$active.ident.updated.simplified2[grepl("CD8 T cells",dat.cca.tmp.sct$active.ident.updated.simplified2)] <- "CD8 T cells"
dat.cca.tmp.sct$active.ident.updated.simplified2[grepl("CD4 T cells",dat.cca.tmp.sct$active.ident.updated.simplified2)] <- "CD4 T cells"
dat.cca.tmp.sct$active.ident.updated.simplified2[dat.cca.tmp.sct$active.ident.updated.simplified2 == "CD4+ CD8+ T cells (unknown)"] <- "CD4+ CD8+ T cells (doublets?)"

categories <- as.data.frame(dat.cca.tmp.sct$active.ident.updated.simplified2)
categories <- data.frame(barcode = rownames(categories), id=categories[,1])
head(categories)
new_names <- names(table(apply(str_split_fixed(categories$barcode,"_",3)[,c(1,3)],1,function(x){paste0(x,collapse="-")})))
new_names <- new_names[order(str_split_fixed(new_names,"-",2)[,2],decreasing = F)]
new_names <- as.data.frame(cbind(new_names,1:length(new_names)))
new_names$snum <- str_split_fixed(new_names[,1],"-",2)[,1]
new_names$pnum <- str_split_fixed(new_names[,1],"-",2)[,2]
colnames(new_names)[1] <- "pair"
colnames(new_names)[2] <- "new_name"

categories$pair <- apply(str_split_fixed(categories$barcode,"_",3)[,c(1,3)],1,function(x){paste0(x,collapse="-")})

categories <- merge(categories,new_names,by="pair",all=T)
head(categories)
categories$barcode_new <- paste0(str_split_fixed(categories$barcode,"_",3)[,2],"-",categories$new_name)
head(categories)

all(categories$barcode %in% rownames(dat.cca.tmp.sct@meta.data))
all(rownames(dat.cca.tmp.sct@meta.data) %in% categories$barcode)
rownames(categories) <- categories$barcode
categories <- categories[rownames(dat.cca.tmp.sct@meta.data),]
identical(categories$barcode, rownames(dat.cca.tmp.sct@meta.data))

map.seurat <- table(as.numeric(categories$new_name),dat.cca.tmp.sct$Sample.ID)

loupe <- read.csv("~/Desktop/SampleID.csv")
loupe$new_name <- str_split_fixed(loupe$Barcode,"-",2)[,2]
map.loupe <- table(as.numeric(loupe$new_name),loupe$SampleID)

colnames(map.loupe) <- str_split_fixed(colnames(map.loupe),"_",2)[,2]
colnames(map.seurat)

map.loupe <- map.loupe[,colnames(map.seurat)]
map.loupe[map.loupe>0] <- 1

map.seurat[map.seurat>0] <- 1

identical(map.seurat,map.loupe)

categories <- categories[,c("barcode_new","id")]
write.table(categories,file="~/proj/um_ss/Pipelines/10x/results/gex_aggr/seurat_categories.csv",sep = ",",quote = F,col.names = T,row.names = F)