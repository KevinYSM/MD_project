# Export UMAP coordinates for Loupe --------------------------------------------

sample_id_map <- as.data.frame(read.table("~/Desktop/SampleID.csv",sep=",",header=T),
                               stringsAsFactors=F)



get_tsne_loupe <- function(s.t,sample_id_map){
  df <- as.data.frame(s.t@reductions$tsne@cell.embeddings)
  df$Barcode <- rownames(df)
  df$Barcode <- gsub("[A-Z][0-9]-GEX_","",df$Barcode)
  df$SampleID <- rownames(df)
  df$SampleID <- gsub("_[0-9]_[A-Z]+_[0-9]$","",df$SampleID)
  
  df$SampleID <- paste0("GWA-JN-388_",df$SampleID)
  
  tab <- table(str_split_fixed(sample_id_map$Barcode,"-",2)[,2],sample_id_map$SampleID)
  tab <- as.data.frame.matrix(tab)
  tab[tab>0] <- 1
  tab <- tab*as.numeric(rownames(tab))
  
  tab <- colSums(tab)
  sample_id_map2 <- as.data.frame(tab)
  colnames(sample_id_map2) <- "sample_nr"
  sample_id_map2$sample_id <- rownames(sample_id_map2)
  rownames(sample_id_map2) <- NULL
  
  df$Barcode2 <- str_split_fixed(df$Barcode,"_",3)[,2]
  df2 <- merge(df,sample_id_map2,by.x="SampleID",by.y="sample_id")
  df2$Barcode2 <- paste0(df2$Barcode2,"-",df2$sample_nr)
  df2 <- df2[,c("Barcode2","tSNE_1","tSNE_2")]
  colnames(df2)[colnames(df2)=="Barcode2"] <- "Barcode"
  colnames(df2)[colnames(df2)=="tSNE_1"] <- "X Coordinate"
  colnames(df2)[colnames(df2)=="tSNE_2"] <- "Y Coordinate"
  return(df2)
}

tsne_loupe.s.t <- get_tsne_loupe(s.t,sample_id_map)
tsne_loupe.s <- get_tsne_loupe(s,sample_id_map)

write.table(tsne_loupe.s.t,file = "~/proj/um_ss/Manuscript/loupe_umap_cd8.csv",
            sep=",",col.names = T,row.names = F,quote = F)
write.table(tsne_loupe.s,file = "~/proj/um_ss/Manuscript/loupe_umap_all.csv",
            sep=",",col.names = T,row.names = F,quote = F)

get_categories_loupe <- function(s.t,sample_id_map){
  df <- data.frame(Barcode=names(s.t@active.ident),Category=as.character(s.t@active.ident))
  
  df$SampleID <- df$Barcode
  df$Barcode <- gsub("[A-Z][0-9]-GEX_","",df$Barcode)
  df$SampleID <- gsub("_[0-9]_[A-Z]+_[0-9]$","",df$SampleID)
  df$SampleID <- paste0("GWA-JN-388_",df$SampleID)
  
  tab <- table(str_split_fixed(sample_id_map$Barcode,"-",2)[,2],sample_id_map$SampleID)
  tab <- as.data.frame.matrix(tab)
  tab[tab>0] <- 1
  tab <- tab*as.numeric(rownames(tab))
  
  tab <- colSums(tab)
  sample_id_map2 <- as.data.frame(tab)
  colnames(sample_id_map2) <- "sample_nr"
  sample_id_map2$sample_id <- rownames(sample_id_map2)
  rownames(sample_id_map2) <- NULL
  
  df$Barcode2 <- str_split_fixed(df$Barcode,"_",3)[,2]
  df2 <- merge(df,sample_id_map2,by.x="SampleID",by.y="sample_id")
  df2$Barcode2 <- paste0(df2$Barcode2,"-",df2$sample_nr)
  df2 <- df2[,c("Barcode2","Category")]
  colnames(df2)[colnames(df2)=="Barcode2"] <- "Barcode"
  df2$Category <- gsub(",",";",df2$Category)
  return(df2)
}

categories_loupe.s.t <- get_categories_loupe(s.t,sample_id_map)
categories_loupe.s <- get_categories_loupe(s,sample_id_map)

write.table(categories_loupe.s.t,file = "~/proj/um_ss/Manuscript/loupe_categories_cd8.csv",
            sep=",",col.names = T,row.names = F,quote = F)
write.table(categories_loupe.s,file = "~/proj/um_ss/Manuscript/loupe_categories_all.csv",
            sep=",",col.names = T,row.names = F,quote = F)
