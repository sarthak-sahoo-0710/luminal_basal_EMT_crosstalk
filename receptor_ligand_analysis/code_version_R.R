library(magrittr)
library(liana)
library(Seurat)
library(AUCell)
library(Biobase)
library(GSEABase)
library(data.table)
library(plotly)
library(GEOquery)
library(Rmagic)
library(readr)
library(ggplot2)
library(viridis)
library(nichenetr)
library(phateR)
library(org.Hs.eg.db)
library("biomaRt")

clist = c('CID44971','CID44991','CID4513','CID4515','CID4523','CID4530N','CID4535')#'CID3941','CID3946','CID3948','CID3963','CID4040','CID4067','CID4290A','CID4398','CID44041','CID4461','CID4463','CID4465','CID4471','CID4495','CID44971','CID44991','CID4513','CID4515','CID4523','CID4530N','CID4535')
for (cn in clist) {
  data <- Read10X(data.dir = paste("./patient_data/", cn ,"/", sep = ""),gene.column=1)
  metadata <- read.csv(paste("./patient_data/", cn ,"/metadata.csv", sep = ""))
  rownames(metadata) = metadata$X
  seurat_object <- CreateSeuratObject(counts = data, meta.data = metadata)
  seurat_object <- NormalizeData(seurat_object)
  liana_test <- liana_wrap(seurat_object,idents_col = "celltype_minor")
  liana_test <- liana_test %>% liana_aggregate()
  write.table(liana_test,paste("./patient_data/", cn ,"/output_data_liana.tsv", sep = ""),sep='\t',quote = FALSE)
  
  x <- seurat_object@assays$RNA@counts
  geneSets <- getGmt("breast_cancer_specific_gene_lists.gmt")
  geneSets <- subsetGeneSets(geneSets, rownames(x)) 
  cbind(nGenes(geneSets))
  
  cells_rankings <- AUCell_buildRankings(x, plotStats=TRUE)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  AUCmatrix <-  as.data.frame(cells_AUC@assays@data$AUC)
  write.table(AUCmatrix,paste("./patient_data/", cn ,"/output_AUCellscores_not_imputed.tsv",sep=""),sep="\t",row.names = TRUE, quote = FALSE)
  
  exprMatrix_T<- t(seurat_object@assays$RNA@counts)
  exprMatrix_T<- library.size.normalize(exprMatrix_T)
  exprMatrix_T <- sqrt(exprMatrix_T)
  data_MAGIC <- magic(exprMatrix_T, solver='approximate')
  selected_imputed_data <- data_MAGIC$result
  write.table(round(selected_imputed_data,digits = 4),paste("./patient_data/", cn ,"/output_imputed_values.tsv",sep=""),sep="\t",row.names = TRUE, quote = FALSE)
  
  geneSets <- getGmt("breast_cancer_specific_gene_lists.gmt")
  geneSets <- subsetGeneSets(geneSets, rownames(t(selected_imputed_data))) 
  cbind(nGenes(geneSets))
  
  cells_rankings <- AUCell_buildRankings(t(selected_imputed_data), plotStats=TRUE)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  AUCmatrix <-  as.data.frame(cells_AUC@assays@data$AUC)
  write.table(AUCmatrix,paste("./patient_data/", cn ,"/output_AUCellscores_imputed.tsv",sep=""),sep="\t",row.names = TRUE, quote = FALSE)
}


            
