library(AUCell)
library(Biobase)
library(GSEABase)
library(data.table)
library(plotly)
library(GEOquery)
library(Seurat)
library(Rmagic)
library(readr)
library(ggplot2)
library(viridis)
library(phateR)
library(org.Hs.eg.db)
library("biomaRt")

setwd("E:/projects/breast_luminal_basal/single cell/sarthak/celllines_atlas/")

cellline_data <- Read10X(data.dir = "./data/",gene.column = 1)

celllines <- CreateSeuratObject(counts = cellline_data, project = "breast_celllines", min.cells = 5, min.features = 500)

# code to run umap analysis
celllines <- NormalizeData(celllines)
celllines <- FindVariableFeatures(celllines, selection.method = "vst", nfeatures = 2000)

celllines <- ScaleData(celllines, verbose = FALSE)
celllines <- RunPCA(celllines, npcs = 30, verbose = FALSE)
celllines <- RunUMAP(celllines, reduction = "pca", dims = 1:30)
celllines <- FindNeighbors(celllines, reduction = "pca", dims = 1:30)
celllines <- FindClusters(celllines, resolution = 0.5)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=celllines@assays$RNA@counts@Dimnames[[1]],
  uniqueRows=TRUE)

annotLookup <- data.frame(
  celllines@assays$RNA@counts@Dimnames[[1]][match(annotLookup$ensembl_gene_id, celllines@assays$RNA@counts@Dimnames[[1]])],
  annotLookup)

colnames(annotLookup) <- c(
  "original_id",
  c("ensembl_gene_id", "gene_biotype", "external_gene_name"))

p2 <- DimPlot(celllines, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(celllines, features = c("ENSG00000091831", "ENSG00000082175", "ENSG00000160182"))

# pathway score calculation

x <- celllines@assays$RNA@counts

geneSets <- getGmt("signatures_geneids.txt")
geneSets <- subsetGeneSets(geneSets, rownames(x)) 
cbind(nGenes(geneSets))

cells_rankings <- AUCell_buildRankings(x, nCores=1, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

AUCmatrix <-  as.data.frame(cells_AUC@assays@data$AUC)
write.table(AUCmatrix,"AUCell_scores_3.txt",sep="\t",row.names = TRUE, quote = FALSE)

rm(AUCmatrix, cellline_data, celllines, cells_AUC, cells_AUC, cells_rankings)

# magic for specific genes

exprMatrix_T<- t(celllines@assays$RNA@counts)
exprMatrix_T<- library.size.normalize(exprMatrix_T)
exprMatrix_T <- sqrt(exprMatrix_T)

data_MAGIC <- magic(exprMatrix_T, genes=c("ENSG00000115594","ENSG00000164362","ENSG00000270141","ENSG00000115008"), solver='approximate')
genes_of_interest = read.delim(file.choose(), header = FALSE)
selected_imputed_data <- data_MAGIC$result[c("ENSG00000115594","ENSG00000164362")]
rm(exprMatrix_T)
expr_Magic <- t(data_MAGIC$result)
rm(data_MAGIC)
write.table(round(selected_imputed_data[,],digits = 4),"selected_imputed_genes_shantanu.txt",sep="\t",row.names = TRUE, quote = FALSE)
