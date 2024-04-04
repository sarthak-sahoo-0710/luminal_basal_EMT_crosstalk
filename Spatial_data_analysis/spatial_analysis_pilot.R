library(Rmagic)
library(Seurat)
library(AUCell)
library(Biobase)
library(GSEABase)
library(GEOquery)

file_tisspos <- "./spatial/1160920F_spatial"
file_features <- "./filtered_count_matrices//1160920F_filtered_count_matrix"

X <- Read10X_Image(
  file_tisspos,
  image.name = "tissue_lowres_image.png",
  filter.matrix = TRUE
)

#  gene.column = 2 by default and that was causing an error, update to 1 for last version with gziped
my_object <- CreateSeuratObject(
  counts = Read10X( data.dir = file_features ,  gene.column = 1),
  assay = 'Spatial'
)

my_image <- Read10X_Image( image.dir = file_tisspos)

image <- my_image[Cells(x = my_object)]

DefaultAssay(object = image) <- 'Spatial'

my_object[['Slice1']] <- image

# running magic on the data

exprMatrix_T<- t(my_object@assays$Spatial@counts)
exprMatrix_T<- library.size.normalize(exprMatrix_T)
exprMatrix_T <- sqrt(exprMatrix_T)
data_MAGIC <- magic(exprMatrix_T, genes="all_genes", solver='approximate')
expr_Magic <- t(data_MAGIC$result)
rm(exprMatrix_T,data_MAGIC)

geneSets <- getGmt("pathways_HT.gmt")
geneSets <- subsetGeneSets(geneSets, rownames(expr_Magic)) 
cbind(nGenes(geneSets))

cells_rankings <- AUCell_buildRankings(expr_Magic, nCores=1, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
AUCmatrix <-  as.data.frame(cells_AUC@assays@data$AUC)
my_object <- AddMetaData(my_object, t(AUCmatrix), col.name = rownames(AUCmatrix))
rm(cells_rankings,cells_AUC,AUCmatrix)

my_object <- FindVariableFeatures(my_object, selection.method = "vst", nfeatures = 2000)
my_object <- ScaleData(my_object)
my_object <- RunPCA(my_object, assay = "Spatial", verbose = FALSE, features = rownames(my_object))
my_object <- FindNeighbors(my_object, reduction = "pca", dims = 1:30)
my_object <- FindClusters(my_object, verbose = FALSE, resolution = 0.5)
my_object <- RunUMAP(my_object, reduction = "pca", dims = 1:30)

#write.table(my_object@meta.data,"1142243F_metadata.txt",sep="\t",row.names = TRUE, quote = FALSE)
#write.table(as.matrix(my_object@reductions$umap@cell.embeddings),"1142243F_UMAP_coords.txt",sep="\t",row.names = TRUE, quote = FALSE)

SpatialFeaturePlot(my_object, features = c("Epithelial_tumour_signature"))
SpatialFeaturePlot(my_object, features = c("Mesenchymal_tumour_signature"))
SpatialFeaturePlot(my_object, features = c("pEMT"))

SpatialFeaturePlot(my_object, features = c("luminal_breast_cancer"))
SpatialFeaturePlot(my_object, features = c("basal_breast_cancer"))

SpatialFeaturePlot(my_object, features = c("HALLMARK_E2F_TARGETS"))
SpatialFeaturePlot(my_object, features = c("HALLMARK_ESTROGEN_RESPONSE_EARLY"))
SpatialFeaturePlot(my_object, features = c("HALLMARK_ESTROGEN_RESPONSE_LATE"))

SpatialFeaturePlot(my_object, features = c("IL1R1"))
SpatialFeaturePlot(my_object, features = c("GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN"))
SpatialFeaturePlot(my_object, features = c("GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN"))
SpatialFeaturePlot(my_object, features = c("basal_breast_cancer"))
SpatialFeaturePlot(my_object, features = c("top_il1r1_correlated_genes_HT1080"))



my_object <- SCTransform(my_object, assay = "Spatial", verbose = FALSE)

my_object <- FindVariableFeatures(my_object, selection.method = "vst", nfeatures = 2000)
my_object <- ScaleData(my_object)

my_object <- RunPCA(my_object, assay = "Spatial", verbose = FALSE, features = rownames(my_object))
my_object <- FindNeighbors(my_object, reduction = "pca", dims = 1:30)
my_object <- FindClusters(my_object, verbose = FALSE, resolution = 0.5)
my_object <- RunUMAP(my_object, reduction = "pca", dims = 1:30)

p1 <- DimPlot(my_object, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(my_object, label = TRUE, label.size = 3,pt.size.factor = 2,alpha = c(0.1, 1))
p1 + p2

SpatialDimPlot(my_object, cells.highlight = CellsByIdentities(object = my_object, idents = c(1)), facet.highlight = TRUE, ncol = 3)
