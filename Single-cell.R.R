library(Seurat)
library(dplyr)
library(ggplot2)
library(celldex)
library(SingleR)


pbmc.data2 <- Read10X(data.dir = "SC2/")
pbmc.data3 <- Read10X(data.dir = "SC3/")
srat_3p   <- CreateSeuratObject(pbmc.data2,project = "SC2")
srat_5p   <- CreateSeuratObject(pbmc.data3,project = "SC3")
#srat_3p
#An object of class Seurat 
#36601 features across 6484 samples within 1 assay 
#Active assay: RNA (36601 features, 0 variable features)
#srat_5p
#An object of class Seurat 
#36601 features across 5648 samples within 1 assay 
#Active assay: RNA (36601 features, 0 variable features)
rm(pbmc.data2)
rm(pbmc.data3)
# QC
meta2 <- srat_3p@meta.data
dim(meta2)
head(meta2)
summary(meta2$nCount_RNA)
summary(meta2$nFeature_RNA)

meta3 <- srat_5p@meta.data
dim(meta3)
head(meta3)
summary(meta3$nCount_RNA)
summary(meta3$nFeature_RNA)

#Let’s calculate the fractions of mitochondrial genes and ribosomal proteins, and do quick-and-dirty filtering of the datasets:

srat_3p[["percent.mt"]]  <- PercentageFeatureSet(srat_3p, pattern = "^MT-")
srat_3p[["percent.rbp"]] <- PercentageFeatureSet(srat_3p, pattern = "^RP[SL]")
srat_5p[["percent.mt"]]  <- PercentageFeatureSet(srat_5p, pattern = "^MT-")
srat_5p[["percent.rbp"]] <- PercentageFeatureSet(srat_5p, pattern = "^RP[SL]")
VlnPlot(srat_3p, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

VlnPlot(srat_5p, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
#смотрим, где заканчивается нормальное распределение #так как количество генов коррелирует с количеством каунтов, то и каунты почистим
srat_3p <- subset(srat_3p, subset = nFeature_RNA > 200 & nFeature_RNA <
                3500 & percent.mt < 15)
srat_5p <- subset(srat_5p, subset = nFeature_RNA > 200 & nFeature_RNA <
                3000 & percent.mt < 20)
VlnPlot(srat_3p, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(srat_5p, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

# Now, let’s follow Seurat vignette for integration. 
# To do this we need to make a simple R list of the two objects, and normalize using SCTransform/find variable genes for each:

pbmc_list <- list()
pbmc_list[["SC2"]] <- srat_3p
pbmc_list[["SC3"]] <- srat_5p
pbmc_list <- lapply(X = pbmc_list, FUN = SCTransform)#трансформация данных для качественной интеграции
features <- SelectIntegrationFeatures(object.list = pbmc_list, nfeatures = 3000)#выбираем фичи, с помощью которых будем интегрировать данные
pbmc_list <- PrepSCTIntegration(object.list = pbmc_list, anchor.features = features)#подготавливаем объект к интеграции
pbmc_anchors    <- FindIntegrationAnchors(object.list = pbmc_list, dims = 1:30)
pbmc_seurat     <- IntegrateData(anchorset = pbmc_anchors, dims = 1:30)
#rm(pbmc_list)
#rm(pbmc_anchors)

# Seurat integration creates a unified object that contains both original data (‘RNA’ assay) as well as integrated data (‘integrated’ assay). 
# Let’s set the assay to RNA and visualize the datasets before integration.

DefaultAssay(pbmc_seurat) <- "RNA"

# Let’s do normalization, HVG finding, scaling, PCA, and UMAP on the un-integrated (RNA) assay:

pbmc_seurat <- NormalizeData(pbmc_seurat, verbose = F)
pbmc_seurat <- FindVariableFeatures(pbmc_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
pbmc_seurat <- ScaleData(pbmc_seurat, verbose = F)
pbmc_seurat <- RunPCA(pbmc_seurat, npcs = 30, verbose = F)
pbmc_seurat <- RunUMAP(pbmc_seurat, reduction = "pca", dims = 1:30, verbose = F)
# UMAP plot of the datasets before integration shows clear separation. 
DimPlot(pbmc_seurat,reduction = "umap")
# The data are visibly very nicely integrated. 
# Let’s try a split plot, which should make the comparison easier:

n_cells <- FetchData(pbmc_seurat, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident)
#Draw a stacked barplot # Save 6x3
ggplot(n_cells, aes(x=orig.ident, y=n, fill=ident)) + 
  geom_bar(position="fill", stat="identity")



