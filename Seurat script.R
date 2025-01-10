# Load libraries
library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(DoubletFinder)


# set seed for reproducibility
set.seed(123)

# Set the working directory
getwd()

# get data location
dirs <- list.dirs(path = '/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/CD72 Adila scRNAseq/RAW/', recursive = F, full.names = F)
dirs

# load the data and assign them a SeuratObject
for (x in dirs) {
   cts <- ReadMtx(
   mtx = paste0('/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/CD72 Adila scRNAseq/RAW/', x, '/matrix.mtx'),
   features = paste0('/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/CD72 Adila scRNAseq/RAW/', x, '/features.tsv'),
   cells = paste0('/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/CD72 Adila scRNAseq/RAW/', x, '/barcodes.tsv'),
   feature.column = 1
   )
# cells = # create seurat objects
 assign(x, CreateSeuratObject(counts = cts))
}  

# merge datasets
ls()

merged_seurat <- merge(EmptyCAR, y = c(H24, nbD4_13),
                       add.cell.ids = ls()[c(3:5)],
                       project = 'CD72 Adila scRNAseq')

merged_seurat

# QC & filtering 
view(merged_seurat@meta.data)

# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('CART', 'Barcode'), 
                                    sep = '_')

# Checking the datasets are being merged properly
unique(merged_seurat@meta.data$CART)
unique(merged_seurat@meta.data$Barcode)

# Calculate percentage of mitochondrial genes
merged_seurat[["mitoPercent"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   nFeature_RNA <= 50000 &
                                   nCount_RNA >= 1000 &
                                   nCount_RNA <= 10000 &
                                   mitoPercent < 10)

merged_seurat_filtered
merged_seurat

# Normalize data 
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)

# standard workflow steps
merged_seurat_filtered <- FindVariableFeatures(merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20, reduction = 'pca')

# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'CART', raster = FALSE) + ggtitle("UMAP by CAR-T")
p1

#saveRDS(merged_seurat_filtered, file = '/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/CD72 Adila scRNAseq/CD72_merged_seurat_filtered.rds')
#merged_seurat_filtered<-readRDS('/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/CD72 Adila scRNAseq/CD72_merged_seurat_filtered.rds')
unique(merged_seurat_filtered@meta.data$CART)

########################
##### run Harmony ######
########################

seurat.harmony <- merged_seurat_filtered %>%
   RunHarmony(group.by.vars = 'CART', plot_convergence = FALSE)

seurat.harmony@reductions

seurat.harmony.embed <- Embeddings(seurat.harmony, "harmony")
seurat.harmony.embed[1:10,1:10]

# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
seurat.harmony <- seurat.harmony %>%
   RunUMAP(reduction = 'harmony', dims = 1:20) %>%
   FindNeighbors(reduction = "harmony", dims = 1:20) %>%
   FindClusters(resolution = 0.5)

# visualize 
p2 <- DimPlot(seurat.harmony, reduction = 'umap', group.by = 'CART', raster=FALSE)+ ggtitle("UMAP by CAR-T after batch reduction using Harmony")
p2

p3 <- DimPlot(seurat.harmony, reduction = 'umap', label = TRUE, raster=FALSE)
p3

grid.arrange(p1, p2, ncol = 2, nrow = 1)
grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

#saveRDS(seurat.harmony, file = "CD72_scRNAseq_Harmony_integrated.rds")
seurat.harmony
