library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(presto)
library(biomaRt)


# Loading RDS
seurat.harmony <-readRDS("CD72_scRNAseq_Harmony_integrated.rds")
view(seurat.harmony@meta.data)

# Run UMAP to see clusters
p3 <- DimPlot(seurat.harmony, reduction = 'umap', label = TRUE, raster=FALSE)
p3

p4<-DimPlot(seurat.harmony, reduction = 'umap', label = TRUE, raster=FALSE, split.by = "CART")
p4

# Connect to Ensembl database using biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 109)


# Retrieve gene annotations (ENSG to gene symbols)
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(seurat.harmony),
  mart = ensembl
)

# Map Ensembl IDs to gene symbols for the RNA assay
current_features <- rownames(seurat.harmony@assays$RNA@counts)
mapped_features <- gene_annotations$hgnc_symbol[
  match(current_features, gene_annotations$ensembl_gene_id)
]

# Replace Ensembl IDs with gene symbols (keep unmapped as is)
rownames(seurat.harmony@assays$RNA@counts) <- ifelse(
  is.na(mapped_features),
  current_features,
  mapped_features
)

rownames(seurat.harmony@assays$RNA@data) <- rownames(seurat.harmony@assays$RNA@counts)


# Find all markers
#seurat.harmony.markers <- FindAllMarkers(seurat.harmony, only.pos = TRUE)
seurat.harmony.markers <- FindAllMarkers(seurat.harmony,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Filter markers by avg_log2FC > 1 and group by cluster
filtered_markers <- seurat.harmony.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5)

# save the filtered markers
#saveRDS(seurat.harmony.markers, file = "CD72_seurat_harmony_markers.rds")
#seurat.harmony.markers<-readRDS("CD72_seurat_harmony_markers.rds")

seurat.harmony.markers%>%filter()
seurat.harmony.markers%>%filter(cluster == 0)

#General plots
Idents(seurat.harmony) <- "CART"

# Mantle cell lymphoma cells: CD19, CD20, CD22, CD79A, PAX5, CD5, BCL2, FMC7, CyclinD1, SOX11 - Negative for CD10, BCL6, CD23, IgMD
VlnPlot(seurat.harmony, features = c("CD19", "MS4A1", "CD22", "CD79A", "PAX5", "CD5", "BCL2", "CCND1", "SOX11", "CD72"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("CD19", "MS4A1", "CD79A", "CCND1", "CD72"))
DotPlot(seurat.harmony, features = c("CD19", "MS4A1", "CD22", "CD79A", "PAX5", "CD5", "BCL2", "CCND1", "SOX11", "CD72")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("CD79B", "MS4A1", "CD79A"))
VlnPlot(seurat.harmony, features = c("CD72", "CD19", "MS4A1", "CD79A", "PAX5", "CCND1"), pt.size = 0, log=F)


# CD8+ T-cells: CD8A, CD8B, CD7, CD3E# CD8+ T-cells: CD8A, CD8B, CD7, CD3E - CD4+ T-cells:"CD4", "IL7R", "CD7", "CD3E"
VlnPlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD4", "IL7R", "CD3E"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD3E", "CD4", "IL7R"))
DotPlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD3E", "CD4", "IL7R")) + RotatedAxis()


# Split by CAR-T type

# Mantle cell lymphoma cells: CD19, CD20, CD22, CD79A, PAX5, CD5, BCL2, FMC7, CyclinD1, SOX11 - Negative for CD10, BCL6, CD23, IgMD
VlnPlot(seurat.harmony, features = c("ELOB", "CD72"), split.by = "CART", pt.size = 0, log=F) + 
  theme(legend.position = "right")  
FeaturePlot(seurat.harmony, features = c("CD19", "MS4A1", "CD79A", "CCND1", "CD72"), split.by = "CART")
DotPlot(seurat.harmony, features = c("CD19", "MS4A1", "CD79A", "CCND1", "CD72"), cols = c("orange", "red", "blue"), split.by = "CART") + RotatedAxis()


# CD8+ T-cells: CD8A, CD8B, CD7, CD3E# CD8+ T-cells: CD8A, CD8B, CD7, CD3E - CD4+ T-cells:"CD4", "IL7R", "CD7", "CD3E"
VlnPlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD4", "IL7R", "CD3E"),  split.by = "CART", pt.size = 0, log=F)+ 
  theme(legend.position = "right")  
FeaturePlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD3E", "CD4", "IL7R"), split.by = "CART")
DotPlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD3E", "CD4", "IL7R"), cols = c("orange", "red", "blue"), split.by = "CART") + RotatedAxis()


# Compare featuresplot between CAR-Ts
FeaturePlot(seurat.harmony, features = c("CD19", "MS4A1", "CD79A", "CCND1", "CD72"), split.by = "CART", min.cutoff = 0, max.cutoff = 6)
FeaturePlot(seurat.harmony, features = c("CD3E", "CD8A", "CD4", "ILR7"), split.by = "CART", min.cutoff = 0, max.cutoff = 6)

FeaturePlot(seurat.harmony, features = c("SDC1", "TNFRSF17"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)
FeaturePlot(seurat.harmony, features = c("GPRC5D", "TNFRSF17"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)

FeaturePlot(seurat.harmony, features = c("CD4", "CD25", "FOXP3"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)
FeaturePlot(seurat.harmony, features = c("CD72", "CD19"), split.by = "Patient")


