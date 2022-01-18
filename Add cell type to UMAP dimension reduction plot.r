# Add cell type to UMAP dimension reduction plot
# demo data = pbmc3k, Seurat

rm(list = ls())

# load pbmc3k
library(SeuratData)
library(Seurat)
InstallData('pbmc3k')
data('pbmc3k')
pbmc3k <- pbmc3k.final

# generate a dimension reduction  plot
DimPlot(pbmc3k, reduction = "umap", label = TRUE)
# retrieve cell type
unique(Idents(pbmc3k))
pbmc3k$celltype = Idents(pbmc3k)

# log transformation
pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize",
                        scale.factor = 1e4)
# local polynomial regress, log(variance) vs log(mean)
pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst",
                               nfeatures = 2000)
pbmc3k <- ScaleData(pbmc3k)
# principal component analysis
pbmc3k <- RunPCA(pbmc3k, pc.genes = VariableFeatures(pbmc3k))
# nearest neighbor construction (shared nearest neighbor) and clusters
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:15)
pbmc3k <- FindClusters(pbmc3k, resolution = 0.8)
set.seed(123)
# get real coordinate for rSNE and UMAP
pbmc3k <- RunTSNE(pbmc3k, dims = 1:15, do.fast = FALSE)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:15, do.fast = FALSE)

DimPlot(pbmc3k, reduction = "umap", label = TRUE)

library(dplyr)
glimpse(pbmc3k)

# define groups
Idents(pbmc3k) <- pbmc3k$celltype

# find significant expressed genes for each subgroup
markers <- FindAllMarkers(pbmc3k,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          thresh.use = 0.25)
# output object
save(pbmc3k, file = "pbmc3k.rdata")

# extract expression matrix and dimension reduction coordinates
pbmc3k.all <- pbmc3k
gene = 'CD4'
p1 <- DimPlot(pbmc3k.all, reduction = 'umap', label = T, repel = T, group.by = 'celltype')
pos <- pbmc3k.all@reductions$umap@cell.embeddings
pos <- pos[pbmc3k.all@assays$RNA@counts[gene, ] > 1, ]

# plot
library(ggplot2)
head(pos)
p1 + geom_point(aes(x = UMAP_1, y = UMAP_2),
                shape = 21, colour = "black",
                fill = "blue", size = 0.35,
                data = as.data.frame(pos)) +
  ggtitle("Cell type")

# Result: https://github.com/liuchen37/Pics/blob/main/UMPA%2Bcelltype.png
