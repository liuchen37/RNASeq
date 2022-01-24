# Generating a good looking cell proportion plot
# demo data = pbmc3k, Seurat

rm(list = ls())

file.choose()
setwd("directory")

# load pbmc3k
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gginnards)
library(ggstatsplot)
library(SeuratData)
library(Seurat)
data('pbmc3k')
sce <- pbmc3k.final

# UMAP
DimPlot(sce, reduction = 'umap', label = TRUE)
# determine cell types
unique(Idents(sce))
sce$celltype = Idents(sce)

# log transformation
sce <- NormalizeData(sce, normalization.method =  "LogNormalize",  
                     scale.factor = 1e4)
# local polynomial regress, log(variance) vs log(mean)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce)
# principal component analysis
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
# nearest neighbor construction (shared nearest neighbor) and clusters
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.8) 
set.seed(123)
# tSNE and UMAP optimization
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE) 
sce <- RunUMAP(object = sce, dims = 1:15, do.fast = TRUE)

DimPlot(object = sce, reduction = "umap",label = TRUE) 

library(dplyr)
glimpse(sce)

# define groups
Idents(sce) = sce$celltype

# find significant expressed genes for each subgroup
sce.markers <- FindAllMarkers(object = sce, 
                              only.pos = TRUE,
                              min.pct = 0.25, 
                              thresh.use = 0.25 )
# output object
save(sce, sce.markers, file = 'tmp.Rdata')

# define a PropPlot function
PropPlot <- function(object, groupBy){
  # retrieve data
  plot_data = object@meta.data %>% 
    dplyr::select(orig.ident, {{groupBy}}) %>% 
    dplyr::rename(group = as.name(groupBy))
  
  # plot
  figure = ggbarstats(data = plot_data, 
                      x = group, y = orig.ident,
                      package = 'ggsci',
                      palette = 'category20c_d3',
                      results.subtitle = FALSE,
                      bf.message = FALSE,
                      proportion.test = FALSE,
                      label.args = list(size = 2, 
                                        fill = 'white', 
                                        alpha = 0.85,
                                        family = 'Arial',
                                        fontface = 'bold'),
                      perc.k = 2,
                      title = '',
                      xlab = '',
                      legend.title = 'Seurat Cluster',
                      ggtheme = ggpubr::theme_pubclean()) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = 'black', lineend = 'round'),
          legend.position = 'right',
          axis.text.x = element_text(size = 15, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 15, color = 'black', family = 'Arial'),
          legend.text = element_text(family = 'Arial', size = 10, color = 'black'),
          legend.title = element_text(family = 'Arial', size = 13, color = 'black')) 
  
  # remove GeomText
  gginnards::delete_layers(x = figure, match_type = 'GeomText')
}

# plot
source('PropPlot.R')
table(sce$seurat_clusters, sce$celltype)
PropPlot(object = sce, groupBy = 'celltype') + 
  PropPlot(object = sce, groupBy = 'seurat_clusters')

# Result: https://github.com/liuchen37/Pics/blob/main/PropPlot.png
