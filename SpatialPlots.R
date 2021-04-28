
library(Seurat)
library(Matrix)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

brain = Load10X_Spatial('/home/sargis', 'V1_Adult_Mouse_Brain_raw_feature_bc_matrix.h5')
brain


plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# rerun normalization to store sctransform residuals for all genes
brain <- SCTransform(brain, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
# also run standard log normalization for comparison
brain <- NormalizeData(brain, verbose = FALSE, assay = "Spatial")


#plot respect to varianc
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)

ElbowPlot(brain, ndims = 50)+ theme(axis.text=element_text(size=20),
                                   axis.title=element_text(size=20,face="bold"),
                                   legend.text = element_text(size=20))


brain <- FindNeighbors(brain, reduction = "tsne", dims = 1:12)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "tsne", dims = 1:12)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2


