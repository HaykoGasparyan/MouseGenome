library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

library(ggpubr)

library(Seurat)
library(Matrix)
library(ape)

#install.packages('BiocManager')
#BiocManager::install('limma')

my_cols <- c('3'='#F68282','15'='#31C53F','5'='#1FA195','1'='#B95FBB','13'='#D4D915',
             '14'='#28CECA','9'='#ff9a36','8'='#2FF18B','11'='#aeadb3','6'='#faf4cf',
             '2'='#CCB1F1','12'='#25aff5','7'='#A4DFF2','4'='#4B4BF7','16'='#AC8F14',
             '10'='#E6C122')

#d <- read.csv("/home/sargis/project/genCountTable.csv")


#d <- read.table("/home/sargis/project/GSE114000_Counts_Batiuk_Martirosyan_Supplementary_Data3 (1).tsv", sep='\t')
#rownames(d) = d[, 1]
#colnames(d) = d[1, ]
#
#d_new = d[2:49661, 2:1024]
#d = d_new
#pbmc.data = d

pbmc.data = read.csv("/home/sargis/project/CleanedData.csv")
#pbmc.data
pbmc.data = pbmc.data[!duplicated(pbmc.data[, "cell_id"]),]
rownames(pbmc.data) = pbmc.data[,2]
pbmc.data = pbmc.data[, 3:3005]

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "/home/sargis/project/GenCount/")

#QC Preproces
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)




pbmc <- NormalizeData(pbmc)



pbmc <- FindVariableFeatures(pbmc)

pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

####P-value

pbmc <- JackStraw(pbmc, dims = 50)
pbmc <- ScoreJackStraw(pbmc, dims = 1:50)

####JackStrawPlot(pbmc, dims = 1:50)

#plot respect to varianc
ElbowPlot(pbmc, ndims = 50)+ theme(axis.text=element_text(size=20),
                                   axis.title=element_text(size=20,face="bold"),
                                   legend.text = element_text(size=20))


#Dim heatMap

DimHeatmap(pbmc, dims = 1:11, cells = 500, balanced = TRUE) + theme(axis.text=element_text(size=20),
                                                                   axis.title=element_text(size=40,face="bold"))


#clustering
pbmc <- FindNeighbors(pbmc, dims = 1:11)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Build a classification hierarchy that places transcriptionally similar clusters adjacent on a tree

#UMAP
pbmc <- RunUMAP(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- RunTSNE(pbmc, features =  VariableFeatures(object = pbmc))
new.cluster.ids <- c("0", "1", "2", "3", "4", "5","6", "7", "8", "0", "0", "0")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE,label.size = 12, pt.size = 2) + theme(axis.text=element_text(size=20),
                                                                                     axis.title=element_text(size=20,face="bold"),
                                                                                     legend.text = element_text(size=20))

FeatureScatter(pbmc, feature1 = "PC_2", feature2 = "PC_1")

DimHeatmap(pbmc, dims = 1, cells = 100, balanced = TRUE)

mat = matrix(, nrow = 6, ncol = 6)
for (i in c(1,2,3,4,5,6)) {
  for (j in c(1,2,3,4,5,6)) {
    pbmc.markers <- FindMarkers(pbmc, ident.1 = i-1, ident.2 = j-1, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.7)
    mat[i, j] = nrow(pbmc.markers)
  }
}
mat



#Top markers in clasters
top = pbmc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
FeaturePlot(pbmc, reduction = 'tsne',features = top$gene) 

deatMap = DoHeatmap(pbmc, features = top$gene)  + theme(axis.text=element_text(size=15),
                                              axis.title=element_text(size=15,face="bold"),
                                              legend.text = element_text(size=18,face="bold"))

ggarrange(dimPlot,deatMap,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)


astrocytes <- c("GFAP","Atp1b2","Fgfr3","Gjb6", "Gluf","Aqp4", "Slc1a3","Slc1a2", "S100B","Aldh1L1")
microglia <- c("Aif1","ADGRE1", "CD68", "CD40", "cx3cr1", "itgam")
oligodendrocytes <- c("Olig1", "Olig2", "Olig3", "Mbp","Mog","Mag", "Cldn11", "Sox10")
neuron <- c("Rbfox3", "Snap25","Celf4","syt1","Syp","l1cam", "Map2", "Dlg4","Slc17a7", "Grin1", "Slc17a6", "Grin2b", "Slc6a1", "Gad2", "Gad1", "Foxa2", "Kcnj6", "Nr4a2", "Lmx1b")
muralClass <- c("Pdgfrb", "Cspg4", "Pid1", "Des", "Plce1", "Col1a2", "Itga8")

endhothelium <- c("Flt1", "Cd34", "Pecam1", "icam2", "tie1")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.4, min.diff.pct = 0.3, logfc.threshold = 0.5)
clus0 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==0))
clus1 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==1))
clus2 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==2))
clus3 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==3))
clus4 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==4))
clus5 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==5))
clus6 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==6))
clus7 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==7))
clus8 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==8))
clus9 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==9))
clus10 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==10))
clus11 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==11))
clus12 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==12))
clus13 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==13))
clus14 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==14))
clus15 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==15))
clus16 <- rownames(subset(pbmc.markers, pbmc.markers$cluster==16))

for (i in seq.int(from=0,to = length(astrocytes))) {
  astrocytes[i] = tolower(astrocytes[i])
}

for (i in seq.int(from=0,to = length(microglia))) {
  microglia[i] = tolower(microglia[i])
}

for (i in seq.int(from=0,to = length(oligodendrocytes))) {
  oligodendrocytes[i] = tolower(oligodendrocytes[i])
}

for (i in seq.int(from=0,to = length(neuron))) {
  neuron[i] = tolower(neuron[i])
}

for (i in seq.int(from=0,to = length(muralClass))) {
  muralClass[i] = tolower(muralClass[i])
}

for (i in seq.int(from=0,to = length(endhothelium))) {
  endhothelium[i] = tolower(endhothelium[i])
}





for (i in seq.int(from=0,to = length(clus0))) {
  clus0[i] = tolower(clus0[i])
}

for (i in seq.int(from=0,to = length(clus1))) {
  clus1[i] = tolower(clus1[i])
}

for (i in seq.int(from=0,to = length(clus2))) {
  clus2[i] = tolower(clus2[i])
}

for (i in seq.int(from=0,to = length(clus3))) {
  clus3[i] = tolower(clus3[i])
}

for (i in seq.int(from=0,to = length(clus4))) {
  clus4[i] = tolower(clus4[i])
}

for (i in seq.int(from=0,to = length(clus5))) {
  clus5[i] = tolower(clus5[i])
}

for (i in seq.int(from=0,to = length(clus6))) {
  clus6[i] = tolower(clus6[i])
}

for (i in seq.int(from=0,to = length(clus7))) {
  clus7[i] = tolower(clus7[i])
}

for (i in seq.int(from=0,to = length(clus8))) {
  clus8[i] = tolower(clus8[i])
}

length(intersect(clus0, astrocytes))
length(intersect(clus0, microglia))
length(intersect(clus0, oligodendrocytes))
length(intersect(clus0, neuron))
length(intersect(clus0, muralClass))
length(intersect(clus0, endhothelium))
print("----------------")

length(intersect(clus1, astrocytes))
length(intersect(clus1, microglia))
length(intersect(clus1, oligodendrocytes))
length(intersect(clus1, neuron))
length(intersect(clus1, muralClass))
length(intersect(clus1, endhothelium))
print("----------------")

length(intersect(clus2, astrocytes))
length(intersect(clus2, microglia))
length(intersect(clus2, oligodendrocytes))
length(intersect(clus2, neuron))
length(intersect(clus2, muralClass))
length(intersect(clus2, endhothelium))
print("----------------")

length(intersect(clus3, astrocytes))
length(intersect(clus3, microglia))
length(intersect(clus3, oligodendrocytes))
length(intersect(clus3, neuron))
length(intersect(clus3, muralClass))
length(intersect(clus3, endhothelium))
print("----------------")

length(intersect(clus4, astrocytes))
length(intersect(clus4, microglia))
length(intersect(clus4, oligodendrocytes))
length(intersect(clus4, neuron))
length(intersect(clus4, muralClass))
length(intersect(clus4, endhothelium))
print("----------------")

length(intersect(clus5, astrocytes))
length(intersect(clus5, microglia))
length(intersect(clus5, oligodendrocytes))
length(intersect(clus5, neuron))
length(intersect(clus5, muralClass))
length(intersect(clus5, endhothelium))
print("----------------")

length(intersect(clus6, astrocytes))
length(intersect(clus6, microglia))
length(intersect(clus6, oligodendrocytes))
length(intersect(clus6, neuron))
length(intersect(clus6, muralClass))
length(intersect(clus6, endhothelium))
print("----------------")

length(intersect(clus7, astrocytes))
length(intersect(clus7, microglia))
length(intersect(clus7, oligodendrocytes))
length(intersect(clus7, neuron))
length(intersect(clus7, muralClass))
length(intersect(clus7, endhothelium))
print("----------------")

length(intersect(clus8, astrocytes))
length(intersect(clus8, microglia))
length(intersect(clus8, oligodendrocytes))
length(intersect(clus8, neuron))
length(intersect(clus8, muralClass))
length(intersect(clus8, endhothelium))
print("----------------")
new.cluster.ids <- c("Oligodendrocytes", "Neuron", "Unknown", "Neuron", "Unknown", "Microglia","Unknown", "Astrocytes", "Microglia")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE,label.size = 10, pt.size = 1.5) + theme(axis.text=element_text(size=20),
                                                                                     axis.title=element_text(size=20,face="bold"),
                                                                                     legend.text = element_text(size=20))



