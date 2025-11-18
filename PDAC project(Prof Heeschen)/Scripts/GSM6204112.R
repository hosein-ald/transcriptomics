path6204112="C:/Users/abbas/Documents/GSE205013/GSM6204112/GSM6204112_RAW"


# Update seurat object
GSM6204112=Read10X(path6204112)
dim(GSM6204112)

#Step create srobj
GSM6204112obj=CreateSeuratObject(GSM6204112, min.cells = 3, min.features = 100)
dim(GSM6204112obj)

GSM6204112obj[["MTpercent"]]= PercentageFeatureSet(GSM6204112obj, pattern = "^MT-")
VlnPlot(GSM6204112obj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)
plot1= FeatureScatter(GSM6204112obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
GSM6204112obj= subset(GSM6204112obj, subset = nFeature_RNA >200 & nFeature_RNA<8000 & MTpercent <25)
plot1

# step Normalization
GSM6204112obj = NormalizeData(GSM6204112obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Step HighestVariableGenes
GSM6204112obj= FindVariableFeatures(GSM6204112obj, selection.method = "vst", nfeatures = 3000)
topVariableFeatures= head(VariableFeatures(GSM6204112obj), 100)
plot1= VariableFeaturePlot(GSM6204112obj)
plot2= LabelPoints(plot = plot1, points = topVariableFeatures, repel = T)
plot2

# Step Linear Reduction
GSM6204112obj= ScaleData(GSM6204112obj, features = rownames(GSM6204112obj))
GSM6204112obj= RunPCA(GSM6204112obj, features = VariableFeatures(GSM6204112obj))
print(GSM6204112obj[["pca"]], dims= 1:50, nfeatures= 5)
VizDimLoadings(GSM6204112obj, dims = 1:5 , reduction = "pca")
DimPlot(GSM6204112obj, reduction = "pca")
DimHeatmap(GSM6204112obj, dims = 1:6, cells = 100)
ElbowPlot(GSM6204112obj)

# Step Jackstraw and PCA
GSM6204112obj= JackStraw(GSM6204112obj, num.replicate = 100, dims = 50)
GSM6204112obj= ScoreJackStraw(GSM6204112obj, dims = 1:50)
JackStrawPlot(GSM6204112obj, dims = 1:50)
ElbowPlot(GSM6204112obj)
saveRDS(gsm8452584p1obj, file = "C:/Users/abbas/Documents/gsm8452584p1obj.rds")

# Step FindNeighbors
GSM6204112obj= FindNeighbors(GSM6204112obj, dims = 1:50)
GSM6204112obj= FindClusters(GSM6204112obj, resolution = 1.3)
length(levels(Idents(gsm6735860)))
table(Idents(GSM6204112obj))
# Step T-SNE
GSM6204112obj= RunTSNE(GSM6204112obj, dims = 1:13)
DimPlot(GSM6204112obj, reduction = "tsne", label = T)

# Step UMAP
GSM6204112obj= RunUMAP(GSM6204112obj, dims = 1:50)
DimPlot(GSM6204112obj, reduction = "umap",repel = T ,label = T )

GSM6204112obj2=GSM6204112obj
markers_of_GSM6204112 =FindAllMarkers(GSM6204112obj, min.pct = 0.25, logfc.threshold = 0.3)
table(markers_of_gsm452585$cluster)
doubeldmarkersofgsm452585= markers_of_gsm452585[which(markers_of_gsm452585$avg_log2FC>2),]
dim(doubeldmarkersofgsm452585)
table(Idents(GSM6204112obj2))

# ASSIGNING Rename idents
linagesofPDAC=c("NK Cells","Mast Cells","Macrophages","Ductal Cells","CD4+ T Cells","CD8+ T Cells","trNK",
                "Hepatocytes","Fibroblasts","Endothelial", "B Cells", "pDC","Plasma Cells")
GSM6204112objcelltype= c("pDC","B Cells", "Hepatocytes", "Mast Cells", "trNK","Treg",
                         "CD8+ T Cells","CD4+ T Cells","Fibroblasts", "Ductal Cells", "CD4+ T Cells", "Endothelial",
                         "CD8+ T Cells", "Ductal Cells","Fibroblasts", "Muscle Cells","Ductal Cells", "Macrophages",
                         "Macrophages","NK Cells","Fibroblasts", "CD4+ T Cells","B Cells","Ductal Cells",
                         "Endothelial", "Acinar Cells","Ductal Cells")
newLevels6204112=c("Muscle Cells","NK Cells","trNK","Hepatocytes","pDC","Macrophages","CD8+ T Cells","Mast Cells",
            "CD4+ T Cells","Ductal Cells","Treg",
            "Endothelial","Acinar Cells","Fibroblasts","B Cells","Plasma Cells")
length(GSM6204112objcelltype)
length(newLevels6204112)
Idents(GSM6204112obj2) <- factor(Idents(GSM6204112obj2), levels = newLevels6204112)

topGenesOfGSM6204112objClusterCD8TCells=topGenesOfGSM6204112objCluster00
topGenesOfGSM6204112objClusterCD4TCells=topGenesOfGSM6204112objCluster01
topGenesOfGSM6204112objClusterFibroblasts=topGenesOfGSM6204112objCluster02
topGenesOfGSM6204112objClusterMuscleCells=topGenesOfGSM6204112objCluster09
topGenesOfGSM6204112objClusterDuctal=topGenesOfGSM6204112objCluster03
topGenesOfGSM6204112objClusterEndothelial=topGenesOfGSM6204112objCluster05
topGenesOfGSM6204112objClusterMacrophages=topGenesOfGSM6204112objCluster12
topGenesOfGSM6204112objClusterNkCell=topGenesOfGSM6204112objCluster13
topGenesOfGSM6204112objClusterBCell=topGenesOfGSM6204112objCluster16
topGenesOfGSM6204112objClusterAcinar=topGenesOfGSM6204112objCluster19

topOfGSM6204112Heatmap=c(topGenesOfGSM6204112objClusterCD8TCells,topGenesOfGSM6204112objClusterCD4TCells,
                         topGenesOfGSM6204112objClusterFibroblasts,topGenesOfGSM6204112objClusterMuscleCells,
                         topGenesOfGSM6204112objClusterDuctal,topGenesOfGSM6204112objClusterMacrophages,
                         topGenesOfGSM6204112objClusterNkCell,topGenesOfGSM6204112objClusterBCell,
                         topGenesOfGSM6204112objClusterAcinar,topGenesClusterPlasmaMarkersGSM6204112vec,
                         topGenesOfGSM6204112objClusterEndothelial)

topGenesOfGSM6204112objCluster00=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==0),])[1:10]
topGenesOfGSM6204112objCluster01=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==1),])[1:10]
topGenesOfGSM6204112objCluster02=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==2),])[1:10]
topGenesOfGSM6204112objCluster03=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==3),])[1:10]
topGenesOfGSM6204112objCluster04=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==4),])[1:10]
topGenesOfGSM6204112objCluster05=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==5),])[1:10]
topGenesOfGSM6204112objCluster06=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==6),])[1:10]
topGenesOfGSM6204112objCluster07=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==7),])[1:10]
topGenesOfGSM6204112objCluster08=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==8),])[1:10]
topGenesOfGSM6204112objCluster09=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==9),])[1:10]
topGenesOfGSM6204112objCluster10=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==10),])[1:10]
topGenesOfGSM6204112objCluster11=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==11),])[1:10]
topGenesOfGSM6204112objCluster12=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==12),])[1:10]
topGenesOfGSM6204112objCluster13=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==13),])[1:10]
topGenesOfGSM6204112objCluster14=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==14),])[1:10]
topGenesOfGSM6204112objCluster15=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==15),])[1:10]
topGenesOfGSM6204112objCluster16=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==16),])[1:10]
topGenesOfGSM6204112objCluster17=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==17),])[1:10]
topGenesOfGSM6204112objCluster18=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==18),])[1:10]
topGenesOfGSM6204112objCluster19=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==19),])[1:10]
topGenesOfGSM6204112objCluster20=rownames(markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==20),])[1:10]

DoHeatmap(GSM6204112obj2 ,features = topOfGSM6204112Heatmap,  angle = 60) 



length(topGenesOfGSM6204112objCluster00)
  
  
topOfGSM6204112Heatmap=c(topGenesOfGSM6204112objClusterDuctal,topGenesOfGSM6204112objClusterFibroblasts,
                         topGenesOfGSM6204112objClusterEndothelial,topGenesOfGSM6204112objClusterMacrophages,
                         topGenesOfGSM6204112objClusterCD8TCells,topGenesOfGSM6204112objClusterCD4TCells,
                         topGenesOfGSM6204112objClusterNkCell,
                         topGenesOfGSM6204112objClusterMuscleCells,
                         topGenesOfGSM6204112objClusterBCell,topGenesClusterPlasmaMarkersGSM6204112vec,
                         topGenesOfGSM6204112objClusterAcinar)

newLevelsGSM320412=c("Ductal Cells","Fibroblasts","Endothelial","Macrophages","CD8+ T Cells","CD4+ T Cells",
                     "NK Cells","trNK","Muscle Cells","Hepatocytes","pDC","Mast Cells","Treg","B Cells","Plasma Cells", "Acinar Cells")

length(newLevelsGSM320412)
#Idents(GSM6204112obj2) <- factor(Idents(GSM6204112obj2), levels = newLevelsGSM320412)
Idents(GSM6204112obj2) <- factor(Idents(GSM6204112obj2), levels = newLevelsGSM320412)

length(GSM6204112objcelltype)
names(x= GSM6204112objcelltype)= levels(x= GSM6204112obj2)
GSM6204112obj2= RenameIdents(object = GSM6204112obj2, GSM6204112objcelltype)
DimPlot(GSM6204112obj2, reduction = "umap",repel = T ,label = T )
length(levels(x= GSM6204112obj2))
VlnPlot(GSM6204112obj2, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)

# Dimplot Of Empty
DimPlot(GSM6204112obj2, reduction = "umap", label = T)
table(Idents(GSM6204112obj2))

cluster16_markers <- FindMarkers(GSM6204112obj, ident.1 = 16)
cluster20_markers <- FindMarkers(GSM6204112obj, ident.1 = 20)
cluster0markers <- FindMarkers(GSM6204112obj, ident.1 = 0, ident.2 = c(10,13, 3))
cluster3markers <- FindMarkers(GSM6204112obj, ident.1 = 3, ident.2 = c(10, 13,0))

head(cluster13.markers)
View(cluster13.markers)
rownames(cluster16_markers[1:30,])
rownames(cluster20_markers[1:30,])
rownames(cluster0.markers[1:25,])
rownames(cluster3.markers[1:25,])
intersect_of_10_13=intersect(rownames(cluster13.markers[1:30,]),rownames(cluster10.markers[1:30,]))
length(intersect_of_10_13)

#ductal_markers 
FeaturePlot(GSM6204112obj, features = "KRT19", label = T, reduction = "umap") #
FeaturePlot(GSM6204112obj, features = "CDH1", label = T, reduction = "umap") #
FeaturePlot(GSM6204112obj, features = "CLDN4", label = T, reduction = "umap") #
FeaturePlot(GSM6204112obj, features = "AQP1", label = T, reduction = "umap") #
FeaturePlot(GSM6204112obj, features = "SOX9", label = T, reduction = "umap") #
FeaturePlot(GSM6204112obj, features = "PAX6", label = T, reduction = "umap") #
FeaturePlot(GSM6204112obj, features = "TP63", label = T, reduction = "umap") #
FeaturePlot(GSM6204112obj, features = "TP53", label = T, reduction = "umap") #
FeaturePlot(GSM6204112obj, features = "AKT2", label = T, reduction = "umap") #


#acinar cell markers
FeaturePlot(GSM6204112obj, features = "AMY2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CPA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CLPS", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PRSS1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CELA2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MUC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "KRT19", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PSMB5", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PLIN2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "NQO1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "RAB37", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "AQP8", label = T, reduction = "umap")#

#Endocrine cells
FeaturePlot(GSM6204112obj, features = "INS", label = T, reduction = "umap")  # Insulin, Beta cell marker
FeaturePlot(GSM6204112obj, features = "GCG", label = T, reduction = "umap")  # Glucagon, Alpha cell marker
FeaturePlot(GSM6204112obj, features = "SST", label = T, reduction = "umap")   # Somatostatin, Delta cell marker
FeaturePlot(GSM6204112obj, features = "PPY", label = T, reduction = "umap")  # Pancreatic Polypeptide, PP cell marker
FeaturePlot(GSM6204112obj, features = "NKX6-1", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204112obj, features = "MAFA", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204112obj, features = "PAX4", label = T, reduction = "umap")  # Beta cell development marker
FeaturePlot(GSM6204112obj, features = "CDH1", label = T, reduction = "umap")  # Cell adhesion, pan-endocrine marker

FeaturePlot(GSM6204112obj, features = "GHRL", label = T, reduction = "umap")  # Ghrelin gene, marker for epsilon (ε) cells
FeaturePlot(GSM6204112obj, features = "PDX1", label = T, reduction = "umap")  # Pancreatic and duodenal homeobox 1, critical for β cell function and development
FeaturePlot(GSM6204112obj, features = "ARX", label = T, reduction = "umap")  # Marker for alpha (α) cell identity (glucagon-secreting)
FeaturePlot(GSM6204112obj, features = "HHEX", label = T, reduction = "umap") # Marker associated with delta (δ) cells (somatostatin-secreting)


#CD4
FeaturePlot(GSM6204112obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD4", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IL7R", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CCR7", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "SELL", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "FOXP3", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IL2RA", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IL17A", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "RORC", label = T, reduction = "umap")#

#CD8
FeaturePlot(GSM6204112obj, features = "CD8A", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD8B", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "ITGA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "ZNF683", label = T, reduction = "umap")#
#Treg
FeaturePlot(GSM6204112obj, features = "FOXP3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "IL2RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "CTLA4", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "IKZF2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "TIGIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "IL10", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "AREG", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "IL1RL1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "ITGAE", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "CD69", label = TRUE, reduction = "umap")


#nk
FeaturePlot(GSM6204112obj, features = "NCAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "NKG7", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "GNLY", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "KLRD1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "EOMES", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "TBX21", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "FCGR3A", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "FGFBP2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CX3CR1", label = T, reduction = "umap")#

#B cells
FeaturePlot(GSM6204112obj, features = "CD19", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MS4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD79A", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD79B", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PAX5", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD22", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IGHM", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IGHG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MZB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "JCHAIN", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PRDM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "XBP1", label = T, reduction = "umap")#

#Macrophage
FeaturePlot(GSM6204112obj, features = "CD163", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MRC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MSR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "ARG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CCL18", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CHI3L1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "VEGFA", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "TGFB1", label = T, reduction = "umap")#

#M1 Macrophages
FeaturePlot(GSM6204112obj, features = "CD80", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD86", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IL1B", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "TNF", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IL6", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "NOS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CXCL9", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CXCL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "STAT1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IRF5", label = T, reduction = "umap")#

#cDC
FeaturePlot(GSM6204112obj, features = "CLEC9A", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "XCR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "BATF3", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "IRF8", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CADM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "THBD", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "HLA-DQB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD207", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD74", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CLNK", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD14", label = T, reduction = "umap")#

#fibroblast
FeaturePlot(GSM6204112obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "DCN", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "LUM", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PDGFRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "VIM", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "FAP", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "FN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "SPP1", label = T, reduction = "umap")#

#myofibroblast
FeaturePlot(GSM6204112obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "TAGLN", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CNN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MYH11", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "FN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "POSTN", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PDGFRA", label = T, reduction = "umap")#

#Erythroblast umap
FeaturePlot(GSM6204112obj,features = "BPGM", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "TRIM58", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "XPO7", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "HBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "HBA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "ALAS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "GYPA", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "AHSP", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "KLF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "SLC4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "HEMGN", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "CA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "TFRC", label = T, reduction = "umap")#

#Muscle cell
FeaturePlot(GSM6204112obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "TAGLN", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MYH11", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CNN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "DES", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "TPM2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "LMOD1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MYLK", label = T, reduction = "umap")#

#pDC
FeaturePlot(GSM6204112obj, features = "IL3RA", label = T, reduction = "umap") # CD123, canonical pDC marker
FeaturePlot(GSM6204112obj, features = "CLEC4C", label = T, reduction = "umap") # BDCA-2, highly specific to pDC
FeaturePlot(GSM6204112obj, features = "LILRA4", label = T, reduction = "umap") # ILT7, inhibitory receptor, pDC hallmark
FeaturePlot(GSM6204112obj, features = "TCF4", label = T, reduction = "umap")   # Lineage-defining TF
FeaturePlot(GSM6204112obj, features = "IRF7", label = T, reduction = "umap")   # High in pDC, type I IFN amplifier
FeaturePlot(GSM6204112obj, features = "TLR7", label = T, reduction = "umap")   # Viral RNA sensor
FeaturePlot(GSM6204112obj, features = "TLR9", label = T, reduction = "umap")   # CpG DNA sensor
FeaturePlot(GSM6204112obj, features = "BST2", label = T, reduction = "umap")   # IFN-inducible, en

#Endothelial umap
FeaturePlot(GSM6204112obj,features = "PECAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "VWF", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "CDH5", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "KDR", label = T, reduction = "umap")#

#Mast cell
FeaturePlot(GSM6204112obj,features = "HDC", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "HPGDS", label = T, reduction = "umap")#

#Pericyte umap
FeaturePlot(GSM6204112obj,features = "PDGFRB", label = T, reduction = "umap")#24 7
FeaturePlot(GSM6204112obj,features = "CSPG4", label = T, reduction = "umap")#21 7
FeaturePlot(GSM6204112obj,features = "RGS5", label = T, reduction = "umap")#24
FeaturePlot(GSM6204112obj,features = "MCAM", label = T, reduction = "umap")#24
FeaturePlot(GSM6204112obj,features = "TAGLN", label = T, reduction = "umap")#24
#Adipocyte umap
FeaturePlot(GSM6204112obj,features = "PPARG", label = T, reduction = "umap")#12 21
FeaturePlot(GSM6204112obj,features = "CEBPA", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "ADIPOQ", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "APOE", label = T, reduction = "umap")#7 20 14
FeaturePlot(GSM6204112obj,features = "LPL", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "CFD", label = T, reduction = "umap")#7

#Neutrophil cell
FeaturePlot(GSM6204112obj,features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj,features = "CYBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "S100A8", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "S100A9", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "ELANE", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "AZU1", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "PRTN3", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CD177", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CXCR2", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "FCGR3B", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "CEACAM8", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "LTF", label = T, reduction = "umap")#
FeaturePlot(GSM6204112obj, features = "DEFA1", label = T, reduction = "umap")#
#CMP
FeaturePlot(GSM6204112obj, features = "CD34", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "KIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "SPI1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "CEBPA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "GATA2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "FLT3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "IL3RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "CSF1R", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "LY86", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204112obj, features = "KIT", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204112obj, features = "CD38", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "CD14", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204112obj, features = "KIT", label = TRUE, reduction = "umap") 

#Hepatocyte
FeaturePlot(GSM6204112obj, features = "ALB", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "ASGR1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "HNF4A", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "FABP1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "APOA1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "CYP3A4", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "GPC3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "AFP", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "KRT8", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204112obj, features = "KRT18", label = TRUE, reduction = "umap")

table(Idents(GSM6204112obj2))
GSM6204112obj2=GSM6204112obj
plotcdsevenGSM6204112obj2= FeaturePlot(GSM6204112obj2, features = "MZB1", label = T, reduction = "umap")
HoverLocator(plot = plotcdsevenGSM6204112obj2, information = FetchData(GSM6204112obj2, vars = c("ident", "PC_1", "nFeature_RNA")))

#Endothelial assigning
selected.cellsAsEndothelialcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsEndothelialcell_CellGSM6204112obj2) <- "Endothelial"


#Ductal assigning
selected.cellsAsDuctalcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsDuctalcell_CellGSM6204112obj2) <- "Ductal Cells"


#Plasma Cells assigning
selected.cellsAsPlasmaCell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsPlasmaCell_CellGSM6204112obj2) <- "Plasma Cells"

#pDC assigning
selected.cellsAspDCcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAspDCcell_CellGSM6204112obj2) <- "pDC"

#Erythroblast assigning
selected.cellsAsNKcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsNKcell_CellGSM6204112obj2) <- "Erythroblast"

#Fibroblast assigning
selected.cellsAsfibroblastcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
length(selected.cellsAsfibroblastcell_CellGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsfibroblastcell_CellGSM6204112obj2) <- "Fibroblast"

#NK Cells assigning
selected.cellsAsNKcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsNKcell_CellGSM6204112obj2) <- "NK Cells"

#trNK Cells assigning
selected.cellsAsTrNKcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsTrNKcell_CellGSM6204112obj2) <- "trNK"

#Treg Cells assigning
selectedcellsAsTregCell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
Idents(GSM6204112obj2, cells = selectedcellsAsTregCell_CellGSM6204112obj2) <- "Treg"

#Mast assigning
selected.cellsAsMastcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
length(selected.cellsAsMastcell_CellGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsMastcell_CellGSM6204112obj2) <- "Mast Cells"


#Hepatocyte assigning
selected.cellsAsHepatocytecell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
length(selected.cellsAsMastcell_CellGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsHepatocytecell_CellGSM6204112obj2) <- "Hepatocyte"
#B cell assigning
selected.cellsAsBcell_CellGSM6204112obj2=CellSelector(plot = plotcdsevenGSM6204112obj2)
length(selected.cellsAsBcell_CellGSM6204112obj2)
Idents(GSM6204112obj2, cells = selected.cellsAsBcell_CellGSM6204112obj2) <- "B Cells"

#Hepatocyte assigning
Idents(GSM6204112obj2, cells = plasmalike) <- "Plasma Cells"


fibrolikes=intersect(selected.cellsAsfibroblastcell_CellGSM6204112obj2,
                     selected.cellsAsEndothelialcell_CellGSM6204112obj2)
selected.cellsAsMastcell_CellGSM6204112obj2[1:2]
length(fibrolikes)
plasmalike=c("GACCCTTGTGATGGCA-1","AGGGTGAAGGTAGCAC-1","GACCCTTGTGATGGCA-1","CTGATCCGTAACTTCG-1")

#=================================================
#=================================================
macrophage_markers6204112=c("IL1B", "CD86", "CLEC10A", "IL23A","CXCL5","TLR2","TLR4","STAT1",
                        "CD163", "CD68", "MMP9", "MSR1", "MRC1","CCL18", "CSF1R")
m1featureForFeaturePlot=c("IL1B", "CD86", "CLEC10A", "IL23A","CXCL5","TLR2","TLR4","STAT1", "IRF4")
m2featureForFeaturePlot=c("CD163", "CD68", "MMP9", "MSR1", "MRC1","CCL18", "CSF1R","SIGLEC1","IL10")

FeaturePlot(subsetMacrophages_GSM6204112, features = m1featureForFeaturePlot , label = F, reduction = "umap")
FeaturePlot(subsetMacrophages_GSM6204112, features = m2featureForFeaturePlot , label = F, reduction = "umap")

#=================================================
#=================================================

rownames( markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==17 & markers_of_GSM6204112$avg_log2FC>9),])[1:25]
View( markers_of_GSM6204112[which(markers_of_GSM6204112$cluster==17),])



# Assays
# Get expression data for M1 cells
# Method 2: Create separate Seurat objects for M1 and M2

subsetMacrophages_GSM6204112=subset(GSM6204112obj2, idents = "Macrophages")
m1featureForFeaturePlot=
m2featureForFeaturePlot=
m1_obj_GSM6204112 <- subset(subsetMacrophages_GSM6204112, idents ="M1 Macrophage" )
m2_obj_GSM6204112 <- subset(subsetMacrophages_GSM6204112, idents = "M2 Macrophage" )
TransitionalM1_obj_GSM6204112<-subset(subsetMacrophages_GSM6204112,  idents = "Phenotype-Switching Macrophages")


m1_expressionGSM6204112 <- GetAssayData(m1_obj_GSM6204112, slot = "data", assay = "RNA")
m2_expressionGSM6204112 <- GetAssayData(m2_obj_GSM6204112, slot = "data", assay = "RNA")
TransitionalM1_expressionGSM6204112 <- GetAssayData(TransitionalM1_obj_GSM6204112, slot = "data",assay = "RNA")
c("M1 Macrophage", "M2 Macrophage","Phenotype-Switching Macrophages")
dim(m1_expressionGSM6204112)

m1_expression[1:5,1:5]
write.csv(m1_expressionGSM6204112,file = "C:/Users/abbas/Documents/GSE205013/GSM6204112/GSM6204112_RAW/Macrophages_M1_GSM6204112_data.csv" )
write.csv(m2_expressionGSM6204112,file = "C:/Users/abbas/Documents/GSE205013/GSM6204112/GSM6204112_RAW/Macrophages_M2_GSM6204112_data.csv"  )
write.csv(TransitionalM1_expressionGSM6204112,file = "C:/Users/abbas/Documents/GSE205013/GSM6204112/GSM6204112_RAW/TransitionalM1_expressionGSM6204112_data.csv")

phynotyeSwitchingMacrophagesOfGSM6204112=read.csv("C:/Users/abbas/Documents/GSE205013/GSM6204112/GSM6204112_RAW/Macrophages_M1_GSM6204112_data.csv")
length(phynotyeSwitchingMacrophagesOfGSM6204112)