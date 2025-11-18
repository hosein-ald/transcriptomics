pathGSM6204110="C:/Users/abbas/Documents/GSE205013/GSM6204110/GSM6204110_RAW"


# Update seurat object
GSM6204110=Read10X(pathGSM6204110)
dim(GSM6204110)

#Step create srobj
GSM6204110obj=CreateSeuratObject(GSM6204110, min.cells = 3, min.features = 100)
dim(GSM6204110obj)

GSM6204110obj[["MTpercent"]]= PercentageFeatureSet(GSM6204110obj, pattern = "^MT-")
VlnPlot(GSM6204110obj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)
plot1= FeatureScatter(GSM6204110obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
GSM6204110obj= subset(GSM6204110obj, subset = nFeature_RNA >200 & nFeature_RNA<7500 & MTpercent <35 & nCount_RNA<35000)
plot1

# step Normalization
GSM6204110obj = NormalizeData(GSM6204110obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Step HighestVariableGenes
GSM6204110obj= FindVariableFeatures(GSM6204110obj, selection.method = "vst", nfeatures = 3000)
topVariableFeatures= head(VariableFeatures(GSM6204110obj), 40)
plot1= VariableFeaturePlot(GSM6204110obj)
plot2= LabelPoints(plot = plot1, points = topVariableFeatures, repel = T)
plot2

# Step Linear Reduction
GSM6204110obj= ScaleData(GSM6204110obj, features = rownames(GSM6204110obj))
GSM6204110obj= RunPCA(GSM6204110obj, features = VariableFeatures(GSM6204110obj))
print(GSM6204110obj[["pca"]], dims= 1:50, nfeatures= 5)
VizDimLoadings(GSM6204110obj, dims = 1:5 , reduction = "pca")
DimPlot(GSM6204110obj, reduction = "pca")
DimHeatmap(GSM6204110obj, dims = 1:6, cells = 100)
ElbowPlot(GSM6204110obj)

# Step Jackstraw and PCA
GSM6204110obj= JackStraw(GSM6204110obj, num.replicate = 100, dims = 50)
GSM6204110obj= ScoreJackStraw(GSM6204110obj, dims = 1:50)
JackStrawPlot(GSM6204110obj, dims = 1:50)
ElbowPlot(GSM6204110obj)

# Step FindNeighbors
GSM6204110obj= FindNeighbors(GSM6204110obj, dims = 1:50)
GSM6204110obj= FindClusters(GSM6204110obj, resolution = 1.3)
length(levels(Idents(gsm6735860)))
table(Idents(GSM6204110obj))
# Step T-SNE
GSM6204110obj= RunTSNE(GSM6204110obj, dims = 1:13)
DimPlot(GSM6204110obj, reduction = "tsne", label = T)
# Step UMAP
GSM6204110obj= RunUMAP(GSM6204110obj, dims = 1:50)
DimPlot(GSM6204110obj, reduction = "umap",repel = T ,label = T )

GSM6204110obj2=GSM6204110obj
markers_of_GSM6204110 =FindAllMarkers(GSM6204110obj, min.pct = 0.25, logfc.threshold = 0.3)
table(markers_of_gsm452585$cluster)
doubeldmarkersofgsm452585= markers_of_gsm452585[which(markers_of_gsm452585$avg_log2FC>2),]
dim(doubeldmarkersofgsm452585)

# ASSIGNING Rename idents

linagesofPDAC=c("NK Cells","Mast Cells","Macrophages","Ductal Cells","CD4+ T Cells","CD8+ T Cells","trNK",
                "Hepatocytes","Fibroblasts","Endothelial", "B Cells", "pDC","Plasma Cells")

GSM6204110objcelltype= c("pDC", "Hepatocytes", "Endothelial", "Fibroblasts","Mast Cells",
                         "CD4+ T Cells","Ductal Cells","CD4+ T Cells", "CD8+ T Cells", "NK Cells", "Ductal Cells",
                         "Treg", "trNK","CD8+ T Cells", "NK Cells","Ductal Cells", "Macrophages",
                         "Ductal Cells","Macrophages","Macrophages", "trNK","B Cells","Ductal Cells",
                         "Plasma Cells")

length(GSM6204110objcelltype)
newLevelsGSM6204110=c("NK Cells","trNK","Hepatocytes","pDC","Macrophages","CD8+ T Cells","Mast Cells","CD4+ T Cells","Ductal Cells","Treg",
            "Endothelial","Fibroblasts","B Cells","Plasma Cells")
length(newLevelsGSM6204110)
Idents(GSM6204110obj2) <- factor(Idents(GSM6204110obj2), levels = newLevels)
#Idents(GSM6204110obj2) <- factor(Idents(GSM6204110obj2), levels = newLevelsGSM320412)
table(Idents(GSM6204110obj2))

length(GSM6204110objcelltype)
names(x= GSM6204110objcelltype)= levels(x= GSM6204110obj2)
GSM6204110obj2= RenameIdents(object = GSM6204110obj2, GSM6204110objcelltype)
DimPlot(GSM6204110obj2, reduction = "umap",repel = T ,label = T )
length(levels(x= GSM6204110obj2))

# Dimplot Of Empty
DimPlot(GSM6204110obj2, reduction = "umap", label = T)

cluster17.markers <- FindMarkers(GSM6204110obj, ident.1 = 17, ident.2 = c(9, 3,14,4))
cluster10_markers <- FindMarkers(GSM6204110obj, ident.1 = 10, ident.2 = c(0, 3,13))
cluster0.markers <- FindMarkers(GSM6204110obj, ident.1 = 0, ident.2 = c(10,13, 3))
cluster3.markers <- FindMarkers(GSM6204110obj, ident.1 = 3, ident.2 = c(10, 13,0))

clusterPlasmaMarkersGSM6204110 <- FindMarkers(GSM6204110obj2, ident.1 = "Plasma Cells")
topGenesClusterPlasmaMarkersGSM6204110vec=rownames(clusterPlasmaMarkersGSM6204110)[1:10]
head(cluster13.markers)
View(cluster13.markers)
rownames(cluster17.markers[1:30,])
rownames(cluster10.markers[1:30,])
rownames(cluster0.markers[1:25,])
rownames(cluster3.markers[1:25,])
intersect_of_10_13=intersect(rownames(cluster13.markers[1:30,]),rownames(cluster10.markers[1:30,]))
length(intersect_of_10_13)

#ductal_markers 
FeaturePlot(GSM6204110obj, features = "KRT19", label = T, reduction = "umap") #
FeaturePlot(GSM6204110obj, features = "CDH1", label = T, reduction = "umap") #
FeaturePlot(GSM6204110obj, features = "CLDN4", label = T, reduction = "umap") #
FeaturePlot(GSM6204110obj, features = "AQP1", label = T, reduction = "umap") #
FeaturePlot(GSM6204110obj, features = "SOX9", label = T, reduction = "umap") #
FeaturePlot(GSM6204110obj, features = "PAX6", label = T, reduction = "umap") #
FeaturePlot(GSM6204110obj, features = "TP63", label = T, reduction = "umap") #
FeaturePlot(GSM6204110obj, features = "TP53", label = T, reduction = "umap") #
FeaturePlot(GSM6204110obj, features = "AKT2", label = T, reduction = "umap") #


#acinar cell markers
FeaturePlot(GSM6204110obj, features = "AMY2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CPA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CLPS", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PRSS1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CELA2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "MUC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "KRT19", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PSMB5", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PLIN2", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "NQO1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "RAB37", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "AQP8", label = T, reduction = "umap")#

#Endocrine cells
FeaturePlot(GSM6204110obj, features = "INS", label = T, reduction = "umap")  # Insulin, Beta cell marker
FeaturePlot(GSM6204110obj, features = "GCG", label = T, reduction = "umap")  # Glucagon, Alpha cell marker
FeaturePlot(GSM6204110obj, features = "SST", label = T, reduction = "umap")   # Somatostatin, Delta cell marker
FeaturePlot(GSM6204110obj, features = "PPY", label = T, reduction = "umap")  # Pancreatic Polypeptide, PP cell marker
FeaturePlot(GSM6204110obj, features = "NKX6-1", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204110obj, features = "MAFA", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204110obj, features = "PAX4", label = T, reduction = "umap")  # Beta cell development marker
FeaturePlot(GSM6204110obj, features = "CDH1", label = T, reduction = "umap")  # Cell adhesion, pan-endocrine marker

FeaturePlot(GSM6204110obj, features = "GHRL", label = T, reduction = "umap")  # Ghrelin gene, marker for epsilon (ε) cells
FeaturePlot(GSM6204110obj, features = "PDX1", label = T, reduction = "umap")  # Pancreatic and duodenal homeobox 1, critical for β cell function and development
FeaturePlot(GSM6204110obj, features = "ARX", label = T, reduction = "umap")  # Marker for alpha (α) cell identity (glucagon-secreting)
FeaturePlot(GSM6204110obj, features = "HHEX", label = T, reduction = "umap") # Marker associated with delta (δ) cells (somatostatin-secreting)


#CD4
FeaturePlot(GSM6204110obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD4", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IL7R", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CCR7", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "SELL", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "FOXP3", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IL2RA", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IL17A", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "RORC", label = T, reduction = "umap")#

#CD8
FeaturePlot(GSM6204110obj, features = "CD8A", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD8B", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "ITGA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "ZNF683", label = T, reduction = "umap")#
#Treg
FeaturePlot(GSM6204110obj, features = "FOXP3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "IL2RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "CTLA4", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "IKZF2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "TIGIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "IL10", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "AREG", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "IL1RL1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "ITGAE", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "CD69", label = TRUE, reduction = "umap")


#nk
FeaturePlot(GSM6204110obj, features = "NCAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "NKG7", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "GNLY", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "KLRD1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "EOMES", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "TBX21", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "FCGR3A", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "FGFBP2", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CX3CR1", label = T, reduction = "umap")#

#B cells
FeaturePlot(GSM6204110obj, features = "CD19", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "MS4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD79A", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD79B", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PAX5", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD22", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IGHM", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IGHG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "MZB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "JCHAIN", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PRDM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "XBP1", label = T, reduction = "umap")#

#Macrophage
FeaturePlot(GSM6204110obj, features = "CD163", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "MRC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "MSR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "ARG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CCL18", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CHI3L1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "VEGFA", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "TGFB1", label = T, reduction = "umap")#

#M1 Macrophages
FeaturePlot(GSM6204110obj, features = "CD80", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD86", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IL1B", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "TNF", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IL6", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "NOS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CXCL9", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CXCL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "STAT1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IRF5", label = T, reduction = "umap")#

#cDC
FeaturePlot(GSM6204110obj, features = "CLEC9A", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "XCR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "BATF3", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "IRF8", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CADM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "THBD", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "HLA-DQB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "HLA-DRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD74", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CLNK", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD207", label = T, reduction = "umap")#

#fibroblast
FeaturePlot(GSM6204110obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "DCN", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "LUM", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PDGFRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "VIM", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "FAP", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "FN1", label = T, reduction = "umap")#

#myofibroblast
FeaturePlot(GSM6204110obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "TAGLN", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CNN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "MYH11", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "FN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "POSTN", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PDGFRA", label = T, reduction = "umap")#

#Erythroblast umap
FeaturePlot(GSM6204110obj,features = "BPGM", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "TRIM58", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "XPO7", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "HBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "HBA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "ALAS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "GYPA", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "AHSP", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "KLF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "SLC4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "HEMGN", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "CA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "TFRC", label = T, reduction = "umap")#

#pDC
FeaturePlot(GSM6204110obj, features = "IL3RA", label = T, reduction = "umap") # CD123, canonical pDC marker
FeaturePlot(GSM6204110obj, features = "CLEC4C", label = T, reduction = "umap") # BDCA-2, highly specific to pDC
FeaturePlot(GSM6204110obj, features = "LILRA4", label = T, reduction = "umap") # ILT7, inhibitory receptor, pDC hallmark
FeaturePlot(GSM6204110obj, features = "TCF4", label = T, reduction = "umap")   # Lineage-defining TF
FeaturePlot(GSM6204110obj, features = "IRF7", label = T, reduction = "umap")   # High in pDC, type I IFN amplifier
FeaturePlot(GSM6204110obj, features = "TLR7", label = T, reduction = "umap")   # Viral RNA sensor
FeaturePlot(GSM6204110obj, features = "TLR9", label = T, reduction = "umap")   # CpG DNA sensor
FeaturePlot(GSM6204110obj, features = "BST2", label = T, reduction = "umap")   # IFN-inducible, en

#Endothelial umap
FeaturePlot(GSM6204110obj,features = "PECAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "VWF", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "CDH5", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "KDR", label = T, reduction = "umap")#

#Mast cell
FeaturePlot(GSM6204110obj,features = "HDC", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "HPGDS", label = T, reduction = "umap")#


#Neutrophil cell
FeaturePlot(GSM6204110obj,features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj,features = "CYBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "S100A8", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "S100A9", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "ELANE", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "AZU1", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "PRTN3", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CD177", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CXCR2", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "FCGR3B", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "CEACAM8", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "LTF", label = T, reduction = "umap")#
FeaturePlot(GSM6204110obj, features = "DEFA1", label = T, reduction = "umap")#
#CMP
FeaturePlot(GSM6204110obj, features = "CD34", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "KIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "SPI1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "CEBPA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "GATA2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "FLT3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "IL3RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "CSF1R", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "LY86", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204110obj, features = "KIT", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204110obj, features = "CD38", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "CD14", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204110obj, features = "KIT", label = TRUE, reduction = "umap") 

# Hepatocyte
FeaturePlot(GSM6204110obj, features = "ALB", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "ASGR1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "HNF4A", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "FABP1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "APOA1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "CYP3A4", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "GPC3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "AFP", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "KRT8", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204110obj, features = "KRT18", label = TRUE, reduction = "umap")

table(Idents(GSM6204110obj2))
GSM6204110obj2=GSM6204110obj
plotcdsevenGSM6204110obj2= FeaturePlot(GSM6204110obj2, features = "ASGR1", label = T, reduction = "umap")
HoverLocator(plot = plotcdsevenGSM6204110obj2, information = FetchData(GSM6204110obj2, vars = c("ident", "PC_1", "nFeature_RNA")))

#Endothelial assigning
selected.cellsAsEndothelialcell_CellGSM6204110obj2=CellSelector(plot = plotcdsevenGSM6204110obj2)
Idents(GSM6204110obj2, cells = selected.cellsAsEndothelialcell_CellGSM6204110obj2) <- "Endothelial"

#Endothelial assigning
selected.cellsAspDCcell_CellGSM6204110obj2=CellSelector(plot = plotcdsevenGSM6204110obj2)
Idents(GSM6204110obj2, cells = selected.cellsAspDCcell_CellGSM6204110obj2) <- "pDC"

#Erythroblast assigning
selected.cellsAsNKcell_CellGSM6204110obj2=CellSelector(plot = plotcdsevenGSM6204110obj2)
Idents(GSM6204110obj2, cells = selected.cellsAsNKcell_CellGSM6204110obj2) <- "Erythroblast"

#Fibroblast assigning
selected.cellsAsfibroblastcell_CellGSM6204110obj2=CellSelector(plot = plotcdsevenGSM6204110obj2)
length(selected.cellsAsfibroblastcell_CellGSM6204110obj2)
Idents(GSM6204110obj2, cells = selected.cellsAsfibroblastcell_CellGSM6204110obj2) <- "Fibroblast"

#NK Cells assigning
selected.cellsAsNKcell_CellGSM6204110obj2=CellSelector(plot = plotcdsevenGSM6204110obj2)
Idents(GSM6204110obj2, cells = selected.cellsAsNKcell_CellGSM6204110obj2) <- "NK Cells"

#Mast assigning
selected.cellsAsMastcell_CellGSM6204110obj2=CellSelector(plot = plotcdsevenGSM6204110obj2)
length(selected.cellsAsMastcell_CellGSM6204110obj2)
Idents(GSM6204110obj2, cells = selected.cellsAsMastcell_CellGSM6204110obj2) <- "Mast Cells"

#Hepatocyte assigning
selected.cellsAsHepatocytescell_CellGSM6204110obj2=CellSelector(plot = plotcdsevenGSM6204110obj2)

Idents(GSM6204110obj2, cells = selected.cellsAsHepatocytescell_CellGSM6204110obj2) <- "Hepatocytes"


fibrolikes=intersect(selected.cellsAsfibroblastcell_CellGSM6204110obj2,
                     selected.cellsAsEndothelialcell_CellGSM6204110obj2)
selected.cellsAsMastcell_CellGSM6204110obj2[1:2]
length(fibrolikes)
hepatolike=c("ACTTATCGTGACGTCC-1","CCGAACGTCTCCGCAT-1")
#=================================================
#=================================================


rownames( markers_of_GSM6204110[which(markers_of_GSM6204110$cluster==17 & markers_of_GSM6204110$avg_log2FC>9),])[1:25]
View( markers_of_GSM6204110[which(markers_of_GSM6204110$cluster==17),])


# Assays
# Get expression data for M1 cells
# Method 2: Create separate Seurat objects for M1 and M2

subsetMacrophages_GSM6204110=subset(GSM6204110obj2, idents = "Macrophages")

m1_obj_GSM6204110 <- subset(subsetMacrophages_GSM6204110, idents ="M1 Macrophage" )
m2_obj_GSM6204110 <- subset(subsetMacrophages_GSM6204110, idents = "M2 Macrophage" )
TransitionalM1_obj_GSM6204110<-subset(subsetMacrophages_GSM6204110,  idents = "Phenotype-Switching Macrophages")


m1_expressionGSM6204110 <- GetAssayData(m1_obj_GSM6204110, slot = "data", assay = "RNA")
m2_expressionGSM6204110 <- GetAssayData(m2_obj_GSM6204110, slot = "data", assay = "RNA")
TransitionalM1_expressionGSM6204110 <- GetAssayData(TransitionalM1_obj_GSM6204110, slot = "data",assay = "RNA")
c("M1 Macrophage", "M2 Macrophage","Phenotype-Switching Macrophages")
dim(m1_expressionGSM6204110)

m1_expression[1:5,1:5]
write.csv(m1_expressionGSM6204110,file = "C:/Users/abbas/Documents/GSE205013/GSM6204111/GSM6204111_RAW/Macrophages_M1_GSM6204110_data.csv" )
write.csv(m2_expressionGSM6204110,file = "C:/Users/abbas/Documents/GSE205013/GSM6204111/GSM6204111_RAW/Macrophages_M2_GSM6204110_data.csv"  )
write.csv(TransitionalM1_expressionGSM6204110,file = "C:/Users/abbas/Documents/GSE205013/GSM6204111/GSM6204111_RAW/TransitionalM1_expressionGSM6204110_data.csv")

phynotyeSwitchingMacrophagesOfGSM6204110=read.csv("C:/Users/abbas/Documents/GSE205013/GSM6204111/GSM6204111_RAW/Macrophages_M1_GSM6204110_data.csv")
length(phynotyeSwitchingMacrophagesOfGSM6204110)

