pathGSM6204109="C:/Users/abbas/Documents/GSE205013/GSM6204109/GSM6204109_RAW"
install.packages('devtools')
library(scales)
library(ggplot2)

library(ggrepel)
# Update seurat object
GSM6204109=Read10X(pathGSM6204109)
dim(GSM6204109)
#Step create srobj
GSM6204109obj=CreateSeuratObject(GSM6204109, min.cells = 3, min.features = 100)
dim(GSM6204109obj)

GSM6204109obj[["MTpercent"]]= PercentageFeatureSet(GSM6204109obj, pattern = "^MT-")
VlnPlot(GSM6204109obj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)
plot1= FeatureScatter(GSM6204109obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
GSM6204109obj= subset(GSM6204109obj, subset = nFeature_RNA >200 & nFeature_RNA<7500 & MTpercent <35)
plot1

# step Normalization
GSM6204109obj = NormalizeData(GSM6204109obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Step HighestVariableGenes
GSM6204109obj= FindVariableFeatures(GSM6204109obj, selection.method = "vst", nfeatures = 3000)
topVariableFeatures= head(VariableFeatures(GSM6204109obj), 40)
plot1= VariableFeaturePlot(GSM6204109obj)
plot2= LabelPoints(plot = plot1, points = topVariableFeatures, repel = T)
plot2

# Step Linear Reduction
GSM6204109obj= ScaleData(GSM6204109obj, features = rownames(GSM6204109obj))
GSM6204109obj= RunPCA(GSM6204109obj, features = VariableFeatures(GSM6204109obj))
print(GSM6204109obj[["pca"]], dims= 1:50, nfeatures= 5)
VizDimLoadings(GSM6204109obj, dims = 1:5 , reduction = "pca")
DimPlot(GSM6204109obj, reduction = "pca")
DimHeatmap(GSM6204109obj, dims = 1:6, cells = 100)
ElbowPlot(GSM6204109obj)

# Step Jackstraw and PCA
GSM6204109obj= JackStraw(GSM6204109obj, num.replicate = 100, dims = 50)
GSM6204109obj= ScoreJackStraw(GSM6204109obj, dims = 1:50)
JackStrawPlot(GSM6204109obj, dims = 1:50)
ElbowPlot(GSM6204109obj)

# Step FindNeighbors
GSM6204109obj= FindNeighbors(GSM6204109obj, dims = 1:50)
GSM6204109obj= FindClusters(GSM6204109obj, resolution = 0.8)
length(levels(Idents(gsm6735860)))
table(Idents(GSM6204109obj2))
# Step T-SNE
GSM6204109obj= RunTSNE(GSM6204109obj, dims = 1:13)
DimPlot(GSM6204109obj, reduction = "tsne", label = T)
# Step UMAP
GSM6204109obj= RunUMAP(GSM6204109obj, dims = 1:50)
DimPlot(GSM6204109obj, reduction = "umap",repel = T ,label = T )
View(markers_of_GSM6204109[which(markers_of_GSM6204109$cluster==17),])

GSM6204109obj2=GSM6204109obj

markers_of_GSM6204109 =FindAllMarkers(GSM6204109obj, min.pct = 0.25, logfc.threshold = 0.3)
table(markers_of_gsm452585$cluster)
doubeldmarkersofgsm452585= markers_of_gsm452585[which(markers_of_gsm452585$avg_log2FC>2),]
dim(doubeldmarkersofgsm452585)

# ASSIGNING Rename idents
GSM6204109objcelltype= c("Macrophages","pDC","Mast Cells","NK Cells",
                         "Macrophages", "Ductal Cells", "Ductal Cells", "CD4+ T Cells", "Ductal Cells","CD8+ T Cells",
                          "Macrophages", "Ductal Cells", "Ductal Cells", "Macrophages" , "trNK", "CD4+ T Cells",
                          "Endothelial", "trNK","CD4+ T Cells", "Macrophages","B Cells", "Hepatocytes",
                          "Fibroblasts", "CD4+ T Cells")

linagesofPDAC=c("NK Cells","Mast Cells","Macrophages","Ductal Cells","CD4+ T Cells","CD8+ T Cells",
                "trNK","Hepatocytes","Fibroblasts","Endothelial", "B Cells", "pDC")
length(GSM6204109objcelltype)
names(x= GSM6204109objcelltype)= levels(x= GSM6204109obj2)
GSM6204109obj2= RenameIdents(object = GSM6204109obj2, GSM6204109objcelltype)
DimPlot(GSM6204109obj2, reduction = "umap",repel = T ,label = T )

# Dimplot Of Empty
DimPlot(GSM6204109obj2, reduction = "umap", label = T)

cluster13.markers <- FindMarkers(GSM6204109obj, ident.1 = 13, ident.2 = c(0, 3))
cluster10_markers <- FindMarkers(GSM6204109obj, ident.1 = 10, ident.2 = c(0, 3,13))
cluster0.markers <- FindMarkers(GSM6204109obj, ident.1 = 0, ident.2 = c(10,13, 3))
cluster3.markers <- FindMarkers(GSM6204109obj, ident.1 = 3, ident.2 = c(10, 13,0))

head(cluster13.markers)
View(cluster13.markers)
rownames(cluster13.markers[1:30,])
rownames(cluster10.markers[1:30,])
rownames(cluster0.markers[1:25,])
rownames(cluster3.markers[1:25,])
intersect_of_10_13=intersect(rownames(cluster13.markers[1:30,]),rownames(cluster10.markers[1:30,]))
length(intersect_of_10_13)

#ductal_markers 
FeaturePlot(GSM6204109obj, features = "KRT19", label = T, reduction = "umap") #
FeaturePlot(GSM6204109obj, features = "CDH1", label = T, reduction = "umap") #
FeaturePlot(GSM6204109obj, features = "CLDN4", label = T, reduction = "umap") #
FeaturePlot(GSM6204109obj, features = "AQP1", label = T, reduction = "umap") #
FeaturePlot(GSM6204109obj, features = "SOX9", label = T, reduction = "umap") #
FeaturePlot(GSM6204109obj, features = "PAX6", label = T, reduction = "umap") #
FeaturePlot(GSM6204109obj, features = "TP63", label = T, reduction = "umap") #
FeaturePlot(GSM6204109obj, features = "TP53", label = T, reduction = "umap") #
FeaturePlot(GSM6204109obj, features = "AKT2", label = T, reduction = "umap") #

#pDC
FeaturePlot(GSM6204109obj, features = "IL3RA", label = T, reduction = "umap") # CD123, canonical pDC marker
FeaturePlot(GSM6204109obj, features = "CLEC4C", label = T, reduction = "umap") # BDCA-2, highly specific to pDC
FeaturePlot(GSM6204109obj, features = "LILRA4", label = T, reduction = "umap") # ILT7, inhibitory receptor, pDC hallmark
FeaturePlot(GSM6204109obj, features = "TCF4", label = T, reduction = "umap")   # Lineage-defining TF
FeaturePlot(GSM6204109obj, features = "IRF7", label = T, reduction = "umap")   # High in pDC, type I IFN amplifier
FeaturePlot(GSM6204109obj, features = "TLR7", label = T, reduction = "umap")   # Viral RNA sensor
FeaturePlot(GSM6204109obj, features = "TLR9", label = T, reduction = "umap")   # CpG DNA sensor
FeaturePlot(GSM6204109obj, features = "BST2", label = T, reduction = "umap")   # IFN-inducible, en


#acinar cell markers
FeaturePlot(GSM6204109obj, features = "AMY2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CPA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CLPS", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PRSS1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CELA2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "MUC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "KRT19", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PSMB5", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PLIN2", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "NQO1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "RAB37", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "AQP8", label = T, reduction = "umap")#

#Endocrine cells
FeaturePlot(GSM6204109obj, features = "INS", label = T, reduction = "umap")  # Insulin, Beta cell marker
FeaturePlot(GSM6204109obj, features = "GCG", label = T, reduction = "umap")  # Glucagon, Alpha cell marker
FeaturePlot(GSM6204109obj, features = "SST", label = T, reduction = "umap")   # Somatostatin, Delta cell marker
FeaturePlot(GSM6204109obj, features = "PPY", label = T, reduction = "umap")  # Pancreatic Polypeptide, PP cell marker
FeaturePlot(GSM6204109obj, features = "NKX6-1", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204109obj, features = "MAFA", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204109obj, features = "PAX4", label = T, reduction = "umap")  # Beta cell development marker
FeaturePlot(GSM6204109obj, features = "CDH1", label = T, reduction = "umap")  # Cell adhesion, pan-endocrine marker

FeaturePlot(GSM6204109obj, features = "GHRL", label = T, reduction = "umap")  # Ghrelin gene, marker for epsilon (ε) cells
FeaturePlot(GSM6204109obj, features = "PDX1", label = T, reduction = "umap")  # Pancreatic and duodenal homeobox 1, critical for β cell function and development
FeaturePlot(GSM6204109obj, features = "ARX", label = T, reduction = "umap")  # Marker for alpha (α) cell identity (glucagon-secreting)
FeaturePlot(GSM6204109obj, features = "HHEX", label = T, reduction = "umap") # Marker associated with delta (δ) cells (somatostatin-secreting)


#CD4
FeaturePlot(GSM6204109obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD4", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IL7R", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CCR7", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "SELL", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "FOXP3", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IL2RA", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IL17A", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "RORC", label = T, reduction = "umap")#

#CD8
FeaturePlot(GSM6204109obj, features = "CD8A", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD8B", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "ITGA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "ZNF683", label = T, reduction = "umap")#


#nk
FeaturePlot(GSM6204109obj, features = "NCAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "NKG7", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "GNLY", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "KLRD1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "EOMES", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "TBX21", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "FCGR3A", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "FGFBP2", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CX3CR1", label = T, reduction = "umap")#

#B cells
FeaturePlot(GSM6204109obj, features = "CD19", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "MS4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD79A", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD79B", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PAX5", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD22", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IGHM", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IGHG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "MZB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "JCHAIN", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PRDM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "XBP1", label = T, reduction = "umap")#

#Macrophage
FeaturePlot(GSM6204109obj, features = "CD163", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "MRC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "MSR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "ARG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CCL18", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CHI3L1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "VEGFA", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "TGFB1", label = T, reduction = "umap")#

plottii=FeaturePlot(GSM6204109obj, features =c("IL1R2", "MMP9", "CLEC10A","CSF1R", "TGFB1"), label = T, reduction = "umap")#
plottii

#M1 Macrophages
FeaturePlot(GSM6204109obj, features = "CD80", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD86", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IL1B", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "TNF", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IL6", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "NOS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CXCL9", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CXCL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "STAT1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IRF5", label = T, reduction = "umap")#

#cDC
FeaturePlot(GSM6204109obj, features = "CLEC9A", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "XCR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "BATF3", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IRF8", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CADM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "THBD", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "HLA-DQB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "HLA-DRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD74", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CLNK", label = T, reduction = "umap")#

#fibroblast
FeaturePlot(GSM6204109obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "DCN", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "LUM", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PDGFRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "VIM", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "FAP", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "FN1", label = T, reduction = "umap")#

#myofibroblast
FeaturePlot(GSM6204109obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "TAGLN", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CNN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "MYH11", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "FN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "POSTN", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PDGFRA", label = T, reduction = "umap")#


#CMP
FeaturePlot(GSM6204109obj, features = "CD34", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "KIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "SPI1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "CEBPA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "GATA2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "FLT3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "IL3RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "CSF1R", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "LY86", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204109obj, features = "KIT", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204109obj, features = "CD38", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204109obj, features = "CD14", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204109obj, features = "KIT", label = TRUE, reduction = "umap") 

#Erythroblast umap
FeaturePlot(GSM6204109obj,features = "BPGM", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "TRIM58", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "XPO7", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "HBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "HBA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "ALAS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "GYPA", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "AHSP", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "KLF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "SLC4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "HEMGN", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "CA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "TFRC", label = T, reduction = "umap")#

#cDC1
FeaturePlot(GSM6204109obj, features = "CLEC9A", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "XCR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "BATF3", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "IRF8", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CADM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "THBD", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "HLA-DQB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "HLA-DRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD74", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CLNK", label = T, reduction = "umap")#
#Endothelial umap
FeaturePlot(GSM6204109obj,features = "PECAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "VWF", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "CDH5", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "KDR", label = T, reduction = "umap")#
#Mast cell
FeaturePlot(GSM6204109obj,features = "HDC", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "HPGDS", label = T, reduction = "umap")#
#Neutrophil cell
FeaturePlot(GSM6204109obj,features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj,features = "CYBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "S100A8", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "S100A9", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "ELANE", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "AZU1", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "PRTN3", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CD177", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CXCR2", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "FCGR3B", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "CEACAM8", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "LTF", label = T, reduction = "umap")#
FeaturePlot(GSM6204109obj, features = "DEFA1", label = T, reduction = "umap")#
#pDC
FeaturePlot(GSM6204109obj, features = "IL3RA", label = T, reduction = "umap") # CD123, canonical pDC marker
FeaturePlot(GSM6204109obj, features = "CLEC4C", label = T, reduction = "umap") # BDCA-2, highly specific to pDC
FeaturePlot(GSM6204109obj, features = "LILRA4", label = T, reduction = "umap") # ILT7, inhibitory receptor, pDC hallmark
FeaturePlot(GSM6204109obj, features = "TCF4", label = T, reduction = "umap")   # Lineage-defining TF
FeaturePlot(GSM6204109obj, features = "IRF7", label = T, reduction = "umap")   # High in pDC, type I IFN amplifier
FeaturePlot(GSM6204109obj, features = "TLR7", label = T, reduction = "umap")   # Viral RNA sensor
FeaturePlot(GSM6204109obj, features = "TLR9", label = T, reduction = "umap")   # CpG DNA sensor
FeaturePlot(GSM6204109obj, features = "BST2", label = T, reduction = "umap")   # IFN-inducible, en

GSM6204109obj2=GSM6204109obj
plotcdsevenGSM6204109obj2= FeaturePlot(GSM6204109obj2, features = "HDC", label = T, reduction = "umap")

#Erythroblast assigning
selected.cellsAsNKcell_CellGSM6204109obj2=CellSelector(plot = plotcdsevenGSM6204109obj2)
Idents(GSM6204109obj2, cells = selected.cellsAsNKcell_CellGSM6204109obj2) <- "Erythroblast"

#NK assigning
selected.cellsAsNKcell_CellGSM6204109obj2=CellSelector(plot = plotcdsevenGSM6204109obj2)
Idents(GSM6204109obj2, cells = selected.cellsAsNKcell_CellGSM6204109obj2) <- "NK Cells"


#trNK Cells assigning
selected.cellsAsTrNKcell_CellGSM6204109obj2=CellSelector(plot = plotcdsevenGSM6204109obj2)
Idents(GSM6204109obj2, cells = selected.cellsAsTrNKcell_CellGSM6204109obj2) <- "NK Cells"

#Mast assigning
selected.cellsAsMastcell_CellGSM6204109obj2=CellSelector(plot = plotcdsevenGSM6204109obj2)
length(selected.cellsAsMastcell_CellGSM6204109obj2)
Idents(GSM6204109obj2, cells = selected.cellsAsMastcell_CellGSM6204109obj2) <- "Mast Cells"

#Macrophages assigning
selected.cellsAsMacrophagesCell_CellGSM6204109obj2=CellSelector(plot = plotcdsevenGSM6204109obj2)
Idents(GSM6204109obj2, cells = selected.cellsAsMacrophagesCell_CellGSM6204109obj2) <- "Macrophages"


#Plasma Cells assigning
selected.cellsAsPlasmaCell_CellGSM6204109obj2=CellSelector(plot = plotcdsevenGSM6204109obj2)
Idents(GSM6204109obj2, cells = selected.cellsAsPlasmaCell_CellGSM6204109obj2) <- "Plasma Cells"

#Erythroblast assigning
selected.cellsAspDCcell_CellGSM6204109obj2=CellSelector(plot = plotcdsevenGSM6204109obj2)
length(selected.cellsAspDCcell_CellGSM6204109obj2)
Idents(GSM6204109obj2, cells = selected.cellsAspDCcell_CellGSM6204109obj2) <- "pDC"

#=================================================
#=================================================


rownames( markers_of_GSM6204109[which(markers_of_GSM6204109$cluster==17 & markers_of_GSM6204109$avg_log2FC>9),])[1:25]
View( markers_of_GSM6204109[which(markers_of_GSM6204109$cluster==17),])
