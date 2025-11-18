path6204113="C:/Users/abbas/Documents/GSE205013/GSM6204113/GSM6204113_RAW"


# Update seurat object
GSM6204113=Read10X(path6204113)
dim(GSM6204113)

#Step create srobj
GSM6204113obj=CreateSeuratObject(GSM6204113, min.cells = 3, min.features = 100)
dim(GSM6204113obj)

GSM6204113obj[["MTpercent"]]= PercentageFeatureSet(GSM6204113obj, pattern = "^MT-")
VlnPlot(GSM6204113obj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)
plot1= FeatureScatter(GSM6204113obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
GSM6204113obj= subset(GSM6204113obj, subset = nFeature_RNA >200 & nFeature_RNA<6000 & MTpercent <30)
plot1

# step Normalization
GSM6204113obj = NormalizeData(GSM6204113obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Step HighestVariableGenes
GSM6204113obj= FindVariableFeatures(GSM6204113obj, selection.method = "vst", nfeatures = 3000)
topVariableFeatures= head(VariableFeatures(GSM6204113obj), 80)
plot1= VariableFeaturePlot(GSM6204113obj)
plot2= LabelPoints(plot = plot1, points = topVariableFeatures, repel = T)
plot2

# Step Linear Reduction
GSM6204113obj= ScaleData(GSM6204113obj, features = rownames(GSM6204113obj))
GSM6204113obj= RunPCA(GSM6204113obj, features = VariableFeatures(GSM6204113obj))
print(GSM6204113obj[["pca"]], dims= 1:50, nfeatures= 5)
VizDimLoadings(GSM6204113obj, dims = 1:5 , reduction = "pca")
DimPlot(GSM6204113obj, reduction = "pca")
DimHeatmap(GSM6204113obj, dims = 1:6, cells = 100)
ElbowPlot(GSM6204113obj)

# Step Jackstraw and PCA
GSM6204113obj= JackStraw(GSM6204113obj, num.replicate = 100, dims = 50)
GSM6204113obj= ScoreJackStraw(GSM6204113obj, dims = 1:50)
JackStrawPlot(GSM6204113obj, dims = 1:50)
ElbowPlot(GSM6204113obj)
saveRDS(gsm8452584p1obj, file = "C:/Users/abbas/Documents/gsm8452584p1obj.rds")

# Step FindNeighbors
GSM6204113obj= FindNeighbors(GSM6204113obj, dims = 1:50)
GSM6204113obj= FindClusters(GSM6204113obj, resolution = 2)
length(levels(Idents(gsm6735860)))
table(Idents(GSM6204113obj))
# Step T-SNE
GSM6204113obj= RunTSNE(GSM6204113obj, dims = 1:13)
DimPlot(GSM6204113obj, reduction = "tsne", label = T)
# Step UMAP
GSM6204113obj= RunUMAP(GSM6204113obj, dims = 1:50)
DimPlot(GSM6204113obj, reduction = "umap",repel = T ,label = T )

GSM6204113obj2=GSM6204113obj
markers_of_GSM6204110 =FindAllMarkers(GSM6204113obj, min.pct = 0.25, logfc.threshold = 0.3)
table(markers_of_gsm452585$cluster)
doubeldmarkersofgsm452585= markers_of_gsm452585[which(markers_of_gsm452585$avg_log2FC>2),]
dim(doubeldmarkersofgsm452585)

# ASSIGNING Rename idents

linagesofPDAC=c("NK Cells","Mast Cells","Macrophages","Ductal Cells","CD4+ T Cells","CD8+ T Cells","trNK",
                "Hepatocytes","Fibroblasts","Endothelial", "B Cells", "pDC","Plasma Cells")

GSM6204113objcelltype= c("NK Cells","CD8+ T Cells", "trNK", "B Cells", "Plasma Cells",
                         "Ductal Cells","Ductal Cells","Ductal Cells","Ductal Cells","Ductal Cells","Ductal Cells",
                         "Ductal Cells","Ductal Cells","Ductal Cells","Ductal Cells","Ductal Cells","Ductal Cells",
                         "Ductal Cells","Ductal Cells","Ductal Cells","Macrophages","Ductal Cells","Ductal Cells",
                         "Ductal Cells", "Macrophages","CD8+ T Cells","Ductal Cells","Fibroblasts","CD4+ T Cells",
                         "Endothelial","Treg","Mast Cells","CD4+ T Cells","Muscle Cells","pDC",
                         "Ductal Cells","Acinar Cells")

newLevels6204112=c("Muscle Cells","NK Cells","trNK","Hepatocytes","pDC","Macrophages","CD8+ T Cells","Mast Cells",
                   "CD4+ T Cells","Ductal Cells","Treg",
                   "Endothelial","Acinar Cells","Fibroblasts","B Cells","Plasma Cells")

length(GSM6204113objcelltype)
length(newLevels6204112)
Idents(GSM6204113obj2) <- factor(Idents(GSM6204113obj2), levels = newLevels6204112)

length(GSM6204113objcelltype)
names(x= GSM6204113objcelltype)= levels(x= GSM6204113obj2)
GSM6204113obj2= RenameIdents(object = GSM6204113obj2, GSM6204113objcelltype)
DimPlot(GSM6204113obj2, reduction = "umap",repel = T ,label = T )
length(levels(x= GSM6204113obj2))
VlnPlot(GSM6204113obj2, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)
table(Idents(GSM6204113obj2))

# Dimplot Of Empty
DimPlot(GSM6204113obj2, reduction = "umap", label = T)

cluster16_markers <- FindMarkers(GSM6204113obj, ident.1 = 16)
cluster20_markers <- FindMarkers(GSM6204113obj, ident.1 = 20)
cluster0markers <- FindMarkers(GSM6204113obj, ident.1 = 0, ident.2 = c(10,13, 3))
cluster3markers <- FindMarkers(GSM6204113obj, ident.1 = 3, ident.2 = c(10, 13,0))

head(cluster13.markers)
View(cluster13.markers)
rownames(cluster16_markers[1:30,])
rownames(cluster20_markers[1:30,])
rownames(cluster0.markers[1:25,])
rownames(cluster3.markers[1:25,])
intersect_of_10_13=intersect(rownames(cluster13.markers[1:30,]),rownames(cluster10.markers[1:30,]))
length(intersect_of_10_13)

#ductal_markers 
FeaturePlot(GSM6204113obj, features = "KRT19", label = T, reduction = "umap") #
FeaturePlot(GSM6204113obj, features = "CDH1", label = T, reduction = "umap") #
FeaturePlot(GSM6204113obj, features = "CLDN4", label = T, reduction = "umap") #
FeaturePlot(GSM6204113obj, features = "AQP1", label = T, reduction = "umap") #
FeaturePlot(GSM6204113obj, features = "SOX9", label = T, reduction = "umap") #
FeaturePlot(GSM6204113obj, features = "PAX6", label = T, reduction = "umap") #
FeaturePlot(GSM6204113obj, features = "TP63", label = T, reduction = "umap") #
FeaturePlot(GSM6204113obj, features = "TP53", label = T, reduction = "umap") #
FeaturePlot(GSM6204113obj, features = "AKT2", label = T, reduction = "umap") #


#acinar cell markers
FeaturePlot(GSM6204113obj, features = "AMY2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CPA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CLPS", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PRSS1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CELA2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MUC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "KRT19", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PSMB5", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PLIN2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "NQO1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "RAB37", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "AQP8", label = T, reduction = "umap")#

#Endocrine cells
FeaturePlot(GSM6204113obj, features = "INS", label = T, reduction = "umap")  # Insulin, Beta cell marker
FeaturePlot(GSM6204113obj, features = "GCG", label = T, reduction = "umap")  # Glucagon, Alpha cell marker
FeaturePlot(GSM6204113obj, features = "SST", label = T, reduction = "umap")   # Somatostatin, Delta cell marker
FeaturePlot(GSM6204113obj, features = "PPY", label = T, reduction = "umap")  # Pancreatic Polypeptide, PP cell marker
FeaturePlot(GSM6204113obj, features = "NKX6-1", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204113obj, features = "MAFA", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204113obj, features = "PAX4", label = T, reduction = "umap")  # Beta cell development marker
FeaturePlot(GSM6204113obj, features = "CDH1", label = T, reduction = "umap")  # Cell adhesion, pan-endocrine marker

FeaturePlot(GSM6204113obj, features = "GHRL", label = T, reduction = "umap")  # Ghrelin gene, marker for epsilon (ε) cells
FeaturePlot(GSM6204113obj, features = "PDX1", label = T, reduction = "umap")  # Pancreatic and duodenal homeobox 1, critical for β cell function and development
FeaturePlot(GSM6204113obj, features = "ARX", label = T, reduction = "umap")  # Marker for alpha (α) cell identity (glucagon-secreting)
FeaturePlot(GSM6204113obj, features = "HHEX", label = T, reduction = "umap") # Marker associated with delta (δ) cells (somatostatin-secreting)


#CD4
FeaturePlot(GSM6204113obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD4", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IL7R", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CCR7", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "SELL", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "FOXP3", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IL2RA", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IL17A", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "RORC", label = T, reduction = "umap")#

#CD8
FeaturePlot(GSM6204113obj, features = "CD8A", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD8B", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ITGA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ZNF683", label = T, reduction = "umap")#
#Treg
FeaturePlot(GSM6204113obj, features = "FOXP3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "IL2RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "CTLA4", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "IKZF2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "TIGIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "IL10", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "AREG", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "IL1RL1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "ITGAE", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "CD69", label = TRUE, reduction = "umap")


#nk
FeaturePlot(GSM6204113obj, features = "NCAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "NKG7", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "GNLY", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "KLRD1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "EOMES", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "TBX21", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "FCGR3A", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "FGFBP2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CX3CR1", label = T, reduction = "umap")#

#B cells
FeaturePlot(GSM6204113obj, features = "CD19", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MS4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD79A", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD79B", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PAX5", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD22", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IGHM", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IGHG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MZB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "JCHAIN", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PRDM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "XBP1", label = T, reduction = "umap")#

#Macrophage
FeaturePlot(GSM6204113obj, features = "CD163", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MRC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MSR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ARG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CCL18", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CHI3L1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "VEGFA", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "TGFB1", label = T, reduction = "umap")#

#M1 Macrophages
FeaturePlot(GSM6204113obj, features = "CD80", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD86", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IL1B", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "TNF", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IL6", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "NOS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CXCL9", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CXCL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "STAT1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IRF5", label = T, reduction = "umap")#

#cDC
FeaturePlot(GSM6204113obj, features = "CLEC9A", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "XCR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "BATF3", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "IRF8", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CADM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "THBD", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "HLA-DQB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD207", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD74", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CLNK", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD14", label = T, reduction = "umap")#

#fibroblast
FeaturePlot(GSM6204113obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "DCN", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "LUM", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PDGFRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "VIM", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "FAP", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "FN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "SPP1", label = T, reduction = "umap")#

#myofibroblast
FeaturePlot(GSM6204113obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "TAGLN", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CNN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MYH11", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "FN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "POSTN", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PDGFRA", label = T, reduction = "umap")#

#Erythroblast umap
FeaturePlot(GSM6204113obj,features = "BPGM", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "TRIM58", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "XPO7", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "HBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "HBA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "ALAS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "GYPA", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "AHSP", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "KLF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "SLC4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "HEMGN", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "CA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "TFRC", label = T, reduction = "umap")#

#Muscle cell
FeaturePlot(GSM6204113obj, features = "LMOD1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MYLK", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CNN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "TAGLN", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MYH11", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "DES", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "TPM2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MYOCD", label = T, reduction = "umap")#

#pDC
FeaturePlot(GSM6204113obj, features = "IL3RA", label = T, reduction = "umap") # CD123, canonical pDC marker
FeaturePlot(GSM6204113obj, features = "CLEC4C", label = T, reduction = "umap") # BDCA-2, highly specific to pDC
FeaturePlot(GSM6204113obj, features = "LILRA4", label = T, reduction = "umap") # ILT7, inhibitory receptor, pDC hallmark
FeaturePlot(GSM6204113obj, features = "TCF4", label = T, reduction = "umap")   # Lineage-defining TF
FeaturePlot(GSM6204113obj, features = "IRF7", label = T, reduction = "umap")   # High in pDC, type I IFN amplifier
FeaturePlot(GSM6204113obj, features = "TLR7", label = T, reduction = "umap")   # Viral RNA sensor
FeaturePlot(GSM6204113obj, features = "TLR9", label = T, reduction = "umap")   # CpG DNA sensor
FeaturePlot(GSM6204113obj, features = "BST2", label = T, reduction = "umap")   # IFN-inducible, en

#Endothelial umap
FeaturePlot(GSM6204113obj,features = "PECAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "VWF", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "CDH5", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "KDR", label = T, reduction = "umap")#

#Mast cell
FeaturePlot(GSM6204113obj,features = "HDC", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "HPGDS", label = T, reduction = "umap")#

#Pericyte umap
FeaturePlot(GSM6204113obj,features = "PDGFRB", label = T, reduction = "umap")#24 7
FeaturePlot(GSM6204113obj,features = "CSPG4", label = T, reduction = "umap")#21 7
FeaturePlot(GSM6204113obj,features = "RGS5", label = T, reduction = "umap")#24
FeaturePlot(GSM6204113obj,features = "MCAM", label = T, reduction = "umap")#24
FeaturePlot(GSM6204113obj,features = "TAGLN", label = T, reduction = "umap")#24
#Adipocyte umap
FeaturePlot(GSM6204113obj,features = "PPARG", label = T, reduction = "umap")#12 21
FeaturePlot(GSM6204113obj,features = "CEBPA", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "ADIPOQ", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "APOE", label = T, reduction = "umap")#7 20 14
FeaturePlot(GSM6204113obj,features = "LPL", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "CFD", label = T, reduction = "umap")#7

#Neutrophil cell
FeaturePlot(GSM6204113obj,features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj,features = "CYBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "S100A8", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "S100A9", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "ELANE", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "AZU1", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "PRTN3", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CD177", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CXCR2", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "FCGR3B", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "CEACAM8", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "LTF", label = T, reduction = "umap")#
FeaturePlot(GSM6204113obj, features = "DEFA1", label = T, reduction = "umap")#
#CMP
FeaturePlot(GSM6204113obj, features = "CD34", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "KIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "SPI1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "CEBPA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "GATA2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "FLT3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "IL3RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "CSF1R", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "LY86", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204113obj, features = "KIT", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204113obj, features = "CD38", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "CD14", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204113obj, features = "KIT", label = TRUE, reduction = "umap") 

#Hepatocyte
FeaturePlot(GSM6204113obj, features = "ALB", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "ASGR1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "HNF4A", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "FABP1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "APOA1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "CYP3A4", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "GPC3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "AFP", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "KRT8", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204113obj, features = "KRT18", label = TRUE, reduction = "umap")

table(Idents(GSM6204113obj2))
GSM6204113obj2=GSM6204113obj
plotcdsevenGSM6204113obj2= FeaturePlot(GSM6204113obj, features = "MZB1", label = T, reduction = "umap")
HoverLocator(plot = plotcdsevenGSM6204113obj2, information = FetchData(GSM6204113obj2, vars = c("ident", "PC_1", "nFeature_RNA")))

#Endothelial assigning
selected.cellsAsEndothelialcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsEndothelialcell_CellGSM6204113obj2) <- "Endothelial"


#Ductal assigning
selected.cellsAsDuctalcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsDuctalcell_CellGSM6204113obj2) <- "Ductal Cells"


#Plasma Cells assigning
selected.cellsAsPlasmaCell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsPlasmaCell_CellGSM6204113obj2) <- "Plasma Cells"

#pDC assigning
selected.cellsAspDCcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAspDCcell_CellGSM6204113obj2) <- "pDC"

#Erythroblast assigning
selected.cellsAsNKcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsNKcell_CellGSM6204113obj2) <- "Erythroblast"

#Fibroblast assigning
selected.cellsAsfibroblastcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
length(selected.cellsAsfibroblastcell_CellGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsfibroblastcell_CellGSM6204113obj2) <- "Fibroblast"

#NK Cells assigning
selected.cellsAsNKcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsNKcell_CellGSM6204113obj2) <- "NK Cells"

#trNK Cells assigning
selected.cellsAsTrNKcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsTrNKcell_CellGSM6204113obj2) <- "trNK"

#Treg Cells assigning
selectedcellsAsTregCell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selectedcellsAsTregCell_CellGSM6204113obj2) <- "Treg"

#T Cells assigning
selectedcellsAsTregCell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
Idents(GSM6204113obj2, cells = selectedcellsAsTregCell_CellGSM6204113obj2) <- "CD8+ T Cells"

#Mast assigning
selected.cellsAsMastcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
length(selected.cellsAsMastcell_CellGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsMastcell_CellGSM6204113obj2) <- "Mast Cells"


#Hepatocyte assigning
selected.cellsAsHepatocytecell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
length(selected.cellsAsMastcell_CellGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsHepatocytecell_CellGSM6204113obj2) <- "Hepatocyte"
#B cell assigning
selected.cellsAsBcell_CellGSM6204113obj2=CellSelector(plot = plotcdsevenGSM6204113obj2)
length(selected.cellsAsBcell_CellGSM6204113obj2)
Idents(GSM6204113obj2, cells = selected.cellsAsBcell_CellGSM6204113obj2) <- "B Cells"

#Hepatocyte assigning
Idents(GSM6204113obj2, cells = plasmalike) <- "Plasma Cells"


fibrolikes=intersect(selected.cellsAsfibroblastcell_CellGSM6204113obj2,
                     selected.cellsAsEndothelialcell_CellGSM6204113obj2)
selected.cellsAsMastcell_CellGSM6204113obj2[1:2]
length(fibrolikes)
plasmalike=c("GACCCTTGTGATGGCA-1","AGGGTGAAGGTAGCAC-1","GACCCTTGTGATGGCA-1","CTGATCCGTAACTTCG-1")

#=================================================
#=================================================
#=================================================
#=================================================
cell_labels <- as.character(Idents(GSM6204113obj2))


counts_df <- as.data.frame(table(Idents(GSM6204113obj2)))
colnames(counts_df) <- c("celltype", "n")
counts_df <- counts_df %>%
  mutate(prop = n / sum(n),
         label = ifelse(n > 0, paste0(celltype, "\n", percent(prop)), paste0(celltype, "\n0%")))

# pie chart (ggplot2)
ggplot(counts_df, aes(x = "", y = prop, fill = celltype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Cell type composition (proportion)", x = NULL, y = NULL, fill = NULL) +
  theme_void() +
  theme(legend.position = "none")
cell = cell_labels
#=================================================


# --- 17-color vector (distinct hex colors) ---
cell_colors <- c(
  "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
  "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
  "#c7c7c7", "#bcbd22"
)

# ---- FUNCTIONAL: build proportions from a Seurat object (or metadata column) ----
# seurat_obj: your Seurat object
# id_col: optional string name of metadata column that stores identities.
#         If NULL (default), uses Idents(seurat_obj)
plot_celltype_pie <- function(GSM6204113obj2, id_col = NULL, title = "Cell type proportions") {
  # Get identity vector
  if (!is.null(id_col)) {
    id_vec <- as.character(GSM6204113obj2@meta.data[[id_col]])
  } else {
    id_vec <- as.character(Idents(GSM6204113obj2))
  }
  
  df_counts <- tibble(CellType = id_vec) %>%
    count(CellType, name = "n")
  
  # ensure all levels appear (even with zero counts)
  df_all <- tibble(CellType = newLevels6204111) %>%
    left_join(df_counts, by = "CellType") %>%
    replace_na(list(n = 0)) %>%
    mutate(
      CellType = factor(CellType, levels = newLevels6204111),
      frac = n / sum(n),
      pct = frac * 100,
      pct_label = sprintf("%.2f%%", pct)
    )
  
  # base pie chart via bar + coord_polar
  ggplot(df_all, aes(x = 1, y = frac, fill = CellType)) +
    geom_col(width = 1, color = "white", size = 0.2) +
    coord_polar(theta = "y") +
    # percent labels placed inside slices; check_overlap avoids collisions
    geom_text(aes(label = ifelse(n > 0, pct_label, "")), 
              position = position_stack(vjust = 0.5), size = 3, color = "white", check_overlap = TRUE) +
    scale_fill_manual(values = cell_colors, breaks = newLevels6204111) +
    theme_void() +
    labs(title = title, fill = "Cell type") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 9)
    )
}

# ---- Example usage ----
# Replace `seurat_obj` with your actual Seurat object variable.
# If you stored cell type calls in seurat_obj$celltype (example), call:
# plot_celltype_pie(seurat_obj, id_col = "celltype")
#
# Otherwise (default use Idents):
# plot_celltype_pie(seurat_obj)

# Example (uncomment and run with real object):
p <- plot_celltype_pie(GSM6204113obj2)
print(p)
#=================================================


rownames( markers_of_GSM6204110[which(markers_of_GSM6204110$cluster==17 & markers_of_GSM6204110$avg_log2FC>9),])[1:25]
View( markers_of_GSM6204110[which(markers_of_GSM6204110$cluster==17),])