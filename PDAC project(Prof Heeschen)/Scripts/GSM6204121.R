path21="C:/Users/abbas/Documents/GSE205013/GSM6204121/GSM6204121_RAW"

# Update seurat object
GSM6204121=Read10X(path21)
dim(GSM6204121)

GSM6204121obj2=readRDS("C:/Users/abbas/Documents/GSE205013/GSM6204121/GSM6204121obj2.rds")
saveRDS(GSM6204121obj2, file="C:/Users/abbas/Documents/GSE205013/GSM6204121/GSM6204121obj2.rds")

GSM6204114obj2@meta.data$mutation="TP53_mutated"

GSM6204121obj2=CreateSeuratObject(GSM6204121, min.cells = 3, min.features = 100)
GSM6204121obj2[["MTpercent"]]= PercentageFeatureSet(GSM6204121obj2, pattern = "^MT-")
rm(GSM6204121)
GSM6204121obj2= subset(GSM6204121obj2, subset = nFeature_RNA >200 & nFeature_RNA<8000 & MTpercent <30& nCount_RNA<50000)
dim(GSM6204121obj2)
GSM6204121obj2 <- NormalizeData(object = GSM6204121obj2, normalization.method = "LogNormalize", scale.factor = 10000)
GSM6204121obj2 <- FindVariableFeatures(object = GSM6204121obj2, selection.method = "vst", nfeatures = 3000)
GSM6204121obj2 <- ScaleData(object = GSM6204121obj2, features = VariableFeatures(GSM6204121obj2))
GSM6204121obj2 <- RunPCA(object = GSM6204121obj2, features = VariableFeatures(GSM6204121obj2))
GSM6204121obj2 <- FindNeighbors(object = GSM6204121obj2, dims = 1:30)
GSM6204121obj2 <- FindClusters(object = GSM6204121obj2)
GSM6204121obj2 <- RunUMAP(object = GSM6204121obj2, dims = 1:30)
DimPlot(object = GSM6204121obj2, reduction = "umap")


anotation_for_GSM6204121obj2=read.table("/Users/abbas/Documents/GSE205013/GSM6204121/anotaionGSM6204121.csv", sep = "," )
dim(anotation_for_GSM6204121obj2)
anotation_for_GSM6204121obj2[1:6,]
colnames(anotation_for_GSM6204121obj2)=anotation_for_GSM6204121obj2[1,]
anotation_for_GSM6204121obj2=anotation_for_GSM6204121obj2[-1,]
unique(anotation_for_GSM6204121obj2$ident)
length(anotatedAsFibroblastGSM6204121obj2)
anotatedAsMast_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Mast Cells"),2]
anotatedAsNK_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "NK Cells"),2]
anotatedAsNeutrophil_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Neutrophil Cells"),2]
anotatedAsTreg_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Treg"),2]
anotatedAsEndothelial_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Endothelial"),2]
anotatedAsB_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Plasma Cells"),2]
anotatedAsDuctalCells_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Ductal Cells"),2]
anotatedAsCD4_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "CD4+ T Cells"),2]
anotatedAsMac_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Macrophages"),2]
anotatedAsFibroblastGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Fibroblasts"),2]
anotatedAsCD8_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "CD8+ T Cells"),2]



length(anotatedAsPlasma_CellsGSM6204121obj2)
# Assigning Obj2
Idents(GSM6204121obj2, cells = anotatedAsMast_CellsGSM6204121obj2) <- "Mast Cells"
Idents(GSM6204121obj2, cells = anotatedAsNK_CellsGSM6204121obj2) <- "NK Cells"
Idents(GSM6204121obj2, cells = anotatedAsNeutrophil_CellsGSM6204121obj2) <- "Neutrophil Cells"
Idents(GSM6204121obj2, cells = anotatedAsTreg_CellsGSM6204121obj2) <- "Treg"
Idents(GSM6204121obj2, cells = anotatedAsEndothelial_CellsGSM6204121obj2) <- "Endothelial"
Idents(GSM6204121obj2, cells = anotatedAsB_CellsGSM6204121obj2) <- "B Cells"
Idents(GSM6204121obj2, cells = anotatedAsPlasma_CellsGSM6204121obj2) <- "Plasma Cells"
Idents(GSM6204121obj2, cells = anotatedAsDuctalCells_CellsGSM6204121obj2) <- "Ductal Cells"
Idents(GSM6204121obj2, cells = anotatedAsCD4_CellsGSM6204121obj2) <- "CD4+ T Cells"
Idents(GSM6204121obj2, cells = anotatedAsMac_CellsGSM6204121obj2) <- "Macrophages"
Idents(GSM6204121obj2, cells = anotatedAsFibroblastGSM6204121obj2) <- "Fibroblasts"
Idents(GSM6204121obj2, cells = anotatedAsCD8_CellsGSM6204121obj2) <- "CD8+ T Cells"

DimPlot(GSM6204121obj2, reduction = "umap",repel = T ,label = T )

table(Idents(GSM6204121obj2))

#Step create srobj
GSM6204121obj=CreateSeuratObject(GSM6204121, min.cells = 3, min.features = 100)
dim(GSM6204121obj)

GSM6204121obj[["MTpercent"]]= PercentageFeatureSet(GSM6204121obj, pattern = "^MT-")
VlnPlot(GSM6204121obj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)
plot1= FeatureScatter(GSM6204121obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
GSM6204121obj= subset(GSM6204121obj, subset = nFeature_RNA >200 & nFeature_RNA<8000 & MTpercent <30 & nCount_RNA<50000)#& nCount_RNA<50000
plot1

# step Normalization
GSM6204121obj = NormalizeData(GSM6204121obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Step HighestVariableGenes
GSM6204121obj= FindVariableFeatures(GSM6204121obj, selection.method = "vst", nfeatures = 3000)
topVariableFeatures= head(VariableFeatures(GSM6204121obj), 40)
plot1= VariableFeaturePlot(GSM6204121obj)
plot2= LabelPoints(plot = plot1, points = topVariableFeatures, repel = T)
plot2
GSM6204121obj2=GSM6204121obj
GSM6204121obj3= subset(GSM6204121obj3, subset = nFeature_RNA >200 & nFeature_RNA<7000 & MTpercent <30& nCount_RNA<50000)
dim(GSM6204121obj3)
GSM6204121obj3 <- NormalizeData(object = GSM6204121obj3)
GSM6204121obj3 <- FindVariableFeatures(object = GSM6204121obj3)
GSM6204121obj3 <- ScaleData(object = GSM6204121obj3)
GSM6204121obj3 <- RunPCA(object = GSM6204121obj3)
GSM6204121obj3 <- FindNeighbors(object = GSM6204121obj3, dims = 1:30)
GSM6204121obj3 <- FindClusters(object = GSM6204121obj3)
GSM6204121obj3 <- RunUMAP(object = GSM6204121obj3, dims = 1:30)
DimPlot(object = GSM6204121obj3, reduction = "umap")


# Step Linear Reduction
GSM6204121obj= ScaleData(GSM6204121obj, features = rownames(GSM6204121obj))
GSM6204121obj= RunPCA(GSM6204121obj, features = VariableFeatures(GSM6204121obj))
print(GSM6204121obj[["pca"]], dims= 1:50, nfeatures= 5)
VizDimLoadings(GSM6204121obj, dims = 1:5 , reduction = "pca")
DimPlot(GSM6204121obj, reduction = "pca")
DimHeatmap(GSM6204121obj, dims = 1:6, cells = 100)
ElbowPlot(GSM6204121obj)

# Step Jackstraw and PCA
GSM6204121obj= JackStraw(GSM6204121obj, num.replicate = 100, dims = 50)
GSM6204121obj= ScoreJackStraw(GSM6204121obj, dims = 1:50)
JackStrawPlot(GSM6204121obj, dims = 1:50)
ElbowPlot(GSM6204121obj)

# Step FindNeighbors
GSM6204121obj= FindNeighbors(GSM6204121obj, dims = 1:50)
GSM6204121obj= FindClusters(GSM6204121obj, resolution = 1)

length(levels(Idents(gsm6735860)))
table(Idents(GSM6204121obj))
# Step T-SNE
GSM6204121obj= RunTSNE(GSM6204121obj, dims = 1:50)
DimPlot(GSM6204121obj, reduction = "tsne", label = T)
# Step UMAP
GSM6204121obj= RunUMAP(GSM6204121obj, dims = 1:50)
DimPlot(GSM6204121obj, reduction = "umap",repel = T ,label = T )


markers_of_GSM6204121obj =FindAllMarkers(GSM6204121obj, min.pct = 0.25, logfc.threshold = 0.5)
summary(markers_of_gsm452585)
table(markers_of_gsm452585$cluster)
doubeldmarkersofgsm452585= markers_of_gsm452585[which(markers_of_gsm452585$avg_log2FC>2),]
dim(doubeldmarkersofgsm452585)
GSE452585p1t2ocelltype= LETTERS[1:13]
names(x= GSE452585p1t2ocelltype)= levels(x= GSM6204121obj)
GSM6204121obj2= RenameIdents(object = GSM6204121obj, GSE452585p1t2ocelltype)
DimPlot(GSM6204121obj2, reduction = "umap",repel = T ,label = T )
DimPlot(GSM6204121obj2, reduction = "tsne", label = T)


GSM6204121obj2=GSM6204121obj
anotation_for_GSM6204121obj2=read.table("C:/Users/abbas/Documents/GSE188711/GSM5688710/anotaionGSM5688710.csv", sep = "," )
dim(anotation_for_GSM6204121obj2)
colnames(anotation_for_GSM6204121obj2)=anotation_for_GSM6204121obj2[1,]
anotation_for_GSM6204121obj2=anotation_for_GSM6204121obj2[-1,]
length(unique(anotation_for_GSM6204121obj2$ident))
length(anotatedAsFibroblastGSM6204121obj2)
anotatedAsFibroblastGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Fibroblast"),2]
anotatedAsNeutrophil_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Neutrophil Cells"),2]
anotatedAsMacrophages_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Macrophages"),2]
anotatedAsEndothelial_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Endothelial"),2]
anotatedAsCD8_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "CD4+ T Cells"),2]
anotatedAsPlasma_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Plasma Cells"),2]
anotatedAsB_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "B Cells"),2]
anotatedAsTreg_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Treg"),2]
anotatedAsMast_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Mast Cells"),2]
anotatedAsNK_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "NK Cells"),2]
anotatedAsDC_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "DC"),2]
anotatedAsDuctal_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Ductal Cells"),2]
anotatedAstrNK_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "trNK"),2]
anotatedAsMonocytes_CellsGSM6204121obj2=anotation_for_GSM6204121obj2[which(anotation_for_GSM6204121obj2$ident == "Monocytes"),2]
length(anotatedAsPlasma_CellsGSM6204121obj2)

# Assigning Obj2
Idents(GSM6204121obj2, cells = anotatedAsNK_CellsGSM6204121obj2) <- "NK Cells"
Idents(GSM6204121obj2, cells = anotatedAstrNK_CellsGSM6204121obj2) <- "trNK"
Idents(GSM6204121obj2, cells = anotatedAsCD8_CellsGSM6204121obj2) <- "CD8+ T Cells"
Idents(GSM6204121obj2, cells = anotatedAsCD4_CellsGSM6204121obj2) <- "CD4+ T Cells"
Idents(GSM6204121obj2, cells = anotatedAsTreg_CellsGSM6204121obj2) <- "Treg"
Idents(GSM6204121obj2, cells = anotatedAsNeutrophil_CellsGSM6204121obj2) <- "Neutrophil Cells"
Idents(GSM6204121obj2, cells = anotatedAsPlasma_CellsGSM6204121obj2) <- "Plasma Cells"
Idents(GSM6204121obj2, cells = anotatedAsB_CellsGSM6204121obj2) <- "B Cells"
Idents(GSM6204121obj2, cells = anotatedAsMacrophages_CellsGSM6204121obj2) <- "Macrophages"
Idents(GSM6204121obj2, cells = anotatedAsFibroblastGSM6204121obj2) <- "Fibroblast"
Idents(GSM6204121obj2, cells = anotatedAsEndothelial_CellsGSM6204121obj2) <- "Endothelial"
Idents(GSM6204121obj2, cells = anotatedAsDC_CellsGSM6204121obj2) <- "DC"
Idents(GSM6204121obj2, cells = anotatedAsDuctal_CellsGSM6204121obj2) <- "Ductal Cells"
Idents(GSM6204121obj2, cells = anotatedAsMast_CellsGSM6204121obj2) <- "Mast Cells"
Idents(GSM6204121obj2, cells = anotatedAsMonocytes_CellsGSM6204121obj2) <- "Monocytes"



markers_of_GSM6204121 =FindAllMarkers(GSM6204121obj2, min.pct = 0.25, logfc.threshold = 0.3)
table(markers_of_gsm452585$cluster)
doubelGSM6204121clusterAmarker= markers_of_GSM6204121[which(markers_of_GSM6204121$avg_log2FC>1.4 
                                                            & markers_of_GSM6204121$cluster=="A"),]
doubelGSM6204121clusterBmarker= markers_of_GSM6204121[which(markers_of_GSM6204121$avg_log2FC>1.4 
                                                            & markers_of_GSM6204121$cluster=="B"),]
doubelGSM6204121clusterCmarker= markers_of_GSM6204121[which(markers_of_GSM6204121$avg_log2FC>1.4 
                                                            & markers_of_GSM6204121$cluster=="C"),]
doubelGSM6204121clusterDmarker= markers_of_GSM6204121[which(markers_of_GSM6204121$avg_log2FC>1.4 
                                                            & markers_of_GSM6204121$cluster=="D"),]
doubelGSM6204121clusterEmarker= markers_of_GSM6204121[which(markers_of_GSM6204121$avg_log2FC>1.4 
                                                            & markers_of_GSM6204121$cluster=="E"),]
doubelGSM6204121clusterFmarker= markers_of_GSM6204121[which(markers_of_GSM6204121$avg_log2FC>1.4 
                                                            & markers_of_GSM6204121$cluster=="F"),]



dim(doubelGSM6204121clusterDmarker)
View(doubelGSM6204121clusterBmarker)
# ASSIGNING Rename idents

linagesofPDAC=c("NK Cells","Mast Cells","Macrophages","Ductal Cells","CD4+ T Cells","CD8+ T Cells","trNK",
                "Hepatocytes","Fibroblasts","Endothelial", "B Cells", "pDC","Plasma Cells")

GSM6204121objcelltype= c("Mast Cells","NK Cells","Neutrophil Cells","Treg","Endothelial","B Cells","Plasma Cells",
                         "Ductal Cells","CD4+ T Cells",
                         "Macrophages","Fibroblasts",
                          "CD8+ T Cells")
                         c("Ductal Cells","Ductal Cells","Macrophages","Ductal Cells","Ductal Cells","Ductal Cells",
                         "CD4+ T Cells","NK Cells","CD8+ T Cells","Ductal Cells","Ductal Cells","Macrophages",
                         "Fibroblasts","CD8+ T Cells","Neutrophil Cells","Ductal Cells","Macrophages","Ductal Cells","Ductal Cells")

GSM6204121objcelltype= c("Plasma Cells","Mast Cells","Macrophages","Endothelial",
                         "Hepatocytes","Fibroblasts","B Cells","NK Cells","Ductal Cells","CD4+ T Cells",
                         "CD8+ T Cells","Endocrine Cells","Endothelial","Treg","Ductal Cells")


length(GSM6204121objcelltype)
newLevelsGSM6204121=c("NK Cells","trNK","CD8+ T Cells","CD4+ T Cells","Treg",
                      "B Cells","Plasma Cells","","Macrophages","DC","cDC","Monocytes",
                      "Ductal Cells","Acinar Cells","Endothelial","Fibroblasts","","","")

length(newLevelsGSM6204121)
Idents(GSM6204121obj2) <- factor(Idents(GSM6204121obj2), levels = newLevelsGSM6204121)

length(GSM6204121objcelltype)
names(x= GSM6204121objcelltype)= levels(x= GSM6204121obj2)
GSM6204121obj2= RenameIdents(object = GSM6204121obj2, GSM6204121objcelltype)
DimPlot(GSM6204121obj2, reduction = "umap",repel = T ,label = T )
table(Idents(GSM6204121obj2))

length(levels(x= GSM6204121obj2))

# Dimplot Of Empty
DimPlot(GSM6204121obj2, reduction = "umap", label = T, repel = T)

clusterA.markersGSM6204121obj2 <- FindMarkers(GSM6204121obj2, ident.1 = "A")
clusterB.markersGSM6204121obj2 <- FindMarkers(GSM6204121obj2, ident.1 = "B")
clusterC_markersGSM6204121obj2 <- FindMarkers(GSM6204121obj2, ident.1 = "C")
clusterE.markersGSM6204121obj2 <- FindMarkers(GSM6204121obj2, ident.1 = "E")
clusterD.markersGSM6204121obj2 <- FindMarkers(GSM6204121obj2, ident.1 = "D")

write.csv(doubelGSM6204121clusterAmarker, file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/clusterA.markersGSM6204121.csv")
write.csv(doubelGSM6204121clusterBmarker, file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/clusterB.markersGSM6204121.csv")
write.csv(doubelGSM6204121clusterCmarker, file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/clusterC.markersGSM6204121.csv")
write.csv(doubelGSM6204121clusterDmarker, file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/clusterD.markersGSM6204121.csv")
write.csv(doubelGSM6204121clusterEmarker, file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/clusterE.markersGSM6204121.csv")
write.csv(doubelGSM6204121clusterFmarker, file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/clusterF.markersGSM6204121.csv")
View(clusterB.markersGSM6204121obj2)
head(cluster11.markers)
View(cluster11.markers)
rownames(cluster17.markers[1:30,])
rownames(cluster10.markers[1:30,])
rownames(cluster0.markers[1:25,])
rownames(cluster3.markers[1:25,])
intersect_of_10_13=intersect(rownames(cluster13.markers[1:30,]),rownames(cluster10.markers[1:30,]))
length(intersect_of_10_13)
VlnPlot(GSM6204121obj2, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)

DoHeatmap(GSM6204121obj2,features =  topOfGSM6204112Heatmap, angle = 60)
newLevelsGSM6204121=c("Ductal Cells","Fibroblasts","Endothelial","Macrophages","CD8+ T Cells","CD4+ T Cells",
                      "NK Cells","trNK","Muscle Cells","Hepatocytes","pDC","Mast Cells","Treg","B Cells","Plasma Cells",
                      "Acinar Cells","Endocrine Cells")
length(newLevelsGSM620411)
multixIdents = c("Ductal Cells", "Macrophages", "Fibroblasts","Endothelial","CD8+ T Cells","CD4+ T Cells", "Treg")
allTcellsIdents=c("Treg","CD8+ T Cells","CD4+ T Cells" )
subsetMacrophages_GSM6204121=subset(GSM6204121obj2, idents = "Macrophages")

Idents(GSM6204121obj2) <- factor(Idents(GSM6204121obj2), levels = newLevelsGSM620411)
GSM6204121obj3=GSM6204121obj2

#Cell Cycling 
GSM6204121obj3 <- CellCycleScoring(GSM6204121obj3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
RidgePlot(GSM6204121obj3, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
GSM6204121obj3 <- RunPCA(GSM6204121obj3, features = c(s.genes, g2m.genes))



FeaturePlot(GSM6204121obj2, features = "MTpercent", label = T)
FeaturePlot(GSM6204121obj3, features = g2m.genes)
FeaturePlot(object = GSM6204121obj3, features = "G2M.Score", label = F)
FeaturePlot(object = GSM6204121obj3, features = "S.Score")
table(Idents(GSM6204121obj2))
# view cell cycle scores and phase assignments
#ductal_markers 
FeaturePlot(GSM6204121obj, features = "KRT19", label = T, reduction = "umap") #
FeaturePlot(GSM6204121obj, features = "CDH1", label = T, reduction = "umap") #
FeaturePlot(GSM6204121obj, features = "CLDN4", label = T, reduction = "umap") #
FeaturePlot(GSM6204121obj, features = "AQP1", label = T, reduction = "umap") #
FeaturePlot(GSM6204121obj, features = "SOX9", label = T, reduction = "umap") #
FeaturePlot(GSM6204121obj, features = "PAX6", label = T, reduction = "umap") #
FeaturePlot(GSM6204121obj, features = "TP63", label = T, reduction = "umap") #
FeaturePlot(GSM6204121obj, features = "TP53", label = T, reduction = "umap") #
FeaturePlot(GSM6204121obj, features = "AKT2", label = T, reduction = "umap") #

#acinar cell markers
FeaturePlot(GSM6204121obj, features = "AMY2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CPA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CLPS", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PRSS1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CELA2A", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MUC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "KRT19", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PSMB5", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PLIN2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "NQO1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "RAB37", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "AQP8", label = T, reduction = "umap")#

#Endocrine cells
FeaturePlot(GSM6204121obj, features = "INS", label = T, reduction = "umap")  # Insulin, Beta cell marker
FeaturePlot(GSM6204121obj, features = "GCG", label = T, reduction = "umap")  # Glucagon, Alpha cell marker
FeaturePlot(GSM6204121obj, features = "SST", label = T, reduction = "umap")   # Somatostatin, Delta cell marker
FeaturePlot(GSM6204121obj, features = "PPY", label = T, reduction = "umap")  # Pancreatic Polypeptide, PP cell marker
FeaturePlot(GSM6204121obj, features = "NKX6-1", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204121obj, features = "MAFA", label = T, reduction = "umap")  # Beta cell transcription factor
FeaturePlot(GSM6204121obj, features = "PAX4", label = T, reduction = "umap")  # Beta cell development marker
FeaturePlot(GSM6204121obj, features = "CDH1", label = T, reduction = "umap")  # Cell adhesion, pan-endocrine marker
FeaturePlot(GSM6204121obj, features = "CHGB", label = T, reduction = "umap")#

FeaturePlot(GSM6204121obj, features = "GHRL", label = T, reduction = "umap")  # Ghrelin gene, marker for epsilon (ε) cells
FeaturePlot(GSM6204121obj, features = "PDX1", label = T, reduction = "umap")  # Pancreatic and duodenal homeobox 1, critical for β cell function and development
FeaturePlot(GSM6204121obj, features = "ARX", label = T, reduction = "umap")  # Marker for alpha (α) cell identity (glucagon-secreting)
FeaturePlot(GSM6204121obj, features = "HHEX", label = T, reduction = "umap") # Marker associated with delta (δ) cells (somatostatin-secreting)


#CD4
FeaturePlot(GSM6204121obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD4", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IL7R", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CCR7", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "SELL", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "FOXP3", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IL2RA", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IL17A", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "RORC", label = T, reduction = "umap")#

#CD8
FeaturePlot(GSM6204121obj, features = "CD8A", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD8B", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ITGA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IFNG", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ZNF683", label = T, reduction = "umap")#

#Treg
FeaturePlot(GSM6204121obj, features = "FOXP3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "IL2RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "CTLA4", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "IKZF2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "TIGIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "IL10", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "AREG", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "IL1RL1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "ITGAE", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "CD69", label = TRUE, reduction = "umap")


#nk
FeaturePlot(GSM6204121obj, features = "GNLY", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "EOMES", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "NKG7", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "NCAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "FCGR3A", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "GZMB", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PRF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "KLRD1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CXCR6", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ITGAE", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "TBX21", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "FGFBP2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CX3CR1", label = T, reduction = "umap")#

#B cells
FeaturePlot(GSM6204121obj, features = "CD19", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MS4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD79A", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD79B", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PAX5", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD22", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IGHM", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IGHG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MZB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "JCHAIN", label = F, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PRDM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "XBP1", label = T, reduction = "umap")#

#Macrophage
FeaturePlot(GSM6204121obj, features = "FCN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "C1QA", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "C1QB", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CHI3L1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "LYZ", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD14", label = T, reduction = "umap")#

#Macrophage M2
FeaturePlot(GSM6204121obj, features = "CD163", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MRC1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MSR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ARG1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CCL18", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CHI3L1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "VEGFA", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "TGFB1", label = T, reduction = "umap")#

#M1 Macrophages
FeaturePlot(GSM6204121obj, features = "CD80", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD86", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IL1B", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "TNF", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IL6", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "NOS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CXCL9", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CXCL10", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "STAT1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IRF5", label = T, reduction = "umap")#

#cDC
FeaturePlot(GSM6204121obj, features = "CLEC9A", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "XCR1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "BATF3", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "IRF8", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CADM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "THBD", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "HLA-DQB1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "HLA-DRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD74", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CLNK", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD207", label = T, reduction = "umap")#


#Erythroblast umap
FeaturePlot(GSM6204121obj,features = "BPGM", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "TRIM58", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "XPO7", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "HBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "HBA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "ALAS2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "GYPA", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "AHSP", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "KLF1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "SLC4A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "HEMGN", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "CA1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "TFRC", label = T, reduction = "umap")#

#pDC
FeaturePlot(GSM6204121obj, features = "IL3RA", label = T, reduction = "umap") # CD123, canonical pDC marker
FeaturePlot(GSM6204121obj, features = "CLEC4C", label = T, reduction = "umap") # BDCA-2, highly specific to pDC
FeaturePlot(GSM6204121obj, features = "LILRA4", label = T, reduction = "umap") # ILT7, inhibitory receptor, pDC hallmark
FeaturePlot(GSM6204121obj, features = "TCF4", label = T, reduction = "umap")   # Lineage-defining TF
FeaturePlot(GSM6204121obj, features = "IRF7", label = T, reduction = "umap")   # High in pDC, type I IFN amplifier
FeaturePlot(GSM6204121obj, features = "TLR7", label = T, reduction = "umap")   # Viral RNA sensor
FeaturePlot(GSM6204121obj, features = "TLR9", label = T, reduction = "umap")   # CpG DNA sensor
FeaturePlot(GSM6204121obj, features = "BST2", label = T, reduction = "umap")   # IFN-inducible, en

#Mast cell
FeaturePlot(GSM6204121obj,features = "HDC", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "HPGDS", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "TPSB2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "TPSAB1", label = T, reduction = "umap")#


#Neutrophil cell
FeaturePlot(GSM6204121obj,features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "CYBB", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "S100A8", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "S100A9", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ELANE", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "AZU1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PRTN3", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CD177", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CXCR2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "FCGR3B", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CEACAM8", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "LTF", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "DEFA1", label = T, reduction = "umap")#

#Fibroblast
FeaturePlot(GSM6204121obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "DCN", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "LUM", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PDGFRA", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "VIM", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "FAP", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "FN1", label = T, reduction = "umap")#

#myofibroblast
FeaturePlot(GSM6204121obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "TAGLN", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CNN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MYH11", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "FN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "COL3A1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "POSTN", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "THY1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "PDGFRA", label = T, reduction = "umap")#

#Muscle cell
FeaturePlot(GSM6204121obj, features = "MYH11", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MYOCD", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "DES", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "LMOD1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "ACTA2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "TAGLN", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "CNN1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "TPM2", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj, features = "MYLK", label = T, reduction = "umap")#

#Adipocyte umap
FeaturePlot(GSM6204121obj,features = "PPARG", label = T, reduction = "umap")#12 21
FeaturePlot(GSM6204121obj,features = "CEBPA", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "ADIPOQ", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "APOE", label = T, reduction = "umap")#7 20 14
FeaturePlot(GSM6204121obj,features = "LPL", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "CFD", label = T, reduction = "umap")#7

#Endothelial umap
FeaturePlot(GSM6204121obj,features = "PECAM1", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "VWF", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "CDH5", label = T, reduction = "umap")#
FeaturePlot(GSM6204121obj,features = "KDR", label = T, reduction = "umap")#


#CMP
FeaturePlot(GSM6204121obj, features = "CD34", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "KIT", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "SPI1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "CEBPA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "GATA2", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "FLT3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "IL3RA", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "CSF1R", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "LY86", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204121obj, features = "KIT", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204121obj, features = "CD38", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "CD14", label = TRUE, reduction = "umap")  
FeaturePlot(GSM6204121obj, features = "KIT", label = TRUE, reduction = "umap") 

# Hepatocyte
FeaturePlot(GSM6204121obj, features = "ALB", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "ASGR1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "HNF4A", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "FABP1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "APOA1", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "CYP3A4", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "GPC3", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "AFP", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "KRT8", label = TRUE, reduction = "umap")
FeaturePlot(GSM6204121obj, features = "KRT18", label = TRUE, reduction = "umap")

table(Idents(GSM6204121obj2))
GSM6204121obj2=GSM6204121obj
FeaturePlot(GSM6204121obj, features = "GZMB", label = T, reduction = "tsne")

plotcdsevenGSM6204121obj2= FeaturePlot(GSM6204121obj, features = "GZMB", label = T, reduction = "umap")
HoverLocator(plot = plotcdsevenGSM6204121obj2, information = FetchData(GSM6204121obj2, vars = c("ident", "PC_1", "nFeature_RNA")))

#Endothelial assigning
selected.cellsAsEndothelialcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsEndothelialcell_CellGSM6204121obj2) <- "Endothelial"

#pDC assigning
selected.cellsAspDCcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAspDCcell_CellGSM6204121obj2) <- "DC"

#Macrophages assigning
selected.cellsAsMacrophagescell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsMacrophagescell_CellGSM6204121obj2) <- "Macrophages"

#Treg assigning
selected.cellsAsTregcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsTregcell_CellGSM6204121obj2) <- "Treg"

#plasma assigning
selected.cellsAsPlasmacell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsPlasmacell_CellGSM6204121obj2) <- "Plasma Cells"

#plasma assigning
selected.cellsAsBcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsBcell_CellGSM6204121obj2) <- "B Cells"
#plasma assigning
selected.cellsAsAcinarcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsAcinarcell_CellGSM6204121obj2) <- "Acinar Cells"

#Ductal assigning
selected.cellsAsDuctalcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsDuctalcell_CellGSM6204121obj2) <- "Ductal Cells"
#Erythroblast assigning
selected.cellsAsMusclecell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsMusclecell_CellGSM6204121obj2) <- "Muscle Cells"

#Fibroblast assigning
selected.cellsAsfibroblastcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
length(selected.cellsAsfibroblastcell_CellGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsfibroblastcell_CellGSM6204121obj2) <- "Fibroblast"

#NK Cells assigning
selected.cellsAsNKcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsNKcell_CellGSM6204121obj2) <- "NK Cells"

#Mast assigning
selected.cellsAsMastcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
length(selected.cellsAsMastcell_CellGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsMastcell_CellGSM6204121obj2) <- "Mast Cells"


#Hepatocytes assigning
selected.cellsAsHepatocytesCell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
length(selected.cellsAsHepatocytesCell_CellGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsHepatocytesCell_CellGSM6204121obj2) <- "Hepatocytes"

#trNK Cells assigning
selected.cellsAsTrNKcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsTrNKcell_CellGSM6204121obj2) <- "trNK"

#cd4 Cells assigning
selected.cellsAscd4cell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAscd4cell_CellGSM6204121obj2) <- "CD4+ T Cells"

#cd8 Cells assigning
selected.cellsAscd8cell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAscd8cell_CellGSM6204121obj2) <- "CD8+ T Cells"

#pDC assigning
selected.cellsacell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsacell_CellGSM6204121obj2) <- "A"

#pDC assigning
selected.cellsbcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsbcell_CellGSM6204121obj2) <- "B"

#pDC assigning
selected.cellsAsCcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsCcell_CellGSM6204121obj2) <- "C"
#pDC assigning
selected.cellsAsDcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsDcell_CellGSM6204121obj2) <- "D"
#pDC assigning
selected.cellsAsEcell_CellGSM6204121obj2=CellSelector(plot = plotcdsevenGSM6204121obj2)
Idents(GSM6204121obj2, cells = selected.cellsAsEcell_CellGSM6204121obj2) <- "E"
#Hepatocyte assigning
Idents(GSM6204121obj2, cells = hepatolike) <- "Hepatocytes"


fibrolikes=intersect(selected.cellsAsfibroblastcell_CellGSM6204121obj2,
                     selected.cellsAsEndothelialcell_CellGSM6204121obj2)
selected.cellsAsMastcell_CellGSM6204121obj2[1:2]
length(fibrolikes)
hepatolike=c("ACTTATCGTGACGTCC-1","CCGAACGTCTCCGCAT-1")
#=================================================
#=================================================
DoHeatmap(multix_obj_GSM6204121, features = differenHeatmapVector)
DimPlot(ductal_obj_GSM6204121, )
DotPlot(multix_obj_GSM6204121, features = pdac_proliferative_genes) + RotatedAxis()
DotPlot(multix_obj_GSM6204121, features =c( pdac_tumor_suppressor_genes,pdac_steady_genes)) + RotatedAxis()

RidgePlot(multix_obj_GSM6204121, features = pdac_proliferative_genes, ncol = 2)

cell_labels <- as.character(Idents(GSM6204121obj2))


counts_df <- as.data.frame(table(Idents(GSM6204121obj2)))
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
plot_celltype_pie <- function(GSM6204121obj2, id_col = NULL, title = "Cell type proportions") {
  # Get identity vector
  if (!is.null(id_col)) {
    id_vec <- as.character(GSM6204121obj2@meta.data[[id_col]])
  } else {
    id_vec <- as.character(Idents(GSM6204121obj2))
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
p <- plot_celltype_pie(GSM6204121obj2)
print(p)
#=================================================
rownames( markers_of_GSM6204121[which(markers_of_GSM6204121$cluster==17 & markers_of_GSM6204121$avg_log2FC>9),])[1:25]
View( markers_of_GSM6204121[which(markers_of_GSM6204121$cluster==17),])

macrophage_markers6204111=c("IL1B", "CD86", "CLEC10A", "IL23A","CXCL5","TLR2","TLR4","STAT1",
                            "CD163", "CD68", "MMP9", "MSR1", "MRC1","CCL18", "CSF1R")
m1featureForFeaturePlot=c("IL1B", "CD86", "CLEC10A", "IL23A","CXCL5","TLR2","TLR4","STAT1", "IRF4")
m2featureForFeaturePlotA=c("CD163", "CD68", "MMP9", "MSR1")
m2featureForFeaturePlotB=c("MRC1","CCL18", "CSF1R","SIGLEC1")#"IL10")
anti_apoptosis_module <- c("BCL2","MCL1","BCL2L1","BIRC5","BIRC2","XIAP")#,"BIRC3"
death_receptor_axis <- c("FAS","FASLG","TNFRSF10B","CASP8","BIRC3","CASP3","BID","CFLAR")
FeaturePlot(multix_obj_GSM6204121, features = anti_apoptosis_module, label = F, reduction = "umap")#

FeaturePlot(subsetMacrophages_GSM6204121, features = exhaustion_core_module , label = F, reduction = "umap")
FeaturePlot(subsetMacrophages_GSM6204121, features = m2featureForFeaturePlotB , label = F, reduction = "umap")
FeaturePlot(subsetMacrophages_GSM6204121, features = "ENO2", label = F, reduction = "umap")#
FeaturePlot(tcells_obj_GSM6204121, features = "ENO2", label = F, reduction = "umap")#

#=================================================
#=================================================
DotPlot(multix_obj_GSM6204121, features = eno1_module)+ RotatedAxis()
RidgePlot(multix_obj_GSM6204121, features = pilot_features)

#=================================================
#=================================================
tgfb_treg_module <- c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","FOXP3","IL2RA","CTLA4")
# H) Exhaustion / checkpoint readout (primarily on T cells)
exhaustion_core_module <- c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX","NR4A1","NR4A2","NR4A3")
# I) Phagocytosis checkpoint (“don’t-eat-me”)
cd47_sirpa_module <- c("CD47","SIRPA")
cholesterol_handling_module <- c("SREBF2","SOAT1","HMGCR","LDLR","MSMO1")
gcn2_isr_module <- c("EIF2AK4","ATF4","ASNS","DDIT3","PPP1R15A")
er_stress_xbp1_module <- c("XBP1","ERN1","ATF6","HSPA5","HSP90B1","DDIT3")
suppressive_secretome_module <- c("TGFB1","IL10","CSF1","CCL2","CCL22","LGALS9","VSIR" )
adenosine_module <- c("ENTPD1","NT5E","ADORA2A","ADORA2B","ADA")
ahr_kyn_module <- c(
  "IDO1","TDO2","IL4I1",           # ligand production
  "AHR","AHRR","CYP1A1","CYP1B1"   # sensor/targets
)
eno1_module <- c(
  "ENO1", "ENO2", "ALDOA", "GAPDH", "PGK1", "PGAM1", "TPI1",
  "PKM", "HK2", "PFKP", "PFKFB3", "PFKFB4", "LDHA",
  "SLC2A1", "SLC16A1", "SLC16A3",        # GLUT1, MCT1, MCT4 (transporters)
  "HIF1A", "MYC", "STAT3", "SP1", "E2F1", "RELA", # transcriptional regulators
  "PDK1", "PDK3",                          # modulate glycolysis/PDH
  "ANXA2", "PLAU", "PLAUR", "PLG", "SERPINE1",  # plasminogen activation / invasion
  "YBX1", "HSPA5"                           # chaperone/regulator links (stress/glycolysis)
)

heatmapENOFeature=c(eno1_module, ahr_kyn_module, suppressive_secretome_module,adenosine_module, gcn2_isr_module,
                    cholesterol_handling_module,cd47_sirpa_module,exhaustion_core_module,tgfb_treg_module)
DoHeatmap(multix_obj_GSM6204121, features = heatmapENOFeature)

#=================================================
#=================================================
#=================================================
pilot_features= c("IL1A", "NLRP3", "PYCARD", "IRF1", "MYD88",
                  "IL6ST", "JAK1", "STAT3", "SOCS3", "NFKB2", "ADAM17", "CEBPB",
                  "TNF", "LITAF", "ADAM17", "TLR4", "NFKB1", "NFKBIA", "RELA", "MAP3K7", "MAPK14", "FOS", "JUN", "CREB1",
                  "CXCL8")
length(pilot_features)
immunosuppressive_TME_genes <- c(
  # CAF (Cancer-Associated Fibroblast) markers:
  "ACTA2",    # alpha-SMA, general CAF marker
  "FAP",      # fibroblast activation protein
  "LRRC15",   # myCAF marker
  "PDPN",     # iCAF marker
  "POSTN",    # periostin, myCAF
  "DCN",      # decorin, general/ECM CAF
  "COL1A1",   # collagen type I
  "MMP2",     # matrix metallopeptidase 2
  "TAGLN",    # myCAF marker
  "PDGFRB",   # general CAF marker
  "SPARC",    # secreted protein acidic and rich in cysteine
  
  # Treg (Regulatory T Cell) markers:
  "FOXP3",    # canonical Treg marker
  "CTLA4",    # immune checkpoint upregulated in Tregs
  "IL2RA",    # CD25, Treg activation marker
  "ENTPD1",   # CD39, immunosuppressive Treg
  "TNFRSF4"  # OX40, Treg activation marker
  #"LAYN"      # exhaustion/Treg marker
)
immusFeaturePlotPolot=c("FAP","MMP2","CTLA4","TNFRSF4")



# -------------------------
# Aging & senescence pathway vectors (Up / Down)
# -------------------------

# 1) Telomere attrition / replicative exhaustion
telomere_attrition_up   <- c("CDKN2A","CDKN1A","TP53","TP53BP1","TERF1","TERF2")
telomere_attrition_down <- c("TERT")

# 2) DNA damage response (genotoxic stress, chemo/radiation)
dna_damage_up   <- c("ATM","ATR","CHEK1","CHEK2","TP53","TP53BP1","H2AFX","GADD45A")
dna_damage_down <- c("MKI67","PCNA")

# 3) Oncogene-induced senescence (OIS)
ois_up   <- c("CDKN2A","CDKN1A","IL1A","NFKB1","RELA")
ois_down <- c("MKI67","PCNA")  # proliferation genes typically down in OIS

# 4) Mitochondrial dysfunction / metabolic stress
mito_dysfunction_up   <- c("SOD2","GPX1","HMOX1","NQO1","PINK1","PARK2")
mito_dysfunction_down <- c("PPARGC1A","NDUFS1","COX4I1","SDHA","ATP5F1A") # OXPHOS/biogenesis markers (context-dependent)

# 5) Epigenetic alterations (chromatin remodel / methylation drift)
epigenetic_up   <- c("HMGA1","HMGA2","EZH2","SUZ12","EED")  # remodeling / PRC components (context-specific)
epigenetic_down <- c("DNMT1")  # regulator often altered (global/focal changes depend on context)

# 6) Proteostasis collapse / UPR and autophagy decline
proteostasis_up   <- c("HSPA5","DDIT3","HSP90AA1","HSP90AB1","UBB","PSMB5")
proteostasis_down <- c("BECN1","ATG5","ATG7")  # autophagy regulators (variable)

# 7) Oxidative stress & chronic inflammation (inflammaging / SASP)
inflammaging_up <- c("IL6","CXCL8","IL1A","IL1B","CCL2","CXCL1","CXCL2","CXCL3",
                     "MMP1","MMP3","SERPINE1")

# 8) Anti-apoptotic shift (senescent cell apoptosis resistance)
anti_apoptotic_up   <- c("BCL2","BCL2L1","BCL2L2","MCL1","XIAP","AKT1","PIK3CA")
anti_apoptotic_down <- c("BAX","BBC3")  # BBC3 = PUMA

# 1) Telomere attrition / replicative exhaustion
# 2) DNA damage response
# 3) Oncogene-induced senescence (OIS)
# 4) Mitochondrial dysfunction / metabolic stress
# 5) Epigenetic alterations (chromatin remodel / methylation drift)
# 6) Proteostasis collapse / UPR
# 7) Oxidative stress & chronic inflammation (inflammaging / SASP)
# 8) Anti-apoptotic shift


# -------------------------
# Pack into a single list for easier programmatic use
# -------------------------
aging_senescence_modules <- list(
  Telomere = list(up = telomere_attrition_up, down = telomere_attrition_down),
  DNA_Damage = list(up = dna_damage_up, down = dna_damage_down),
  OIS = list(up = ois_up, down = ois_down),
  Mitochondria = list(up = mito_dysfunction_up, down = mito_dysfunction_down),
  Epigenetic = list(up = epigenetic_up, down = epigenetic_down),
  Proteostasis = list(up = proteostasis_up, down = proteostasis_down),
  Inflammaging = list(up = inflammaging_up, down = character(0)),
  AntiApoptotic = list(up = anti_apoptotic_up, down = anti_apoptotic_down)
)

# -------------------------
# Flattened aggregated vectors (if you want a single "up" or "down" module)
# -------------------------
all_upregulated   <- unique(c(telomere_attrition_up, dna_damage_up, ois_up,
                              mito_dysfunction_up, epigenetic_up, proteostasis_up,
                              inflammaging_up, anti_apoptotic_up))

all_downregulated <- unique(c(telomere_attrition_down, dna_damage_down, ois_down,
                              mito_dysfunction_down, epigenetic_down, proteostasis_down,
                              anti_apoptotic_down))

alteredFeature=c(all_upregulated, all_downregulated)

tcellFeaturesOfActivation=c("CD3E", "CD8B", "GZMB", "PDCD1")
FeaturePlot(tcells_obj_GSM6204121, features = tcellFeaturesOfActivation, reduction = "umap")
xhs

DoHeatmap(multix_obj_GSM6204121, features = all_upregulated)
DoHeatmap(multix_obj_GSM6204121, features = all_downregulated)
DotPlot(multix_obj_GSM6204121, features = alteredFeature) +RotatedAxis()

RidgePlot(multix_obj_GSM6204121,features =  immunosuppressive_TME_genes)
DotPlot(multix_obj_GSM6204121,features =  immunosuppressive_TME_genes) +RotatedAxis()
VlnPlot(multix_obj_GSM6204121, features =pdac_proliferative_genes)
DoHeatmap(multix_obj_GSM6204121, features =immunosuppressive_TME_genes)
FeaturePlot(multix_obj_GSM6204121, features = immusFeaturePlotPolot, reduction = "umap")








#=================================================
#=================================================
#=================================================
# Assays
# Get expression data for M1 cells
# Method 2: Create separate Seurat objects for M1 and M2






dim(subsetMacrophages_GSM6204121)
m1_obj_GSM6204121 <- subset(subsetMacrophages_GSM6204121, idents ="M1 Macrophage" )
m2_obj_GSM6204121 <- subset(subsetMacrophages_GSM6204121, idents = "M2 Macrophage" )
TransitionalM1_obj_GSM6204121<-subset(subsetMacrophages_GSM6204121,  idents = "Phenotype-Switching Macrophages")


macrophages_expressionGSM6204121 <- GetAssayData(subsetMacrophages_GSM6204121, slot = "data", assay = "RNA")
m2_expressionGSM6204121 <- GetAssayData(m2_obj_GSM6204121, slot = "data", assay = "RNA")
TransitionalM1_expressionGSM6204121 <- GetAssayData(TransitionalM1_obj_GSM6204121, slot = "data",assay = "RNA")
ductal_obj_GSM6204121 <- subset(GSM6204121obj2, idents = "Ductal Cells" )
multix_obj_GSM6204121 <- subset(GSM6204121obj2, idents = multixIdents )
tcells_obj_GSM6204121 <- subset(GSM6204121obj2, idents = allTcellsIdents )
newLevelsGSM6204121=newLevels6204117
expressionIdentsGSM6204121obj2_list <- list()
dim(ductal_obj_GSM6204121)
newLevelsGSM6204121=c(newLevelsGSM6204121, "Treg")
length(newLevelsGSM6204121)
xhsGSM6204121=data.frame()
expressionIdentsGSMGSM6204121obj2_list <- list()
newLevelsGSM6204121=unique(GSM6204121objcelltype)
for(celltype in newLevelsGSM6204121){
  subset_obj <- subset(GSM6204121obj2, idents = celltype)
  expr_data <- GetAssayData(subset_obj, slot = "data", assay = "RNA")
  expressionIdentsGSMGSM6204121obj2_list[[celltype]] <- expr_data
  xhsGSM6204121=rbind(xhsGSM6204121,data.frame(barcodes=colnames(expr_data),ident=celltype))
}

write.csv(xhsGSM6204121, file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/anotaionGSM6204121.csv")
dim(xhsGSM6204121)

differenHeatmapVector=c(pdac_proliferative_genes, pdac_tumor_suppressor_genes,pdac_steady_genes)
dim(m1_expressionGSM6204121)
pdac_proliferative_genes <- c(
  "KRAS",      # Oncogene frequently mutated in PDAC
  "MYC",       # Regulates cell cycle & proliferation
  "CCND1",     # Cyclin D1, cell cycle progression
  "CDK4",      # Cell cycle dependent kinase
  "PCNA",      # Proliferating cell nuclear antigen
  "MKI67",     # Marker of proliferation
  "EGFR",      # Epidermal growth factor receptor
  "FGFR1",     # Fibroblast growth factor receptor 1
  "IGFL2"      # Promotes cell survival and proliferation
)
pdac_tumor_suppressor_genes <- c(
  "TP53",      # Most commonly mutated tumor suppressor
  "CDKN2A",    # Encodes p16, regulates cell cycle arrest
  "SMAD4",     # TGF-beta signaling, lost in many PDAC
  "PTEN",      # PI3K pathway inhibitor, often deleted
  "WWOX"       # Loss accelerates PDAC progression
)
pdac_steady_genes <- c(
  "ACTB", "GAPDH", "TUBB", "B2M", "UBC", "RPLP0"
)


#=================================================
#=================================================
newLevelsGSM6204121=c("NK Cells","trNK","CD8+ T Cells","CD4+ T Cells","Treg",
                      "Macrophages","DC","Mast Cells","B Cells","Plasma Cells",
                      "Ductal Cells",
                      "Acinar Cells","Endothelial","Fibroblasts")
xhsGSM6204121=data.frame()
expressionIdentsGSMGSM6204121obj2_list <- list()
dim(ductal_obj_GSM5688710)
for(celltype in newLevelsGSM6204121){
  subset_obj <- subset(GSM6204121obj2, idents = celltype)
  expr_data <- GetAssayData(subset_obj, slot = "data", assay = "RNA")
  expressionIdentsGSMGSM6204121obj2_list[[celltype]] <- expr_data
  xhsGSM6204121=rbind(xhsGSM6204121,data.frame(barcodes=colnames(expr_data),ident=celltype))
}
write.csv(xhsGSM6204121, file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/anotaionGSM6204121.csv")
dim(xhsGSM6204121)


expressionIdentsGSM6204121obj2_list$
  m1_expression[1:5,1:5]

write.csv(m1_expressionGSM6204121,file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/GSM6204121_RAW/Macrophages_M1_GSM6204121_data.csv" )
write.csv(m2_expressionGSM6204121,file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/GSM6204121_RAW/Macrophages_M2_GSM6204121_data.csv")
write.csv(TransitionalM1_expressionGSM6204121,file = "C:/Users/abbas/Documents/GSE205013/GSM6204121/GSM6204121_RAW/TransitionalM1_expressionGSM6204121_data.csv")

phynotyeSwitchingMacrophagesOfGSM6204121=read.csv("C:/Users/abbas/Documents/GSE205013/GSM6204121/GSM6204121_RAW/Macrophages_M1_GSM6204121_data.csv")
length(phynotyeSwitchingMacrophagesOfGSM6204121)
