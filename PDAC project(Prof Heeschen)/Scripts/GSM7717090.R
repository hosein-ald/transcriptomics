require(ggplot2)
require(Seurat)
require(reshape2)

GSM7717090_RAW=Read10X("C:/Users/abbas/Documents/GSE241132/GSM7717090/GSM7717090_PWH28D30/PWH28D30/")
dim(GSM7717090_RAW)
GSM7717090_RAW[1:5,1:5]


#read obj file
GSM7717090_obj2=readRDS("C:/Users/abbas/Documents/GSE241132/GSM7717090/GSM7717090_obj.rds")
saveRDS(GSM7717090_obj2, file="C:/Users/abbas/Documents/GSE241132/GSM7717090/GSM7717090_obj.rds")


GSM7717090_obj2= CreateSeuratObject(GSM7717090_RAW, min.cells =3, min.features = 100 )
GSM7717090_obj2[["MTpercent"]]=PercentageFeatureSet(GSM7717090_obj2, pattern = "^MT-")
GSM7717090_obj2=subset(GSM7717090_obj2, subset = nFeature_RNA<7500 &  MTpercent<25 & nCount_RNA<90000)
GSM7717090_obj2= NormalizeData(GSM7717090_obj2, normalization.method = "LogNormalize", scale.factor = 10000)
GSM7717090_obj2=FindVariableFeatures(GSM7717090_obj2,selection.method = "vst", nfeatures = 2000)
GSM7717090_obj2=ScaleData(GSM7717090_obj2, features = rownames(GSM7717090_obj2))
GSM7717090_obj2=RunPCA(GSM7717090_obj2,features = VariableFeatures(GSM7717090_obj2))
GSM7717090_obj2=RunUMAP(GSM7717090_obj2, dims = 1:50)
DimPlot(GSM7717090_obj2,reduction = "umap", label.size =4 ,repel = T, label = T)


anotation_for_GSM7717090_obj2=read.table("C:/Users/abbas/Documents/GSE241132/GSM7717090/anotaionGSM7717090.csv", sep = "," )
dim(anotation_for_GSM7717090_obj2)
anotation_for_GSM7717090_obj2[1:6,]
colnames(anotation_for_GSM7717090_obj2)=anotation_for_GSM7717090_obj2[1,]
anotation_for_GSM7717090_obj2=anotation_for_GSM7717090_obj2[-1,]
length(unique(anotation_for_GSM7717090_obj2$ident))
length(anotatedAsFibroblastGSM7717090_obj2)
anotatedAsFibroblastGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Fibroblasts"),2]
anotatedAsKeratinocytes_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Keratinocytes"),2]
anotatedAsNK_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "NK Cells"),2]
anotatedAsLangerhans_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Langerhans"),2]
anotatedAsCD8_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "CD4+ T Cells"),2]
anotatedAscDC_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "cDC"),2]
anotatedAsm1mac_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "M1 Macrophages"),2]
anotatedAsm2Macrophages_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "M2 Macrophages"),2]
anotatedAsm3_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Phenotype-Switching Macrophages"),2]
anotatedAsMast_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Mast Cells"),2]
anotatedAsB_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "B Cells"),2]
anotatedAsEndothelial_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Endothelial"),2]
anotatedAsPericytes_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Pericytes"),2]
anotatedAsMelanocytes_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Melanocytes"),2]
anotatedAsPlasma_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Plasma Cells"),2]
anotatedAsTreg_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Treg"),2]
anotatedAsNeutrophil_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Neuronal Cells"),2]

anotatedAsEpithelial_CellsGSM7717090_obj2=anotation_for_GSM7717090_obj2[which(anotation_for_GSM7717090_obj2$ident == "Epithelial"),2]
length(anotatedAsPlasma_CellsGSM7717090_obj2)
# Assigning Obj2
Idents(GSM7717090_obj2, cells = anotatedAsFibroblastGSM7717090_obj2) <- "Fibroblasts"
Idents(GSM7717090_obj2, cells = anotatedAsKeratinocytes_CellsGSM7717090_obj2) <- "Keratinocytes"
Idents(GSM7717090_obj2, cells = anotatedAsNK_CellsGSM7717090_obj2) <- "NK Cells"
Idents(GSM7717090_obj2, cells = anotatedAsLangerhans_CellsGSM7717090_obj2) <- "Langerhans"
Idents(GSM7717090_obj2, cells = anotatedAsCD8_CellsGSM7717090_obj2) <- "CD8+ T Cells"
Idents(GSM7717090_obj2, cells = anotatedAsCD4_CellsGSM7717090_obj2) <- "CD4+ T Cells"
Idents(GSM7717090_obj2, cells = anotatedAscDC_CellsGSM7717090_obj2) <- "cDC"
Idents(GSM7717090_obj2, cells = anotatedAsm1mac_CellsGSM7717090_obj2) <- "M1 Macrophages"
Idents(GSM7717090_obj2, cells = anotatedAsm2Macrophages_CellsGSM7717090_obj2) <- "M2 Macrophages"
Idents(GSM7717090_obj2, cells = anotatedAsm3_CellsGSM7717090_obj2) <- "Phenotype-Switching Macrophages"
Idents(GSM7717090_obj2, cells = anotatedAsMast_CellsGSM7717090_obj2) <- "Mast Cells"
Idents(GSM7717090_obj2, cells = anotatedAsB_CellsGSM7717090_obj2) <- "B Cells"
Idents(GSM7717090_obj2, cells = anotatedAsEndothelial_CellsGSM7717090_obj2) <- "Endothelial"
Idents(GSM7717090_obj2, cells = anotatedAsPericytes_CellsGSM7717090_obj2) <- "Pericytes"
Idents(GSM7717090_obj2, cells = anotatedAsMelanocytes_CellsGSM7717090_obj2) <- "Melanocytes"
Idents(GSM7717090_obj2, cells = anotatedAsPlasma_CellsGSM7717090_obj2) <- "Plasma Cells"
Idents(GSM7717090_obj2, cells = anotatedAsTreg_CellsGSM7717090_obj2) <- "Treg"
Idents(GSM7717090_obj2, cells = anotatedAsNeutrophil_CellsGSM7717090_obj2) <- "Neuronal Cells"

DimPlot(GSM7717090_obj2, reduction = "umap", label = T,repel = T)

Idents(GSM7717090_obj2, cells = anotatedAsEpithelial_CellsGSM7717090_obj2) <- "Epithelial"



#creating SuratObject
GSM7717090_obj= CreateSeuratObject(GSM7717090_RAW, min.cells =3, min.features = 100 )
dim(GSM7717090_obj)

GSM7717090_obj= readRDS("/Users/hossein.allahdadi/Downloads/GSE241132/GSM7717090/GSM7717090_obj.rds")
dim(GSM7717090_obj)

cell_types <- c("Basal_epi", "Suprabasal_epi", "Wound_edges", "Hair_follicle", "Papillary_dermis",
                "Sebaceous_gland", "Sweat_gland", "Sweatgland_FB", "FB1", "FB2", "Dermis1",
                "Dermis2", "Smooth_muscle", "Immune", "Immune_endo", "LE", "Mast_cell")

#MT count
GSM7717090_obj[["MTpercent"]]=PercentageFeatureSet(GSM7717090_obj, pattern = "^MT-")
VlnPlot(GSM7717090_obj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)
plot1= FeatureScatter(GSM7717090_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#filtering
GSM7717090_obj=subset(GSM7717090_obj, subset = nFeature_RNA<7500 &  MTpercent<25 & nCount_RNA<90000)

#normalization
GSM7717090_obj= NormalizeData(GSM7717090_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#highest variable
GSM7717090_obj=FindVariableFeatures(GSM7717090_obj,selection.method = "vst", nfeatures = 2000)
topVariableFeatures_GSM7717090=head(VariableFeatures(GSM7717090_obj),80)
plot_variableFeature_GSM7717090=VariableFeaturePlot(GSM7717090_obj)
plot_variableFeatureGSM7717090__tagged=LabelPoints(plot = plot_variableFeature_GSM7717090, points = topVariableFeatures_GSM7717090, repel = TRUE)
plot_variableFeatureGSM7717090__tagged

#scaling for PCA
GSM7717090_obj=ScaleData(GSM7717090_obj, features = rownames(GSM7717090_obj))
GSM7717090_obj=RunPCA(GSM7717090_obj,features = VariableFeatures(GSM7717090_obj))

#plot PCA
VizDimLoadings(GSM7717090_obj, dims = 1:5,  reduction = "pca")
ElbowPlot(GSM7717090_obj)
DimPlot(GSM7717090_obj, reduction = "pca")
DimHeatmap(GSM7717090_obj, dims = 1:6, cells = 200)


#Jak
GSM7717090_obj=JackStraw(GSM7717090_obj,num.replicate = 100, dims = 50)
GSM7717090_obj=ScoreJackStraw(GSM7717090_obj, dims = 1:50)
JackStrawPlot(GSM7717090_obj, 1:50)

# find clusters
GSM7717090_obj= FindNeighbors(GSM7717090_obj, dims = 1:50)
GSM7717090_obj= FindClusters(GSM7717090_obj, resolution = 3.1)

#tsne
GSM7717090_obj= RunTSNE(GSM7717090_obj, dims = 1:50)
DimPlot(GSM7717090_obj, reduction = "tsne",label.size =4 ,repel = T, label = T)

#umap
GSM7717090_obj=RunUMAP(GSM7717090_obj, dims = 1:50)
DimPlot(GSM7717090_obj,reduction = "umap", label.size =4 ,repel = T, label = T)

#markers:
GSM7717090_allMarkers=FindAllMarkers(GSM7717090_obj,min.pct = 0.25, logfc.threshold = 0.5 )
doubel_markers_of_GSM7717090_obj= GSE241136_allMarkers[which(GSM7717090_allMarkers$avg_log2FC>3),]
View(doubel_markers_of_GSM7717090_obj)
dim(doubel_markers_of_GSM7717090_obj)

markersOfClusterAGSM7717086=FindMarkers(GSM7717090_obj2, ident.1 = "A")
View(markersOfClusterAGSM7717086)
markersOfClusterBGSM7717086=FindMarkers(GSM7717090_obj2, ident.1 = "B")
View(markersOfClusterBGSM7717086)
markersOfClusterCGSM7717086=FindMarkers(GSM7717090_obj2, ident.1 = "C")
View(markersOfClusterCGSM7717086)
markersOfClusterDGSM7717086=FindMarkers(GSM7717090_obj2, ident.1 = "D")
View(markersOfClusterDGSM7717086)
write.csv(markersOfClusterAGSM7717086, file = "C:/Users/abbas/Documents/GSE241132/GSM7717090/GSM7717090_AclusterMarkers.csv")
# at the end of your analysis session, save the object
saveRDS(GSM7717090_obj, file = "/Users/hossein.allahdadi/Downloads/GSM7717090/GSM7717090_obj.rds")

# cell_type annotation
Cell_types_GSM7717090=c("Keratinocyte","Fibroblast","Keratinocyte","Keratinocyte","Mast Cells","Naive T cell",
                        "Effector T cell","T cell","Endothelial","Keratinocyte","Fibroblast","Pericyte",
                        "Fibroblast","Keratinocyte","Keratinocyte","Keratinocyte","Keratinocyte Stem cells",
                        "Endothelial","Th1","DC","Pericyte","Melanocyte")
length(Cell_types_GSM7717090)


Cell_types_GSM7717090_obj=c("Plasma Cells","Phenotype-Switching Macrophages","NK Cells","CD8+ T Cells",
                             "Treg","CD4+ T Cells","Keratinocytes","Fibroblasts",
                             "Langerhans",
                             "M1 Macrophages","M2 Macrophages","cDC",
                             "Mast Cells","B Cells",
                             "Endothelial","Pericytes","Melanocytes","Neuronal Cells"
)
length(Cell_types_GSM7717090_obj)
GSM7717090_obj2=GSM7717090_obj
names(x=Cell_types_GSM7717090_obj)=levels(x=GSM7717090_obj2)
GSM7717090_obj2=RenameIdents(object = GSM7717090_obj2, Cell_types_GSM7717090_obj
                             )

DimPlot(GSM7717090_obj2, reduction = "umap", label = T)
DimPlot(GSM7717090_obj2, reduction = "umap", label = T,repel = T)
Idents(GSM7717090_obj2) <- factor(Idents(GSM7717090_obj2), levels = newLevelsGSM7717090)


newLevelsGSM7717090=unique(Cell_types_GSM7717090)
newLevelsGSM7717090=c("Keratinocytes","Fibroblasts","NK Cells","Langerhans","CD8+ T Cells","CD4+ T Cells",
                      "Treg", "M1 Macrophages","M2 Macrophages","cDC","Phenotype-Switching Macrophages",
                      "Mast Cells","B Cells","Plasma Cells",
                      "Endothelial","Pericytes","Melanocytes","Neuronal Cells"
                     )
length(newLevelsGSM7717090)
xhsGSM7717090=data.frame()
expressionIdentsGSMGSM7717090obj2_list <- list()
for(celltype in newLevelsGSM7717090){
  subset_obj <- subset(GSM7717090_obj2, idents = celltype)
  expr_data <- GetAssayData(subset_obj, slot = "data", assay = "RNA")
  expressionIdentsGSMGSM7717090obj2_list[[celltype]] <- expr_data
  xhsGSM7717090=rbind(xhsGSM7717090,data.frame(barcodes=colnames(expr_data),ident=celltype))
}
write.csv(xhsGSM7717090, file = "C:/Users/abbas/Documents/GSE241132/GSM7717090/anotaionGSM7717090.csv")
dim(xhsGSM7717090)

#Keratinocytes umap
FeaturePlot(GSM7717090_obj,features = "KRT14", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "KRT5", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "KRT1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "KRT10", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CDH3", label = T, reduction = "umap")#
#basal layer
FeaturePlot(GSM7717090_obj,features = "IVL", label = T, reduction = "umap")#
#Fibroblasts
FeaturePlot(GSM7717090_obj,features = "COL1A1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "COL1A2", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "COL5A1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "LUM", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "FBLN1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "FBLN2", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "PDGFRA", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "DCN", label = T, reduction = "umap")#

#Merkels
#FeaturePlot(GSM7717090_obj,features = "NCAM1", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "ATOH1", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "ISL1", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "KRT20", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "KRT8", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "KRT18", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "SOX2", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "LUM", label = T, reduction = "umap")#

#langerhans
FeaturePlot(GSM7717090_obj,features = "CD207", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD1A", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "HLA‑DRB1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD74", label = T, reduction = "umap")#

#Melanocytes umap
FeaturePlot(GSM7717090_obj,features = "MITF", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "PMEL", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "MLANA", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "TYR", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "TYRP1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "DCT", label = T, reduction = "umap")#
#Endothelial umap
FeaturePlot(GSM7717090_obj,features = "PECAM1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "VWF", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CDH5", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "KDR", label = T, reduction = "umap")#
#Pericyte umap
FeaturePlot(GSM7717090_obj,features = "PDGFRB", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CSPG4", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "RGS5", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "MCAM", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "TAGLN", label = T, reduction = "umap")#
#Adipocyte umap
FeaturePlot(GSM7717090_obj,features = "PPARG", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CEBPA", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "ADIPOQ", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "APOE", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "LPL", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "FABP4", label = T, reduction = "umap")#


#Bcell
FeaturePlot(GSM7717090_obj,features = "CD19", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "PAX5", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "CD79A", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD79B", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "MS4A1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD22", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "JCHAIN", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "MZB1", label = T, reduction = "umap")#

#nk
FeaturePlot(GSM7717090_obj,features = "FCGR3A", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "NCAM1", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "CD14", label = T, reduction = "umap")#5,19
FeaturePlot(GSM7717090_obj,features = "KLRD1", label = T, reduction = "umap")#15
FeaturePlot(GSM7717090_obj,features = "KLRB1", label = T, reduction = "umap")#8,15

# ================================================================
#Tcell umap
FeaturePlot(GSM7717090_obj,features = "CD3D", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD1C", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD69", label = T, reduction = "umap")#

FeaturePlot(GSM7717090_obj,features = "CD4", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "IL7R", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "CCR7", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "SELL", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "FOXP3", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "IL17A", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "CTLA4", label = T, reduction = "umap")

FeaturePlot(GSM7717090_obj,features = "CD8B", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD8A", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "TBX21", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "ITGAE", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "ITGA1", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "PDCD1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "GZMA", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "GZMB", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "PRF1", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "IFNG", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "ZNF683", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj, features = "CCR7", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CXCR6", label = T, reduction = "umap")#

#Effector Tcell
FeaturePlot(GSM7717090_obj,features = "ICOS", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD62L", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "CCR6", label = T, reduction = "umap")#

#CD4
FeaturePlot(GSM7717090_obj,features = "FHIT", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "TSHZ2", label = T, reduction = "umap")#
#previuosly V.0.4 Tcell
FeaturePlot(GSM7717090_obj,features = "CD3E", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD3D", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD4", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD8A", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "ITGAE", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "CD69", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "PDCD1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "TBX21", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD1C", label = T, reduction = "umap")#

#Mast cell
FeaturePlot(GSM7717090_obj,features = "HDC", label = T, reduction = "umap")#18
FeaturePlot(GSM7717090_obj,features = "HPGDS", label = T, reduction = "umap")#18


#Macro umap #
FeaturePlot(GSM7717090_obj,features = "MARCO", label = T, reduction = "umap")#5
FeaturePlot(GSM7717090_obj,features = "CSF1R", label = T, reduction = "umap")#5
FeaturePlot(GSM7717090_obj,features = "CD68", label = T, reduction = "umap")#5
FeaturePlot(GSM7717090_obj,features = "CD86", label = T, reduction = "umap")#5,6
FeaturePlot(GSM7717090_obj,features = "CLEC10A", label = T, reduction = "umap")#6
FeaturePlot(GSM7717090_obj,features = "CD80", label = T, reduction = "umap")#6
FeaturePlot(GSM7717090_obj, features = "CD163" , reduction = "umap", label = T) #5
FeaturePlot(GSM7717090_obj, features = "AIF1" , reduction = "umap", label = T) #5,6
FeaturePlot(GSM7717090_obj, features = "ITGAM", reduction = "umap", label = T ) #5
FeaturePlot(GSM7717090_obj, features = "PTPRC" , reduction = "umap", label = T) #5,6,8,15
FeaturePlot(GSM7717090_obj, features = "CX3CR1" , reduction = "umap", label = T) #6
FeaturePlot(GSM7717090_obj, features = "MSR1" , reduction = "umap", label = T) #5
FeaturePlot(GSM7717090_obj, features = "CCL18" , reduction = "umap", label = T) #5
FeaturePlot(GSM7717090_obj, features = "CD14" , reduction = "umap", label = T) #5
FeaturePlot(GSM7717090_obj, features = "SIGLEC1" , reduction = "umap", label = T)#5
FeaturePlot(GSM7717090_obj, features = "MERTK" , reduction = "umap", label = T)#5

FeaturePlot(GSM7717090_obj, features = "SNHG26" , reduction = "umap", label = T)#5

#M1
FeaturePlot(GSM7717090_obj, features = "CD86", min.cutoff = 0.2,   label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M1, can appear on activated M2 in disease.
FeaturePlot(GSM7717090_obj, features = "CD80",   label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M1; activation-associated.
FeaturePlot(GSM7717090_obj, features = "CD68", min.cutoff = 0.2,  label = TRUE, reduction = "umap") # Specificity: low (pan-macrophage); Expression: high in both subtypes.
FeaturePlot(GSM7717090_obj, features = "CD40",min.cutoff = 0.2, label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "HLA-DRA",label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M1 (IFN-γ driven), mid in M2.
FeaturePlot(GSM7717090_obj, features = "NOS2", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "STAT1",min.cutoff = 0.01, label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "STAT2",min.cutoff = 0.01, label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IRF4", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "TLR4",min.cutoff = 0.1, label = TRUE, reduction = "umap") # Specificity: highish; Expression: ↑M1 (LPS pathway hallmark), low in M2.
FeaturePlot(GSM7717090_obj, features = "TLR2",   label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M1, low–mid in M2.
FeaturePlot(GSM7717090_obj, features = "NLRP3",  label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M1 (feedback brake in inflammatory signaling).
FeaturePlot(GSM7717090_obj, features = "TNF", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IL1B", min.cutoff = 1.0, label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M1 (feedback brake in inflammatory signaling).
FeaturePlot(GSM7717090_obj, features = "IL1R2",  label = TRUE, min.cutoff = 0.5, reduction = "umap") # Specificity: moderate; Expression: ↑M1 (signaling receptor), M2 tends to ↑IL1R2.
FeaturePlot(GSM7717090_obj, features = "IL6", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IL12A", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IL12B", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IL23A", min.cutoff = 0.1, label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IL1R2", min.cutoff = 0.3,  label = TRUE,  reduction = "umap") # Specificity: moderate; Expression: ↑M1 (signaling receptor), M2 tends to ↑IL1R2.
FeaturePlot(GSM7717090_obj, features = "CLEC10A", min.cutoff = 0.1, label = F, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CXCL8", min.cutoff = 3.2, label = TRUE, reduction = "umap") # Specificity: high (M1); Expression: ↑M1; low/undetectable in human skin scRNA.
FeaturePlot(GSM7717090_obj, features = "CXCL10", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CXCL5", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IRF5", label = T, reduction = "umap")#

FeaturePlot(GSM7717090_obj, features = "IRF4", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "STAT6", label = T, reduction = "umap")#
#Macrophage M2
FeaturePlot(GSM7717090_obj, features = "CD163",min.cutoff = 1.2,  label = TRUE, reduction = "umap") # Specificity: high; Expression: ↑M2, low in M1.
FeaturePlot(GSM7717090_obj, features = "MSR1",min.cutoff = 0.8, label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "MMP9",min.cutoff = 1.0, label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "MMP9",    label = TRUE,min.cutoff = 0.5, reduction = "umap") # Specificity: moderate; Expression: ↑M2/TAM remodeling programs.
FeaturePlot(GSM7717090_obj, features = "MRC1" ,min.cutoff = 0.9, label = TRUE, reduction = "umap") # Specificity: high; Expression: ↑M2 (signature), low in M1.
FeaturePlot(GSM7717090_obj, features = "TGFB1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IL10", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CCL18", min.cutoff = 0.1, label = TRUE, reduction = "umap") # Mouse-only M2; Expression: ↑M2 (human may show CHI3L1 instead).
FeaturePlot(GSM7717090_obj, features = "IL1R2", min.cutoff = 0.1,   label = TRUE, reduction = "umap") # Specificity: highish; Expression: ↑M2 decoy receptor; minimal in M1.
FeaturePlot(GSM7717090_obj, features = "PPARG", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "ARG1",    label = TRUE, reduction = "umap") # Specificity: high in mouse; Expression: ↑M2; often low in human scRNA.
FeaturePlot(GSM7717090_obj, features = "CHI3L1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CD200R1", label = TRUE, reduction = "umap") # Specificity: moderate–high; Expression: ↑M2/regulatory; low M1.
FeaturePlot(GSM7717090_obj, features = "TLR1",    label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M2 relative to M1.
FeaturePlot(GSM7717090_obj, features = "TLR8",    label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M2 in some tissues.
FeaturePlot(GSM7717090_obj, features = "VEGFA",   label = TRUE, reduction = "umap") # Specificity: moderate; Expression: ↑M2 (angiogenesis/wound repair); variable.
FeaturePlot(GSM7717090_obj, features = "TGM2",    label = TRUE, reduction = "umap") # Specificity: highish; Expression: ↑M2; low M1.
FeaturePlot(GSM7717090_obj, features = "RETNLA",  label = TRUE, reduction = "umap") # Mouse-only M2; Expression: ↑M2 (do not expect in human).
FeaturePlot(GSM7717090_obj, features = "CHI3L1",  label = TRUE, reduction = "umap") # Mouse-only M2; Expression: ↑M2 (human may show CHI3L1 instead).
FeaturePlot(GSM7717090_obj, features = "CD14" ,  min.cutoff = 0.3, reduction = "umap") #
FeaturePlot(GSM7717090_obj, features = "CSF1R",  min.cutoff = 0.2, label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "SIGLEC1" , reduction = "umap", label = T)#
FeaturePlot(GSM7717090_obj, features = "MERTK" , reduction = "umap", label = T)#
FeaturePlot(GSM7717090_obj, features = "PTPRC" , reduction = "umap", label = T) #
FeaturePlot(GSM7717090_obj, features = "CX3CR1" , reduction = "umap", label = T) #
FeaturePlot(GSM7717090_obj, features = "CCL18" ,min.cutoff = 0.4, reduction = "umap", label = T) #
FeaturePlot(GSM7717090_obj, features = "CCR2", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "S100A8", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "S100A9", label = T, reduction = "umap")#
# ================================================================
# ================================================================
#Neut umap
FeaturePlot(GSM7717090_obj,features = "MPO", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CEBPE", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "FCGR3B", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj,features = "CD177", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "ELANE", label = T, reduction = "umap")
FeaturePlot(GSM7717090_obj,features = "CXCR2", label = T, reduction = "umap")


# ================================================================
#cDC1
FeaturePlot(GSM7717090_obj, features = "CLEC9A", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "XCR1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "BATF3", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IRF8", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CADM1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "THBD", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "HLA-DQB1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "HLA-DRA", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CD74", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CLNK", label = T, reduction = "umap")#
#cDC2
FeaturePlot(GSM7717090_obj, features = "CD1C", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "FCER1A", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CLEC10A", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "SIRPA", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "IRF4", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "KLF4", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "HLA-DQB1", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "TNFSF11", label = T, reduction = "umap")#
FeaturePlot(GSM7717090_obj, features = "CCR6", label = T, reduction = "umap")#
#pDC
FeaturePlot(GSM7717090_obj, features = "IL3RA", label = T, reduction = "umap") # CD123, canonical pDC marker
FeaturePlot(GSM7717090_obj, features = "CLEC4C", label = T, reduction = "umap") # BDCA-2, highly specific to pDC
FeaturePlot(GSM7717090_obj, features = "LILRA4", label = T, reduction = "umap") # ILT7, inhibitory receptor, pDC hallmark
FeaturePlot(GSM7717090_obj, features = "TCF4", label = T, reduction = "umap")   # Lineage-defining TF
FeaturePlot(GSM7717090_obj, features = "IRF7", label = T, reduction = "umap")   # High in pDC, type I IFN amplifier
FeaturePlot(GSM7717090_obj, features = "TLR7", label = T, reduction = "umap")   # Viral RNA sensor
FeaturePlot(GSM7717090_obj, features = "TLR9", label = T, reduction = "umap")   # CpG DNA sensor
FeaturePlot(GSM7717090_obj, features = "BST2", label = T, reduction = "umap")   # IFN-inducible, en
#Langerhans
FeaturePlot(GSM7717090_obj, features = "CD207", label = T, reduction = "umap") # Exclusive LC marker (Langerin)
FeaturePlot(GSM7717090_obj, features = "CD1A", label = T, reduction = "umap") # Strong LC marker in epidermis
FeaturePlot(GSM7717090_obj, features = "HLA-DRA", label = T, reduction = "umap") # MHC-II, shared APC signature
FeaturePlot(GSM7717090_obj, features = "CD80", label = T, reduction = "umap") # Costimulatory, upregulated on activation
FeaturePlot(GSM7717090_obj, features = "CD86", label = T, reduction = "umap") # Costimulatory, LC maturation
FeaturePlot(GSM7717090_obj, features = "CCR7", label = T, reduction = "umap") # Migratory LC marker after activation
FeaturePlot(GSM7717090_obj, features = "EPCAM", label = T, reduction = "umap") # Epidermal retention, LC residency
FeaturePlot(GSM7717090_obj, features = "IL15RA", label = T, reduction = "umap") # Survival in epidermal niche
FeaturePlot(GSM7717090_obj, features = "CSF1R", label = T, reduction = "umap") # Myeloid lineage marker
FeaturePlot(GSM7717090_obj, features = "CLEC4E", label = T, reduction = "umap") # Pathogen recognition receptor


#Pan-neuronal structural markers
FeaturePlot(GSM7717090_obj, features = "TUBB3", label = T, reduction = "umap")  # βIII-Tubulin, core neuronal cytoskeleton
FeaturePlot(GSM7717090_obj, features = "RBFOX3", label = T, reduction = "umap") # NeuN, nuclear neuronal identity
FeaturePlot(GSM7717090_obj, features = "UCHL1", label = T, reduction = "umap")  # PGP9.5, pan-neuronal cytoplasmic marker
FeaturePlot(GSM7717090_obj, features = "PRPH",  label = T, reduction = "umap")  # Peripherin, peripheral sensory axons
#Sensory neuron receptors & ion channels
FeaturePlot(GSM7717090_obj, features = "NTRK1", label = T, reduction = "umap")  # TrkA, hallmark nociceptor receptor (NGF)
FeaturePlot(GSM7717090_obj, features = "SCN9A", label = T, reduction = "umap")  # Nav1.7 sodium channel, pain signaling
FeaturePlot(GSM7717090_obj, features = "SCN10A",label = T, reduction = "umap")  # Nav1.8 sodium channel, nociceptors
FeaturePlot(GSM7717090_obj, features = "PIEZO2",label = T, reduction = "umap")  # Mechanosensation, touch/pressure
#Neuropeptides & neurotransmitters
FeaturePlot(GSM7717090_obj, features = "CALCA", label = T, reduction = "umap")  # CGRP, nociceptive neuropeptide
FeaturePlot(GSM7717090_obj, features = "TAC1",  label = T, reduction = "umap")  # Substance P, pain & neurogenic inflammation
FeaturePlot(GSM7717090_obj, features = "CHAT",  label = T, reduction = "umap")  # Choline acetyltransferase, cholinergic subset
#Synaptic machinery
FeaturePlot(GSM7717090_obj, features = "SYN1",  label = T, reduction = "umap")  # Synapsin I, presynaptic vesicle regulation
FeaturePlot(GSM7717090_obj, features = "SNAP25",label = T, reduction = "umap")  # Synaptic vesicle fusion machinery
FeaturePlot(GSM7717090_obj, features = "SYT4",  label = T, reduction = "umap")  # Synaptotagmin-4, Ca²⁺-regulated vesicle release
#Subtype & functional markers
FeaturePlot(GSM7717090_obj, features = "TH",    label = T, reduction = "umap")  # Tyrosine hydroxylase, sympathetic catecholaminergic
FeaturePlot(GSM7717090_obj, features = "PVALB", label = T, reduction = "umap")  # Parvalbumin, proprioceptive mechanosensory neurons
FeaturePlot(GSM7717090_obj, features = "GFRA2", label = T, reduction = "umap")  # GDNF receptor, trophic sensory subset

FeaturePlot(GSM7717090_obj, features = "KRT20", label = T, reduction = "umap") # Exclusive Merkel marker
FeaturePlot(GSM7717090_obj, features = "SOX2", label = T, reduction = "umap")  # Merkel progenitor & identity
FeaturePlot(GSM7717090_obj, features = "INSM1", label = T, reduction = "umap") # Neuroendocrine Merkel identity
FeaturePlot(GSM7717090_obj, features = "PIEZO2", label = T, reduction = "umap")# Touch-sensing channel
FeaturePlot(GSM7717090_obj, features = "SYT4", label = T, reduction = "umap")  # Synaptic vesicle machinery
FeaturePlot(GSM7717090_obj, features = "CHRNB2", label = T, reduction = "umap")# Cholinergic receptor
FeaturePlot(GSM7717090_obj, features = "NCAM1", label = T, reduction = "umap") # Neural adhesion
FeaturePlot(GSM7717090_obj, features = "TRPV4", label = T, reduction = "umap") # Mechanosensitive Ca2+ channel
# ================================================================
# ================================================================
# ================================================================

plotcdsevenGSM7717090=FeaturePlot(GSM7717090_obj2, features = "CD207", label = T, reduction = "umap")
selected.cellsasLangerhans086obj2 <- CellSelector(plot = plotcdsevenGSM7717090086)
length(selected.cellsasLangerhans086obj2)
Idents(GSM7717090_obj2, cells = selected.cellsasLangerhans) <- "Langerhans"
table(Idents(GSM7717090_obj2))
GSM7717090_obj2=GSM7717090_obj

DimPlot(GSM7717090_obj2, reduction = "umap", label = T,repel = T)
#NK cell assigning
selected.cellsAsNKCell086=CellSelector(plot = plotcdsevenGSM7717090)
length(selected.cellsAsNKCell086)
Idents(GSM7717090_obj2, cells = selected.cellsAsNKCell086) <- "NK cell"

#Tcell assigning CD4
selected.cellsAsCD4T_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsCD4T_Cell086) <- "CD4+ T cell"

#m1 macrophages assigning CD4
selected.cellsAsm1macrophages_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsm1macrophages_Cell086) <- "M1 Macrophages"

#m2 macrophages assigning CD4
selected.cellsAsm2macrophages_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsm2macrophages_Cell086) <- "M2 Macrophages"

#m3 macrophages assigning CD4
selected.cellsAsm3macrophages_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsm3macrophages_Cell086) <- "Phenotype-Switching Macrophages"

#Tcell assigning CD8
selected.cellsAsCD8T_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsCD8T_Cell086) <- "CD8+ T cell"
#melanocytes assigning CD8
selected.cellsAsmelanocytes_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsmelanocytes_Cell086) <- "Melanocytes"

#Tcell assigning CD4
selected.cellsAsCD4T_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsCD4T_Cell086) <- "CD4+ T cell"

#Keratinocytes assigning CD4
selected.cellsAsKeratinocytes_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsKeratinocytes_Cell086) <- "Keratinocytes"

#pericyte assigning CD8
selected.cellsAsPericytes_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsPericytes_Cell7717090) <- "Pericytes"

#Fibroblasts assigning CD4
selected.cellsAsFibroblasts_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsFibroblasts_Cell086) <- "Fibroblasts"

#Endothelial assigning CD4
selected.cellsAsEndothelial_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsEndothelial_Cell086) <- "Endothelial"

#Tcell assigning CD8
selected.cellsAsCD8T_Cell086=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsCD8T_Cell086) <- "CD8+ T cell"

#pericyte assigning CD8
selected.cellsAsmast_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsmast_Cell7717090) <- "Mast Cells"

#Malanocytes assigning CD8
selected.cellsAsMelanocytes_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsMelanocytes_Cell7717090) <- "Malanocytes"


#B cell assigning CD8
selected.cellsAsBCells_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsBCells_Cell7717090) <- "B Cells"

#B cell assigning CD8
selected.cellsAsplasmaCells_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsplasmaCells_Cell7717090) <- "Plasma Cells"


#B cell assigning CD8
selected.cellsAsTreg_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsTreg_Cell7717090) <- "Treg"

#A assigning CD8
selected.cellsAsA_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsA_Cell7717090) <- "A"
Idents(GSM7717090_obj2, cells = selected.cellsAsA_Cell7717090) <- "Neuronal Cells"

#A assigning CD8
selected.cellsAsB_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsB_Cell7717090) <- "B"
Idents(GSM7717090_obj2, cells = selected.cellsAsB_Cell7717090) <- "Langerhans"

#A assigning CD8
selected.cellsAsC_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsC_Cell7717090) <- "C"
Idents(GSM7717090_obj2, cells = selected.cellsAsC_Cell7717090) <- "cDC"

#A assigning CD8
selected.cellsAsD_Cell7717090=CellSelector(plot = plotcdsevenGSM7717090)
Idents(GSM7717090_obj2, cells = selected.cellsAsD_Cell7717090) <- "D"
Idents(GSM7717090_obj2, cells = selected.cellsAsD_Cell7717090) <- "Endothelial"
# ================================================================
# ========================






table(Idents(GSM7717090_obj3))

macrophage_markers090 = c("IL1B","TLR2","TLR4","CXCL5","STAT1","STAT2","CD86",
                         "CD163","CD68","IRF4","IL10","CCL18","MRC1","MSR1", "MMP9")
# ================================================================
target_clusters090 <- c("14", "19","40","8")  # Replace with your actual cluster names/numbers
subset_90macros <- subset(GSM7717090_obj, idents = target_clusters090)
dim(subset_90macros)
table(Idents(subset_90macros))

subset90CellType= c("M2 Macrophage","M1 Macrophage","Phenotype-Switching Macrophages","Phenotype-Switching Macrophages")
names(x=subset90CellType)=levels(x=subset_90macros)
subset_90macros=RenameIdents(object = subset_90macros, subset90CellType)
# Create classification column
subset_90macros$macrophage_type_v2 <- "Unclassified"
subset_88macros$macrophage_type_v2[colnames(subset_88macros) %in% m1pure] <- "M1"
subset_88macros$macrophage_type_v2[colnames(subset_88macros) %in% m2pure] <- "M2"
subset_88macros$macrophage_type_v2[colnames(subset_88macros) %in% transitional_cd163_posit_ccl18_neg] <- "TransitionallM1"

subset_90macros$macrophage_id <- Idents(subset_90macros)

df90macrophages <- FetchData(subset_90macros, vars = c(macrophage_markers090, "macrophage_id"))

df_long090 <- df90macrophages %>%
  tidyr::pivot_longer(cols = all_of(macrophage_markers090),
                      names_to = "gene", values_to = "expr")

# Visualization of results
# ========================
table(Idents(subset_88macros))
# Plot scores
p1 <- FeatureScatter(subset_90macros, feature1 = "M1_score", feature2 = "M2_score", group.by = "macrophage_type")
# UMAP plots
p2 <- DimPlot(subset_90macros, group.by = "macrophage_id", reduction = "umap")
p3 <- DimPlot(subset_90macros, group.by = "macrophage_type_v2", reduction = "umap")
# Feature plots with custom thresholds
p5 <- FeaturePlot(subset_88macros, features = "IL1B", min.cutoff = il1b_threshold)
p6 <- FeaturePlot(subset_88macros, features = "MRC1", min.cutoff = mrc1_threshold)
FeaturePlot(subset_88macros,features = "CD86", label = T, reduction = "umap")#
featuresDimPlotMacro
dotplot090=DotPlot(subset_90macros, features = macrophage_markers090) + RotatedAxis()
dotplot090+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom",   # put legends below plot
    legend.box = "vertical"     # stack them horizontally
  )
DotPlot(subset_90macros,features = macrophage_markers090 )
# Visualization DotPlot
# ========================
p090 = DotPlot(subset_90macros, features = macrophage_markers090)+ RotatedAxis()+
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",   # put legends below plot
    legend.box = "horizontal"     # stack them horizontally
  ) +
  guides(
    size = guide_legend(title = "Fraction of cells in group (%)",
                        title.position = "top", direction = "horizontal"),
    color = guide_colorbar(title = "Mean expression in group",
                           title.position = "top", barwidth = 10, barheight = 0.5)
  )
p090 + coord_fixed(ratio = 0.3)



# Pull expression (from the current assay’s data slot; usually log-normalized)
df <- FetchData(obj_mm, vars = c(macrophage_markers, "macrophage_type_v2"))
obj_mm <- subset(subset_88macros, idents = c("M1","M2"))

# Long format for ggplot
df_long090 <- df %>%
  pivot_longer(cols = all_of(macrophage_markers),
               names_to = "gene", values_to = "expr") %>%
  mutate(macrophage_type_v2 = factor(macrophage_type_v2, levels = c("M1","M2")))


#Dark orange (#B75A1D to #A64B00)
#Light orange (#F6C25B to #F5B942)
#Pale yellow (#FAF4C6 to #F5EDAD)
#Medium cyan/light blue (#6DC1E2 to #5AB5DB)
#Darker blue (#3A90B5 to #317FA3)
cols <- c(M1 = "#44B4FF",  # light blue
          M2 = "#9B1C1C")  # dark red
boxCols=c( "#9B1C1C", "#44B4FF" ,"#1F9D7A")
p_box090 <- ggplot(df_long090, aes(x = macrophage_id, y = expr, fill = macrophage_id)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  scale_fill_manual(values = boxCols, guide = "none") +
  facet_wrap(~ gene, ncol = 3, scales = "free_y") +
  labs(x = NULL, y = "Expression (data slot)") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank()) + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

p_box090
greens= c("#3BA99C","#1F9D7A","#20C997","#00A087")
