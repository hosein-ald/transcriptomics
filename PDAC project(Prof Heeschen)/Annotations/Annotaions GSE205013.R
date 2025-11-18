#GSM6204109

anotation_for_GSM6204109_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204109.csv", sep = "," )
dim(anotation_for_GSM6204109_obj2)
anotation_for_GSM6204109_obj2[1:6,]
colnames(anotation_for_GSM6204109_obj2)=anotation_for_GSM6204109_obj2[1,]
anotation_for_GSM6204109_obj2=anotation_for_GSM6204109_obj2[-1,]
length(unique(anotation_for_GSM6204109_obj2$ident))
length(anotatedAsFibroblastGSM6204109_obj2)
anotatedAsNK_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "NK Cells"),2]
anotatedAstrNK_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "trNK"),2]
anotatedAsCD8_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "CD4+ T Cells"),2]
anotatedAsTreg_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "Treg"),2]
anotatedAsMast_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "Mast Cells"),2]
anotatedAsB_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "Plasma Cells"),2]
anotatedAsMac_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "Macrophages"),2]
anotatedAsDC_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "DC"),2]
anotatedAsDuctalCells_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "Ductal Cells"),2]
anotatedAsFibroblastGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "Fibroblasts"),2]
anotatedAsEndothelial_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "Endothelial"),2]
anotatedAsHepatocytes_CellsGSM6204109_obj2=anotation_for_GSM6204109_obj2[which(anotation_for_GSM6204109_obj2$ident == "Hepatocytes"),2]


length(anotatedAsPlasma_CellsGSM6204109_obj2)
# Assigning Obj2
Idents(GSM6204109_obj2, cells = anotatedAsNK_CellsGSM6204109_obj2) <- "NK Cells"
Idents(GSM6204109_obj2, cells = anotatedAstrNK_CellsGSM6204109_obj2) <- "trNK Cells"
Idents(GSM6204109_obj2, cells = anotatedAsCD8_CellsGSM6204109_obj2) <- "CD8+ T Cells"
Idents(GSM6204109_obj2, cells = anotatedAsCD4_CellsGSM6204109_obj2) <- "CD4+ T Cells"
Idents(GSM6204109_obj2, cells = anotatedAsTreg_CellsGSM6204109_obj2) <- "Treg"
Idents(GSM6204109_obj2, cells = anotatedAsMast_CellsGSM6204109_obj2) <- "Mast Cells"
Idents(GSM6204109_obj2, cells = anotatedAsB_CellsGSM6204109_obj2) <- "B Cells"
Idents(GSM6204109_obj2, cells = anotatedAsPlasma_CellsGSM6204109_obj2) <- "Plasma Cells"
Idents(GSM6204109_obj2, cells = anotatedAsMac_CellsGSM6204109_obj2) <- "Macrophages"
Idents(GSM6204109_obj2, cells = anotatedAsDC_CellsGSM6204109_obj2) <- "DC"
Idents(GSM6204109_obj2, cells = anotatedAsDuctalCells_CellsGSM6204109_obj2) <- "Ductal Cells"
Idents(GSM6204109_obj2, cells = anotatedAsFibroblastGSM6204109_obj2) <- "Fibroblasts"
Idents(GSM6204109_obj2, cells = anotatedAsEndothelial_CellsGSM6204109_obj2) <- "Endothelial"
Idents(GSM6204109_obj2, cells = anotatedAsHepatocytes_CellsGSM6204109_obj2) <- "Hepatocytes"


DimPlot(GSM6204109_obj2, reduction = "umap", label = T,repel = T)



#GSM6204110

anotation_for_GSM6204110_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204110.csv", sep = "," )
dim(anotation_for_GSM6204110_obj2)
anotation_for_GSM6204110_obj2[1:6,]
colnames(anotation_for_GSM6204110_obj2)=anotation_for_GSM6204110_obj2[1,]
anotation_for_GSM6204110_obj2=anotation_for_GSM6204110_obj2[-1,]
length(unique(anotation_for_GSM6204110_obj2$ident))
length(anotatedAsFibroblastGSM6204110_obj2)
anotatedAsNK_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "NK Cells"),2]
anotatedAstrNK_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "trNK"),2]
anotatedAsCD8_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "CD4+ T Cells"),2]
anotatedAsTreg_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "Treg"),2]
anotatedAsMast_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "Mast Cells"),2]
anotatedAsB_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "Plasma Cells"),2]
anotatedAsMac_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "Macrophages"),2]
anotatedAsDC_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "DC"),2]
anotatedAsDuctalCells_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "Ductal Cells"),2]
anotatedAsFibroblastGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "Fibroblasts"),2]
anotatedAsEndothelial_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "Endothelial"),2]
anotatedAsHepatocytes_CellsGSM6204110_obj2=anotation_for_GSM6204110_obj2[which(anotation_for_GSM6204110_obj2$ident == "Hepatocytes"),2]


length(anotatedAsPlasma_CellsGSM6204110_obj2)
# Assigning Obj2
Idents(GSM6204110_obj2, cells = anotatedAsNK_CellsGSM6204110_obj2) <- "NK Cells"
Idents(GSM6204110_obj2, cells = anotatedAstrNK_CellsGSM6204110_obj2) <- "trNK Cells"
Idents(GSM6204110_obj2, cells = anotatedAsCD8_CellsGSM6204110_obj2) <- "CD8+ T Cells"
Idents(GSM6204110_obj2, cells = anotatedAsCD4_CellsGSM6204110_obj2) <- "CD4+ T Cells"
Idents(GSM6204110_obj2, cells = anotatedAsTreg_CellsGSM6204110_obj2) <- "Treg"
Idents(GSM6204110_obj2, cells = anotatedAsMast_CellsGSM6204110_obj2) <- "Mast Cells"
Idents(GSM6204110_obj2, cells = anotatedAsB_CellsGSM6204110_obj2) <- "B Cells"
Idents(GSM6204110_obj2, cells = anotatedAsPlasma_CellsGSM6204110_obj2) <- "Plasma Cells"
Idents(GSM6204110_obj2, cells = anotatedAsMac_CellsGSM6204110_obj2) <- "Macrophages"
Idents(GSM6204110_obj2, cells = anotatedAsDC_CellsGSM6204110_obj2) <- "DC"
Idents(GSM6204110_obj2, cells = anotatedAsDuctalCells_CellsGSM6204110_obj2) <- "Ductal Cells"
Idents(GSM6204110_obj2, cells = anotatedAsFibroblastGSM6204110_obj2) <- "Fibroblasts"
Idents(GSM6204110_obj2, cells = anotatedAsEndothelial_CellsGSM6204110_obj2) <- "Endothelial"
Idents(GSM6204110_obj2, cells = anotatedAsHepatocytes_CellsGSM6204110_obj2) <- "Hepatocytes"


DimPlot(GSM6204110_obj2, reduction = "umap", label = T,repel = T)


#GSM6204111

anotation_for_GSM6204111_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204111.csv", sep = "," )
dim(anotation_for_GSM6204111_obj2)
anotation_for_GSM6204111_obj2[1:6,]
colnames(anotation_for_GSM6204111_obj2)=anotation_for_GSM6204111_obj2[1,]
anotation_for_GSM6204111_obj2=anotation_for_GSM6204111_obj2[-1,]
length(unique(anotation_for_GSM6204111_obj2$ident))
length(anotatedAsFibroblastGSM6204111_obj2)
anotatedAsDuctalCells_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Ductal Cells"),2]
anotatedAsFibroblastGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Fibroblasts"),2]
anotatedAsEndothelial_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Endothelial"),2]
anotatedAsMac_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Macrophages"),2]
anotatedAsCD8_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "CD4+ T Cells"),2]
anotatedAsNK_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "NK Cells"),2]
anotatedAstrNK_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "trNK"),2]
anotatedAsMuscle_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Muscle Cells"),2]
anotatedAsHepatocytes_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Hepatocytes"),2]
anotatedAsDC_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "DC"),2]
anotatedAsMast_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Mast Cells"),2]
anotatedAsTreg_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Treg"),2]
anotatedAsB_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Plasma Cells"),2]
anotatedAsAcinar_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Acinar Cells"),2]
anotatedAsEndocrine_CellsGSM6204111_obj2=anotation_for_GSM6204111_obj2[which(anotation_for_GSM6204111_obj2$ident == "Endocrine Cells"),2]


length(anotatedAsPlasma_CellsGSM6204111_obj2)
# Assigning Obj2
Idents(GSM6204111_obj2, cells = anotatedAsDuctalCells_CellsGSM6204111_obj2) <- "Ductal Cells"
Idents(GSM6204111_obj2, cells = anotatedAsFibroblastGSM6204111_obj2) <- "Fibroblasts"
Idents(GSM6204111_obj2, cells = anotatedAsEndothelial_CellsGSM6204111_obj2) <- "Endothelial"
Idents(GSM6204111_obj2, cells = anotatedAsMac_CellsGSM6204111_obj2) <- "Macrophages"
Idents(GSM6204111_obj2, cells = anotatedAsCD8_CellsGSM6204111_obj2) <- "CD8+ T Cells"
Idents(GSM6204111_obj2, cells = anotatedAsCD4_CellsGSM6204111_obj2) <- "CD4+ T Cells"
Idents(GSM6204111_obj2, cells = anotatedAsNK_CellsGSM6204111_obj2) <- "NK Cells"
Idents(GSM6204111_obj2, cells = anotatedAstrNK_CellsGSM6204111_obj2) <- "trNK Cells"
Idents(GSM6204111_obj2, cells = anotatedAsMuscle_CellsGSM6204111_obj2) <- "Muscle Cells"
Idents(GSM6204111_obj2, cells = anotatedAsHepatocytes_CellsGSM6204111_obj2) <- "Hepatocytes"
Idents(GSM6204111_obj2, cells = anotatedAsDC_CellsGSM6204111_obj2) <- "DC"
Idents(GSM6204111_obj2, cells = anotatedAsMast_CellsGSM6204111_obj2) <- "Mast Cells"
Idents(GSM6204111_obj2, cells = anotatedAsTreg_CellsGSM6204111_obj2) <- "Treg"
Idents(GSM6204111_obj2, cells = anotatedAsB_CellsGSM6204111_obj2) <- "B Cells"
Idents(GSM6204111_obj2, cells = anotatedAsPlasma_CellsGSM6204111_obj2) <- "Plasma Cells"
Idents(GSM6204111_obj2, cells = anotatedAsAcinar_CellsGSM6204111_obj2) <- "Acinar Cells"
Idents(GSM6204111_obj2, cells = anotatedAsEndocrine_CellsGSM6204111_obj2) <- "Endocrine Cells"


DimPlot(GSM6204111_obj2, reduction = "umap", label = T,repel = T)



#GSM6204112

anotation_for_GSM6204112_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204112.csv", sep = "," )
dim(anotation_for_GSM6204112_obj2)
anotation_for_GSM6204112_obj2[1:6,]
colnames(anotation_for_GSM6204112_obj2)=anotation_for_GSM6204112_obj2[1,]
anotation_for_GSM6204112_obj2=anotation_for_GSM6204112_obj2[-1,]
length(unique(anotation_for_GSM6204112_obj2$ident))
length(anotatedAsFibroblastGSM6204112_obj2)
anotatedAsNK_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "NK Cells"),2]
anotatedAstrNK_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "trNK"),2]
anotatedAsCD8_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "CD4+ T Cells"),2]
anotatedAsTreg_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Treg"),2]
anotatedAsMac_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Macrophages"),2]
anotatedAsMast_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Mast Cells"),2]
anotatedAsDC_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "DC"),2]
anotatedAsB_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Plasma Cells"),2]
anotatedAsNeutrophil_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Neutrophil Cells"),2]
anotatedAsDuctalCells_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Ductal Cells"),2]
anotatedAsEntrocytes_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Entrocytes"),2]
anotatedAsFibroblastGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Fibroblasts"),2]
anotatedAsEndothelial_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Endothelial"),2]
anotatedAsAcinar_CellsGSM6204112_obj2=anotation_for_GSM6204112_obj2[which(anotation_for_GSM6204112_obj2$ident == "Acinar Cells"),2]


length(anotatedAsPlasma_CellsGSM6204112_obj2)
# Assigning Obj2
Idents(GSM6204112_obj2, cells = anotatedAsNK_CellsGSM6204112_obj2) <- "NK Cells"
Idents(GSM6204112_obj2, cells = anotatedAstrNK_CellsGSM6204112_obj2) <- "trNK Cells"
Idents(GSM6204112_obj2, cells = anotatedAsCD8_CellsGSM6204112_obj2) <- "CD8+ T Cells"
Idents(GSM6204112_obj2, cells = anotatedAsCD4_CellsGSM6204112_obj2) <- "CD4+ T Cells"
Idents(GSM6204112_obj2, cells = anotatedAsTreg_CellsGSM6204112_obj2) <- "Treg"
Idents(GSM6204112_obj2, cells = anotatedAsMac_CellsGSM6204112_obj2) <- "Macrophages"
Idents(GSM6204112_obj2, cells = anotatedAsMast_CellsGSM6204112_obj2) <- "Mast Cells"
Idents(GSM6204112_obj2, cells = anotatedAsDC_CellsGSM6204112_obj2) <- "DC"
Idents(GSM6204112_obj2, cells = anotatedAsB_CellsGSM6204112_obj2) <- "B Cells"
Idents(GSM6204112_obj2, cells = anotatedAsPlasma_CellsGSM6204112_obj2) <- "Plasma Cells"
Idents(GSM6204112_obj2, cells = anotatedAsNeutrophil_CellsGSM6204112_obj2) <- "Neutrophil Cells"
Idents(GSM6204112_obj2, cells = anotatedAsDuctalCells_CellsGSM6204112_obj2) <- "Ductal Cells"
Idents(GSM6204112_obj2, cells = anotatedAsEntrocytes_CellsGSM6204112_obj2) <- "Enterocytes"
Idents(GSM6204112_obj2, cells = anotatedAsFibroblastGSM6204112_obj2) <- "Fibroblasts"
Idents(GSM6204112_obj2, cells = anotatedAsEndothelial_CellsGSM6204112_obj2) <- "Endothelial"
Idents(GSM6204112_obj2, cells = anotatedAsAcinar_CellsGSM6204112_obj2) <- "Acinar Cells"


DimPlot(GSM6204112_obj2, reduction = "umap", label = T,repel = T)


#GSM6204113

anotation_for_GSM6204113_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204113.csv", sep = "," )
dim(anotation_for_GSM6204113_obj2)
anotation_for_GSM6204113_obj2[1:6,]
colnames(anotation_for_GSM6204113_obj2)=anotation_for_GSM6204113_obj2[1,]
anotation_for_GSM6204113_obj2=anotation_for_GSM6204113_obj2[-1,]
length(unique(anotation_for_GSM6204113_obj2$ident))
length(anotatedAsFibroblastGSM6204113_obj2)
anotatedAsNK_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "NK Cells"),2]
anotatedAstrNK_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "trNK"),2]
anotatedAsCD8_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "CD4+ T Cells"),2]
anotatedAsTreg_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "Treg"),2]
anotatedAsMac_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "Macrophages"),2]
anotatedAsDC_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "DC"),2]
anotatedAsMast_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "Mast Cells"),2]
anotatedAsB_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "Plasma Cells"),2]
anotatedAsDuctalCells_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "Ductal Cells"),2]
anotatedAsAcinar_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "Acinar Cells"),2]
anotatedAsEndothelial_CellsGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "Endothelial"),2]
anotatedAsFibroblastGSM6204113_obj2=anotation_for_GSM6204113_obj2[which(anotation_for_GSM6204113_obj2$ident == "Fibroblasts"),2]



length(anotatedAsPlasma_CellsGSM6204113_obj2)
# Assigning Obj2
Idents(GSM6204113_obj2, cells = anotatedAsNK_CellsGSM6204113_obj2) <- "NK Cells"
Idents(GSM6204113_obj2, cells = anotatedAstrNK_CellsGSM6204113_obj2) <- "trNK Cells"
Idents(GSM6204113_obj2, cells = anotatedAsCD8_CellsGSM6204113_obj2) <- "CD8+ T Cells"
Idents(GSM6204113_obj2, cells = anotatedAsCD4_CellsGSM6204113_obj2) <- "CD4+ T Cells"
Idents(GSM6204113_obj2, cells = anotatedAsTreg_CellsGSM6204113_obj2) <- "Treg"
Idents(GSM6204113_obj2, cells = anotatedAsMac_CellsGSM6204113_obj2) <- "Macrophages"
Idents(GSM6204113_obj2, cells = anotatedAsDC_CellsGSM6204113_obj2) <- "DC"
Idents(GSM6204113_obj2, cells = anotatedAsMast_CellsGSM6204113_obj2) <- "Mast Cells"
Idents(GSM6204113_obj2, cells = anotatedAsB_CellsGSM6204113_obj2) <- "B Cells"
Idents(GSM6204113_obj2, cells = anotatedAsPlasma_CellsGSM6204113_obj2) <- "Plasma Cells"
Idents(GSM6204113_obj2, cells = anotatedAsDuctalCells_CellsGSM6204113_obj2) <- "Ductal Cells"
Idents(GSM6204113_obj2, cells = anotatedAsAcinar_CellsGSM6204113_obj2) <- "Acinar Cells"
Idents(GSM6204113_obj2, cells = anotatedAsEndothelial_CellsGSM6204113_obj2) <- "Endothelial"
Idents(GSM6204113_obj2, cells = anotatedAsFibroblastGSM6204113_obj2) <- "Fibroblasts"



DimPlot(GSM6204113_obj2, reduction = "umap", label = T,repel = T)


#GSM6204114

anotation_for_GSM6204114_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204114.csv", sep = "," )
dim(anotation_for_GSM6204114_obj2)
anotation_for_GSM6204114_obj2[1:6,]
colnames(anotation_for_GSM6204114_obj2)=anotation_for_GSM6204114_obj2[1,]
anotation_for_GSM6204114_obj2=anotation_for_GSM6204114_obj2[-1,]
length(unique(anotation_for_GSM6204114_obj2$ident))
length(anotatedAsFibroblastGSM6204114_obj2)
anotatedAsNK_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "NK Cells"),2]
anotatedAstrNK_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "trNK"),2]
anotatedAsCD8_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "CD4+ T Cells"),2]
anotatedAsTreg_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Treg"),2]
anotatedAsB_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Plasma Cells"),2]
anotatedAsMast_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Mast Cells"),2]
anotatedAsMac_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Macrophages"),2]
anotatedAsDC_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "DC"),2]
anotatedAscDC_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "cDC"),2]
anotatedAsMonocytes_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Monocytes"),2]
anotatedAsDuctalCells_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Ductal Cells"),2]
anotatedAsAcinar_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Acinar Cells"),2]
anotatedAsEndothelial_CellsGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Endothelial"),2]
anotatedAsFibroblastGSM6204114_obj2=anotation_for_GSM6204114_obj2[which(anotation_for_GSM6204114_obj2$ident == "Fibroblasts"),2]


length(anotatedAsPlasma_CellsGSM6204114_obj2)
# Assigning Obj2
Idents(GSM6204114_obj2, cells = anotatedAsNK_CellsGSM6204114_obj2) <- "NK Cells"
Idents(GSM6204114_obj2, cells = anotatedAstrNK_CellsGSM6204114_obj2) <- "trNK Cells"
Idents(GSM6204114_obj2, cells = anotatedAsCD8_CellsGSM6204114_obj2) <- "CD8+ T Cells"
Idents(GSM6204114_obj2, cells = anotatedAsCD4_CellsGSM6204114_obj2) <- "CD4+ T Cells"
Idents(GSM6204114_obj2, cells = anotatedAsTreg_CellsGSM6204114_obj2) <- "Treg"
Idents(GSM6204114_obj2, cells = anotatedAsB_CellsGSM6204114_obj2) <- "B Cells"
Idents(GSM6204114_obj2, cells = anotatedAsPlasma_CellsGSM6204114_obj2) <- "Plasma Cells"
Idents(GSM6204114_obj2, cells = anotatedAsMast_CellsGSM6204114_obj2) <- "Mast Cells"
Idents(GSM6204114_obj2, cells = anotatedAsMac_CellsGSM6204114_obj2) <- "Macrophages"
Idents(GSM6204114_obj2, cells = anotatedAsDC_CellsGSM6204114_obj2) <- "DC"
Idents(GSM6204114_obj2, cells = anotatedAscDC_CellsGSM6204114_obj2) <- "cDC"
Idents(GSM6204114_obj2, cells = anotatedAsMonocytes_CellsGSM6204114_obj2) <- "Monocytes"
Idents(GSM6204114_obj2, cells = anotatedAsDuctalCells_CellsGSM6204114_obj2) <- "Ductal Cells"
Idents(GSM6204114_obj2, cells = anotatedAsAcinar_CellsGSM6204114_obj2) <- "Acinar Cells"
Idents(GSM6204114_obj2, cells = anotatedAsEndothelial_CellsGSM6204114_obj2) <- "Endothelial"
Idents(GSM6204114_obj2, cells = anotatedAsFibroblastGSM6204114_obj2) <- "Fibroblasts"


DimPlot(GSM6204114_obj2, reduction = "umap", label = T,repel = T)



#GSM6204115

anotation_for_GSM6204115_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204115.csv", sep = "," )
dim(anotation_for_GSM6204115_obj2)
anotation_for_GSM6204115_obj2[1:6,]
colnames(anotation_for_GSM6204115_obj2)=anotation_for_GSM6204115_obj2[1,]
anotation_for_GSM6204115_obj2=anotation_for_GSM6204115_obj2[-1,]
length(unique(anotation_for_GSM6204115_obj2$ident))
length(anotatedAsFibroblastGSM6204115_obj2)
anotatedAsDuctalCells_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "Ductal Cells"),2]
anotatedAsFibroblastGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "Fibroblasts"),2]
anotatedAsEndothelial_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "Endothelial"),2]
anotatedAsMac_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "Macrophages"),2]
anotatedAsCD8_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "CD4+ T Cells"),2]
anotatedAsNK_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "NK Cells"),2]
anotatedAstrNK_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "trNK"),2]
anotatedAsMast_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "Mast Cells"),2]
anotatedAsTreg_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "Treg"),2]
anotatedAsB_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "Plasma Cells"),2]
anotatedAsPericyte_CellsGSM6204115_obj2=anotation_for_GSM6204115_obj2[which(anotation_for_GSM6204115_obj2$ident == "Pericyte"),2]



length(anotatedAsPlasma_CellsGSM6204115_obj2)
# Assigning Obj2
Idents(GSM6204115_obj2, cells = anotatedAsDuctalCells_CellsGSM6204115_obj2) <- "Ductal Cells"
Idents(GSM6204115_obj2, cells = anotatedAsFibroblastGSM6204115_obj2) <- "Fibroblasts"
Idents(GSM6204115_obj2, cells = anotatedAsEndothelial_CellsGSM6204115_obj2) <- "Endothelial"
Idents(GSM6204115_obj2, cells = anotatedAsMac_CellsGSM6204115_obj2) <- "Macrophages"
Idents(GSM6204115_obj2, cells = anotatedAsCD8_CellsGSM6204115_obj2) <- "CD8+ T Cells"
Idents(GSM6204115_obj2, cells = anotatedAsCD4_CellsGSM6204115_obj2) <- "CD4+ T Cells"
Idents(GSM6204115_obj2, cells = anotatedAsNK_CellsGSM6204115_obj2) <- "NK Cells"
Idents(GSM6204115_obj2, cells = anotatedAstrNK_CellsGSM6204115_obj2) <- "trNK Cells"
Idents(GSM6204115_obj2, cells = anotatedAsMast_CellsGSM6204115_obj2) <- "Mast Cells"
Idents(GSM6204115_obj2, cells = anotatedAsTreg_CellsGSM6204115_obj2) <- "Treg"
Idents(GSM6204115_obj2, cells = anotatedAsB_CellsGSM6204115_obj2) <- "B Cells"
Idents(GSM6204115_obj2, cells = anotatedAsPlasma_CellsGSM6204115_obj2) <- "Plasma Cells"
Idents(GSM6204115_obj2, cells = anotatedAsPericyte_CellsGSM6204115_obj2) <- "Pericytes"



DimPlot(GSM6204115_obj2, reduction = "umap", label = T,repel = T)



#GSM6204117

anotation_for_GSM6204117_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204117.csv", sep = "," )
dim(anotation_for_GSM6204117_obj2)
anotation_for_GSM6204117_obj2[1:6,]
colnames(anotation_for_GSM6204117_obj2)=anotation_for_GSM6204117_obj2[1,]
anotation_for_GSM6204117_obj2=anotation_for_GSM6204117_obj2[-1,]
length(unique(anotation_for_GSM6204117_obj2$ident))
length(anotatedAsFibroblastGSM6204117_obj2)
anotatedAsMast_CellsGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "Mast Cells"),2]
anotatedAsFibroblastGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "Fibroblasts"),2]
anotatedAsTreg_CellsGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "Treg"),2]
anotatedAsDuctalCells_CellsGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "Ductal Cells"),2]
anotatedAsCD4_CellsGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "CD4+ T Cells"),2]
anotatedAsMac_CellsGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "Macrophages"),2]
anotatedAsCD8_CellsGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "CD8+ T Cells"),2]
anotatedAsNeutrophil_CellsGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "Neutrophil Cells"),2]
anotatedAsEndothelial_CellsGSM6204117_obj2=anotation_for_GSM6204117_obj2[which(anotation_for_GSM6204117_obj2$ident == "Endothelial"),2]




length(anotatedAsPlasma_CellsGSM6204117_obj2)
# Assigning Obj2
Idents(GSM6204117_obj2, cells = anotatedAsMast_CellsGSM6204117_obj2) <- "Mast Cells"
Idents(GSM6204117_obj2, cells = anotatedAsFibroblastGSM6204117_obj2) <- "Fibroblasts"
Idents(GSM6204117_obj2, cells = anotatedAsTreg_CellsGSM6204117_obj2) <- "Treg"
Idents(GSM6204117_obj2, cells = anotatedAsDuctalCells_CellsGSM6204117_obj2) <- "Ductal Cells"
Idents(GSM6204117_obj2, cells = anotatedAsCD4_CellsGSM6204117_obj2) <- "CD4+ T Cells"
Idents(GSM6204117_obj2, cells = anotatedAsMac_CellsGSM6204117_obj2) <- "Macrophages"
Idents(GSM6204117_obj2, cells = anotatedAsCD8_CellsGSM6204117_obj2) <- "CD8+ T Cells"
Idents(GSM6204117_obj2, cells = anotatedAsNeutrophil_CellsGSM6204117_obj2) <- "Neutrophil Cells"
Idents(GSM6204117_obj2, cells = anotatedAsEndothelial_CellsGSM6204117_obj2) <- "Endothelial"



DimPlot(GSM6204117_obj2, reduction = "umap", label = T,repel = T)


#GSM6204118

anotation_for_GSM6204118_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204118.csv", sep = "," )
dim(anotation_for_GSM6204118_obj2)
anotation_for_GSM6204118_obj2[1:6,]
colnames(anotation_for_GSM6204118_obj2)=anotation_for_GSM6204118_obj2[1,]
anotation_for_GSM6204118_obj2=anotation_for_GSM6204118_obj2[-1,]
length(unique(anotation_for_GSM6204118_obj2$ident))
length(anotatedAsFibroblastGSM6204118_obj2)
anotatedAsDC_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "DC"),2]
anotatedAsPlasma_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "Plasma Cells"),2]
anotatedAsCD8_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "CD4+ T Cells"),2]
anotatedAsNK_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "NK Cells"),2]
anotatedAsTreg_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "Treg"),2]
anotatedAsDuctalCells_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "Ductal Cells"),2]
anotatedAsFibroblastGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "Fibroblasts"),2]
anotatedAsMac_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "Macrophages"),2]
anotatedAsEndothelial_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "Endothelial"),2]
anotatedAstrNK_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "trNK"),2]
anotatedAsB_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "B Cells"),2]
anotatedAsEndocrine_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "Endocrine Cells"),2]
anotatedAsMast_CellsGSM6204118_obj2=anotation_for_GSM6204118_obj2[which(anotation_for_GSM6204118_obj2$ident == "Mast Cells"),2]


length(anotatedAsPlasma_CellsGSM6204118_obj2)
# Assigning Obj2
Idents(GSM6204118_obj2, cells = anotatedAsDC_CellsGSM6204118_obj2) <- "DC"
Idents(GSM6204118_obj2, cells = anotatedAsPlasma_CellsGSM6204118_obj2) <- "Plasma Cells"
Idents(GSM6204118_obj2, cells = anotatedAsCD8_CellsGSM6204118_obj2) <- "CD8+ T Cells"
Idents(GSM6204118_obj2, cells = anotatedAsCD4_CellsGSM6204118_obj2) <- "CD4+ T Cells"
Idents(GSM6204118_obj2, cells = anotatedAsNK_CellsGSM6204118_obj2) <- "NK Cells"
Idents(GSM6204118_obj2, cells = anotatedAsTreg_CellsGSM6204118_obj2) <- "Treg"
Idents(GSM6204118_obj2, cells = anotatedAsDuctalCells_CellsGSM6204118_obj2) <- "Ductal Cells"
Idents(GSM6204118_obj2, cells = anotatedAsFibroblastGSM6204118_obj2) <- "Fibroblasts"
Idents(GSM6204118_obj2, cells = anotatedAsMac_CellsGSM6204118_obj2) <- "Macrophages"
Idents(GSM6204118_obj2, cells = anotatedAsEndothelial_CellsGSM6204118_obj2) <- "Endothelial"
Idents(GSM6204118_obj2, cells = anotatedAstrNK_CellsGSM6204118_obj2) <- "trNK Cells"
Idents(GSM6204118_obj2, cells = anotatedAsB_CellsGSM6204118_obj2) <- "B Cells"
Idents(GSM6204118_obj2, cells = anotatedAsEndocrine_CellsGSM6204118_obj2) <- "Endocrine Cells"
Idents(GSM6204118_obj2, cells = anotatedAsMast_CellsGSM6204118_obj2) <- "Mast Cells"



DimPlot(GSM6204118_obj2, reduction = "umap", label = T,repel = T)



#GSM6204119

anotation_for_GSM6204119_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204119.csv", sep = "," )
dim(anotation_for_GSM6204119_obj2)
anotation_for_GSM6204119_obj2[1:6,]
colnames(anotation_for_GSM6204119_obj2)=anotation_for_GSM6204119_obj2[1,]
anotation_for_GSM6204119_obj2=anotation_for_GSM6204119_obj2[-1,]
length(unique(anotation_for_GSM6204119_obj2$ident))
length(anotatedAsFibroblastGSM6204119_obj2)
anotatedAsPlasma_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Plasma Cells"),2]
anotatedAsMast_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Mast Cells"),2]
anotatedAsMac_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Macrophages"),2]
anotatedAsEndothelial_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Endothelial"),2]
anotatedAsHepatocytes_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Hepatocytes"),2]
anotatedAsFibroblastGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Fibroblasts"),2]
anotatedAsB_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "B Cells"),2]
anotatedAsNK_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "NK Cells"),2]
anotatedAsDuctalCells_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Ductal Cells"),2]
anotatedAsCD4_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "CD4+ T Cells"),2]
anotatedAsCD8_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "CD8+ T Cells"),2]
anotatedAsEndocrine_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Endocrine Cells"),2]
anotatedAsTreg_CellsGSM6204119_obj2=anotation_for_GSM6204119_obj2[which(anotation_for_GSM6204119_obj2$ident == "Treg"),2]



length(anotatedAsPlasma_CellsGSM6204119_obj2)
# Assigning Obj2
Idents(GSM6204119_obj2, cells = anotatedAsPlasma_CellsGSM6204119_obj2) <- "Plasma Cells"
Idents(GSM6204119_obj2, cells = anotatedAsMast_CellsGSM6204119_obj2) <- "Mast Cells"
Idents(GSM6204119_obj2, cells = anotatedAsMac_CellsGSM6204119_obj2) <- "Macrophages"
Idents(GSM6204119_obj2, cells = anotatedAsEndothelial_CellsGSM6204119_obj2) <- "Endothelial"
Idents(GSM6204119_obj2, cells = anotatedAsHepatocytes_CellsGSM6204119_obj2) <- "Hepatocytes"
Idents(GSM6204119_obj2, cells = anotatedAsFibroblastGSM6204119_obj2) <- "Fibroblasts"
Idents(GSM6204119_obj2, cells = anotatedAsB_CellsGSM6204119_obj2) <- "B Cells"
Idents(GSM6204119_obj2, cells = anotatedAsNK_CellsGSM6204119_obj2) <- "NK Cells"
Idents(GSM6204119_obj2, cells = anotatedAsDuctalCells_CellsGSM6204119_obj2) <- "Ductal Cells"
Idents(GSM6204119_obj2, cells = anotatedAsCD4_CellsGSM6204119_obj2) <- "CD4+ T Cells"
Idents(GSM6204119_obj2, cells = anotatedAsCD8_CellsGSM6204119_obj2) <- "CD8+ T Cells"
Idents(GSM6204119_obj2, cells = anotatedAsEndocrine_CellsGSM6204119_obj2) <- "Endocrine Cells"
Idents(GSM6204119_obj2, cells = anotatedAsTreg_CellsGSM6204119_obj2) <- "Treg"


#GSM6204120

anotation_for_GSM6204120_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204120.csv", sep = "," )
dim(anotation_for_GSM6204120_obj2)
anotation_for_GSM6204120_obj2[1:6,]
colnames(anotation_for_GSM6204120_obj2)=anotation_for_GSM6204120_obj2[1,]
anotation_for_GSM6204120_obj2=anotation_for_GSM6204120_obj2[-1,]
length(unique(anotation_for_GSM6204120_obj2$ident))
length(anotatedAsFibroblastGSM6204120_obj2)
anotatedAsEndothelial_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "Endothelial"),2]
anotatedAsMast_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "Mast Cells"),2]
anotatedAsB_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "Plasma Cells"),2]
anotatedAsDuctalCells_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "Ductal Cells"),2]
anotatedAsMac_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "Macrophages"),2]
anotatedAsCD4_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "CD4+ T Cells"),2]
anotatedAsNK_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "NK Cells"),2]
anotatedAsCD8_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "CD8+ T Cells"),2]
anotatedAsFibroblastGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "Fibroblasts"),2]
anotatedAsNeutrophil_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "Neutrophil Cells"),2]
anotatedAsTreg_CellsGSM6204120_obj2=anotation_for_GSM6204120_obj2[which(anotation_for_GSM6204120_obj2$ident == "Treg"),2]




length(anotatedAsPlasma_CellsGSM6204120_obj2)
# Assigning Obj2
Idents(GSM6204120_obj2, cells = anotatedAsEndothelial_CellsGSM6204120_obj2) <- "Endothelial"
Idents(GSM6204120_obj2, cells = anotatedAsMast_CellsGSM6204120_obj2) <- "Mast Cells"
Idents(GSM6204120_obj2, cells = anotatedAsB_CellsGSM6204120_obj2) <- "B Cells"
Idents(GSM6204120_obj2, cells = anotatedAsPlasma_CellsGSM6204120_obj2) <- "Plasma Cells"
Idents(GSM6204120_obj2, cells = anotatedAsDuctalCells_CellsGSM6204120_obj2) <- "Ductal Cells"
Idents(GSM6204120_obj2, cells = anotatedAsMac_CellsGSM6204120_obj2) <- "Macrophages"
Idents(GSM6204120_obj2, cells = anotatedAsCD4_CellsGSM6204120_obj2) <- "CD4+ T Cells"
Idents(GSM6204120_obj2, cells = anotatedAsNK_CellsGSM6204120_obj2) <- "NK Cells"
Idents(GSM6204120_obj2, cells = anotatedAsCD8_CellsGSM6204120_obj2) <- "CD8+ T Cells"
Idents(GSM6204120_obj2, cells = anotatedAsFibroblastGSM6204120_obj2) <- "Fibroblasts"
Idents(GSM6204120_obj2, cells = anotatedAsNeutrophil_CellsGSM6204120_obj2) <- "Neutrophil Cells"
Idents(GSM6204120_obj2, cells = anotatedAsTreg_CellsGSM6204120_obj2) <- "Treg"





#GSM6204121

anotation_for_GSM6204121_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204121.csv", sep = "," )
dim(anotation_for_GSM6204121_obj2)
anotation_for_GSM6204121_obj2[1:6,]
colnames(anotation_for_GSM6204121_obj2)=anotation_for_GSM6204121_obj2[1,]
anotation_for_GSM6204121_obj2=anotation_for_GSM6204121_obj2[-1,]
unique(anotation_for_GSM6204121_obj2$ident)
length(anotatedAsFibroblastGSM6204121_obj2)
anotatedAsMast_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "Mast Cells"),2]
anotatedAsNK_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "NK Cells"),2]
anotatedAsNeutrophil_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "Neutrophil Cells"),2]
anotatedAsTreg_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "Treg"),2]
anotatedAsEndothelial_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "Endothelial"),2]
anotatedAsB_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "Plasma Cells"),2]
anotatedAsDuctalCells_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "Ductal Cells"),2]
anotatedAsCD4_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "CD4+ T Cells"),2]
anotatedAsMac_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "Macrophages"),2]
anotatedAsFibroblastGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "Fibroblasts"),2]
anotatedAsCD8_CellsGSM6204121_obj2=anotation_for_GSM6204121_obj2[which(anotation_for_GSM6204121_obj2$ident == "CD8+ T Cells"),2]




length(anotatedAsPlasma_CellsGSM6204121_obj2)
# Assigning Obj2
Idents(GSM6204121_obj2, cells = anotatedAsMast_CellsGSM6204121_obj2) <- "Mast Cells"
Idents(GSM6204121_obj2, cells = anotatedAsNK_CellsGSM6204121_obj2) <- "NK Cells"
Idents(GSM6204121_obj2, cells = anotatedAsNeutrophil_CellsGSM6204121_obj2) <- "Neutrophil Cells"
Idents(GSM6204121_obj2, cells = anotatedAsTreg_CellsGSM6204121_obj2) <- "Treg"
Idents(GSM6204121_obj2, cells = anotatedAsEndothelial_CellsGSM6204121_obj2) <- "Endothelial"
Idents(GSM6204121_obj2, cells = anotatedAsB_CellsGSM6204121_obj2) <- "B Cells"
Idents(GSM6204121_obj2, cells = anotatedAsPlasma_CellsGSM6204121_obj2) <- "Plasma Cells"
Idents(GSM6204121_obj2, cells = anotatedAsDuctalCells_CellsGSM6204121_obj2) <- "Ductal Cells"
Idents(GSM6204121_obj2, cells = anotatedAsCD4_CellsGSM6204121_obj2) <- "CD4+ T Cells"
Idents(GSM6204121_obj2, cells = anotatedAsMac_CellsGSM6204121_obj2) <- "Macrophages"
Idents(GSM6204121_obj2, cells = anotatedAsFibroblastGSM6204121_obj2) <- "Fibroblasts"
Idents(GSM6204121_obj2, cells = anotatedAsCD8_CellsGSM6204121_obj2) <- "CD8+ T Cells"



#GSM6204122

anotation_for_GSM6204122_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/anotaionGSM6204122.csv", sep = "," )
dim(anotation_for_GSM6204122_obj2)
anotation_for_GSM6204122_obj2[1:6,]
colnames(anotation_for_GSM6204122_obj2)=anotation_for_GSM6204122_obj2[1,]
anotation_for_GSM6204122_obj2=anotation_for_GSM6204122_obj2[-1,]
unique(anotation_for_GSM6204122_obj2$ident)
length(anotatedAsFibroblastGSM6204122_obj2)
anotatedAsFibroblastGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Fibroblasts"),2]
anotatedAsCD8_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "CD4+ T Cells"),2]
anotatedAsEndothelial_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Endothelial"),2]
anotatedAsMac_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Macrophages"),2]
anotatedAsDuctalCells_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Ductal Cells"),2]
anotatedAsB_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "B Cells"),2]
anotatedAsErythroblast_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Erythroblast"),2]
anotatedAsPlasma_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Plasma Cells"),2]
anotatedAsNK_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "NK Cells"),2]
anotatedAsNeutrophil_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Neutrophil Cells"),2]
anotatedAsMast_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Mast Cells"),2]
anotatedAsDC_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "DC"),2]
anotatedAsAcinarCells_CellsGSM6204122_obj2=anotation_for_GSM6204122_obj2[which(anotation_for_GSM6204122_obj2$ident == "Acinar Cells"),2]





length(anotatedAsPlasma_CellsGSM6204122_obj2)
# Assigning Obj2
Idents(GSM6204122_obj2, cells = anotatedAsFibroblastGSM6204122_obj2) <- "Fibroblasts"
Idents(GSM6204122_obj2, cells = anotatedAsCD8_CellsGSM6204122_obj2) <- "CD8+ T Cells"
Idents(GSM6204122_obj2, cells = anotatedAsCD4_CellsGSM6204122_obj2) <- "CD4+ T Cells"
Idents(GSM6204122_obj2, cells = anotatedAsEndothelial_CellsGSM6204122_obj2) <- "Endothelial"
Idents(GSM6204122_obj2, cells = anotatedAsMac_CellsGSM6204122_obj2) <- "Macrophages"
Idents(GSM6204122_obj2, cells = anotatedAsDuctalCells_CellsGSM6204122_obj2) <- "Ductal Cells"
Idents(GSM6204122_obj2, cells = anotatedAsB_CellsGSM6204122_obj2) <- "B Cells"
Idents(GSM6204122_obj2, cells = anotatedAsErythroblast_CellsGSM6204122_obj2) <- "Erythroblast"
Idents(GSM6204122_obj2, cells = anotatedAsPlasma_CellsGSM6204122_obj2) <- "Plasma Cells"
Idents(GSM6204122_obj2, cells = anotatedAsNK_CellsGSM6204122_obj2) <- "NK Cells"
Idents(GSM6204122_obj2, cells = anotatedAsNeutrophil_CellsGSM6204122_obj2) <- "Neutrophil Cells"
Idents(GSM6204122_obj2, cells = anotatedAsMast_CellsGSM6204122_obj2) <- "Mast Cells"
Idents(GSM6204122_obj2, cells = anotatedAsAcinar_CellsGSM6204122_obj2) <- "Acinar Cells"
Idents(GSM6204122_obj2, cells = anotatedAsDC_CellsGSM6204122_obj2) <- "DC"




#GSM6204123

anotation_for_GSM6204123_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204123.csv", sep = "," )
dim(anotation_for_GSM6204123_obj2)
anotation_for_GSM6204123_obj2[1:6,]
colnames(anotation_for_GSM6204123_obj2)=anotation_for_GSM6204123_obj2[1,]
anotation_for_GSM6204123_obj2=anotation_for_GSM6204123_obj2[-1,]
unique(anotation_for_GSM6204123_obj2$ident)
length(anotatedAsFibroblastGSM6204123_obj2)
anotatedAstrNK_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "trNK"),2]
anotatedAsNK_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "NK Cells"),2]
anotatedAsTreg_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Treg"),2]
anotatedAsFibroblastGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Fibroblasts"),2]
anotatedAsCD8_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "CD4+ T Cells"),2]
anotatedAsEndothelial_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Endothelial"),2]
anotatedAsMac_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Macrophages"),2]
anotatedAsDuctalCells_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Ductal Cells"),2]
anotatedAsB_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "B Cells"),2]
anotatedAsErythroblast_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Erythroblast"),2]
anotatedAsPlasma_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Plasma Cells"),2]
anotatedAsNeutrophil_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Neutrophil Cells"),2]
anotatedAsMast_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Mast Cells"),2]
anotatedAsDC_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "DC"),2]
anotatedAsAcinar_CellsGSM6204123_obj2=anotation_for_GSM6204123_obj2[which(anotation_for_GSM6204123_obj2$ident == "Acinar Cells"),2]





length(anotatedAsPlasma_CellsGSM6204123_obj2)
# Assigning Obj2
Idents(GSM6204123_obj2, cells = anotatedAsFibroblastGSM6204123_obj2) <- "Fibroblasts"
Idents(GSM6204123_obj2, cells = anotatedAsCD8_CellsGSM6204123_obj2) <- "CD8+ T Cells"
Idents(GSM6204123_obj2, cells = anotatedAsCD4_CellsGSM6204123_obj2) <- "CD4+ T Cells"
Idents(GSM6204123_obj2, cells = anotatedAsEndothelial_CellsGSM6204123_obj2) <- "Endothelial"
Idents(GSM6204123_obj2, cells = anotatedAsMac_CellsGSM6204123_obj2) <- "Macrophages"
Idents(GSM6204123_obj2, cells = anotatedAsDuctalCells_CellsGSM6204123_obj2) <- "Ductal Cells"
Idents(GSM6204123_obj2, cells = anotatedAsB_CellsGSM6204123_obj2) <- "B Cells"
Idents(GSM6204123_obj2, cells = anotatedAsErythroblast_CellsGSM6204123_obj2) <- "Erythroblast"
Idents(GSM6204123_obj2, cells = anotatedAsPlasma_CellsGSM6204123_obj2) <- "Plasma Cells"
Idents(GSM6204123_obj2, cells = anotatedAsNK_CellsGSM6204123_obj2) <- "NK Cells"
Idents(GSM6204123_obj2, cells = anotatedAsNeutrophil_CellsGSM6204123_obj2) <- "Neutrophil Cells"
Idents(GSM6204123_obj2, cells = anotatedAsMast_CellsGSM6204123_obj2) <- "Mast Cells"
Idents(GSM6204123_obj2, cells = anotatedAsAcinar_CellsGSM6204123_obj2) <- "Acinar Cells"
Idents(GSM6204123_obj2, cells = anotatedAsDC_CellsGSM6204123_obj2) <- "DC"



#GSM6204124

anotation_for_GSM6204124_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204124.csv", sep = "," )
dim(anotation_for_GSM6204124_obj2)
anotation_for_GSM6204124_obj2[1:6,]
colnames(anotation_for_GSM6204124_obj2)=anotation_for_GSM6204124_obj2[1,]
anotation_for_GSM6204124_obj2=anotation_for_GSM6204124_obj2[-1,]
unique(anotation_for_GSM6204124_obj2$ident)
length(anotatedAsFibroblastGSM6204124_obj2)
anotatedAsTreg_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Treg"),2]
anotatedAsFibroblastGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Fibroblasts"),2]
anotatedAstrNK_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "trNK"),2]
anotatedAsMast_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Mast Cells"),2]
anotatedAsCD4_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "CD4+ T Cells"),2]
anotatedAsCD8_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "CD8+ T Cells"),2]
anotatedAsDuctalCells_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Ductal Cells"),2]
anotatedAsNK_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "NK Cells"),2]
anotatedAsMac_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Macrophages"),2]
anotatedAsNeutrophil_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Neutrophil Cells"),2]
anotatedAsEndothelial_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Endothelial"),2]
anotatedAsPlasma_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Plasma Cells"),2]
anotatedAsHepatocytes_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "Hepatocytes"),2]
anotatedAsB_CellsGSM6204124_obj2=anotation_for_GSM6204124_obj2[which(anotation_for_GSM6204124_obj2$ident == "B Cells"),2]






length(anotatedAsPlasma_CellsGSM6204124_obj2)
# Assigning Obj2
Idents(GSM6204124_obj2, cells = anotatedAsTreg_CellsGSM6204124_obj2) <- "Treg"
Idents(GSM6204124_obj2, cells = anotatedAsFibroblastGSM6204124_obj2) <- "Fibroblasts"
Idents(GSM6204124_obj2, cells = anotatedAstrNK_CellsGSM6204124_obj2) <- "trNK Cells"
Idents(GSM6204124_obj2, cells = anotatedAsMast_CellsGSM6204124_obj2) <- "Mast Cells"
Idents(GSM6204124_obj2, cells = anotatedAsCD4_CellsGSM6204124_obj2) <- "CD4+ T Cells"
Idents(GSM6204124_obj2, cells = anotatedAsCD8_CellsGSM6204124_obj2) <- "CD8+ T Cells"
Idents(GSM6204124_obj2, cells = anotatedAsDuctalCells_CellsGSM6204124_obj2) <- "Ductal Cells"
Idents(GSM6204124_obj2, cells = anotatedAsNK_CellsGSM6204124_obj2) <- "NK Cells"
Idents(GSM6204124_obj2, cells = anotatedAsMac_CellsGSM6204124_obj2) <- "Macrophages"
Idents(GSM6204124_obj2, cells = anotatedAsNeutrophil_CellsGSM6204124_obj2) <- "Neutrophil Cells"
Idents(GSM6204124_obj2, cells = anotatedAsEndothelial_CellsGSM6204124_obj2) <- "Endothelial"
Idents(GSM6204124_obj2, cells = anotatedAsPlasma_CellsGSM6204124_obj2) <- "Plasma Cells"
Idents(GSM6204124_obj2, cells = anotatedAsHepatocytes_CellsGSM6204124_obj2) <- "Hepatocytes"
Idents(GSM6204124_obj2, cells = anotatedAsB_CellsGSM6204124_obj2) <- "B Cells"



#GSM6204125

anotation_for_GSM6204125_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204125.csv", sep = "," )
dim(anotation_for_GSM6204125_obj2)
anotation_for_GSM6204125_obj2[1:6,]
colnames(anotation_for_GSM6204125_obj2)=anotation_for_GSM6204125_obj2[1,]
anotation_for_GSM6204125_obj2=anotation_for_GSM6204125_obj2[-1,]
unique(anotation_for_GSM6204125_obj2$ident)
length(anotatedAsFibroblastGSM6204125_obj2)
anotatedAsTreg_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "Treg"),2]
anotatedAsB_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "B Cells"),2]
anotatedAsDuctalCells_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "Ductal Cells"),2]
anotatedAsNeutrophil_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "Neutrophil Cells"),2]
anotatedAsCD4_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "CD4+ T Cells"),2]
anotatedAsCD8_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "CD8+ T Cells"),2]
anotatedAsMac_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "Macrophages"),2]
anotatedAsNK_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "NK Cells"),2]
anotatedAsFibroblastGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "Fibroblasts"),2]
anotatedAsHepatocytes_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "Hepatocytes"),2]
anotatedAsDC_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "DC"),2]
anotatedAstrNK_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "trNK"),2]
anotatedAsEndothelial_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "Endothelial"),2]
anotatedAsPlasma_CellsGSM6204125_obj2=anotation_for_GSM6204125_obj2[which(anotation_for_GSM6204125_obj2$ident == "Plasma Cells"),2]






length(anotatedAsPlasma_CellsGSM6204125_obj2)
# Assigning Obj2
Idents(GSM6204125_obj2, cells = anotatedAsTreg_CellsGSM6204125_obj2) <- "Treg"
Idents(GSM6204125_obj2, cells = anotatedAsB_CellsGSM6204125_obj2) <- "B Cells"
Idents(GSM6204125_obj2, cells = anotatedAsDuctalCells_CellsGSM6204125_obj2) <- "Ductal Cells"
Idents(GSM6204125_obj2, cells = anotatedAsNeutrophil_CellsGSM6204125_obj2) <- "Neutrophil Cells"
Idents(GSM6204125_obj2, cells = anotatedAsCD4_CellsGSM6204125_obj2) <- "CD4+ T Cells"
Idents(GSM6204125_obj2, cells = anotatedAsCD8_CellsGSM6204125_obj2) <- "CD8+ T Cells"
Idents(GSM6204125_obj2, cells = anotatedAsMac_CellsGSM6204125_obj2) <- "Macrophages"
Idents(GSM6204125_obj2, cells = anotatedAsNK_CellsGSM6204125_obj2) <- "NK Cells"
Idents(GSM6204125_obj2, cells = anotatedAsFibroblastGSM6204125_obj2) <- "Fibroblasts"
Idents(GSM6204125_obj2, cells = anotatedAsHepatocytes_CellsGSM6204125_obj2) <- "Hepatocytes"
Idents(GSM6204125_obj2, cells = anotatedAsDC_CellsGSM6204125_obj2) <- "DC"
Idents(GSM6204125_obj2, cells = anotatedAstrNK_CellsGSM6204125_obj2) <- "trNK Cells"
Idents(GSM6204125_obj2, cells = anotatedAsEndothelial_CellsGSM6204125_obj2) <- "Endothelial"
Idents(GSM6204125_obj2, cells = anotatedAsPlasma_CellsGSM6204125_obj2) <- "Plasma Cells"




#GSM6204127

anotation_for_GSM6204127_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204127.csv", sep = "," )
dim(anotation_for_GSM6204127_obj2)
anotation_for_GSM6204127_obj2[1:6,]
colnames(anotation_for_GSM6204127_obj2)=anotation_for_GSM6204127_obj2[1,]
anotation_for_GSM6204127_obj2=anotation_for_GSM6204127_obj2[-1,]
unique(anotation_for_GSM6204127_obj2$ident)
length(anotatedAsFibroblastGSM6204127_obj2)
anotatedAstrNK_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "trNK"),2]
anotatedAsTreg_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "Treg"),2]
anotatedAsCD4_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "CD4+ T Cells"),2]
anotatedAsB_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "B Cells"),2]
anotatedAsCD8_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "CD8+ T Cells"),2]
anotatedAsDuctalCells_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "Ductal Cells"),2]
anotatedAsMac_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "Macrophages"),2]
anotatedAsNK_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "NK Cells"),2]
anotatedAsFibroblastGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "Fibroblasts"),2]
anotatedAsPlasma_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "Plasma Cells"),2]
anotatedAsEndothelial_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "Endothelial"),2]
anotatedAsDC_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "DC"),2]
anotatedAsAcinar_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "Acinar Cells"),2]
anotatedAsMast_CellsGSM6204127_obj2=anotation_for_GSM6204127_obj2[which(anotation_for_GSM6204127_obj2$ident == "Mast Cells"),2]






length(anotatedAsPlasma_CellsGSM6204127_obj2)
# Assigning Obj2
Idents(GSM6204127_obj2, cells = anotatedAstrNK_CellsGSM6204127_obj2) <- "trNK Cells"
Idents(GSM6204127_obj2, cells = anotatedAsTreg_CellsGSM6204127_obj2) <- "Treg"
Idents(GSM6204127_obj2, cells = anotatedAsCD4_CellsGSM6204127_obj2) <- "CD4+ T Cells"
Idents(GSM6204127_obj2, cells = anotatedAsB_CellsGSM6204127_obj2) <- "B Cells"
Idents(GSM6204127_obj2, cells = anotatedAsCD8_CellsGSM6204127_obj2) <- "CD8+ T Cells"
Idents(GSM6204127_obj2, cells = anotatedAsDuctalCells_CellsGSM6204127_obj2) <- "Ductal Cells"
Idents(GSM6204127_obj2, cells = anotatedAsMac_CellsGSM6204127_obj2) <- "Macrophages"
Idents(GSM6204127_obj2, cells = anotatedAsNK_CellsGSM6204127_obj2) <- "NK Cells"
Idents(GSM6204127_obj2, cells = anotatedAsFibroblastGSM6204127_obj2) <- "Fibroblasts"
Idents(GSM6204127_obj2, cells = anotatedAsPlasma_CellsGSM6204127_obj2) <- "Plasma Cells"
Idents(GSM6204127_obj2, cells = anotatedAsEndothelial_CellsGSM6204127_obj2) <- "Endothelial"
Idents(GSM6204127_obj2, cells = anotatedAsDC_CellsGSM6204127_obj2) <- "DC"
Idents(GSM6204127_obj2, cells = anotatedAsAcinar_CellsGSM6204127_obj2) <- "Acinar Cells"
Idents(GSM6204127_obj2, cells = anotatedAsMast_CellsGSM6204127_obj2) <- "Mast Cells"


#GSM6204128

anotation_for_GSM6204128_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204128.csv", sep = "," )
dim(anotation_for_GSM6204128_obj2)
anotation_for_GSM6204128_obj2[1:6,]
colnames(anotation_for_GSM6204128_obj2)=anotation_for_GSM6204128_obj2[1,]
anotation_for_GSM6204128_obj2=anotation_for_GSM6204128_obj2[-1,]
unique(anotation_for_GSM6204128_obj2$ident)
length(anotatedAsFibroblastGSM6204128_obj2)
anotatedAstrNK_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "trNK"),2]
anotatedAsTreg_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Treg"),2]
anotatedAsMast_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Mast Cells"),2]
anotatedAsCD4_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "CD4+ T Cells"),2]
anotatedAsCD8_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "CD8+ T Cells"),2]
anotatedAsDuctalCells_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Ductal Cells"),2]
anotatedAsMac_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Macrophages"),2]
anotatedAsNK_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "NK Cells"),2]
anotatedAsB_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "B Cells"),2]
anotatedAsAcinar_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Acinar Cells"),2]
anotatedAsFibroblastGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Fibroblasts"),2]
anotatedAsEndothelial_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Endothelial"),2]
anotatedAsPlasma_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Plasma Cells"),2]
anotatedAsGastricEpithelial_CellsGSM6204128_obj2=anotation_for_GSM6204128_obj2[which(anotation_for_GSM6204128_obj2$ident == "Gastric Epithelial"),2]






length(anotatedAsPlasma_CellsGSM6204128_obj2)
# Assigning Obj2
Idents(GSM6204128_obj2, cells = anotatedAstrNK_CellsGSM6204128_obj2) <- "trNK Cells"
Idents(GSM6204128_obj2, cells = anotatedAsTreg_CellsGSM6204128_obj2) <- "Treg"
Idents(GSM6204128_obj2, cells = anotatedAsMast_CellsGSM6204128_obj2) <- "Mast Cells"
Idents(GSM6204128_obj2, cells = anotatedAsCD4_CellsGSM6204128_obj2) <- "CD4+ T Cells"
Idents(GSM6204128_obj2, cells = anotatedAsCD8_CellsGSM6204128_obj2) <- "CD8+ T Cells"
Idents(GSM6204128_obj2, cells = anotatedAsDuctalCells_CellsGSM6204128_obj2) <- "Ductal Cells"
Idents(GSM6204128_obj2, cells = anotatedAsMac_CellsGSM6204128_obj2) <- "Macrophages"
Idents(GSM6204128_obj2, cells = anotatedAsNK_CellsGSM6204128_obj2) <- "NK Cells"
Idents(GSM6204128_obj2, cells = anotatedAsB_CellsGSM6204128_obj2) <- "B Cells"
Idents(GSM6204128_obj2, cells = anotatedAsAcinar_CellsGSM6204128_obj2) <- "Acinar Cells"
Idents(GSM6204128_obj2, cells = anotatedAsFibroblastGSM6204128_obj2) <- "Fibroblasts"
Idents(GSM6204128_obj2, cells = anotatedAsEndothelial_CellsGSM6204128_obj2) <- "Endothelial"
Idents(GSM6204128_obj2, cells = anotatedAsPlasma_CellsGSM6204128_obj2) <- "Plasma Cells"
Idents(GSM6204128_obj2, cells = anotatedAsGastricEpithelial_CellsGSM6204128_obj2) <- "Gastric Epithelial"



#GSM6204129

anotation_for_GSM6204129_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204129.csv", sep = "," )
dim(anotation_for_GSM6204129_obj2)
anotation_for_GSM6204129_obj2[1:6,]
colnames(anotation_for_GSM6204129_obj2)=anotation_for_GSM6204129_obj2[1,]
anotation_for_GSM6204129_obj2=anotation_for_GSM6204129_obj2[-1,]
unique(anotation_for_GSM6204129_obj2$ident)
length(anotatedAsFibroblastGSM6204129_obj2)
anotatedAsTreg_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "Treg"),2]
anotatedAsEndothelial_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "Endothelial"),2]
anotatedAsHepatocytes_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "Hepatocytes"),2]
anotatedAsDuctalCells_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "Ductal Cells"),2]
anotatedAsNeutrophil_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "Neutrophil Cells"),2]
anotatedAsMac_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "Macrophages"),2]
anotatedAsCD4_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "CD4+ T Cells"),2]
anotatedAstrNK_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "trNK"),2]
anotatedAsCD8_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "CD8+ T Cells"),2]
anotatedAsNK_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "NK Cells"),2]
anotatedAsGastricEpithelial_CellsGSM6204129_obj2=anotation_for_GSM6204129_obj2[which(anotation_for_GSM6204129_obj2$ident == "Gastric Epithelial"),2]





length(anotatedAsPlasma_CellsGSM6204129_obj2)
# Assigning Obj2
Idents(GSM6204129_obj2, cells = anotatedAsTreg_CellsGSM6204129_obj2) <- "Treg"
Idents(GSM6204129_obj2, cells = anotatedAsEndothelial_CellsGSM6204129_obj2) <- "Endothelial"
Idents(GSM6204129_obj2, cells = anotatedAsHepatocytes_CellsGSM6204129_obj2) <- "Hepatocytes"
Idents(GSM6204129_obj2, cells = anotatedAsDuctalCells_CellsGSM6204129_obj2) <- "Ductal Cells"
Idents(GSM6204129_obj2, cells = anotatedAsNeutrophil_CellsGSM6204129_obj2) <- "Neutrophil Cells"
Idents(GSM6204129_obj2, cells = anotatedAsMac_CellsGSM6204129_obj2) <- "Macrophages"
Idents(GSM6204129_obj2, cells = anotatedAsCD4_CellsGSM6204129_obj2) <- "CD4+ T Cells"
Idents(GSM6204129_obj2, cells = anotatedAstrNK_CellsGSM6204129_obj2) <- "trNK Cells"
Idents(GSM6204129_obj2, cells = anotatedAsCD8_CellsGSM6204129_obj2) <- "CD8+ T Cells"
Idents(GSM6204129_obj2, cells = anotatedAsNK_CellsGSM6204129_obj2) <- "NK Cells"
Idents(GSM6204129_obj2, cells = anotatedAsGastricEpithelial_CellsGSM6204129_obj2) <- "Gastric Epithelial"




#GSM6204130

anotation_for_GSM6204130_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204130.csv", sep = "," )
dim(anotation_for_GSM6204130_obj2)
anotation_for_GSM6204130_obj2[1:6,]
colnames(anotation_for_GSM6204130_obj2)=anotation_for_GSM6204130_obj2[1,]
anotation_for_GSM6204130_obj2=anotation_for_GSM6204130_obj2[-1,]
unique(anotation_for_GSM6204130_obj2$ident)
length(anotatedAsFibroblastGSM6204130_obj2)
anotatedAstrNK_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "trNK"),2]
anotatedAsPlasma_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "Plasma Cells"),2]
anotatedAsTreg_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "Treg"),2]
anotatedAsCD4_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "CD4+ T Cells"),2]
anotatedAsNK_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "NK Cells"),2]
anotatedAsDuctalCells_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "Ductal Cells"),2]
anotatedAsMac_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "Macrophages"),2]
anotatedAsCD8_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "CD8+ T Cells"),2]
anotatedAsDC_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "DC"),2]
anotatedAsB_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "B Cells"),2]
anotatedAsNeutrophil_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "Neutrophils"),2]
anotatedAsMast_CellsGSM6204130_obj2=anotation_for_GSM6204130_obj2[which(anotation_for_GSM6204130_obj2$ident == "Mast Cells"),2]



length(anotatedAsPlasma_CellsGSM6204130_obj2)
# Assigning Obj2
Idents(GSM6204130_obj2, cells = anotatedAstrNK_CellsGSM6204130_obj2) <- "trNK Cells"
Idents(GSM6204130_obj2, cells = anotatedAsPlasma_CellsGSM6204130_obj2) <- "Plasma Cells"
Idents(GSM6204130_obj2, cells = anotatedAsTreg_CellsGSM6204130_obj2) <- "Treg"
Idents(GSM6204130_obj2, cells = anotatedAsCD4_CellsGSM6204130_obj2) <- "CD4+ T Cells"
Idents(GSM6204130_obj2, cells = anotatedAsNK_CellsGSM6204130_obj2) <- "NK Cells"
Idents(GSM6204130_obj2, cells = anotatedAsDuctalCells_CellsGSM6204130_obj2) <- "Ductal Cells"
Idents(GSM6204130_obj2, cells = anotatedAsMac_CellsGSM6204130_obj2) <- "Macrophages"
Idents(GSM6204130_obj2, cells = anotatedAsCD8_CellsGSM6204130_obj2) <- "CD8+ T Cells"
Idents(GSM6204130_obj2, cells = anotatedAsDC_CellsGSM6204130_obj2) <- "DC"
Idents(GSM6204130_obj2, cells = anotatedAsB_CellsGSM6204130_obj2) <- "B Cells"
Idents(GSM6204130_obj2, cells = anotatedAsNeutrophil_CellsGSM6204130_obj2) <- "Neutrophil Cells"
Idents(GSM6204130_obj2, cells = anotatedAsMast_CellsGSM6204130_obj2) <- "Mast Cells"



#GSM6204131

anotation_for_GSM6204131_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204131.csv", sep = "," )
dim(anotation_for_GSM6204131_obj2)
anotation_for_GSM6204131_obj2[1:6,]
colnames(anotation_for_GSM6204131_obj2)=anotation_for_GSM6204131_obj2[1,]
anotation_for_GSM6204131_obj2=anotation_for_GSM6204131_obj2[-1,]
unique(anotation_for_GSM6204131_obj2$ident)
length(anotatedAsFibroblastGSM6204131_obj2)
anotatedAsCD8_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "CD8+ T Cells"),2]
anotatedAsTreg_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Treg"),2]
anotatedAsNK_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "NK Cells"),2]
anotatedAstrNK_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "trNK"),2]
anotatedAsDuctalCells_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Ductal Cells"),2]
anotatedAsMac_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Macrophages"),2]
anotatedAsCD4_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "CD4+ T Cells"),2]
anotatedAsEndothelial_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Endothelial"),2]
anotatedAsB_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "B Cells"),2]
anotatedAsFibroblasts_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Fibroblasts"),2]
anotatedAsPlasma_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Plasma Cells"),2]
anotatedAsMuscle_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Muscle Cells"),2]
anotatedAsEndocrine_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Endocrine Cells"),2]
anotatedAsMast_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Mast Cells"),2]
anotatedAsAcinar_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "Acinar Cells"),2]
anotatedAsDC_CellsGSM6204131_obj2=anotation_for_GSM6204131_obj2[which(anotation_for_GSM6204131_obj2$ident == "DC"),2]



length(anotatedAsPlasma_CellsGSM6204131_obj2)
# Assigning Obj2
Idents(GSM6204131_obj2, cells = anotatedAsCD8_CellsGSM6204131_obj2) <- "CD8+ T Cells"
Idents(GSM6204131_obj2, cells = anotatedAsTreg_CellsGSM6204131_obj2) <- "Treg"
Idents(GSM6204131_obj2, cells = anotatedAsNK_CellsGSM6204131_obj2) <- "NK Cells"
Idents(GSM6204131_obj2, cells = anotatedAstrNK_CellsGSM6204131_obj2) <- "trNK Cells"
Idents(GSM6204131_obj2, cells = anotatedAsDuctalCells_CellsGSM6204131_obj2) <- "Ductal Cells"
Idents(GSM6204131_obj2, cells = anotatedAsMac_CellsGSM6204131_obj2) <- "Macrophages"
Idents(GSM6204131_obj2, cells = anotatedAsCD4_CellsGSM6204131_obj2) <- "CD4+ T Cells"
Idents(GSM6204131_obj2, cells = anotatedAsEndothelial_CellsGSM6204131_obj2) <- "Endothelial"
Idents(GSM6204131_obj2, cells = anotatedAsB_CellsGSM6204131_obj2) <- "B Cells"
Idents(GSM6204131_obj2, cells = anotatedAsFibroblasts_CellsGSM6204131_obj2) <- "Fibroblasts"
Idents(GSM6204131_obj2, cells = anotatedAsPlasma_CellsGSM6204131_obj2) <- "Plasma Cells"
Idents(GSM6204131_obj2, cells = anotatedAsMuscle_CellsGSM6204131_obj2) <- "Muscle Cells"
Idents(GSM6204131_obj2, cells = anotatedAsEndocrine_CellsGSM6204131_obj2) <- "Endocrine Cells"
Idents(GSM6204131_obj2, cells = anotatedAsMast_CellsGSM6204131_obj2) <- "Mast Cells"
Idents(GSM6204131_obj2, cells = anotatedAsAcinar_CellsGSM6204131_obj2) <- "Acinar Cells"
Idents(GSM6204131_obj2, cells = anotatedAsDC_CellsGSM6204131_obj2) <- "DC"


#GSM6204132

anotation_for_GSM6204132_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204132.csv", sep = "," )
dim(anotation_for_GSM6204132_obj2)
anotation_for_GSM6204132_obj2[1:6,]
colnames(anotation_for_GSM6204132_obj2)=anotation_for_GSM6204132_obj2[1,]
anotation_for_GSM6204132_obj2=anotation_for_GSM6204132_obj2[-1,]
unique(anotation_for_GSM6204132_obj2$ident)
length(anotatedAsFibroblastGSM6204132_obj2)
anotatedAsEndothelial_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Endothelial"),2]
anotatedAsMast_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Mast Cells"),2]
anotatedAstrNK_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "trNK"),2]
anotatedAsDC_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "DC"),2]
anotatedAsCD4_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "CD4+ T Cells"),2]
anotatedAsCD8_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "CD8+ T Cells"),2]
anotatedAsDuctalCells_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Ductal Cells"),2]
anotatedAsNK_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "NK Cells"),2]
anotatedAsHepatocytes_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Hepatocytes"),2]
anotatedAsFibroblasts_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Fibroblasts"),2]
anotatedAsMac_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Macrophages"),2]
anotatedAsNeutrophil_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Neutrophil Cells"),2]
anotatedAsTreg_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Treg"),2]
anotatedAsPlasma_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "Plasma Cells"),2]
anotatedAsB_CellsGSM6204132_obj2=anotation_for_GSM6204132_obj2[which(anotation_for_GSM6204132_obj2$ident == "B Cells"),2]




length(anotatedAsPlasma_CellsGSM6204132_obj2)
# Assigning Obj2
Idents(GSM6204132_obj2, cells = anotatedAsEndothelial_CellsGSM6204132_obj2) <- "Endothelial"
Idents(GSM6204132_obj2, cells = anotatedAsMast_CellsGSM6204132_obj2) <- "Mast Cells"
Idents(GSM6204132_obj2, cells = anotatedAstrNK_CellsGSM6204132_obj2) <- "trNK Cells"
Idents(GSM6204132_obj2, cells = anotatedAsDC_CellsGSM6204132_obj2) <- "DC"
Idents(GSM6204132_obj2, cells = anotatedAsCD4_CellsGSM6204132_obj2) <- "CD4+ T Cells"
Idents(GSM6204132_obj2, cells = anotatedAsCD8_CellsGSM6204132_obj2) <- "CD8+ T Cells"
Idents(GSM6204132_obj2, cells = anotatedAsDuctalCells_CellsGSM6204132_obj2) <- "Ductal Cells"
Idents(GSM6204132_obj2, cells = anotatedAsNK_CellsGSM6204132_obj2) <- "NK Cells"
Idents(GSM6204132_obj2, cells = anotatedAsHepatocytes_CellsGSM6204132_obj2) <- "Hepatocytes"
Idents(GSM6204132_obj2, cells = anotatedAsFibroblasts_CellsGSM6204132_obj2) <- "Fibroblasts"
Idents(GSM6204132_obj2, cells = anotatedAsMac_CellsGSM6204132_obj2) <- "Macrophages"
Idents(GSM6204132_obj2, cells = anotatedAsNeutrophil_CellsGSM6204132_obj2) <- "Neutrophil Cells"
Idents(GSM6204132_obj2, cells = anotatedAsTreg_CellsGSM6204132_obj2) <- "Treg"
Idents(GSM6204132_obj2, cells = anotatedAsPlasma_CellsGSM6204132_obj2) <- "Plasma Cells"
Idents(GSM6204132_obj2, cells = anotatedAsB_CellsGSM6204132_obj2) <- "B Cells"



#GSM6204133

anotation_for_GSM6204133_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204133.csv", sep = "," )
dim(anotation_for_GSM6204133_obj2)
anotation_for_GSM6204133_obj2[1:6,]
colnames(anotation_for_GSM6204133_obj2)=anotation_for_GSM6204133_obj2[1,]
anotation_for_GSM6204133_obj2=anotation_for_GSM6204133_obj2[-1,]
unique(anotation_for_GSM6204133_obj2$ident)
length(anotatedAsFibroblastGSM6204133_obj2)
anotatedAstrNK_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "trNK"),2]
anotatedAsFibroblasts_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "Fibroblasts"),2]
anotatedAsMast_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "Mast Cells"),2]
anotatedAsHepatocytes_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "Hepatocytes"),2]
anotatedAsCD8_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "CD8+ T Cells"),2]
anotatedAsCD4_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "CD4+ T Cells"),2]
anotatedAsDuctalCells_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "Ductal Cells"),2]
anotatedAsMac_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "Macrophages"),2]
anotatedAsNK_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "NK Cells"),2]
anotatedAsB_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "B Cells"),2]
anotatedAsPlasma_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "Plasma Cells"),2]
anotatedAsEndothelial_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "Endothelial"),2]
anotatedAsTreg_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "Treg"),2]
anotatedAsDC_CellsGSM6204133_obj2=anotation_for_GSM6204133_obj2[which(anotation_for_GSM6204133_obj2$ident == "DC"),2]





length(anotatedAsPlasma_CellsGSM6204133_obj2)
# Assigning Obj2
Idents(GSM6204133_obj2, cells = anotatedAstrNK_CellsGSM6204133_obj2) <- "trNK Cells"
Idents(GSM6204133_obj2, cells = anotatedAsFibroblasts_CellsGSM6204133_obj2) <- "Fibroblasts"
Idents(GSM6204133_obj2, cells = anotatedAsMast_CellsGSM6204133_obj2) <- "Mast Cells"
Idents(GSM6204133_obj2, cells = anotatedAsHepatocytes_CellsGSM6204133_obj2) <- "Hepatocytes"
Idents(GSM6204133_obj2, cells = anotatedAsCD8_CellsGSM6204133_obj2) <- "CD8+ T Cells"
Idents(GSM6204133_obj2, cells = anotatedAsCD4_CellsGSM6204133_obj2) <- "CD4+ T Cells"
Idents(GSM6204133_obj2, cells = anotatedAsDuctalCells_CellsGSM6204133_obj2) <- "Ductal Cells"
Idents(GSM6204133_obj2, cells = anotatedAsMac_CellsGSM6204133_obj2) <- "Macrophages"
Idents(GSM6204133_obj2, cells = anotatedAsNK_CellsGSM6204133_obj2) <- "NK Cells"
Idents(GSM6204133_obj2, cells = anotatedAsB_CellsGSM6204133_obj2) <- "B Cells"
Idents(GSM6204133_obj2, cells = anotatedAsPlasma_CellsGSM6204133_obj2) <- "Plasma Cells"
Idents(GSM6204133_obj2, cells = anotatedAsEndothelial_CellsGSM6204133_obj2) <- "Endothelial"
Idents(GSM6204133_obj2, cells = anotatedAsTreg_CellsGSM6204133_obj2) <- "Treg"
Idents(GSM6204133_obj2, cells = anotatedAsDC_CellsGSM6204133_obj2) <- "DC"



#GSM6204134

anotation_for_GSM6204134_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204134.csv", sep = "," )
dim(anotation_for_GSM6204134_obj2)
anotation_for_GSM6204134_obj2[1:6,]
colnames(anotation_for_GSM6204134_obj2)=anotation_for_GSM6204134_obj2[1,]
anotation_for_GSM6204134_obj2=anotation_for_GSM6204134_obj2[-1,]
unique(anotation_for_GSM6204134_obj2$ident)
length(anotatedAsFibroblastGSM6204134_obj2)
anotatedAsDC_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "DC"),2]
anotatedAsEndothelial_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "Endothelial"),2]
anotatedAstrNK_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "trNK"),2]
anotatedAsPlasma_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "Plasma Cells"),2]
anotatedAsB_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "B Cells"),2]
anotatedAsDuctalCells_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "Ductal Cells"),2]
anotatedAsCD4_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "CD4+ T Cells"),2]
anotatedAsMac_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "Macrophages"),2]
anotatedAsNK_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "NK Cells"),2]
anotatedAsFibroblasts_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "Fibroblasts"),2]
anotatedAsCD8_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "CD8+ T Cells"),2]
anotatedAsTreg_CellsGSM6204134_obj2=anotation_for_GSM6204134_obj2[which(anotation_for_GSM6204134_obj2$ident == "Treg"),2]





length(anotatedAsPlasma_CellsGSM6204134_obj2)
# Assigning Obj2
Idents(GSM6204134_obj2, cells = anotatedAsDC_CellsGSM6204134_obj2) <- "DC"
Idents(GSM6204134_obj2, cells = anotatedAsEndothelial_CellsGSM6204134_obj2) <- "Endothelial"
Idents(GSM6204134_obj2, cells = anotatedAstrNK_CellsGSM6204134_obj2) <- "trNK Cells"
Idents(GSM6204134_obj2, cells = anotatedAsPlasma_CellsGSM6204134_obj2) <- "Plasma Cells"
Idents(GSM6204134_obj2, cells = anotatedAsB_CellsGSM6204134_obj2) <- "B Cells"
Idents(GSM6204134_obj2, cells = anotatedAsDuctalCells_CellsGSM6204134_obj2) <- "Ductal Cells"
Idents(GSM6204134_obj2, cells = anotatedAsCD4_CellsGSM6204134_obj2) <- "CD4+ T Cells"
Idents(GSM6204134_obj2, cells = anotatedAsMac_CellsGSM6204134_obj2) <- "Macrophages"
Idents(GSM6204134_obj2, cells = anotatedAsNK_CellsGSM6204134_obj2) <- "NK Cells"
Idents(GSM6204134_obj2, cells = anotatedAsFibroblasts_CellsGSM6204134_obj2) <- "Fibroblasts"
Idents(GSM6204134_obj2, cells = anotatedAsCD8_CellsGSM6204134_obj2) <- "CD8+ T Cells"
Idents(GSM6204134_obj2, cells = anotatedAsTreg_CellsGSM6204134_obj2) <- "Treg"



#GSM6204135

anotation_for_GSM6204135_obj2=read.table("/Users/hossein.allahdadi/Downloads/annotations/CellAnnotation_GSM6204135.csv", sep = "," )
dim(anotation_for_GSM6204135_obj2)
anotation_for_GSM6204135_obj2[1:6,]
colnames(anotation_for_GSM6204135_obj2)=anotation_for_GSM6204135_obj2[1,]
anotation_for_GSM6204135_obj2=anotation_for_GSM6204135_obj2[-1,]
unique(anotation_for_GSM6204135_obj2$ident)
length(anotatedAsFibroblastGSM6204135_obj2)
anotatedAsDC_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "DC"),2]
anotatedAsTreg_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "Treg"),2]
anotatedAstrNK_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "trNK"),2]
anotatedAsNK_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "NK Cells"),2]
anotatedAsDuctalCells_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "Ductal Cells"),2]
anotatedAsCD4_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "CD4+ T Cells"),2]
anotatedAsMac_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "Macrophages"),2]
anotatedAsCD8_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "CD8+ T Cells"),2]
anotatedAsPlasma_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "Plasma Cells"),2]
anotatedAsEndothelial_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "Endothelial"),2]
anotatedAsFibroblasts_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "Fibroblasts"),2]
anotatedAsB_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "B Cells"),2]
anotatedAsHepatocytes_CellsGSM6204135_obj2=anotation_for_GSM6204135_obj2[which(anotation_for_GSM6204135_obj2$ident == "Hepatocytes"),2]





length(anotatedAsPlasma_CellsGSM6204135_obj2)
# Assigning Obj2
Idents(GSM6204135_obj2, cells = anotatedAsDC_CellsGSM6204135_obj2) <- "DC"
Idents(GSM6204135_obj2, cells = anotatedAsTreg_CellsGSM6204135_obj2) <- "Treg"
Idents(GSM6204135_obj2, cells = anotatedAstrNK_CellsGSM6204135_obj2) <- "trNK Cells"
Idents(GSM6204135_obj2, cells = anotatedAsNK_CellsGSM6204135_obj2) <- "NK Cells"
Idents(GSM6204135_obj2, cells = anotatedAsDuctalCells_CellsGSM6204135_obj2) <- "Ductal Cells"
Idents(GSM6204135_obj2, cells = anotatedAsCD4_CellsGSM6204135_obj2) <- "CD4+ T Cells"
Idents(GSM6204135_obj2, cells = anotatedAsMac_CellsGSM6204135_obj2) <- "Macrophages"
Idents(GSM6204135_obj2, cells = anotatedAsCD8_CellsGSM6204135_obj2) <- "CD8+ T Cells"
Idents(GSM6204135_obj2, cells = anotatedAsPlasma_CellsGSM6204135_obj2) <- "Plasma Cells"
Idents(GSM6204135_obj2, cells = anotatedAsEndothelial_CellsGSM6204135_obj2) <- "Endothelial"
Idents(GSM6204135_obj2, cells = anotatedAsFibroblasts_CellsGSM6204135_obj2) <- "Fibroblasts"
Idents(GSM6204135_obj2, cells = anotatedAsB_CellsGSM6204135_obj2) <- "B Cells"
Idents(GSM6204135_obj2, cells = anotatedAsHepatocytes_CellsGSM6204135_obj2) <- "Hepatocytes"