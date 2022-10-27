---
Title: "Dental pulp data integration"
Author: Natnicha Jiravejchakul
Date: '2022-07-11'
---


library(Seurat)
library(SeuratObject)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(viridis)


###


DTP_01[["RNA"]]@counts
pulp1[["RNA"]]@counts
pulp2[["RNA"]]@counts
pulp3[["RNA"]]@counts
pulp4[["RNA"]]@counts
pulp5[["RNA"]]@counts
DTP_MH1[["RNA"]]@counts
DTP_MH2[["RNA"]]@counts
DTP_MH3[["RNA"]]@counts

DTP_01@meta.data$dataset <- "DTP1"
pulp1@meta.data$dataset <- "DTP2"
pulp2@meta.data$dataset <- "DTP3"
pulp3@meta.data$dataset <- "DTP4"
pulp4@meta.data$dataset <- "DTP5"
pulp5@meta.data$dataset <- "DTP6"
DTP_MH1@meta.data$dataset <- "DTP7"
DTP_MH2@meta.data$dataset <- "DTP8"
DTP_MH3@meta.data$dataset <- "DTP9"

DTP_01@meta.data$tissue <- "DTP"
pulp1@meta.data$tissue <- "DTP"
pulp2@meta.data$tissue <- "DTP"
pulp3@meta.data$tissue <- "DTP"
pulp4@meta.data$tissue <- "DTP"
pulp5@meta.data$tissue <- "DTP"
DTP_MH1@meta.data$tissue <- "DTP"
DTP_MH2@meta.data$tissue <- "DTP"
DTP_MH3@meta.data$tissue <- "DTP"



================================================= SEURAT INTEGRATION ==============================================
  

library(Seurat)
library(SeuratData)
library(patchwork)

Integ_list_pulp <- list(DTP_01, pulp1, pulp2, pulp3, pulp4, pulp5, DTP_MH1, DTP_MH2, DTP_MH3)

Pulp_data_integ_seurat <- lapply(X = Integ_list_pulp, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Pulp_data_integ_seurat)
Integ.anchors <- FindIntegrationAnchors(object.list = Pulp_data_integ_seurat, anchor.features = features)
Pulp_data_integ_seurat <- IntegrateData(anchorset = Integ.anchors)

DefaultAssay(Pulp_data_integ_seurat) <- "integrated"

Pulp_data_integ_seurat <- ScaleData(Pulp_data_integ_seurat, verbose = FALSE)
Pulp_data_integ_seurat <- RunPCA(Pulp_data_integ_seurat, npcs = 30, verbose = FALSE)
Pulp_data_integ_seurat <- RunUMAP(Pulp_data_integ_seurat, reduction = "pca", dims = 1:30)
Pulp_data_integ_seurat <- FindNeighbors(Pulp_data_integ_seurat, reduction = "pca", dims = 1:30)
Pulp_data_integ_seurat <- FindClusters(Pulp_data_integ_seurat, resolution = 0.6)
Pulp_data_integ_seurat

# Visualization - FIGURE 3A
DimPlot(Pulp_data_integ_seurat, reduction = "umap", label = TRUE, repel = TRUE, cols = c("#749B58B2", "#F0E685B2", "#466983B2", "#BA6338B2", "#5DB1DDB2", "#802268B2", "#D595A7B2",
"#837B8DB2", "#C75127B2", "#D58F5CB2", "#7A65A5B2", "#E4AF69B2", "#3B1B53B2", "#CDDEB7B2", 
"#612A79B2", "#AE1F63B2", "#E7C76FB2", "#5A655EB2", "#CC9900B2", "#99CC00B2", "#A9A9A9B2")) # igv palette

Pulp_data_integ_seurat # 30972 features across 48659 samples
levels(Pulp_data_integ_seurat)

# Save Data
saveRDS(Pulp_data_integ_seurat, file = "~/path_to_directory/object_name.rds")



================================================= CLUSTER ANNOTATION ================================================


# Expression of marker genes by clusters - FIGURE 3B

Pulp_integ_newCR <- readRDS(file = "~/path_to_directory")

Pulp_integ_newCR
DimPlot(Pulp_integ_newCR)

DotPlot(Pulp_integ_newCR, assay = "RNA", 
        features = c("DMP1",
                     "COL1A1", "COL1A2", "COL3A1", "COL21A1", "DCN",
                     "SFN", "FDCSP", "ODAM",                            
                     "KRT5", "KRT14",
                     "ACTA2", "FRZB", "NOTCH3", "MYH11",                
                     "CSPG4", "PDGFRB",                                           
                     "TAGLN", "TPM2",
                     "SOX10", "GFRA3", "MBP", 
                     "SELE", "PECAM1", "EMCN",  
                     "PTPRC",                                          
                     "CD3D", "CD3E", "NKG7", "GZMA", "PRF1",
                     "C1QA", "LYZ", "LST1", "MS4A7", "CD14",
                     "CD79A", "JCHAIN", "IGKC",                 
                     "HBA1", "HBB"
                     ),                                    
        cols = c("lightgrey", "blue"), 
        cluster.idents = FALSE) + #RotatedAxis() + 
        coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option = "cividis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))



# Rename clusters
cluster_ids <- c("Fibroblast_1", "Fibroblast_2", "Non-myelinated Schw", "Endothelial_1", "MSCs_1", "Endothelial_2",
                 "NK/T cells", "Fibroblast_3", "Fibroblast_4", "MSCs_2", "Myelinated Schw",
                 "Monocytes/Macrophages", "Endothelial_3", "Fibroblast_5", "Fibroblast_6", "Odontoblast",  
                 "Epithelial", "Erythrocytes", "B cells")

names(cluster_ids) <- levels(Pulp_integ_newCR)
Pulp_integ_newCR_id <- RenameIdents(Pulp_integ_newCR, cluster_ids)
DimPlot(Pulp_integ_newCR_id, reduction = "umap", label= TRUE, repel = TRUE, label.size = 3,  cols = c("#749B58B2", "#466983B2", 
        "#BA6338B2", "#5DB1DDB2", "#802268B2", "#D595A7B2", "#837B8DB2", "#C75127B2", "#D58F5CB2", "#7A65A5B2", "#E4AF69B2", 
        "#3B1B53B2", "#CDDEB7B2", "#612A79B2", "#AE1F63B2", "#E7C76FB2", "#5A655EB2", "#CC9900B2", "#99CC00B2", "#A9A9A9B2"))
Pulp_integ_newCR_id

saveRDS(Pulp_data_integ_seurat, file = "~/path_to_directory/object_name.rds")



# Retrieve cell number/percentage from the UMAP - FIGURE 3C


library(ggplot2)
library(tidyr)
library(tibble)

# Number of cells by datasets, by clusters  
df <- data.frame(rbind(table(Idents(Pulp_integ_newCR_id), Pulp_integ_newCR_id$dataset)))
df

# calulate percentage 
for (RN in 1:length(df)) {df[,RN] <- df[,RN] / sum(df[,RN])}
rm(RN)
df_p <- round(x = df * 100 , digits = 2)
df_p



# Cell number
data <- data.frame(
  class = c("Fibroblast_1", "Fibroblast_2", "Non-myelinated Schw", "Endothelial_1", "MSCs_1", "Endothelial_2", "NK/T cells", "Fibroblast_3", "Fibroblast_4", "MSCs_2", "Myelinated Schw", "Monocytes/Macrophages", "Endothelial_3", "Fibroblast_5", "Fibroblast_6", "Odontoblast", "Epithelial", "Fibroblast_7", "B cells"),
  cell_number = c(6954,  6291,  3908, 3813, 3706, 3549, 3403, 3336, 3330, 2888, 2018, 1379, 1271, 541,  541,  292,  221,  196, 185))

ggplot() + geom_bar(aes(y = cell_number, x= "", fill = class), stat="identity" , data = data ) + scale_fill_manual(values=c("#99CC00B2", "#5DB1DDB2", "#D595A7B2",  
"#CDDEB7B2", "#5A655EB2", "#749B58B2", "#466983B2", "#C75127B2", "#D58F5CB2", "#612A79B2", "#AE1F63B2", "#CC9900B2", "#3B1B53B2", "#802268B2", "#7A65A5B2", "#BA6338B2",   "#837B8DB2", "#E4AF69B2", "#CC9900B2")) + theme(axis.text = element_text(size = 12 , face="bold"), 
        axis.text.x = element_text(size=12, angle=45 , hjust = 1), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(color="black", size=15, face="bold"), 
        strip.background = element_rect(colour = "white", fill = FALSE), 
        legend.text = element_text(size = 12) ) + 
        ylab("Percentage") + theme_minimal()


# Percentage - FIGURE 3C
data <- data.frame(
  class = c("Fibroblast_1", "Fibroblast_2", "Non-myelinated Schw", "Endothelial_1", "MSCs_1", "Endothelial_2", "NK/T cells", "Fibroblast_3", "Fibroblast_4", "MSCs_2", "Myelinated Schw", "Monocytes/Macrophages", "Endothelial_3", "Fibroblast_5", "Fibroblast_6", "Odontoblast", "Epithelial", "Fibroblast_7", "B cells"),
  cell_percentage = c(14.54, 13.16, 8.17, 7.97, 7.75, 7.42, 7.1,  7,  6.96, 6.04, 4.22, 2.88, 2.66, 1.13, 1.13, 0.61, 0.46, 0.4, 0.4))

ggplot() + geom_bar(aes(y = cell_percentage, x= "", fill = class), stat="identity" , data = data ) + scale_fill_manual(values=c("#99CC00B2", "#5DB1DDB2", "#D595A7B2",  
"#CDDEB7B2", "#5A655EB2", "#749B58B2", "#466983B2", "#C75127B2", "#D58F5CB2", "#612A79B2", "#AE1F63B2", "#CC9900B2", "#3B1B53B2", "#802268B2", "#7A65A5B2", "#BA6338B2",   "#837B8DB2", "#E4AF69B2", "#CC9900B2")) + theme(axis.text = element_text(size = 12 , face="bold"), 
        axis.text.x = element_text(size=12, angle=45 , hjust = 1), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(color="black", size=15, face="bold"), 
        strip.background = element_rect(colour = "white", fill = FALSE), 
        legend.text = element_text(size = 12) ) + 
        ylab("Percentage") + theme_minimal()

