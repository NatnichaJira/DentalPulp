---
Title: "Data Integration - 19 datasets, 6 tissues"
Author: Natnicha Jiravejchakul
Date: '2022-07-11'
---

ALL DATASETS USED IN THIS PROJECT WERE DOWLOADED AS FASTQ FILES 
PROCESSED WITH CELL RANGER 7.0.0 - REFERENCE = HUMAN GRCH38-2020-A


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
New_pbmc1k[["RNA"]]@counts
New_pbmc8k[["RNA"]]@counts
BM_CD34_1[["RNA"]]@counts
BM_CD34_2[["RNA"]]@counts
ADP_L1[["RNA"]]@counts
ADP_L2[["RNA"]]@counts
Lung_CR7_1[["RNA"]]@counts
Lung_CR7_2[["RNA"]]@counts
Skin_1[["RNA"]]@counts
Skin_2[["RNA"]]@counts


DTP_01@meta.data$dataset <- "DTP1"
pulp1@meta.data$dataset <- "DTP2"
pulp2@meta.data$dataset <- "DTP3"
pulp3@meta.data$dataset <- "DTP4"
pulp4@meta.data$dataset <- "DTP5"
pulp5@meta.data$dataset <- "DTP6"
DTP_MH1@meta.data$dataset <- "DTP7"
DTP_MH2@meta.data$dataset <- "DTP8"
DTP_MH3@meta.data$dataset <- "DTP9"
New_pbmc1k@meta.data$dataset <- "PBMC1"
New_pbmc8k@meta.data$dataset <- "PBMC2"
BM_CD34_1@meta.data$dataset <- "BM1"
BM_CD34_2@meta.data$dataset <- "BM2"
ADP_L1@meta.data$dataset <- "ADP1"
ADP_L2@meta.data$dataset <- "ADP2"
Lung_CR7_1@meta.data$dataset <- "LUNG1"
Lung_CR7_2@meta.data$dataset <- "LUNG2"
Skin_1@meta.data$dataset <- "SKIN1"
Skin_2@meta.data$dataset <- "SKIN2"


DTP_01@meta.data$tissue <- "DTP"
pulp1@meta.data$tissue <- "DTP"
pulp2@meta.data$tissue <- "DTP"
pulp3@meta.data$tissue <- "DTP"
pulp4@meta.data$tissue <- "DTP"
pulp5@meta.data$tissue <- "DTP"
DTP_MH1@meta.data$tissue <- "DTP"
DTP_MH2@meta.data$tissue <- "DTP"
DTP_MH3@meta.data$tissue <- "DTP"
New_pbmc1k@meta.data$tissue <- "PBMC"
New_pbmc8k@meta.data$tissue <- "PBMC"
BM_CD34_1@meta.data$tissue <- "BM"
BM_CD34_2@meta.data$tissue <- "BM"
ADP_L1@meta.data$tissue <- "ADP"
ADP_L2@meta.data$tissue <- "ADP"
Lung_CR7_1@meta.data$tissue <- "LUNG"
Lung_CR7_2@meta.data$tissue <- "LUNG"
Skin_1@meta.data$tissue <- "SKIN"
Skin_2@meta.data$tissue <- "SKIN"


-------------------------------------- SEURAT INTEGRATION ---------------------------------------

# FIGURE 1A - 1B --- INTEGRATION OF 19 DATASETS FROM 6 TISSUES


library(Seurat)
library(patchwork)
library(ggplot2)

# Create a list of seurat objects
Integ_list_new_anno <- list(DTP_01, pulp1, pulp2, pulp3, pulp4, pulp5, DTP_MH1, DTP_MH2, DTP_MH3, 
                            New_pbmc1k, New_pbmc8k, 
                            BM_CD34_1, BM_CD34_2,
                            ADP_L1, ADP_L2,
                            Lung_CR7_1, Lung_CR7_2,
                            Skin_1, Skin_2)

Integ_list_new_anno

Integ_list_new_anno_seurat <- lapply(X = Integ_list_new_anno, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features_new_anno_seurat <- SelectIntegrationFeatures(object.list = Integ_list_new_anno_seurat)
In_anchors_new_anno_seurat <- FindIntegrationAnchors(object.list = Integ_list_new_anno_seurat, anchor.features = features_new_anno_seurat)

Integ_combined_new_anno_seurat <- IntegrateData(anchorset = In_anchors_new_anno_seurat)
DefaultAssay(Integ_combined_new_anno_seurat) <- "integrated"

Integ_combined_new_anno_seurat <- ScaleData(Integ_combined_new_anno_seurat, verbose = FALSE)
Integ_combined_new_anno_seurat <- RunPCA(Integ_combined_new_anno_seurat, npcs = 30, verbose = FALSE)
Integ_combined_new_anno_seurat <- RunUMAP(Integ_combined_new_anno_seurat, reduction = "pca", dims = 1:30)
Integ_combined_new_anno_seurat <- FindNeighbors(Integ_combined_new_anno_seurat, reduction = "pca", dims = 1:30)
Integ_combined_new_anno_seurat <- FindClusters(Integ_combined_new_anno_seurat, resolution = 0.5)

DimPlot(Integ_combined_new_anno_seurat, reduction = "umap", group.by = "tissue", 
                              cols = c("#5A655EB2","#466983B2", "#CC9900B2" , "#802268B2", "#BA6338B2", "#749B58B2"))

DimPlot(Integ_combined_new_anno_seurat, reduction = "umap", cols = c("#749B58B2", "#466983B2",
        "#BA6338B2", "#5DB1DDB2", "#802268B2", "#D595A7B2", "#837B8DB2", "#C75127B2", "#D58F5CB2", "#7A65A5B2", "#E4AF69B2", 
        "#3B1B53B2", "#CDDEB7B2", "#612A79B2", "#AE1F63B2", "#E7C76FB2", "#5A655EB2", "#CC9900B2", "#99CC00B2", "#A9A9A9B2"))

DimPlot(Integ_combined_new_anno_seurat, reduction = "umap", group.by = "dataset")
DimPlot(Integ_combined_new_anno_seurat, reduction = "umap", split.by = "tissue")

# Save the integrated data
saveRDS(Integ_combined_new_anno_seurat, file = "~/path_to_directory/object_name.rds")



---------------------------------- CLUSTER ANNOTATION ----------------------------------

# Gene markers used for the annotation (from literatures) - Supplementary figure 1

library(ggsci)
library(ggplot2)
library(gridExtra)
library(viridis)

Integ_new_CR7_seurat <- readRDS(file = "~/path_to_directory/object_name.rds")
Integ_new_CR7_seurat 

DotPlot(Integ_new_CR7_seurat, assay = "RNA", 
        features = c("PTPRC",                                           # CD45-Immune cells
                     "CD19", "CD79A", "JCHAIN", "IGKC",                 # B cells
                     "CD3D", "CD3E", "NCR1", "NKG7", "GZMA", "PRF1",    # T and NK
                     "C1QA", "LYZ", "LST1", "MS4A7", "CD14",            # Mono- and Macrophages
                     "DMP1",                                            # Odontoblast 
                     "COL1A1", "COL1A2", "COL3A1", "DCN",               # Fibroblast
                     "SFN", "FDCSP", "ODAM",                            # Epithelial cells
                     "KRT1", "KRT5", "KRT14", "KRT15",                  # Keratinocytes
                     "SOX10", "MBP", "GFRA3",                           # Schwann cells
                     "SELE", "PECAM1", "EMCN",                          # Endothelial
                     "ACTA2", "FRZB", "NOTCH3", "MYH11",                # MSC
                     "CSPG4",                                           # Pericytes 
                     "TAGLN", "TPM2",                                   # SMC
                     "MKI67", "TOP2A",                                  # Proliferating
                     "HBA1", "HBB"                                      # RBC
                    ),                                                             
        cols = c("lightgrey", "blue"), 
        cluster.idents = TRUE) + RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option = "cividis") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


# Rename Clusters

cluster_ids <- c("Fibroblast_1", "NK/T cells_1", "Fibroblast_2", "Fibroblast_3", 
                 "Unknown", "MSC-like_1", "Endothelial_1", "Fibroblast_4", "Macrophages", 
                 "Non-myelinated Schw", "MSC-like_2", "Myelinated Schw", "Endothelial_2",
                 "Monocytes", "NK/T cells_2", "B cells", "Epithelial",
                 "Odontoblast", "Pro_Fibroblast")

names(cluster_ids) <- levels(Integ_new_CR7_seurat)
Integ_new_CR7_seurat_ClusterAnno <- RenameIdents(Integ_new_CR7_seurat, cluster_ids)

DimPlot(Integ_new_CR7_seurat_ClusterAnno, reduction = "umap", repel = TRUE, label = TRUE, label.size = 3) + NoLegend()

saveRDS(Integ_new_CR7_seurat_ClusterAnno, file = "~/path_to_directory/object_name.rds")


# FIGURE 1D - Cell type Fraction

ggplot(Integ_new_CR7_seurat[[]], aes(seurat_clusters, fill=tissue)) + geom_bar(position="fill") +
  scale_fill_manual(values=c("#5A655EB2","#466983B2", "#CC9900B2" , "#802268B2", "#BA6338B2", "#749B58B2"))



========================================================================================================

FIBROBLASTS

# Subset the fibroblasts: clusters 0, 2, 3, and 7 for further downstream analysis -- FIGURE 2E, 2F


Idents(Integ_new_CR7_seurat) <- "seurat_clusters"
levels(Integ_new_CR7_seurat)

Integ_new_CR7_fibroblast <- subset(x = Integ_new_CR7_seurat, idents = c("0", "2", "3", "7"))
Integ_new_CR7_fibroblast

DefaultAssay(Integ_new_CR7_fibroblast) <- "integrated"

Integ_new_CR7_fibroblast <- NormalizeData(Integ_new_CR7_fibroblast, normalization.method = "LogNormalize")
Integ_new_CR7_fibroblast <- FindVariableFeatures(Integ_new_CR7_fibroblast, selection.method = "vst", nfeatures = 2000)
Integ_new_CR7_fibroblast <- ScaleData(Integ_new_CR7_fibroblast)
Integ_new_CR7_fibroblast <- RunPCA(Integ_new_CR7_fibroblast, npcs = 30, features = VariableFeatures(object = Integ_new_CR7_fibroblast))
Integ_new_CR7_fibroblast <- FindNeighbors(Integ_new_CR7_fibroblast, dims = 1:10)
Integ_new_CR7_fibroblast <- FindClusters(Integ_new_CR7_fibroblast, resolution = 0.5)

Integ_new_CR7_fibroblast <- RunUMAP(Integ_new_CR7_fibroblast, dims = 1:10)
DimPlot(Integ_new_CR7_fibroblast, reduction = "umap", label = TRUE)
DimPlot(Integ_new_CR7_fibroblast, reduction = "umap", group.by = "dataset")
DimPlot(Integ_new_CR7_fibroblast, reduction = "umap", group.by = "tissue", 
        cols = c("#5A655EB2","#466983B2", "#CC9900B2" , "#802268B2", "#BA6338B2", "#749B58B2"))

saveRDS(Integ_new_CR7_fibroblast, file = "~/path_to_directory/object_name.rds")


# Cell fraction : FIGURE 2G

ggplot(Integ_new_CR7_fibroblast[[]], aes(seurat_clusters, fill=tissue)) + geom_bar(position="fill") + 
  scale_fill_manual(values=c("#5A655EB2","#466983B2", "#CC9900B2" , "#802268B2", "#BA6338B2", "#749B58B2"))



# Expression of PTN and MDK in DP fibroblast

# FIGURE 2D
Idents(Integ_new_CR7_fibroblast) <- "dataset"
DotPlot(Integ_new_CR7_seurat, features = c("PTN", "MDK"), assay = "RNA") + RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option = "cividis", direction = 1) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

# FIGURE 2H
Idents(Integ_new_CR7_fibroblast) <- "seurat_clusters"
DotPlot(Integ_new_CR7_seurat, features = c("PTN", "MDK"), assay = "RNA") + RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option = "cividis", direction = 1) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))




