---
Title: "CPDB analysis of pulp data"
Author: Natnicha Jiravejchakul
Date: '2022-08-1'
---


library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(viridis)

# Import integrated result 
Seurat_obj_integPulp <- readRDS(file ="~/path_to_directory")
Seurat_obj_integPulp


--------------------------------------------------------------------------------------------------------------------------------
  

1. Prepared input data for CPDB


levels(Seurat_obj_integPulp)
Seurat_obj_integPulp <- NormalizeData(object = Seurat_obj_integPulp)

writeMM(Seurat_obj_integPulp@assays$RNA@data, file = '~/matrix.mtx')
write(x = rownames(Seurat_obj_integPulp@assays$RNA@data), file = "~/features.tsv")
write(x = colnames(Seurat_obj_integPulp@assays$RNA@data), file = "~/barcodes.tsv")

table(Seurat_obj_integPulp@meta.data$cell_type)
Seurat_obj_integPulp@meta.data$Cell = rownames(Seurat_obj_integPulp@meta.data)
df = Seurat_obj_integPulp@meta.data[, c('Cell', 'cell_type')]
write.table(df, file ='~/path_to_directory/object_name.tsv', sep = '\t', quote = F, row.names = F)

## Extract DEGs for each cell type
DefaultAssay(Seurat_obj_integPulp) <- "RNA"
DEGs <- FindAllMarkers(Seurat_obj_integPulp, 
                       test.use = 'LR', 
                       verbose = F, 
                       only.pos = T, 
                       random.seed = 1, 
                       logfc.threshold = 0.2, 
                       min.pct = 0.1, 
                       return.thresh = 0.05)

library(limma)
'DKK1' %in% rownames(Seurat_obj_integPulp@assays$RNA@counts)

fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC > 0.1)

# 1st column = cluster; 2nd column = gene 
fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
write.table(fDEGs, file ='~/path_to_directory/object_name.tsv', sep = '\t', quote = F, row.names = F)

head(fDEGs)

## Data was further processed by CellPhoneDB DEG-based Methods 


--------------------------------------------------------------------------------------------------------------------------------

  
2. Vizualized CPDB result

library(scmisc)
library(here)

CPDB <- read_cpdb_out("~/path_to_directiry_CPDBoutput")
plot_cpdb_dotplot(CPDB)

# Select specific ligands/receptors to plot
PTN_gene <- grep("PTN", CPDB$proteins, value=TRUE)
plot_cpdb_dotplot(CPDB, 
                  cells=c("Non-myelinated Schw|Fibroblast_1", "Non-myelinated Schw|Fibroblast_2", "Non-myelinated Schw|Fibroblast_3",  
                          "Non-myelinated Schw|Fibroblast_4", "Non-myelinated Schw|Fibroblast_5", "Non-myelinated Schw|Fibroblast_6",
                          "Myelinated Schw|Fibroblast_1", "Myelinated Schw|Fibroblast_2", "Myelinated Schw|Fibroblast_3",
                          "Myelinated Schw|Fibroblast_4", "Myelinated Schw|Fibroblast_5", "Myelinated Schw|Fibroblast_6",
                          "Monocytes/Macrophages|Fibroblast_1","Monocytes/Macrophages|Fibroblast_2",
                          "Monocytes/Macrophages|Fibroblast_3","Monocytes/Macrophages|Fibroblast_4",
                          "Monocytes/Macrophages|Fibroblast_5","Monocytes/Macrophages|Fibroblast_6",
                          "Odontoblast|Fibroblast_1","Odontoblast|Fibroblast_2","Odontoblast|Fibroblast_3",
                          "Odontoblast|Fibroblast_4","Odontoblast|Fibroblast_5","Odontoblast|Fibroblast_6",
                          "Endothelial_1|Fibroblast_1",
                          "Endothelial_1|Fibroblast_2",
                          "Endothelial_1|Fibroblast_3",
                          "Endothelial_1|Fibroblast_4",
                          "Endothelial_1|Fibroblast_5",
                          "Endothelial_1|Fibroblast_6",
                          "Endothelial_2|Fibroblast_1",
                          "Endothelial_2|Fibroblast_2",
                          "Endothelial_2|Fibroblast_3",
                          "Endothelial_2|Fibroblast_4",
                          "Endothelial_2|Fibroblast_5",
                          "Endothelial_2|Fibroblast_6",
                          "Endothelial_3|Fibroblast_1",
                          "Endothelial_3|Fibroblast_2",
                          "Endothelial_3|Fibroblast_3",
                          "Endothelial_3|Fibroblast_4",
                          "Endothelial_3|Fibroblast_5",
                          "Endothelial_3|Fibroblast_6",
                          "Fibroblast_1|MSCs_1","Fibroblast_1|MSCs_2",
                          "Fibroblast_2|MSCs_1","Fibroblast_2|MSCs_2",
                          "Fibroblast_3|MSCs_1","Fibroblast_3|MSCs_2",
                          "Fibroblast_4|MSCs_1","Fibroblast_4|MSCs_2",
                          "Fibroblast_5|MSCs_1","Fibroblast_5|MSCs_2",
                          "Fibroblast_6|MSCs_1","Fibroblast_6|MSCs_2",
                          "Fibroblast_1|Odontoblast",
                          "Fibroblast_2|Odontoblast",
                          "Fibroblast_3|Odontoblast",
                          "Fibroblast_4|Odontoblast",
                          "Fibroblast_5|Odontoblast",
                          "Fibroblast_6|Odontoblast",
                          "Fibroblast_1|Monocytes/Macrophages",
                          "Fibroblast_2|Monocytes/Macrophages",
                          "Fibroblast_3|Monocytes/Macrophages",
                          "Fibroblast_4|Monocytes/Macrophages",
                          "Fibroblast_5|Monocytes/Macrophages",
                          "Fibroblast_6|Monocytes/Macrophages",
                          "Fibroblast_1|Non-myelinated Schw",
                          "Fibroblast_2|Non-myelinated Schw",
                          "Fibroblast_3|Non-myelinated Schw",
                          "Fibroblast_4|Non-myelinated Schw",
                          "Fibroblast_5|Non-myelinated Schw",
                          "Fibroblast_6|Non-myelinated Schw",
                          "Fibroblast_1|Myelinated Schw",
                          "Fibroblast_2|Myelinated Schw",
                          "Fibroblast_3|Myelinated Schw",
                          "Fibroblast_4|Myelinated Schw",
                          "Fibroblast_5|Myelinated Schw",
                          "Fibroblast_6|Myelinated Schw",
                          "Fibroblast_1|Endothelial_1",
                          "Fibroblast_2|Endothelial_1",
                          "Fibroblast_3|Endothelial_1",
                          "Fibroblast_4|Endothelial_1",
                          "Fibroblast_5|Endothelial_1",
                          "Fibroblast_6|Endothelial_1",
                          "Fibroblast_1|Endothelial_2",
                          "Fibroblast_2|Endothelial_2",
                          "Fibroblast_3|Endothelial_2",
                          "Fibroblast_4|Endothelial_2",
                          "Fibroblast_5|Endothelial_2",
                          "Fibroblast_6|Endothelial_2",
                          "Fibroblast_1|Endothelial_3",
                          "Fibroblast_2|Endothelial_3",
                          "Fibroblast_3|Endothelial_3",
                          "Fibroblast_4|Endothelial_3",
                          "Fibroblast_5|Endothelial_3",
                          "Fibroblast_6|Endothelial_3"
                          ), 
                  proteins = PTN_gene, cluster_proteins=TRUE) + scale_fill_viridis() #+ coord_flip()

