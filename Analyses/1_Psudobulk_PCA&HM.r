---
Title: "Pseudobulk_PCA_HM"
Author: Natnicha Jiravejchakul
Date: '2022-07-11'
---

  
library(tidyverse)
library(reshape2)
library(ggplot2)
library(parallelDist)
library(dplyr)
library(ComplexHeatmap)
library(viridis)


Integ_new_CR7_seurat <- readRDS(file = "~/path_to_directory")
Integ_new_CR7_seurat



============================ PSEUDOBULK PCA AND HEATMAP OF EACH CLUSTERS ============================
                              FIGURE 1E, 1F AND SUPPLEMENTARY FIGURE 5

  
  
1. PCA

DefaultAssay(Integ_new_CR7_seurat) <- "RNA"
Idents(Integ_new_CR7_seurat) <- "seurat_clusters"
levels(Integ_new_CR7_seurat)

Integ_new_CR7_Cluster0 <- subset(x = Integ_new_CR7_seurat, idents = "0") # Change the cluster of interest here 
# In figure 1D, we combined cluster 0,2,3,7 (Fibroblasts) together using idents = c("0", "2", "3", "7")

# Ckeck which datasets were composed in this cluster
Integ_new_CR7_Cluster0
DimPlot(Integ_new_CR7_Cluster0, split.by = "dataset")

Idents(Integ_new_CR7_Cluster0) <- "dataset" 
Average_data_Integ_new_CR7_Cluster0 <- AverageExpression(Integ_new_CR7_Cluster0 , assays = "RNA" , slot = "data")

# PCA Calculation
pca <- prcomp(t(Average_data_Integ_new_CR7_Cluster0$RNA))
pca_perc <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
df_pca_Integ_new_CR7_Cluster0 <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], 
                                 sample = colnames(Average_data_Integ_new_CR7_Cluster0$RNA))

# Add and order the data sets
df_pca_Integ_new_CR7_Cluster0$dataset <- c("DTP1", "DTP2", "DTP3", "DTP4", "DTP5", "DTP6", "DTP7", "DTP8", "DTP9",
                                           "PBMC1", "PBMC2", "BM1", "BM2", 
                                           "ADP1", "ADP2", "LUNG1", "LUNG2", 
                                           "SKIN1", "SKIN2") 
df_pca_Integ_new_CR7_Cluster0$dataset <- factor(df_pca_Integ_new_CR7_Cluster0$dataset , 
                                         levels = c("DTP1", "DTP2", "DTP3", "DTP4", "DTP5", "DTP6", "DTP7", "DTP8", "DTP9",
                                                    "PBMC1", "PBMC2", "BM1", "BM2", 
                                                    "ADP1", "ADP2", "LUNG1", "LUNG2", 
                                                    "SKIN1", "SKIN2")) 

# Plot PCA
pca_plot_Integ_new_CR7_Cluster0 <-  ggplot(df_pca_Integ_new_CR7_Cluster0, aes(PC1,PC2, color = dataset)) + 
  geom_point(colour = c("#CC9900B2", "#CC9900B2", "#CC9900B2", "#CC9900B2",                 # DP
                        "#CC9900B2", "#CC9900B2", "#CC9900B2", "#CC9900B2", "#CC9900B2",                                                                               "#BA6338B2", "#BA6338B2",                                           # PBMC
                        "#466983B2", "#466983B2",                                           # BM
                        "#5A655EB2", "#5A655EB2",                                           # ADP
                        "#802268B2", "#802268B2",                                           # LUNG
                        "#749B58B2", "#749B58B2"), size = 5) +                              # SKIN
  labs(x=paste0("PC1 (",pca_perc[1],"%)"), y=paste0("PC2 (",pca_perc[2],"%)")) +  
  theme(axis.text = element_text(size = 13, face="bold", colour = "black"), 
        axis.title.y = element_text(color="black", size=13, face="bold"),  
        axis.title.x  = element_text(color="black", size=13, face="bold"),  
        legend.title = element_text(face = "bold" , size = 13),  
        legend.text = element_text(size = 13), legend.key.size = unit(0.5, "cm"), 
        legend.key.width = unit(0.5,"cm"), legend.key = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(pca_plot_Integ_new_CR7_Cluster0)


----------------------------------------------------------------------------------------------------
  

2. Heat Map 
   # HM was constructed using list of top 500 genes that highly contributed to the PC1 and PC2 of the PCA


# Retrieve list of 500 genes (loading) from PC1 and PC2 
PC1_Genes <- data.frame(sort(abs(pca$rotation[,"PC1"]), decreasing=TRUE)[1:500])
PC2_Genes <- data.frame(sort(abs(pca$rotation[,"PC2"]), decreasing=TRUE)[1:500])

# Union of the PC1 genes and PC2 genes
union_genes_Integ_new_CR7_seurat <- union(rownames(PC1_Genes) , rownames(PC2_Genes))


# Calculate Z-score of the those genes using the average expression level

input_heatmap_Integ_new_CR7_seurat <- Average_data_Integ_new_CR7_Cluster0$RNA[rownames(Average_data_Integ_new_CR7_Cluster0$RNA) %in% union_genes_Integ_new_CR7_seurat,]

input_heatmap_Integ_new_CR7_seurat <- t(apply((input_heatmap_Integ_new_CR7_seurat), 1, function(x){
  mean <- mean(x)
  SD <- sd(x)
  Z_score <- (x-mean)/SD
  Z_score
}))

input_heatmap_Integ_new_CR7_seurat <- as.data.frame(input_heatmap_Integ_new_CR7_seurat)


# Construct the HM

parallelDist_dtw <- hclust(parDist(x = as.matrix(input_heatmap_Integ_new_CR7_seurat), method = "dtw") , method="ward")
parallelDist_dtw_log2_de_or <-input_heatmap_Integ_new_CR7_seurat[parallelDist_dtw$order,] 

ggplot(parallelDist_dtw$height %>% as.tibble() %>% add_column(groups = length(parallelDist_dtw$height):1), 
       aes(x=groups, y=value)) +geom_point() + geom_line() + xlim(0,20)

row_cutree <- data.frame(cutree(parallelDist_dtw, k = 4)) 
# Adjust number of k based on data 
# Number of k correlates with number of gene clusters showed on the heatmap

colnames(row_cutree) <- "cluster"

parallelDist_dtw_log2_de_or <- parallelDist_dtw_log2_de_or %>%  rownames_to_column %>%  
                               right_join( row_cutree %>%  rownames_to_column %>%  
                               select(cluster ,rowname ) , by = "rowname")
parallelDist_dtw_log2_de_or <- parallelDist_dtw_log2_de_or[parallelDist_dtw$order,]%>% 
                               cbind.data.frame(index=seq(1:nrow(.)))

out_heatmap_Integ_new_CR7_seurat <- Heatmap(parallelDist_dtw_log2_de_or %>% remove_rownames() %>% 
                                    column_to_rownames() %>% select(-cluster,-index),  
                                    col = circlize::colorRamp2(c(-3,0,3), c("blue","black","orange")), 
                                    name = "Average Gene Expression (Z-score)", 
                                    cluster_columns = F, show_column_dend = F, 
                                    cluster_rows = F, show_row_names = T, 
                                    show_column_names = T, na_col = "grey" , 
                                    row_names_gp = gpar(fontsize = 1))

draw(out_heatmap_Integ_new_CR7_seurat, row_split = parallelDist_dtw_log2_de_or$cluster)


----------------------------------------------------------------------------------------------------
  

3. Retrieve the set of genes from the heatmap 

HeapmapCluster_Gene_list <-list(cluster1=row_cutree %>% rownames_to_column() %>%     
                                    filter(cluster=="1")%>%select(rowname)%>%unlist%>%as.character(), 
                                cluster2=row_cutree %>% rownames_to_column() %>% 
                                    filter(cluster=="2")%>%select(rowname)%>%unlist%>%as.character(), 
                                cluster3=row_cutree %>% rownames_to_column() %>% 
                                    filter(cluster=="3")%>%select(rowname)%>%unlist%>%as.character(), 
                                cluster4=row_cutree %>% rownames_to_column() %>% 
                                    filter(cluster=="4")%>%select(rowname)%>%unlist%>%as.character())

HeapmapCluster_Gene_list[1] # The number inside [] indicate the gene clusters of interest


----------------------------------------------------------------------------------------------------
  

4. Dot plot shows expression of selected genes (Highly express in DP) retrieved from the heat map - FIGURE 1G


Idents(Integ_new_CR7_fibroblast) <- "dataset"
DotPlot(Integ_new_CR7_fibroblast, features = c( "IFI6",         "NES",      "KCNT2",      "KLHL29",   "CYP1B1",     "NRXN1",   
                                                "THSD7B",       "LRP1B",    "CHN1",       "SATB2",    "CNTN4",      "CCK",    
                                                "ADAMTS9-AS2",  "MAGI1",    "LSAMP",      "KALRN",    "TF",         "NSG1",  
                                                "LIMCH1",       "SLC4A4",   "TMEM150C",   "BMPR1B",   "GALNTL6",    "SCRG1",    
                                                "SORBS2",       "SEMA5A",   "CDH12",      "MSX2",     "RUNX2",      "COL21A1",  
                                                "AL445250.1",   "DYNC1I1",  "DLX6-AS1",   "PTN",      "SLC20A2",    "CPA6",  
                                                "PTPRD",        "MOB3B",    "TNC",        "KIAA1217", "BAMBI",      "DKK3",  
                                                "INSC",         "MDK",      "FAT3",       "PDGFD",    "AC007368.1", "PCDH9",  
                                                "MYO16",        "NRXN3",    "UACA",       "CRABP1",   "CDH11",      "TIAM1",  
                                                "GPM6B" ), # Selected highly expressed genes in DP from gene cluster [1] of the heat map
  assay = "RNA")  + RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option = "cividis", direction = 1) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
    



