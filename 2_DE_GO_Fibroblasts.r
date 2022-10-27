---
Title: "Analysis of dental pulp dataset"
Author: "Diego Diez"
date: "2022-09-26"
---


library(here)
library(tidyverse)
library(patchwork)
library(Seurat)
library(scmisc)

theme_set(cowplot::theme_cowplot())


# Load data
x <- read_rds(here("data/SeuratInteg_CR7_Ref2020-A_6tis_2nd.rds"))
x


DefaultAssay(x) <- "RNA"
x


ggplot(x[[]], aes(dataset)) +
  geom_bar() + RotatedAxis()

ggplot(x[[]], aes(dataset, fill=seurat_clusters)) +
  geom_bar(position="fill") + RotatedAxis()

ggplot(x[[]], aes(seurat_clusters, fill=dataset)) +
  geom_bar(position="fill") + RotatedAxis()

plot_coord2(x, expand="dataset", label=TRUE) & NoAxes() + NoLegend()

plot_coord2(x, expand="tissue", label=TRUE) & NoAxes() + NoLegend()

VlnPlot(x, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), pt.size=.1, ncol=2, assay="RNA")

FeaturePlot(x, "percent.mt", label=TRUE) + NoLegend() + NoAxes()



x <- CellCycleScoring(x, cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)

FeaturePlot(x, c("S.Score", "G2M.Score"), order=TRUE, blend=TRUE, label=TRUE)[[3]] + NoAxes()

DimPlot(x, label=TRUE) + NoLegend() + NoAxes()

plot_coord2(x, expand="seurat_clusters", label=TRUE) & NoAxes() + NoLegend()

ggplot(x[[]], aes(seurat_clusters)) +
  geom_bar() + RotatedAxis()


--------------------------------------------------------------------------------------------------------
  
# Find DEG between clusters

deg <- FindAllMarkers(x, max.cells.per.ident=500, verbose=FALSE, assay="RNA")
head(deg)

plot_volcano(deg, max.overlaps=Inf)

top_genes <- top_deg(deg, fdr=0.05, lfc=1, n=20)
top_genes

y <- sample_cells(x, "seurat_clusters", n_max=50)

plot_heatmap(y[top_genes$gene, ], column_split=y$seurat_clusters)


## Enrichment analysis

db <- org.Hs.eg.db::org.Hs.eg.db
deg$entrezgene <- AnnotationDbi::mapIds(db, keys=deg$gene, column="ENTREZID", keytype="SYMBOL")
head(deg)

res_kegg <- run_enrichment(deg, type="kegg", org="Hs")
head(res_kegg)

lapply(unique(res_kegg$cluster), function(cluster) {
  plot_enrichment_barplot(res_kegg |> filter(cluster==!!cluster)) +
    labs(title=paste0("Cluster: ", cluster))
}) |> wrap_plots()

res_go <- scmisc::run_enrichment(deg, type="go", org="Hs")
head(res_go)

lapply(unique(res_go$cluster), function(cluster) {
  plot_enrichment_barplot(res_go |> filter(cluster==!!cluster, Ont=="BP")) +
    labs(title=paste0("Cluster: ", cluster, "; Ontology: BP"))
}) |> wrap_plots()

lapply(unique(res_go$cluster), function(cluster) {
  plot_enrichment_barplot(res_go |> filter(cluster==!!cluster, Ont=="MF"))+
    labs(title=paste0("Cluster: ", cluster, "; Ontology: MF"))
}) |> wrap_plots()

lapply(unique(res_go$cluster), function(cluster) {
  plot_enrichment_barplot(res_go |> filter(cluster==!!cluster, Ont=="CC"))+
    labs(title=paste0("Cluster: ", cluster, "; Ontology: CC"))
}) |> wrap_plots()


--------------------------------------------------------------------------------------------------------
  
FIGURE 2A and 2B

# Find DEG between DP and rest in fibroblasts

deg_dtp <- FindMarkers(x, ident.1="DTP", group.by="tissue", subset.ident=c("0", "2", "3", "7"), max.cells.per.ident=500, verbose=FALSE, assay="RNA")
deg_dtp <- deg_dtp |> 
  rownames_to_column("gene") |>
  mutate(cluster="DTP vs. Others (Fibroblasts)")
head(deg_dtp)

plot_volcano(deg_dtp, n=20, max.overlaps=Inf)

y <- x[, x$seurat_clusters %in% c("0", "2", "3", "7")]
y

y$tissue_simple <- "Other"
y$tissue_simple[y$tissue == "DTP"] <- "DTP"

deg_dtp

VlnPlot(y, c("PTN", "MDK"), split.by="tissue_simple", assay="RNA")

VlnPlot(y, head(deg_dtp$gene, 6), split.by="tissue_simple", assay="RNA")


## Enrichment analysis

db <- org.Hs.eg.db::org.Hs.eg.db
deg_dtp$entrezgene <- AnnotationDbi::mapIds(db, keys=deg_dtp$gene, column="ENTREZID", keytype="SYMBOL")
head(deg_dtp)

res_kegg <- run_enrichment(deg_dtp, type="kegg", org="Hs", FDR=0.05)
limma::topKEGG(res_kegg)

plot_enrichment_barplot(res_kegg)

res_go <- scmisc::run_enrichment(deg_dtp, type="go", org="Hs", FDR=0.05)
head(res_go)

plot_enrichment_barplot(res_go |> filter(Ont == "BP"))
plot_enrichment_barplot(res_go |> filter(Ont == "MF"))
plot_enrichment_barplot(res_go |> filter(Ont == "CC"))

