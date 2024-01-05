# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the data --------------------------------------------------------
# mat <- fread("../data/exprMatrix.tsv.gz")
mat <- read.table("data/exprMatrix.tsv", header=T, sep="\t", as.is=T, row.names=1)

meta <- read.table("data/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
# fix the rownames of the metadata compared to the genenames of the matrix
rownames(meta) <- str_replace(rownames(meta),pattern = "-",replacement = "\\.")

# confirm all the barcodes 
sum(!(colnames(mat)==rownames(meta)))

# create the object -------------------------------------------------------
so <- CreateSeuratObject(counts = mat, project = "Tanaka_organoids", meta.data=meta)

# save the seurat object
saveRDS(object = so,file = "out/object/Tanaka_organoids_init.rds")

so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so) 
so <- RunPCA(so) 
so <- FindNeighbors(so,dims = 1:30) 
so <- FindClusters(so) 
so <- RunUMAP(so,dims = 1:30) 

DimPlot(so)
DimPlot(so,group.by = "Cluster_name",label = T)
ggsave("out/image/UMAP_Tanaka_organoids_full_reanalyzed.pdf",width = 7,height = 6)
saveRDS(object = so,file = "out/object/Tanaka_organoids_full_reanalyzed.rds")

# add the plot following martina's suggestion -----------------------------
# add the plot with no label on the cluster
DimPlot(so,group.by = "Cluster_name",label = T)
ggsave("out/image/UMAP_Tanaka_organoids_full_reanalyzed_ClusterName_Label.pdf",width = 7,height = 6)

DimPlot(so,group.by = "Cluster_name")
ggsave("out/image/UMAP_Tanaka_organoids_full_reanalyzed_ClusterName_NoLabel.pdf",width = 7,height = 6)

DimPlot(so,group.by = "seurat_clusters",label = T)
ggsave("out/image/UMAP_Tanaka_organoids_full_reanalyzed_SeuratClusters_Label.pdf",width = 7,height = 6)

DimPlot(so,group.by = "seurat_clusters")
ggsave("out/image/UMAP_Tanaka_organoids_full_reanalyzed_SeuratClusters_NoLabel.pdf",width = 7,height = 6)

# split the sample by the factor martina is interested in 
table(so@meta.data$orig.ident)
so$split_martina <- case_when(so@meta.data$orig.ident == "Fetal"~"Fetal",
          T ~"Organoids")
DimPlot(so,group.by = "Cluster_name",split.by = "split_martina")
ggsave("out/image/UMAP_Tanaka_organoids_full_reanalyzed_ClusterName_noLabels_split.pdf",width = 14,height = 6)

