# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(monocle)

# read in the object ------------------------------------------------------
# the input dataset is the one out of the BS_drop_RR25CTRL
data.combined <- readRDS("data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
Idents(data.combined) <- "orig.ident"
# test_CTRL <- subset(data.combined,idents = c("pool_MOCK_CTRL","RR16_SOX10_CTRL"))
test_CTRL <- subset(data.combined,subset = orig.ident %in% c("pool_MOCK_CTRL","RR16_SOX10_CTRL") & seurat_clusters !=12 & seurat_clusters != 13)
DimPlot(test_CTRL)
DefaultAssay(test_CTRL) <- "RNA"
DimPlot(test_CTRL,split.by="clone",group.by = "seurat_clusters",label = T)
ggsave("out/image/UMAP_POOLandRR16_alt_split.pdf",height = 5,width = 9)

# convert he object ot CellDataSet
monocle_df <- as.CellDataSet(test_CTRL)
monocle_df <- estimateSizeFactors(monocle_df)
monocle_df <- estimateDispersions(monocle_df)

# Trajectory step 1: choose genes that define a cell's progress -----------
# read in the markers genes
ordering_genes <- c("CLU", "SPARCL1","B2M","GFAP", "VIM","SLC1A2","AQP4",
                    "EGR1","TOP2A", "MKI67", "PBK", "NUSAP1", "CENPU",
                    "SEZ6L2", "NSG2", "SYT4", "SERPINI1", "STMN2", "BASP1", "MAP1B","GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1","NEUROD6", "NEUROD1", "NEUROD2", "NHLH1", "BHLHE22",
                    "PLP1","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG","S100B")

# check the marker genes
Idents(test_CTRL) <- "seurat_clusters"
DotPlot(test_CTRL,features = ordering_genes)+theme(axis.text.x = element_text(hjust = 1,angle = 45))

setdiff(ordering_genes2,ordering_genes)
setdiff(ordering_genes,ordering_genes2)

monocle_df <- setOrderingFilter(monocle_df, ordering_genes)
plot_ordering_genes(monocle_df)

# Trajectory step 2: reduce data dimensionality ---------------------------
# Next, we will reduce the space down to one with two dimensions, which we will be able to easily visualize and interpret while Monocle is ordering the cells.
monocle_df <- reduceDimension(monocle_df, max_components = 3,method = 'DDRTree')

# Trajectory step 3: order cells along the trajectory ---------------------
# Now that the space is reduced, it's time to order the cells using the orderCells function as shown below.
# monocle_df <- orderCells(monocle_df)
monocle_df <- orderCells(monocle_df)

plot_cell_trajectory(monocle_df, color_by = "seurat_clusters")
plot_cell_trajectory(monocle_df, color_by = "State")
plot_cell_trajectory(monocle_df, color_by = "Pseudotime")

# change the root state. the criterion for the root is the identification of the positino of the proliferating cells. in this case the cluster 7 cells
test <- plot_cell_trajectory(monocle_df, color_by = "seurat_clusters")
test$data %>% 
  ggplot(aes(x=data_dim_1,y=data_dim_2))+geom_point()+facet_wrap(~seurat_clusters)

# the source of cluster 7 cells seems to be located into the state. therefore I will put the root there. this will affect the ralative value of the pseudotime
monocle_df <- orderCells(monocle_df, root_state = 2)

plot_cell_trajectory(monocle_df, color_by = "Pseudotime")
plot_cell_trajectory(monocle_df, color_by = "seurat_clusters")
plot_cell_trajectory(monocle_df, color_by = "State")

# saveRDS(monocle_df,"out/monocle_df_RR16_both.rds")
saveRDS(monocle_df,"out/object/monocle_df_POOLandRR16_CTRL.rds")

# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
plot_cell_trajectory(monocle_df, color_by = "orig.ident")

plot_cell_trajectory(monocle_df, color_by = "seurat_clusters") +
  facet_wrap(~seurat_clusters, nrow = 3)

plot_cell_trajectory(monocle_df, color_by = "seurat_clusters") +
  facet_wrap(~orig.ident, nrow = 1)


# -------------------------------------------------------------------------
marker_genes <- row.names(subset(fData(monocle_df),
                                 gene_short_name %in% ordering_genes))

diff_test_res <- differentialGeneTest(monocle_df[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))

pdf("out/image/heatmap_differentialGeneTest_monocle_df_POOLandRR16_CTRL.pdf",width = 5,height = 10)
plot_pseudotime_heatmap(monocle_df[sig_gene_names,],
                        num_clusters = 5,
                        cores = 2,
                        show_rownames = T)
dev.off()
