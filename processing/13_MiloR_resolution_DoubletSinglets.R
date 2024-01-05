# libraries ---------------------------------------------------------------
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

# 01 Introduction ---------------------------------------------------------
# Milo is a tool for analysis of complex single cell datasets generated from replicated multi-condition experiments, which detects changes in composition between conditions. While differential abundance (DA) is commonly quantified in discrete cell clusters, Milo uses partially overlapping neighbourhoods of cells on a KNN graph. Starting from a graph that faithfully recapitulates the biology of the cell population, Milo analysis consists of 3 steps:
  
# Sampling of representative neighbourhoods
# Testing for differential abundance of conditions in all neighbourhoods
# Accounting for multiple hypothesis testing using a weighted FDR procedure that accounts for the overlap of neighbourhoods
# In this vignette we will elaborate on how these steps are implemented in the miloR package.

# 02 Load data ------------------------------------------------------------
data.combined <- readRDS(file = "data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
sce <- as.SingleCellExperiment(data.combined)
milo <- Milo(sce)
# add the graph to the object
miloR::graph(milo) <- miloR::graph(buildFromAdjacency(data.combined@graphs$RNA_snn, k=10))

# 03 Pre-processing -------------------------------------------------------
# For DA analysis we need to construct an undirected KNN graph of single-cells. Standard single-cell analysis pipelines usually do this from distances in PCA. We normalize and calculate principal components using scater. I also run UMAP for visualization purposes.
# the dataset is already preprocessed
plotUMAP(milo)
# get some more metadata to plot
colData(milo)

plotUMAP(milo,colour_by="seurat_clusters")

# 06 1. Defining representative neighbourhoods ----------------------------
# We define the neighbourhood of a cell, the index, as the group of cells connected by an edge in the KNN graph to the index cell. For efficiency, we don’t test for DA in the neighbourhood of every cell, but we sample as indices a subset of representative cells, using a KNN sampling algorithm used by Gut et al. 2015.

# For sampling you need to define a few parameters:
  
# prop: the proportion of cells to randomly sample to start with (usually 0.1 - 0.2 is sufficient)
# k: the k to use for KNN refinement (we recommend using the same k used for KNN graph building)
# d: the number of reduced dimensions to use for KNN refinement (we recommend using the same d used for KNN graph building)
# refined indicated whether you want to use the sampling refinement algorithm, or just pick cells at random. The default and recommended way to go is to use refinement. The only situation in which you might consider using random instead, is if you have batch corrected your data with a graph based correction algorithm, such as BBKNN, but the results of DA testing will be suboptimal.
milo <- makeNhoods(milo, prop = 0.1, k = 10, d=30, refined = TRUE)

# Once we have defined neighbourhoods, it’s good to take a look at how big the neighbourhoods are (i.e. how many cells form each neighbourhood). This affects the power of DA testing. We can check this out using the plotNhoodSizeHist function. Empirically, we found it’s best to have a distribution peaking between 50 and 100. Otherwise you might consider rerunning makeNhoods increasing k and/or prop (here the distribution looks ludicrous because it’s a small dataset).
plotNhoodSizeHist(milo)

# 07 Counting cells in neighbourhoods -------------------------------------
# Now we have to count how many cells from each sample are in each neighbourhood. We need to use the cell metadata and specify which column contains the sample information.
milo <- countCells(milo, meta.data = data.frame(colData(milo)), samples="sample")

# This adds to the Milo object a n \times m matrix, where n is the number of neighbourhoods and m is the number of experimental samples. Values indicate the number of cells from each sample counted in a neighbourhood. This count matrix will be used for DA testing.
head(nhoodCounts(milo))

# 08 Differential abundance testing ---------------------------------------
# Now we are all set to test for differential abundance in neighbourhoods. We implement this hypothesis testing in a generalized linear model (GLM) framework, specifically using the Negative Binomial GLM implementation in edgeR.

# We first need to think about our experimental design. The design matrix should match samples to a condition of interest. In this case the doxy is the covariate we are going to test for.
traj_design <- data.frame(colData(milo))[,c("sample", "doxy")]
# comparison tested
table(colData(milo)$doxy)
# comparision in the original dataset
table(data.combined$doxy)
# potential alternative comparison with more than one level
table(data.combined$group_martina)

traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample
## Reorder rownames to match columns of nhoodCounts(milo)
traj_design <- traj_design[colnames(nhoodCounts(milo)), , drop=FALSE]

traj_design

# Milo uses an adaptation of the Spatial FDR correction introduced by cydar, which accounts for the overlap between neighbourhoods. Specifically, each hypothesis test P-value is weighted by the reciprocal of the kth nearest neighbour distance. To use this statistic we first need to store the distances between nearest neighbors in the Milo object.
milo <- calcNhoodDistance(milo, d=30)

# Now we can do the test, explicitly defining our experimental design.
rownames(traj_design) <- traj_design$sample
da_results <- testNhoods(milo, design = ~doxy, design.df = traj_design,fdr.weighting = "graph-overlap")

# This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates whether there is significant differential abundance between conditions.
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

# 09 Visualize neighbourhoods displaying DA -------------------------------
# To visualize DA results, we build an abstracted graph of neighbourhoods that we can superimpose on the single-cell embedding.
milo <- buildNhoodGraph(milo)

plotUMAP(milo,colour_by="Phase") + plotNhoodGraphDA(milo, da_results, alpha=1) +
  plot_layout(guides="collect")
ggsave("out/image/miloR_test.pdf",width = 13,height = 6)

## -----------------------------------------------------------------------------
saveRDS(milo,"out/object/miloR.rds")
saveRDS(da_results,"out/object/da_results.rds")
sessionInfo()

# change color of the plot ------------------------------------------------
milo <- readRDS("out/object/miloR.rds")
da_results <- readRDS("out/object/da_results.rds")
#
test <- plotNhoodGraphDA(milo, da_results, alpha=1)
test+scale_fill_viridis_c(option = "turbo")

test + 
  scale_fill_gradient2(
    name="log2FC",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill")

ggsave("out/image/miloR_test_update.pdf",width = 7,height = 6)

test$data %>% 
  ggplot(aes(x=x,y=y))+
  geom_point(aes(fill=colour_by,size = size),shape=21,alpha=0.5)+scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )+theme_bw()
ggsave("out/image/miloR_test_update2.pdf",width = 7,height = 6)

# simulate a volcano for the differential aboundance of the neighbours
da_results %>% 
  ggplot(aes(x=logFC,y=-log(FDR)))+geom_point(alpha=0.5)+theme_bw()+geom_hline(yintercept = -log(0.05),linetype="dashed",col="gray")
ggsave("out/image/miloR_volcano_neighbour_FDR.pdf",width = 4,height = 4)

da_results %>% 
  ggplot(aes(x=logFC,y=-log(SpatialFDR)))+geom_point(alpha=0.5)+theme_bw()+geom_hline(yintercept = -log(0.05),linetype="dashed",col="gray")
ggsave("out/image/miloR_volcano_neighbour_spatialFDR.pdf",width = 4,height = 4)

da_results %>% 
  ggplot(aes(x=logFC,y=-log(PValue)))+geom_point(alpha=0.5)+theme_bw()+geom_hline(yintercept = -log(0.05),linetype="dashed",col="gray")
ggsave("out/image/miloR_volcano_neighbour_pvalue.pdf",width = 4,height = 4)
