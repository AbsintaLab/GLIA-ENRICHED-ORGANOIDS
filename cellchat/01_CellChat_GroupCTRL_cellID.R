# Load the required libraries ---------------------------------------------
library(patchwork)
library(tidyverse)
library(CellChat)
library(Matrix)
library(NMF)
library(ggalluvial)
library(Seurat)
options(stringsAsFactors = FALSE)

# set up the objects ------------------------------------------------------
# read in the seurat object
seurat_GroupCTRL <- readRDS("../../out/object/subset_GroupCCC_CTRL_data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
# add a cellID classification following Martina's suggestion
ref_meta <- seurat_GroupCTRL@meta.data %>%
  rownames_to_column("barcode") %>% 
  mutate(cellID = case_when(seurat_clusters %in% c(16)~"OLIGO",
                            seurat_clusters %in% c(2)~"OPC",
                            seurat_clusters %in% c(4,6,9)~"ASTRO",
                            seurat_clusters %in% c(7)~"CYCLING",
                            seurat_clusters %in% c(12,13)~"iMIC",
                            seurat_clusters %in% c(0,1,3,5,10,11,15)~"NEU",
                            seurat_clusters %in% c(8,14)~"GLIA"))

# update the meta also in Martina's object
seurat_GroupCTRL$cellID <- ref_meta$cellID
# update the Idents
Idents(seurat_GroupCTRL) <- "cellID"

#
DimPlot(seurat_GroupCTRL,label = T)+guides(colour = guide_legend(override.aes = list(size=5),ncol=1))
ggsave("../../out/image/UMAP_seurat_GroupCTRLCellID.pdf",width = 5,height = 3)

DimPlot(seurat_GroupCTRL,label = T,split.by = "orig.ident")+guides(colour = guide_legend(override.aes = list(size=5),ncol=1))
ggsave("../../out/image/UMAP_seurat_GroupCTRLCellID_split.pdf",width = 10,height = 3)

df_summary <- seurat_GroupCTRL@meta.data %>%
  # mutate(cluster_fix = paste0("clu_",str_pad(seurat_clusters,width = 2,side = "left",pad = 0))) %>%
  group_by(cellID,orig.ident) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  mutate(prop = n/tot) %>%
  mutate(dataset = "CTRL")

df_summary %>%
  write_tsv("../../out/table/prop_seurat_GroupCTRLCellID.tsv")

# fix the cluster name
seurat_GroupCTRL$cluster_fix <- paste0("clu_",str_pad(seurat_GroupCTRL$seurat_clusters,width = 2,side = "left",pad = "0"))

# Create a CellChat object ------------------------------------------------
cellchat_GroupCTRL <- createCellChat(object = seurat_GroupCTRL, group.by = "cellID",assay = "RNA")

# number of cells in each cell group
groupSize <- table(cellchat_GroupCTRL@idents)

# Set the ligand-receptor interaction database ----------------------------
# use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# ggsave("../../out/image/CellChatDB.pdf",width = 5,height = 5)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB$interaction %>%
  group_by(annotation) %>%
  summarise()

# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat_GroupCTRL@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication an --------
# To infer the cell state-specific communications, we identify over-expressed ligands or receptors in one cell group and then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
cellchat_GroupCTRL <- subsetData(cellchat_GroupCTRL)

cellchat_GroupCTRL <- identifyOverExpressedGenes(cellchat_GroupCTRL)
cellchat_GroupCTRL <- identifyOverExpressedInteractions(cellchat_GroupCTRL)
# project gene expression data onto PPI network (optional)
cellchat_GroupCTRL <- projectData(cellchat_GroupCTRL, PPI.human)

# Part II: Inference of cell-cell communication network -------------------
# CellChat infers the biologically significant cell-cell communication by assigning each interaction with a probability value and peforming a permutation test. CellChat models the probability of cell-cell communication by integrating gene expression with prior known knowledge of the interactions between signaling ligands, receptors and their cofactors using the law of mass action.

# Compute the communication probability and infer cellular communi --------
# if there is an error in this step try to check:
# 1) Use the correct CellChatDB (mouse or human)
# 2) the input data matrix cellchat@data is present and non 0 values
# 3) the data matrix cellchat@data.signaling is present and non 0 values
# 4) the cell group information is correct: unique(cellchat@idents): 1 to 11.
# 5) try to change the name of the idents
# population.size: whether consider the proportion of cells in each group across all sequenced cells. Set population.size = FALSE if analyzing sorting-enriched single cells, to remove the potential artifact of population size. Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, with the reason that abundant cell populations tend to send collectively stronger signals than the rare cell populations.
cellchat_GroupCTRL <- computeCommunProb(cellchat_GroupCTRL,population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_GroupCTRL <- filterCommunication(cellchat_GroupCTRL, min.cells = 10)

# Extract the inferred cellular communication network as a data fr --------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,

# Infer the cell-cell communication at a signaling pathway level ----------
# CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.
cellchat_GroupCTRL <- computeCommunProbPathway(cellchat_GroupCTRL)
cellchat_GroupCTRL@net
cellchat_GroupCTRL@netP

# Calculate the aggregated cell-cell communication network ----------------
# We can calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.

cellchat_GroupCTRL <- aggregateNet(cellchat_GroupCTRL)
# save the object up to this poitn
# saveRDS(cellchat,file = "cellchat.rds")

# try to plot some communications
groupSize <- as.numeric(table(cellchat_GroupCTRL@idents))
# par(mfrow = c(1,2), xpd=TRUE)

pdf("../../out/image/01_number_of_interactions_cellchat_GroupCTRLCellID.pdf",width = 10,height = 10)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat_GroupCTRL@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("../../out/image/01_strength_of_interactions_cellchat_GroupCTRLCellID.pdf",width = 10,height = 10)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat_GroupCTRL@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# focus on the IMMUNE cells communication
mat <- cellchat_GroupCTRL@net$weight

pdf("../../out/image/01_grid_interactions_cellchat_GroupCTRLCellID.pdf",width = 15,height = 15)
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = F, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# partIII: visualization of ccc networks ----------------------------------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
cellchat_GroupCTRL <- netAnalysis_computeCentrality(cellchat_GroupCTRL, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_GroupCTRL, pattern = "outgoing",height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_GroupCTRL, pattern = "incoming",height = 15)

pdf("../../out/image/02_contribution_cellchat_GroupCTRL.pdf",width = 11,height = 11)
ht1 + ht2
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_GroupCTRL)
#> Signaling role analysis on the cell-cell communication network from user's input
gg1
ggsave("../../out/image/02_contribution_all_cellchat_GroupCTRLCellID.pdf",width = 5,height = 5)

# -------------------------------------------------------------------------
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2_SPP1 <- netAnalysis_signalingRole_scatter(cellchat_GroupCTRL, signaling = c("SPP1"))

gg2_SPP1
ggsave("../../out/image/02_contribution_subset_SPP1_cellchat_GroupCTRLCellID.pdf",width = 5,height = 5)
# gg2_BMP
# ggsave("out/image/02_contribution_subset_BMP_cellchat_GroupCTRLCellID.pdf",width = 5,height = 5)
# confirm the evdence by looking at the average expression
computeAveExpr(cellchat_GroupCTRL, features = c("SPP1","ITGAV","ITGB1","ITGA4","ITGB5"),group.by = c("cluster_fix"))

#
pathways.show <- c("SPP1")

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,3) # a numeric vector. 
pdf("../../out/image/02_contribution_hierarchy_SPP1_cellchat_GroupCTRLCellID.pdf",width = 10,height = 7)
lapply(pathways.show, function(x){
  netVisual_aggregate(cellchat_GroupCTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "hierarchy")
})
dev.off()
# also other approach in terms of aesthtuchs
vertex.receiver = seq(1,3) # a numeric vector. 
pdf("../../out/image/02_contribution_chord_SPP1_cellchat_GroupCTRLCellID.pdf",width = 10,height = 7)
lapply(pathways.show, function(x){
  netVisual_aggregate(cellchat_GroupCTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "chord")
})
dev.off()
#
vertex.receiver = seq(1,3) # a numeric vector. 
pdf("../../out/image/02_contribution_circle_SPP1_cellchat_GroupCTRLCellID.pdf",width = 10,height = 7)
lapply(pathways.show, function(x){
  netVisual_aggregate(cellchat_GroupCTRL, signaling = x,  vertex.receiver = vertex.receiver,layout = "circle")
})
dev.off()

# Heatmap
pdf("../../out/image/02_contribution_heatmap_SPP1_cellchat_GroupCTRLCellID.pdf",width = 10,height = 7)
lapply(pathways.show, function(x){
  netVisual_heatmap(cellchat_GroupCTRL, signaling = x, color.heatmap = "Reds")
})
dev.off()
# par(mfrow=c(1,1))
#> Do heatmap based on a single object

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat_GroupCTRL, signaling = pathways.show)
ggsave("../../out/image/02_contribution_braplotLR_SPP1_cellchat_GroupCTRLCellID.pdf",width = 5,height = 5)

#
pairLR <- extractEnrichedLR(cellchat_GroupCTRL, signaling = pathways.show, geneLR.return = FALSE)

# Hierarchy plot
vertex.receiver = seq(1,3) # a numeric vector
pdf("../../out/image/02_contribution_hierarcyLR_SPP1_cellchat_GroupCTRLCellID.pdf",width = 10,height = 7)
lapply(pairLR$interaction_name, function(x){
  netVisual_individual(cellchat_GroupCTRL, signaling = pathways.show,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "hierarchy")
})
dev.off()
#
pdf("../../out/image/02_contribution_circleLR_SPP1_cellchat_GroupCTRLCellID.pdf",width = 10,height = 7)
lapply(pairLR$interaction_name, function(x){
  netVisual_individual(cellchat_GroupCTRL, signaling = pathways.show,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "circle")
})
dev.off()
#
pdf("../../out/image/02_contribution_chordLR_SPP1_cellchat_GroupCTRLCellID.pdf",width = 10,height = 7)
lapply(pairLR$interaction_name, function(x){
  netVisual_individual(cellchat_GroupCTRL, signaling = pathways.show,  pairLR.use = x, vertex.receiver = vertex.receiver,layout = "chord")
})
dev.off()

# Part IV: Systems analysis of cell-cell communication network ------------
# To facilitate the interpretation of the complex intercellular communication networks, CellChat quantitively measures networks through methods abstracted from graph theory, pattern recognition and manifold learning.

# Compute and visualize the network centrality scores
# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# cellchat_GroupCTRLCellID <- netAnalysis_computeCentrality(cellchat_GroupCTRLCellID, slot.name = "netP")

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pdf("../../out/image/04_contribution_role_SPP1_cellchat_GroupCTRLCellID.pdf",width = 5,height = 5)
netAnalysis_signalingRole_network(cellchat_GroupCTRL, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()

# Identify signals contributing most to outgoing or incoming signa --------
# We can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups.

# Identify global communication patterns to explore how multiple c --------
# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
# In addition to exploring detailed communications for individual pathways, an important question is how multiple cell groups and signaling pathways coordinate to function. CellChat employs a pattern recognition method to identify the global communication patterns.

# Identify and visualize outgoing communication pattern of secreti --------
# Outgoing patterns reveal how the sender cells (i.e. cells as signal source) coordinate with each other as well as how they coordinate with certain signaling pathways to drive communication.

# library(NMF)
# library(ggalluvial)
# Here we run selectK to infer the number of patterns.
selectK(cellchat_GroupCTRL, pattern = "outgoing")
ggsave("../../out/image/selectK_outgoing_cellchat_GroupCTRLCellID.pdf",width = 8,height = 4)
dev.off()

# Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 8.
nPatterns <- 5
pdf("../../out/image/patterns_outgoing_cellchat_GroupCTRLCellID.pdf",width = 15,height = 10)
cellchat_GroupCTRL <- identifyCommunicationPatterns(cellchat_GroupCTRL, pattern = "outgoing", k = nPatterns,height = 15,width = 15)
dev.off()
# river plot
netAnalysis_river(cellchat_GroupCTRL, pattern = "outgoing")
ggsave("../../out/image/river_outgoing_cellchat_GroupCTRLCellID.pdf",width = 8,height = 10)

# Identify and visualize incoming communication pattern of target  --------
# Incoming patterns show how the target cells (i.e. cells as signal receivers) coordinate with each other as well as how they coordinate with certain signaling pathways to respond to incoming signals.
selectK(cellchat_GroupCTRL, pattern = "incoming")
ggsave("../../out/image/selectK_incoming_cellchat_GroupCTRLCellID.pdf",width = 8,height = 4)
dev.off()
# siluette values begin to drop when the number of incoming patterns is 4.
nPatterns <- 5
pdf("../../out/image/patterns_incoming_cellchat_GroupCTRLCellID.pdf",width = 15,height = 10)
cellchat_GroupCTRL <- identifyCommunicationPatterns(cellchat_GroupCTRL, pattern = "incoming", k = nPatterns,height = 15,width = 15)
dev.off()
# river plot
netAnalysis_river(cellchat_GroupCTRL, pattern = "incoming")
ggsave("../../out/image/river_incoming_cellchat_GroupCTRLCellID.pdf",width = 8,height = 10)

# Manifold and classification learning analysis of signaling netwo --------
# Manifold and classification learning analysis of signaling networks

# Identify signaling groups based on structure similarity
cellchat_GroupCTRL <- computeNetSimilarity(cellchat_GroupCTRL, type = "functional")
cellchat_GroupCTRL <- netEmbedding(cellchat_GroupCTRL, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat_GroupCTRL <- netClustering(cellchat_GroupCTRL, type = "functional",do.parallel = F)
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat_GroupCTRL, type = "functional", label.size = 3.5)
ggsave("../../out/image/netVisual_functional_cellchat_GroupCTRLCellID.pdf",width = 5,height = 5)

#
cellchat_GroupCTRL <- computeNetSimilarity(cellchat_GroupCTRL, type = "structural")
cellchat_GroupCTRL <- netEmbedding(cellchat_GroupCTRL, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat_GroupCTRL <- netClustering(cellchat_GroupCTRL, type = "structural",do.parallel = F)
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat_GroupCTRL, type = "structural", label.size = 3.5)
ggsave("../../out/image/netVisual_structural_cellchat_GroupCTRLCellID.pdf",width = 5,height = 5)

# Part V: Save the CellChat object ----------------------------------------
saveRDS(cellchat_GroupCTRL, file = "../../out/object/cellchat_GroupCTRLCellID_full.rds")
