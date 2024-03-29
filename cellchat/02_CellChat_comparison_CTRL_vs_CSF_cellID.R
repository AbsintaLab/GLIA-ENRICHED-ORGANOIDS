# This vignette shows how to apply CellChat to identify major signaling changes as well as conserved and context-specific signaling by joint manifCSF learning and quantitative contrasts of multiple cell-cell communication networks. We showcase CellChat’s diverse functionalities by applying it to a scRNA-seq data on cells from two biological conditions: nonlesional (NL, normal) and lesional (LS, diseased) human skin from patients with atopic dermatitis. These two datasets (conditions) have the same cell population compositions after joint clustering. If there are slightly or vastly different cell population compositions between different datasets, please check out another related tutorial (Comparison analysis of multiple datasets with different cell type compositions).

# CellChat employs a top-down approach, i.e., starting with the big picture and then refining it in a greater detail on the signaling mechanisms, to identify signaling changes at different levels, including both general principles of cell-cell communication and dysfunctional cell populations/signaling pathways/ligand-receptors.

# Load the required libraries ---------------------------------------------
library(CellChat)
library(patchwork)
library(circlize)

# Create a directory to save figures --------------------------------------
# Load CellChat object of each dataset and then merge together ------------
# USERS need to run CellChat on each dataset seperately and then merge different CellChat objects together. Please do updateCellChat if you have CellChat objects that are obtained using the earlier version (< 0.5.0).

cellchat.CTRL <- readRDS("../../out/object/cellchat_GroupCTRLCellID_full.rds")
cellchat.CSF <- readRDS("../../out/object/cellchat_GroupCSFCellID_full.rds")

# explore the entity of the rds objects
cellchat.CTRL
cellchat.CSF
# these are cellchat objects that can be generated usign:

# put the objects in a list
object.list <- list(CTRL = cellchat.CTRL, CSF = cellchat.CSF)

# merge the objects
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# print the content of the object
cellchat

# Part I: Predict general principles of cell-cell communication -----------
# CellChat starts with the big picture to predict general principles of cell-cell communication. When comparing cell-cell communication among multiple biological conditions, it can answer the following biological questions:
# Whether the cell-cell communication is enhanced or not
# The interaction between which cell types is significantly changed
# How the major sources and targets change from one condition to another

# Compare the total number of interactions and interaction strength -------
# To answer on question on whether the cell-cell communication is enhanced or not, CellChat compares the the total number of interactions and interaction strength of the inferred cell-cell communication networks from different biological conditions.

# measure: "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

gg1 + gg2
ggsave("../../out/image/comparison_01_barplot_interaction_cellID.pdf",width = 10,height = 5)

# Compare the number of interactions and interaction strength amon --------
# Compare the number of interactions and interaction strength among different cell populations
# To identify the interaction between which cell populations showing significant changes, CellChat compares the number of interactions and interaction strength among different cell populations.

# Differential number of interactions or interaction strength amon --------
# Differential number of interactions or interaction strength among different cell populations
# The differential number of interactions or interaction strength in the cell-cell communication network between two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.

pdf("../../out/image/comparison_02_differential_interaction_cellID.pdf",width = 20,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
list(netVisual_diffInteraction_fix(cellchat, weight.scale = T),
     netVisual_diffInteraction_fix(cellchat, weight.scale = T, measure = "weight"))
dev.off()

pdf("../../out/image/comparison_02_differential_interaction_UP_cellID.pdf",width = 20,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
list(fun_test(cellchat, weight.scale = T,filter_edges = "up"),
     fun_test(cellchat, weight.scale = T,filter_edges = "up", measure = "weight"))
dev.off()

pdf("../../out/image/comparison_02_differential_interaction_DOWN_cellID.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
list(fun_test(cellchat, weight.scale = T,filter_edges = "down"),
     fun_test(cellchat, weight.scale = T,filter_edges = "down", measure = "weight"))
dev.off()

# We can also show differential number of interactions or interaction strength in a greater details using a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.

gg1 <- netVisual_heatmap(cellchat,height = 5,comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",height = 5)

pdf("../../out/image/comparison_02_differential_interaction_heatmap_cellID.pdf",width = 10,height = 5)
gg1 + gg2
dev.off()

# To better control the node size and edge weights of the inferred networks across different datasets, we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
weight.max

pdf("../../out/image/comparison_03_interaction_individual_number_cellID.pdf",width = 20,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
lapply(1:length(object.list),function(i){
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale = T,
                   label.edge= F,
                   edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
})
dev.off()

# same as above but for strength
pdf("../../out/image/comparison_03_interaction_individual_strenght_cellID.pdf",width = 20,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
lapply(1:length(object.list),function(i){
  netVisual_circle(object.list[[i]]@net$weight,
                   weight.scale = T,
                   label.edge= F,
                   # edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = paste0("Strenght of interactions - ", names(object.list)[i]))
})
dev.off()

# Compare the major sources and targets in 2D space -----------------------
# Comparing the outgoing and incoming interaction strength in 2D space allows ready identification of the cell populations with significant changes in sending or receiving signals between different datasets.

num.link <- sapply(object.list, function(x){
  rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)
})
num.link

# control the dot size in the different datasets
weight.MinMax <- c(min(num.link), max(num.link))
weight.MinMax

gg <- lapply(1:length(object.list),function(i){
  netAnalysis_signalingRole_scatter(object.list[[i]],
                                    title = names(object.list)[i],
                                    weight.MinMax = weight.MinMax)
})

pdf("../../out/image/comparison_04_compare_source_target_cellID.pdf",width = 10,height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

# From the scatter plot, we can see that Inflam.DC and cDC1 emerge as one of the major source and targets in LS compared to NL. Fibroblast populations also become the major sources in LS.

# Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS. ## Identify signaling changes associated with one cell group

# Visualizing differential outgoing and incoming signaling changes from NL to LS
id <- c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC")

list_gg <- lapply(id, function(x){
  netAnalysis_signalingChanges_scatter(cellchat, idents.use = x)
})

pdf("../../out/image/comparison_05_compare_signalling_cellID.pdf",width = 30,height = 15)
patchwork::wrap_plots(plots = list_gg)
dev.off()

# Part II: Identify the conserved and context-specific signaling p --------
# CellChat then can identify signaling networks with larger (or less) difference, signaling groups, and the conserved and context-specific signaling pathways based on their cell-cell communication networks among multiple biological conditions.

# Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
# CellChat performs joint manifCSF learning and classification of the inferred communication networks based on their functional and topological similarity. NB: Such analysis is applicable to more than two datasets.

# Functional similarity: High degree of functional similarity indicates major senders and receivers are similar, and it can be interpreted as the two signaling pathways or two ligand-receptor pairs exhibit similar and/or redundant roles. NB: Functional similarity analysis is not applicable to multiple datsets with different cell type composition.

# Structural similarity: A structural similarity was used to compare their signaling network structure, without considering the similarity of senders and receivers. NB: Structural similarity analysis is applicable to multiple datsets with the same cell type composition or the vastly different cell type composition.

#  Here we can run the manifCSF and classification learning analysis based on the functional similarity because the two datasets have the the same cell type composition.

# Identify signaling groups based on their functional similarity

# Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural",do.parallel = F)

# Visualization in 2D-space
pdf("../../out/image/comparison_05_structural_similarity_cellID.pdf",width = 10,height = 10)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()

pdf("../../out/image/comparison_05_structural_similarity_grid_cellID.pdf",width = 10,height = 10)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

# Compute and visualize the pathway distance in the learned joint  --------

# Compute and visualize the pathway distance in the learned joint manifCSF
# We can identify the signaling networks with larger (or less) difference based on their Euclidean distance in the shared two-dimensions space. Larger distance implies larger difference of the communication networks between two datasets in terms of either functional or structure similarity. NB: We only compute the distance of overlapped signaling pathways between two datasets. Those signaling pathways that are only identified in one dataset are not considered here. If there are more than three datasets, one can do pairwise comparisons by defining comparison in the function rankSimilarity.

pdf("../../out/image/comparison_06_rank_similarity_cellID.pdf",width = 5,height = 5)
rankSimilarity(cellchat, type = "structural")
dev.off()

# Identify and visualize the conserved and context-specific signal --------

# Identify and visualize the conserved and context-specific signaling pathways
# By comparing the information flow/interaction strengh of each signaling pathway, we can identify signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase, by change their information flow at one condition as compared to another condition.

# Compare the overall information flow of each signaling pathway
# We can identify the conserved and context-specific signaling pathways by simply comparing the information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).

# This bar graph can be plotted in a stacked mode or not. Significant signaling pathways were ranked based on differences in the overall information flow within the inferred networks between NL and LS skin. The top signaling pathways colored red are enriched in NL skin, and these colored green were enriched in the LS skin.

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

pdf("../../out/image/comparison_07_signaling_pathway_comparison_cellID.pdf",width = 10,height = 8)
gg1 + gg2
dev.off()

# Compare outgoing (or incoming) signaling associated with each ce --------

# Compare outgoing (or incoming) signaling associated with each cell population
# The above analysis summarize the information from the outgoing and incoming signaling together. We can also compare the outgoing (or incoming) signaling pattern between two datasets, allowing to identify signaling pathways/ligand-receptors that exhibit different signaling patterns.

# We can combine all the identified signaling pathways from different datasets and thus compare them side by side, including outgoing signaling, incoming signaling and overall signaling by aggregating outgoing and incoming signaling together. NB: rankNet also shows the comparison of overall signaling, but it does not show the signaling strength in specific cell populations.

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pathway.union

names(object.list[1])
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                         pattern = "outgoing",
                                         color.heatmap = "GnBu",
                                         signaling = pathway.union,
                                         title = names(object.list)[i], width = 6, height = 18)

ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                         pattern = "outgoing",
                                         color.heatmap = "GnBu",
                                         signaling = pathway.union,
                                         title = names(object.list)[i+1], width = 6, height = 18)

pdf("../../out/image/comparison_08_signaling_pathway_comparison_outgoing_cellID.pdf",width = 10,height = 10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()

#
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                         pattern = "incoming",
                                         signaling = pathway.union,
                                         color.heatmap = "GnBu",
                                         title = names(object.list)[i], width = 6, height = 18)

ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                         pattern = "incoming",
                                         color.heatmap = "GnBu",
                                         signaling = pathway.union,
                                         title = names(object.list)[i+1], width = 6, height = 18)

pdf("../../out/image/comparison_08_signaling_pathway_comparison_incoming_cellID.pdf",width = 10,height = 10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()

# overall
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
                                         pattern = "all",
                                         signaling = pathway.union,
                                         color.heatmap = "GnBu",
                                         title = names(object.list)[i], width = 6, height = 18)

ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                         pattern = "all",
                                         color.heatmap = "GnBu",
                                         signaling = pathway.union,
                                         title = names(object.list)[i+1], width = 6, height = 18)

pdf("../../out/image/comparison_08_signaling_pathway_comparison_overall_cellID.pdf",width = 10,height = 10)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()

# Part III: Identify the upgulated and down-regulated signaling li --------

# Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
# Identify dysfunctional signaling by comparing the communication probabities
# We can compare the communication probabilities mediated by ligand-receptor pairs from some cell groups to other cell groups. This can be done by setting comparison in the function netVisual_bubble.

# check the position of the levels
levels(cellchat@idents$joint)

gg1 <- netVisual_bubble(cellchat,
                        sources.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
                        targets.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
                        comparison = c(1, 2),
                        max.dataset = 2,
                        title.name = "Increased signaling in CSF", angle.x = 45, remove.isolate = T)

gg2 <- netVisual_bubble(cellchat,
                        sources.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
                        targets.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
                        comparison = c(1, 2),
                        max.dataset = 1,
                        title.name = "Decreased signaling in CSF", angle.x = 45, remove.isolate = T)

pdf("../../out/image/comparison_09_L_R_pair_significant_ALL_cellID.pdf",width = 30,height = 15)
gg1 + gg2
dev.off()
#
# NB: The ligand-receptor pairs shown in the bubble plot can be accessed via signaling.LSIncreased = gg1$data.
gg1$data

# Identify dysfunctional signaling by using differential expressio --------

# Identify dysfunctional signaling by using differential expression analysis
# The above method for identifying the upgulated and down-regulated signaling is perfomed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups. Alternative, we can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis. Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, and then obtain the upgulated and down-regulated signaling based on the fCSF change of ligands in the sender cells and receptors in the receiver cells. Such analysis can be done as follows.

# define a positive dataset, i.e., the dataset with positive fCSF change against the other dataset
pos.dataset = "CSF"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat,
                                       group.dataset = "datasets",
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

# net %>% filter(ligand == "Vegfa",receptor == "Kdr")

# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "CSF",ligand.logFC = 0.2, receptor.logFC = NULL)
net.up %>% 
  write_tsv("../../out/table/net.up_cellID.tsv")

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = -0.1, receptor.logFC = -0.1)

# Since the signaling genes in the net.up and net.down might be complex with multi-subunits, we can do further deconvolution to obtain the individual signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.up,
                        sources.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
                        targets.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.down,
                        sources.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
                        targets.use = c("ASTRO","CYCLING","GLIA","iMIC","NEU","OLIGO","OPC"),
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

pdf("../../out/image/comparison_09_L_R_pair_significant_DEG_cellID.pdf",width = 35,height = 8)
gg1 + gg2
dev.off()

# Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram

# Chord diagram
levels(cellchat@idents$joint)

# -------------------------------------------------------------------------
# whole dataset
# plot all the pathway that are goind up in CCA
pdf("../../out/image/comparison_09_pathway_significant_CSF_chord_cellID.pdf",width = 15,height = 15)
list(
  netVisual_chord_gene2(object.list[[2]],
                        slot.name = 'netP',
                        net = net.up, lab.cex = 0.8, small.gap = 3.5,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2])),
  netVisual_chord_gene2(object.list[[1]],
                        slot.name = 'netP',
                        net = net.down, lab.cex = 0.8, small.gap = 3.5,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
)
dev.off()

# check the position of the levels
levels(cellchat@idents$joint)

# plot all the LR
pdf("../../out/image/comparison_09_LR_significant_CSF_chord_cellID.pdf",width = 15,height = 15)
list(
  netVisual_chord_gene2(object.list[[2]],
                        slot.name = 'net',
                        net = net.up, lab.cex = 0.8, small.gap = 3.5,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2])),
  netVisual_chord_gene2(object.list[[1]],
                        slot.name = 'net',
                        net = net.down, lab.cex = 0.8, small.gap = 3.5,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
)
dev.off()

# Save the merged CellChat object -----------------------------------------
saveRDS(cellchat, file = "../../out/object/cellchat_comparison_CTRL_vs_CSF_cellID.rds")
