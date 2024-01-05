# libraries ---------------------------------------------------------------
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(ggrepel)

# read in the data --------------------------------------------------------
data.combined <- readRDS(file = "data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
Idents(data.combined)<-"seurat_clusters"

# save the table of normalized expression
write_tsv(data.combined@assays$RNA@data %>%
            as.matrix() %>%
            data.frame() %>%
            rownames_to_column(var = "gene"),"out/table/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina_normCounts.tsv")

# plots -------------------------------------------------------------------
# main umap
plot03 <- DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters",label = T)
ggsave(plot=plot03,"out/image/UMAP_harmonyMartina_cluster_norm_fix_label_resolution_DoubletSinglet.pdf",width = 5,height = 4)

# main umap cell cycle
plot04 <- DimPlot(data.combined, reduction = "umap", group.by = "Phase")
ggsave(plot=plot04,"out/image/UMAP_harmonyMartina_cluster_norm_fix_phase_resolution_DoubletSinglet.pdf",width = 5,height = 4)

# split by sample
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcode")
df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

data2 <- left_join(df_UMAP,df_meta,"barcode")
data2_avg <- data2 %>% group_by(seurat_clusters) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# single plot
plot05 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(plot=plot05,"out/image/UMAP_harmonyMartina_cluster_norm_fix_label_resolution_DoubletSinglet_ggplot.pdf",width = 5,height = 3)

# save the same with theme cowplot
plot05_alt <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot()

ggsave(plot=plot05_alt,"out/image/UMAP_harmonyMartina_cluster_norm_fix_label_resolution_DoubletSinglet_ggplot_noGrid.pdf",width = 5,height = 3)

# single plot
plot051 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.ribo) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=percent.ribo),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  scale_color_viridis_c(option = "turbo")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(plot = plot051,"out/image/UMAP_harmonyMartina_cluster_norm_fix_label_resolution_DoubletSinglet_ggplot_RIBO.pdf",width = 5,height = 3)

# single plot
plot052 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.mt) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=percent.mt),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  scale_color_viridis_c(option = "turbo")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(plot = plot052,"out/image/UMAP_harmonyMartina_cluster_norm_fix_label_resolution_DoubletSinglet_ggplot_MITO.pdf",width = 5,height = 3)

# single plot
plot053 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.mt) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=nCount_RNA),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  scale_color_viridis_c(option = "turbo")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(plot = plot053,"out/image/UMAP_harmonyMartina_cluster_norm_fix_label_resolution_DoubletSinglet_ggplot_nCount.pdf",width = 5,height = 3)

# split by sample
plot06 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  facet_grid(~orig.ident)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(plot=plot06,"out/image/UMAP_harmonyMartina_cluster_norm_fix_label_split_resolution_DoubletSinglet.pdf",width = 20,height = 3)

# split martina
data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  facet_grid(~group_martina)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/UMAP_harmonyMartina_cluster_norm_fix_label_splitMartina_resolution_DoubletSinglet.pdf",width = 8,height = 3)

# split martina no grid
data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  facet_grid(~group_martina)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot()+
  theme(strip.background = element_blank())
ggsave("out/image/UMAP_harmonyMartina_cluster_norm_fix_label_splitMartina_resolution_DoubletSinglet_alt.pdf",width = 8,height = 3)

data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  filter(sample == "RR25_SOX10_CSF") %>% 
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  facet_grid(~clone)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot()+
  theme(strip.background = element_blank())
ggsave("out/image/UMAP_RR25.pdf",width = 4,height = 3)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(orig.ident,seurat_clusters) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary,"out/table/summary_data.combined_NOT_annotated_norm_fix_resolution_DoubletSinglet_harmonyMartina.tsv")

# define a convenient palette of colors
col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
show_col(col_pal)

plot07<-df_summary %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  ggplot(aes(x=seurat_clusters,y=prop,fill=orig.ident))+geom_col()+theme_bw()+scale_fill_manual(values = col_pal)
ggsave(plot=plot07,"out/image/summary_stacked_data.combined_NOT_annotated_norm_fix_resolution_DoubletSinglet_harmonyMartina.pdf",width = 5,height = 4)

plot08<-df_summary %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  ggplot(aes(x=seurat_clusters,y=prop,fill=orig.ident))+geom_col(position = "dodge")+theme_bw()+scale_fill_manual(values = col_pal)
ggsave(plot=plot08,"out/image//summary_dodge_data.combined_NOT_annotated_norm_fix_resolution_DoubletSinglet_harmonyMartina.pdf",width = 8,height = 4)

shortlist_features_list <- list(
  IMMUNE = c("AIF1","TYROBP","FTL","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
  # Neu = c("SYT1")
  NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
)

# martina asked to remove FTL
shortlist_features_list2 <- list(
  IMMUNE = c("AIF1","TYROBP","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
  # Neu = c("SYT1")
  NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
)

shortlist_features_list_long <- list(
  IMMUNE = c("IGHG1","CD38","CD8A","CD2","SKAP1","LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"))

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_short <- DotPlot(data.combined, features = shortlist_features_list, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot = test_short,"out/image/Dotplot_clusters_NOT_annotated_short_panle_integration_norm_fix_long_resolution_DoubletSinglet_harmonyMartina.pdf",width = 13,height = 5)

test_short2 <- DotPlot(data.combined, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot = test_short2,"out/image/Dotplot_clusters_NOT_annotated_short_panle_integration_norm_fix_long_resolution_DoubletSinglet_harmonyMartina2.pdf",width = 13,height = 5)

test_long <- DotPlot(data.combined, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()

ggsave(plot=test_long,"out/image/Dotplot_clusters_NOT_annotated_big_panle_integration_norm_fix_long_resolution_DoubletSinglet_harmonyMartina.pdf",width = 25,height = 5)

df_meta %>%
  group_by(orig.ident) %>%
  summarise(n=n())


# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(data.combined) <- "RNA"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("out/table/FindAllMarkers_harmonyMartina_data.combined_fix_DoubletSinglet.tsv")

