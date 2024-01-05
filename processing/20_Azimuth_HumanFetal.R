# libraries ---------------------------------------------------------------
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(ComplexHeatmap)

# available datasets ------------------------------------------------------
available_data <- AvailableData()

# run azimuth -------------------------------------------------------------
query <- LoadFileInput(path = "out/object/diet_data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
head(query@meta.data)

# The RunAzimuth function can take a Seurat object as input
query <- RunAzimuth(query, reference = "fetusref")
saveRDS(query,"out/object/Azimuth_humanfetusref_harmonyMartina.rds")

# the metadata have all the info of interest
head(query@meta.data)

DimPlot(query, group.by = "predicted.annotation.l1", label = TRUE, label.size = 3) + NoLegend()
DimPlot(query, group.by = "predicted.annotation.l1", label = TRUE, label.size = 3,reduction = "ref.umap") + NoLegend()
DimPlot(query, group.by = "predicted.annotation.l1", label = TRUE, label.size = 3,split.by = "orig.ident") + NoLegend()

# save the UMAP coordinates of the new reference
query@reductions$ref.umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>% 
  write_tsv("out/table/data.combined_NOT_annotated_fix_normalized_CoordUMAP_AzimuthHumanFetalRef_resolution_DoubletSinglet_harmonyMartina.tsv")

# save the new annotation from azimuth
query@meta.data %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>% 
  write_tsv("out/table/data.combined_NOT_annotated_fix_normalized_meta_AzimuthHumanFetalRef_resolution_DoubletSinglet_harmonyMartina.tsv")

# costume plots threshold 0.75 --------------------------------------------
# costume visualization
UMAP_query <- read_tsv("out/table/data.combined_NOT_annotated_fix_normalized_CoordUMAP_resolution_DoubletSinglet_harmonyMartina.tsv")
meta_query <- read_tsv("out/table/data.combined_NOT_annotated_fix_normalized_meta_AzimuthHumanFetalRef_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(robust_score_organ = case_when(predicted.organ.score>0.75&mapping.score>0.75~predicted.organ,
                                           T~"uncertain"),
         robust_score_annol1 = case_when(predicted.annotation.l1.score>0.75&mapping.score>0.75~predicted.annotation.l1,
                                        T~"uncertain"),
         robust_score_annol2 = case_when(predicted.annotation.l2.score>0.75&mapping.score>0.75~predicted.annotation.l2,
                                          T~"uncertain"))

# define the levels of the cluster variable
level_organ <- meta_query %>%
  group_by(predicted.organ) %>%
  summarise(med = median(predicted.organ.score)) %>%
  mutate(predicted.organ = fct_reorder(predicted.organ, med,.desc = T)) %>%
  pull(predicted.organ) %>%
  levels()

# for each assignment what is the distribution of the scores
meta_query %>%
  mutate(predicted.organ = factor(predicted.organ, levels = level_organ)) %>%
  ggplot(aes(x=predicted.organ,y=predicted.organ.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_organ_score_075_harmonyMartina.pdf",height = 4,width = 4)

meta_query %>%
  mutate(predicted.organ = factor(predicted.organ, levels = level_organ)) %>%
  ggplot(aes(y=predicted.organ,x=predicted.organ.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_organ_score_ridges_075_harmonyMartina.pdf",height = 5,width = 4)

# -------------------------------------------------------------------------
# define the levels of the cluster variable
level_annol1 <- meta_query %>%
  group_by(predicted.annotation.l1) %>%
  summarise(med = median(predicted.annotation.l1.score)) %>%
  mutate(predicted.annotation.l1 = fct_reorder(predicted.annotation.l1, med,.desc = T)) %>%
  pull(predicted.annotation.l1) %>%
  levels()

# for each assignment what is the distribution of the scores
meta_query %>%
  mutate(predicted.annotation.l1 = factor(predicted.annotation.l1, levels = level_annol1)) %>%
  ggplot(aes(x=predicted.annotation.l1,y=predicted.annotation.l1.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol1_score_075_harmonyMartina.pdf",height = 4,width = 16)

meta_query %>%
  mutate(predicted.annotation.l1 = factor(predicted.annotation.l1, levels = level_annol1)) %>%
  ggplot(aes(y=predicted.annotation.l1,x=predicted.annotation.l1.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol1_score_ridges_075_harmonyMartina.pdf",height = 8,width = 5)

# define the anno l2
level_annol2 <- meta_query %>%
  group_by(predicted.annotation.l2) %>%
  summarise(med = median(predicted.annotation.l2.score)) %>%
  mutate(predicted.annotation.l2 = fct_reorder(predicted.annotation.l2, med,.desc = T)) %>%
  pull(predicted.annotation.l2) %>%
  levels()

# for each assignment what is the distribution of the scores
meta_query %>%
  mutate(predicted.annotation.l2 = factor(predicted.annotation.l2, levels = level_annol2)) %>%
  ggplot(aes(x=predicted.annotation.l2,y=predicted.annotation.l2.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol2_score_075_harmonyMartina.pdf",height = 4,width = 25)

meta_query %>%
  mutate(predicted.annotation.l2 = factor(predicted.annotation.l2, levels = level_annol2)) %>%
  ggplot(aes(y=predicted.annotation.l2,x=predicted.annotation.l2.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol2_score_ridges_075_harmonyMartina.pdf",height = 16,width = 5)

# identifyt he most likely assignment for each seurat cluster
# first using all the subcluster annotation, not filtering for threshold of scores

# level annotation organ
prop_table_organ <- meta_query %>%
  group_by(seurat_clusters,predicted.organ) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.organ,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.organ") %>% 
  as.matrix()

pdf("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_organ_heatmapAll_harmonyMartina.pdf",height = 6,width = 12)
Heatmap(prop_table_organ,
        name = "prop", 
        column_title = "organ score all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# first using all the subcluster annotation, filtering for threshold of scores
prop_table_organ_filter <- meta_query %>%
  # filter(robust_score_subclass != "uncertain") %>%
  filter(robust_score_organ != "uncertain") %>%
  group_by(seurat_clusters,predicted.organ) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.organ,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.organ") %>%
  as.matrix()

pdf("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_organ_heatmapFilter_075_harmonyMartina.pdf",height = 4,width = 12)
Heatmap(prop_table_organ_filter,
        name = "prop", 
        column_title = "organ score high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# level annotation l1
prop_table_annol1 <- meta_query %>%
  group_by(seurat_clusters,predicted.annotation.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.annotation.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.annotation.l1") %>% 
  as.matrix()

pdf("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol1_heatmapAll_harmonyMartina.pdf",height = 12,width = 12)
Heatmap(prop_table_annol1,
        name = "prop", 
        column_title = "annol1 score all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# first using all the subcluster annotation, filtering for threshold of scores
prop_table_annol1_filter <- meta_query %>%
  # filter(robust_score_subclass != "uncertain") %>%
  filter(robust_score_annol1 != "uncertain") %>%
  group_by(seurat_clusters,predicted.annotation.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.annotation.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.annotation.l1") %>%
  as.matrix()

pdf("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol1_heatmapFilter_075_harmonyMartina.pdf",height = 6,width = 12)
Heatmap(prop_table_annol1_filter,
        name = "prop", 
        column_title = "annol1 score high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# level annotation l2
prop_table_annol2 <- meta_query %>%
  group_by(seurat_clusters,predicted.annotation.l2) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.annotation.l2,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.annotation.l2") %>% 
  as.matrix()

pdf("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol2_heatmapAll_harmonyMartina.pdf",height = 20,width = 12)
Heatmap(prop_table_annol2,
        name = "prop", 
        column_title = "annol2 score all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# first using all the subcluster annotation, filtering for threshold of scores
prop_table_annol2_filter <- meta_query %>%
  # filter(robust_score_subclass != "uncertain") %>%
  filter(robust_score_annol2 != "uncertain") %>%
  group_by(seurat_clusters,predicted.annotation.l2) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.annotation.l2,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.annotation.l2") %>%
  as.matrix()

pdf("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol2_heatmapFilter_075_harmonyMartina.pdf",height = 6,width = 12)
Heatmap(prop_table_annol2_filter,
        name = "prop", 
        column_title = "annol2 score high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# add the meta to the coordinates
data_query <- left_join(UMAP_query,meta_query,"barcode")

# divide the dataset into uncertain and not
data_query_unc <- data_query %>%
  filter(robust_score_annol1 == "uncertain")
#
data_query_defined <- data_query %>%
  filter(robust_score_annol1 != "uncertain")

# average the position of the clusters
data_query_avg <- data_query_defined %>% group_by(robust_score_annol1) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)
data_query_avg2 <- data_query_defined %>% group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_annol1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_annol1),col="black")+theme_bw()
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol1_UMAP_075_harmonyMartina.pdf",width = 5,height = 3)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_annol1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_annol1),col="black")+theme_bw()+facet_grid(~orig.ident)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol1_split_UMAP_075_harmonyMartina.pdf",width = 20,height = 3)

ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_annol1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_annol1),col="black")+theme_bw()+facet_grid(~group_martina)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_fetalv1.0.0_annol1_split_UMAP_075_queryMartina_harmonyMartina.pdf",width = 8,height = 3)

# costume visualization
UMAP_ref <- read_tsv("out/table/azimuth_fetusref_CoordUMAP.tsv")
UMAP_ref2 <- read_tsv("out/table/data.combined_NOT_annotated_fix_normalized_CoordUMAP_AzimuthHumanFetalRef_resolution_DoubletSinglet_harmonyMartina.tsv")
meta_ref <- read_tsv("out/table/azimuth_fetusref_Metadata.tsv")
meta_query <- read_tsv("out/table/data.combined_NOT_annotated_fix_normalized_meta_AzimuthHumanFetalRef_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(robust_score_organ = case_when(predicted.organ.score>0.75&mapping.score>0.75~predicted.organ,
                                        T~"uncertain"),
         robust_score_annol1 = case_when(predicted.annotation.l1.score>0.75&mapping.score>0.75~predicted.annotation.l1,
                                         T~"uncertain"),
         robust_score_annol2 = case_when(predicted.annotation.l2.score>0.75&mapping.score>0.75~predicted.annotation.l2,
                                         T~"uncertain"))

# add the meta to the coordinates
data_ref <- left_join(UMAP_ref,meta_ref,"barcode")
data_ref2 <- left_join(UMAP_ref2,meta_query,"barcode") %>%
  mutate(seurat_clusters=factor(seurat_clusters))

# average the position of the clusters
data_ref_avg <- data_ref %>% group_by(annotation.l1) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)

# average the position of the clusters
data_ref_avg2 <- data_ref2 %>%
  filter(robust_score_annol1 != "uncertain") %>%
  group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_annol1 != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black",max.overlaps = 12)+theme_bw()
ggsave("out/image/human_fetalv1.0.0_annol1_data.combined_annotated_norm_fix_resolution_UMAP_075_harmonyMartina.pdf",width = 9,height = 6)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_annol1 != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black",max.overlaps = 15)+theme_bw()+facet_wrap(~orig.ident)+theme(strip.background = element_blank())
ggsave("out/image/human_fetalv1.0.0_annol1_data.combined_annotated_norm_fix_resolution_UMAP_075_split_harmonyMartina.pdf",width = 20,height = 22)

# average the position of the clusters
data_ref_avg2annol1 <- data_ref2 %>%
  filter(robust_score_annol1 != "uncertain") %>%
  group_by(predicted.annotation.l1) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_annol1 != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.annotation.l1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2annol1,aes(x = UMAP_1,y = UMAP_2,label = predicted.annotation.l1),col="black")+theme_bw()
ggsave("out/image/human_fetalv1.0.0_annol1_data.combined_annotated_norm_fix_resolution_UMAP_075_annol1_harmonyMartina.pdf",width = 8,height = 6)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_annol1 != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.annotation.l1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2annol1,aes(x = UMAP_1,y = UMAP_2,label = predicted.annotation.l1),col="black")+theme_bw()+facet_wrap(~orig.ident)+theme(strip.background = element_blank())
ggsave("out/image/human_fetalv1.0.0_annol1_data.combined_annotated_norm_fix_resolution_UMAP_075_annol1_split_harmonyMartina.pdf",width = 20,height = 21)
