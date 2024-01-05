# libraries ---------------------------------------------------------------
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)

# available datasets ------------------------------------------------------
available_data <- AvailableData()

# run azimuth -------------------------------------------------------------
query <- readRDS("out/object/Azimuth_humancortexref_Tanaka_organoids_full_reanalyzed.rds")

# the metadata have all the info of interest
head(query@meta.data)

DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3) + NoLegend()
DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3,reduction = "ref.umap") + NoLegend()
DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3,split.by = "orig.ident") + NoLegend()

# save the UMAP coordinates of the new reference
query@reductions$ref.umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>% 
  write_tsv("out/table/Tanaka_CoordUMAP_AzimuthHumanCortexRef.tsv")

# save the new annotation from azimuth
query@meta.data %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>% 
  write_tsv("out/table/Tanaka_meta_AzimuthHumanCortexRef.tsv")

# costume plots threshold 0.75 --------------------------------------------
# costume visualization
UMAP_query <- read_tsv("out/table/Tanaka_organoids_full_reanalyzed_CoordUMAP.tsv")
meta_query <- read_tsv("out/table/Tanaka_meta_AzimuthHumanCortexRef.tsv") %>%
  mutate(robust_score_subclass = case_when(predicted.subclass.score>0.75&mapping.score>0.75~predicted.subclass,
                                           T~"uncertain"),
         robust_score_class = case_when(predicted.class.score>0.75&mapping.score>0.75~predicted.class,
                                        T~"uncertain"),
         robust_score_cluster = case_when(predicted.cluster.score>0.75&mapping.score>0.75~predicted.cluster,
                                          T~"uncertain"))

# define the levels of the cluster variable
level_subclass <- meta_query %>%
  group_by(predicted.subclass) %>%
  summarise(med = median(predicted.subclass.score)) %>%
  mutate(predicted.subclass = fct_reorder(predicted.subclass, med,.desc = T)) %>%
  pull(predicted.subclass) %>%
  levels()

# for each assignment what is the distribution of the scores
meta_query %>%
  mutate(predicted.subclass = factor(predicted.subclass, levels = level_subclass)) %>%
  ggplot(aes(x=predicted.subclass,y=predicted.subclass.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("out/image/Tanaka_human_motorcortexv1.0.0_subclass_score_075.pdf",height = 4,width = 4)

meta_query %>%
  mutate(predicted.subclass = factor(predicted.subclass, levels = level_subclass)) %>%
  ggplot(aes(y=predicted.subclass,x=predicted.subclass.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
ggsave("out/image/Tanaka_human_motorcortexv1.0.0_subclass_score_ridges_075.pdf",height = 5,width = 4)

# define the levels of the cluster variable
level_class <- meta_query %>%
  group_by(predicted.class) %>%
  summarise(med = median(predicted.class.score)) %>%
  mutate(predicted.class = fct_reorder(predicted.class, med,.desc = T)) %>%
  pull(predicted.class) %>%
  levels()

# plot the score varaible for the cluster
meta_query %>%
  mutate(predicted.class = factor(predicted.class, levels = level_class)) %>%
  ggplot(aes(x=predicted.class,y=predicted.class.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("out/image/Tanaka_human_motorcortexv1.0.0_class_score_075.pdf",height = 4,width = 2)

# define the levels of the cluster variable
level_cluster <- meta_query %>%
  group_by(predicted.cluster) %>%
  summarise(med = median(predicted.cluster.score)) %>%
  mutate(predicted.cluster = fct_reorder(predicted.cluster, med,.desc = T)) %>%
  pull(predicted.cluster) %>%
  levels()

# plot the score varaible for the cluster
meta_query %>%
  mutate(predicted.cluster = factor(predicted.cluster, levels = level_cluster)) %>%
  ggplot(aes(x=predicted.cluster,y=predicted.cluster.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45),
        plot.margin = margin(0.5, 0.5, 0.5, 2, "cm")) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave("out/image/Tanaka_human_motorcortexv1.0.0_cluster_score_075.pdf",height = 4,width = 6)

# identifyt he most likely assignment for each seurat cluster
# first using all the subcluster annotation, not filtering for threshold of scores
prop_table_subclass <- meta_query %>%
  group_by(seurat_clusters,predicted.subclass) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.subclass,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.subclass")

pdf("out/image/Tanaka_human_motorcortexv1.0.0_subclass_heatmapAll.pdf",height = 3,width = 6)
Heatmap(prop_table_subclass,
        name = "prop", 
        column_title = "subclass score all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

meta_query %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.subclass) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# first using all the subcluster annotation, filtering for threshold of scores
prop_table_subclass_filter <- meta_query %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.subclass) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.subclass,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.subclass") %>%
  as.matrix()

pdf("out/image/Tanaka_human_motorcortexv1.0.0_subclass_heatmapFilter_075.pdf",height = 3,width = 5)
Heatmap(prop_table_subclass_filter,
        name = "prop", 
        column_title = "subclass score high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

meta_query %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.subclass) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.subclass,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.subclass")

# add the meta to the coordinates
data_query <- left_join(UMAP_query,meta_query,by = c("barcodes"="barcode"))

# divide the dataset into uncertain and not
data_query_unc <- data_query %>%
  filter(robust_score_subclass == "uncertain")
#
data_query_defined <- data_query %>%
  filter(robust_score_subclass != "uncertain")

# average the position of the clusters
data_query_avg <- data_query_defined %>% group_by(robust_score_subclass) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)
data_query_avg2 <- data_query_defined %>% group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_subclass),col="black")+theme_bw()
ggsave("out/image/Tanaka_human_motorcortexv1.0.0_subclass_UMAP_075.pdf",width = 5,height = 3)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_subclass),col="black")+theme_bw()+facet_grid(~orig.ident)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("out/image/Tanaka_human_motorcortexv1.0.0_subclass_split_UMAP_075.pdf",width = 20,height = 3)

# costume visualization
UMAP_ref <- read_tsv("out/table/azimuth_humancortex_CoordUMAP.tsv")
UMAP_ref2 <- read_tsv("out/table/Tanaka_CoordUMAP_AzimuthHumanCortexRef.tsv")
meta_ref <- read_tsv("out/table/azimuth_humancortex_Metadata.tsv")
meta_query <- read_tsv("out/table/Tanaka_meta_AzimuthHumanCortexRef.tsv") %>%
  mutate(robust_score_subclass = case_when(predicted.subclass.score>0.75&mapping.score>0.75~predicted.subclass,
                                           T~"uncertain"),
         robust_score_class = case_when(predicted.class.score>0.75&mapping.score>0.75~predicted.class,
                                        T~"uncertain"),
         robust_score_cluster = case_when(predicted.cluster.score>0.75&mapping.score>0.75~predicted.cluster,
                                          T~"uncertain"))

# add the meta to the coordinates
data_ref <- left_join(UMAP_ref,meta_ref,"barcode")
data_ref2 <- left_join(UMAP_ref2,meta_query,"barcode") %>%
  mutate(seurat_clusters=factor(seurat_clusters))

# average the position of the clusters
data_ref_avg <- data_ref %>% group_by(subclass) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)

# average the position of the clusters
data_ref_avg2 <- data_ref2 %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+theme_bw()
ggsave("out/image/human_motorcortexv1.0.0_SeuratCluster_Tanaka_UMAP_075.pdf",width = 9,height = 6)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+theme_bw()+facet_wrap(~orig.ident)+theme(strip.background = element_blank())
ggsave("out/image/human_motorcortexv1.0.0_SeuratCluster_Tanaka_UMAP_075_split.pdf",width = 20,height = 18)

# average the position of the clusters
data_ref_avg2subclass <- data_ref2 %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(predicted.subclass) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.subclass),col="black")+theme_bw()
ggsave("out/image/human_motorcortexv1.0.0_subclass_Tanaka_UMAP_075.pdf",width = 8,height = 6)

# same as above no grid
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.subclass),col="black")+theme_cowplot()
ggsave("out/image/human_motorcortexv1.0.0_subclass_Tanaka_UMAP_075_alt.pdf",width = 8,height = 6)

# same as above, no text
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  theme_cowplot()
ggsave("out/image/human_motorcortexv1.0.0_subclass_Tanaka_UMAP_075_alt_noLabel.pdf",width = 8,height = 6)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  theme_cowplot()
ggsave("out/image/human_motorcortexv1.0.0_subclass_Tanaka_UMAP_075_alt_noLabel2.pdf",width = 6,height = 4)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.subclass),col="black")+theme_bw()+facet_wrap(~orig.ident)+theme(strip.background = element_blank())
ggsave("out/image/human_motorcortexv1.0.0_subclass_Tanaka_UMAP_075_split.pdf",width = 19,height = 18)

