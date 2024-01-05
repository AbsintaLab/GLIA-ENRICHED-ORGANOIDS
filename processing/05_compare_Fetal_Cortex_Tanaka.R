# libraries ---------------------------------------------------------------
library(tidyverse)

# read in the data --------------------------------------------------------
UMAP_ref_fetal <- read_tsv("out/table/azimuth_humanfetus_CoordUMAP.tsv")
UMAP_ref2_fetal <- read_tsv("out/table/Tanaka_CoordUMAP_AzimuthHumanFetusRef.tsv")

meta_ref_fetal <- read_tsv("out/table/azimuth_humanfetus_Metadata.tsv")
meta_query_fetal <- read_tsv("out/table/Tanaka_meta_AzimuthHumanFetusRef.tsv") %>%
  mutate(robust_score_organ = case_when(predicted.organ.score>0.75&mapping.score>0.75~predicted.organ,
                                        T~"uncertain"),
         robust_score_annol1 = case_when(predicted.annotation.l1.score>0.75&mapping.score>0.75~predicted.annotation.l1,
                                         T~"uncertain"),
         robust_score_annol2 = case_when(predicted.annotation.l2.score>0.75&mapping.score>0.75~predicted.annotation.l2,
                                         T~"uncertain"))

data_ref_fetal <- left_join(UMAP_ref_fetal,meta_ref_fetal,"barcode")
data_ref2_fetal <- left_join(UMAP_ref2_fetal,meta_query_fetal,"barcode") %>%
  mutate(seurat_clusters=factor(seurat_clusters))

# average the position of the clusters
data_ref_avg_fetal <- data_ref_fetal %>% group_by(annotation.l1) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)

# average the position of the clusters
data_ref_avg2_fetal <- data_ref2_fetal %>%
  filter(robust_score_annol1 != "uncertain") %>%
  group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

meta_ref_fetal %>% 
  pull(annotation.l1) %>%
  table()

meta_query_fetal %>%
  filter(predicted.annotation.l1 %in% c("Astrocytes","ENS neurons","Excitatory neurons","Granule neurons","Inhibitory interneurons","Inhibitory neurons","Limbic system neurons","Myeloid cells","Oligodendrocytes","Purkinje neurons","Visceral neurons")) %>% 
  ggplot(aes(y=predicted.annotation.l1,x=predicted.annotation.l1.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")

# do the same for the Cortex ----------------------------------------------
# costume visualization
UMAP_ref_cortex <- read_tsv("out/table/azimuth_humancortex_CoordUMAP.tsv")
UMAP_ref2_cortex <- read_tsv("out/table/Tanaka_CoordUMAP_AzimuthHumanCortexRef.tsv")
meta_ref_cortex <- read_tsv("out/table/azimuth_humancortex_Metadata.tsv")
meta_query_cortex <- read_tsv("out/table/Tanaka_meta_AzimuthHumanCortexRef.tsv") %>%
  mutate(robust_score_subclass = case_when(predicted.subclass.score>0.50&mapping.score>0.50~predicted.subclass,
                                           T~"uncertain"),
         robust_score_class = case_when(predicted.class.score>0.50&mapping.score>0.50~predicted.class,
                                        T~"uncertain"),
         robust_score_cluster = case_when(predicted.cluster.score>0.50&mapping.score>0.50~predicted.cluster,
                                          T~"uncertain"))

# add the meta to the coordinates
data_ref_cortex <- left_join(UMAP_ref_cortex,meta_ref_cortex,"barcode")
data_ref2_cortex <- left_join(UMAP_ref2_cortex,meta_query_cortex,"barcode") %>%
  mutate(seurat_clusters=factor(seurat_clusters))

# average the position of the clusters
data_ref_avg_cortex <- data_ref_cortex %>% group_by(subclass) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)

# average the position of the clusters
data_ref_avg2_cortex <- data_ref2_cortex %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

meta_ref_cortex %>% 
  pull(subclass) %>%
  table()

meta_query_cortex %>% 
  ggplot(aes(x=predicted.subclass,y=predicted.subclass.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")

meta_query_cortex %>%
  filter(predicted.subclass %in% c("Astro","Oligo","Micro-PVM")) %>% 
  ggplot(aes(y=predicted.subclass,x=predicted.subclass.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")

# harmonyze the annotations -----------------------------------------------
# try to harmonyze that annotations and put the two datasets together

df_cortex <- meta_query_cortex %>%
  filter(predicted.subclass %in% c("Astro","Oligo","Micro-PVM","OPC","L2/3 IT","L5 ET","L5 IT","L5/6 NP","L6 CT","L6 IT","L6 IT Car3","L6b","Lamp5","Pvalb","Sncg","Sst","Sst Chodl","Vip")) %>% 
  dplyr::select(annotation = predicted.subclass,score = predicted.subclass.score,seurat_clusters) %>% 
  mutate(ref = "adult") %>% 
  # harmonyze the annotation
  mutate(annotation_harmony = case_when(annotation == "Astro"~"ASTRO",
                                        annotation == "Oligo"~"OLIGO",
                                        annotation == "Micro-PVM"~"MG",
                                        annotation == "OPC"~"OPC",
                                        T~"NEU"))

df_fetal <- meta_query_fetal %>%
  filter(predicted.annotation.l1 %in% c("Astrocytes","Microglia","Oligodendrocytes","Inhibitory interneurons","Inhibitory neurons")) %>% 
  dplyr::select(annotation = predicted.annotation.l1,score = predicted.annotation.l1.score,seurat_clusters) %>% 
  mutate(ref = "fetal") %>% 
  # harmonyze the annotation
  mutate(annotation_harmony = case_when(annotation == "Astrocytes"~"ASTRO",
                                        annotation == "Oligodendrocytes"~"OLIGO",
                                        annotation == "Microglia"~"MG",
                                        T~"NEU"))

bind_rows(list(df_cortex,df_fetal)) %>% 
  ggplot(aes(y=annotation_harmony,x=score,fill=ref))+
  ggridges::geom_density_ridges(alpha=0.6)+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
ggsave("out/image/ridges_fetal_vs_adult_tanaka.pdf",width = 6,height = 3)

# try the breakdown for cluster identity on the astro assignament
bind_rows(list(df_cortex,df_fetal)) %>% 
  filter(annotation_harmony=="ASTRO") %>% 
  mutate(seurat_clusters = factor(seurat_clusters)) %>% 
  ggplot(aes(y=seurat_clusters,x=score))+
  ggridges::geom_density_ridges(alpha=0.6)+
  theme_bw() + facet_wrap(~ref) +theme(strip.background = element_blank())+
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
ggsave("out/image/ridges_fetal_vs_adult2_tanaka.pdf",width = 6,height = 5)

# save the table version of the image
df_summary <- bind_rows(list(df_cortex,df_fetal)) %>% 
  mutate(thr_cat = case_when(score >= 0.75~"high",
                             T~"low")) %>% 
  group_by(annotation_harmony,ref,thr_cat) %>% 
  summarise(n=n())

df_summary %>% 
  write_tsv("out/table/df_summary_BS_fetal_vs_adult_tanaka.tsv")
