# libraries ---------------------------------------------------------------
library(tidyverse)
library(monocle)
library(cowplot)

# read in the data --------------------------------------------------------
monocle_df <- readRDS("out/object/monocle_df_POOLandRR16_CTRL.rds")

# -------------------------------------------------------------------------
# 
test <- plot_cell_trajectory(monocle_df, color_by = "seurat_clusters")

test$layers[[1]]$aes_params$alpha <- 0.5
test$layers[[3]]$aes_params$alpha <- 0.5

test$data %>%
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=seurat_clusters),alpha=0.2)+
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol = 2))+theme_bw()
ggsave("out/image/trajectory_mono_POOLandRR16.pdf",height = 3,width = 5)

# single plot
test$data %>%
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=seurat_clusters),alpha=0.2)+
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol = 2))+theme_cowplot()
ggsave("out/image/trajectory_mono_POOLandRR16_alt.pdf",height = 3,width = 5)

# show the split over the whole sample
test$data %>%
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=seurat_clusters),alpha=0.2)+
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol = 2))+theme_cowplot()+facet_wrap(~clone)+theme(strip.background = element_blank())
ggsave("out/image/trajectory_mono_POOLandRR16_alt_split.pdf",height = 3,width = 9)

# -------------------------------------------------------------------------
# split plot by seurat clusters
test <- plot_cell_trajectory(monocle_df, color_by = "Pseudotime")

test$layers[[1]]$aes_params$alpha <- 0.5
test$layers[[3]]$aes_params$alpha <- 0.5

test$data %>%
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=Pseudotime),alpha=0.2) +
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  theme_bw() +  
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) + 
  viridis::scale_color_viridis(option = "H")+
  facet_wrap(~seurat_clusters)
ggsave("out/image/pseuditime_trajectory_POOLandRR16.pdf",height = 10,width = 14)

test$data %>%
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=Pseudotime),alpha=0.2) +
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  theme_bw() +  
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) + 
  viridis::scale_color_viridis(option = "H")+
  facet_grid(seurat_clusters~clone)
ggsave("out/image/pseuditime_trajectory_POOLandRR16_2.pdf",height = 20,width = 6)

# updated plots -----------------------------------------------------------
# split plot by seurat clusters
test <- plot_cell_trajectory(monocle_df, color_by = "Pseudotime")
saveRDS(test,file = "out/object/test_updated_POOLandRR16.rds")

test <- readRDS(file = "out/object/test_updated_POOLandRR16.rds")

test$layers[[1]]$aes_params$alpha <- 0.5
test$layers[[3]]$aes_params$alpha <- 0.5

test$data %>%
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=Pseudotime),alpha=0.2) +
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  theme_bw() +  
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) + 
  viridis::scale_color_viridis(option = "H")+
  facet_wrap(~seurat_clusters,scales = "free")
ggsave("out/image/pseuditime_trajectory_POOLandRR16_new.pdf",height = 10,width = 14)

# same as above no grid
test$data %>%
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=Pseudotime),alpha=0.2) +
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  theme(panel.border = element_rect(colour = "black", fill = NA)) + 
  viridis::scale_color_viridis(option = "H")+
  facet_wrap(~seurat_clusters,scales = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave("out/image/pseuditime_trajectory_POOLandRR16_new_alt.pdf",height = 10,width = 14)

# option panel 2
test$data %>%
  filter(seurat_clusters %in% c(3,10,4,6,16,7)) %>% 
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=Pseudotime),alpha=0.2) +
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  theme_bw() +  
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) + 
  viridis::scale_color_viridis(option = "H")+
  # facet_rep_grid(~seurat_clusters,repeat.tick.labels = T)
  facet_wrap(~seurat_clusters,ncol = 2,scales = "free")
ggsave("out/image/pseuditime_trajectory_POOLandRR16_new_2.pdf",height = 9,width = 7)

# same as above no Grid
test$data %>%
  filter(seurat_clusters %in% c(3,10,4,6,16,7)) %>% 
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=Pseudotime),alpha=0.2) +
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  theme(panel.border = element_rect(colour = "black", fill = NA)) + 
  viridis::scale_color_viridis(option = "H")+
  # facet_rep_grid(~seurat_clusters,repeat.tick.labels = T)
  facet_wrap(~seurat_clusters,ncol = 2,scales = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave("out/image/pseuditime_trajectory_POOLandRR16_new_2_alt.pdf",height = 9,width = 7)

test$data %>%
  filter(seurat_clusters %in% c(3,10,4,6,16,7)) %>% 
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=Pseudotime),alpha=0.2) +
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  theme(panel.border = element_rect(colour = "black", fill = NA)) + 
  viridis::scale_color_viridis(option = "H")+
  # facet_rep_grid(~seurat_clusters,repeat.tick.labels = T)
  facet_grid(seurat_clusters~clone,scales = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave("out/image/pseuditime_trajectory_POOLandRR16_new_2_alt2.pdf",height = 10,width = 7)

# matina asked to group the cells by seurat cluster
test$data %>%
  filter(seurat_clusters %in% c(3,10,4,6,16,7)) %>%
  mutate(cell_ID = case_when(seurat_clusters %in% c(3,10)~"NEU",
                             seurat_clusters %in% c(4,6)~"ASTRO",
                             seurat_clusters %in% c(16)~"OLIGO",
                             seurat_clusters %in% c(7)~"CYCLING")) %>% 
  ggplot()+
  geom_point(aes(x=data_dim_1,y=data_dim_2,col=Pseudotime),alpha=0.2) +
  test$layers[1]+
  test$layers[3]+
  test$layers[4]+
  theme(panel.border = element_rect(colour = "black", fill = NA)) + 
  viridis::scale_color_viridis(option = "H")+
  # facet_rep_grid(~seurat_clusters,repeat.tick.labels = T)
  facet_grid(cell_ID~clone,scales = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave("out/image/pseuditime_trajectory_POOLandRR16_new_2_alt3.pdf",height = 8,width = 7)
