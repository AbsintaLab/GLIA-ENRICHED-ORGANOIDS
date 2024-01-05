# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(limma)
library(multtest)
library(metap)
library(ggbreak)
library(patchwork)

# read in the dataset -----------------------------------------------------
# after the meetign we decided to fix the analysis of harmony martina to run the rest of the analysis
data.combined <- readRDS(file = "data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
DimPlot(data.combined,label = T)

# summary for martina
data.combined@meta.data %>% 
  group_by(sample) %>% 
  summarise(n=n(),avg_count = mean(nCount_RNA),avg_feature = mean(nFeature_RNA),avg_mito_prc=mean(percent.mt))

# wrangling ---------------------------------------------------------------
# add the treatmend varibale to the metadata
meta <- data.combined@meta.data %>%
  rownames_to_column("barcodes")

# find markers ------------------------------------------------------------
# perform the DE using the cluster identity
head(data.combined)

# confirm the correct ident
Idents(data.combined) <- "seurat_clusters"

# make sure the dafault assay is RNA and make sure the data have been scaled also on the RNA slot
DefaultAssay(data.combined)
dim(data.combined@assays$RNA@scale.data)
data.combined@assays$RNA@scale.data[1:10,1:10]

#
data.combined.markers <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
data.combined.markers %>%
  write_tsv("out/table/FindAllMarkers_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_seuratClusters.tsv")

# top 100 per cluster
data.combined.markers %>%
  group_by(cluster) %>%
  mutate(rank = rank(order(p_val_adj, -abs(avg_log2FC)), ties.method='first')) %>%
  arrange(cluster,rank) %>%
  filter(rank < 101) %>%
  write_tsv("out/table/FindAllMarkers_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_seuratClusters_top100.tsv")

# find conserved markers --------------------------------------------------
# try the other approach using the
# the issue is that some of the clusters have very few cells in the old condition that is throwing an error

# run the identification of the markers across all the clusters
cluster_id <- sort(unique(data.combined@meta.data$seurat_clusters))

# explore the number of cells per orig.ident
data.combined@meta.data %>%
  group_by(seurat_clusters,orig.ident) %>%
  summarise(n = n()) %>%
  # print(n=50)
  pivot_wider(names_from = seurat_clusters,values_from = n)

# unfortunetaly the conserved genes per orig.ident might be problematic
# explore the number of cells per clone
data.combined@meta.data %>%
  group_by(seurat_clusters,clone) %>%
  summarise(n = n()) %>%
  # print(n=50)
  pivot_wider(names_from = seurat_clusters,values_from = n)

# explore the number of cells per treatment
data.combined@meta.data %>%
  group_by(seurat_clusters,treat) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = seurat_clusters,values_from = n)

# define the grouping for the test
Idents(data.combined) <- "seurat_clusters"

# try to keep all the cluster despite the splitting for orig.ident
cluster_id2 <- cluster_id

# test the conserved markers across all the clones
list_conserved_markers_sample_01 <- lapply(cluster_id2, function(x){
  FindConservedMarkers(data.combined, ident.1 = x, grouping.var = "orig.ident", verbose = T,min.cells.group = 0)
}) %>%
  setNames(cluster_id2) %>%
  bind_rows(.id = "cluster_number")
# 
list_conserved_markers_sample_01 %>%
  rownames_to_column("gene") %>%
  separate(gene,into = c("gene","id"),sep = "\\.\\.") %>%
  write_tsv("out/table/FindConservedMarkers_OrigIdent_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_seuratClusters.tsv")
# 
list_conserved_markers_sample_01 %>%
  rownames_to_column("gene") %>%
  separate(gene,into = c("gene","id"),sep = "\\.\\.") %>%
  group_by(cluster_number,gene) %>%
  mutate(avg_log2FC = mean(c(CTRL8_WT_CSF_avg_log2FC,CTRL4_SOX10_CTRL_avg_log2FC,CTRL8_WT_CTRL_avg_log2FC,CTRL8_SOX10_CSF_avg_log2FC,RR25_SOX10_CSF_avg_log2FC,RR16_SOX10_CTRL_avg_log2FC,pool_MOCK_CTRL_avg_log2FC,CTRL8_SOX10_CTRL_avg_log2FC,RR16_SOX10_CSF_avg_log2FC))) %>%
  filter(avg_log2FC>0) %>%
  arrange(cluster_number,max_pval,desc(avg_log2FC)) %>%
  write_tsv("out/table/FindConservedMarkers_OrigIdent_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_seuratClusters_reordered.tsv")
# 

GOI_01 <- read_tsv("out/table/FindAllMarkers_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_seuratClusters.tsv") %>%
  filter(str_detect(gene,pattern = "MT-|-AS1|\\.",negate = T)) %>%
  mutate(delta_prc = pct.1-pct.2) %>% 
  select(gene,cluster,p_val,avg_log2FC,delta_prc) %>%
  group_by(cluster) %>%
  arrange(cluster,p_val,desc(avg_log2FC),desc(delta_prc)) %>%
  # mutate(rank = rank(pval)) %>%
  # top_n(15,wt = -rank) %>%
  dplyr::slice(1:20) %>% 
  pull(gene) %>%
  unique()

read_tsv("out/table/FindAllMarkers_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_seuratClusters.tsv") %>%
  filter(str_detect(gene,pattern = "MT-|-AS1|\\.",negate = T)) %>%
  mutate(delta_prc = pct.1-pct.2) %>% 
  select(gene,cluster,p_val,avg_log2FC,delta_prc) %>%
  group_by(cluster) %>%
  arrange(cluster,p_val,desc(avg_log2FC),desc(delta_prc)) %>%
  # mutate(rank = rank(pval)) %>%
  # top_n(15,wt = -rank) %>%
  dplyr::slice(1:20) %>% 
  write_tsv("out/table/GOI_01.tsv")

# merge the list of markers
GOI <- c(GOI_01)

#produce the heatmap
test1 <- DoHeatmap(data.combined,group.by = "seurat_clusters", features = GOI) + NoLegend()
# # beware! the Identity column is not correct
data_heatmap_long <- test1$data %>%
  drop_na() %>%
  select(-Identity)

# add the original meta
data_heatmap_long_fix <- left_join(data_heatmap_long,meta,by = c("Cell"="barcodes"))

mat <- data_heatmap_long_fix %>%
  arrange(seurat_clusters,treat) %>%
  dplyr::select(Feature,Cell,Expression) %>%
  pivot_wider(names_from = Cell,values_from = Expression) %>%
  column_to_rownames("Feature")

ref_meta <- meta %>%
  dplyr::slice(match(colnames(mat),meta$barcodes))

id_clusters <- names(table(meta$seurat_clusters))
id_colors <- hue_pal()(length(id_clusters))

#
names(id_colors) <- id_clusters
column_ha <- HeatmapAnnotation(cluster = ref_meta$seurat_clusters,
                               sample = ref_meta$orig.ident,
                               col = list(sample = c("CTRL4_SOX10_CTRL" = "#466BE3FF",
                                                     "CTRL8_WT_CTRL"= "#25ECA7FF",
                                                     "CTRL8_WT_CSF"= "#62FC6BFF",
                                                     "CTRL8_SOX10_CSF"= "#28BBECFF", 
                                                     "CTRL8_SOX10_CTRL"= "#31F299FF",
                                                     "pool_MOCK_CTRL"= "#30123BFF",
                                                     "RR16_SOX10_CSF"= "#D23105FF",
                                                     "RR16_SOX10_CTRL"= "#7A0403FF",
                                                     "RR25_SOX10_CSF"= "#FB8022FF"),
                                          cluster = id_colors))

id_lab <- which(rownames(mat) %in% c("GFAP", "S100b", "AQP4", "VIM", "CLU", "SPARCL1","B2M","RPS6", "RPS27", "EGR1","TOP2A", "MKI67", "PBK", "NUSAP1", "CENPU","SEZ6L2", "NSG2", "SYT4", "SERPINI1","NEFL", "STMN2", "BASP1", "MAP1B","NEUROD6", "NEUROD1", "NEUROD2", "NHLH1", "BHLHE22","PLP1","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG","AQP4", "GFAP", "VIM","SLC1A2","S100B","GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"))
name_lab <- rownames(mat)[id_lab]

ha = rowAnnotation(link = anno_mark(at = id_lab, labels = name_lab))

ht <- Heatmap(mat,
              name = "z-score",
              # col = rainbow(20),
              # col = viridis::turbo(20),
              col = viridis::inferno(20),
              # col = viridis::magma(20),
              # col = fun_col,
              use_raster = F,
              show_row_names = F,
              show_column_names = F,
              cluster_columns = F,
              right_annotation = ha,
              column_split = ref_meta$seurat_clusters,
              # row_km = 8,
              top_annotation = column_ha)

pdf("out/image/heatmap_panel_expression_FindConservedMarkers_seurat_clusters_harmonyMartina.pdf",width = 15,height = 12)
draw(ht,heatmap_legend_side = "left",annotation_legend_side = "left")
dev.off()

# compare CSF vs CTRL -----------------------------------------------------
data.combined
# run it on RR16
# make suer the correct defalult dataset il loaded should be RNA
DefaultAssay(data.combined)

# check that in the RNA slot the data object is indeed loaded with normalized values
data.combined@assays$RNA@data[1:10,1:10]

# define the grouping variables for the comparison of CSF vs CTRL
head(data.combined@meta.data)

# define the new grouping
data.combined$treat.cluster <- paste(data.combined$treat,data.combined$seurat_clusters, sep = "_")
head(data.combined@meta.data)

data.combined$treat.cluster.test <- paste(data.combined$treat,
                                          data.combined$seurat_clusters,
                                          case_when(data.combined$orig.ident %in% c("CTRL4_SOX10_CTRL","CTRL8_SOX10_CSF","CTRL8_SOX10_CTRL","pool_MOCK_CTRL","RR25_SOX10_CSF")~"IN",
                                                    T~"OUT"),sep = "_")

# update the idents of the object
Idents(data.combined) <- "treat.cluster"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# run the differential expression over the WT
clusters_id <- as.character(sort(unique(data.combined$seurat_clusters)))

# the option below is neede to avoid errors for max size
options(future.globals.maxSize = 8000 * 1024^2)

table(data.combined@meta.data$seurat_clusters,data.combined@meta.data$treat)

# comparison 1
# skip RPE and IMM_14 form the compariosn as there are too few cells
clusters_id2 <- clusters_id
# str_subset(pattern = "IMM_14|RPE",negate = T)

list_CSF_vs_CTRL <- lapply(clusters_id2,function(x){
  id_1 <- paste0("CSF_",x)
  id_2 <- paste0("CTRL_",x)
  response <- FindMarkers(data.combined,
                          ident.1 = id_1,
                          ident.2 = id_2,
                          verbose = T,
                          logfc.threshold = 0,
                          min.pct = 0.1)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_CSF_vs_CTRL %>%
  setNames(clusters_id2) %>%
  bind_rows(.id = "cell_id") %>%
  write_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv")

# update the idents of the object
Idents(data.combined) <- "seurat_clusters"
DimPlot(data.combined,split.by = "orig.ident")

FeaturePlot(data.combined,features = "SPP1",split.by = "orig.ident")

# the reason why the SPP1 is not popping up in the CSF vs CTRL anymore is due to the presence of the cells from RR16 and CTRL8 WT
Idents(data.combined) <- "treat.cluster.test"

test <- FindMarkers(data.combined,
                        ident.1 = "CSF_12_IN",
                        ident.2 = "CTRL_12_IN",
                        verbose = T,
                        logfc.threshold = 0,
                        min.pct = 0.1)

test["SPP1",]

test2 <- FindMarkers(data.combined,
                    ident.1 = "CSF_13_IN",
                    ident.2 = "CTRL_13_IN",
                    verbose = T,
                    logfc.threshold = 0,
                    min.pct = 0.1)

test2["SPP1",]

test3 <- FindMarkers(data.combined,
                     ident.1 = "CSF_12_OUT",
                     ident.2 = "CTRL_12_OUT",
                     verbose = T,
                     logfc.threshold = 0,
                     min.pct = 0.1)

test3["SPP1",]

test4 <- FindMarkers(data.combined,
                     ident.1 = "CSF_13_OUT",
                     ident.2 = "CTRL_13_OUT",
                     verbose = T,
                     logfc.threshold = 0,
                     min.pct = 0.1)

test4["SPP1",]

# update the idents of the object
table(data.combined@meta.data$seurat_clusters,data.combined@meta.data$orig.ident)

# explore outputs comparison ----------------------------------------------
# check the position of genes francesca submitted
read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  split(f = .$cluster) %>%
  lapply(.,function(x){
    x %>%
      mutate(dPct = pct.1 - pct.2) %>%
      filter(avg_log2FC > 0.5,
             p_val_adj < 0.01)
  })

# explore the list of DE gens
test <- read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  split(f = .$cluster) %>%
  lapply(.,function(x){
    x %>%
      mutate(dPct = pct.1 - pct.2) %>%
      filter(avg_log2FC > 0.5,
             p_val_adj < 0.01)
  })

read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  filter(gene %in% c("SPP1"))

read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(DE_cat != "no") %>%
  mutate(cluster=factor(cluster)) %>%
  group_by(cluster,DE_cat) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = fct_reorder(cluster,tot,.desc = T)) %>%
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  scale_fill_manual(values = c("blue","red"))+
  theme_bw()
ggsave("out/image/DE_barplot_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.pdf",width = 8,height = 4)

# try to plot woth breaks
read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(DE_cat != "no") %>%
  mutate(cluster=factor(cluster)) %>%
  group_by(cluster,DE_cat) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = fct_reorder(cluster,tot,.desc = T)) %>%
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  scale_fill_manual(values = c("blue","red"))+
  theme_bw()+
  scale_y_break(breaks = c(35,120),scales = 1)+
  scale_y_break(breaks = c(230,400),scales = 0.2)
ggsave("out/image/DE_barplot_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_break.pdf",width = 8,height = 4)

# try to plot woth breaks
read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(DE_cat != "no") %>%
  mutate(cluster=factor(cluster)) %>%
  group_by(cluster,DE_cat) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = fct_reorder(cluster,tot,.desc = T)) %>%
  mutate(split_tab = case_when(cluster%in%c(13,12)~"high",
                               T~"low")) %>% 
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  scale_fill_manual(values = c("blue","red"))+
  theme_bw()+
  facet_wrap(~split_tab,scales = "free")+
  theme(strip.background = element_blank(),strip.text=element_blank())

# make the plot more balanced
list_plot <- read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(DE_cat != "no") %>%
  mutate(cluster=factor(cluster)) %>%
  group_by(cluster,DE_cat) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(cluster = fct_reorder(cluster,tot,.desc = T)) %>%
  mutate(split_tab = case_when(cluster%in%c(13,12)~"high",
                               T~"low")) %>% 
  split(.$split_tab)
  
p1 <- list_plot$high %>% 
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  scale_fill_manual(values = c("blue","red"))+
  theme_bw()+
  facet_wrap(~split_tab,scales = "free")+
  theme(strip.background = element_blank(),strip.text=element_blank(),legend.position = "none")

p2 <- list_plot$low %>% 
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  scale_fill_manual(values = c("blue","red"))+
  theme_bw()+
  facet_wrap(~split_tab,scales = "free")+
  theme(strip.background = element_blank(),strip.text=element_blank())

p1+p2+plot_layout(widths = c(1,3))

ggsave("out/image/DE_barplot_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_split.pdf",width = 8,height = 4)

# same as above without the grid
p1_alt <- list_plot$high %>% 
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  scale_fill_manual(values = c("blue","red"))+
  theme_cowplot()+
  facet_wrap(~split_tab,scales = "free")+
  theme(strip.background = element_blank(),strip.text=element_blank(),legend.position = "none")

p2_alt <- list_plot$low %>% 
  ggplot(aes(x=cluster,y=n,fill=DE_cat)) + geom_col() +
  scale_fill_manual(values = c("blue","red"))+
  theme_cowplot()+
  facet_wrap(~split_tab,scales = "free")+
  theme(strip.background = element_blank(),strip.text=element_blank())

p1_alt+p2_alt+plot_layout(widths = c(1,3))

ggsave("out/image/DE_barplot_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina_split_alt.pdf",width = 8,height = 4)


# plot the volcano for the cluster 12 DE
volcano_tot <- read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  filter(cluster %in% c(13,12))

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+theme_bw()+facet_wrap(~cluster)+theme(strip.background = element_blank())
ggsave("out/image/volcano_data.combined_annotated_norm_fix_cell_resolution_12_13_CSFvsCTRL_harmonyMartina.pdf",width = 12,height = 9)

# plot the volcano for every cluster
volcano_tot <- read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no"))

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+theme_bw()+facet_wrap(~cluster,scales = "free")+theme(strip.background = element_blank())
ggsave("out/image/volcano_data.combined_annotated_norm_fix_cell_resolution_ALL_CSFvsCTRL_harmonyMartina.pdf",width = 20,height = 16)

