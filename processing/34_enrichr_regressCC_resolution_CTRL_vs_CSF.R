# libraries ---------------------------------------------------------------
library(tidyverse)
library(enrichR)
library(scales)
library(patchwork)

# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Azim"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","Human_Gene_Atlas","Azimuth_Cell_Types_2021")

# GENE SELECTION ----------------------------------------------------------
res <- read_tsv("out/table/response_CSF_vs_CTRL_data.combined_fix_resolution_DoubletSinglet_harmonyMartina.tsv") %>%
  mutate(cluster = paste0("cluster_",cluster))

list_df <- split(res,f = res$cluster)

list_res_tot <- lapply(list_df, function(x){
  x %>%
    filter(p_val_adj<0.05,avg_log2FC>0.5) %>%
    pull(gene)
})

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs
list_res_tot_filter <- list_res_tot[lengths(list_res_tot)>5]

list <- lapply(list_res_tot_filter,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db)
  #
  out_enrich %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list %>%
  write_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")

plot_list <- list %>%
  split(f = .$comparison)

list_plot <- pmap(list(plot_list,names(plot_list)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:20) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
})

wrap_plots(list_plot,nrow = 1)
ggsave("out/image/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.pdf",width = 25,height = 15,limitsize = FALSE)

# make the plot following martina's requests
list_plot2 <- pmap(list(plot_list,names(plot_list)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:20) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    ungroup() %>% 
    mutate(Term = fct_reorder(Term, P.value,.desc = T)) %>%
    ggplot(aes(y=Term,col=annotation,size = Odds.Ratio,x = P.value)) + geom_point() + 
    # facet_wrap(~annotation,scales = "free",ncol = 1)+
    theme_bw() +
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
})

wrap_plots(list_plot2,nrow = 1)
ggsave("out/image/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina2.pdf",width = 25,height = 15,limitsize = FALSE)

# make the plot following martina's requests
list_plot3 <- pmap(list(plot_list,names(plot_list)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:20) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    ungroup() %>% 
    mutate(Term = fct_reorder(Term, P.value,.desc = T)) %>%
    filter(annotation %in% c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016")) %>% 
    ggplot(aes(y=Term,col=annotation,size = Odds.Ratio,x = P.value)) + geom_point() + 
    theme_bw() +
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA),legend.position = "bottom",legend.box="vertical")+
    guides(color=guide_legend(nrow=3, byrow=TRUE))+
    ggtitle(y)
})

wrap_plots(list_plot3,nrow = 1)
ggsave("out/image/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina3.pdf",width = 20,height = 10,limitsize = FALSE)

