# libraries ---------------------------------------------------------------
library(monocle)
library(RColorBrewer)

# read in the data --------------------------------------------------------
monocle_df <- readRDS("out/object/monocle_df_RR16_CTRL.rds")

# wrangling ---------------------------------------------------------------
BEAM_res_monocle <- BEAM(monocle_df, branch_point = 1, cores = 12)
BEAM_res_monocle <- BEAM_res_monocle[order(BEAM_res_monocle$qval),]
BEAM_res_monocle <- BEAM_res_monocle[,c("gene_short_name", "pval", "qval")]

saveRDS(BEAM_res_monocle,"out/object/BEAM_res_monocle_branch01_RR16.rds")

plot_genes_branched_heatmap(monocle_df[row.names(subset(BEAM_res_monocle,
                                                        qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)

pdf("out/image/heatmap_branch_01_RR16.pdf",width = 4,height = 9)
plot_genes_branched_heatmap(monocle_df[row.names(subset(BEAM_res_monocle,
                                                        qval == 0)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T,
                            hmcols = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 5, name ="RdYlBu")))(60))
dev.off()

# -------------------------------------------------------------------------
BEAM_res_monocle2 <- BEAM(monocle_df, branch_point = 2, cores = 12)
BEAM_res_monocle2 <- BEAM_res_monocle2[order(BEAM_res_monocle2$qval),]
BEAM_res_monocle2 <- BEAM_res_monocle2[,c("gene_short_name", "pval", "qval")]

saveRDS(BEAM_res_monocle2,"out/object/BEAM_res_monocle_branch02_RR16.rds")

pdf("out/image/heatmap_branch_02_RR16.pdf",width = 4,height = 6)
plot_genes_branched_heatmap(monocle_df[row.names(subset(BEAM_res_monocle2,qval < 1e-170)),],
                            branch_point = 2,
                            num_clusters = 3,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T,
                            hmcols = colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(60))
dev.off()

scales::show_col(colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(20))
