library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(openxlsx)
library(ReactomePA)
library(org.Hs.eg.db)
library(pathfindR)
library(cowplot)
library(DESeq2)
library(fastcluster)
library(parallel)
library(foreach)
library(doParallel)

# Import the full pathway results
# load("Resources_and_environments/Pathways.RData")

# ashr_TR_vs_TNR #####
nonzero_res = which(lapply(pathfindR_outputs_ashr_TR_vs_TNR, nrow) != 0)
enrichment_dotplots_ashr_TR_vs_TNR = list()
cluster_enrichment_dotplots_ashr_TR_vs_TNR = list()

# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  # unclustered results
  enrichment_dotplots_ashr_TR_vs_TNR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[names]],
                                                                 top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    scale_size(range = c(0.1, 2)) +
    theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 3.5),
          axis.title.x = element_text(size = 4.5, face = "bold"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 3.5, face = "bold"),
          legend.text = element_text(size = 2.75))+
    labs(title = paste0("Top 10 ", names,
                        " terms enrichment dotplot - (", comparisons["ashr_TR_vs_TNR"], ")"))
  print(enrichment_dotplots_ashr_TR_vs_TNR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.png"),
         path = paste0("Reports/png_pathways/ashr_TR_vs_TNR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_TR_vs_TNR)) {
    cluster_enrichment_dotplots_ashr_TR_vs_TNR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[names]][clustered_results_ashr_TR_vs_TNR[[names]]$Status
                                                                                                                                                   == "Representative", ][1:10,],
                                                                           top_terms = NULL,
                                                                           plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      scale_size(range = c(0.1, 2)) +
      theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
            axis.text.y = element_text(color = "black", size = 3.5),
            axis.text.x = element_text(color = "black", size = 3.5),
            axis.title.x = element_text(size = 4.5, face = "bold"),
            legend.key.size = unit(1.5, units = "mm"),
            legend.spacing.y = unit(0.5, units = "mm"),
            legend.spacing.x = unit(0.5, units = "mm"),
            legend.title = element_text(size = 3.5, face = "bold"),
            legend.text = element_text(size = 2.75))+
      labs(title = paste0("Top 10 clustered ", names,
                          " terms enrichment dotplot - (", comparisons["ashr_TR_vs_TNR"], ")"))
    print(cluster_enrichment_dotplots_ashr_TR_vs_TNR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.png"),
           path = paste0("Reports/png_pathways/ashr_TR_vs_TNR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
}

term_gene_heatmaps_ashr_TR_vs_TNR = list()
term_gene_graphs_ashr_TR_vs_TNR = list()

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  # term-gene heatmaps
  term_gene_heatmaps_ashr_TR_vs_TNR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[names]],
                                                                 genes_df = pathf_input_ashr_TR_vs_TNR,
                                                                 num_terms = 10,
                                                                 use_description = TRUE,
                                                                 low = "darkblue",
                                                                 high = "red",
                                                                 mid = "white",
                                                                 pin_name_path = 'Biogrid')+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 4),
          axis.text.x = element_text(color = "black", size = 4,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names, 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_TR_vs_TNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  #print(ggdraw(align_legend(term_gene_heatmaps_ashr_TR_vs_TNR[[names]], hjust = 0.5)))
  print(term_gene_heatmaps_ashr_TR_vs_TNR[[names]])
  ggsave(filename = paste0(names,
                           "_top10_term_gene_heatmap.png"),
         path = paste0("Reports/png_pathways/ashr_TR_vs_TNR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_TR_vs_TNR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_TR_vs_TNR[[names]],
                                                             num_terms = 5,
                                                             use_description = TRUE,
                                                             node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 2, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names, 
                        " term - gene graph (",
                        comparisons["ashr_TR_vs_TNR"], ")"))
  term_gene_graphs_ashr_TR_vs_TNR[[names]]$layers = list(term_gene_graphs_ashr_TR_vs_TNR[[names]]$layers[[1]], 
                                                         term_gene_graphs_ashr_TR_vs_TNR[[names]]$layers[[2]], 
                                                         term_gene_graphs_ashr_TR_vs_TNR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_TR_vs_TNR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.png"),
         path = paste0("Reports/png_pathways/ashr_TR_vs_TNR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_TR_vs_TNR = list()
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  # UpSet plot
  UpSet_plots_ashr_TR_vs_TNR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[names]],
                                                   genes_df = pathf_input_ashr_TR_vs_TNR,
                                                   num_terms = 5,
                                                   use_description = TRUE,
                                                   low = "darkgreen",
                                                   high = "darkred",
                                                   mid = "white")+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 3),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                              combmatrix.panel.point.size = 0.8,
                              combmatrix.panel.line.size = 0.5)+
    labs(fill = "logFC",
         title = paste0("Top 5 ", names, 
                        " term - UpSet plot (",
                        comparisons["ashr_TR_vs_TNR"], ")"))
  tiff(paste0("Reports/png_pathways/ashr_TR_vs_TNR/", names, "/",
              names, "_top5_UpSet_plot.png"), 
       width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_TR_vs_TNR[[names]])
  dev.off()
}

# ashr_LTR_vs_ETR #####
nonzero_res = which(lapply(pathfindR_outputs_ashr_LTR_vs_ETR, nrow) != 0)
enrichment_dotplots_ashr_LTR_vs_ETR = list()
cluster_enrichment_dotplots_ashr_LTR_vs_ETR = list()

# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  # unclustered results
  enrichment_dotplots_ashr_LTR_vs_ETR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[names]],
                                                                  top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    scale_size(range = c(0.1, 2)) +
    theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 3.5),
          axis.title.x = element_text(size = 4.5, face = "bold"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 3.5, face = "bold"),
          legend.text = element_text(size = 2.75))+
    labs(title = paste0("Top 10 ", names,
                        " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_ETR"], ")"))
  print(enrichment_dotplots_ashr_LTR_vs_ETR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.png"),
         path = paste0("Reports/png_pathways/ashr_LTR_vs_ETR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_LTR_vs_ETR)) {
    cluster_enrichment_dotplots_ashr_LTR_vs_ETR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[names]][clustered_results_ashr_LTR_vs_ETR[[names]]$Status
                                                                                                                                                     == "Representative", ][1:10,],
                                                                            top_terms = NULL,
                                                                            plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      scale_size(range = c(0.1, 2)) +
      theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
            axis.text.y = element_text(color = "black", size = 3.5),
            axis.text.x = element_text(color = "black", size = 3.5),
            axis.title.x = element_text(size = 4.5, face = "bold"),
            legend.key.size = unit(1.5, units = "mm"),
            legend.spacing.y = unit(0.5, units = "mm"),
            legend.spacing.x = unit(0.5, units = "mm"),
            legend.title = element_text(size = 3.5, face = "bold"),
            legend.text = element_text(size = 2.75))+
      labs(title = paste0("Top 10 clustered ", names,
                          " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_ETR"], ")"))
    print(cluster_enrichment_dotplots_ashr_LTR_vs_ETR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.png"),
           path = paste0("Reports/png_pathways/ashr_LTR_vs_ETR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
}

term_gene_heatmaps_ashr_LTR_vs_ETR = list()
term_gene_graphs_ashr_LTR_vs_ETR = list()

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  # term-gene heatmaps
  term_gene_heatmaps_ashr_LTR_vs_ETR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[names]],
                                                                  genes_df = pathf_input_ashr_LTR_vs_ETR,
                                                                  num_terms = 10,
                                                                  use_description = TRUE,
                                                                  low = "darkblue",
                                                                  high = "red",
                                                                  mid = "white",
                                                                  pin_name_path = 'Biogrid')+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 4),
          axis.text.x = element_text(color = "black", size = 4,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names, 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_LTR_vs_ETR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  #print(ggdraw(align_legend(term_gene_heatmaps_ashr_LTR_vs_ETR[[names]], hjust = 0.5)))
  print(term_gene_heatmaps_ashr_LTR_vs_ETR[[names]])
  ggsave(filename = paste0(names,
                           "_top10_term_gene_heatmap.png"),
         path = paste0("Reports/png_pathways/ashr_LTR_vs_ETR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_LTR_vs_ETR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_LTR_vs_ETR[[names]],
                                                              num_terms = 5,
                                                              use_description = TRUE,
                                                              node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 2, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names, 
                        " term - gene graph (",
                        comparisons["ashr_LTR_vs_ETR"], ")"))
  term_gene_graphs_ashr_LTR_vs_ETR[[names]]$layers = list(term_gene_graphs_ashr_LTR_vs_ETR[[names]]$layers[[1]], 
                                                          term_gene_graphs_ashr_LTR_vs_ETR[[names]]$layers[[2]], 
                                                          term_gene_graphs_ashr_LTR_vs_ETR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_LTR_vs_ETR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.png"),
         path = paste0("Reports/png_pathways/ashr_LTR_vs_ETR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_LTR_vs_ETR = list()
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  # UpSet plot
  UpSet_plots_ashr_LTR_vs_ETR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[names]],
                                                    genes_df = pathf_input_ashr_LTR_vs_ETR,
                                                    num_terms = 5,
                                                    use_description = TRUE,
                                                    low = "darkgreen",
                                                    high = "darkred",
                                                    mid = "white")+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 3),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                              combmatrix.panel.point.size = 0.8,
                              combmatrix.panel.line.size = 0.5)+
    labs(fill = "logFC",
         title = paste0("Top 5 ", names, 
                        " term - UpSet plot (",
                        comparisons["ashr_LTR_vs_ETR"], ")"))
  tiff(paste0("Reports/png_pathways/ashr_LTR_vs_ETR/", names, "/",
              names, "_top5_UpSet_plot.png"), 
       width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_LTR_vs_ETR[[names]])
  dev.off()
}

# ashr_LTR_vs_LTNR #####
nonzero_res = which(lapply(pathfindR_outputs_ashr_LTR_vs_LTNR, nrow) != 0)
enrichment_dotplots_ashr_LTR_vs_LTNR = list()
cluster_enrichment_dotplots_ashr_LTR_vs_LTNR = list()

# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  # unclustered results
  enrichment_dotplots_ashr_LTR_vs_LTNR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]],
                                                                   top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    scale_size(range = c(0.1, 2)) +
    theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 3.5),
          axis.title.x = element_text(size = 4.5, face = "bold"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 3.5, face = "bold"),
          legend.text = element_text(size = 2.75))+
    labs(title = paste0("Top 10 ", names,
                        " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_LTNR"], ")"))
  print(enrichment_dotplots_ashr_LTR_vs_LTNR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.png"),
         path = paste0("Reports/png_pathways/ashr_LTR_vs_LTNR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_LTR_vs_LTNR)) {
    cluster_enrichment_dotplots_ashr_LTR_vs_LTNR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]][clustered_results_ashr_LTR_vs_LTNR[[names]]$Status
                                                                                                                                                       == "Representative", ][1:10,],
                                                                             top_terms = NULL,
                                                                             plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      scale_size(range = c(0.1, 2)) +
      theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
            axis.text.y = element_text(color = "black", size = 3.5),
            axis.text.x = element_text(color = "black", size = 3.5),
            axis.title.x = element_text(size = 4.5, face = "bold"),
            legend.key.size = unit(1.5, units = "mm"),
            legend.spacing.y = unit(0.5, units = "mm"),
            legend.spacing.x = unit(0.5, units = "mm"),
            legend.title = element_text(size = 3.5, face = "bold"),
            legend.text = element_text(size = 2.75))+
      labs(title = paste0("Top 10 clustered ", names,
                          " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_LTNR"], ")"))
    print(cluster_enrichment_dotplots_ashr_LTR_vs_LTNR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.png"),
           path = paste0("Reports/png_pathways/ashr_LTR_vs_LTNR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
}

term_gene_heatmaps_ashr_LTR_vs_LTNR = list()
term_gene_graphs_ashr_LTR_vs_LTNR = list()

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  # term-gene heatmaps
  term_gene_heatmaps_ashr_LTR_vs_LTNR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]],
                                                                   genes_df = pathf_input_ashr_LTR_vs_LTNR,
                                                                   num_terms = 10,
                                                                   use_description = TRUE,
                                                                   low = "darkblue",
                                                                   high = "red",
                                                                   mid = "white",
                                                                   pin_name_path = 'Biogrid')+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 4),
          axis.text.x = element_text(color = "black", size = 4,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names, 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  #print(ggdraw(align_legend(term_gene_heatmaps_ashr_LTR_vs_LTNR[[names]], hjust = 0.5)))
  print(term_gene_heatmaps_ashr_LTR_vs_LTNR[[names]])
  ggsave(filename = paste0(names,
                           "_top10_term_gene_heatmap.png"),
         path = paste0("Reports/png_pathways/ashr_LTR_vs_LTNR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_LTR_vs_LTNR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_LTR_vs_LTNR[[names]],
                                                               num_terms = 5,
                                                               use_description = TRUE,
                                                               node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 2, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names, 
                        " term - gene graph (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"))
  term_gene_graphs_ashr_LTR_vs_LTNR[[names]]$layers = list(term_gene_graphs_ashr_LTR_vs_LTNR[[names]]$layers[[1]], 
                                                           term_gene_graphs_ashr_LTR_vs_LTNR[[names]]$layers[[2]], 
                                                           term_gene_graphs_ashr_LTR_vs_LTNR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_LTR_vs_LTNR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.png"),
         path = paste0("Reports/png_pathways/ashr_LTR_vs_LTNR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_LTR_vs_LTNR = list()
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  # UpSet plot
  UpSet_plots_ashr_LTR_vs_LTNR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]],
                                                     genes_df = pathf_input_ashr_LTR_vs_LTNR,
                                                     num_terms = 5,
                                                     use_description = TRUE,
                                                     low = "darkgreen",
                                                     high = "darkred",
                                                     mid = "white")+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 3),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                              combmatrix.panel.point.size = 0.8,
                              combmatrix.panel.line.size = 0.5)+
    labs(fill = "logFC",
         title = paste0("Top 5 ", names, 
                        " term - UpSet plot (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"))
  tiff(paste0("Reports/png_pathways/ashr_LTR_vs_LTNR/", names, "/",
              names, "_top5_UpSet_plot.png"), 
       width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_LTR_vs_LTNR[[names]])
  dev.off()
}

# ashr_ETR_vs_ETNR #####
nonzero_res = which(lapply(pathfindR_outputs_ashr_ETR_vs_ETNR, nrow) != 0)
enrichment_dotplots_ashr_ETR_vs_ETNR = list()
cluster_enrichment_dotplots_ashr_ETR_vs_ETNR = list()

# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  # unclustered results
  enrichment_dotplots_ashr_ETR_vs_ETNR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]],
                                                                   top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    scale_size(range = c(0.1, 2)) +
    theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 3.5),
          axis.title.x = element_text(size = 4.5, face = "bold"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 3.5, face = "bold"),
          legend.text = element_text(size = 2.75))+
    labs(title = paste0("Top 10 ", names,
                        " terms enrichment dotplot - (", comparisons["ashr_ETR_vs_ETNR"], ")"))
  print(enrichment_dotplots_ashr_ETR_vs_ETNR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.png"),
         path = paste0("Reports/png_pathways/ashr_ETR_vs_ETNR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_ETR_vs_ETNR)) {
    cluster_enrichment_dotplots_ashr_ETR_vs_ETNR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]][clustered_results_ashr_ETR_vs_ETNR[[names]]$Status
                                                                                                                                                       == "Representative", ][1:10,],
                                                                             top_terms = NULL,
                                                                             plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      scale_size(range = c(0.1, 2)) +
      theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
            axis.text.y = element_text(color = "black", size = 3.5),
            axis.text.x = element_text(color = "black", size = 3.5),
            axis.title.x = element_text(size = 4.5, face = "bold"),
            legend.key.size = unit(1.5, units = "mm"),
            legend.spacing.y = unit(0.5, units = "mm"),
            legend.spacing.x = unit(0.5, units = "mm"),
            legend.title = element_text(size = 3.5, face = "bold"),
            legend.text = element_text(size = 2.75))+
      labs(title = paste0("Top 10 clustered ", names,
                          " terms enrichment dotplot - (", comparisons["ashr_ETR_vs_ETNR"], ")"))
    print(cluster_enrichment_dotplots_ashr_ETR_vs_ETNR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.png"),
           path = paste0("Reports/png_pathways/ashr_ETR_vs_ETNR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
}

term_gene_heatmaps_ashr_ETR_vs_ETNR = list()
term_gene_graphs_ashr_ETR_vs_ETNR = list()

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  # term-gene heatmaps
  term_gene_heatmaps_ashr_ETR_vs_ETNR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]],
                                                                   genes_df = pathf_input_ashr_ETR_vs_ETNR,
                                                                   num_terms = 10,
                                                                   use_description = TRUE,
                                                                   low = "darkblue",
                                                                   high = "red",
                                                                   mid = "white",
                                                                   pin_name_path = 'Biogrid')+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 4),
          axis.text.x = element_text(color = "black", size = 4,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names, 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  #print(ggdraw(align_legend(term_gene_heatmaps_ashr_ETR_vs_ETNR[[names]], hjust = 0.5)))
  print(term_gene_heatmaps_ashr_ETR_vs_ETNR[[names]])
  ggsave(filename = paste0(names,
                           "_top10_term_gene_heatmap.png"),
         path = paste0("Reports/png_pathways/ashr_ETR_vs_ETNR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_ETR_vs_ETNR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_ETR_vs_ETNR[[names]],
                                                               num_terms = 5,
                                                               use_description = TRUE,
                                                               node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 2, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names, 
                        " term - gene graph (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"))
  term_gene_graphs_ashr_ETR_vs_ETNR[[names]]$layers = list(term_gene_graphs_ashr_ETR_vs_ETNR[[names]]$layers[[1]], 
                                                           term_gene_graphs_ashr_ETR_vs_ETNR[[names]]$layers[[2]], 
                                                           term_gene_graphs_ashr_ETR_vs_ETNR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_ETR_vs_ETNR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.png"),
         path = paste0("Reports/png_pathways/ashr_ETR_vs_ETNR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_ETR_vs_ETNR = list()
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  # UpSet plot
  UpSet_plots_ashr_ETR_vs_ETNR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]],
                                                     genes_df = pathf_input_ashr_ETR_vs_ETNR,
                                                     num_terms = 5,
                                                     use_description = TRUE,
                                                     low = "darkgreen",
                                                     high = "darkred",
                                                     mid = "white")+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 3),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                              combmatrix.panel.point.size = 0.8,
                              combmatrix.panel.line.size = 0.5)+
    labs(fill = "logFC",
         title = paste0("Top 5 ", names, 
                        " term - UpSet plot (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"))
  tiff(paste0("Reports/png_pathways/ashr_ETR_vs_ETNR/", names, "/",
              names, "_top5_UpSet_plot.png"), 
       width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_ETR_vs_ETNR[[names]])
  dev.off()
}