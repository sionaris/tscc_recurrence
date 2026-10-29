# Exploratory plots #####

# Libraries
library(dplyr)
library(DESeq2)
library(tximport)
library(tximeta)
library(readr)
library(stringr)
library(edgeR)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(rcartocolor)
library(viridisLite)
library(vsn)
library(RColorBrewer)
library(tibble)
library(ggrepel)
library(glmpca)
library(openxlsx)
library(PoiClaClu)
library(mclust)
library(aricode)
library(EnhancedVolcano)
library(genefilter)
library(grid)

# Import environment
load("Resources_and_environments/Preprocessing_and_exploratory_analysis.RData")

png("Reports/png_exploratory/dimred_prior_to_sum.png", width = 2980, height = 1380, res = 700)
ggarrange(gpcaplot, poisson_mdsplot,
          ncol = 2, nrow = 1, common.legend = TRUE, 
          labels = c("A", "B"), legend = "right",
          font.label = list(size = 4, color = "black", face = "bold"))
dev.off()

png("Reports/png_exploratory/MDS.png", res = 700, width = 2880, height = 1080)
fig
annotate_figure(fig, left = textGrob(bquote(bold("MDS2")), rot = 90, vjust = 1, gp = gpar(cex = 0.45)))
dev.off()

fig = ggarrange(plot1, plot2, plot3, 
                ncol = 3, nrow = 1, 
                common.legend = TRUE, legend = "bottom",
                labels = c("A", "B", "C"),
                align = "hv",
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
png("Reports/png_exploratory/meanSdPlot.png", res = 300, width = 4500, height = 1600)
fig
annotate_figure(fig, left = textGrob(bquote(bold("sd")), rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob(bquote(bold("rank (mean)")), gp = gpar(cex = 1)))
dev.off()

create_density_curve(log2.norm.counts.plus1) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_log2normcounts_density_curves.png",
       path = "Reports/png_exploratory",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

create_density_plot_color(log2.norm.counts.plus1) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_log2normcounts_color_densities.png",
       path = "Reports/png_exploratory",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

create_density_curve(counts(dds)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_read_counts_density_curves.png",
       path = "Reports/png_exploratory",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

create_density_curve(counts(dds)) +
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 10)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_read_counts_density_curves_xlim30.png",
       path = "Reports/png_exploratory",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

create_density_plot_color(counts(dds)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_read_counts_color_densities.png",
       path = "Reports/png_exploratory",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

create_density_plot_color(counts(dds))+
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 10)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_read_counts_color_densities_xlim30.png",
       path = "Reports/png_exploratory",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

png("Reports/png_exploratory/PCA.png", width = 2880, height = 1380, res = 700)
ggarrange(pcaplot, gpcaplot, 
          ncol = 2, nrow = 1, common.legend = TRUE, 
          labels = c("A", "B"), legend = "right",
          font.label = list(size = 4, color = "black", face = "bold"))
dev.off()

png("Reports/png_exploratory/sample_distance_heatmap_vst.png", 
     res = 700, width = 6760, 
     height = 4760)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleVSTdists,
         clustering_distance_cols = sampleVSTdists,
         col = colors,
         fontsize = 8,
         main = "Sample distance heatmap",
         annotation_col = df,
         annotation_colors = ann_colors)
dev.off()

png("Reports/png_exploratory/sample_poisson_distance_heatmap.png", res = 700, width = 6760, 
     height = 4760)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         fontsize = 8,
         main = "Sample Poisson distance heatmap",
         annotation_col = df,
         annotation_colors = ann_colors)
dev.off()

png("Reports/png_exploratory/top50_heatmap_log2(norm.counts+1).png",
    res = 200, width = 1420, 
    height = 1080)
pheatmap(log2.norm.counts.plus1[select1, rownames(df)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "log2(norm.counts + 1)",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

png("Reports/png_exploratory/top50_heatmap_rlog.png",
    res = 200, width = 1420, 
    height = 1080)
pheatmap(assay(rld)[select2, rownames(df)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "rlog",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

png("Reports/png_exploratory/top50_heatmap_vst.png", res = 200, width = 1420, 
    height = 1080)
pheatmap(assay(vsd)[select3, rownames(df)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "vst",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

rm(list = ls()); gc()

# DGEA #####

# Libraries #####
library(dplyr)
library(DESeq2)
library(tximport)
library(tximeta)
library(readr)
library(stringr)
library(edgeR)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(rcartocolor)
library(viridisLite)
library(vsn)
library(RColorBrewer)
library(tibble)
library(ggrepel)
library(glmpca)
library(openxlsx)
library(PoiClaClu)
library(EnhancedVolcano)
library(genefilter)
library(grid)
library(gridExtra)

# Import environment
load("Resources_and_environments/DGEA.RData")

# Heatmaps for ashr-shrunk LFCs
for (i in 1:length(ashr_shrunk_results)) {
  # Construct the data frame
  genes = rownames(ashr_shrunk_results[[i]])
  res = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  topDEGs = res[which(res$padj < 0.05),]
  
  # Gene clustering with raw counts and Poisson distance
  topDEGs = topDEGs[1:n_degs, ]
  matcount = assay(dds_Wald_filt)[topDEGs$Gene.Symbol,]
  DEGpoisd = PoissonDistance(t(counts(dds_Wald_filt[topDEGs$Gene.Symbol,])) + .Machine$double.eps)
  DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
  rownames(DEGsamplePoisDistMatrix) = dds_Wald_filt$Patient.ID
  colnames(DEGsamplePoisDistMatrix) = dds_Wald_filt$Patient.ID
  
  # Same for genes
  DEGpoisd_genes = PoissonDistance(counts(dds_Wald_filt[topDEGs$Gene.Symbol,]) + .Machine$double.eps)
  DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
  rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  colnames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  
  # Log-transformed matrix
  logmat = log2(matcount + 1)
  
  # The same matrix but order samples by recurrence status and genes by logFC (raw)
  genes_order = topDEGs %>% dplyr::filter(padj < 0.05) %>%
    arrange(log2FoldChange)
  gaps_row = length(which(genes_order$log2FoldChange < 0))
  genes_order = genes_order$Gene.Symbol
  
  sample_order = anno %>% dplyr::arrange(Recurrence)
  gaps_col = length(which(sample_order$Recurrence == "No recurrence"))
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[genes_order, sample_order]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[genes_order, sample_order]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  png(paste0("Reports/png_dgea/",
              names(ashr_shrunk_results)[i], "_DEG_Poisson_clustered_samples_and_genes.png"),
       res = 700, width = 7260, 
       height = 5760)
  pheatmap(logmat_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = DEGpoisd$dd,
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 4,
           fontsize_col = 6,
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  
  # Ordered genes and samples
  png(paste0("Reports/png_dgea/",
              names(ashr_shrunk_results)[i], "_DEG_samples_ordered.png"),
       res = 700, width = 7260, 
       height = 5760)
  pheatmap(logmat_ordered_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           col = colors2,
           gaps_row = gaps_row,
           gaps_col = gaps_col,
           cutree_cols = 2,
           fontsize = 8,
           fontsize_row = 4,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Filtered versions
  png(paste0("Reports/png_dgea/filtered_",
              names(ashr_shrunk_results)[i], "_DEG_Poisson_clustered_samples_and_genes.png"),
       res = 700, width = 7260, 
       height = 5760)
  pheatmap(logmat_znorm[, filtered_samples[[i]]], annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = as.dist(DEGsamplePoisDistMatrix[filtered_samples[[i]], filtered_samples[[i]]]),
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 4,
           fontsize_col = 6,
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Ordered genes and samples
  column_of_interest = columns_of_interest[i]
  indices = colnames(logmat_ordered_znorm) %in% filtered_samples[[i]]
  arrangement = pheno_filt[pheno_filt$Patient.ID %in% colnames(logmat_ordered_znorm)[indices], ]
  arrangement = arrangement[order(arrangement[[column_of_interest]]), ]
  gaps_col = as.numeric(cumsum(table(arrangement[[column_of_interest]])))
  gaps_col = gaps_col[-length(gaps_col)]
  new_indices = arrangement$Patient.ID
  png(paste0("Reports/png_dgea/filtered_",
              names(results)[i], "_DEG_samples_ordered.png"),
       res = 700, width = 7260, 
       height = 5760)
  pheatmap(logmat_ordered_znorm[, new_indices], 
           annotation_col = anno[new_indices, ],
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           col = colors2,
           gaps_row = gaps_row,
           gaps_col = gaps_col,
           cutree_cols = 2,
           fontsize = 8,
           fontsize_row = 4,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}

for (i in 1:length(results)) {
  # Extract data from results
  res_df <- as.data.frame(results[[i]])
  
  # Get log2.norm.counts and compute mean normalized counts for each gene
  mean_norm_counts <- rowMeans(log2.norm.counts.plus1)
  
  # Ensure that the genes are matched correctly
  # Subset mean_norm_counts to include only genes in res_df
  common_genes <- intersect(rownames(res_df), names(mean_norm_counts))
  res_df <- res_df[common_genes, ]
  mean_norm_counts <- mean_norm_counts[common_genes]
  
  # Add mean_norm_counts to res_df
  res_df$mean_counts <- mean_norm_counts
  
  # Remove rows with NA in 'log2FoldChange' or 'mean_counts'
  res_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$mean_counts), ]
  
  # Create 'Color' column based on your logic
  res_df$Color <- ifelse(res_df$padj < 0.05, "blue", "grey")
  
  # Repeat similar steps for ashr_shrunk_results
  res_shrunk_df <- as.data.frame(ashr_shrunk_results[[i]])
  
  # Ensure genes are matched
  res_shrunk_df <- res_shrunk_df[common_genes, ]
  
  # Add mean_counts to res_shrunk_df
  res_shrunk_df$mean_counts <- mean_norm_counts
  
  # Remove rows with NA in 'log2FoldChange' or 'mean_counts'
  res_shrunk_df <- res_shrunk_df[!is.na(res_shrunk_df$log2FoldChange) & !is.na(res_shrunk_df$mean_counts), ]
  
  # Ensure 'Note' column matches the filtered res_df
  res_shrunk_df$Note <- res_df$Note[match(rownames(res_shrunk_df), rownames(res_df))]
  
  # Create 'Color' column for shrunk results
  res_shrunk_df$Color <- ifelse(res_shrunk_df$padj < 0.05, "blue", "grey")
  
  # Create MA plots
  p1 <- plot_ma(res_df, "MA plot: raw LFCs")
  p2 <- plot_ma(res_shrunk_df, "MA plot: ashr shrunk LFCs", show_legend = TRUE)
  
  # Save the plots side by side
  png(filename = paste0("Reports/png_dgea/MA_plot_", names(results)[i], ".png"),
       res = 300, width = 4000, height = 2000)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
}

for (i in 1:length(ashr_shrunk_results)) {
  # Construct the data frame
  genes = rownames(ashr_shrunk_results[[i]])
  res = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  # create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
  keyvals.colour <- ifelse(
    res$log2FoldChange < -1 & res$padj < 0.05, 'royalblue',
    ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red4',
           ifelse(abs(res$log2FoldChange) < 1 & res$padj < 0.05, 'pink',
                  'grey')))
  
  # keyvals.colour[is.na(keyvals.colour)] <- 'black'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
  names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
  names(keyvals.colour)[keyvals.colour == 'pink'] <- '|log2FC| < 1'
  names(keyvals.colour)[keyvals.colour == 'grey'] <- 'p.adj > 0.05'
  
  res$aes = keyvals.colour
  
  volcano = EnhancedVolcano(res,
                            lab = res[, "Gene.Symbol"],
                            caption = NULL,
                            x = 'log2FoldChange',
                            y = 'padj',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            cutoffLineType = "dashed",
                            cutoffLineWidth = 0.3,
                            cutoffLineCol = "black",
                            colCustom = keyvals.colour,
                            colAlpha = 0.7,
                            ylab = bquote(bold(-log[10]("BH adj. p-value"))),
                            xlab = bquote(bold(log[2]("Fold Change"))),
                            pointSize = 1,
                            axisLabSize = 7,
                            subtitle = NULL,
                            labSize = 2,
                            selectLab = res[1:20, "Gene.Symbol"],
                            legendLabSize = 8,
                            legendIconSize = 4,
                            title = paste0("Volcano plot: contrast = ",
                                           names(results)[i]),
                            labFace = "bold",
                            boxedLabels = TRUE,
                            drawConnectors = TRUE,
                            typeConnectors = "closed",
                            arrowheads = FALSE,
                            widthConnectors = 0.3,
                            max.overlaps = Inf) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 10, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(face = "bold", size = 8),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.ticks.length = unit(1, units = "mm"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "mm"),
          legend.spacing.y = unit(1, units = "mm"),
          legend.spacing.x = unit(1, units = "mm"),
          legend.background = element_blank(),
          legend.key.size = unit(0.5, units = "mm"),
          legend.text = element_text(size = 6),
          legend.box.background = element_rect(colour = "black"))
  
  volcano
  ggsave(filename = paste0("Volcano_ashr_shrunk_LFCs_", names(ashr_shrunk_results)[i], ".png"),
         path = "Reports/png_dgea", 
         width = 152, height = 142, device = 'png', units = "mm",
         dpi = 700)
}

rm(list = ls()); gc()

# Pathways #####

# Libraries
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
library(crayon)

# Import environment
load("Resources_and_environments/Pathways.RData")

# Create subdirectories
imgdirs = names(comparisons)

for (imgdir in imgdirs){
  if (!dir.exists(paste0("Reports/png_pathways/", imgdir))) {
    dir.create(paste0("Reports/png_pathways/", imgdir))
  }
}

subimgdirs = pathfindR_sheets_names

for (imgdir in imgdirs) {
  for (subimgdir in subimgdirs) {
    if (!dir.exists(paste0("Reports/png_pathways/", imgdir, "/", subimgdir))) {
      dir.create(paste0("Reports/png_pathways/", imgdir, "/", subimgdir))
    }
  }
}

nonzero_res = which(lapply(pathfindR_outputs_ashr_TR_vs_TNR, nrow) != 0)
# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  cat(blue("Trying enrichment dotplot for", names, "\n"))
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
    cat(blue("Trying clustered enrichment dotplot for", names, "\n"))
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

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  cat(blue("Trying term gene heatmap for", names, "\n"))
  # term-gene heatmaps
  if (nrow(wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[names]]) > 1) {
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
  }
  
  # term-gene graphs
  cat(blue("Trying term gene graph for", names, "\n"))
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

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  cat(blue("Trying UpSet plot for", names, "\n"))
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
  png(paste0("Reports/png_pathways/ashr_TR_vs_TNR/", names, "/",
              names, "_top5_UpSet_plot.png"), 
       width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_TR_vs_TNR[[names]])
  dev.off()
}

nonzero_res = which(lapply(pathfindR_outputs_ashr_LTR_vs_ETR, nrow) != 0)
# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  cat(blue("Trying enrichment dotplot for", names, "\n"))
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
    cat(blue("Trying clustered enrichment dotplot for", names, "\n"))
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

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  cat(blue("Trying term gene heatmap for", names, "\n"))
  # term-gene heatmaps
  if (nrow(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[names]]) > 1) {
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
  }
  
  # term-gene graphs
  cat(blue("Trying term gene graph for", names, "\n"))
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

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  cat(blue("Trying UpSet plot for", names, "\n"))
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
  png(paste0("Reports/png_pathways/ashr_LTR_vs_ETR/", names, "/",
             names, "_top5_UpSet_plot.png"), 
      width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_LTR_vs_ETR[[names]])
  dev.off()
}

nonzero_res = which(lapply(pathfindR_outputs_ashr_ETR_vs_ETNR, nrow) != 0)
# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  cat(blue("Trying enrichment dotplot for", names, "\n"))
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
    cat(blue("Trying clustered enrichment dotplot for", names, "\n"))
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

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  cat(blue("Trying term gene heatmap for", names, "\n"))
  # term-gene heatmaps
  if (nrow(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]]) > 1) {
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
  }
  
  # term-gene graphs
  cat(blue("Trying term gene graph for", names, "\n"))
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

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  cat(blue("Trying UpSet plot for", names, "\n"))
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
  png(paste0("Reports/png_pathways/ashr_ETR_vs_ETNR/", names, "/",
             names, "_top5_UpSet_plot.png"), 
      width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_ETR_vs_ETNR[[names]])
  dev.off()
}

nonzero_res = which(lapply(pathfindR_outputs_ashr_LTR_vs_LTNR, nrow) != 0)
# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  cat(blue("Trying enrichment dotplot for", names, "\n"))
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
    cat(blue("Trying clustered enrichment dotplot for", names, "\n"))
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

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  cat(blue("Trying term gene heatmap for", names, "\n"))
  # term-gene heatmaps
  if (nrow(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]]) > 1) {
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
  }
  
  # term-gene graphs
  cat(blue("Trying term gene graph for", names, "\n"))
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

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  cat(blue("Trying UpSet plot for", names, "\n"))
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
  png(paste0("Reports/png_pathways/ashr_LTR_vs_LTNR/", names, "/",
             names, "_top5_UpSet_plot.png"), 
      width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_LTR_vs_LTNR[[names]])
  dev.off()
}

nonzero_res = which(lapply(pathfindR_outputs_ashr_L_vs_E, nrow) != 0)
# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_L_vs_E)[i]
  cat(blue("Trying enrichment dotplot for", names, "\n"))
  # unclustered results
  enrichment_dotplots_ashr_L_vs_E[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_L_vs_E[[names]],
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
                        " terms enrichment dotplot - (", comparisons["ashr_L_vs_E"], ")"))
  print(enrichment_dotplots_ashr_L_vs_E[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.png"),
         path = paste0("Reports/png_pathways/ashr_L_vs_E/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_L_vs_E)) {
    cat(blue("Trying clustered enrichment dotplot for", names, "\n"))
    cluster_enrichment_dotplots_ashr_L_vs_E[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_L_vs_E[[names]][clustered_results_ashr_L_vs_E[[names]]$Status
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
                          " terms enrichment dotplot - (", comparisons["ashr_L_vs_E"], ")"))
    print(cluster_enrichment_dotplots_ashr_L_vs_E[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.png"),
           path = paste0("Reports/png_pathways/ashr_L_vs_E/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
}

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_L_vs_E)[i]
  cat(blue("Trying term gene heatmap for", names, "\n"))
  # term-gene heatmaps
  if (nrow(wrapped_pathfindR_outputs_ashr_L_vs_E[[names]]) > 1) {
    term_gene_heatmaps_ashr_L_vs_E[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_L_vs_E[[names]],
                                                                genes_df = pathf_input_ashr_L_vs_E,
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
                          comparisons["ashr_L_vs_E"], ")"),
           fill = expression(log[2] ~ "FoldChange"))
    #print(ggdraw(align_legend(term_gene_heatmaps_ashr_L_vs_E[[names]], hjust = 0.5)))
    print(term_gene_heatmaps_ashr_L_vs_E[[names]])
    ggsave(filename = paste0(names,
                             "_top10_term_gene_heatmap.png"),
           path = paste0("Reports/png_pathways/ashr_L_vs_E/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
  
  # term-gene graphs
  cat(blue("Trying term gene graph for", names, "\n"))
  term_gene_graphs_ashr_L_vs_E[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_L_vs_E[[names]],
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
                        comparisons["ashr_L_vs_E"], ")"))
  term_gene_graphs_ashr_L_vs_E[[names]]$layers = list(term_gene_graphs_ashr_L_vs_E[[names]]$layers[[1]], 
                                                      term_gene_graphs_ashr_L_vs_E[[names]]$layers[[2]], 
                                                      term_gene_graphs_ashr_L_vs_E[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_L_vs_E[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.png"),
         path = paste0("Reports/png_pathways/ashr_L_vs_E/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_L_vs_E)[i]
  cat(blue("Trying UpSet plot for", names, "\n"))
  # UpSet plot
  UpSet_plots_ashr_L_vs_E[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_L_vs_E[[names]],
                                                genes_df = pathf_input_ashr_L_vs_E,
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
                        comparisons["ashr_L_vs_E"], ")"))
  png(paste0("Reports/png_pathways/ashr_L_vs_E/", names, "/",
             names, "_top5_UpSet_plot.png"), 
      width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_L_vs_E[[names]])
  dev.off()
}

nonzero_res = which(lapply(pathfindR_outputs_ashr_LTNR_vs_ETNR, nrow) != 0)
# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTNR_vs_ETNR)[i]
  cat(blue("Trying enrichment dotplot for", names, "\n"))
  # unclustered results
  enrichment_dotplots_ashr_LTNR_vs_ETNR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_LTNR_vs_ETNR[[names]],
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
                        " terms enrichment dotplot - (", comparisons["ashr_LTNR_vs_ETNR"], ")"))
  print(enrichment_dotplots_ashr_LTNR_vs_ETNR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.png"),
         path = paste0("Reports/png_pathways/ashr_LTNR_vs_ETNR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_LTNR_vs_ETNR)) {
    cat(blue("Trying clustered enrichment dotplot for", names, "\n"))
    cluster_enrichment_dotplots_ashr_LTNR_vs_ETNR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_LTNR_vs_ETNR[[names]][clustered_results_ashr_LTNR_vs_ETNR[[names]]$Status
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
                          " terms enrichment dotplot - (", comparisons["ashr_LTNR_vs_ETNR"], ")"))
    print(cluster_enrichment_dotplots_ashr_LTNR_vs_ETNR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.png"),
           path = paste0("Reports/png_pathways/ashr_LTNR_vs_ETNR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
}

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTNR_vs_ETNR)[i]
  cat(blue("Trying term gene heatmap for", names, "\n"))
  # term-gene heatmaps
  if (nrow(wrapped_pathfindR_outputs_ashr_LTNR_vs_ETNR[[names]]) > 1) {
    term_gene_heatmaps_ashr_LTNR_vs_ETNR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_LTNR_vs_ETNR[[names]],
                                                                      genes_df = pathf_input_ashr_LTNR_vs_ETNR,
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
                          comparisons["ashr_LTNR_vs_ETNR"], ")"),
           fill = expression(log[2] ~ "FoldChange"))
    #print(ggdraw(align_legend(term_gene_heatmaps_ashr_LTNR_vs_ETNR[[names]], hjust = 0.5)))
    print(term_gene_heatmaps_ashr_LTNR_vs_ETNR[[names]])
    ggsave(filename = paste0(names,
                             "_top10_term_gene_heatmap.png"),
           path = paste0("Reports/png_pathways/ashr_LTNR_vs_ETNR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
  
  # term-gene graphs
  cat(blue("Trying term gene graph for", names, "\n"))
  term_gene_graphs_ashr_LTNR_vs_ETNR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_LTNR_vs_ETNR[[names]],
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
                        comparisons["ashr_LTNR_vs_ETNR"], ")"))
  term_gene_graphs_ashr_LTNR_vs_ETNR[[names]]$layers = list(term_gene_graphs_ashr_LTNR_vs_ETNR[[names]]$layers[[1]], 
                                                            term_gene_graphs_ashr_LTNR_vs_ETNR[[names]]$layers[[2]], 
                                                            term_gene_graphs_ashr_LTNR_vs_ETNR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_LTNR_vs_ETNR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.png"),
         path = paste0("Reports/png_pathways/ashr_LTNR_vs_ETNR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTNR_vs_ETNR)[i]
  cat(blue("Trying UpSet plot for", names, "\n"))
  # UpSet plot
  UpSet_plots_ashr_LTNR_vs_ETNR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_LTNR_vs_ETNR[[names]],
                                                      genes_df = pathf_input_ashr_LTNR_vs_ETNR,
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
                        comparisons["ashr_LTNR_vs_ETNR"], ")"))
  png(paste0("Reports/png_pathways/ashr_LTNR_vs_ETNR/", names, "/",
             names, "_top5_UpSet_plot.png"), 
      width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_LTNR_vs_ETNR[[names]])
  dev.off()
}

nonzero_res = which(lapply(pathfindR_outputs_ashr_LT_vs_ET, nrow) != 0)
# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LT_vs_ET)[i]
  cat(blue("Trying enrichment dotplot for", names, "\n"))
  # unclustered results
  enrichment_dotplots_ashr_LT_vs_ET[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_LT_vs_ET[[names]],
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
                        " terms enrichment dotplot - (", comparisons["ashr_LT_vs_ET"], ")"))
  print(enrichment_dotplots_ashr_LT_vs_ET[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.png"),
         path = paste0("Reports/png_pathways/ashr_LT_vs_ET/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'png', units = "px",
         dpi = 700)
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_LT_vs_ET)) {
    cat(blue("Trying clustered enrichment dotplot for", names, "\n"))
    cluster_enrichment_dotplots_ashr_LT_vs_ET[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_LT_vs_ET[[names]][clustered_results_ashr_LT_vs_ET[[names]]$Status
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
                          " terms enrichment dotplot - (", comparisons["ashr_LT_vs_ET"], ")"))
    print(cluster_enrichment_dotplots_ashr_LT_vs_ET[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.png"),
           path = paste0("Reports/png_pathways/ashr_LT_vs_ET/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
}

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LT_vs_ET)[i]
  cat(blue("Trying term gene heatmap for", names, "\n"))
  # term-gene heatmaps
  if (nrow(wrapped_pathfindR_outputs_ashr_LT_vs_ET[[names]]) > 1) {
    term_gene_heatmaps_ashr_LT_vs_ET[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_LT_vs_ET[[names]],
                                                                  genes_df = pathf_input_ashr_LT_vs_ET,
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
                          comparisons["ashr_LT_vs_ET"], ")"),
           fill = expression(log[2] ~ "FoldChange"))
    #print(ggdraw(align_legend(term_gene_heatmaps_ashr_LT_vs_ET[[names]], hjust = 0.5)))
    print(term_gene_heatmaps_ashr_LT_vs_ET[[names]])
    ggsave(filename = paste0(names,
                             "_top10_term_gene_heatmap.png"),
           path = paste0("Reports/png_pathways/ashr_LT_vs_ET/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'png', units = "px",
           dpi = 700)
    dev.off()
  }
  
  # term-gene graphs
  cat(blue("Trying term gene graph for", names, "\n"))
  term_gene_graphs_ashr_LT_vs_ET[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_LT_vs_ET[[names]],
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
                        comparisons["ashr_LT_vs_ET"], ")"))
  term_gene_graphs_ashr_LT_vs_ET[[names]]$layers = list(term_gene_graphs_ashr_LT_vs_ET[[names]]$layers[[1]], 
                                                        term_gene_graphs_ashr_LT_vs_ET[[names]]$layers[[2]], 
                                                        term_gene_graphs_ashr_LT_vs_ET[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_LT_vs_ET[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.png"),
         path = paste0("Reports/png_pathways/ashr_LT_vs_ET/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'png', units = "px",
         dpi = 700)
  dev.off()
}

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LT_vs_ET)[i]
  cat(blue("Trying UpSet plot for", names, "\n"))
  # UpSet plot
  UpSet_plots_ashr_LT_vs_ET[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_LT_vs_ET[[names]],
                                                  genes_df = pathf_input_ashr_LT_vs_ET,
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
                        comparisons["ashr_LT_vs_ET"], ")"))
  png(paste0("Reports/png_pathways/ashr_LT_vs_ET/", names, "/",
             names, "_top5_UpSet_plot.png"), 
      width = 2880, height = 1620*4, res = 700)
  print(UpSet_plots_ashr_LT_vs_ET[[names]])
  dev.off()
}

rm(list = ls()); gc()