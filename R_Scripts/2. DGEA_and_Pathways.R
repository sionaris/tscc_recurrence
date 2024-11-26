# In this script we carry out differential gene expression analysis

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

# Import environment
load("Resources_and_environments/Preprocessing_and_exploratory_analysis.RData")

# Differential Expression #####
# We are going to use a LRT as described in the DESeq2 vignette:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
# and the RNAseq workflow by Michael Love:
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments

dds = DESeq(dds)
resultsNames(dds) # Names of the coefficients

# Rename for easy contrast call
contrast_tumor_vs_normal = "Pheno.Sample.typeTumor"
contrast_late_vs_early = "Pheno.StageLate"
contrast_recurrence_vs_norec = "Pheno.RecurrenceRecurrence"
contrast_tumor_rec = "Pheno.Sample.typeTumor.Pheno.RecurrenceRecurrence"
contrast_late_rec = "Pheno.StageLate.Pheno.RecurrenceRecurrence"

# Create contrast list
contrasts = list(Rec_vs_norec = list(c(contrast_recurrence_vs_norec,
                                       contrast_tumor_rec)),
                 )

contrasts = list(HMBA24h_vs_Control48h = list("treatment_HMBA_vs_control",
                                              "timepoint_48h_vs_24h"),
                 HMBA48h_vs_Control48h = list(c("treatment_HMBA_vs_control", 
                                                "treatmentHMBA.timepoint48h")),
                 HMBA72h_vs_Control48h = list(c("treatment_HMBA_vs_control",
                                                "timepoint_72h_vs_24h",
                                                "treatmentHMBA.timepoint72h"),
                                              "timepoint_48h_vs_24h"),
                 HMBA48h_vs_HMBA24h = list(c("timepoint_48h_vs_24h",
                                             "treatmentHMBA.timepoint48h")),
                 HMBA72h_vs_HMBA24h = list(c("treatmentHMBA.timepoint72h",
                                             "timepoint_72h_vs_24h")),
                 HMBA72h_vs_HMBA48h = list(c("treatmentHMBA.timepoint72h",
                                             "timepoint_72h_vs_24h"), 
                                           c("timepoint_48h_vs_24h",
                                             "treatmentHMBA.timepoint48h")))

# hmba24h vs control48h
res_HMBA24h_vs_control48h <- results(dds, contrast=contrasts[["HMBA24h_vs_Control48h"]], 
                                     listValues=c(1, -1), test = "Wald")

# hmba48h vs control48h
res_HMBA48h_vs_control48h <- results(dds, 
                                     contrast=contrasts[["HMBA48h_vs_Control48h"]], 
                                     listValues=c(1, -1), test = "Wald")

# hmba72h vs control48h
res_HMBA72h_vs_control48h <- results(dds, contrast=contrasts[["HMBA72h_vs_Control48h"]],
                                     listValues=c(1, -1), test = "Wald")

# HMBA48h vs HMBA24h
res_HMBA48h_vs_HMBA24h <- results(dds, 
                                  contrast=contrasts[["HMBA48h_vs_HMBA24h"]],
                                  listValues=c(1, -1), test = "Wald")

# HMBA72h vs HMBA24h
res_HMBA72h_vs_HMBA24h <- results(dds, 
                                  contrast=contrasts[["HMBA72h_vs_HMBA24h"]], 
                                  listValues=c(1, -1), test = "Wald")

# HMBA72h vs HMBA48h
res_HMBA72h_vs_HMBA48h <- results(dds, contrast=contrasts[["HMBA72h_vs_HMBA48h"]],
                                  listValues=c(1, -1), test = "Wald")

results = list(res_HMBA24h_vs_control48h, res_HMBA48h_vs_control48h, 
               res_HMBA72h_vs_control48h, res_HMBA48h_vs_HMBA24h,
               res_HMBA72h_vs_HMBA24h, res_HMBA72h_vs_HMBA48h)
names(results) = names(contrasts)

# Shrinking logFCs
# apeglm_shrunk_results = list() # apeglm does not work with contrasts
# normal_shrunk_results = list() # normal cannot be implemented with interactions designs
ashr_shrunk_results = list()

for (i in 1:length(results)) {
  # apeglm_shrunk_results[[i]] = lfcShrink(dds, contrast = contrasts[[i]], type = "apeglm")
  # normal_shrunk_results[[i]] = lfcShrink(dds, contrast = contrasts[[i]], type = "normal")
  ashr_shrunk_results[[i]] = lfcShrink(dds, res = results[[i]], type = "ashr")
}
names(ashr_shrunk_results) = names(results)

# Create a workbook and save data there to export later
library(openxlsx)
dgea_wb = createWorkbook()
dgea_wb_shrunk = createWorkbook()

# Create header and body styles for the output Excel Sheet
header_style = createStyle(textDecoration = "bold")
body_style = createStyle(halign = "left", valign = "center")

# Choose if you want the colored gene symbols or not (color coding is time-consuming)
color_handle = FALSE # Change that to FALSE if you don't want colored symbols

if (color_handle) {
  blue_style = createStyle(fontColour = "#279AE5", textDecoration = "bold")  # Blue and bold
  red_style = createStyle(fontColour = "#B91B60", textDecoration = "bold")  # Red and bold
  grey_style = createStyle(fontColour = "#818383", textDecoration = "bold")  # Grey and bold
  
  # Define a function to determine the color style based on the conditions
  getStyle <- function(log2foldchange, padj) {
    if (is.na(log2foldchange) | is.na(padj)) {
      return(grey_style)
    } else if (log2foldchange < 0 & padj < 0.05) {
      return(blue_style)
    } else if (log2foldchange > 0 & padj < 0.05) {
      return(red_style)
    } else {
      return(grey_style)
    }
  }
}

# Create data frames that will be exported into .xlsx files
for (i in 1:length(results)) {
  # Raw LFCs
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(Gene.Symbol, everything())
  
  addWorksheet(dgea_wb, names(results)[i])
  writeData(dgea_wb, names(results)[i], res)
  
  # Apply the styles
  addStyle(dgea_wb, names(results)[i], style = header_style, rows = 1, 
           cols = 1:ncol(res), gridExpand = TRUE) # first row only
  addStyle(dgea_wb, names(results)[i], style = body_style, rows = 2:nrow(res), 
           cols = 1:ncol(res), gridExpand = TRUE)
  
  if (color_handle) {
    color_styles = mapply(getStyle, res$log2FoldChange, res$padj)
    mapply(addStyle, MoreArgs = list(wb = dgea_wb, sheet = names(results)[i], 
                                     cols = which(colnames(res) == "Gene.Symbol")), 
           style = color_styles, rows = 2:(nrow(res)+1))
  }
  
  # ASHR LFCs
  genes = rownames(ashr_shrunk_results[[i]])
  res_ashr = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(Gene.Symbol, everything())
  
  addWorksheet(dgea_wb_shrunk, names(ashr_shrunk_results)[i])
  writeData(dgea_wb_shrunk, names(ashr_shrunk_results)[i], res_ashr)
  
  # Apply the styles
  addStyle(dgea_wb_shrunk, names(ashr_shrunk_results)[i], style = header_style, 
           rows = 1, cols = 1:ncol(res_ashr), gridExpand = TRUE) # first row only
  addStyle(dgea_wb_shrunk, names(ashr_shrunk_results)[i], style = body_style, 
           rows = 2:nrow(res_ashr), cols = 1:ncol(res_ashr), gridExpand = TRUE)
  
  if (color_handle) {
    color_styles = mapply(getStyle, res_ashr$log2FoldChange, res_ashr$padj)
    mapply(addStyle, MoreArgs = list(wb = dgea_wb_shrunk, 
                                     sheet = names(ashr_shrunk_results)[i], 
                                     cols = which(colnames(res_ashr) == "Gene.Symbol")), 
           style = color_styles, rows = 2:(nrow(res_ashr)+1))
  }
}

# Create DGEA folder
if (!dir.exists("DGEA")) {
  dir.create("DGEA")
}

saveWorkbook(dgea_wb, "DGEA/DGEA_raw.xlsx",
             overwrite = TRUE)
saveWorkbook(dgea_wb_shrunk, "DGEA/DGEA_ashr.xlsx",
             overwrite = TRUE)
rm(genes, res, res_ashr, i); gc()

# MA-plots
if (!dir.exists("DGEA/MA_plots")) {
  dir.create("DGEA/MA_plots")
}

for (i in 1:length(results)) {
  tiff(paste0("DGEA/MA_plots/MA_plot_",
              names(results)[i], ".tiff"), res = 700,
       width = 8000, height = 4000, compression = "lzw")
  par(mfrow=c(1,2), mar=c(4,4,2,1))
  print(DESeq2::plotMA(results[[i]], main = "MA plot: raw LFCs",
                       xlab = substitute(paste(bold('Mean of normalized counts'))),
                       ylab = substitute(paste(bold('log fold change'))),
                       alpha = 0.05))
  print(DESeq2::plotMA(ashr_shrunk_results[[i]], main = "MA plot: ashr shrunk LFCs",
                       xlab = substitute(paste(bold('Mean of normalized counts'))),
                       ylab = substitute(paste(bold('log fold change'))),
                       alpha = 0.05))
  dev.off()
}

# Volcano plots
library(EnhancedVolcano)

if (!dir.exists("DGEA/Volcano_plots")) {
  dir.create("DGEA/Volcano_plots")
}

# Volcano loop (raw LFCs)
for (i in 1:length(results)) {
  # Construct the data frame
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
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
  ggsave(filename = paste0("Volcano_raw_LFCs_", names(results)[i], ".tiff"),
         path = "DGEA/Volcano_plots", 
         width = 152, height = 142, device = 'tiff', units = "mm",
         dpi = 700, compression = "lzw")
}
rm(genes, i, res, keyvals.colour, volcano)

# Volcano loop (ashr-shrunk LFCs)
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
  ggsave(filename = paste0("Volcano_ashr_shrunk_LFCs_", names(ashr_shrunk_results)[i], ".tiff"),
         path = "DGEA/Volcano_plots", 
         width = 152, height = 142, device = 'tiff', units = "mm",
         dpi = 700, compression = "lzw")
}
rm(genes, i, res, keyvals.colour, volcano)

# Gene clustering with mean-centered values from the VST data
library(genefilter)

# Create output directories
if (!dir.exists("DGEA/Post_DGEA_heatmaps")) {
  dir.create("DGEA/Post_DGEA_heatmaps")
}

# Subdirectories
if (!dir.exists("DGEA/Post_DGEA_heatmaps/raw_lfcs")) {
  dir.create("DGEA/Post_DGEA_heatmaps/raw_lfcs")
}

if (!dir.exists("DGEA/Post_DGEA_heatmaps/ashr_shrunk")) {
  dir.create("DGEA/Post_DGEA_heatmaps/ashr_shrunk")
}

# Number of features for plots
n_degs = 20

# Top n_degs variable genes fromt the VST object
topVarGenes = head(order(rowVars(assay(vsd)), decreasing = TRUE), n_degs)
matvar  = assay(vsd)[topVarGenes, ]
matvar  = matvar - rowMeans(matvar)
colnames(matvar) = colnames(assay(vsd))
anno = as.data.frame(colData(vsd)[, c("batch", "treatment", "timepoint")])

# Create a custom palette with 20*14 = 280 color breaks (one for each tile)
colors2 = colorRampPalette(c(viridisLite::mako(280)[140:280],
                             rev(viridisLite::rocket(280)[140:280])))(280)

tiff(paste0("DGEA/Post_DGEA_heatmaps/VST_top20vargenes_samples_heatmap.tiff"),
     res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(matvar, annotation_col = anno,
         col = colors2,
         fontsize = 8,
         cutree_col = 2,
         main = "Top 20 variable VST genes (centered) heatmap",
         annotation_colors = ann_colors,
         cluster_rows = FALSE)
dev.off()

# Heatmaps for raw LFCs
for (i in 1:length(results)) {
  # Construct the data frame
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  topDEGs = res[which(res$padj < 0.05),]
  
  # Gene clustering with raw counts and Poisson distance
  topDEGs = topDEGs[1:n_degs, ]
  matcount = assay(dds)[topDEGs$Gene.Symbol,]
  DEGpoisd = PoissonDistance(t(counts(dds[topDEGs$Gene.Symbol,])))
  DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
  rownames(DEGsamplePoisDistMatrix) = dds$cell_line
  colnames(DEGsamplePoisDistMatrix) = NULL
  
  # Same for genes
  DEGpoisd_genes = PoissonDistance(counts(dds[topDEGs$Gene.Symbol,]))
  DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
  rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  colnames(DEGgenesPoisDistMatrix) = NULL
  
  # Log-transformed matrix
  logmat = log2(matcount + 1)
  
  # The same matrix but order samples by metastatic status and genes by logFC (raw)
  genes_order = topDEGs %>% dplyr::filter(padj < 0.05) %>%
    arrange(log2FoldChange)
  gaps_row = length(which(genes_order$log2FoldChange < 0))
  genes_order = genes_order$Gene.Symbol
  
  sample_order = anno %>% dplyr::arrange(treatment)
  gaps_col = length(which(sample_order$treatment == "control"))
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[genes_order, sample_order]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[genes_order, sample_order]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  tiff(paste0("DGEA/Post_DGEA_heatmaps/raw_lfcs/",
              names(results)[i], "_DEG_Poisson_clustered_samples_and_genes.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = DEGpoisd$dd,
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  
  # Ordered genes and samples
  tiff(paste0("DGEA/Post_DGEA_heatmaps/raw_lfcs/",
              names(results)[i], "_DEG_samples_ordered.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
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
           fontsize_row = 8,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}
rm(genes, i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col)

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
  matcount = assay(dds)[topDEGs$Gene.Symbol,]
  DEGpoisd = PoissonDistance(t(counts(dds[topDEGs$Gene.Symbol,])))
  DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
  rownames(DEGsamplePoisDistMatrix) = dds$cell_line
  colnames(DEGsamplePoisDistMatrix) = NULL
  
  # Same for genes
  DEGpoisd_genes = PoissonDistance(counts(dds[topDEGs$Gene.Symbol,]))
  DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
  rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  colnames(DEGgenesPoisDistMatrix) = NULL
  
  # Log-transformed matrix
  logmat = log2(matcount + 1)
  
  # The same matrix but order samples by metastatic status and genes by logFC (raw)
  genes_order = topDEGs %>% dplyr::filter(padj < 0.05) %>%
    arrange(log2FoldChange)
  gaps_row = length(which(genes_order$log2FoldChange < 0))
  genes_order = genes_order$Gene.Symbol
  
  sample_order = anno %>% dplyr::arrange(treatment)
  gaps_col = length(which(sample_order$treatment == "control"))
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[genes_order, sample_order]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[genes_order, sample_order]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  tiff(paste0("DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
              names(ashr_shrunk_results)[i], "_DEG_Poisson_clustered_samples_and_genes.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = DEGpoisd$dd,
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  
  # Ordered genes and samples
  tiff(paste0("DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
              names(ashr_shrunk_results)[i], "_DEG_samples_ordered.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
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
           fontsize_row = 8,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}

rm(genes, i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col)

gc()

# Session info
sessionInfo()