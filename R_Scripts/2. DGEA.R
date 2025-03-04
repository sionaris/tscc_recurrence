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
library(gridExtra)

# Import environment
load("Resources_and_environments/Preprocessing_and_exploratory_analysis.RData")

# Differential Expression #####
# We are going to use a LRT as described in the DESeq2 vignette:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
# and the RNAseq workflow by Michael Love:
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
# in order to determine which interactions term add to the model

dds_lrt_3way = DESeq(dds, test = "LRT", reduced = ~ Recurrence + Stage + Sample.type +
                       Recurrence:Stage +
                       Recurrence:Sample.type +
                       Stage:Sample.type)
res_lrt_3way = results(dds_lrt_3way) %>% as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(padj)

# using pre-existing size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# final dispersion estimates
# fitting model and testing
# 43 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# -- replacing outliers and refitting for 19567 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
# 27 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]

# Only 40 significant genes. The three-way interaction will be discarded.
# Testing two-way interactions now with a new full model
dds_2way = DESeqDataSetFromMatrix(countData = assay(dds),
                                  colData = coldata,
                                  design = ~ Recurrence + Stage + Sample.type +
                                    Recurrence:Sample.type +
                                    Stage:Sample.type +
                                    Recurrence:Stage)

# Recurrence:Stage checks
dds_lrt_2way_a = DESeq(dds_2way, test = "LRT", reduced = ~ Recurrence + Stage + Sample.type +
                         Recurrence:Sample.type +
                         Stage:Sample.type)
res_lrt_2way_a = results(dds_lrt_2way_a) %>% as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(padj)

# estimating size factors
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# final dispersion estimates
# fitting model and testing
# 2 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# -- replacing outliers and refitting for 8731 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]

# 2327 significant genes. We keep it.

# Recurrence:Sample.type checks
dds_lrt_2way_b = DESeq(dds_2way, test = "LRT", reduced = ~ Recurrence + Stage + Sample.type +
                         Stage:Sample.type +
                         Recurrence:Stage)
res_lrt_2way_b = results(dds_lrt_2way_b) %>% as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(padj)

# estimating size factors
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# final dispersion estimates
# fitting model and testing
# 2 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# -- replacing outliers and refitting for 8731 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]

# Only 1 significant gene. We can discard this term
dds_2way_upd = DESeqDataSetFromMatrix(countData = assay(dds),
                                      colData = coldata,
                                      design = ~ Recurrence + Stage + Sample.type +
                                        Stage:Sample.type +
                                        Recurrence:Stage)

# Stage:Sample.type checks
dds_lrt_2way_c = DESeq(dds_2way_upd, test = "LRT", reduced = ~ Recurrence + Stage + Sample.type +
                         # Recurrence:Sample.type +
                         Recurrence:Stage)
res_lrt_2way_c = results(dds_lrt_2way_c) %>% as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(padj)

# estimating size factors
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# final dispersion estimates
# fitting model and testing
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# -- replacing outliers and refitting for 8812 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]

# 528 significant genes. We keep it

dds_Wald = DESeq(dds_2way_upd, test = "Wald")

# estimating size factors
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# final dispersion estimates
# fitting model and testing
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# -- replacing outliers and refitting for 8812 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]

# which(!mcols(dds_Wald)$betaConv)

# We discard these rows from dds_Wald (some of them may end up as DEGs otherwise)
dds_Wald_filt = dds_Wald[which(mcols(dds_Wald)$betaConv), ]
resultsNames(dds_Wald_filt) # Names of the coefficients

# Intercept: Baseline expression for the reference levels (Normal samples, Early stage, No recurrence).
# Sample.type_Tumor_vs_Normal: Effect of tumor samples compared to normal.
# Stage_Late_vs_Early: Effect of late stage compared to early.
# Recurrence_Recurrence_vs_No.recurrence: Effect of recurrence compared to no recurrence in normal samples.
# StageLate.Sample.typeTumor: Additional effect in tumor samples at late stage.
# RecurrenceRecurrence.StageLate: Additional effect of recurrence in late-stage samples.

# Create contrast list
contrasts <- list(
  # 1. Recurrent Tumors vs. Non-Recurrent Tumors Overall
  tumor_rec_vs_tumor_no_rec = list(
    contrast = c(
      "Recurrence_Recurrence_vs_No.recurrence",
      "RecurrenceRecurrence.StageLate"
    ),
    listValues = c(1, 0.5)
  ),
  
  # 2. Recurrent Late-Stage Tumors vs. Recurrent Early-Stage Tumors
  late_tumor_rec_vs_early_tumor_rec = list(
    contrast = c(
      "Stage_Late_vs_Early",
      "StageLate.Sample.typeTumor",
      "RecurrenceRecurrence.StageLate"
    ),
    listValues = c(1, 1, 1)
  ),
  
  # 3. Recurrent Early-Stage Tumors vs. Non-Recurrent Early-Stage Tumors
  early_tumor_rec_vs_early_tumor_no_rec = list(
    contrast = c(
      "Recurrence_Recurrence_vs_No.recurrence"
    ),
    listValues = c(1)
  ),
  
  # 4. Recurrent Late-Stage Tumors vs. Non-Recurrent Late-Stage Tumors
  late_tumor_rec_vs_late_tumor_no_rec = list(
    contrast = c(
      "Recurrence_Recurrence_vs_No.recurrence",
      "RecurrenceRecurrence.StageLate"
    ),
    listValues = c(1, 1)
  )
)

results <- list()
coefNames <- resultsNames(dds_Wald_filt)  # Get the names of the coefficients in the model

for (i in 1:length(contrasts)) {
  contrastNames <- contrasts[[i]][[1]]
  contrastValues <- contrasts[[i]][[2]]
  
  # Initialize a contrast vector with zeros
  contrastVector <- setNames(rep(0, length(coefNames)), coefNames)
  
  # Assign the specified values to the corresponding coefficients
  contrastVector[contrastNames] <- contrastValues
  
  # Pass the named numeric vector to the contrast argument
  results[[i]] <- results(dds_Wald_filt, contrast = contrastVector, test = "Wald")
}
names(results) <- paste0("res_", names(contrasts))
gc()

# Independent filtering
metadata(results[[1]])$filterThreshold
# 25.20408% of genes do not pass the threshold
# 1.207857

metadata(results[[2]])$filterThreshold
# 31.02041%% of genes do not pass the threshold
# 2.470962

metadata(results[[3]])$filterThreshold
# 31.02041%% of genes do not pass the threshold
# 2.470962

metadata(results[[4]])$filterThreshold
# 25.20408% of genes do not pass the threshold
# 1.207857

# Shrinking logFCs
# apeglm_shrunk_results = list() # apeglm does not work with contrasts
# normal_shrunk_results = list() # normal cannot be implemented with interactions designs
ashr_shrunk_results = list()

for (i in 1:length(results)) {
  # normal_shrunk_results[[i]] = lfcShrink(dds_Wald, contrast = contrasts[[i]], type = "normal")
  ashr_shrunk_results[[i]] = lfcShrink(dds_Wald_filt, res = results[[i]], type = "ashr")
  # ashr_shrunk_results[[i]]@listData$betaConv <- beta_conv
  # ashr_shrunk_results[[i]]@listData$Note <- ifelse(ashr_shrunk_results[[i]]@listData$betaConv == FALSE, 
  #                                                  "Non-converged", "")
}
names(ashr_shrunk_results) = names(results)

# Create a workbook and save data there to export later
dgea_wb <- createWorkbook()
dgea_wb_shrunk <- createWorkbook()

# Create header and body styles
header_style <- createStyle(textDecoration = "bold")
body_style <- createStyle(halign = "left", valign = "center")

# Set color_handle to TRUE to enable coloring
color_handle <- FALSE  # Change this to TRUE to enable colored symbols

if (color_handle) {
  blue_style <- createStyle(fontColour = "#279AE5", textDecoration = "bold")  # Blue and bold
  red_style <- createStyle(fontColour = "#B91B60", textDecoration = "bold")   # Red and bold
  grey_style <- createStyle(fontColour = "#818383", textDecoration = "bold")  # Grey and bold
  
  # Function to determine the color style based on conditions
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

# Create short sheet names
short_sheet_names <- c("TR_vs_TNR", "LTR_vs_ETR", "ETR_vs_ETNR", "LTR_vs_LTNR")
short_sheet_names_ashr <- paste0("ashr_", short_sheet_names)

# Loop through results
for (i in seq_along(results)) {
  # Raw LFCs
  genes <- rownames(results[[i]])
  res <- results[[i]]@listData %>%
    as.data.frame() %>%
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(Gene.Symbol, everything())
  
  # Add worksheet and write data
  addWorksheet(dgea_wb, short_sheet_names[i])
  writeData(dgea_wb, short_sheet_names[i], res)
  
  # Apply header style
  addStyle(dgea_wb, sheet = short_sheet_names[i], style = header_style,
           rows = 1, cols = 1:ncol(res), gridExpand = TRUE)
  
  # Apply body style
  addStyle(dgea_wb, sheet = short_sheet_names[i], style = body_style,
           rows = 2:(nrow(res) + 1), cols = 1:ncol(res), gridExpand = TRUE)
  
  if (color_handle) {
    # Generate color styles for each gene
    color_styles <- mapply(getStyle, res$log2FoldChange, res$padj, SIMPLIFY = FALSE)
    
    # Apply color styles to the 'Gene.Symbol' column
    for (j in seq_along(color_styles)) {
      addStyle(dgea_wb, sheet = short_sheet_names[i],
               style = color_styles[[j]],
               rows = j + 1,  # Adjust for header row
               cols = which(colnames(res) == "Gene.Symbol"),
               gridExpand = FALSE,
               stack = TRUE)
    }
  }
  
  # ASHR LFCs
  genes_ashr <- rownames(ashr_shrunk_results[[i]])
  res_ashr <- ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>%
    mutate(Gene.Symbol = genes_ashr) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(Gene.Symbol, everything())
  
  # Add worksheet and write data
  addWorksheet(dgea_wb_shrunk, short_sheet_names_ashr[i])
  writeData(dgea_wb_shrunk, short_sheet_names_ashr[i], res_ashr)
  
  # Apply header style
  addStyle(dgea_wb_shrunk, sheet = short_sheet_names_ashr[i], style = header_style,
           rows = 1, cols = 1:ncol(res_ashr), gridExpand = TRUE)
  
  # Apply body style
  addStyle(dgea_wb_shrunk, sheet = short_sheet_names_ashr[i], style = body_style,
           rows = 2:(nrow(res_ashr) + 1), cols = 1:ncol(res_ashr), gridExpand = TRUE)
  
  if (color_handle) {
    # Generate color styles for each gene
    color_styles_ashr <- mapply(getStyle, res_ashr$log2FoldChange, res_ashr$padj, SIMPLIFY = FALSE)
    
    # Apply color styles to the 'Gene.Symbol' column
    for (j in seq_along(color_styles_ashr)) {
      addStyle(dgea_wb_shrunk, sheet = short_sheet_names_ashr[i],
               style = color_styles_ashr[[j]],
               rows = j + 1,  # Adjust for header row
               cols = which(colnames(res_ashr) == "Gene.Symbol"),
               gridExpand = FALSE,
               stack = TRUE)
    }
  }
}

# Save workbooks
saveWorkbook(dgea_wb, "Results/DGEA/DGEA_raw.xlsx", overwrite = TRUE)
saveWorkbook(dgea_wb_shrunk, "Results/DGEA/DGEA_ashr.xlsx", overwrite = TRUE)

# Clean up
rm(genes, res, res_ashr, i, dgea_wb, dgea_wb_shrunk); gc()

# MA-plots #####
plot_ma <- function(data, title, show_legend=FALSE) {
  # Create a named vector for colors and labels
  color_values <- c("blue" = "blue", "grey" = "grey")
  color_labels <- c("blue" = "Significant", "grey" = "Not significant")
  
  # Start building the plot
  p <- ggplot(data, aes(x = mean_counts, y = log2FoldChange)) +
    geom_point(aes(color = Color), size = 0.5, alpha = 0.5)
  
  # Add scale_color_manual with or without guide based on show_legend
  if (show_legend) {
    p <- p + scale_color_manual(values = color_values,
                                labels = color_labels,
                                name = "Gene status")
  } else {
    p <- p + scale_color_manual(values = color_values,
                                labels = color_labels,
                                guide = "none")
  }
  
  # Continue building the plot
  p <- p +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    labs(title = title,
         x = expression(bold("Mean of Log"[2]*" (normalized counts + 1)")),
         y = expression(bold("Log"[2]*" fold change"))) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Position the legend on the right if show_legend is TRUE
  if (show_legend) {
    p <- p + theme(legend.position = "right")
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
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
  tiff(filename = paste0("Results/DGEA/MA_plots/MA_plot_", names(results)[i], ".tiff"),
       res = 300, width = 4000, height = 2000, compression = "lzw")
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
}

rm(res_df, res_shrunk_df, mean_norm_counts, common_genes, p1, p2); gc()

# Volcano plots #####
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
         path = "Results/DGEA/Volcano_plots", 
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
         path = "Results/DGEA/Volcano_plots", 
         width = 152, height = 142, device = 'tiff', units = "mm",
         dpi = 700, compression = "lzw")
}
rm(genes, i, res, keyvals.colour, volcano)

# Number of features for plots
n_degs = 100

# Top n_degs variable genes fromt the VST object
topVarGenes = head(order(rowVars(assay(vsd)), decreasing = TRUE), n_degs)
matvar  = assay(vsd)[topVarGenes, ]
matvar  = matvar - rowMeans(matvar)
colnames(matvar) = colnames(assay(vsd))
anno = as.data.frame(colData(vsd)[, c("Sample.type", "Stage", "Recurrence")])

# Create a custom palette with 20*14 = 280 color breaks (one for each tile)
colors2 = colorRampPalette(c(viridisLite::mako(280)[140:280],
                             rev(viridisLite::rocket(280)[140:280])))(280)

tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/VST_top", n_degs, "vargenes_samples_heatmap.tiff"),
     res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(matvar, annotation_col = anno,
         col = colors2,
         fontsize = 8,
         fontsize_row = 3,
         cutree_col = 2,
         main = paste0("Top ", n_degs, " variable VST genes (centered) heatmap"),
         annotation_colors = ann_colors,
         cluster_rows = FALSE)
dev.off()

# Filtered cohorts for filtered heatmaps
filtered_samples1 = Pheno$Patient.ID[Pheno$Sample.type == "Tumor"]
filtered_samples2 = Pheno$Patient.ID[Pheno$Sample.type == "Tumor" &
                                       Pheno$Recurrence == "Recurrence"]
filtered_samples3 = Pheno$Patient.ID[Pheno$Sample.type == "Tumor" &
                                       Pheno$Stage == "Early"]
filtered_samples4 = Pheno$Patient.ID[Pheno$Sample.type == "Tumor" &
                                       Pheno$Stage == "Late"]
filtered_samples = list(filtered_samples1, filtered_samples2,
                        filtered_samples3, filtered_samples4)
names(filtered_samples) = names(ashr_shrunk_results)
rm(filtered_samples1, filtered_samples2, filtered_samples3, filtered_samples4); gc()
columns_of_interest = c("Recurrence", "Stage", "Recurrence", "Recurrence")

# Heatmaps for raw LFCs
for (i in 1:length(results)) {
  # Construct the data frame
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj) %>%
    na.omit()
  
  topDEGs = res[which(res$padj < 0.05),]
  
  # Gene clustering with raw counts and Poisson distance
  topDEGs = topDEGs[1:n_degs, ]
  matcount = assay(dds_Wald_filt)[topDEGs$Gene.Symbol,]
  DEGpoisd = PoissonDistance(t(counts(dds_Wald_filt[topDEGs$Gene.Symbol,])) + .Machine$double.eps) # to avoid NA values
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
  tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/raw_lfcs/",
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
           fontsize_row = 4,
           fontsize_col = 6,
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Ordered genes and samples
  tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/raw_lfcs/",
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
           fontsize_row = 4,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Filtered versions
  tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/raw_lfcs/filtered_",
              names(results)[i], "_DEG_Poisson_clustered_samples_and_genes.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
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
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Ordered genes and samples
  column_of_interest = columns_of_interest[i]
  indices = colnames(logmat_ordered_znorm) %in% filtered_samples[[i]]
  arrangement = Pheno[Pheno$Patient.ID %in% colnames(logmat_ordered_znorm)[indices], ]
  arrangement = arrangement[order(arrangement[[column_of_interest]]), ]
  gaps_col = as.numeric(cumsum(table(arrangement[[column_of_interest]])))
  gaps_col = gaps_col[-length(gaps_col)]
  new_indices = arrangement$Patient.ID
  tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/raw_lfcs/filtered_",
              names(results)[i], "_DEG_samples_ordered.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
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
rm(genes, i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col, indices, arrangement,
   new_indices)

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
  tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
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
           fontsize_row = 4,
           fontsize_col = 6,
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  
  # Ordered genes and samples
  tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
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
           fontsize_row = 4,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Filtered versions
  tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/ashr_shrunk/filtered_",
              names(ashr_shrunk_results)[i], "_DEG_Poisson_clustered_samples_and_genes.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
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
  arrangement = Pheno[Pheno$Patient.ID %in% colnames(logmat_ordered_znorm)[indices], ]
  arrangement = arrangement[order(arrangement[[column_of_interest]]), ]
  gaps_col = as.numeric(cumsum(table(arrangement[[column_of_interest]])))
  gaps_col = gaps_col[-length(gaps_col)]
  new_indices = arrangement$Patient.ID
  tiff(paste0("Results/DGEA/Post_DGEA_heatmaps/ashr_shrunk/filtered_",
              names(results)[i], "_DEG_samples_ordered.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
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
rm(genes, i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col, indices, arrangement,
   new_indices)

gc()

# Export a session info text file for reproducibility
writeLines(capture.output(sessionInfo()), 
           paste0("Resources_and_environments/DGEA_sessionInfo.txt"))

# Export dds_Wald_filt object
saveRDS(dds_Wald_filt, "Resources_and_environments/final_dds.rds")

# Save environment for easy loading for DGEA later
save.image("Resources_and_environments/DGEA.RData")