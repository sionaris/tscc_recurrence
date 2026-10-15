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
library(IHW)
library(ashr)
library(patchwork)

# Import environment
load("Resources_and_environments/Preprocessing_and_exploratory_analysis.RData")

# Differential Expression #####
# Tumor-only model for contrasts of interest
dds_Wald_filt <- dds[, dds$Sample.type == "Tumor"]
dds_Wald_filt$Stage      <- droplevels(dds_Wald_filt$Stage)
dds_Wald_filt$Recurrence <- droplevels(dds_Wald_filt$Recurrence)

dds_Wald_filt <- estimateSizeFactors(dds_Wald_filt)     # re-estimate on subset
design(dds_Wald_filt) <- ~ Stage * Recurrence
dds_Wald_filt <- DESeq(dds_Wald_filt, test = "Wald")

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
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]
# -- replacing outliers and refitting for 9295 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing
# Note: levels of factors in the design contain characters other than
# letters, numbers, '_' and '.'. It is recommended (but not required) to use
# only letters, numbers, and delimiters '_' or '.', as these are safe characters
# for column names in R. [This is a message, not a warning or an error]

coefNames <- resultsNames(dds_Wald_filt)
coefNames

# [1] "Intercept"                              "Stage_Late_vs_Early"                   
# [3] "Recurrence_Recurrence_vs_No.recurrence" "StageLate.RecurrenceRecurrence"   

# Define contrasts (use the actual interaction name "StageLate.RecurrenceRecurrence")
contrasts <- list(
  tumor_rec_vs_tumor_no_rec = list(
    contrast   = c("Recurrence_Recurrence_vs_No.recurrence", "StageLate.RecurrenceRecurrence"),
    listValues = c(1, 0.5)   # TR vs TNR, equal stage weighting
  ),
  
  late_tumor_rec_vs_early_tumor_rec = list(
    contrast   = c("Stage_Late_vs_Early", "StageLate.RecurrenceRecurrence"),
    listValues = c(1, 1)     # LTR vs ETR
  ),
  
  early_tumor_rec_vs_early_tumor_no_rec = list(
    contrast   = c("Recurrence_Recurrence_vs_No.recurrence"),
    listValues = c(1)        # ETR vs ETNR
  ),
  
  late_tumor_rec_vs_late_tumor_no_rec = list(
    contrast   = c("Recurrence_Recurrence_vs_No.recurrence", "StageLate.RecurrenceRecurrence"),
    listValues = c(1, 1)     # LTR vs LTNR
  )
)

# Get results and perform IHW
results_raw  <- list()   # DESeqResults objects
results_tbls <- list()   # tidy data.frames with padj_ihw
ihw_objs     <- list()   # IHW fit objects (optional to keep)

for (nm in names(contrasts)) {
  cn   <- contrasts[[nm]]$contrast
  vals <- contrasts[[nm]]$listValues
  
  # Build a named numeric contrast vector over all coefficients
  cv <- setNames(rep(0, length(coefNames)), coefNames)
  
  # Sanity check: all requested coef names exist
  miss <- setdiff(cn, names(cv))
  if (length(miss) > 0) stop(sprintf("Contrast '%s' has unknown coefficients: %s",
                                     nm, paste(miss, collapse = ", ")))
  
  cv[cn] <- vals
  
  # Get DESeq2 results WITHOUT independent filtering,
  # since we will do IHW ourselves on pvalue ~ baseMean
  res <- results(dds_Wald_filt,
                 contrast = cv,
                 independentFiltering = FALSE,
                 cooksCutoff = FALSE)
  
  # Tidy table
  deRes <- as.data.frame(res)
  deRes$gene <- rownames(deRes)
  deRes <- deRes[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
  
  # IHW on (pvalue ~ baseMean). Remove NAs first.
  keep <- is.finite(deRes$pvalue) & !is.na(deRes$baseMean)
  ihwRes <- ihw(pvalue ~ baseMean, data = deRes[keep, c("pvalue","baseMean")], alpha = 0.05)
  
  deRes$padj_ihw <- NA_real_
  deRes$padj_ihw[keep] <- adj_pvalues(ihwRes)
  
  # Store
  results_raw[[paste0("res_", nm)]]  <- res
  results_tbls[[paste0("res_", nm)]] <- deRes
  ihw_objs[[paste0("ihw_", nm)]]     <- ihwRes
}

# Rejected null hypotheses in IHW results
lapply(ihw_objs, rejections)
gc()

# IHW diagnostics #####

# Create output dirs for diagnostics
diag_root <- "Results/DGEA/IHW_diagnostics"
if (!dir.exists(diag_root)) dir.create(diag_root, recursive = TRUE)

for (nm in names(ihw_objs)) {
  ihw_obj <- ihw_objs[[nm]]
  res_name <- sub("^ihw_", "res_", nm)
  deRes <- results_tbls[[res_name]]
  out_dir <- file.path(diag_root, nm)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 1) Official IHW diagnostic plot (weights vs covariate strata, stratified p-value histograms)
  pdf_file <- file.path(out_dir, paste0(nm, "_ihw_diagnostic_plot.pdf"))
  obj <- plot(ihw_objs[[nm]])  # build plots
  
  grDevices::pdf(pdf_file, width = 8, height = 5, onefile = TRUE)
  if (is.list(obj)) {
    for (g in obj) {
      if (inherits(g, "ggplot")) {
        print(g)
      } else {
        grid::grid.newpage(); grid::grid.draw(g)
      }
    }
  } else {
    if (inherits(obj, "ggplot")) {
      print(obj)
    } else {
      grid::grid.newpage(); grid::grid.draw(obj)
    }
  }
  dev.off()
  
  # 2) Stratified p-value histograms by baseMean quintiles (null-independence visual check)
  dfp <- deRes[is.finite(deRes$pvalue) & !is.na(deRes$baseMean), c("pvalue","baseMean")]
  if (nrow(dfp) > 0) {
    dfp$bin <- cut(rank(dfp$baseMean)/nrow(dfp),
                   breaks = seq(0, 1, by = 0.2),
                   include.lowest = TRUE,
                   labels = paste0("Q", 1:5))
    p_hist <- ggplot2::ggplot(dfp, aes(x = pvalue)) +
      geom_histogram(bins = 40, color = "grey20", fill = "grey70") +
      facet_wrap(~bin, nrow = 1) +
      labs(title = paste0(nm, " — p-value histograms by baseMean quintiles"),
           x = "p-value", y = "Count") +
      theme_minimal(base_size = 10)
    ggplot2::ggsave(filename = file.path(out_dir, paste0(nm, "_pval_hist_by_quintile.png")),
                    plot = p_hist, width = 260, height = 70, units = "mm",
                    dpi = 300)
  }
  
  # 3) Rejection summaries (BH vs IHW) overall and by baseMean quintile
  deRes$padj_bh <- p.adjust(deRes$pvalue, method = "BH")
  deRes$ihw_reject <- !is.na(deRes$padj_ihw) & deRes$padj_ihw < 0.05
  deRes$bh_reject  <- !is.na(deRes$padj_bh)  & deRes$padj_bh  < 0.05
  
  # overall counts
  overall <- data.frame(
    contrast   = nm,
    bh_rejects = sum(deRes$bh_reject,  na.rm = TRUE),
    ihw_rejects= sum(deRes$ihw_reject, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # by quintile
  dfq <- deRes[is.finite(deRes$pvalue) & !is.na(deRes$baseMean), ]
  if (nrow(dfq) > 0) {
    dfq$bin <- cut(rank(dfq$baseMean)/nrow(dfq),
                   breaks = seq(0, 1, by = 0.2),
                   include.lowest = TRUE,
                   labels = paste0("Q", 1:5))
    byq <- dfq %>%
      dplyr::group_by(bin) %>%
      dplyr::summarise(
        n = dplyr::n(),
        bh_rejects  = sum(bh_reject,  na.rm = TRUE),
        ihw_rejects = sum(ihw_reject, na.rm = TRUE),
        prop_bh_rejects  = bh_rejects / n,
        prop_ihw_rejects = ihw_rejects / n,
        # uniformity sanity check for large p-values within strata
        prop_p_gt_0.5 = mean(pvalue > 0.5, na.rm = TRUE),
        .groups = "drop"
      )
    readr::write_csv(overall, file.path(out_dir, paste0(nm, "_overall_rejections.csv")))
    readr::write_csv(byq,     file.path(out_dir, paste0(nm, "_rejections_by_quintile.csv")))
  } else {
    readr::write_csv(overall, file.path(out_dir, paste0(nm, "_overall_rejections.csv")))
  }
  
  # 4) Quick console summary for go/no-go
  cat("\n========== IHW diagnostics for", nm, "==========\n")
  cat("BH rejects:", overall$bh_rejects, " | IHW rejects:", overall$ihw_rejects, "\n")
  if (exists("byq")) {
    print(byq)
    cat("Heuristic check: prop_p_gt_0.5 should be ~0.5 in each bin under null-uniformity.\n")
    rm(byq)
  }
}

rm(diag_root, out_dir, ihw_obj, res_name, deRes, dfp, dfq, overall, p_hist, nm)
gc()

# Shrinking logFCs #####
# apeglm_shrunk_results = list() # apeglm does not work with contrasts
# normal_shrunk_results = list() # normal cannot be implemented with interactions designs
ashr_shrunk_results_raw <- list()  # DESeqResults with shrunken LFCs
ashr_shrunk_tbls       <- list()  # data.frames with log2FC_shrunk + padj_ihw retained

for (nm in names(contrasts)) {
  res_nm <- paste0("res_", nm)
  
  # 1) Shrink using the *same* fitted contrast by passing the DESeqResults
  shr <- lfcShrink(dds_Wald_filt, res = results_raw[[res_nm]], type = "ashr")
  
  # 2) Tidy + merge with your IHW table so we preserve pvalue/padj/padj_ihw
  shr_df <- as.data.frame(shr)
  shr_df$gene <- rownames(shr_df)
  shr_df <- shr_df[, c("gene","baseMean","log2FoldChange","lfcSE")]
  
  # keep original names and add a clear column for shrunk LFC
  colnames(shr_df)[colnames(shr_df) == "log2FoldChange"] <- "log2FC_shrunk"
  colnames(shr_df)[colnames(shr_df) == "lfcSE"]          <- "lfcSE_shrunk"
  
  # merge with the IHW-augmented table you already made
  merged <- merge(
    results_tbls[[res_nm]],  # has gene, baseMean, log2FoldChange (MLE), lfcSE, stat, pvalue, padj, padj_ihw
    shr_df[, c("gene","log2FC_shrunk","lfcSE_shrunk")],
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  
  # store
  ashr_shrunk_results_raw[[res_nm]] <- shr
  ashr_shrunk_tbls[[res_nm]]        <- merged
}

## ---- Workbooks, styles ----
dgea_wb        <- createWorkbook()
dgea_wb_shrunk <- createWorkbook()

header_style <- createStyle(textDecoration = "bold")
body_style   <- createStyle(halign = "left", valign = "center")

color_handle <- FALSE  # set TRUE to enable coloured gene symbols

if (color_handle) {
  blue_style <- createStyle(fontColour = "#279AE5", textDecoration = "bold")  # down
  red_style  <- createStyle(fontColour = "#B91B60", textDecoration = "bold")  # up
  grey_style <- createStyle(fontColour = "#818383", textDecoration = "bold")  # ns
  
  getStyle <- function(lfc, padj) {
    if (is.na(lfc) || is.na(padj))        return(grey_style)
    if (lfc < 0  && padj < 0.05)          return(blue_style)
    if (lfc > 0  && padj < 0.05)          return(red_style)
    return(grey_style)
  }
}

## ---- Map contrast names -> short sheet names ----
short_names <- c(
  tumor_rec_vs_tumor_no_rec             = "TR_vs_TNR",
  late_tumor_rec_vs_early_tumor_rec     = "LTR_vs_ETR",
  early_tumor_rec_vs_early_tumor_no_rec = "ETR_vs_ETNR",
  late_tumor_rec_vs_late_tumor_no_rec   = "LTR_vs_LTNR"
)

## ---- Helper to write one sheet ----
.write_sheet <- function(wb, sheet, df, lfc_col, padj_col, color_handle=FALSE) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df)
  addStyle(wb, sheet=sheet, style=header_style, rows=1, cols=1:ncol(df), gridExpand=TRUE)
  if (nrow(df) > 0) {
    addStyle(wb, sheet=sheet, style=body_style, rows=2:(nrow(df)+1), cols=1:ncol(df), gridExpand=TRUE)
  }
  if (color_handle && nrow(df) > 0) {
    color_styles <- mapply(getStyle, df[[lfc_col]], df[[padj_col]], SIMPLIFY = FALSE)
    gs_col <- which(colnames(df) == "Gene.Symbol")
    for (j in seq_along(color_styles)) {
      addStyle(wb, sheet=sheet, style=color_styles[[j]], rows=j+1, cols=gs_col, gridExpand=FALSE, stack=TRUE)
    }
  }
}

## ---- RAW workbook: results_tbls (uses padj_ihw if present) ----
for (nm in names(short_names)) {
  res_nm <- paste0("res_", nm)
  if (!res_nm %in% names(results_tbls)) next
  
  df <- results_tbls[[res_nm]]
  # ensure a consistent column set and sorting
  if (!"gene" %in% names(df)) df$gene <- rownames(df)
  df$Gene.Symbol <- df$gene
  
  # pick adj p column: prefer padj_ihw, else padj
  padj_col <- if ("padj_ihw" %in% names(df)) "padj_ihw" else "padj"
  lfc_col  <- "log2FoldChange"
  
  keep_cols <- c("Gene.Symbol","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj", 
                 if ("padj_ihw" %in% names(df)) "padj_ihw" else NULL)
  keep_cols <- intersect(keep_cols, names(df))
  
  df_out <- df[, keep_cols]
  # sort by chosen adjusted p
  ord <- order(ifelse(is.na(df_out[[padj_col]]), Inf, df_out[[padj_col]]), na.last = TRUE)
  df_out <- df_out[ord, , drop=FALSE]
  
  .write_sheet(dgea_wb, short_names[[nm]], df_out, lfc_col=lfc_col, padj_col=padj_col, color_handle=color_handle)
}

## ---- ASHR-shrunk workbook: ashr_shrunk_tbls (includes shrunk LFCs) ----
for (nm in names(short_names)) {
  res_nm <- paste0("res_", nm)
  if (!res_nm %in% names(ashr_shrunk_tbls)) next
  
  df <- ashr_shrunk_tbls[[res_nm]]
  if (!"gene" %in% names(df)) df$gene <- rownames(df)
  df$Gene.Symbol <- df$gene
  
  padj_col <- if ("padj_ihw" %in% names(df)) "padj_ihw" else "padj"
  lfc_col  <- if ("log2FC_shrunk" %in% names(df)) "log2FC_shrunk" else "log2FoldChange"
  
  keep_cols <- c("Gene.Symbol","baseMean",
                 # show both MLE and shrunk if available
                 if ("log2FoldChange" %in% names(df)) "log2FoldChange" else NULL,
                 if ("log2FC_shrunk" %in% names(df)) "log2FC_shrunk" else NULL,
                 if ("lfcSE" %in% names(df)) "lfcSE" else NULL,
                 if ("lfcSE_shrunk" %in% names(df)) "lfcSE_shrunk" else NULL,
                 "stat","pvalue","padj",
                 if ("padj_ihw" %in% names(df)) "padj_ihw" else NULL)
  keep_cols <- intersect(keep_cols, names(df))
  
  df_out <- df[, keep_cols]
  ord <- order(ifelse(is.na(df_out[[padj_col]]), Inf, df_out[[padj_col]]), na.last = TRUE)
  df_out <- df_out[ord, , drop=FALSE]
  
  .write_sheet(dgea_wb_shrunk, paste0("ashr_", short_names[[nm]]), df_out,
               lfc_col=lfc_col, padj_col=padj_col, color_handle=color_handle)
}

## ---- Save ----
saveWorkbook(dgea_wb,        "Results/DGEA/DGEA_raw.xlsx",   overwrite = TRUE)
saveWorkbook(dgea_wb_shrunk, "Results/DGEA/DGEA_ashr.xlsx",  overwrite = TRUE)

# Clean up
rm(dgea_wb, dgea_wb_shrunk, df, df_out, res_nm, nm, ord, short_names); gc()

# MA-plots #####
plot_ma <- function(data, title, show_legend=FALSE) {
  color_values <- c("blue" = "blue", "grey" = "grey")
  color_labels <- c("blue" = "Significant", "grey" = "Not significant")
  
  p <- ggplot(data, aes(x = mean_counts, y = log2FoldChange)) +
    geom_point(aes(color = Color), size = 0.5, alpha = 0.5)
  
  if (show_legend) {
    p <- p + scale_color_manual(values = color_values,
                                labels = color_labels,
                                name = "Gene status")
  } else {
    p <- p + scale_color_manual(values = color_values,
                                labels = color_labels,
                                guide = "none")
  }
  
  p <- p +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    labs(title = title,
         x = expression(bold("Mean of Log"[2]*" (normalized counts + 1)")),
         y = expression(bold("Log"[2]*" fold change"))) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(legend.position = if (show_legend) "right" else "none")
  
  return(p)
}

# objects adapted to new pipeline
plot_list_raw    <- results_tbls
plot_list_shrunk <- ashr_shrunk_tbls

for (i in 1:length(plot_list_raw)) {
  # Construct raw results data.frame
  res_df <- as.data.frame(plot_list_raw[[i]])
  if ("gene" %in% colnames(res_df)) rownames(res_df) <- res_df$gene
  if ("Gene.Symbol" %in% colnames(res_df)) rownames(res_df) <- res_df$Gene.Symbol
  
  # Mean of log2(normalized counts + 1)
  mean_norm_counts <- rowMeans(log2.norm.counts.plus1)
  
  # Align genes
  common_genes <- intersect(rownames(res_df), names(mean_norm_counts))
  res_df <- res_df[common_genes, , drop = FALSE]
  mean_norm_counts <- mean_norm_counts[common_genes]
  
  # Add mean_counts and p-value column to use (prefer IHW)
  res_df$mean_counts <- mean_norm_counts
  res_df$padj_use <- if ("padj_ihw" %in% colnames(res_df) && any(!is.na(res_df$padj_ihw))) res_df$padj_ihw else res_df$padj
  
  # Keep rows with finite LFC, mean_counts, and padj_use
  res_df <- res_df[is.finite(res_df$log2FoldChange) & is.finite(res_df$mean_counts), , drop = FALSE]
  
  # Color by significance (avoid name-matching that produced NAs)
  res_df$Color <- ifelse(!is.na(res_df$padj_use) & res_df$padj_use < 0.05, "blue", "grey")
  res_df$Color <- factor(res_df$Color, levels = c("grey", "blue"))
  
  # Shrunk results
  res_shrunk_df <- as.data.frame(plot_list_shrunk[[i]])
  if ("gene" %in% colnames(res_shrunk_df)) rownames(res_shrunk_df) <- res_shrunk_df$gene
  if ("Gene.Symbol" %in% colnames(res_shrunk_df)) rownames(res_shrunk_df) <- res_shrunk_df$Gene.Symbol
  
  # Align to the same common genes
  res_shrunk_df <- res_shrunk_df[common_genes, , drop = FALSE]
  
  # Use shrunk LFCs if present
  if ("log2FC_shrunk" %in% colnames(res_shrunk_df)) {
    res_shrunk_df$log2FoldChange <- res_shrunk_df$log2FC_shrunk
  }
  
  # Add mean_counts and p-value column to use (prefer IHW)
  res_shrunk_df$mean_counts <- mean_norm_counts
  res_shrunk_df$padj_use <- if ("padj_ihw" %in% colnames(res_shrunk_df) && any(!is.na(res_shrunk_df$padj_ihw))) res_shrunk_df$padj_ihw else res_shrunk_df$padj
  
  # Keep rows with finite LFC, mean_counts
  res_shrunk_df <- res_shrunk_df[is.finite(res_shrunk_df$log2FoldChange) & is.finite(res_shrunk_df$mean_counts), , drop = FALSE]
  
  # Color by significance
  res_shrunk_df$Color <- ifelse(!is.na(res_shrunk_df$padj_use) & res_shrunk_df$padj_use < 0.05, "blue", "grey")
  res_shrunk_df$Color <- factor(res_shrunk_df$Color, levels = c("grey", "blue"))
  
  # MA plots
  p1 <- plot_ma(res_df, "MA plot: raw LFCs")
  p2 <- plot_ma(res_shrunk_df, "MA plot: ashr shrunk LFCs", show_legend = TRUE)
  
  # Save side-by-side
  png(filename = paste0("Results/DGEA/MA_plots/MA_plot_", names(plot_list_raw)[i], ".png"),
      res = 300, width = 4000, height = 2000)
  grid.arrange(p1, p2, ncol = 2)
  dev.off()
}

rm(res_df, res_shrunk_df, mean_norm_counts, common_genes, p1, p2); gc()

# Volcano plots #####
contrast_nice_names = c(
  res_tumor_rec_vs_tumor_no_rec             = "Recurrent tumors vs. non-recurrent tumors (overall)",
  res_late_tumor_rec_vs_early_tumor_rec     = "Late-stage recurrent vs. Early-stage recurrent tumors",
  res_early_tumor_rec_vs_early_tumor_no_rec = "Early-stage recurrent vs. Early-stage non-recurrent tumors",
  res_late_tumor_rec_vs_late_tumor_no_rec   = "Late-stage recurrent vs. Late-stage non-recurrent tumors"
)
legend_cols <- c(
  "Down-regulated"            = "royalblue",
  "Up-regulated"              = "red4",
  "|log2FC| < 1"              = "pink",
  "Not significant / missing" = "grey"
)

# Volcano loop (raw LFCs)
volcanos_raw = list()
for (i in 1:length(results_tbls)) {
  # Construct the data frame
  genes = if ("gene" %in% colnames(results_tbls[[i]])) results_tbls[[i]]$gene else rownames(results_tbls[[i]])
  res = results_tbls[[i]] %>%
    as.data.frame() %>%
    mutate(Gene.Symbol = if ("Gene.Symbol" %in% colnames(.)) Gene.Symbol else genes)
  
  # prefer IHW-adjusted p where available
  if ("padj_ihw" %in% colnames(res)) {
    res$padj_bh <- res$padj
    res$padj    <- res$padj_ihw
  }
  res <- dplyr::arrange(res, padj)
  
  # create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
  keyvals.colour <- ifelse(
    res$log2FoldChange < -1 & res$padj < 0.05, 'royalblue',
    ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red4',
           ifelse(abs(res$log2FoldChange) < 1 & res$padj < 0.05, 'pink',
                  'grey')))
  
  keyvals.colour[is.na(keyvals.colour)] <- 'grey'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
  names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
  names(keyvals.colour)[keyvals.colour == 'pink'] <- '|log2FC| < 1'
  names(keyvals.colour)[keyvals.colour == "grey"] <- "Not significant / missing"
  
  res$aes = keyvals.colour
  
  volcanos_raw[[names(results_tbls)[i]]] = EnhancedVolcano(res,
                                                           lab = res[, "Gene.Symbol"],
                                                           caption = NULL,
                                                           x = 'log2FoldChange',
                                                           y = 'padj',
                                                           pCutoff = 0.05,
                                                           FCcutoff = 1,
                                                           cutoffLineType = "dashed",
                                                           cutoffLineWidth = 0.7,
                                                           cutoffLineCol = "black",
                                                           colCustom = keyvals.colour,
                                                           colAlpha = 0.7,
                                                           ylab = bquote(bold(-log[10]("IHW adj. p-value"))),
                                                           xlab = bquote(bold("Shrunk"~log[2]("Fold Change"))),
                                                           pointSize = 3,
                                                           axisLabSize = 12,
                                                           subtitle = NULL,
                                                           labSize = 3.5,
                                                           selectLab = res[1:20, "Gene.Symbol"],
                                                           legendIconSize = 6,
                                                           title = paste0(contrast_nice_names[i]),
                                                           labFace = "bold",
                                                           boxedLabels = TRUE,
                                                           drawConnectors = TRUE,
                                                           typeConnectors = "closed",
                                                           arrowheads = FALSE,
                                                           widthConnectors = 0.3,
                                                           max.overlaps = Inf) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(face = "bold", size = 14),
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
          legend.text = element_text(size = 10),
          legend.box.background = element_rect(colour = "black"))
  
  print(volcanos_raw[[names(results_tbls)[i]]])
  ggsave(filename = paste0("Volcano_raw_LFCs_", names(results_tbls)[i], ".png"),
         path = "Results/DGEA/Volcano_plots", 
         width = 162, height = 142, device = 'png', units = "mm",
         dpi = 700)
}
rm(genes, i, res, keyvals.colour, volcano)

# Volcano loop (ashr-shrunk LFCs)
volcanos_ashr = list()
for (i in 1:length(ashr_shrunk_tbls)) {
  # Construct the data frame
  genes = if ("gene" %in% colnames(ashr_shrunk_tbls[[i]])) ashr_shrunk_tbls[[i]]$gene else rownames(ashr_shrunk_tbls[[i]])
  res = ashr_shrunk_tbls[[i]] %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = if ("Gene.Symbol" %in% colnames(.)) Gene.Symbol else genes)
  
  # prefer IHW-adjusted p where available
  if ("padj_ihw" %in% colnames(res)) {
    res$padj_bh <- res$padj
    res$padj    <- res$padj_ihw
  }
  
  # use shrunk LFCs if available
  if ("log2FC_shrunk" %in% colnames(res)) {
    res$log2FoldChange <- res$log2FC_shrunk
  }
  res <- dplyr::arrange(res, padj)
  
  # create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
  keyvals.colour <- ifelse(
    res$log2FoldChange < -1 & res$padj < 0.05, 'royalblue',
    ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red4',
           ifelse(abs(res$log2FoldChange) < 1 & res$padj < 0.05, 'pink',
                  'grey')))
  
  keyvals.colour[is.na(keyvals.colour)] <- 'grey'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
  names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
  names(keyvals.colour)[keyvals.colour == 'pink'] <- '|log2FC| < 1'
  names(keyvals.colour)[keyvals.colour == "grey"] <- "Not significant / missing"
  
  res$aes = keyvals.colour
  
  volcanos_ashr[[names(ashr_shrunk_tbls)[i]]] = EnhancedVolcano(res,
                            lab = res[, "Gene.Symbol"],
                            caption = NULL,
                            x = 'log2FoldChange',
                            y = 'padj',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            cutoffLineType = "dashed",
                            cutoffLineWidth = 0.7,
                            cutoffLineCol = "black",
                            colCustom = keyvals.colour,
                            colAlpha = 0.7,
                            ylab = bquote(bold(-log[10]("IHW adj. p-value"))),
                            xlab = bquote(bold("Shrunk"~log[2]("Fold Change"))),
                            pointSize = 3,
                            axisLabSize = 12,
                            subtitle = NULL,
                            labSize = 3.5,
                            selectLab = res[1:20, "Gene.Symbol"],
                            legendIconSize = 6,
                            title = paste0(contrast_nice_names[i]),
                            labFace = "bold",
                            boxedLabels = TRUE,
                            drawConnectors = TRUE,
                            typeConnectors = "closed",
                            arrowheads = FALSE,
                            widthConnectors = 0.3,
                            max.overlaps = Inf) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(face = "bold", size = 14),
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
          legend.text = element_text(size = 10),
          legend.box.background = element_rect(colour = "black"))
  
  print(volcanos_ashr[[names(ashr_shrunk_tbls)[i]]])
  ggsave(filename = paste0("Volcano_ashr_shrunk_LFCs_", names(ashr_shrunk_tbls)[i], ".png"),
         path = "Results/DGEA/Volcano_plots", 
         width = 162, height = 142, device = 'png', units = "mm",
         dpi = 700)
}
rm(genes, i, res, keyvals.colour)

# Combine with patchwork
# Standardise the legends across all four plots
volcanos_raw_fixed <- lapply(volcanos_raw[1:4], function(p) {
  
  # Ensure unused categories have visible legend symbols
  p$layers[[1]]$show.legend <- TRUE
  
  p +
    scale_colour_manual(
      values = legend_cols,
      limits = names(legend_cols),
      breaks = names(legend_cols),
      drop = FALSE,
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        override.aes = list(size = 7)
      )
    )
})

# Combine into 2x2 grid
combined_volcano_raw <- wrap_plots(
  volcanos_raw_fixed,
  ncol = 2,
  byrow = TRUE,
  guides = "collect"
) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 14),
    plot.tag = element_text(face = "bold", size = 20)
  )

print(combined_volcano_raw)

ggsave(
  filename = "Volcano_raw_combined.png",
  plot = combined_volcano_raw,
  path = "Results/DGEA/Volcano_plots",
  width = 300,
  height = 350,
  units = "mm",
  dpi = 700
)

# Standardise the legends across all four plots
volcanos_ashr_fixed <- lapply(volcanos_ashr[1:4], function(p) {
  
  # Ensure unused categories have visible legend symbols
  p$layers[[1]]$show.legend <- TRUE
  
  p +
    scale_colour_manual(
      values = legend_cols,
      limits = names(legend_cols),
      breaks = names(legend_cols),
      drop = FALSE,
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        override.aes = list(size = 7)
      )
    )
})

# Combine into 2x2 grid
combined_volcano_ashr <- wrap_plots(
  volcanos_ashr_fixed,
  ncol = 2,
  byrow = TRUE,
  guides = "collect"
) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 14),
    plot.tag = element_text(face = "bold", size = 20)
  )

print(combined_volcano_ashr)

ggsave(
  filename = "Volcano_ashr_combined.png",
  plot = combined_volcano_ashr,
  path = "Results/DGEA/Volcano_plots",
  width = 300,
  height = 350,
  units = "mm",
  dpi = 700
)

# Post-DGEA heatmaps #####
# Number of features for plots
n_degs = 50

# VST transformation for tumor samples
vsd_tumor = varianceStabilizingTransformation(dds_Wald_filt[, dds_Wald_filt$Sample.type == "Tumor"],
                                                        blind = TRUE)

# Top n_degs variable genes fromt the VST object
topVarGenes = head(order(rowVars(assay(vsd_tumor)), decreasing = TRUE), n_degs)
matvar  = assay(vsd_tumor)[topVarGenes, ]
matvar  = matvar - rowMeans(matvar)
colnames(matvar) = colnames(assay(vsd_tumor))
anno = as.data.frame(colData(vsd_tumor)[, c("Sample.type", "Stage", "Recurrence")])

# Create a custom palette with 20*14 = 280 color breaks (one for each tile)
colors2 = colorRampPalette(c(viridisLite::mako(280)[140:280],
                             rev(viridisLite::rocket(280)[140:280])))(280)

png(paste0("Results/DGEA/Post_DGEA_heatmaps/VST_top", n_degs, "vargenes_samples_heatmap.png"),
     res = 700, width = 6760, 
     height = 4760)
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
pheno_filt = Pheno[Pheno$Patient.ID %in% colnames(dds_Wald_filt), ]
filtered_samples1 = pheno_filt$Patient.ID[pheno_filt$Sample.type == "Tumor"]
filtered_samples2 = pheno_filt$Patient.ID[pheno_filt$Sample.type == "Tumor" &
                                            pheno_filt$Recurrence == "Recurrence"]
filtered_samples3 = pheno_filt$Patient.ID[pheno_filt$Sample.type == "Tumor" &
                                            pheno_filt$Stage == "Early"]
filtered_samples4 = pheno_filt$Patient.ID[pheno_filt$Sample.type == "Tumor" &
                                            pheno_filt$Stage == "Late"]
filtered_samples = list(filtered_samples1, filtered_samples2,
                        filtered_samples3, filtered_samples4)
names(filtered_samples) = names(results_tbls)
rm(filtered_samples1, filtered_samples2, filtered_samples3, filtered_samples4); gc()
columns_of_interest = c("Recurrence", "Stage", "Recurrence", "Recurrence")

# Export standardized normalized counts
# Firstly, transform the normalized counts with a data-driven pseudocount
# https://www.biorxiv.org/content/10.1101/404962v1.full -> Process description

# We need a grouping for the OSCC samples. We are going to go for 
# Recurrence vs No recurrence
rec_samples = pheno_filt$Patient.ID[pheno_filt$Recurrence == "Recurrence"]
norec_samples = pheno_filt$Patient.ID[pheno_filt$Recurrence == "No recurrence"]

# Get size factors
# Mock sample conditions
sample_conditions <- data.frame(
  row.names = colnames(counts(dds_Wald_filt)),
  condition = factor(rep("condition", ncol(counts(dds_Wald_filt))))
)

# Create DESeqDataSet object
dds_null <- DESeqDataSetFromMatrix(countData = counts(dds_Wald_filt), 
                                   colData = sample_conditions, 
                                   design = ~ 1) 

# Estimate size factors
dds_null <- estimateSizeFactors(dds_null)

# Save size factors for later if needed
size_factors = dds_null$sizeFactor
names(size_factors) = dds_null@colData@rownames

# Get the normalized counts
normalized_counts <- counts(dds_null, normalized = TRUE)

# Assign proper row and column names to the normalized counts
rownames(normalized_counts) = rownames(counts(dds_Wald_filt))
colnames(normalized_counts) = colnames(counts(dds_Wald_filt))

# Estimate averages
avg_rec_size_factor = mean(size_factors[rec_samples], na.rm = TRUE)
avg_norec_size_factor = mean(size_factors[norec_samples], na.rm = TRUE)
sf = c(avg_rec_size_factor, avg_norec_size_factor)

# Lun et al. suggest pseudocount = max{1, r|1/smin - 1/smax|}, where r = 1 (suggestion)
pseudocount_RNA = max(c(1, abs(1/min(sf) - 1/max(sf)))) # 1

# We proceed with the log2(norm.counts + 1) transformation
log2.norm.counts.plus1 = log2(normalized_counts + pseudocount_RNA)

# Export
saveRDS(log2.norm.counts.plus1, "Resources_and_environments/log2.norm.counts.plus1.rds")

rm(rec_samples, norec_samples, sample_conditions, normalized_counts, dds_null, 
   size_factors, avg_norec_size_factor, avg_rec_size_factor, sf, pseudocount_RNA); gc()

# Heatmaps for raw LFCs
anno = anno[pheno_filt$Patient.ID, ]
for (i in 1:length(results_tbls)) {
  # Construct the data frame
  res = results_tbls[[i]] %>%
    as.data.frame() %>%
    mutate(Gene.Symbol = if ("Gene.Symbol" %in% colnames(.)) Gene.Symbol else
      if ("gene" %in% colnames(.)) gene else rownames(.)) %>%
    dplyr::arrange(if ("padj_ihw" %in% colnames(.)) padj_ihw else padj) %>%
    na.omit()
  
  # choose adj-p for filtering (prefer IHW)
  padj_use <- if ("padj_ihw" %in% colnames(res)) res$padj_ihw else res$padj
  
  topDEGs = res[which(padj_use < 0.05),]
  if (nrow(topDEGs) == 0) next
  topDEGs = topDEGs[seq_len(min(n_degs, nrow(topDEGs))), ]
  
  # Gene clustering with raw counts and Poisson distance
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
  genes_order = res %>% dplyr::filter((if ("padj_ihw" %in% colnames(res)) padj_ihw else padj) < 0.05) %>%
    arrange(log2FoldChange)
  # order after intersecting with the available rows
  ordered_genes = intersect(genes_order$Gene.Symbol, rownames(logmat))
  
  # compute gaps_row on the ordered subset; must be strictly inside (1..nrow-1)
  neg_flags = genes_order$log2FoldChange[match(ordered_genes, genes_order$Gene.Symbol)] < 0
  gaps_row = sum(neg_flags, na.rm = TRUE)
  if (is.na(gaps_row) || gaps_row <= 0 || gaps_row >= length(ordered_genes)) gaps_row = NULL
  
  sample_order = anno %>% dplyr::arrange(Recurrence)
  n_cols_full = nrow(sample_order)
  gaps_col = length(which(sample_order$Recurrence == "No recurrence"))
  if (gaps_col <= 0 || gaps_col >= n_cols_full) gaps_col = NULL
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[ordered_genes, sample_order, drop = FALSE]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[ordered_genes, sample_order, drop = FALSE]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  png(paste0("Results/DGEA/Post_DGEA_heatmaps/raw_lfcs/",
              names(results_tbls)[i], "_DEG_Poisson_clustered_samples_and_genes.png"),
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
           main = paste0(names(results_tbls)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Ordered genes and samples
  png(paste0("Results/DGEA/Post_DGEA_heatmaps/raw_lfcs/",
              names(results_tbls)[i], "_DEG_samples_ordered.png"),
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
           main = paste0(names(results_tbls)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Filtered versions
  filt_cols = intersect(filtered_samples[[i]], colnames(logmat_znorm))
  png(paste0("Results/DGEA/Post_DGEA_heatmaps/raw_lfcs/filtered_",
              names(results_tbls)[i], "_DEG_Poisson_clustered_samples_and_genes.png"),
       res = 700, width = 7260, 
       height = 5760)
  pheatmap(logmat_znorm[, filt_cols, drop = FALSE], annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = as.dist(DEGsamplePoisDistMatrix[filt_cols, filt_cols]),
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 4,
           fontsize_col = 6,
           main = paste0(names(results_tbls)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Ordered genes and samples (filtered)
  column_of_interest = columns_of_interest[i]
  indices = colnames(logmat_ordered_znorm) %in% filt_cols
  arrangement = pheno_filt[pheno_filt$Patient.ID %in% colnames(logmat_ordered_znorm)[indices], ]
  arrangement = arrangement[order(arrangement[[column_of_interest]]), ]
  if (nrow(arrangement) > 1) {
    gcuts = as.numeric(cumsum(table(arrangement[[column_of_interest]])))
    gcuts = gcuts[gcuts > 0 & gcuts < nrow(arrangement)]
    gaps_col_filt = if (length(gcuts)) gcuts else NULL
  } else {
    gaps_col_filt = NULL
  }
  new_indices = arrangement$Patient.ID
  
  png(paste0("Results/DGEA/Post_DGEA_heatmaps/raw_lfcs/filtered_",
              names(results_tbls)[i], "_DEG_samples_ordered.png"),
       res = 700, width = 7260, 
       height = 5760)
  pheatmap(logmat_ordered_znorm[, new_indices, drop = FALSE], 
           annotation_col = anno[new_indices, , drop = FALSE],
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           col = colors2,
           gaps_row = gaps_row,
           gaps_col = gaps_col_filt,
           cutree_cols = 2,
           fontsize = 8,
           fontsize_row = 4,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(results_tbls)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}
rm(i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col, gaps_row, neg_flags,
   ordered_genes, filt_cols, gaps_col_filt, arrangement, indices); gc()

# Heatmaps for ashr-shrunk LFCs
anno = anno[pheno_filt$Patient.ID, ]
for (i in 1:length(ashr_shrunk_tbls)) {
  # Construct the data frame
  res = ashr_shrunk_tbls[[i]] %>%
    as.data.frame() %>%
    mutate(Gene.Symbol = if ("Gene.Symbol" %in% colnames(.)) Gene.Symbol else
      if ("gene" %in% colnames(.)) gene else rownames(.)) %>%
    dplyr::arrange(if ("padj_ihw" %in% colnames(.)) padj_ihw else padj)
  
  # choose adj-p for filtering (prefer IHW)
  padj_use <- if ("padj_ihw" %in% colnames(res)) res$padj_ihw else res$padj
  topDEGs = res[which(padj_use < 0.05),]
  if (nrow(topDEGs) == 0) next
  topDEGs = topDEGs[seq_len(min(n_degs, nrow(topDEGs))), ]
  
  # Gene clustering with raw counts and Poisson distance
  matcount = assay(dds_Wald_filt)[topDEGs$Gene.Symbol, , drop = FALSE]
  DEGpoisd = PoissonDistance(t(counts(dds_Wald_filt[topDEGs$Gene.Symbol, , drop = FALSE])) + .Machine$double.eps)
  DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
  rownames(DEGsamplePoisDistMatrix) = dds_Wald_filt$Patient.ID
  colnames(DEGsamplePoisDistMatrix) = dds_Wald_filt$Patient.ID
  
  # Same for genes
  DEGpoisd_genes = PoissonDistance(counts(dds_Wald_filt[topDEGs$Gene.Symbol, , drop = FALSE]) + .Machine$double.eps)
  DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
  rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  colnames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  
  # Log-transformed matrix
  logmat = log2(matcount + 1)
  
  # The same matrix but order samples by recurrence status and genes by shrunk logFC
  # prefer log2FC_shrunk if present
  lfc_col <- if ("log2FC_shrunk" %in% colnames(res)) "log2FC_shrunk" else "log2FoldChange"
  genes_order_df = res %>% dplyr::filter(padj_use < 0.05) %>%
    dplyr::select(Gene.Symbol, !!lfc_col := all_of(lfc_col)) %>%
    arrange(.data[[lfc_col]])
  ordered_genes = intersect(genes_order_df$Gene.Symbol, rownames(logmat))
  
  # compute gaps_row on the ordered subset; must be strictly inside (1..nrow-1)
  neg_flags = genes_order_df[[lfc_col]][match(ordered_genes, genes_order_df$Gene.Symbol)] < 0
  gaps_row = sum(neg_flags, na.rm = TRUE)
  if (is.na(gaps_row) || gaps_row <= 0 || gaps_row >= length(ordered_genes)) gaps_row = NULL
  
  sample_order = anno %>% dplyr::arrange(Recurrence)
  n_cols_full = nrow(sample_order)
  gaps_col = length(which(sample_order$Recurrence == "No recurrence"))
  if (gaps_col <= 0 || gaps_col >= n_cols_full) gaps_col = NULL
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[ordered_genes, sample_order, drop = FALSE]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[ordered_genes, sample_order, drop = FALSE]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  png(paste0("Results/DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
              names(ashr_shrunk_tbls)[i], "_DEG_Poisson_clustered_samples_and_genes.png"),
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
           main = paste0(names(ashr_shrunk_tbls)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Ordered genes and samples
  png(paste0("Results/DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
              names(ashr_shrunk_tbls)[i], "_DEG_samples_ordered.png"),
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
           main = paste0(names(ashr_shrunk_tbls)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Filtered versions
  filt_cols = intersect(filtered_samples[[i]], colnames(logmat_znorm))
  png(paste0("Results/DGEA/Post_DGEA_heatmaps/ashr_shrunk/filtered_",
              names(ashr_shrunk_tbls)[i], "_DEG_Poisson_clustered_samples_and_genes.png"),
       res = 700, width = 7260, 
       height = 5760)
  pheatmap(logmat_znorm[, filt_cols, drop = FALSE], annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = as.dist(DEGsamplePoisDistMatrix[filt_cols, filt_cols]),
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 4,
           fontsize_col = 6,
           main = paste0(names(ashr_shrunk_tbls)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  # Ordered genes and samples
  column_of_interest = columns_of_interest[i]
  indices = colnames(logmat_ordered_znorm) %in% filt_cols
  arrangement = pheno_filt[pheno_filt$Patient.ID %in% colnames(logmat_ordered_znorm)[indices], ]
  arrangement = arrangement[order(arrangement[[column_of_interest]]), ]
  if (nrow(arrangement) > 1) {
    gcuts = as.numeric(cumsum(table(arrangement[[column_of_interest]])))
    gcuts = gcuts[gcuts > 0 & gcuts < nrow(arrangement)]
    gaps_col_filt = if (length(gcuts)) gcuts else NULL
  } else {
    gaps_col_filt = NULL
  }
  new_indices = arrangement$Patient.ID
  
  png(paste0("Results/DGEA/Post_DGEA_heatmaps/ashr_shrunk/filtered_",
              names(ashr_shrunk_tbls)[i], "_DEG_samples_ordered.png"),
       res = 700, width = 7260, 
       height = 5760)
  pheatmap(logmat_ordered_znorm[, new_indices, drop = FALSE], 
           annotation_col = anno[new_indices, , drop = FALSE],
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           col = colors2,
           gaps_row = gaps_row,
           gaps_col = gaps_col_filt,
           cutree_cols = 2,
           fontsize = 8,
           fontsize_row = 4,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(ashr_shrunk_tbls)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}
rm(i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order_df, ordered_genes, neg_flags,
   gaps_col, gaps_row, gaps_col_filt, filt_cols, arrangement, indices); gc()

# Export a session info text file for reproducibility
writeLines(capture.output(sessionInfo()), 
           paste0("Resources_and_environments/DGEA_sessionInfo.txt"))

# Export dds_Wald_filt object
saveRDS(dds_Wald_filt, "Resources_and_environments/final_dds.rds")

# Save environment for easy loading for DGEA later
save.image("Resources_and_environments/DGEA.RData")
