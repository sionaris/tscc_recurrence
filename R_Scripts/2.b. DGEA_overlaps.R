# Libraries #####
library(dplyr)
library(readr)
#library(stringr)
library(ggplot2)
#library(tidyr)
#library(ggpubr)
library(cowplot)
#library(pheatmap)
library(rcartocolor)
#library(viridisLite)
#library(RColorBrewer)
#library(tibble)
#library(ggrepel)
library(openxlsx)
#library(grid)
#library(gridExtra)
library(ComplexUpset)
library(patchwork)
library(ggnewscale)

# Import environment
load("Resources_and_environments/DGEA.RData")


# If TRUE, Fisher tests use the pairwise common tested-gene universe.
# However, all contrasts were tested for the same set of genes, so we set to FALSE
use_pairwise_common_tested_universe <- FALSE

# If TRUE, similarity denominators use only the pairwise-eligible significant genes
# (i.e. genes in the pairwise common tested universe). If FALSE, they use the full
# DEG list sizes as extracted from each comparison.Therefore, we also leave to FALSE
use_pairwise_eligible_sizes_for_similarity <- FALSE

# Figure output
out_dir <- "Results/DGEA/Overlaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

figure_prefix <- "DEG_overlap"
figure_width_in  <- 16
figure_height_in <- 22
figure_dpi <- 700

# Helpers for formatting etc
fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f%%"), 100 * x))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

fmt_or <- function(x, digits = 2) {
  ifelse(
    is.na(x), "",
    ifelse(is.infinite(x), "Inf", formatC(x, digits = digits, format = "fg"))
  )
}

fmt_p <- function(x) {
  ifelse(is.na(x), "", format.pval(x, digits = 2, eps = 1e-300))
}

# If duplicates exist, keep the row with the largest absolute effect size.
deduplicate_by_gene <- function(df, gene_col = "gene", effect_col = "log2FoldChange") {
  df %>%
    dplyr::filter(!is.na(.data[[gene_col]])) %>%
    dplyr::mutate(.abs_effect = abs(.data[[effect_col]])) %>%
    dplyr::arrange(dplyr::desc(.abs_effect)) %>%
    dplyr::distinct(.data[[gene_col]], .keep_all = TRUE) %>%
    dplyr::select(-.abs_effect)
}

extract_tested_genes <- function(df, gene_col = "gene", padj_col = "padj", effect_col = "log2FoldChange") {
  df %>%
    dplyr::filter(
      !is.na(.data[[gene_col]]),
      !is.na(.data[[padj_col]]),
      !is.na(.data[[effect_col]])
    ) %>%
    dplyr::pull(dplyr::all_of(gene_col)) %>%
    unique()
}

extract_sig_tbl <- function(df, alpha = 0.05, gene_col = "gene", padj_col = "padj", effect_col = "log2FoldChange") {
  df %>%
    dplyr::filter(
      !is.na(.data[[gene_col]]),
      !is.na(.data[[padj_col]]),
      !is.na(.data[[effect_col]]),
      .data[[padj_col]] < 0.05
    ) %>%
    dplyr::select(
      gene = dplyr::all_of(gene_col),
      log2FoldChange = dplyr::all_of(effect_col)
    ) %>%
    deduplicate_by_gene(gene_col = "gene", effect_col = "log2FoldChange")
}

to_matrix <- function(df, value_col, row_order, col_order) {
  mat <- matrix(
    NA,
    nrow = length(row_order),
    ncol = length(col_order),
    dimnames = list(row_order, col_order)
  )
  df = as.data.frame(df)
  for (r in row_order) {
    for (c in col_order) {
      mat[r, c] <- df[df$row_comp == r & df$col_comp == c, value_col, drop = TRUE]
    }
  }
  mat
}

# Prepare input
# comp_names <- c("Recurrent vs. Non-recurrent tumors",
#                 "Late-stage recurrent vs. Early-stage recurrent tumors",
#                 "Early-stage: recurrent vs. non-recurrent tumors",
#                 "Late-stage: recurrent vs. non-recurrent tumors")

comp_names <- c("R vs NR", "LR vs ER", "ER vs ENR", "LR vs LNR")

if (is.null(comp_names) || any(comp_names == "")) {
  stop("ashr_shrunk_tbls must be a *named* list.")
}

sig_tbls <- lapply(
  ashr_shrunk_tbls,
  extract_sig_tbl,
  alpha = 0.05,
  gene_col = "gene",
  padj_col = "padj_ihw",
  effect_col = "log2FoldChange"
)

tested_sets <- lapply(
  ashr_shrunk_tbls,
  extract_tested_genes,
  gene_col = "gene",
  padj_col = "padj_ihw",
  effect_col = "log2FoldChange"
)

sig_sizes_full <- vapply(sig_tbls, nrow, integer(1))
n_comp <- length(comp_names)
global_common_universe <- Reduce(intersect, tested_sets)

build_membership_df <- function(sig_tbls_named_list) {
  all_genes <- sort(unique(unlist(lapply(sig_tbls_named_list, function(x) x$gene))))
  membership_df <- tibble::tibble(gene = all_genes)
  for (nm in names(sig_tbls_named_list)) {
    membership_df[[nm]] <- membership_df$gene %in% sig_tbls_named_list[[nm]]$gene
  }
  membership_df
}

# Compute pairwise metrics
calc_pair_metrics <- function(i, j) {
  
  nm_i <- comp_names[i]
  nm_j <- comp_names[j]
  
  # Pairwise common tested gene universe
  pair_universe <- if (use_pairwise_common_tested_universe) {
    intersect(tested_sets[[i]], tested_sets[[j]])
  } else {
    global_common_universe
  }
  
  # Significant tables restricted to the pairwise common universe
  sig_i_eligible <- sig_tbls[[i]] %>%
    dplyr::filter(gene %in% pair_universe) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  
  sig_j_eligible <- sig_tbls[[j]] %>%
    dplyr::filter(gene %in% pair_universe) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  
  # Full sizes vs eligible sizes
  n_i_full <- sig_sizes_full[[i]]
  n_j_full <- sig_sizes_full[[j]]
  n_i_eligible <- nrow(sig_i_eligible)
  n_j_eligible <- nrow(sig_j_eligible)
  
  # For similarity denominators, choose full or eligible sizes
  n_i_sim <- if (use_pairwise_eligible_sizes_for_similarity) n_i_eligible else n_i_full
  n_j_sim <- if (use_pairwise_eligible_sizes_for_similarity) n_j_eligible else n_j_full
  
  # Join by gene before comparing logFC signs
  overlap_tbl <- dplyr::inner_join(
    sig_i_eligible %>% dplyr::rename(logFC_i = log2FoldChange),
    sig_j_eligible %>% dplyr::rename(logFC_j = log2FoldChange),
    by = "gene"
  )
  
  overlap_k <- nrow(overlap_tbl)
  
  concordant_k <- if (overlap_k > 0) {
    sum(overlap_tbl$logFC_i * overlap_tbl$logFC_j > 0, na.rm = TRUE)
  } else {
    0L
  }
  
  discordant_k <- overlap_k - concordant_k
  
  union_size <- n_i_sim + n_j_sim - overlap_k
  
  row_recovery <- if (n_i_sim > 0) overlap_k / n_i_sim else NA_real_
  col_recovery <- if (n_j_sim > 0) overlap_k / n_j_sim else NA_real_
  
  row_same_dir_recovery <- if (n_i_sim > 0) concordant_k / n_i_sim else NA_real_
  col_same_dir_recovery <- if (n_j_sim > 0) concordant_k / n_j_sim else NA_real_
  
  overlap_coef <- if (min(n_i_sim, n_j_sim) > 0) {
    overlap_k / min(n_i_sim, n_j_sim)
  } else {
    NA_real_
  }
  
  jaccard <- if (union_size > 0) {
    overlap_k / union_size
  } else {
    NA_real_
  }
  
  sign_concordance <- if (overlap_k > 0) {
    concordant_k / overlap_k
  } else {
    NA_real_
  }
  
  signed_jaccard <- if (union_size > 0) {
    concordant_k / union_size
  } else {
    NA_real_
  }

  # Theoretical maxima
  
  max_possible_overlap <- min(n_i_sim, n_j_sim)
  
  max_possible_jaccard <- if (max(n_i_sim, n_j_sim) > 0) {
    min(n_i_sim, n_j_sim) / max(n_i_sim, n_j_sim)
  } else {
    NA_real_
  }
  
  # Direction-aware maximum for signed Jaccard.
  # This accounts for up/down composition in each DEG list.
  n_i_up <- sum(sig_i_eligible$log2FoldChange > 0, na.rm = TRUE)
  n_i_down <- sum(sig_i_eligible$log2FoldChange < 0, na.rm = TRUE)
  n_j_up <- sum(sig_j_eligible$log2FoldChange > 0, na.rm = TRUE)
  n_j_down <- sum(sig_j_eligible$log2FoldChange < 0, na.rm = TRUE)
  
  max_possible_concordant <- min(n_i_up, n_j_up) + min(n_i_down, n_j_down)
  
  max_possible_signed_jaccard <- if (max(n_i_sim, n_j_sim) > 0) {
    max_possible_concordant / max(n_i_sim, n_j_sim)
  } else {
    NA_real_
  }
  
  # Fisher exact test on the eligible universe
  universe_n <- length(pair_universe)
  
  a <- overlap_k
  b <- n_i_eligible - overlap_k
  c_ <- n_j_eligible - overlap_k
  d <- universe_n - (a + b + c_)
  
  if (i == j || universe_n == 0 || any(c(a, b, c_, d) < 0)) {
    odds_ratio <- NA_real_
    p_value <- NA_real_
  } else {
    ft <- fisher.test(
      matrix(c(a, b, c_, d), nrow = 2, byrow = TRUE),
      alternative = "greater"
    )
    odds_ratio <- unname(ft$estimate)
    p_value <- ft$p.value
  }
  
  tibble::tibble(
    row_index = i,
    col_index = j,
    row_comp = nm_i,
    col_comp = nm_j,
    
    n_row_sig_full = n_i_full,
    n_col_sig_full = n_j_full,
    n_row_sig_eligible = n_i_eligible,
    n_col_sig_eligible = n_j_eligible,
    pairwise_tested_universe = universe_n,
    
    n_row_up = n_i_up,
    n_row_down = n_i_down,
    n_col_up = n_j_up,
    n_col_down = n_j_down,
    
    overlap = overlap_k,
    concordant = concordant_k,
    discordant = discordant_k,
    
    row_recovery = row_recovery,
    col_recovery = col_recovery,
    row_same_dir_recovery = row_same_dir_recovery,
    col_same_dir_recovery = col_same_dir_recovery,
    
    jaccard = jaccard,
    overlap_coef = overlap_coef,
    sign_concordance = sign_concordance,
    signed_jaccard = signed_jaccard,
    
    max_possible_overlap = max_possible_overlap,
    max_possible_jaccard = max_possible_jaccard,
    max_possible_concordant = max_possible_concordant,
    max_possible_signed_jaccard = max_possible_signed_jaccard,
    
    odds_ratio = odds_ratio,
    p_value = p_value
  )
}

pairwise_metrics <- purrr::map_dfr(
  seq_len(n_comp),
  function(i) {
    purrr::map_dfr(seq_len(n_comp), function(j) calc_pair_metrics(i, j))
  }
)

# Correct off-diagonal p-values
upper_unique <- pairwise_metrics %>%
  dplyr::filter(row_index < col_index) %>%
  dplyr::mutate(fdr = p.adjust(p_value, method = "BH"))

fdr_map <- dplyr::bind_rows(
  upper_unique %>% dplyr::select(row_comp, col_comp, fdr),
  upper_unique %>% dplyr::transmute(row_comp = col_comp, col_comp = row_comp, fdr = fdr)
)

pairwise_metrics <- pairwise_metrics %>%
  dplyr::left_join(fdr_map, by = c("row_comp", "col_comp")) %>%
  dplyr::mutate(
    fdr = ifelse(row_index == col_index, NA_real_, fdr),
    neglog10_fdr = ifelse(is.na(fdr), NA_real_, -log10(pmax(fdr, .Machine$double.xmin)))
  ) %>% distinct()

# Export
# Full table
readr::write_csv(
  pairwise_metrics %>% dplyr::arrange(row_index, col_index),
  file.path(out_dir, paste0(figure_prefix, "_pairwise_metrics_full.csv"))
)

# Unique off-diagonal pairs only
readr::write_csv(
  pairwise_metrics %>%
    dplyr::filter(row_index < col_index) %>%
    dplyr::arrange(row_index, col_index),
  file.path(out_dir, paste0(figure_prefix, "_pairwise_metrics_unique_pairs.csv"))
)

# Word-friendly matrix
formatted_matrix <- pairwise_metrics %>%
  dplyr::mutate(
    cell = dplyr::case_when(
      row_index == col_index ~ paste0("n=", n_row_sig_full),
      TRUE ~ paste0(
        overlap, " (", concordant, " same-sign)",
        "; row ", fmt_pct(row_recovery),
        "; col ", fmt_pct(col_recovery),
        "; J=", fmt_num(jaccard, 3),
        "; OC=", fmt_num(overlap_coef, 3),
        "; SC=", fmt_pct(sign_concordance),
        "; p.adj=", ifelse(is.na(fdr), "NA", fmt_p(fdr))
      )
    )
  ) %>%
  dplyr::select(row_comp, col_comp, cell) %>%
  tidyr::pivot_wider(names_from = col_comp, values_from = cell)

readr::write_csv(
  formatted_matrix,
  file.path(out_dir, paste0(figure_prefix, "_formatted_matrix_for_Word.csv"))
)

# Save individual metric matrices
metric_names <- c(
  "overlap", "concordant", "discordant",
  "row_recovery", "col_recovery",
  "row_same_dir_recovery", "col_same_dir_recovery",
  "jaccard", "overlap_coef", "sign_concordance", "signed_jaccard",
  "max_possible_overlap",
  "max_possible_jaccard",
  "max_possible_concordant",
  "max_possible_signed_jaccard",
  "odds_ratio", "p_value", "fdr", "neglog10_fdr"
)

for (m in metric_names) {
  mat <- to_matrix(pairwise_metrics, m, row_order = comp_names, col_order = comp_names)
  write.csv(
    mat,
    file = file.path(out_dir, paste0(figure_prefix, "_", m, "_matrix.csv")),
    row.names = TRUE
  )
}

# Plotting
# Off-diagonal thresholds for text colour
overlap_q75_offdiag <- quantile(
  pairwise_metrics$overlap[pairwise_metrics$row_index != pairwise_metrics$col_index],
  probs = 0.75,
  na.rm = TRUE
)

neglog10_fdr_q75_offdiag <- quantile(
  pairwise_metrics$neglog10_fdr[pairwise_metrics$row_index != pairwise_metrics$col_index],
  probs = 0.75,
  na.rm = TRUE
)

plot_df <- pairwise_metrics %>%
  dplyr::mutate(
    row_comp = factor(row_comp, levels = rev(comp_names)),
    col_comp = factor(col_comp, levels = comp_names),
    
    is_diag = row_index == col_index,
    
    # fill variable for the raw-overlap heatmap.
    # Diagonal is excluded from the colour scale.
    fill_overlap = ifelse(is_diag, NA_real_, overlap),
    
    label_overlap = as.character(overlap),
    
    label_row_recovery = ifelse(
      is.na(row_recovery), "",
      paste0(fmt_pct(row_recovery), "\n(", overlap, "/", n_row_sig_full, ")")
    ),
    
    label_jaccard = ifelse(is.na(jaccard), "", fmt_pct(jaccard)),
    label_overlap_coef = ifelse(is.na(overlap_coef), "", fmt_pct(overlap_coef)),
    
    label_sign_concordance = ifelse(
      overlap == 0 | is.na(sign_concordance), "",
      paste0(fmt_pct(sign_concordance), "\n(", concordant, "/", overlap, ")")
    ),
    
    label_signed_jaccard = ifelse(is.na(signed_jaccard), "", fmt_pct(signed_jaccard)),
    
    label_significance = ifelse(
      is.na(fdr),
      "",
      paste0("OR=", fmt_or(odds_ratio), "\np.adj=", fmt_p(fdr))
    ),
    
    # Text colours
    # For raw overlap:
    # diagonal always black
    # off-diagonal threshold uses only off-diagonal values
    text_col_overlap = dplyr::case_when(
      is_diag ~ "black",
      overlap < overlap_q75_offdiag ~ "white",
      TRUE ~ "black"
    ),
    
    text_col_row_recovery = ifelse(row_recovery < 0.75, "white", "black"),
    text_col_jaccard = ifelse(jaccard < 0.75, "white", "black"),
    text_col_overlap_coef = ifelse(overlap_coef < 0.75, "white", "black"),
    text_col_sign_concordance = ifelse(sign_concordance < 0.75, "white", "black"),
    text_col_signed_jaccard = ifelse(signed_jaccard < 0.75, "white", "black"),
    
    text_col_significance = dplyr::case_when(
      is_diag ~ "black",
      neglog10_fdr < neglog10_fdr_q75_offdiag ~ "white",
      TRUE ~ "black"
    ),
    
    label_max_overlap = as.character(max_possible_overlap),
    
    label_max_jaccard = ifelse(
      is.na(max_possible_jaccard),
      "",
      fmt_pct(max_possible_jaccard)
    ),
    
    label_max_signed_jaccard = ifelse(
      is.na(max_possible_signed_jaccard),
      "",
      fmt_pct(max_possible_signed_jaccard)
    )
  )

# Heatmap using ggplot2
heatmap_plot <- function(df, fill_col, label_col, title,
                         text_col = NULL,
                         fill_limits = NULL,
                         na_fill = "grey95",
                         option = "E",
                         label_size = 3.0,
                         separate_diagonal = FALSE,
                         diagonal_fill = "white",
                         diagonal_outline = "#3B528B",
                         diagonal_outline_width = 1.0,
                         triangle = c("full", "lower", "upper"),
                         include_diagonal = TRUE,
                         axis_levels = NULL,
                         fill_title = NULL,
                         show_upper_max = FALSE,
                         upper_max_label_col = NULL,
                         upper_max_fill = "#009E73",
                         upper_max_legend_label = "theoretical maximum",
                         upper_max_text_colour = "black") {
  
  triangle <- match.arg(triangle)
  
  if (is.null(axis_levels)) {
    axis_levels <- unique(c(as.character(df$row_comp), as.character(df$col_comp)))
  }
  
  df <- df %>%
    dplyr::mutate(
      row_comp = factor(as.character(row_comp), levels = axis_levels),
      col_comp = factor(as.character(col_comp), levels = axis_levels)
    )
  
  # Ensure diagonal flag exists
  if (!"is_diag" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(is_diag = row_index == col_index)
  }
  
  # Base triangle restriction
  df_plot <- df %>%
    dplyr::filter(
      dplyr::case_when(
        triangle == "full"  ~ TRUE,
        triangle == "lower" ~ row_index > col_index | (include_diagonal & row_index == col_index),
        triangle == "upper" ~ row_index < col_index | (include_diagonal & row_index == col_index),
        TRUE ~ TRUE
      )
    )
  
  # If showing upper-triangle maxima, split the data:
  # lower + diagonal = observed scale; upper = constant maximum layer.
  if (show_upper_max) {
    df_observed <- df %>%
      dplyr::filter(row_index > col_index | (include_diagonal & row_index == col_index))
    
    df_upper_max <- df %>%
      dplyr::filter(row_index < col_index) %>%
      dplyr::mutate(.upper_max_legend = upper_max_legend_label)
  } else {
    df_observed <- df_plot
    df_upper_max <- NULL
  }
  
  if (separate_diagonal) {
    df_offdiag <- df_observed %>% dplyr::filter(!is_diag)
    df_diag <- df_observed %>% dplyr::filter(is_diag)
    
    p <- ggplot() +
      geom_tile(
        data = df_offdiag,
        aes(x = col_comp, y = row_comp, fill = .data[[fill_col]]),
        color = "white",
        linewidth = 0.6
      ) +
      geom_tile(
        data = df_diag,
        aes(x = col_comp, y = row_comp),
        fill = diagonal_fill,
        color = diagonal_outline,
        linewidth = diagonal_outline_width
      )
    
  } else {
    p <- ggplot(
      df_observed,
      aes(x = col_comp, y = row_comp, fill = .data[[fill_col]])
    ) +
      geom_tile(color = "white", linewidth = 0.6)
  }
  
  # Observed-value labels: lower triangle + diagonal
  if (is.null(text_col)) {
    p <- p +
      geom_text(
        data = df_observed,
        aes(x = col_comp, y = row_comp, label = .data[[label_col]]),
        size = label_size,
        lineheight = 0.95
      )
  } else {
    p <- p +
      geom_text(
        data = df_observed,
        aes(
          x = col_comp,
          y = row_comp,
          label = .data[[label_col]],
          colour = .data[[text_col]]
        ),
        size = label_size,
        lineheight = 0.95
      ) +
      scale_colour_identity()
  }
  
  # Continuous colour scale for observed values
  p <- p +
    scale_fill_viridis_c(
      name = fill_title,
      option = option,
      limits = fill_limits,
      na.value = na_fill,
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        order = 1
      )
    )
  
  # Upper-triangle theoretical maximum layer
  if (show_upper_max) {
    
    if (is.null(upper_max_label_col)) {
      stop("When show_upper_max = TRUE, please provide upper_max_label_col.")
    }
    
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_tile(
        data = df_upper_max,
        aes(x = col_comp, y = row_comp, fill = .upper_max_legend),
        color = "white",
        linewidth = 0.6
      ) +
      geom_text(
        data = df_upper_max,
        aes(
          x = col_comp,
          y = row_comp,
          label = .data[[upper_max_label_col]]
        ),
        colour = upper_max_text_colour,
        size = label_size,
        lineheight = 0.95
      ) +
      scale_fill_manual(
        name = NULL,
        values = stats::setNames(upper_max_fill, upper_max_legend_label),
        guide = guide_legend(
          order = 2,
          override.aes = list(
            colour = "white",
            linewidth = 0.6
          )
        )
      )
  }
  
  p +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = title) +
    scale_x_discrete(limits = axis_levels, drop = FALSE) +
    scale_y_discrete(limits = rev(axis_levels), drop = FALSE) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title = element_text(
        face = "bold",
        size = 12,
        margin = margin(b = 8)
      ),
      legend.title.align = 0.5,
      legend.text = element_text(size = 10)
    )
}

# Heatmaps
pA <- heatmap_plot(
  plot_df,
  fill_col = "fill_overlap",
  label_col = "label_overlap",
  text_col = "text_col_overlap",
  title = "Raw pairwise overlap",
  fill_limits = c(
    0,
    max(plot_df$overlap[plot_df$row_index > plot_df$col_index], na.rm = TRUE)
  ),
  fill_title = "Number\nof genes",
  option = "E",
  label_size = 8,
  separate_diagonal = TRUE,
  diagonal_fill = "white",
  diagonal_outline = "white",
  diagonal_outline_width = 1.2,
  axis_levels = comp_names,
  show_upper_max = TRUE,
  upper_max_label_col = "label_max_overlap",
  upper_max_fill = paste0(carto_pal(n = 7, name = "Burg")[5], "CC"),
  upper_max_legend_label = "theoretical\nmaximum",
  upper_max_text_colour = "black"
) +
  theme(
    axis.text = element_text(face = "bold", size = 19),
    plot.title = element_text(face = "bold", size = 23),
    legend.text = element_text(size = 15),
    legend.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                margin = ggplot2::margin(b = 10, unit = "pt"))
  )

pB <- heatmap_plot(
  plot_df,
  fill_col = "row_recovery",
  label_col = "label_row_recovery",
  text_col = "text_col_row_recovery",
  title = "Directional overlap (row-normalised)",
  fill_limits = c(0, 1),
  fill_title = "Proportion",
  option = "E",
  label_size = 5.5,
  separate_diagonal = TRUE,
  diagonal_fill = "white",
  diagonal_outline = "white",
  diagonal_outline_width = 1.2,
  triangle = "full",
  include_diagonal = TRUE,
  axis_levels = comp_names
) +
  theme(axis.text = element_text(face = "bold", size = 19),
        plot.title = element_text(face = "bold", size = 23),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                    margin = ggplot2::margin(b = 10, unit = "pt")))

pC <- heatmap_plot(
  plot_df,
  fill_col = "jaccard",
  label_col = "label_jaccard",
  text_col = "text_col_jaccard",
  title = "Jaccard similarity",
  fill_limits = c(0, 1),
  fill_title = "Score",
  option = "E",
  label_size = 7,
  separate_diagonal = TRUE,
  diagonal_fill = "white",
  diagonal_outline = "white",
  diagonal_outline_width = 1.2,
  axis_levels = comp_names,
  show_upper_max = TRUE,
  upper_max_label_col = "label_max_jaccard",
  upper_max_fill = paste0(carto_pal(n = 7, name = "Burg")[5], "CC"),
  upper_max_legend_label = "theoretical\nmaximum",
  upper_max_text_colour = "black"
) +
  theme(axis.text = element_text(face = "bold", size = 19),
        plot.title = element_text(face = "bold", size = 23),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                    margin = ggplot2::margin(b = 10, unit = "pt")))

pF <- heatmap_plot(
  plot_df,
  fill_col = "signed_jaccard",
  label_col = "label_signed_jaccard",
  text_col = "text_col_signed_jaccard",
  title = "Signed Jaccard similarity",
  fill_limits = c(0, 1), # (same direction overlap / union)
  fill_title = "Score",
  option = "E",
  label_size = 7,
  separate_diagonal = TRUE,
  diagonal_fill = "white",
  diagonal_outline = "white",
  diagonal_outline_width = 1.2,
  axis_levels = comp_names,
  show_upper_max = TRUE,
  upper_max_label_col = "label_max_signed_jaccard",
  upper_max_fill = paste0(carto_pal(n = 7, name = "Burg")[5], "CC"),
  upper_max_legend_label = "theoretical\nmaximum",
  upper_max_text_colour = "black"
) +
  theme(axis.text = element_text(face = "bold", size = 19),
        plot.title = element_text(face = "bold", size = 23),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                    margin = ggplot2::margin(b = 10, unit = "pt")))

pD <- heatmap_plot(
  plot_df,
  fill_col = "overlap_coef",
  label_col = "label_overlap_coef",
  text_col = "text_col_overlap_coef",
  title = "Overlap coefficient",
  fill_limits = c(0, 1),
  fill_title = "Score",
  option = "E",
  label_size = 8,
  separate_diagonal = TRUE,
  diagonal_fill = "white",
  diagonal_outline = "white",
  diagonal_outline_width = 1.2,
  triangle = "lower",
  include_diagonal = TRUE,
  axis_levels = comp_names
) +
  theme(axis.text = element_text(face = "bold", size = 19),
        plot.title = element_text(face = "bold", size = 23),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                    margin = ggplot2::margin(b = 10, unit = "pt")))

pE <- heatmap_plot(
  plot_df,
  fill_col = "sign_concordance",
  label_col = "label_sign_concordance",
  text_col = "text_col_sign_concordance",
  title = "Sign concordance among overlaps",
  fill_limits = c(0, 1),
  fill_title = "Score",
  option = "E",
  label_size = 6,
  separate_diagonal = TRUE,
  diagonal_fill = "white",
  diagonal_outline = "white",
  diagonal_outline_width = 1.2,
  triangle = "lower",
  include_diagonal = TRUE,
  axis_levels = comp_names
) +
  theme(axis.text = element_text(face = "bold", size = 19),
        plot.title = element_text(face = "bold", size = 23),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                    margin = ggplot2::margin(b = 10, unit = "pt")))

# Duplicate Fisher results in plot_df
compiler = function (x, y) {
  comps = c()
  for (i in 1:length(x)) {
      pair = sort(c(x[i],y[i]))
      pair = paste(pair, collapse = "_vs_")
      comps = c(comps, pair)
  }
  comps
}
plot_df$comp = compiler(plot_df$row_comp, plot_df$col_comp)

# If a comp entry has Fisher results, copy it to other rows with the same comp entry
stats_df = plot_df %>% filter(!is.na(fdr)) %>% 
  select(comp, fdr, neglog10_fdr, label_significance, text_col_significance) %>%
  na.omit()
for (comp in unique(stats_df$comp)) {
  plot_df$fdr[plot_df$comp == comp] = stats_df$fdr[stats_df$comp == comp]
  plot_df$neglog10_fdr[plot_df$comp == comp] = stats_df$neglog10_fdr[stats_df$comp == comp]
  plot_df$label_significance[plot_df$comp == comp] = stats_df$label_significance[stats_df$comp == comp]
  plot_df$text_col_significance[plot_df$comp == comp] = stats_df$text_col_significance[stats_df$comp == comp]
}

pG <- heatmap_plot(
  plot_df,
  fill_col = "neglog10_fdr",
  label_col = "label_significance",
  text_col = "text_col_significance",
  title = "Overlap enrichment: Fisher exact test",
  fill_limits = c(0, max(plot_df$neglog10_fdr, na.rm = TRUE)),
  fill_title = expression(bold(-log[10] * p["adj"])),
  option = "E",
  label_size = 4.5,
  separate_diagonal = TRUE,
  diagonal_fill = "white",
  diagonal_outline = "white",
  diagonal_outline_width = 1.2,
  triangle = "lower",
  include_diagonal = TRUE,
  axis_levels = comp_names
) +
  theme(axis.text = element_text(face = "bold", size = 19),
        plot.title = element_text(face = "bold", size = 23),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", size = 15, hjust = 0.5,
                                    margin = ggplot2::margin(b = 10, unit = "pt")))

# UpSet plot
membership_all <- build_membership_df(sig_tbls)
colnames(membership_all) = c("gene", comp_names)

pH <- ComplexUpset::upset(
  membership_all,
  intersect = comp_names,
  min_size = 1,
  width_ratio = 0.16,
  sort_sets = FALSE,
  sort_intersections_by = c("degree", "cardinality"),
  
  base_annotations = list(
    "Intersection size" = ComplexUpset::intersection_size(
      text = list(size = 7, vjust = -0.25)
    ) +
      ylab("Intersection size") +
      theme(
        axis.title.y = element_text(size = 22, face = "bold"),
        axis.text.y  = element_text(size = 21)
      )
  ),
  
  set_sizes = (
    ComplexUpset::upset_set_size() +
      ylab("Set size") +
      theme(
        axis.title.x = element_text(
          size = 22,
          face = "bold",
          margin = margin(t = -2, unit = "cm")   # move "Set size" x-title closer
        ),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16)
      )
  ),
  
  themes = ComplexUpset::upset_modify_themes(
    list(
      intersections_matrix = theme(
        axis.title.x = element_blank(),  # removes plain "group"
        axis.text.y = element_text(size = 20)
      )
    )
  )
) +
  labs(
    title = "UpSet plot of DEG intersections",
    subtitle = "All significant DEGs across comparisons"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 19),
    axis.text = element_text(size = 15)
  )

# Multi-figure
combined_plot <-
  ((pA | pB | pC) /
     (pD | pE | pF) /
     (pG | pH)) +
  patchwork::plot_annotation(
    tag_levels = "A") &
    theme(
      plot.tag = element_text(face = "bold", size = 22)
  )

# Save at 700 dpi
ggsave(
  filename = file.path(out_dir, paste0(figure_prefix, "_multifigure_700dpi.png")),
  plot = combined_plot,
  width = figure_width_in*1.5,
  height = figure_height_in,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(out_dir, paste0(figure_prefix, "_multifigure_700dpi.tiff")),
  plot = combined_plot,
  width = figure_width_in*1.5,
  height = figure_height_in,
  units = "in",
  dpi = figure_dpi,
  bg = "white",
  compression = "lzw",
  limitsize = FALSE
)

# Also save the plot object as RDS for later editing
saveRDS(
  combined_plot,
  file = file.path("Resources_and_environments", paste0(figure_prefix, "_multifigure_plot_object.rds"))
)

# Log summary
cat("\nDone.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Detected comparisons:\n")
print(sig_sizes_full)

cat("\nMain figure files written:\n")
cat(" - ", file.path(out_dir, paste0(figure_prefix, "_multifigure_700dpi.png")), "\n", sep = "")
cat(" - ", file.path(out_dir, paste0(figure_prefix, "_multifigure_700dpi.tiff")), "\n", sep = "")

cat("\nMain tables written:\n")
cat(" - ", file.path(out_dir, paste0(figure_prefix, "_pairwise_metrics_full.csv")), "\n", sep = "")
cat(" - ", file.path(out_dir, paste0(figure_prefix, "_pairwise_metrics_unique_pairs.csv")), "\n", sep = "")
cat(" - ", file.path(out_dir, paste0(figure_prefix, "_formatted_matrix_for_Word.csv")), "\n", sep = "")


save.image("C:/Users/as3582/OneDrive - University of Cambridge/Desktop/GitHub/oscc_seq_rec/Resources_and_environments/DGEA_post.RData")

# Export session info as txt
writeLines(capture.output(sessionInfo()), 
           "Resources_and_environments/DGEA_post_session_info.txt")
