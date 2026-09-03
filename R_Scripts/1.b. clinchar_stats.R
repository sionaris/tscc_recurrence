# Libraries
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

# Import clinical characteristics
clinchar = read.csv("clinchar.csv", header = TRUE) %>%
  dplyr::filter(grepl("SAMPLEC", SAMPLE_ID))

# Variables to summarise
categorical_vars_pre <- c(
  "GENDER",
  "AGE_GROUP",
  "SITE",
  "SMOKING",
  "ALCOHOL",
  "PATHOLOGY_N_STAGE",
  "PCT_TUMOR_CELLS_FFPE",
  "GRADE",
  "PERINEURAL_INVASION",
  "VASCULAR_OR_LYMPHATIC_INVASION",
  "EXTRANODAL_EXTENSION_ENE",
  "POSITIVE_MARGINS",
  "DEPTH_OF_INVASION",
  "ADJUVANT_THERAPY",
  "TYPE_OF_RECURRENCE"
)

continuous_vars_pre <- c(
  "NUMBER_PACK_YEARS_SMOKED",
  "NUMBER_OF_LYMPH_NODES",
  "DAYS_TO_DEATH_OR_LAST_FUP",
  "DAYS_TO_RECURRENCE"
)

# Categorical / binary variables
for (v in categorical_vars_pre) {
  cat("\n\n###", v, "###\n")
  print(table(clinchar[[v]], useNA = "ifany"))
  cat("Effective n:", sum(!is.na(clinchar[[v]])), "\n")
}

# Continuous variables
for (v in continuous_vars_pre) {
  cat("\n\n###", v, "###\n")
  
  x <- clinchar[[v]]
  x_complete <- x[!is.na(x)]
  
  print(summary(x_complete))
  cat("Effective n:", length(x_complete), "\n")
  cat("Missing n:", sum(is.na(x)), "\n")
}


# Tests: setup and run
categorical_vars <- c(
  "GENDER",
  "AGE_GROUP",
  "SITE",
  "SMOKING",
  "ALCOHOL",
  "PATHOLOGY_N_STAGE",
  "PCT_TUMOR_CELLS_FFPE",
  "GRADE",
  "PERINEURAL_INVASION",
  "VASCULAR_OR_LYMPHATIC_INVASION",
  "EXTRANODAL_EXTENSION_ENE",
  "POSITIVE_MARGINS",
  "DEPTH_OF_INVASION",
  "ADJUVANT_THERAPY"
)

continuous_vars <- c(
  "NUMBER_OF_LYMPH_NODES"
)

# Clean basic values 

clinchar_test <- clinchar |>
  mutate(
    across(where(is.character), ~ na_if(trimws(.x), "")),
    GROUP = factor(
      GROUP,
      levels = c("T1T2R", "T1T2NR", "T3T4R", "T3T4NR")
    ),
    T_STAGE_GROUP = case_when(
      GROUP %in% c("T1T2R", "T1T2NR") ~ "T1T2",
      GROUP %in% c("T3T4R", "T3T4NR") ~ "T3T4",
      TRUE ~ NA_character_
    ),
    RECURRENCE_GROUP = case_when(
      GROUP %in% c("T1T2R", "T3T4R") ~ "R",
      GROUP %in% c("T1T2NR", "T3T4NR") ~ "NR",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(GROUP))


# Helper: categorical variable

run_fisher <- function(data, variable, grouping, family, comparison) {
  
  dat <- data |>
    select(all_of(c(variable, grouping))) |>
    filter(
      !is.na(.data[[variable]]),
      !is.na(.data[[grouping]])
    )
  
  tab <- table(dat[[variable]], dat[[grouping]])
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(
      tibble(
        family = family,
        comparison = comparison,
        variable = variable,
        variable_type = "Categorical",
        test = "Fisher's exact test",
        n = nrow(dat),
        statistic = NA_real_,
        effect_size = NA_real_,
        p_value = NA_real_,
        note = "Test not estimable: fewer than 2 observed levels"
      )
    )
  }
  
  ft <- fisher.test(tab)
  
  # Fisher provides an odds-ratio estimate only for 2 x 2 tables
  OR <- if (all(dim(tab) == c(2, 2))) {
    unname(ft$estimate)
  } else {
    NA_real_
  }
  
  tibble(
    family = family,
    comparison = comparison,
    variable = variable,
    variable_type = "Categorical",
    test = "Fisher's exact test",
    n = nrow(dat),
    statistic = NA_real_,
    effect_size = OR,
    p_value = ft$p.value,
    note = ifelse(
      all(dim(tab) == c(2, 2)),
      "Effect size = Fisher odds ratio",
      NA_character_
    )
  )
}


# Helper: continuous variable, >2 groups

run_kruskal <- function(data, variable, grouping, family, comparison) {
  
  dat <- data |>
    select(all_of(c(variable, grouping))) |>
    filter(
      !is.na(.data[[variable]]),
      !is.na(.data[[grouping]])
    )
  
  if (
    nrow(dat) < 2 ||
    dplyr::n_distinct(dat[[grouping]]) < 2
  ) {
    return(
      tibble(
        family = family,
        comparison = comparison,
        variable = variable,
        variable_type = "Continuous",
        test = "Kruskal-Wallis",
        n = nrow(dat),
        statistic = NA_real_,
        p_value = NA_real_,
        note = "Test not estimable"
      )
    )
  }
  
  kt <- kruskal.test(
    dat[[variable]] ~ as.factor(dat[[grouping]])
  )
  
  tibble(
    family = family,
    comparison = comparison,
    variable = variable,
    variable_type = "Continuous",
    test = "Kruskal-Wallis",
    n = nrow(dat),
    statistic = unname(kt$statistic),
    p_value = kt$p.value,
    note = NA_character_
  )
}


# Helper: continuous variable, 2 groups 

run_wilcox <- function(data, variable, grouping, family, comparison) {
  
  dat <- data |>
    select(all_of(c(variable, grouping))) |>
    filter(
      !is.na(.data[[variable]]),
      !is.na(.data[[grouping]])
    )
  
  if (
    nrow(dat) < 2 ||
    dplyr::n_distinct(dat[[grouping]]) < 2
  ) {
    return(
      tibble(
        family = family,
        comparison = comparison,
        variable = variable,
        variable_type = "Continuous",
        test = "Wilcoxon rank-sum",
        n = nrow(dat),
        statistic = NA_real_,
        p_value = NA_real_,
        note = "Test not estimable"
      )
    )
  }
  
  wt <- wilcox.test(
    dat[[variable]] ~ as.factor(dat[[grouping]]),
    exact = FALSE
  )
  
  tibble(
    family = family,
    comparison = comparison,
    variable = variable,
    variable_type = "Continuous",
    test = "Wilcoxon rank-sum",
    n = nrow(dat),
    statistic = unname(wt$statistic),
    p_value = wt$p.value,
    note = NA_character_
  )
}


# Family 1: global comparison across all four tumor groups

family1_cat <- map_dfr(
  categorical_vars,
  ~ run_fisher(
    data = clinchar_test,
    variable = .x,
    grouping = "GROUP",
    family = "1. Four-group global",
    comparison = "T1T2R vs T1T2NR vs T3T4R vs T3T4NR"
  )
)

family1_cont <- map_dfr(
  continuous_vars,
  ~ run_kruskal(
    data = clinchar_test,
    variable = .x,
    grouping = "GROUP",
    family = "1. Four-group global",
    comparison = "T1T2R vs T1T2NR vs T3T4R vs T3T4NR"
  )
)


# Family 2: early vs advanced T-stage

family2_cat <- map_dfr(
  categorical_vars,
  ~ run_fisher(
    data = clinchar_test,
    variable = .x,
    grouping = "T_STAGE_GROUP",
    family = "2. T-stage comparison",
    comparison = "T1T2 vs T3T4"
  )
)

family2_cont <- map_dfr(
  continuous_vars,
  ~ run_wilcox(
    data = clinchar_test,
    variable = .x,
    grouping = "T_STAGE_GROUP",
    family = "2. T-stage comparison",
    comparison = "T1T2 vs T3T4"
  )
)


# Family 3A: recurrence within T1-T2

early_df <- clinchar_test |>
  filter(GROUP %in% c("T1T2R", "T1T2NR")) |>
  mutate(
    WITHIN_STAGE_RECURRENCE = if_else(GROUP == "T1T2R", "R", "NR")
  )

family3a_cat <- map_dfr(
  categorical_vars,
  ~ run_fisher(
    data = early_df,
    variable = .x,
    grouping = "WITHIN_STAGE_RECURRENCE",
    family = "3. Within-stage recurrence",
    comparison = "T1T2R vs T1T2NR"
  )
)

family3a_cont <- map_dfr(
  continuous_vars,
  ~ run_wilcox(
    data = early_df,
    variable = .x,
    grouping = "WITHIN_STAGE_RECURRENCE",
    family = "3. Within-stage recurrence",
    comparison = "T1T2R vs T1T2NR"
  )
)


# Family 3B: recurrence within T3-T4

advanced_df <- clinchar_test |>
  filter(GROUP %in% c("T3T4R", "T3T4NR")) |>
  mutate(
    WITHIN_STAGE_RECURRENCE = if_else(GROUP == "T3T4R", "R", "NR")
  )

family3b_cat <- map_dfr(
  categorical_vars,
  ~ run_fisher(
    data = advanced_df,
    variable = .x,
    grouping = "WITHIN_STAGE_RECURRENCE",
    family = "3. Within-stage recurrence",
    comparison = "T3T4R vs T3T4NR"
  )
)

family3b_cont <- map_dfr(
  continuous_vars,
  ~ run_wilcox(
    data = advanced_df,
    variable = .x,
    grouping = "WITHIN_STAGE_RECURRENCE",
    family = "3. Within-stage recurrence",
    comparison = "T3T4R vs T3T4NR"
  )
)


# Family 4: characteristics of recurrence itself 

recurrent_df <- clinchar_test |>
  filter(GROUP %in% c("T1T2R", "T3T4R")) |>
  mutate(
    RECURRENT_STAGE = if_else(GROUP == "T1T2R", "T1T2R", "T3T4R"),
    TYPE_OF_RECURRENCE = na_if(trimws(TYPE_OF_RECURRENCE), "")
  )

# Recurrence type
family4_type <- run_fisher(
  data = recurrent_df,
  variable = "TYPE_OF_RECURRENCE",
  grouping = "RECURRENT_STAGE",
  family = "4. Recurrence characteristics",
  comparison = "T1T2R vs T3T4R"
)

# Time to recurrence
# Only recurrent patients are included, so non-recurrent zero placeholders
# cannot enter this test.
family4_time <- run_wilcox(
  data = recurrent_df,
  variable = "DAYS_TO_RECURRENCE",
  grouping = "RECURRENT_STAGE",
  family = "4. Recurrence characteristics",
  comparison = "T1T2R vs T3T4R"
)


# Combine all inferential results

clinical_test_results <- bind_rows(
  family1_cat,
  family1_cont,
  family2_cat,
  family2_cont,
  family3a_cat,
  family3a_cont,
  family3b_cat,
  family3b_cont,
  family4_type,
  family4_time
) |>
  arrange(family, comparison, p_value)

print(clinical_test_results, n = Inf)


# inspect only nominally significant findings

clinical_test_results |>
  filter(!is.na(p_value), p_value < 0.05) |>
  arrange(p_value)

# BH correction within each analysis family

clinical_test_results_bh <- clinical_test_results |>
  group_by(family, comparison) |>
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) |>
  ungroup() |>
  arrange(family, p_adj_BH, p_value)

print(clinical_test_results_bh, n = Inf)


# Nominal and BH-significant findings side-by-side

clinical_test_results_bh |>
  mutate(
    nominal_significant = !is.na(p_value) & p_value < 0.05,
    BH_significant = !is.na(p_adj_BH) & p_adj_BH < 0.05
  ) |>
  filter(nominal_significant | BH_significant) |>
  select(
    family,
    comparison,
    variable,
    test,
    n,
    p_value,
    p_adj_BH,
    nominal_significant,
    BH_significant
  ) |>
  arrange(family, p_adj_BH)

# Export a session info text file for reproducibility
writeLines(capture.output(sessionInfo()), 
           paste0("Resources_and_environments/Stats_sessionInfo.txt"))

# Save environment for easy loading for DGEA later
save.image("Resources_and_environments/Stats.RData")