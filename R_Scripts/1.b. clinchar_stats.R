# Libraries
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(officer)
library(flextable)

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

############################################################
# TABLE 1: CLINICOPATHOLOGICAL CHARACTERISTICS
############################################################

# clinchar_test already contains the 32 tumor cases only,
# with blank character values converted to NA.
table1_dat <- clinchar_test

table1_groups <- c("T1T2R", "T1T2NR", "T3T4R", "T3T4NR")


# -------------------------------------------------------------------------
# Formatting helpers
# -------------------------------------------------------------------------

# Maximum 1 decimal place; remove trailing ".0"
fmt1 <- function(x) {
  if (length(x) == 0 || is.na(x)) return("")
  
  out <- formatC(
    round(x, 1),
    format = "f",
    digits = 1
  )
  
  sub("\\.0$", "", out)
}


# n (%)
fmt_n_pct <- function(n, denom) {
  if (denom == 0) return("")
  
  paste0(
    n,
    " (",
    fmt1(100 * n / denom),
    "%)"
  )
}


# Median (Q1-Q3)
fmt_med_iqr <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0)
    return("-")
  
  # If only one observation is available, an IQR is not informative
  if (length(x) == 1)
    return(fmt1(x))
  
  q <- quantile(
    x,
    probs = c(0.25, 0.50, 0.75),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
  
  paste0(
    fmt1(q[2]),
    " (",
    fmt1(q[1]),
    "\u2013",
    fmt1(q[3]),
    ")"
  )
}


# Group sample sizes
group_ns <- setNames(
  vapply(
    table1_groups,
    function(g) sum(table1_dat$GROUP == g),
    integer(1)
  ),
  table1_groups
)


# -------------------------------------------------------------------------
# Helper for categorical variables
# -------------------------------------------------------------------------

make_cat_block <- function(
    variable,
    characteristic,
    raw_levels,
    display_levels = raw_levels,
    include_missing = FALSE,
    block_id) {
  
  stopifnot(length(raw_levels) == length(display_levels))
  
  rows <- purrr::map2_dfr(
    raw_levels,
    display_levels,
    function(raw, lab) {
      
      vals <- purrr::map_chr(
        table1_groups,
        function(g) {
          
          d <- table1_dat[
            table1_dat$GROUP == g,
            ,
            drop = FALSE
          ]
          
          x <- d[[variable]]
          
          n <- sum(
            as.character(x) == raw,
            na.rm = TRUE
          )
          
          fmt_n_pct(n, nrow(d))
        }
      )
      
      total_n <- sum(
        as.character(table1_dat[[variable]]) == raw,
        na.rm = TRUE
      )
      
      tibble(
        Characteristic = characteristic,
        Category = lab,
        T1T2R = vals[1],
        T1T2NR = vals[2],
        T3T4R = vals[3],
        T3T4NR = vals[4],
        Total = as.character(total_n),
        Block = block_id
      )
    }
  )
  
  
  # Add n/a row only where missing values actually exist
  if (include_missing) {
    
    miss_vals <- purrr::map_chr(
      table1_groups,
      function(g) {
        
        d <- table1_dat[
          table1_dat$GROUP == g,
          ,
          drop = FALSE
        ]
        
        n_miss <- sum(is.na(d[[variable]]))
        
        if (n_miss == 0)
          ""
        else
          as.character(n_miss)
      }
    )
    
    total_miss <- sum(is.na(table1_dat[[variable]]))
    
    if (total_miss > 0) {
      
      rows <- bind_rows(
        rows,
        tibble(
          Characteristic = characteristic,
          Category = "n/a",
          T1T2R = miss_vals[1],
          T1T2NR = miss_vals[2],
          T3T4R = miss_vals[3],
          T3T4NR = miss_vals[4],
          Total = as.character(total_miss),
          Block = block_id
        )
      )
    }
  }
  
  rows
}


# -------------------------------------------------------------------------
# Helper for continuous variables
# -------------------------------------------------------------------------

make_cont_block <- function(
    variable,
    characteristic,
    block_id,
    recurrent_only = FALSE) {
  
  vals <- purrr::map_chr(
    table1_groups,
    function(g) {
      
      d <- table1_dat[
        table1_dat$GROUP == g,
        ,
        drop = FALSE
      ]
      
      # Time to recurrence only makes sense in recurrent patients
      if (recurrent_only &&
          !g %in% c("T1T2R", "T3T4R")) {
        return("-")
      }
      
      if (recurrent_only) {
        d <- d[
          d$RECURRENCE == 1,
          ,
          drop = FALSE
        ]
      }
      
      fmt_med_iqr(d[[variable]])
    }
  )
  
  tibble(
    Characteristic = characteristic,
    Category = "",
    T1T2R = vals[1],
    T1T2NR = vals[2],
    T3T4R = vals[3],
    T3T4NR = vals[4],
    Total = "",
    Block = block_id
  )
}


# -------------------------------------------------------------------------
# Construct Table 1
# -------------------------------------------------------------------------

table1 <- bind_rows(
  
  # Total patients
  tibble(
    Characteristic = "Total patients",
    Category = "",
    T1T2R = as.character(group_ns["T1T2R"]),
    T1T2NR = as.character(group_ns["T1T2NR"]),
    T3T4R = as.character(group_ns["T3T4R"]),
    T3T4NR = as.character(group_ns["T3T4NR"]),
    Total = as.character(nrow(table1_dat)),
    Block = 1
  ),
  
  
  # Gender
  make_cat_block(
    "GENDER",
    "Gender",
    c("MALE", "FEMALE"),
    c("Male", "Female"),
    block_id = 2
  ),
  
  
  # Age
  make_cat_block(
    "AGE_GROUP",
    "Age (years)",
    c("<40", "40-65", ">65"),
    c("<40", "40-65", ">65"),
    block_id = 3
  ),
  
  
  # Smoking
  make_cat_block(
    "SMOKING",
    "Smoking",
    c("1", "0"),
    c("Yes", "No"),
    include_missing = TRUE,
    block_id = 4
  ),
  
  
  # Pack-years
  make_cont_block(
    "NUMBER_PACK_YEARS_SMOKED",
    "Pack-years",
    block_id = 5
  ),
  
  
  # Smoking onset
  make_cont_block(
    "YEAR_OF_TOBACCO_SMOKING_ONSET",
    "Smoking onset (age)",
    block_id = 6
  ),
  
  
  # Alcohol
  make_cat_block(
    "ALCOHOL",
    "Alcohol consumption",
    c("1", "0"),
    c("Yes", "No"),
    include_missing = TRUE,
    block_id = 7
  ),
  
  
  # Otalgia
  make_cat_block(
    "OTALGIA",
    "Otalgia",
    c("1", "0"),
    c("Yes", "No"),
    include_missing = TRUE,
    block_id = 8
  ),
  
  
  # Clinical T stage
  make_cat_block(
    "CLINICAL_T_STAGE",
    "Clinical T Stage",
    c("T1", "T2", "T3", "T4"),
    c("cT1", "cT2", "cT3", "cT4"),
    block_id = 9
  ),
  
  
  # Pathological T stage
  make_cat_block(
    "PATHOLOGY_T_STAGE",
    "Pathological T Stage",
    c("T1", "T2", "T3", "T4"),
    c("pT1", "pT2", "pT3", "pT4"),
    block_id = 10
  ),
  
  
  # Pathological N stage
  make_cat_block(
    "PATHOLOGY_N_STAGE",
    "Pathological N Stage",
    c("N0", "N1orN2orN3"),
    c("pN0", "pN+"),
    block_id = 11
  ),
  
  
  # Number of lymph nodes
  make_cont_block(
    "NUMBER_OF_LYMPH_NODES",
    "Number of lymph nodes",
    block_id = 12
  ),
  
  
  # Pathological M stage
  make_cat_block(
    "PATHOLOGY_M_STAGE",
    "Pathological M Stage",
    c("M0", "M1", "Mx"),
    c("pM0", "pM1", "pMx"),
    block_id = 13
  ),
  
  
  # Overall pathological stage
  make_cat_block(
    "PATHOLOGIC_STAGE",
    "Pathologic Stage",
    c(
      "STAGE I",
      "STAGE II",
      "STAGE III",
      "STAGE IV"
    ),
    c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    ),
    block_id = 14
  ),
  
  
  # Tumor-cell fraction
  make_cat_block(
    "PCT_TUMOR_CELLS_FFPE",
    "% of Tumor Cells / FFPE (%)",
    c(
      "<40%",
      "41-60%",
      "61-80%",
      "81-100%"
    ),
    c(
      "<40",
      "41-60",
      "61-80",
      "81-100"
    ),
    block_id = 15
  ),
  
  
  # Grade
  make_cat_block(
    "GRADE",
    "Grade",
    c("G1", "G2", "G3"),
    c("G1", "G2", "G3"),
    block_id = 16
  ),
  
  
  # PNI
  make_cat_block(
    "PERINEURAL_INVASION",
    "Perineural Invasion (PNI)",
    c("1", "0"),
    c("Yes", "No"),
    block_id = 17
  ),
  
  
  # LVI
  make_cat_block(
    "VASCULAR_OR_LYMPHATIC_INVASION",
    "Lympho-Vascular Invasion (LVI)",
    c("1", "0"),
    c("Yes", "No"),
    block_id = 18
  ),
  
  
  # ENE
  make_cat_block(
    "EXTRANODAL_EXTENSION_ENE",
    "Extranodal Extension (ENE)",
    c("1", "0"),
    c("Yes", "No"),
    block_id = 19
  ),
  
  
  # Margins
  make_cat_block(
    "POSITIVE_MARGINS",
    "Positive Margins",
    c("1", "0"),
    c("Yes", "No"),
    block_id = 20
  ),
  
  
  # DOI
  make_cat_block(
    "DEPTH_OF_INVASION",
    "Depth of Invasion (DOI) (mm)",
    c("<5", "6-10mm", ">11"),
    c("<5", "6-10", ">11"),
    block_id = 21
  ),
  
  
  # Adjuvant therapy
  make_cat_block(
    "ADJUVANT_THERAPY",
    "Adjuvant Therapy",
    c(
      "RADIOTHERAPY",
      "CHEMOTHERAPY",
      "CHEMO-RADIOTHERAPY",
      "CLOSE FOLLOW UP"
    ),
    c(
      "Radiotherapy R/T",
      "Chemotherapy C/T",
      "Chemo-Radiotherapy",
      "Close Follow up"
    ),
    block_id = 22
  ),
  
  
  # Recurrence
  make_cat_block(
    "RECURRENCE",
    "Recurrence",
    c("1", "0"),
    c("Yes", "No"),
    block_id = 23
  ),
  
  
  # Type of recurrence:
  # only recurrent groups receive percentages
  purrr::map2_dfr(
    c("LOCAL", "NODAL", "DISTANT"),
    c("Local", "Nodal", "Distant"),
    function(raw, lab) {
      
      vals <- purrr::map_chr(
        table1_groups,
        function(g) {
          
          d <- table1_dat[
            table1_dat$GROUP == g,
            ,
            drop = FALSE
          ]
          
          if (!g %in% c("T1T2R", "T3T4R"))
            return("-")
          
          n <- sum(
            d$TYPE_OF_RECURRENCE == raw,
            na.rm = TRUE
          )
          
          fmt_n_pct(n, nrow(d))
        }
      )
      
      tibble(
        Characteristic = "Type of recurrence",
        Category = lab,
        T1T2R = vals[1],
        T1T2NR = vals[2],
        T3T4R = vals[3],
        T3T4NR = vals[4],
        Total = as.character(
          sum(
            table1_dat$TYPE_OF_RECURRENCE == raw,
            na.rm = TRUE
          )
        ),
        Block = 24
      )
    }
  ),
  
  
  # Follow-up
  make_cont_block(
    "DAYS_TO_DEATH_OR_LAST_FUP",
    "Days to death or last follow-up",
    block_id = 25
  ),
  
  
  # Time to recurrence
  make_cont_block(
    "DAYS_TO_RECURRENCE",
    "Days to recurrence",
    block_id = 26,
    recurrent_only = TRUE
  )
)

# -------------------------------------------------------------------------
# Export plain CSV as a reproducibility/checking copy
# -------------------------------------------------------------------------

dir.create(
  "Tables",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  table1 |>
    select(-Block),
  "Tables/Table1_clinicopathological_characteristics.csv",
  row.names = FALSE
)


# -------------------------------------------------------------------------
# Create publication-ready Word table
# -------------------------------------------------------------------------

table1_export <- table1 |>
  select(-Block)

ft <- flextable(
  table1_export
)


# Column headings
ft <- set_header_labels(
  ft,
  
  Characteristic = "Characteristic",
  Category = "Category",
  
  T1T2R = paste0(
    "T1-T2 with recurrence\n",
    "ETR\n",
    "(n=", group_ns["T1T2R"], ")"
  ),
  
  T1T2NR = paste0(
    "T1-T2 without recurrence\n",
    "ETNR\n",
    "(n=", group_ns["T1T2NR"], ")"
  ),
  
  T3T4R = paste0(
    "T3-T4 with recurrence\n",
    "LTR\n",
    "(n=", group_ns["T3T4R"], ")"
  ),
  
  T3T4NR = paste0(
    "T3-T4 without recurrence\n",
    "LTNR\n",
    "(n=", group_ns["T3T4NR"], ")"
  ),
  
  Total = paste0(
    "Total\n",
    "(n=", nrow(table1_dat), ")"
  )
)


# Merge characteristic names across their category rows
ft <- merge_v(
  ft,
  j = "Characteristic"
)


# -------------------------------------------------------------------------
# Minimal manuscript-style Word formatting
# -------------------------------------------------------------------------

# Start with NO borders or shading anywhere
ft <- border_remove(ft)

# Define rules
top_bottom_rule <- officer::fp_border(
  color = "black",
  width = 1.5
)

header_rule <- officer::fp_border(
  color = "black",
  width = 0.75
)


# Top rule: top of header
ft <- hline_top(
  ft,
  border = top_bottom_rule,
  part = "header"
)


# Bottom border under header only
ft <- hline_bottom(
  ft,
  border = header_rule,
  part = "header"
)


# Bottom rule: bottom of final body row
ft <- hline_bottom(
  ft,
  border = top_bottom_rule,
  part = "body"
)


# All header entries bold
ft <- bold(
  ft,
  part = "header"
)


# All entries in column 1 bold
ft <- bold(
  ft,
  j = 1,
  part = "body"
)


# Alignment

# Characteristic column: left + top aligned
ft <- align(
  ft,
  j = "Characteristic",
  align = "left",
  part = "all"
)

ft <- valign(
  ft,
  j = "Characteristic",
  valign = "top",
  part = "all"
)

# Category column: left aligned, vertically centred
ft <- align(
  ft,
  j = "Category",
  align = "left",
  part = "all"
)

ft <- valign(
  ft,
  j = "Category",
  valign = "top",
  part = "all"
)

# Group and total columns: horizontally and vertically centred
ft <- align(
  ft,
  j = c(
    "T1T2R",
    "T1T2NR",
    "T3T4R",
    "T3T4NR",
    "Total"
  ),
  align = "center",
  part = "all"
)

ft <- valign(
  ft,
  j = c(
    "T1T2R",
    "T1T2NR",
    "T3T4R",
    "T3T4NR",
    "Total"
  ),
  valign = "center",
  part = "all"
)


# Font size
ft <- fontsize(
  ft,
  size = 8,
  part = "all"
)

ft <- fontsize(
  ft,
  size = 8.5,
  part = "header"
)


# Column widths
ft <- width(
  ft,
  j = "Characteristic",
  width = 1.55
)

ft <- width(
  ft,
  j = "Category",
  width = 1.05
)

ft <- width(
  ft,
  j = c(
    "T1T2R",
    "T1T2NR",
    "T3T4R",
    "T3T4NR"
  ),
  width = 1.35
)

ft <- width(
  ft,
  j = "Total",
  width = 0.75
)

# -------------------------------------------------------------------------
# Table notes
# -------------------------------------------------------------------------

ft <- add_footer_lines(
  ft,
  values = paste0(
    "Categorical variables are presented as n (%); ",
    "continuous variables as median (IQR). ",
    "Percentages use the full group denominator; ",
    "missing categorical observations are shown as n/a."
  )
)

ft <- add_footer_lines(
  ft,
  values = paste0(
    "Pack-years and smoking-onset age are summarised from available ",
    "observations only. ETR: early-stage tumor with recurrence; ",
    "ETNR: early-stage tumor without recurrence; ",
    "LTR: late-stage tumor with recurrence; ",
    "LTNR: late-stage tumor without recurrence."
  )
)

ft <- fontsize(
  ft,
  size = 7,
  part = "footer"
)

ft <- italic(
  ft,
  part = "footer"
)


# -------------------------------------------------------------------------
# Export standalone landscape DOCX
# -------------------------------------------------------------------------

save_as_docx(
  ft,
  path = "Tables/Table1_clinicopathological_characteristics.docx",
  pr_section = officer::prop_section(
    page_size = officer::page_size(
      orient = "landscape"
    ),
    page_margins = officer::page_mar(
      top = 0.5,
      bottom = 0.5,
      left = 0.5,
      right = 0.5
    )
  )
)


cat(
  "\nCreated:\n",
  "  Tables/Table1_clinicopathological_characteristics.docx\n",
  "  Tables/Table1_clinicopathological_characteristics.csv\n"
)

# Export a session info text file for reproducibility
writeLines(capture.output(sessionInfo()), 
           paste0("Resources_and_environments/Stats_sessionInfo.txt"))

# Save environment for easy loading for DGEA later
save.image("Resources_and_environments/Stats.RData")