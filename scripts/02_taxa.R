# Identify differential taxa using `lefser` package

# 1. Define output directory and functions
# 2. Load data using utils.R
# 3. Execute main function run_analysis() which:
#   * Runs process_combination()
#   * Runs conduct_lefse(), which:
#   * Runs plot_lefse()
#   * Exports contingency tables as out_dir/contingency
#   * Exports qc results as out_dir/qc.csv
#   * Exports lefse results as out_dir/results.csv
#   * Exports plots as out_dir/plots[combination]_lefse.pdf


# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(lefser)
library(ggplot2)
library(dplyr)
library(stringr)
library(rlang)
source("R/utils.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/02_taxa"


# Define helper functions ------------------------------------------------------
# Pre-process a single cohort and extract metadata
process_combination <- function(
    cohort,
    cohort_name,
    condition,
    age_col = "age_category",
    condition_col = "study_condition",
    class_imbalance_mid = 2,
    class_imbalance_high = 5,
    class_imbalance_severe = 10,
    age_imbalance_thresh = 0.25,
    n_status_thresh = 10) {
  
  meta <- as.data.frame(colData(cohort))
  
  # Create contingency table data frame
  contingency <- as.data.frame(
    table(
      age_category = meta[[age_col]],
      condition = meta[[condition_col]],
      useNA = "ifany"))
  
  names(contingency)[names(contingency) == "Freq"] <- "count"
  
  contingency <- contingency |>
    dplyr::mutate(cohort = cohort_name) |>
    dplyr::select(cohort, age_category, condition, count)
  
  # Check class balance
  counts <- table(meta[[condition_col]])
  
  n_min <- min(counts)
  n_max <- max(counts)
  
  class_ratio <- n_max / n_min
  
  class_imbalance <- if (class_ratio > class_imbalance_severe) {
    paste0("severe (>", class_imbalance_severe, ")")
  } else if (class_ratio > class_imbalance_high) {
    paste0("high (>", class_imbalance_high, ")")
  } else if (class_ratio > class_imbalance_mid) {
    paste0("moderate (>", class_imbalance_mid, ")")
  } else {
    "low"
  }
  
  # Check age balance
  age_tab <- table(meta[[age_col]], meta[[condition_col]])
  age_prop <- prop.table(age_tab, margin = 2)
  
  group_diffs <- combn(ncol(age_prop), 2, function(cols) {
    sum(abs(age_prop[, cols[1]] - age_prop[, cols[2]]))
  })
  
  age_ratio <- max(group_diffs)
  
  age_imbalance <- if (age_ratio > age_imbalance_thresh) {
    paste0("high (>", age_imbalance_thresh, ")")
  } else {
    "low"
  }
  
  # Get counts
  counts <- as.list(table(meta[[condition_col]]))
  counts <- as.list(table(meta[[condition_col]]))
  
  n_control <- counts[["control"]] %||% 0
  n_condition <- counts[[condition]] %||% 0
  
  n_status <- character()
  
  if (n_control < n_status_thresh) {
    n_status <- c(n_status, paste0("n controls < ", n_status_thresh))
  }
  
  if (n_condition < n_status_thresh) {
    n_status <- c(n_status, paste0("n ", condition, " < ", n_status_thresh))
  }
  
  n_status <- if (length(n_status) == 0) {
    "OK"
  } else {
    paste(n_status, collapse = " | ")
  }
  
  return(list(
    cohort = cohort,
    metadata = meta,
    contingency = contingency,
    qc = list(
      cohort_name = cohort_name,
      condition = condition,
      combination = paste(cohort_name, condition, sep = "_"),
      class_ratio = class_ratio,
      class_imbalance = class_imbalance,
      age_ratio = age_ratio,
      age_imbalance = age_imbalance,
      n_control = n_control,
      n_condition = n_condition,
      n_status = n_status)))
  }


# Run LEfSe on a single combination
conduct_lefse <- function(
    processed,
    subclass_col = "age_category",
    condition_col = "study_condition",
    control_label = "control",
    seed = 1234) {
  
  # Skip if sample size is insufficient
  if (processed$qc$n_status != "OK") {
    message("  - skipped: insufficient sample size")
    
    processed$qc$lefse_status <- "insufficient sample size"
    processed$qc$n_features <- 0
    processed$lefse <- NULL
    return(processed)
  }
  
  cohort <- processed$cohort
  meta <- processed$metadata
  condition <- processed$qc$condition
  subclass <- NULL
  
  # Enforce ordering of control and disease
  meta$study_condition <- factor(
    meta$study_condition,
    levels = c(control_label, condition))
  
  cohort$study_condition <- meta$study_condition
  
  # Enable/disable subclass: age_category when:
  if (subclass_col %in% colnames(meta)) {
    
    age_tab <- table(meta$age_category, meta$study_condition)
    
    ## there is more than one age group
    if (nrow(age_tab) < 2) {
      message("  - subclass disabled: only one age group")
    
    } else {

      ## and each age group has both control x disease samples
      valid_strata <- apply(age_tab, 1, function(x) all(x > 0))
      
      if (all(valid_strata)) {
        subclass <- subclass_col
        message("  - subclass enabled: age stratification valid")
      } else {
        message("  - subclass disabled: some age groups lack class balance")
      }
    }
  }
  
  message(
    "  - running LEfSe", " (samples = ", ncol(cohort), ", features = ", nrow(cohort), ")"
  )
  
  # Filter terminal node 
  tn <- get_terminal_nodes(rownames(cohort))
  cohort <- cohort[tn, , drop = FALSE]
  
  # Aggregate species to genus
  abund <- SummarizedExperiment::assay(cohort)
  
  genus <- sub(".*\\|g__", "g__", rownames(abund))
  genus <- sub("\\|.*", "", genus)
  
  keep <- grepl("^g__", genus)
  
  abund <- abund[keep, , drop = FALSE]
  genus <- genus[keep]
  
  if (nrow(abund) == 0) {
    message("  - ERROR: no genus-level features found")
    
    processed$qc$lefse_status <- "no genus-level features"
    processed$qc$n_features <- 0
    processed$lefse <- NULL
    return(NULL)
  }
  
  abund_genus <- rowsum(
    abund,
    group = genus)
  
  # Rebuild object
  cohort <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = abund_genus),
    colData = meta)
  
  message(
    "  - features after genus aggregation = ",
    nrow(cohort))
  
  # Transform relative abundance
  cohort <- relativeAb(cohort)
  
  # Run LEfSe
  set.seed(seed)
  
  result <- tryCatch({
    suppressWarnings(
      suppressMessages(
        lefser(
          cohort,
          classCol = condition_col,
          subclassCol = subclass)))

  }, error = function(e) {
    NULL
  })
  
  # Convert result to data frame
  result_df <- tryCatch(
    as.data.frame(result),
    error = function(e) NULL)
  
  # Determine status
  if (is.null(result_df) || nrow(result_df) == 0) {
    
    processed$qc$lefse_status <- "no results"
    processed$qc$n_features <- 0
    processed$qc$plot_status <- "no LEfSe results"
    processed$lefse <- NULL
    
  } else {
    
    processed$qc$lefse_status <- "SUCCESS"
    processed$qc$n_features <- nrow(result_df)
    
    processed$lefse <- cbind(
      data.frame(
        combination = processed$qc$combination,
        cohort = processed$qc$cohort_name,
        condition = processed$qc$condition,
        stringsAsFactors = FALSE),
      result_df)
  }
  
  return(processed)
}


# Build a LEfSe bar plot for a single combination
plot_lefse <- function(
    analysed,
    plot_dir,
    min_features = 2,
    max_features = 10,
    control_label = "control") {
  
  df <- analysed$lefse
  qc <- analysed$qc
  condition <- qc$condition
  
  # Skip if too few features
  if (nrow(df) < min_features) {
    message("  - plot skipped: too few features")
    analysed$qc$plot_status <- "too few features"
    return(analysed)
  }
  
  # Convert condition label for plotting
  plot_condition <- gsub("_", " ", condition)
  plot_condition <- stringr::str_to_title(tolower(condition))
  
  # Extract genus name
  df <- df |>
    dplyr::mutate(
      taxon = dplyr::case_when(
        stringr::str_detect(features, "g__") ~
          sub(".*g__", "g__", features),
        TRUE ~ features),
      taxon = sub("\\|.*", "", taxon),
      taxon = sub("^g__", "", taxon),
      taxon = gsub("_", " ", taxon),
      scores = as.numeric(scores)) |>
    dplyr::filter(!is.na(scores)) |>
    dplyr::slice_max(
      order_by = abs(scores),
      n = max_features,
      with_ties = FALSE) |>
    dplyr::mutate(
      group_label = ifelse(
        scores > 0,
        condition, control_label))
  
  # Skip if no valid scores
  if (nrow(df) == 0) {
    message("  - plot skipped: no valid scores")
    analysed$qc$plot_status <- "no valid scores"
    return(analysed)
  }
  
  # Create plot
  p <- ggplot(
    df,
    aes(
      x = scores,
      y = stats::reorder(taxon, scores, FUN = mean),
      fill = group_label)) +
    
    geom_col(
      width = 0.7,
      colour = "white",
      linewidth = 0.2) +
    
    geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      colour = "grey40") +
    
    scale_fill_manual(
      values = c(setNames("#c5cb93ff", control_label), 
                 setNames("#dfa57eff", condition)),
      labels = c("Control", "Case"),
      name = "Enriched in") +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_discrete(labels = function(x) gsub("unclassified", "*", x)) +
    
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_text(size = 22, face = "italic"),
      axis.title.x = element_text(size = 22),
      axis.text.x = element_text(size = 20, margin = margin (t = 5)),
      
      legend.position = "bottom",
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()) +
    
    labs(x = "LDA score (log10)", y = NULL)
  
  # Save plot
  combination <- qc$combination
  safe_combination <- gsub("[/\\\\]", "", combination)
  
  plot_height <- max(4, nrow(df) * 0.45 + 1.5 )
  
  ggsave(
    filename = file.path(plot_dir, paste0(safe_combination, ".pdf")),
    plot = p,
    width = 7,
    height = plot_height,
    device = "pdf")
  
  analysed$qc$plot_status <- "SUCCESS"
  message("  - plot saved: ", safe_combination, ".pdf")
  
  invisible(p)
  
  return(analysed)
}


analyse_lefse <- function(lefse_res) {
  
  # All features
  all_features <- lefse_res %>%
    group_by(features) %>%
    summarise(
      n_condition = n_distinct(condition),
      conditions = paste(unique(condition), collapse = ", "),
      n_cohort = n_distinct(cohort),
      cohorts = paste(unique(cohort), collapse = ", "),
      scores = mean(scores)
    ) %>%
    arrange(desc(n_condition))
  
  # Find features unique to conditions
  unique_features <- lefse_res %>%
    group_by(features) %>%
    summarise(
      n_condition = n_distinct(condition),
      n_cohort = n_distinct(cohort),
      conditions = paste(unique(condition), collapse = ", "),
      cohorts = paste(unique(cohort), collapse = ", "),
      .groups = "drop"
    ) %>%
    filter(n_condition == 1)
  
  return(list(
    all_features = all_features,
    unique_features = unique_features))
  
}

# Define main function ---------------------------------------------------------
run_analysis <- function(cohorts, out_dir) {
  
  plot_dir <- file.path(out_dir, "plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list()
  
  for (cohort_name in names(cohorts)) {  
    
    # Process and extract study_conditions from cohort
    cohort <- cohorts[[cohort_name]]
    info <- process_conditions(cohort)
    
    message("\n", cohort_name, ": control vs ", info$n_condition, " condition/s")
    
    # Check cohort contains both case and control samples
    validation <- validate_conditions(info, cohort_name)
    if (!is.null(validation)) {next}
    
    for (condition in info$conditions) {
      
      message(" > Condition: ", condition)
      
      combination <- paste(cohort_name, condition, sep = "_")
      meta <- as.data.frame(colData(cohort))
      keep <- meta$study_condition %in% c("control", condition)
      cohort_subset <- cohort[, keep]
      
      # Process combination
      processed <- process_combination(
        cohort = cohort_subset,
        cohort_name = cohort_name,
        condition = condition)
      
      # Conduct LEfSe
      analysed <- conduct_lefse(
        processed = processed)
      
      # Plot figure
      if (!is.null(analysed$lefse)) {
        
        analysed <- plot_lefse(
          analysed = analysed,
          plot_dir = plot_dir)
        }

      results[[combination]] <- analysed
    }
  }
  
  contingency <- dplyr::bind_rows(lapply(results, `[[`, "contingency"))
  qc <- dplyr::bind_rows(lapply(results, `[[`, "qc"))
  lefse <- dplyr::bind_rows(lapply(results, `[[`, "lefse"))
  
  # Analyse LEfSe features
  analysed <- analyse_lefse(lefse)
  
  write.csv(contingency, file.path(out_dir, "contingency.csv"), row.names = FALSE)
  write.csv(qc, file.path(out_dir, "qc.csv"), row.names = FALSE)
  write.csv(lefse, file.path(out_dir, "results.csv"), row.names = FALSE)
  write.csv(analysed$all_features, file.path(out_dir, "all_features.csv"), row.names = FALSE)
  write.csv(analysed$unique_features, file.path(out_dir, "unique_features.csv"), row.names = FALSE)
 
  message("\nResults saved to: ", out_dir)
  
  return(list(
    contingency = contingency,
    qc = qc,
    lefse = results))
}


# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary.rds")


# Execute ----------------------------------------------------------------------
analysis_res <- run_analysis(primary_cohorts, out_dir)
