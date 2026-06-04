# Pipeline:
# 1. Read dataset registry from datasets.R and load cohorts via handler.R
# 2. For each cohort, clean and extract unique study conditions (disease):
#    * ANALYSIS BRANCH 1: healthy x 1 disease  -> 
#       - Run QC (contingency table, check class and age balance, sample size)
#         * log QC results
#       - Perform Lefser
#         * log analysis status and combination
#    * ANALYSIS BRANCH 2: healthy x >1 disease  -> 
#       - Create pair-wise subsets (healthy vs disease) and on each subset:
#         - Run QC (contingency table, check class and age balance, sample size)
#           * log QC results
#         - Perform Lefser
#           * log analysis status and combination
# 3. Bulk save pipeline results as RDS and CSV


# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(patchwork)
library(ggplot2)
source("src/datasets.R")
source("src/handler.R")


# Load data --------------------------------------------------------------------
# Optional: build primary_cohorts for scratch and save
# primary_cohorts <- load_cohorts(primary_cohort_names)
# 
# saveRDS(primary_cohorts, "data/primary_cohorts.rds")

# Load primary_cohorts after saving as rds
primary_cohorts <- readRDS("data/primary_cohorts.rds")

# Define helper functions ------------------------------------------------------

#' Clean cohort by removing invalid or missing study_condition values
#'
#' @param cohort A SummarizedExperiment / TreeSummarizedExperiment object
#' @return A filtered cohort object with invalid samples removed
clean_cohort <- function(cohort) {
  
  meta <- as.data.frame(colData(cohort))
  meta$study_condition <- as.character(meta$study_condition)
  
  # Remove samples where study_condition is NA, empty or whitespace only
  keep <- !is.na(meta$study_condition) &
    meta$study_condition != "" &
    meta$study_condition != " "
  
  cohort <- cohort[, keep]
  
  return(cohort)
}

#' Identify and extract healthy and disease groups from study_condition metadata
#' @param cohort A SummarizedExperiment / TreeSummarizedExperiment object
#' @param condition_col Column name containing sample condition labels
#' @param healthy_label Label used for healthy/control samples
#' @return A list with:
#' \describe{
#'   \item{healthy_present}{Logical, whether healthy samples exist}
#'   \item{diseases}{Character vector of disease labels}
#'   \item{n_diseases}{Number of disease conditions}
#' }
get_disease_groups <- function(cohort,
                               condition_col = "study_condition",
                               healthy_label = "control") {
  
  meta <- as.data.frame(colData(cohort))
  
  # Extract condition vector, remove NAs and return unique conditions
  conditions <- meta[[condition_col]]
  conditions <- conditions[!is.na(conditions)]
  conditions <- unique(conditions)
  
  # Disease labels = every study_condition except healthy
  diseases <- setdiff(conditions, healthy_label)
  
  list(
    healthy_present = healthy_label %in% conditions,
    diseases = diseases,
    n_diseases = length(diseases)
  )
}

#' Run QC checks on a cohort and compute:
#' - contingency table (age_category × study_condition)
#' - class imbalance 
#' - age distribution imbalance
#' - sample sizes per condition
#'
#' @param cohort A SummarizedExperiment / TreeSummarizedExperiment object
#' @param age_col Column name for age grouping variable
#' @param condition_col Column name for study condition
#' @param class_imbalance_mid Threshold for moderate class imbalance
#' @param class_imbalance_high Threshold for high class imbalance
#' @param age_imbalance_thresh Threshold for age distribution imbalance
#'
#' @return A list containing:
#' \describe{
#'   \item{contingency_table}{Age × condition contingency table}
#'   \item{class_balance}{List with class_ratio and class_imbalance}
#'   \item{age_balance}{List with age_ratio and age_status}
#'   \item{sample_size}{Data frame of sample counts per condition}
#' }
run_qc <- function(
    cohort,
    age_col = "age_category",
    condition_col = "study_condition",
    class_imbalance_mid = 2,
    class_imbalance_high = 600,
    age_imbalance_thresh = 0.3
) {
  
  meta <- as.data.frame(colData(cohort))
  
  # Ensure metadata exist
  if (!all(c(age_col, condition_col) %in% colnames(meta))) {
    stop("Required metadata columns not found")
  }
  
  # Create contingency table
  contingency_table <- table(
    meta[[age_col]],
    meta[[condition_col]],
    useNA = "ifany"
  )
  
  # Check class balance
  counts <- table(meta[[condition_col]])
  
  n_min <- min(counts)
  n_max <- max(counts)
  
  class_ratio <- n_max / n_min
  
  class_imbalance <- if (class_ratio > class_imbalance_high) {
    paste0("high (>", class_imbalance_high, ")")
  } else if (class_ratio > class_imbalance_mid) {
    paste0("med (>", class_imbalance_mid, ")")
  } else {
    "low"
  }
  
  class_balance <- list(
    class_ratio = class_ratio,
    class_imbalance = class_imbalance
  )
  
  # Check age balance
  age_tab <- table(meta[[age_col]], meta[[condition_col]])
  age_prop <- prop.table(age_tab, margin = 2)
  
  group_diffs <- combn(ncol(age_prop), 2, function(cols) {
    sum(abs(age_prop[, cols[1]] - age_prop[, cols[2]]))
  })
  
  age_ratio <- max(group_diffs)
  
  age_status <- if (age_ratio > age_imbalance_thresh) {
    paste0("high (>", age_imbalance_thresh, ")")
  } else {
    "low"
  }
  
  age_balance <- list(
    age_ratio = age_ratio,
    age_status = age_status
  )
  
  # Get sample size
  sample_size <- as.data.frame(table(meta[[condition_col]]))
  colnames(sample_size) <- c("condition", "n")
  
  # Return results
  result <- list(
    contingency_table = contingency_table,
    class_balance = class_balance,
    age_balance = age_balance,
    sample_size = sample_size
  )
  
  return(result)
}

#' Run LEfSe analysis on a cohort by performing:
#' - Terminal node filtering
#' - Relative abundance transformation
#' - Subclass stratification only if:
#'    * age_category exists
#'    * at least 2 age groups are present
#'    * all age groups contain both healthy and disease samples
#' - LEfSe differential abundance analysis
#'
#' @param cohort A SummarizedExperiment / TreeSummarizedExperiment object
#' @param comparison Optional name for logging/debugging
#' @param subclassCol Column used for subclass stratification (default: age_category)
#' @param seed Random seed for reproducibility
#'
#' @return LEfSe result object or NULL if analysis fails
run_lefser <- function(cohort,
                       comparison = NULL,
                       subclassCol = "age_category",
                       seed = 1234) {
  
  meta <- as.data.frame(colData(cohort))
  
  subclass <- NULL
  
  # Enable/disable subclass: age_category when:
  if ("age_category" %in% colnames(meta)) {
    
    age_tab <- table(meta$age_category, meta$study_condition)
    
    ## there is more than one age group
    if (nrow(age_tab) < 2) {
      message("Subclass disabled: only one age group")
      
    } else {
      
      ## and each age group has both healthy x disease samples
      valid_strata <- apply(age_tab, 1, function(x) all(x > 0))
      
      if (all(valid_strata)) {
        subclass <- "age_category"
        message("Subclass enabled: age stratification valid")
      } else {
        message("Subclass disabled: some age groups lack class balance")
      }
    }
  }
  
  message(
    "Running LEfSe",
    " | samples = ", ncol(cohort),
    " | features = ", nrow(cohort)
  )
  
  # Filter terminal node 
  tn <- get_terminal_nodes(rownames(cohort))
  cohort <- cohort[tn, , drop = FALSE]
  
  # Transform relative abundance
  cohort <- relativeAb(cohort)
  
  # Run LEfSe
  set.seed(seed)
  
  tryCatch({
    suppressWarnings(
      suppressMessages(
        lefser(
          cohort,
          classCol = "study_condition",
          subclassCol = subclass
        )
      )
    )
  }, error = function(e) {
    message("ERROR in ", comparison, ": ", e$message)
    NULL
  })
}

#' Format QC results into a summary data frame for logging purposes
#'
#' @param qc A list returned by \code{run_qc()}, containing class balance,
#' age balance, and sample size information.
#' @param comparison_name Character string identifying the cohort comparison.
#'
#' @return A one-row data frame containing:
#' \describe{
#'   \item{comparison}{Name of the cohort comparison}
#'   \item{class_ratio}{Ratio of largest to smallest class size}
#'   \item{class_imbalance}{Categorical description of class imbalance}
#'   \item{age_ratio}{Maximum age distribution difference between groups}
#'   \item{age_status}{Categorical description of age imbalance}
#'   \item{sample_size}{Formatted string of sample counts per condition}
#' }
log_qc <- function(qc, comparison_name) {
  
  data.frame(
    comparison = comparison_name,
    class_ratio = qc$class_balance$class_ratio,
    class_imbalance = qc$class_balance$class_imbalance,
    age_ratio = qc$age_balance$age_ratio,
    age_status = qc$age_balance$age_status,
    sample_size = paste(
      paste(
        qc$sample_size$condition,
        qc$sample_size$n,
        sep = "="
      ),
      collapse = "; "
    ),
    stringsAsFactors = FALSE
  )
}

#' Convert named list of results into a combined data frame
#'
#' @param x Named list of objects (e.g., LEfSe results or contingency tables)
#' @param transform_fn Function that converts each element into a data frame
#'
#' @return A combined data frame with a `comparison` column added
bind_list_to_df <- function(x, transform_fn) {
  
  do.call(
    rbind,
    lapply(names(x), function(nm) {
      
      obj <- x[[nm]]
      
      if (is.null(obj))
        return(NULL)
      
      df <- transform_fn(obj)
      
      if (is.null(df) || nrow(df) == 0)
        return(NULL)
      
      df$comparison <- nm
      
      df
    })
  )
}

#' Write multiple data frames to CSV files
#'
#' @param out_dir Output directory path
#' @param data_list Named list of data frames to write
#'
#' @return NULL (called for side effects)
write_csvs <- function(out_dir, data_list) {
  for (nm in names(data_list)) {
    write.csv(
      data_list[[nm]],
      file.path(out_dir, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }
}

# Storage objects --------------------------------------------------------------
pipeline <- list(
  qc_summary = list(),
  contingency_tables = list(),
  analysis_log = list(),
  lefser_results = list()
)

# Analysis ---------------------------------------------------------------------
# Iterate for every cohort in primary_cohorts
for (cohort_name in names(primary_cohorts)) {  
  
  # Clean cohort 
  cohort <- primary_cohorts[[cohort_name]]
  cohort <- clean_cohort(cohort)
  
  # Extract study_conditions
  info <- get_disease_groups(cohort)
  
  # Filter out cohorts with no control (healthy) samples
  if (!info$healthy_present) {
    message("No control (healthy samples) in ", cohort_name, ". Skipping.")
    
    ## log analysis status and combination
    pipeline$analysis_log[[length(pipeline$analysis_log) + 1]] <- data.frame(
      comparison = comparison_name,
      status = "SKIPPED",
      reason = "No healthy samples"
    )
    
    next
  }
  # Filter out cohorts with control (healthy) samples only
  if (info$n_diseases == 0) {
    message("Only healthy samples in ", cohort_name, ". Skipping.")
    
    ## log analysis status and combination
    pipeline$analysis_log[[length(pipeline$analysis_log) + 1]] <- data.frame(
      comparison = comparison_name,
      status = "SKIPPED",
      reason = "Healthy samples only"
    )
    
    
    # ANALYSIS BRANCH 1: healthy x 1 disease
  } else if (info$n_diseases == 1) {
    disease = info$diseases
    
    comparison_name <- paste(cohort_name, disease, sep = "_")
    message("\n", cohort_name, ": healthy vs ", info$diseases)
    
    # Run QC and log results
    qc <- run_qc(cohort)
    
    pipeline$qc_summary[[length(pipeline$qc_summary) + 1]] <- log_qc(qc, comparison_name)
    pipeline$contingency_tables[[comparison_name]] <- qc$contingency_table
    
    # Perform Lefser, log analysis status/combination and save results
    res <- run_lefser(cohort, comparison = paste("control vs", disease))
    
    if (is.null(res)) {
      pipeline$analysis_log[[length(pipeline$analysis_log) + 1]] <- data.frame(
        comparison = comparison_name,
        status = "FAILED",
        reason = "Lefse returned NULL"
      )
      
    } else {
      pipeline$analysis_log[[length(pipeline$analysis_log) + 1]] <- data.frame(
        comparison = comparison_name,
        status = "SUCCESS",
        reason = paste("n_features =", nrow(res))
      )
      
      pipeline$lefser_results[[comparison_name]] <- res
    }
    
    # ANALYSIS BRANCH 2: healthy x >1 disease
  } else {
    
    message("\n", cohort_name, ": healthy vs ", info$n_diseases, " diseases")
    
    for (disease in info$diseases) {
      
      message("\n", cohort_name, " subset: healthy vs ", disease)
      
      # Subset cohort for healthy x 1 disease
      meta <- as.data.frame(colData(cohort))
      keep <- meta$study_condition %in% c("control", disease)
      cohort_subset <- cohort[, keep]
      cohort_subset <- clean_cohort(cohort_subset)
      
      comparison_name <- paste(cohort_name, disease, sep = "_")
      
      # Run QC and log results
      qc <- run_qc(cohort_subset)
      
      pipeline$qc_summary[[length(pipeline$qc_summary) + 1]] <- log_qc(qc, comparison_name)
      pipeline$contingency_tables[[comparison_name]] <- qc$contingency_table
      
      # Perform Lefser, log analysis status/combination and save results
      res <- run_lefser(cohort_subset, comparison = paste("control vs", disease))
      
      if (is.null(res)) {
        pipeline$analysis_log[[length(pipeline$analysis_log) + 1]] <- data.frame(
          comparison = comparison_name,
          status = "FAILED",
          reason = "LEfSe returned NULL"
        )
        
      } else {
        pipeline$analysis_log[[length(pipeline$analysis_log) + 1]] <- data.frame(
          comparison = comparison_name,
          status = "SUCCESS",
          reason = paste("n_features =", nrow(res))
        )
        
        pipeline$lefser_results[[comparison_name]] <- res
      }
    }
  }
}


# Save results -----------------------------------------------------------------
# Create/set output directory
dir.create("results/lefser_analysis", recursive = TRUE, showWarnings = FALSE)
out_dir <- "results/lefser_analysis"

# Bulk save pipeline results 
saveRDS(pipeline, file.path(out_dir, "pipeline.rds"))

# Stack pipeline list of rows to single table 
pipeline$analysis_log <- dplyr::bind_rows(pipeline$analysis_log)
pipeline$qc_summary <- dplyr::bind_rows(pipeline$qc_summary)

# Convert pipeline lists to dataframes
pipeline$lefser_df <- bind_list_to_df(
  pipeline$lefser_results,
  function(res) as.data.frame(res)
)

pipeline$contingency_df <- bind_list_to_df(
  pipeline$contingency_tables,
  function(tbl) {
    df <- as.data.frame(as.table(tbl))
    names(df) <- c("age_category", "condition", "count")
    df
  }
)

# Save pipeline dfs as CSVs
write_csvs(
  out_dir,
  list(
    lefser_results_combined = pipeline$lefser_df,
    contingency_tables = pipeline$contingency_df,
    qc_summary = pipeline$qc_summary,
    analysis_log = pipeline$analysis_log
  )
)
