# Conduct LEfSe analysis using lefser package
# 1. Define output directory and functions
# 2. Load data using utils.R
# 3. Execute main functions:
#   a. run_analysis()
#     * Clean, extract disease and validate conditions
#     * ANALYSIS BRANCH 1 (healthy x 1 disease) -> run_pipeline() executes:
#           - run_qc() 
#           - run_lefser() 
#           - log_qc()
#           - log_analysis()
#           - log_result()
#     * ANALYSIS BRANCH 2 (healthy x >1 disease)  -> 
#           - Create pair-wise subsets
#           - Execute run_pipeline()
#   b. export_analysis()
#     * Bulk save analysis results into:
#           - results/lefse_analysis/qc_summary.csv
#           - results/lefse_analysis/contingency_tables.csv
#           - results/lefse_analysis/analysis_log.csv
#           - results/lefse_analysis/lefse_results.csv
#           - results/lefse_analysis/lefse_results.rds


# Load packages and dependencies -----------------------------------------------
library(rlang)
library(curatedMetagenomicData)
library(lefser)
source("R/utils.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/lefse_analysis"


# Define helper functions ------------------------------------------------------
#' Run QC on a cohort and return a contingency table, class imbalance and 
#' age imbalance ratios and counts
run_qc <- function(
    cohort,
    cohort_name,
    disease,
    age_col = "age_category",
    condition_col = "study_condition",
    class_imbalance_mid = 2,
    class_imbalance_high = 5,
    class_imbalance_severe = 10,
    age_imbalance_thresh = 0.25,
    n_status_thresh = 5) {
  
  meta <- as.data.frame(colData(cohort))
  
  # Ensure metadata exist
  if (!all(c(age_col, condition_col) %in% colnames(meta))) {
    stop("Required metadata columns not found")
  }
  
  # Create contingency table data frame
  contingency_df <- as.data.frame(as.table(
    table(
      factor(meta[[age_col]]),
      factor(meta[[condition_col]]),
      useNA = "ifany")))
  
  names(contingency_df) <- c("age_category", "condition", "count")
  
  contingency_df$cohort <- cohort_name

  contingency_df <- contingency_df[, c(
    "cohort", "age_category", "condition", "count"
  )]
  
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
  
  class_balance <- list(
    class_ratio = class_ratio,
    class_imbalance = class_imbalance)
  
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
  
  age_balance <- list(
    age_ratio = age_ratio,
    age_imbalance = age_imbalance)
  
  # Get counts
  counts <- as.list(table(meta[[condition_col]]))
  counts <- as.list(table(meta[[condition_col]]))
  
  n_controls <- counts[["control"]] %||% 0
  n_conditions <- counts[[disease]] %||% 0
  
  n_status <- character()
  
  if (n_controls < n_status_thresh) {
    n_status <- c(n_status, paste0("n controls < ", n_status_thresh))
  }
  
  if (n_conditions < 5) {
    n_status <- c(n_status, paste0("n ", disease, " < ", n_status_thresh))
  }
  
  n_status <- if (length(n_status) == 0) {
    NA_character_
  } else {
    paste(n_status, collapse = " | ")
  }
  
  
  # Return results
  result <- list(
    contingency_table = contingency_df,
    class_balance = class_balance,
    age_balance = age_balance,
    n_controls = n_controls,
    n_conditions = n_conditions,
    n_status = n_status)
  
  return(result)
}

#' Analyse a cohort by performing:
#' - Terminal node filtering
#' - Relative abundance transformation
#' - Subclass stratification only if:
#'    * age_category exists
#'    * at least 2 age groups are present
#'    * all age groups contain both healthy and disease samples
#' - LEfSe differential abundance analysis
run_lefser <- function(
    cohort,
    comparison = NULL,
    subclassCol = "age_category",
    disease_label = NULL,
    control_label = "control",
    seed = 1234) {
  
  meta <- as.data.frame(colData(cohort))
  
  subclass <- NULL
  
  # Enforce ordering of control and disease
  meta$study_condition <- factor(
    meta$study_condition,
    levels = c(control_label, disease_label))
  
  cohort$study_condition <- meta$study_condition
  
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
          subclassCol = subclass)))
    
  }, error = function(e) {
    message("ERROR in ", comparison, ": ", e$message)
    NULL
  })
}

#' Log QC results from a single comparison including contingency table, class 
#' ratio, class imbalance, age ratio, age imbalance and condition counts, append
#' analysis_res$qc_summary and return analysis_res
log_qc <- function(
    comparison_name, 
    cohort_name, 
    disease, 
    qc, 
    analysis_res) {
  
  # initialise containers if missing
  if (is.null(analysis_res$qc_summary)) {
    analysis_res$qc_summary <- list()
  }
  
  if (is.null(analysis_res$contingency_tables)) {
    analysis_res$contingency_tables <- list()
  }

  row <- data.frame(
    comparison = comparison_name,
    cohort = cohort_name,
    disease = disease,
    
    class_ratio = qc$class_balance$class_ratio,
    class_imbalance = qc$class_balance$class_imbalance,
    
    age_ratio = qc$age_balance$age_ratio,
    age_imbalance = qc$age_balance$age_imbalance,
    
    n_controls = qc$n_controls,
    n_conditions = qc$n_conditions,
    n_status = qc$n_status,
    
    stringsAsFactors = FALSE)
  
  analysis_res$qc_summary[[comparison_name]] <- row
  
  analysis_res$contingency_tables[[comparison_name]] <- 
    qc$contingency_table
  
  return(analysis_res)
}

#' Log LEfSe analysis status for a single comparison, append 
#' analysis_res$analysis_log and return analysis_res
log_analysis <- function(
    comparison_name, 
    cohort_name, 
    disease, 
    status = NULL,
    reason = NULL,
    result, 
    analysis_res) {
  
  # initialise container if missing
  if (is.null(analysis_res$analysis_log)) {
    analysis_res$analysis_log <- list()
  }
  
  create_log <- function(status, reason) {
    data.frame(
      comparison = comparison_name,
      cohort = cohort_name,
      disease = disease,
      status = status,
      reason = reason,
      stringsAsFactors = FALSE)
  }
  
  # Case 1: explicit SKIPPED when condition_invalid == TRUE
  if (!is.null(status) && status == "SKIPPED") {
    analysis_res$analysis_log[[comparison_name]] <-
      create_log("SKIPPED", reason)
    return(analysis_res)
  }
  
  # Case 2: null result
  if (is.null(result)) {
    analysis_res$analysis_log[[comparison_name]] <-
      create_log("FAILED", "LEfSe returned NULL")
    return(analysis_res)
  }
  
  # Case 3: coercion
  result_df <- tryCatch(
    as.data.frame(result),
    error = function(e) NULL)
  
  if (is.null(result_df)) {
    analysis_res$analysis_log[[comparison_name]] <-
      create_log("FAILED", "Could not coerce LEfSe output to data frame")
    return(analysis_res)
  }
  
  # Case 4: success
  analysis_res$analysis_log[[comparison_name]] <-
    create_log("SUCCESS", paste("n_features =", nrow(result_df)))
  
  return(analysis_res)
}

#' Log LEfSe result from a single analysis, append analysis_res$lefse_results
#' and return analysis_res
log_result <- function(
    comparison_name, 
    cohort_name, 
    disease, 
    result, 
    analysis_res) {
  
  # initialise container if missing
  if (is.null(analysis_res$lefse_results)) {
    analysis_res$lefse_results <- list()
  }
  
  # Store metadata
  row <- data.frame(
    comparison = comparison_name,
    cohort = cohort_name,
    disease = disease,
    stringsAsFactors = FALSE)
  
  # No result
  if (is.null(result)) {
    return(analysis_res)
  }
  
  # Safe coercion
  result_df <- tryCatch(
    as.data.frame(result),
    error = function(e) NULL)

  # Failed coercion or empty
  if (is.null(result_df) || nrow(result_df) == 0) {
    return(analysis_res)
  }
  
  # Combine and store data frames
  combined <- cbind(row[rep(1, nrow(result_df)), , drop = FALSE], result_df)
  analysis_res$lefse_results[[comparison_name]] <- combined
  
  return(analysis_res)
}

#' Run analysis pipeline which includes QC checks, LEfSe analysis and logging
run_pipeline <- function(
    cohort,
    cohort_name,
    comparison_name,
    disease,
    analysis_res) {
  
  # Run and log QC
  qc <- run_qc(
    cohort,
    cohort_name = cohort_name, 
    disease = disease)
  
  analysis_res <- log_qc(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    qc = qc,
    analysis_res = analysis_res)
  
  # Run and log LEfSe
  result <- run_lefser(
    cohort,
    comparison = paste("control vs", disease),
    disease_label = disease)
  
  analysis_res <- log_analysis(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    result = result,
    analysis_res = analysis_res)
  
  # Store results in analysis
  analysis_res <- log_result(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    result = result,
    analysis_res = analysis_res)
  
  return(analysis_res)
}


# Define main functions --------------------------------------------------------
#' Run main LEfSe loop for two analysis branches
run_analysis <- function(primary_cohorts) {
  
  # Initialise storage object
  analysis_res <- list()
  
  for (cohort_name in names(primary_cohorts)) {  
    
    # Process and extract study_conditions from cohort
    cohort <- primary_cohorts[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check cohort contains both healthy and disease samples
    condition_invalid <- validate_conditions(info, cohort_name)
    
    if (!is.null(condition_invalid)) {
      analysis_res <- log_analysis(
        comparison_name,
        cohort_name,
        disease,
        status = "SKIPPED",
        reason = condition_invalid,
        analysis_res = analysis_res)
      next
      
    # ANALYSIS BRANCH 1: healthy x 1 disease
    } else if (info$n_diseases == 1) {
      
      disease = info$diseases
      
      comparison_name <- paste(cohort_name, disease, sep = "_")
      message("\n", cohort_name, ": healthy vs ", info$diseases)
      
      analysis_res <- run_pipeline(
        cohort = cohort,
        cohort_name = cohort_name,
        disease = disease,
        comparison_name = comparison_name,
        analysis_res = analysis_res)
      
    # ANALYSIS BRANCH 2: healthy x >1 disease
    } else {
      
      message("\n", cohort_name, ": healthy vs ", info$n_diseases, " diseases")
      
      for (disease in info$diseases) {
        
        message("\n", cohort_name, " subset: healthy vs ", disease)
        
        # Subset cohort for healthy x 1 disease
        meta <- as.data.frame(colData(cohort))
        keep <- meta$study_condition %in% c("control", disease)
        cohort_subset <- cohort[, keep]

        comparison_name <- paste(cohort_name, disease, sep = "_")
        
        analysis_res <- run_pipeline(cohort = cohort_subset,
                                 cohort_name = cohort_name,
                                 disease = disease,
                                 comparison_name = comparison_name,
                                 analysis_res = analysis_res)
        }
      }
  }
  sapply(analysis_res$contingency_tables, function(x) {
    if (is.null(x)) return("NULL")
    if (length(x) == 0) return("EMPTY")
    class(x)
  })
  lapply(analysis_res$contingency_tables, dim)
  
  return(analysis_res)
}


#' Export LEfSe analysis outputs as raw R objects and CSVs
export_analysis <- function(analysis_res, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save full object
  saveRDS(analysis_res, file.path(out_dir, "lefse_results.rds"))
  
  # Bind data frames
  analysis_res$analysis_log <- dplyr::bind_rows(analysis_res$analysis_log)
  analysis_res$qc_summary <- dplyr::bind_rows(analysis_res$qc_summary)
  analysis_res$contingency_df <- dplyr::bind_rows(analysis_res$contingency_tables)
  analysis_res$lefser_df <- dplyr::bind_rows(analysis_res$lefse_results)
  
  # Write CSVs
  write_csvs(
    out_dir,
    list(
      lefse_results = analysis_res$lefser_df,
      contingency_tables = analysis_res$contingency_df,
      qc_summary = analysis_res$qc_summary,
      analysis_log = analysis_res$analysis_log))
}


# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary_cohorts.rds") |>
  create_cohort_objects()


# Execute ----------------------------------------------------------------------
analysis_res <- run_analysis(primary_cohorts)
export_analysis(analysis_res, out_dir)
