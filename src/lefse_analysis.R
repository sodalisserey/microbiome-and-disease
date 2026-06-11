# Conduct LEfSe analysis using lefser package
# 1. Read, extract and process primary cohorts via handler.R
# 2. Define output directory and functions
# 3. Load data and execute main functions:
#   a. run_analysis()
#     * Clean, extract disease and validate conditions
#     * ANALYSIS BRANCH 1 (healthy x 1 disease) -> run_pipeline() executes:
#           - run_qc()
#           - run_lefser()
#           - log_qc()
#           - log_analysis()
#           - log_result()
#     * ANALYSIS BRANCH 2 (healthy x >1 disease)  -> 
#           - Create pair-wise subsets created
#           - Execute run_piipeline()
#   b. export_analysis()
#     * Bulk save analysis results into:
#           - qc_summary.csv
#           - contingency_tables.csv
#           - analysis_log.csv
#           - lefse_results.csv
#           - lefse_results.rds


# Load packages and dependencies -----------------------------------------------
library(rlang)
library(curatedMetagenomicData)
library(lefser)
source("src/handler.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/lefse_analysis"


# Define helper functions ------------------------------------------------------
#' Run QC on a cohort and return a contingency table, class imbalance and 
#' age imbalance ratios and counts
run_qc <- function(cohort,
                   cohort_name,
                   disease,
                   age_col = "age_category",
                   condition_col = "study_condition",
                   class_imbalance_mid = 2,
                   class_imbalance_high = 5,
                   class_imbalance_severe = 10,
                   age_imbalance_thresh = 0.25) {
  
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
      useNA = "ifany"
    )
  ))
  
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
    class_imbalance = class_imbalance
  )
  
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
    age_imbalance = age_imbalance
  )
  
  # Get counts
  counts <- as.list(table(meta[[condition_col]]))
  
  # Return results
  result <- list(
    contingency_table = contingency_df,
    class_balance = class_balance,
    age_balance = age_balance,
    n_controls = counts[["control"]] %||% 0,
    n_conditions = counts[[disease]] %||% 0
    )
  
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
run_lefser <- function(cohort,
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
    levels = c(control_label, disease_label)
  )
  
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
          subclassCol = subclass
        )
      )
    )
  }, error = function(e) {
    message("ERROR in ", comparison, ": ", e$message)
    NULL
  })
}

#' Format QC results into a summary data frame containing: comparison, 
#' class ratio, class imbalance, age ratio, age imbalance and sample size
log_qc <- function(comparison_name, 
                   cohort_name, 
                   disease, 
                   qc, 
                   analysis) {
  
  # Initialize containers if missing
  if (is.null(analysis$qc_summary)) {
    analysis$qc_summary <- list()
  }
  
  if (is.null(analysis$contingency_tables)) {
    analysis$contingency_tables <- list()
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
    
    stringsAsFactors = FALSE
  )
  
  analysis$qc_summary[[comparison_name]] <- row
  
  analysis$contingency_tables[[comparison_name]] <- 
    qc$contingency_table
  
  return(analysis)
}

#' Log LEfSe analysis status for a single comparison
log_analysis <- function(comparison_name, 
                       cohort_name, 
                       disease, 
                       status = NULL,
                       reason,
                       res, 
                       analysis) {
  
  # Initialize container if missing
  if (is.null(analysis$analysis_log)) {
    analysis$analysis_log <- list()
  }
  
  create_log <- function(status, reason) {
    data.frame(
      comparison = comparison_name,
      cohort = cohort_name,
      disease = disease,
      status = status,
      reason = reason,
      stringsAsFactors = FALSE
    )
  }
  
  # Case 1: explicit SKIPPED when condition_invalid == TRUE
  if (!is.null(status) && status == "SKIPPED") {
    analysis$analysis_log[[comparison_name]] <-
      create_log("SKIPPED", reason)
    return(analysis)
  }
  
  # Case 2: null result
  if (is.null(res)) {
    analysis$analysis_log[[comparison_name]] <-
      create_log("FAILED", "LEfSe returned NULL")
    return(analysis)
  }
  
  # Case 3: coercion
  res_df <- tryCatch(
    as.data.frame(res),
    error = function(e) NULL
  )
  
  if (is.null(res_df)) {
    analysis$analysis_log[[comparison_name]] <-
      create_log("FAILED", "Could not coerce LEfSe output to data frame")
    return(analysis)
  }
  
  # Case 4: success
  analysis$analysis_log[[comparison_name]] <-
    create_log("SUCCESS", paste("n_features =", nrow(res_df)))
  
  return(analysis)
  
}

#' Log LEfSe result from a single analysis
log_result <- function(comparison_name, 
                       cohort_name, 
                       disease, 
                       res, 
                       analysis) {
  
  # Initialize container if missing
  if (is.null(analysis$lefse_results)) {
    analysis$lefse_results <- list()
  }
  
  # Store metadata
  row <- data.frame(
    comparison = comparison_name,
    cohort = cohort_name,
    disease = disease,
    stringsAsFactors = FALSE
  )
  
  # No result
  if (is.null(res)) {
    return(analysis)
  }
  
  # Safe coercion
  res_df <- tryCatch(
    as.data.frame(res),
    error = function(e) NULL)

  # Failed coercion or empty
  if (is.null(res_df) || nrow(res_df) == 0) {
    return(analysis)
  }
  
  # Combine and store data frames
  combined <- cbind(row[rep(1, nrow(res_df)), , drop = FALSE], res_df)
  analysis$lefse_results[[comparison_name]] <- combined
  
  return(analysis)
}

#' Run analysis pipeline which includes QC checks, LEfSe analysis and logging
run_pipeline <- function(cohort,
                         cohort_name,
                         comparison_name,
                         disease,
                         analysis) {
  # Run and log QC
  qc <- run_qc(
    cohort,
    cohort_name = cohort_name, 
    disease = disease
  )
  
  analysis <- log_qc(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    qc = qc,
    analysis = analysis
  )
  
  # Run and log LEfSe
  res <- run_lefser(
    cohort,
    comparison = paste("control vs", disease),
    disease_label = disease
  )
  
  analysis <- log_analysis(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    res = res,
    analysis = analysis
  )
  
  # Store results in analysis
  analysis <- log_result(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    res = res,
    analysis = analysis
  )
  
  return(analysis)
}


# Define main functions --------------------------------------------------------
#' Run main LEfSe loop for two analysis branches
run_analysis <- function(primary_cohorts) {
  
  # Initialise storage object
  analysis <- list()
  
  for (cohort_name in names(primary_cohorts)) {  
    
    # Process and extract study_conditions from cohort
    cohort <- primary_cohorts[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check cohort contains both healthy and disease samples
    condition_invalid <- validate_conditions(info, cohort_name)
    
    if (!is.null(condition_invalid)) {
      analysis <- log_analysis(
        comparison_name,
        cohort_name,
        disease,
        status = "SKIPPED",
        reason = condition_invalid,
        analysis = analysis
      )
      next
      
    # ANALYSIS BRANCH 1: healthy x 1 disease
    } else if (info$n_diseases == 1) {
      
      disease = info$diseases
      
      comparison_name <- paste(cohort_name, disease, sep = "_")
      message("\n", cohort_name, ": healthy vs ", info$diseases)
      
      analysis <- run_pipeline(cohort = cohort,
                               cohort_name = cohort_name,
                               disease = disease,
                               comparison_name = comparison_name,
                               analysis = analysis)
      
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
        
        analysis <- run_pipeline(cohort = cohort_subset,
                                 cohort_name = cohort_name,
                                 disease = disease,
                                 comparison_name = comparison_name,
                                 analysis = analysis)
        }
      }
  }
  sapply(analysis$contingency_tables, function(x) {
    if (is.null(x)) return("NULL")
    if (length(x) == 0) return("EMPTY")
    class(x)
  })
  lapply(analysis$contingency_tables, dim)
  
  return(analysis)
}


#' Export LEfSe analysis outputs as raw R objects and CSVs
export_analysis <- function(analysis, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save full object
  saveRDS(analysis, file.path(out_dir, "lefse_results.rds"))
  
  # Bind data frames
  analysis$analysis_log <- dplyr::bind_rows(analysis$analysis_log)
  analysis$qc_summary <- dplyr::bind_rows(analysis$qc_summary)
  analysis$contingency_df <- dplyr::bind_rows(analysis$contingency_tables)
  analysis$lefser_df <- dplyr::bind_rows(analysis$lefse_results)
  
  # Write CSVs
  write_csvs(
    out_dir,
    list(
      lefse_results = analysis$lefser_df,
      contingency_tables = analysis$contingency_df,
      qc_summary = analysis$qc_summary,
      analysis_log = analysis$analysis_log))
}


# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary_cohorts.rds") |>
  create_cohort_objects()


# Execute ----------------------------------------------------------------------
analysis <- run_analysis(primary_cohorts)
export_analysis(analysis, out_dir)
