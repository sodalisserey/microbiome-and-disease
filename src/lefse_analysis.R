# Conduct LEfSe analysis using lefser package
# 1. Read and extract primary cohorts via handler.R
# 2. Run main() which calls:
#   a. Analysis function: run_analysis()
#     * Clean and extract disease from study_conditions
#     * ANALYSIS BRANCH 1: healthy x 1 disease  -> 
#       - Run QC (contingency table, class/age balance, sample size)
#       - Perform LEfSe
#     * ANALYSIS BRANCH 2: healthy x >1 disease  -> 
#       - Create pair-wise subsets (healthy vs disease) and on each subset:
#       - Run QC (contingency table, class/age balance, sample size)
#       - Perform LEfSe
#   b. Export function: export_analysis()
#     * Bulk save analysis/QC log and analysis results

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(lefser)
source("src/handler.R")


# Define output directory ------------------------------------------------------
dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_dir <- "results/lefse_analysis"


# Load and extract data --------------------------------------------------------
primary_cohorts <- readRDS("data/primary_cohorts.rds") |>
  create_cohort_objects()


# Define helper functions ------------------------------------------------------
#' Run QC checks on a cohort and return a contingency table, class imbalance and 
#' age imbalance ratios and sample sizes
run_qc <- function(
    cohort,
    age_col = "age_category",
    condition_col = "study_condition",
    class_imbalance_mid = 2,
    class_imbalance_high = 5,
    class_imbalance_severe = 10,
    age_imbalance_thresh = 0.25
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
  
  class_imbalance <- if (class_ratio > class_imbalance_severe) {
    "severe"
  } else if (class_ratio > class_imbalance_high) {
    "high"
  } else if (class_ratio > class_imbalance_mid) {
    "moderate"
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

#' Format QC results into a summary data frame containing: comparison (health x 
#' disease), class ratio, class imbalance, age ratio, age imbalance and sample size
log_qc <- function(qc, comparison_name) {
  
  data.frame(
    comparison = comparison_name,
    class_ratio = qc$class_balance$class_ratio,
    class_imbalance = qc$class_balance$class_imbalance,
    age_ratio = qc$age_balance$age_ratio,
    age_imbalance = qc$age_balance$age_imbalance,
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

#' Log LEfSe analysis status and reason for a single comparison
log_lefser <- function(analysis, res, comparison_name) {
  
  if (is.null(res)) {
    
    analysis$analysis_log[[length(analysis$analysis_log) + 1]] <- data.frame(
      comparison = comparison_name,
      status = "FAILED",
      reason = "LEfSe returned NULL"
    )
    
  } else {
    
    res_df <- tryCatch({
      as.data.frame(res)
    }, error = function(e) {
      NULL
    })
    
    if (is.null(res_df)) {
      
      analysis$analysis_log[[length(analysis$analysis_log) + 1]] <- data.frame(
        comparison = comparison_name,
        status = "FAILED",
        reason = "Could not coerce LEfSe output to data frame"
      )
      
    } else {
      
      analysis$analysis_log[[length(analysis$analysis_log) + 1]] <- data.frame(
        comparison = comparison_name,
        status = "SUCCESS",
        reason = paste(
          "n_features =", nrow(res_df)
        )
      )
    }
  }
  
  return(analysis)
}

#' Convert named list of results into a combined data frame
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
write_csvs <- function(out_dir, data_list) {
  for (nm in names(data_list)) {
    write.csv(
      data_list[[nm]],
      file.path(out_dir, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }
  message("\nDone. CSVs written to: ", file.path(out_dir))
}


# Define main functions --------------------------------------------------------
#' Run full LEfSe analysis analysis and return list of results by executing the 
#' complete analysis workflow for a list of cohorts, including preprocessing, 
#' cohort validation, quality control, and LEfSe analysis
#'
#' Each cohort is evaluated for the presence of healthy (control) and disease
#' samples, and the appropriate analysis branch is applied:
#'   Skip cohorts with no healthy samples
#'   Skip cohorts with healthy samples only
#'   Run single comparison for cohorts with one disease
#'   Run multiple comparisons for cohorts with multiple diseases
#'
#' For each valid comparison, the function:
#'   Cleans the cohort using `clean_cohort()`
#'   Extracts disease group structure via `get_disease_groups()`
#'   Runs QC using `run_qc()`
#'   Runs LEfSe analysis using `run_lefser()`
#'   Logs results using `log_qc()` and `log_lefser()`
run_analysis <- function(primary_cohorts) {
  
  # Initialise storage objects
  analysis <- list(
    qc_summary = list(),
    contingency_tables = list(),
    analysis_log = list(),
    lefse_results = list()
  )
  
  for (cohort_name in names(primary_cohorts)) {  
    
    # Process and extract study_conditions from cohort
    cohort <- primary_cohorts[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Filter out cohorts with no control (healthy) samples
    if (!info$healthy_present || info$n_diseases == 0) {
      
      reason <- dplyr::case_when(
        !info$healthy_present ~ "No healthy samples",
        info$n_diseases == 0 ~ "Healthy samples only"
      )
      
      message(reason, " in ", cohort_name, ". Skipping.")
      
      analysis$analysis_log[[length(analysis$analysis_log) + 1]] <- data.frame(
        comparison = cohort_name,
        status = "SKIPPED",
        reason = reason
      )
      
      next
      
    # ANALYSIS BRANCH 1: healthy x 1 disease
    } else if (info$n_diseases == 1) {
      
      disease = info$diseases
      
      comparison_name <- paste(cohort_name, disease, sep = "_")
      message("\n", cohort_name, ": healthy vs ", info$diseases)
      
      # Run QC and log results
      qc <- run_qc(cohort)
      
      analysis$qc_summary[[length(analysis$qc_summary) + 1]] <- log_qc(qc, comparison_name)
      analysis$contingency_tables[[comparison_name]] <- qc$contingency_table
      
      # Perform Lefser, log analysis status/combination and save results
      res <- run_lefser(
        cohort,
        comparison = paste("control vs", disease),
        disease_label = disease
      )
      
      analysis <- log_lefser(analysis, res, comparison_name)
      
      analysis$lefse_results[[comparison_name]] <- res
      
      res_df <- as.data.frame(res)
      
      if (nrow(res_df) > 0) {
        res_df$cohort <- cohort_name
        res_df$disease <- disease
        res_df$comparison <- comparison_name
        
        analysis$lefse_results[[comparison_name]] <- res_df
      } else {
        message("No LEfSe hits for ", comparison_name)
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
        #cohort_subset <- clean_cohort(cohort_subset)

        comparison_name <- paste(cohort_name, disease, sep = "_")
        
        # Run QC and log results
        qc <- run_qc(cohort_subset)
        
        analysis$qc_summary[[length(analysis$qc_summary) + 1]] <- log_qc(qc, comparison_name)
        analysis$contingency_tables[[comparison_name]] <- qc$contingency_table
        
        # Perform LEfSe, log analysis status/combination and save results
        res <- run_lefser(
          cohort_subset,
          comparison = paste("control vs", disease),
          disease_label = disease
        )
        
        analysis <- log_lefser(analysis, res, comparison_name)
        
        analysis$lefse_results[[comparison_name]] <- res
        
        res_df <- as.data.frame(res)
        
        if (nrow(res_df) > 0) {
          res_df$cohort <- cohort_name
          res_df$disease <- disease
          res_df$comparison <- comparison_name
          
          analysis$lefse_results[[comparison_name]] <- res_df
        } else {
          message("No LEfSe hits for ", comparison_name)
        }
        }
      }
  }
  return(analysis)
}


#' Export LEfSe analysis analysis outputs as raw R objects and CSVs
export_analysis <- function(analysis, out_dir) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save full object
  saveRDS(analysis, file.path(out_dir, "lefse_analysis.rds"))
  
  # Bind logs
  analysis$analysis_log <- dplyr::bind_rows(analysis$analysis_log)
  analysis$qc_summary <- dplyr::bind_rows(analysis$qc_summary)
  
  # LEfSe results
  analysis$lefser_df <- bind_list_to_df(
    analysis$lefse_results,
    function(res) {
      df <- as.data.frame(res)
      if (nrow(df) == 0) return(NULL)
      
      df$disease <- attr(res, "case")
      df$control <- attr(res, "lclassf")
      df
    }
  )
  
  # Contingency tables
  analysis$contingency_df <- bind_list_to_df(
    analysis$contingency_tables,
    function(tbl) {
      df <- as.data.frame(as.table(tbl))
      names(df) <- c("age_category", "condition", "count")
      df
    }
  )
  
  # Write CSVs
  write_csvs(
    out_dir,
    list(
      lefse_results = analysis$lefser_df,
      contingency_tables = analysis$contingency_df,
      qc_summary = analysis$qc_summary,
      analysis_log = analysis$analysis_log
    )
  )
}


# Execute ----------------------------------------------------------------------
analysis <- run_analysis(primary_cohorts)
export_analysis(analysis, out_dir)
