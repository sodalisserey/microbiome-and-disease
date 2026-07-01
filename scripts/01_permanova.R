# PERMANOVA analysis using vegan package
# 1. Define output directory and functions
# 2. Load data using utils.R
# 3. Execute main function: run_analysis()
#   * Clean, extract disease and validate conditions
#   * Execute run_permanova() and log_permanova()
#   * Export permanova_res as out_dir/[cohort_name]_plot.pdf


# Load packages and dependencies -----------------------------------------------
library(vegan)
source("R/utils.R")

# Define output directory ------------------------------------------------------
out_dir <- "results/01_permanova"
library(readr)


# Define helper functions ------------------------------------------------------
#' Run PERMANOVA on a single cohort and save plot
run_permanova <- function(
    cohort,
    cohort_name,
    condition_col = "study_condition",
    permutations = 9999,
    dispersion_permutations = 999,
    out_dir) {
  
  plot_dir <- file.path(out_dir, "plots")
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  warnings_list <- NA
  error_msg <- NA
  
  result <- tryCatch({
    
    # Extract abundance
    abund <- as.matrix(assay(cohort))
    
    # Extract metadata
    meta <- as.data.frame(colData(cohort))
    
    # Ensure rownames are sample IDs
    if (is.null(rownames(meta))) {
      rownames(meta) <- colnames(cohort)
    }
    
    # Align samples correctly
    common_samples <- intersect(colnames(abund), rownames(meta))
    
    if (length(common_samples) < 3) {
      stop("Too few overlapping samples")
    }
    
    abund <- abund[, common_samples, drop = FALSE]
    meta  <- meta[common_samples, , drop = FALSE]
    
    # Distance
    X <- t(abund)
    bray <- vegdist(X, method = "bray")
    
    # PERMANOVA
    meta[[condition_col]] <- as.factor(meta[[condition_col]])
    
    ad <- adonis2(
      bray ~ meta[[condition_col]],
      permutations = permutations)
    
    # Dispersion
    bd <- betadisper(bray, meta[[condition_col]])
    bd_test <- permutest(bd, permutations = dispersion_permutations)
    
    list(
      adonis = ad,
      betadisper = bd,
      dispersion_test = bd_test,
      meta = meta
    )
    
  }, warning = function(w) {
    warnings_list <<- c(warnings_list, w$message)
    invokeRestart("muffleWarning")
    
  }, error = function(e) {
    error_msg <<- e$message
    return(NULL)
  })
  
  pdf(
    file.path(plot_dir, paste0(cohort_name, "_plot.pdf")),
    width = 8,
    height = 6)
  
  bd <- betadisper(bray, meta[[condition_col]])
  plot(bd, main = cohort_name)
  
  dev.off()
  
  list(
    result = result,
    warnings = warnings_list,
    error = error_msg
  )
}

#' Log PERMANOVA results and append to permanova_res()
log_permanova <- function(
    cohort_name,
    conditions,
    n_disease,
    result,
    permanova_res) {
  
  # If no results
  if (is.null(result$result)) {
    
    row <- data.frame(
      cohort_name = cohort_name,
      conditions = conditions,
      n_disease = n_disease,
      status = "FAILED",
      r2 = NA_real_,
      p_value = NA_real_,
      dispersion_p = NA_real_,
      n = NA_integer_,
      warnings = if (!is.null(result$warnings)) {
        paste(result$warnings, collapse = " | ")
      } else NA,
      error = if (!is.null(result$error)) {
        paste(result$error, collapse = " | ")
      } else NA,
      stringsAsFactors = FALSE)
    
    permanova_res[[cohort_name]] <- row
    return(permanova_res)
  }
  
  # If success
  res <- result$result
  ad <- res$adonis
  
  r2 <- if (!is.null(ad$R2)) ad$R2[1] else NA_real_
  p_value <- if (!is.null(ad$`Pr(>F)`)) ad$`Pr(>F)`[1] else NA_real_
  
  bd_test <- res$dispersion_test
  
  disp_p <- NA_real_
  if (!is.null(bd_test) &&
      !is.null(bd_test$tab) &&
      "Pr(>F)" %in% colnames(bd_test$tab)) {
    disp_p <- bd_test$tab$`Pr(>F)`[1]
  }
  
  n_samples <- if (!is.null(res$meta)) nrow(res$meta) else NA_integer_

  row <- data.frame(
    cohort_name = cohort_name,
    conditions = conditions,
    n_disease = n_disease,
    status = "SUCCESS",
    r2 = r2,
    p_value = p_value,
    dispersion_p = disp_p,
    n = n_samples,
    warnings = if (!is.null(result$warnings)) {
      paste(result$warnings, collapse = " | ")
    } else NA,
    error = result$error,
    stringsAsFactors = FALSE
  )
  
  permanova_res[[cohort_name]] <- row
  
  return(permanova_res)
}

#' Export PERMANOVA results by saving as CSV
export_analysis <- function(permanova_res, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Bind log
  permanova_df <- dplyr::bind_rows(permanova_res)
  
  # Save CSV
  write_csv(permanova_df, file.path(out_dir, "permanova_results.csv"))
}


# Define main function ---------------------------------------------------------
run_analysis <- function(cohorts, out_dir) {
  
  permanova_res <- list()
  
  for (cohort_name in names(cohorts)) {
    
    cohort <- cohorts[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check conditions
    condition_invalid <- validate_conditions(info, cohort_name)
    
    if (!is.null(condition_invalid)) {
      permanova_res <- log_permanova(
        cohort_name = cohort_name,
        conditions = NA_character_,
        n_disease = NA_character_,
        result = NULL,
        permanova_res = permanova_res)
      
      next
    }
    
    message("\nRunning PERMANOVA for cohort: ", cohort_name)
    
    # Run PERMANOVA on filtered cohort
    result <- run_permanova(
      cohort = info$cohort, 
      cohort_name = cohort_name,
      out_dir = out_dir)
    
    message("   Logging results")
    
    conditions = paste(unique(info$conditions), collapse = ", ")
    
    # Log results
    permanova_res <- log_permanova(
      cohort_name = cohort_name,
      conditions = conditions,
      n_disease = info$n_disease,
      result = result,
      permanova_res = permanova_res)
  }
  
  message("\nExporting analysis results")
  
  export_analysis(permanova_res, out_dir)
  message("Done. Results written to: ", out_dir)
  
  return(permanova_res)
}
  
    
# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary.rds") 


# Execute ----------------------------------------------------------------------
results <- run_analysis(primary_cohorts, out_dir)


