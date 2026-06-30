# Process ML models for external validation using CrcBiomeScreen package
# 1. Combine individual logs from 05_ml_train/logs/[cohort_name]_model_log.csv

# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
library(readr)
library(tools)
library(dplyr)
source("R/utils.R")


# Define input/output directories ----------------------------------------------
prep_dir <- "results/04_ml_train_prep"
train_dir <- "results/05_ml_train"
out_dir <- "results/06_ml_validate_prep"


# Define helper functions ------------------------------------------------------
#' 
create_comparison_table <- function(
    compiled_logs) {
  
  # Extract training models (one row per model)
  train_models <- compiled_logs %>%
    select(train_model = comparison, train_cohort = cohort, train_disease = disease)
  
  # Fill in validation sets (one row per cohort-disease pair)
  validation_sets <- compiled_logs %>%
    distinct(test_cohort = cohort, test_disease = disease)
  
  # Compile all pairwise comparisons
  comparison_table <- merge(train_models, validation_sets, by = NULL)
  
  # Add comparison type
  comparison_table <- comparison_table %>%
    mutate(
      comparison_type = case_when(
        train_cohort == test_cohort &
          train_disease != test_disease ~ "same_cohort_diff_disease",
        
        train_cohort != test_cohort &
          train_disease == test_disease ~ "cross_cohort_same_disease",
        
        train_cohort != test_cohort &
          train_disease != test_disease ~ "cross_cohort_diff_disease",
        
        TRUE ~ "internal"))
  
  return(comparison_table)
}

#' 
align_features <- function(
    train_obj, 
    test_obj) {
  
  # Extract matrices
  train_norm <- getNormalizedData(train_obj)
  test_norm  <- getNormalizedData(test_obj)
  
  # Find shared features
  common_features <- intersect(
    colnames(train_norm),
    colnames(test_norm))
  
  # Subset both objects
  train_norm <- train_norm[, common_features, drop = FALSE]
  test_norm  <- test_norm[, common_features, drop = FALSE]
  
  # Reassign
  setNormalizedData(train_obj) <- train_norm
  setNormalizedData(test_obj)  <- test_norm
  
  # Return both objects
  return(list(
    train_obj = train_obj,
    test_obj  = test_obj,
    features  = common_features))
}

log_processing <- function(
    processing_res,
    train_cohort,
    train_disease,
    test_cohort,
    test_disease,
    comparison_type,
    test_obj_qc) {
  
  if (is.null(processing_res$processing_log)) {
    processing_res$processing_log <- list()
  }
  
  outliers <- test_obj_qc$metadata$outliers
  
  row <- data.frame(
    train_cohort = train_cohort,
    train_disease = train_disease,
    test_cohort = test_cohort,
    test_disease = test_disease,
    comparison_type = comparison_type,
    test_outliers = if (length(outliers) == 0)
      NA_character_
    else
      paste(outliers, collapse = "; "),
    stringsAsFactors = FALSE
  )
  
  processing_res$processing_log[[length(processing_res$processing_log) + 1]] <- row
  
  return(processing_res)
}


#' Pipeline for processing a single external validation cohort
run_pipeline <- function(
    test_norm,
    test_cohort,
    cohort_rows,
    train_list,
    out_dir,
    processing_res) {
  
  # Initialise
  aligned_dir <- file.path(out_dir, "aligned")
  if (is.null(processing_res$aligned)) {
    processing_res$aligned <- list()
  }
  dir.create(aligned_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract unique diseases from cohort_rows
  diseases <- unique(cohort_rows$test_disease)
  diseases <- clean_condition(diseases)
  
  # Iterate per disease
  for (disease in diseases) {
    
    message("   Filtering disease: ", disease)
    
    test_obj <- FilterDataSet(
      test_norm,
      label = c("control", disease),
      condition_col = "study_condition")
    
    test_name <- paste0(test_cohort, "_", disease)
    
    # Evaluate all models
    for (model_name in names(train_list)) {
      
      message("      Aligning with: ", model_name)
      
      train_obj <- train_list[[model_name]]
      
      aligned <- align_features(
        train_obj = train_obj,
        test_obj  = test_obj)
      
      # QC test_obj
      test_obj_qc <- qc_obj(
        cohort_name = test_name,
        norm_obj    = aligned$test_obj,
        out_dir     = out_dir)
      
      # Log results
      comparison_row <- cohort_rows[
        cohort_rows$test_disease == disease &
          cohort_rows$train_model == sub("_model$", "", model_name),]
      
      processing_res <- log_processing(
        processing_res = processing_res,
        train_cohort = comparison_row$train_cohort,
        train_disease = comparison_row$train_disease,
        test_cohort = comparison_row$test_cohort,
        test_disease = comparison_row$test_disease,
        comparison_type = comparison_row$comparison_type,
        test_obj_qc = test_obj_qc)
      
      aligned$test_obj <- test_obj_qc$obj
      
      # Assign generalised "disease" tag to disease samples
      aligned$test_obj@SampleData$validation_condition <-
        ifelse(aligned$test_obj@SampleData$study_condition == disease,
          "disease", "control")
      
      aligned$train_obj@SampleData$validation_condition <-
        ifelse(aligned$train_obj@SampleData$study_condition == disease,
               "disease", "control")
      
      # Save aligned objs as RDS
      aligned_file <- file.path(
        aligned_dir,
        paste0(test_name, "_x_", model_name, "_aligned.rds"))
      
      saveRDS(aligned, aligned_file)
    }
  }
  
  return(processing_res)
}

# Define main functions --------------------------------------------------------
#'

# make for different types: external val
run_processing <- function(
    test_list,
    train_list,
    compiled_logs,
    out_dir){
  
  # Initialise storage object and output directory
  processing_res <- list()
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate comparison table and split into list of unique cohorts
  comparison_table <- create_comparison_table(compiled_logs) 
  processing_res$comparison_table <- as.data.frame(comparison_table) 
  
  cohort_split <- split(
    processing_res$comparison_table,
    processing_res$comparison_table$test_cohort)
  
  message("\nComparison table created")
  
  for (test_cohort in names(cohort_split)) {
    
    message("\nProcessing validation cohort: ", test_cohort)
    
    # Index test_list list using unique test cohorts
    key <- paste0(test_cohort, "_normalised")
    test_norm <- test_list[[key]]
    
    # Extract current cohort row
    cohort_rows <- cohort_split[[test_cohort]]
    
    processing_res <- run_pipeline(
      test_norm = test_norm,
      test_cohort = test_cohort,
      cohort_rows = cohort_rows,
      train_list = train_list,
      out_dir = out_dir,
      processing_res = processing_res)
    
    return(processing_res)
    
        
  }
  
  return(processing_res)
}

export_processing <- function(
    processing_res, 
    out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Write processing_res$comparison table
  write.csv(
    processing_res$comparison_table,
    file.path(out_dir, "comparison_table.csv"),
    row.names = FALSE
  )
}
  

# Load data --------------------------------------------------------------------
# test_list <- read_rds_files(file.path(prep_dir, "normalised"))
# train_list <- read_rds_files(file.path(train_dir, "models"))
# compiled_logs <- compile_logs(
#   file.path(train_dir, "logs"),
#   file_name = "compiled_logs.csv")


# Execute ----------------------------------------------------------------------
processing_res <- run_processing(test_list, train_list, compiled_logs, out_dir)
# export_processing(processing_res, out_dir)


aligned <- readRDS(
  file.path(
    out_dir, 
    "aligned/FengQ_2015_adenoma_x_FengQ_2015_adenoma_model_aligned.rds"))
