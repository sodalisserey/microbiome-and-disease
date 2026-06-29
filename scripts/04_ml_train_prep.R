# Process ML models for training using CrcBiomeScreen package
# 1. Define output directory and functions
# 2. Execute main functions:
#   a. run_processing()
#     * Process, extract and validate study_conditions
#     * Execute R/utils.R/normalise_obj(), qc_obj() and log_qc() on each cohort
#     * Execute run_pipeline() calling run_partition() and log_partition()
#   b. export_processing()
#     * Save processed object as out_dir/processed/[comparison]_processed.rds
#     * Save qc_summary and partition_summary as out_dir/[summary_type].csv


# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
source("R/utils.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/04_ml_train_prep_test"


# Define helper functions ------------------------------------------------------
#' Split a prepped cohort into training and testing datasets and return norm_obj
run_partition <- function(
    cohort_name,
    comparison,
    comparison_obj,
    disease,
    partition_ratio = 0.7,
    n_threshold = 20,
    out_dir) {
  
  # Partition into datasets
  suppressWarnings(  
    comparison_obj <- SplitDataSet(
    comparison_obj,
    label = c("control", disease),
    partition = partition_ratio))
  
  # Check class balance in training subset
  suppressMessages(
    balance <- checkClassBalance(
      getModelData(comparison_obj)$TrainLabel, 
      outdir = out_dir, 
      plot = FALSE))

  # Get counts from training subset
  class_counts <- table(balance$class_counts)
  
  n_controls_training <- if (
    !is.null(balance$class_counts[["control"]])) balance$class_counts[["control"]] else 0
  n_disease_training  <- if (
    !is.null(balance$class_counts[[disease]])) balance$class_counts[[disease]] else 0
  
  n_status_training <- character()
  
  if (n_controls_training < n_threshold) {
    n_status_training <- c(n_status_training, paste0("controls < ", n_threshold))
  }
  
  if (n_disease_training < n_threshold) {
    n_status_training <- c(n_status_training, paste0(disease, " < ", n_threshold))
  }
  
  n_status_training <- if (length(n_status_training) == 0) {
    "OK"
  } else {
    paste(n_status_training, collapse = " | ")
  }     
  
  return(list(
    obj = comparison_obj,
    metadata = list(
      cohort_name = cohort_name,
      comparison = comparison,
      disease = disease,
      partition_ratio = partition_ratio,
      n_threshold = n_threshold,
      cl_ratio_training =
        max(balance$class_counts) / min(balance$class_counts),
      cl_imbalance_training = balance$is_imbalanced,
      n_controls_training = n_controls_training,
      n_disease_training = n_disease_training,
      n_status_training = n_status_training)))
}

#' Log partition results containing training dataset class balance and sample 
#' status, export suitable norm_obj.rds, then append processing_res$partition_summary
#'  and return processing_res
log_partition <- function(
    cohort_name, 
    comparison,
    disease,
    part_obj,
    processing_res) {
  
  # initialise containers if missing
  if (is.null(processing_res$partition_summary)) {
    processing_res$partition_summary <- list()
  }
  
  row <- data.frame(
    cohort_name = cohort_name,
    disease = disease,
    
    partition_ratio = part_obj$metadata$partition_ratio,
    
    cl_ratio_training = part_obj$metadata$cl_ratio_training,
    cl_imbalance_training = part_obj$metadata$cl_imbalance_training,
    
    n_controls_training = part_obj$metadata$n_controls_training,
    n_disease_training = part_obj$metadata$n_disease_training,
    
    n_status_training = part_obj$metadata$n_status_training,
    n_threshold = part_obj$metadata$n_threshold,

    
    stringsAsFactors = FALSE)
  
  processing_res$partition_summary[[comparison]] <- row
  
  return(processing_res)
}

#'
run_pipeline <- function(
    processed_obj,
    cohort_name,
    disease,
    out_dir,
    processing_res) {
  
  prep_dir <- file.path(out_dir, "processed")
  dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)
  
  comparison <- paste(cohort_name, disease, sep = "_")
  
  message("Partitioning: ", comparison)
  
  # Filter for single disease
  suppressWarnings(
    comparison_obj <- FilterDataSet(
      processed_obj$obj,
      label = c("control", disease),
      condition_col = "study_condition"))
  
  # Partition into training dataset
  part_obj <- run_partition(
    cohort_name = cohort_name,
    comparison = comparison,
    disease = disease,
    comparison_obj = comparison_obj,
    out_dir = out_dir)
  
  # Log partition
  processing_res <- log_partition(
    cohort_name = cohort_name,
    comparison = comparison,
    disease = disease,
    part_obj = part_obj,
    processing_res = processing_res)
  
  # Save as norm_obj if suitable
  if (part_obj$metadata$n_status_training == "OK") {
    saveRDS(
      part_obj,
      file.path(prep_dir, paste0(comparison, "_processed.rds"))
    )
  }
  return(processing_res)
}
  
# Define main functions --------------------------------------------------------
#' Run main ML normalisedaration for cohorts
run_processing <- function(cohort_list, out_dir) {
  
  # Initialise storage object
  processing_res <- list()
  
  for (cohort_name in names(cohort_list)) {
    
    # Process and extract study_conditions from cohort
    cohort <- cohort_list[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check cohort contains both healthy and disease samples
    if (!is.null(validate_conditions(info, cohort_name))) next
  
    message("\nRunning QC: ", cohort_name)
    
    # Process (normalise and QC) whole cohort
    norm_obj <- normalise_obj(
      cohort = cohort,
      cohort_name = cohort_name)
    
    processed_obj <- qc_obj(
      cohort = cohort,
      cohort_name = cohort_name,
      obj = norm_obj$obj,
      norm_method = norm_obj$norm_method,
      warnings = norm_obj$warnings,
      out_dir = out_dir)
    
    processing_res <- log_qc_obj(
      cohort_name = cohort_name, 
      processed_obj = processed_obj,
      processing_res = processing_res)
    
    diseases <- clean_condition(info$diseases)
    
    # Run partition pipeline for each comparison
    for (disease in diseases) {
      
      processing_res <- run_pipeline(
        processed_obj = processed_obj,
        cohort_name = cohort_name,
        disease = disease,
        out_dir = out_dir,
        processing_res = processing_res)
    }
  }
  return(processing_res)
}

#' Export processing outputs as raw R objects and CSVs
export_processing <- function(processing_res, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Bind data frames
  processing_res$qc_summary <- dplyr::bind_rows(processing_res$qc_summary)
  processing_res$partition_summary <- dplyr::bind_rows(processing_res$partition_summary)
  
  # Write CSV
  write_csvs(out_dir, list(
    qc_summary = processing_res$qc_summary,
    partition_summary = processing_res$partition_summary))
  
  # Write txt file of all processed file paths
  files <- list.files(
    path = file.path(out_dir, "processed"),
    pattern = "_processed\\.rds$",
    full.names = TRUE)
  
  writeLines(
    files, file.path(out_dir,"processed_files.txt"))
}
  

# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary.rds") |>
  create_cohort_obj(cohorts = _)

# Execute ----------------------------------------------------------------------
processing_res <- run_processing(primary_cohorts, out_dir)
export_processing(processing_res, out_dir)
