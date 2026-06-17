# Process ML models for training using CrcBiomeScreen package

# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
source("R/utils.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/04_ml_process"


# Define helper functions ------------------------------------------------------
#' normalisedare cohort by creating CrcBiomeScreenObject, splitting/setting taxa, 
#' normalising data, filtering study_condition, running QC, splitting dataset 
#' into training/testing and checking class balance
run_qc <- function(
    cohort,
    cohort_name,
    out_dir,
    norm_method = "GMPR") {
  
  qc_dir <- file.path(out_dir, "plots")
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  
  warnings <- character()
  
  result <- withCallingHandlers({
    obj <- CreateCrcBiomeScreenObject(
      RelativeAbundance = cohort@assays@data@listData$relative_abundance,
      TaxaData = cohort@rowLinks$nodeLab,
      SampleData = cohort@colData)
    
    # Split taxa and keep genus level
    obj <- SplitTaxas(obj)
    taxa <- getTaxaData(obj)
    colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    setTaxaData(obj) <- taxa
    obj <- KeepTaxonomicLevel(obj, level = "Genus")
    
    # Normalise data
    obj <- NormalizeData(obj, method = "GMPR", level = "Genus")
    
    # Standardise labels
    obj@SampleData$study_condition <-
      clean_label(obj@SampleData$study_condition)
    
    # Perform QC
    suppressMessages(
      obj <- qcByCmdscale(
        obj,
        TaskName = paste0(cohort_name, "_QC"),
        outdir = qc_dir,
        normalize_method = norm_method,
        plot = TRUE))

    file.rename(
      file.path(qc_dir, paste0("cmdscale_", cohort_name, "_QC_", norm_method, ".pdf")),
      file.path(qc_dir, paste0(cohort_name, "_QC_", norm_method, ".pdf")))
    
  },
  warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  
  return(list(
    obj = obj,
    cohort_name = cohort_name,
    norm_method = norm_method,
    warnings = if (length(warnings) == 0) NA_character_ else
      paste(warnings, collapse = " | "),
    
    n_features = ncol(obj@OriginalNormalizedData),
    n_samples = nrow(obj@OriginalNormalizedData),
    
    outliers = obj@OutlierSamples,
    n_features_qc = ncol(obj@NormalizedData),
    n_samples_qc = nrow(obj@NormalizedData)))
}

#' Log QC results from a single cohort_name containing number of features and 
#' samples before and after QC and outliers, then append processing_res$qc_summary 
#' and return processing_res
log_qc <- function(
    cohort_name, 
    qc_obj,
    processing_res) {
  
  # initialise containers if missing
  if (is.null(processing_res$qc_summary)) {
    processing_res$qc_summary <- list()
  }
  
  row <- data.frame(
    cohort_name = cohort_name,
    norm_method = qc_obj$norm_method,
    warnings = if (length(qc_obj$warnings) == 0) NA_character_ else
      paste(qc_obj$warnings, collapse = " | "),
    
    # N features (genus level) and samples pre-normalisation/qc
    n_features = qc_obj$n_features,
    n_samples = qc_obj$n_samples,
    
    # Number of outliers identified by QC
    outliers = if (length(qc_obj$outliers) == 0) NA_character_ else 
      paste(qc_obj$outliers, collapse = "; "),
    
    # N features (genus level) and samples post-normalisation/qc
    n_features_qc = qc_obj$n_features_qc,
    n_samples_qc = qc_obj$n_samples_qc,
    
    stringsAsFactors = FALSE)
  
  processing_res$qc_summary[[cohort_name]] <- row
  
  return(processing_res)
}

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
    cohort_name = cohort_name,
    comparison = comparison,
    disease = disease,
    partition_ratio = partition_ratio,
    n_threshold = n_threshold,
    cl_ratio_training = max(balance$class_counts)/min(balance$class_counts),
    cl_imbalance_training = balance$is_imbalanced,
    n_controls_training = n_controls_training,
    n_disease_training = n_disease_training,
    n_status_training = n_status_training))
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
    
    partition_ratio = part_obj$partition_ratio,
    
    cl_ratio_training = part_obj$cl_ratio_training,
    cl_imbalance_training = part_obj$cl_imbalance_training,
    
    n_controls_training = part_obj$n_controls_training,
    n_disease_training = part_obj$n_disease_training,
    
    n_status_training = part_obj$n_status_training,
    n_threshold = part_obj$n_threshold,

    
    stringsAsFactors = FALSE)
  
  processing_res$partition_summary[[comparison]] <- row
  
  return(processing_res)
}

#'
run_pipeline <- function(
    qc_obj,
    cohort_name,
    disease,
    out_dir,
    processing_res) {
  
  comp_dir <- file.path(out_dir, "processed")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  
  comparison <- paste(cohort_name, disease, sep = "_")
  
  message("Processing: ", comparison)
  
  # Filter for single disease
  suppressWarnings(
    comparison_obj <- FilterDataSet(
      qc_obj$obj,
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
  if (part_obj$n_status_training == "OK") {
    saveRDS(
      part_obj,
      file.path(comp_dir, paste0(comparison, "_processed.rds"))
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
    
    # Run QC on full cohort
    qc_obj <- run_qc(
      cohort = cohort,
      cohort_name = cohort_name,
      out_dir = out_dir)
    
    processing_res <- log_qc(
      cohort_name = cohort_name,
      qc_obj = qc_obj,
      processing_res = processing_res)
    
    diseases <- clean_label(info$diseases)
    
    for (disease in diseases) {
      
      processing_res <- run_pipeline(
        qc_obj = qc_obj,
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
}
  

# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary_cohorts.rds")


# Execute ----------------------------------------------------------------------
processing_res <- run_processing(primary_cohorts, out_dir)
export_processing(processing_res, out_dir)
