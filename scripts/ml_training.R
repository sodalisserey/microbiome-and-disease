# Train ML models using CrcBiomeScreen package
# 1. Define output directory and functions
# 2. Load data from cohorts.R
# 3. Execute main functions:
#   a. run_training()
#     * Clean, extract disease and validate conditions
#     * TRAINING BRANCH 1 (healthy x 1 disease) -> run_pipeline() executes:
#           - run_qc() 
#           - run_partition()
#           - create_model() if n_training_status == "OK"
#           - log_qc()
#           - log_training()
#           - log_result()
#     * TRAINING BRANCH 2 (healthy x >1 disease)  -> 
#           - Create pair-wise subsets
#           - Execute run_pipeline()
#   b. export_training()
#     * Bulk save analysis results into:
#           - results/ml_training/qc_summary.csv
#           - results/ml_training/training_log.csv
#           - results/ml_training/models/[comparison]_RF.rds
#           - results/ml_training/all_models.rds


# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
source("R/utils.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/ml_training"


# Define helper functions ------------------------------------------------------
#' Prepare cohort by creating CrcBiomeScreenObject, splitting/setting taxa, 
#' normalising data, filtering study_condition, running QC, splitting dataset 
#' into training/testing and checking class balance
run_qc <- function(
    cohort,
    comparison_name,
    disease,
    out_dir,
    norm_method = "GMPR") {
  
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
      
      # Filter for control and disease
      obj <- FilterDataSet(
        obj,
        label = c("control", disease),
        condition_col = "study_condition")
      
      # Perform QC
      suppressMessages(
        obj <- qcByCmdscale(
          obj,
          TaskName = paste0(comparison_name, "_QC"),
          outdir = out_dir,
          normalize_method = norm_method,
          plot = TRUE))
      
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  
  return(list(
    obj = obj,
    comparison = comparison_name,
    disease = disease,
    norm_method = norm_method,
    warnings = if (length(warnings) == 0) NA_character_ else
      paste(warnings, collapse = " | "),
    
    n_features = ncol(obj@OriginalNormalizedData),
    n_samples = nrow(obj@OriginalNormalizedData),
    
    outliers = obj@OutlierSamples,
    n_features_qc = ncol(obj@NormalizedData),
    n_samples_qc = nrow(obj@NormalizedData)))
}

#' Split a prepped cohort into training and testing datasets and return results
run_partition <- function(
    prep,
    partition_ratio = 0.7,
    n_threshold = 20,
    out_dir) {
  
  disease <- prep$disease

  # Partition into datasets
  obj <- SplitDataSet(
    prep$obj,
    label = c("control", disease),
    partition = partition_ratio)
  
  # Check class balance in training subset
  suppressMessages(
    balance <- checkClassBalance(
      getModelData(obj)$TrainLabel, 
      outdir = out_dir, 
      plot = FALSE))
  
  # Get counts from training subset
  class_counts <- table(balance$class_counts)
  
  n_controls_training = (balance$class_counts[["control"]] %||% 0)
  n_disease_training = (balance$class_counts[[disease]] %||% 0)
  
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
  
  prep$obj <- obj
  
  prep$partition_ratio <- partition_ratio
  prep$n_threshold <- n_threshold
  
  prep$cl_ratio_training <- max(balance$class_counts)/min(balance$class_counts)
  
  prep$cl_imbalance_training <- balance$is_imbalanced
  
  prep$n_controls_training <- n_controls_training
  
  prep$n_disease_training <- n_disease_training
  
  prep$n_status_training <- n_status_training
  
  prep
}

#' Create an ML model based on the training subset of a cohort
create_model <- function(
    prep,
    disease,
    comparison,
    out_dir,
    cfg) {
  
  safe_comp <- gsub("[/\\\\]", "", comparison)
  
  model_dir <- file.path(out_dir, "models")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Class weights enabled/disabled according to training data class imbalance
  class_weights <- if (prep$cl_imbalance_training == FALSE) {
    message("Class weights disabled")
    FALSE
  } else {
    message("Class weights enabled")
    TRUE
  }
  
  warnings <- character()
  error_msg <- NULL
  
  num_cores <- cfg$num_cores
  n_cv <- cfg$n_cv
  
  rf_model <- suppressMessages(
    withCallingHandlers(
      tryCatch({
        rf_model <- TrainModels(
          prep$obj,
          model_type = "RF",
          TaskName = paste0(safe_comp, "_RF"),
          ClassWeights = class_weights,
          TrueLabel = disease,
          num_cores = num_cores,
          n_cv = n_cv)
        
        rf_model <- EvaluateModel(
          rf_model,
          model_type = "RF",
          outdir = model_dir,
          TaskName = paste0(safe_comp, "_RF_test"),
          TrueLabel = disease,
          PlotAUC = FALSE)
        
        file.rename(
          file.path(model_dir, paste0("CrcBiomeScreenObject_", safe_comp, "_RF_test.rds")),
          file.path(model_dir, paste0(safe_comp, "_RF.rds")))
        
        rf_model
      
      }, error = function(e) {
        error_msg <<- e$message
        message("RF failed: ", e$message)
        NULL
      }),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }))
        
  # TODO XGBoost model
  return(list(
    rf_model = rf_model,
    class_weights = class_weights,
    warnings = unique(warnings),
    error = error_msg,
    status = if (is.null(error_msg)) "SUCCESS" else "FAILED"))
}

#' Log QC results from a single comparison containing number of features and 
#' samples before and after QC, outliers and training dataset class balance,
#' then append training_res$qc_summary and return training_res
log_qc <- function(
    comparison_name, 
    cohort_name,
    prep,
    training_res) {
  
  # initialise container if missing
  if (is.null(training_res$qc_summary)) {
    training_res$qc_summary <- list()
  }
  
  row <- data.frame(
    comparison = comparison_name,
    cohort = cohort_name,
    disease = prep$disease,
    norm_method = prep$norm_method,
    warnings = if (length(prep$warnings) == 0) NA_character_ else
      paste(prep$warnings, collapse = " | "),
    
    # N features (genus level) and samples pre-normalisation/qc
    n_features = prep$n_features,
    n_samples = prep$n_samples,
    
    # Number of outliers identified by QC
    outliers = if (length(prep$outliers) == 0) NA_character_ else 
      paste(prep$outliers, collapse = "; "),
    
    # N features (genus level) and samples post-normalisation/qc
    n_features_qc = prep$n_features_qc,
    n_samples_qc = prep$n_samples_qc,
    
    # Training dataset partitioning
    partition_ratio = prep$partition_ratio,
    cl_ratio_training = prep$cl_ratio_training,
    cl_imbalance_training = prep$cl_imbalance_training,
    n_controls_training = prep$n_controls_training,
    n_disease_training = prep$n_disease_training,
    n_status_training = prep$n_status_training,
    n_threshold = prep$n_threshold,
    
    stringsAsFactors = FALSE)
  
  training_res$qc_summary[[comparison_name]] <- row
  
  return(training_res)
}

#' Log model training status for a single comparison, append 
#' training_res$training_log and return training_res
log_training <- function(
    comparison_name, 
    cohort_name, 
    disease, 
    status = NULL,
    reason = NULL,
    model = list(error = NA, warnings = NA, class_weights = NA),
    training_res) {
  
  # Initialise container if missing
  if (is.null(training_res$training_log)) {
    training_res$training_log <- list()
  }
  
  create_log <- function(status, reason, AUC = NA_real_) {
    data.frame(
      comparison = comparison_name,
      cohort = cohort_name,
      disease = disease,
      status = status,
      reason = reason,
      internal_AUC = AUC,
      error = if (length(model$error) == 0) NA_character_ else 
        paste(model$error, collapse = "| "),
      warnings = if (length(model$warnings) == 0) NA_character_ else
        paste(model$warnings, collapse = " | "),
      class_weights = if (is.null(model$class_weights)) NA else
        model$class_weights,
      stringsAsFactors = FALSE)
  }
  
  # Case 1: explicit SKIPPED when condition_invalid == TRUE
  if (!is.null(status) && status == "SKIPPED") {
    training_res$training_log[[comparison_name]] <-
      create_log("SKIPPED", reason, NA_real_)
    return(training_res)
  }
  
  # Case 2: null result
  # TODO this is not triggered when models aren't created. try with model$rf_model?
  if (is.null(model) || is.null(model$rf_model)) {
    training_res$training_log[[comparison_name]] <-
      create_log("FAILED", "training returned NULL", NA_real_)
    return(training_res)
  
  # Case 3: success
  } else {
    AUC <- tryCatch(
      as.numeric(sub(".*: ", "", model$rf_model@EvaluateResult$RF$AUC)),
      error = function(e) NA_real_
    )
    
    message(paste0("AUC = ", AUC))
    
    training_res$training_log[[comparison_name]] <-
      create_log("SUCCESS", NA_character_, AUC)
  }
  
  return(training_res)
}

#' Log model from a single training, append training_res$all_models and return 
#' training_res
log_result <- function(
    comparison_name, 
    model,
    training_res){
  
  # initialise container if missing
  if (is.null(training_res$all_models)) {
    training_res$all_models <- list()
  }
  
  if (!is.null(model$rf_model)) {
    training_res$all_models[[comparison_name]] <- model$rf_model
  }
  
  return(training_res)
}

#' Run training pipeline which includes preprocessing and QC checks, RF modellng
#' and logging
run_pipeline <- function(
    cohort,
    cohort_name,
    comparison_name,
    disease,
    training_res){
  
  # Prepare cohort for training and log QC results
  prep <- run_qc(
    cohort = cohort,
    comparison_name = comparison_name,
    disease = disease,
    out_dir = out_dir)
  
  prep <- run_partition(
    prep = prep,
    out_dir = out_dir)
  
  training_res <- log_qc(
    comparison_name, 
    cohort_name = cohort_name,
    prep = prep,
    training_res = training_res)
  
  message(paste0("N training samples = ", prep$n_status_training))
  
  # Skip training if n_sample too small
  if (prep$n_status_training != "OK") {
    
    training_res <- log_training(
      comparison_name,
      cohort_name = cohort_name,
      disease = disease,
      status = "SKIPPED",
      reason = prep$n_status_training,
      model = list(
        error = NA_character_,
        warnings = NA_character_,
        class_weights = NA),
      training_res = training_res)
    
    return(training_res)
  }
    
  # Train model on cohort, log training info/status and save results
  cfg <- get_config()
  
  model <- create_model(
    prep,
    disease = disease,
    comparison = comparison_name,
    out_dir = out_dir,
    cfg = cfg)
  
  training_res <- log_training(
    comparison_name, 
    cohort_name = cohort_name, 
    disease = disease, 
    model = model,
    training_res = training_res) 
  
  training_res <- log_result(
    comparison_name,
    model = model,
    training_res = training_res)
  
  return(training_res)
}


# Define main functions --------------------------------------------------------
#' Run main ML modelling loop for two training branches
run_training <- function(cohort_list, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialise storage object
  training_res <- list()
  
  for (cohort_name in names(cohort_list)) {
    
    # Process and extract study_conditions from cohort
    cohort <- cohort_list[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check cohort contains both healthy and disease samples
    condition_invalid <- validate_conditions(info, cohort_name)
    
    if (!is.null(condition_invalid)) {
      training_res <- log_training(
        comparison_name = NA_character_,
        cohort_name,
        disease = NA_character_,
        status = "SKIPPED",
        reason = condition_invalid,
        training_res = training_res)
      next  
    }
      
    # TRAINING BRANCH 1: healthy x 1 disease
    if (info$n_diseases == 1) {
      
      disease = info$diseases
    
      comparison_name <- paste(cohort_name, disease, sep = "_")
      message("\n", cohort_name, ": healthy vs ", info$diseases)
      
      training_res <- run_pipeline(
        cohort = cohort,
        cohort_name = cohort_name,
        disease = disease,
        comparison_name = comparison_name,
        training_res = training_res)
      
    # TRAINING BRANCH 2: healthy x >1 disease
    } else {
      message("\n", cohort_name, ": healthy vs ", info$n_diseases, " diseases")

      for (disease in info$diseases) {
        message("\n", cohort_name, " subset: healthy vs ", disease)

        # Subset cohort for healthy x 1 disease
        meta <- as.data.frame(colData(cohort))
        keep <- meta$study_condition %in% c("control", disease)
        cohort_subset <- cohort[, keep]

        comparison_name <- paste(cohort_name, disease, sep = "_")
        
        training_res <- run_pipeline(
          cohort = cohort_subset,
          cohort_name = cohort_name,
          disease = disease,
          comparison_name = comparison_name,
          training_res = training_res)
      }
    }
  }
  return(training_res)
}

#' Export training outputs as raw R objects and CSVs
export_training <- function(training_res, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Bulk save all models
  saveRDS(training_res$all_models, file.path(out_dir, "all_models.rds"))
  
  # Bind data frames
  training_res$qc_summary <- dplyr::bind_rows(training_res$qc_summary)
  training_res$training_log <- dplyr::bind_rows(training_res$training_log)
  
  # Write CSVs
  write_csvs(
    out_dir,
    list(
      qc_summary = training_res$qc_summary,
      training_log = training_res$training_log))
}


# Load data --------------------------------------------------------------------
test_cohorts <- readRDS("data/test_cohorts.rds")

# primary_cohorts <- readRDS("data/primary_cohorts.rds")


# Execute ----------------------------------------------------------------------
training_res <- run_training(test_cohorts, out_dir)
export_training(training_res, out_dir)
