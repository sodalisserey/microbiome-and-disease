# Train ML model on a single processed comparison using CrcBiomeScreen package
# 1. Define output directory, functions and args
# 3. Execute main functions:
#   a. run_training()
#     * Read processed (normalised and partitioned) comparison object
#     * Validate training dataset sample sizes
#     * Get computing configuration
#     * Execute run_model() and log_model()
#   b. export_training()
#     * Save model as results/ml_training/models/[comparison]_model.rds
#     * Save log as results/ml_training/[comparison]_model_log.csv


# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
library(xgboost)
source("R/utils.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/05_ml_train"


# Define helper functions ------------------------------------------------------
#' Create an ML model based on the training subset of a cohort
run_model <- function(
    part_obj,
    cfg) {
  
  # Extract cfg
  num_cores <- cfg$num_cores
  n_cv <- cfg$n_cv
  
  # Extract metadata
  comparison <- part_obj$comparison
  disease <- part_obj$disease
  
  # Class weights enabled/disabled according to training data class imbalance
  class_weights <- if (part_obj$cl_imbalance_training == FALSE) {
    message("   Class weights disabled")
    FALSE
  } else {
    message("   Class weights enabled")
    TRUE
  }
  
  warnings <- character()
  error_msg <- NULL
  
  message ("   Training RF model")
  mod <- suppressMessages(
    withCallingHandlers(
      tryCatch({
        mod <- TrainModels(
          part_obj$obj,
          model_type = "RF",
          TaskName = paste0(comparison, "_RF"),
          ClassWeights = class_weights,
          TrueLabel = disease,
          num_cores = num_cores,
          n_cv = n_cv)

        mod <- EvaluateModel(
          mod,
          model_type = "RF",
          TaskName = paste0(comparison, "_RF_test"),
          TrueLabel = disease,
          PlotAUC = FALSE)

        mod

      }, error = function(e) {
        error_msg <<- e$message
        message("RF failed: ", e$message)
        NULL
      }),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }))

  message ("   Training XGBoost model")
  mod <- suppressMessages(
    withCallingHandlers(
      tryCatch({
        mod <- TrainModels(
          mod,
          model_type = "XGBoost",
          TaskName = paste0(comparison, "_XGB"),
          ClassWeights = class_weights,
          TrueLabel = disease,
          num_cores = num_cores,
          n_cv = n_cv)

        mod <- EvaluateModel(
          mod,
          model_type = "XGBoost",
          TaskName = paste0(comparison, "_XGB_test"),
          TrueLabel = disease,
          PlotAUC = FALSE)

        mod

      }, error = function(e) {
        error_msg <<- e$message
        message("XGBoost failed: ", e$message)
        NULL
      }),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }))
  
  return(list(
    model = mod,
    cohort_name = part_obj$cohort_name,
    comparison = comparison,
    disease = disease,
    
    partition_ratio = part_obj$partition_ratio,
    n_threshold = part_obj$n_threshold,
    cl_ratio_training = part_obj$cl_ratio_training,
    cl_imbalance_training = part_obj$cl_imbalance_training,
    n_controls_training = part_obj$n_controls_training,
    n_disease_training = part_obj$n_disease_training,
    n_status_training = part_obj$n_status_training,
    
    class_weights = class_weights,
    warnings = unique(warnings),
    error = error_msg,
    status = if (is.null(error_msg)) "SUCCESS" else "FAILED"))
}

#' Log model status for a single comparison, append raining_res$training_log 
#' and return training_res
log_model <- function(
    result) {
  
  create_log <- function(status, reason, AUC_RF = NA_real_, AUC_XGB = NA_real_) {
    data.frame(
      comparison = result$comparison,
      cohort = result$cohort_name,
      disease = result$disease,
      status = status,
      reason = reason,
      AUC_RF = AUC_RF,
      AUC_XGB = AUC_XGB,
      error = if (length(result$error) == 0) NA_character_ else 
        paste(result$error, collapse = "| "),
      warnings = if (length(result$warnings) == 0) NA_character_ else
        paste(result$warnings, collapse = " | "),
      class_weights = if (is.null(result$class_weights)) NA else
        result$class_weights,
      stringsAsFactors = FALSE)
  }
  
  # If null result
  if (is.null(result) || is.null(result$model)) {
    return(create_log("FAILED", "model is NULL"))
  } 
  
  # If success
  AUC_RF <- extract_auc(result$model@EvaluateResult$RF$AUC)
  AUC_XGB <- extract_auc(result$model@EvaluateResult$XGBoost$AUC)
  
  message(paste0("   > RF AUC = ", AUC_RF))
  message(paste0("   > XGB AUC = ", AUC_XGB))
  
  create_log("SUCCESS", NA_character_, AUC_RF, AUC_XGB)
}

# Define main functions --------------------------------------------------------
#' Run training for a single processed object
run_training <- function(file, out_dir) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Loading: ", basename(file))
  
  # Read partitioned object
  part_obj <- readRDS(file)
  
  # Check n_training_status
  if (part_obj$n_status_training != "OK") {
    message("Skipping ", part_obj$comparison, ": ", part_obj$n_status_training)
    return(NULL)
  }
  
  # Get configuration
  cfg <- get_config()  
  
  # Create and log model
  result <- run_model(
    part_obj = part_obj,
    cfg = cfg)
  
  log_df <- log_model(result)
  
  # Save model
  mod_dir <- file.path(out_dir, "models")
  dir.create(mod_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (result$status == "SUCCESS") {
    saveRDS(
      result$model,
      file.path(mod_dir, paste0(result$comparison, "_model.rds")))
  }
  
  # Save log
  log_dir <- file.path(out_dir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    log_df,
    file.path(log_dir, paste0(result$comparison, "_model_log.csv")),
    row.names = FALSE)
  
  message("\n   Done. Model written to: ", mod_dir)
}


# Define args ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

if (length(args) != 1) {
  stop("Usage: Rscript ml_train.R <processed.rds>")
}


# Execute ----------------------------------------------------------------------
run_training(file, out_dir)
