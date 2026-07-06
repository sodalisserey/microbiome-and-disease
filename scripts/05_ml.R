# Train ML models from prepped comparisons using CrcBiomeScreen package
# 1. Define output directory and functions
# 2. Load comparison_table.csv and define command line args
# 3. Execute main function run_modelling() which runs:
#     a. get_configuration()
#         * returns num_cores and n_cv according to computing configuration
#     b. partition_training()
#         * splits cohort object for training and internal validation
#         * checks for sufficient training sample size (>= 20) and class balance
#         * returns partitioned training object
#     c. train_model() 
#         * enables/disables class weights
#         * models partitioned training object using RF and XGBoost
#         * runs extract_auc() and prints values
#         * saves as out_dir/models/[comparison]_model.rds
#     d. validate_model() 
#         * filters for current validation rows from comparison_table
#         * for each row, load validation objects from in_dir/final
#         * externally validates trained objects using RF and XGBoost
#         * runs extract_auc() and prints values
#         * saves as out_dir/final/[comparison]_evaluated.rds
#     e. export_logs()
#         * combines results from partition_training() and train_model() and save
#           as out_dir/training_log/[comparison].csv
#         * combines results from validate_model() and save as out_dir/validation_
#           log/[comparison].csv


# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
library(dplyr)
source("R/utils.R")


# Define input/output directories ----------------------------------------------
in_dir <- "results/04_pre_ml"
out_dir <- "results/05_ml_test"


# Define helper functions ------------------------------------------------------
# Get computing configuration and return num_cores and n_cv
get_configuration <- function() {
  
  is_slurm <- nzchar(Sys.getenv("SLURM_JOB_ID"))
  
  if (is_slurm) {
    cfg <- list(
      num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
      n_cv = 10,
      mode = "slurm")
    message("\nConfiguration: Slurm (num_cores = ", 
            cfg$num_cores, " | n_cv = ", cfg$n_cv, ")")
    
  } else {
    cfg <- list(
      num_cores = 1,
      n_cv = 2,
      mode = "local")
    message("\nConfiguration: local (num_cores = ", 
            cfg$num_cores, " | n_cv = ", cfg$n_cv, ")")
    
  }
  return(cfg)
}

# Split cohort object for training and internal validation, extract class
# balance and class counts and check for sufficient sample size (>=20)
partition_training <- function(
    train_obj,
    partition_ratio = 0.7,
    partition_threshold = 20) {
  
  message("\nPartitioning for training")
  
  obj <- train_obj$obj
  meta <- train_obj$metadata

  # Partition for internal validation
  suppressWarnings(  
    obj <- SplitDataSet(
      obj,
      label = c("control", "disease"),
      partition = partition_ratio,
      condition_col = "validation_condition"))
  
  # Check class balance of training subset
  suppressMessages(
    balance <- checkClassBalance(
      getModelData(obj)$TrainLabel, 
      plot = FALSE))
  
  class_weights <- balance$is_imbalanced
  
  message("   > class imbalance = ", class_weights)
  
  # Extract counts and check sample size
  class_counts <- table(balance$class_counts)
  
  part_n_control <- if (
    !is.null(balance$class_counts[["control"]])) balance$class_counts[["control"]] else 0
  part_n_disease  <- if (
    !is.null(balance$class_counts[["disease"]])) balance$class_counts[["disease"]] else 0
  
  n_status <- character()
  
  if (part_n_control < partition_threshold) {
    n_status <- c(n_status, paste0("controls < ", partition_threshold))
  }
  
  if (part_n_disease < partition_threshold) {
    n_status <- c(n_status, paste0(train_disease, " < ", partition_threshold))
  }
  
  n_status <- if (length(n_status) == 0) {
    "OK"
  } else {
    paste(n_status, collapse = " | ")
  }     
  
  message("   > sample size = ", n_status, 
          " (controls = ", part_n_control, " | disease = ", part_n_disease, ")")
  
  # Compile results
  partition <- list(
    partition_ratio = partition_ratio,
    partition_threshold = partition_threshold,
    n_status = n_status,
    part_n_control = part_n_control,
    part_n_disease = part_n_disease,
    class_weights = class_weights)
  
  return(list(
    obj = obj,
    metadata = meta,
    partition = partition))
}

# Extract, print and return numeric AUC value from trained/validated object
extract_auc <- function(auc_val) {
  
  auc <- suppressWarnings(as.numeric(sub(".*: ", "", auc_val)))
  if (is.na(auc))
    auc <- NA_real_
  message("   AUC = ", auc)
  auc
}

# Train partitioned object using RF and XGBoost, save as RDS and return object
train_model <- function(
    train_obj,
    cfg) {
  
  # Define model path
  mod_dir <- file.path(out_dir, "models")
  dir.create(mod_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract obj items
  obj <- train_obj$obj
  meta <- train_obj$metadata
  part <- train_obj$partition
  comparison <- meta$comparison
  disease <- meta$disease
  
  # Confirm if model already exists
  mod_file <- file.path(
    mod_dir,
    paste0(comparison, "_model.rds"))
  
  # Skip if model already exists
  if (file.exists(mod_file)) {
    message("\nSkipping model training")
    message("   > already exists at: ", mod_file)
    
    train_obj <- readRDS(mod_file)

    return(train_obj)
  }
  
  message("\nTraining partitioned subset")
  
  # Extract cfg
  num_cores <- cfg$num_cores
  n_cv <- cfg$n_cv
  
  # Class weights enabled/disabled
  class_weights <- if (part$class_weights == FALSE) {
    message("   > class weights disabled")
    FALSE
  } else {
    message("   > class weights enabled")
    TRUE
  }
  
  warnings <- character()
  error_msg <- character()
  
  message ("   > RF model training")
  mod <- suppressMessages(
    withCallingHandlers(
      tryCatch({
        mod <- TrainModels(
          obj,
          model_type = "RF",
          TaskName = paste0(comparison, "_RF"),
          ClassWeights = class_weights,
          TrueLabel = "disease",
          num_cores = num_cores,
          n_cv = n_cv)
        
        mod <- EvaluateModel(
          mod,
          model_type = "RF",
          TaskName = paste0(comparison, "_RF_test"),
          TrueLabel = "disease",
          PlotAUC = FALSE)
        
        mod
        
      }, error = function(e) {
        error_msg <<- e$message
        message("RF training failed: ", e$message)
        NULL
      }),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }))
  
  AUC_RF <- if (length(mod@EvaluateResult$RF$AUC) == 0) { NA_real_
  } else {
    extract_auc(mod@EvaluateResult$RF$AUC)
  }
  
  message ("   > XGBoost model training")
  mod <- suppressMessages(
    withCallingHandlers(
      tryCatch({
        mod <- TrainModels(
          mod,
          model_type = "XGBoost",
          TaskName = paste0(comparison, "_XGB"),
          ClassWeights = class_weights,
          TrueLabel = "disease",
          num_cores = num_cores,
          n_cv = n_cv)
        
        mod <- EvaluateModel(
          mod,
          model_type = "XGBoost",
          TaskName = paste0(comparison, "_XGB_test"),
          TrueLabel = "disease",
          PlotAUC = FALSE)
        
        mod
        
      }, error = function(e) {
        error_msg <<- e$message
        message("XGB training failed: ", e$message)
        NULL
      }),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }))
  
  AUC_XGB <- if (length(mod@EvaluateResult$XGBoost$AUC) == 0) { NA_real_
  } else {
    extract_auc(mod@EvaluateResult$XGBoost$AUC)
  }
  
  # Compile results
  training <- list(
    AUC_RF = AUC_RF,
    AUC_XGB = AUC_XGB,
    warnings = if (length(warnings) == 0) NA_character_ else
      paste(unique(warnings), collapse = " | "),
    error = if (length(error_msg) == 0) NA_character_ else error_msg,
    status = if (length(error_msg) == 0) "SUCCESS" else "FAILED")
  
  train_obj <- list(
    obj = mod,
    metadata = meta,
    partition = part,
    training = training)
  
  if (train_obj$training$status == "SUCCESS") {
    saveRDS(
      train_obj,
      file.path(mod_dir, paste0(comparison, "_model.rds")))
    
    message("   > model saved to: ", mod_dir)
  }
  
  return(train_obj)
}
  
# Filter for current validation rows from comparison_table, load validation
# object, externally validate trained object using RF and XGBoost, save as RDS 
# and return final evaluated object
validate_model <- function(
    train_obj,
    cfg,
    comparisons, 
    in_dir) {
  
  message("\nRunning external validation")
  
  # Define results path
  fin_dir <- file.path(out_dir, "final")
  dir.create(fin_dir, recursive = TRUE, showWarnings = FALSE)
  validation <- list()
  
  # Extract obj items
  result <- train_obj$obj
  meta <- train_obj$metadata
  part <- train_obj$partition
  train <- train_obj$train
  
  for (i in seq_len(nrow(comparisons))) {
    
    row <- comparisons[i, ]
    val_name <- row$val_name
    comparison_type <- row$comparison_type
    
    if (comparison_type == "internal") {
      message(" [", i, "] ", val_name, " (internal)\n")
      next
    }
    
    message(" [", i, "] ", val_name)
    
    # Load validation cohort
    val_rds <- readRDS(
      file.path(in_dir, paste0("final/", val_name, "_final.rds")))
    
    val_obj <- val_rds$obj
    
    # Validate model
    warnings <- character()
    error_msg <- character()
    
    message ("   > RF model validation")
    tmp <- suppressMessages(
      withCallingHandlers(
        tryCatch({
          ValidateModelOnData(
            result,
            ValidationData = val_obj,
            model_type = "RF",
            TaskName = paste0(val_name),
            TrueLabel = "disease",
            condition_col = "validation_condition",
            PlotAUC = FALSE)
          
        }, error = function(e) {
          error_msg <<- e$message
          message("RF validation failed: ", e$message)
          NULL
        }),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }))
    
    if (!is.null(tmp))
      result <- tmp
    
    # Extract results
    rf_res <- result@PredictResult$RF[[val_name]]
    
    AUC_RF <- if (is.null(rf_res)) NA_real_ else {
      extract_auc(rf_res$AUC)
    }

    message ("   > XGBoost model validation")
    tmp <- suppressMessages(
      withCallingHandlers(
        tryCatch({
          ValidateModelOnData(
            result,
            ValidationData = val_obj,
            model_type = "XGBoost",
            TaskName = paste0(val_name),
            TrueLabel = "disease",
            condition_col = "validation_condition",
            PlotAUC = FALSE)
          
        }, error = function(e) {
          error_msg <<- e$message
          message("XGB validation failed: ", e$message)
          NULL
        }),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }))
    
    if (!is.null(tmp))
      result <- tmp
    
    # Extract results
    xgb_res <- result@PredictResult$XGB[[val_name]]
    
    AUC_XGB <- if (is.null(xgb_res)) NA_real_ else {
      extract_auc(xgb_res$AUC)
    }
    
    message("")
    
    # Extract metrics
    validation[[i]] <- data.frame(
      train_name = meta$comparison,
      train_cohort = meta$cohort,
      train_disease = meta$disease,
      val_name = val_name,
      val_cohort = row$val_cohort,
      val_disease = row$val_disease,
      comparison_type = row$comparison_type,
      status = if (length(error_msg) == 0) "SUCCESS" else "FAILED",
      AUC_RF = AUC_RF,
      AUC_XGB = AUC_XGB,
      warnings = if (length(warnings) == 0) NA_character_ else
        paste(unique(warnings), collapse = " | "),
      error = if (length(error_msg) == 0) NA_character_ else error_msg,
      stringsAsFactors = FALSE)
    
  }
  final_obj <- list(
    obj = result,
    metadata = meta,
    partition = part,
    training = train,
    validation = bind_rows(validation))
  
  saveRDS(
      final_obj,
      file.path(fin_dir, paste0(meta$comparison, "_evaluated.rds")))

  message("Evaluated model saved as: ", fin_dir, "/", meta$comparison, "_evaluated.rds")

  return(final_obj)
}
  
# Export training and validation results from final object and save as CSV
export_logs <- function(final_obj, out_dir) {
  
  training_log <- file.path(out_dir, "training_log")
  dir.create(training_log, recursive = TRUE, showWarnings = FALSE)
  
  validation_log <- file.path(out_dir, "validation_log")
  dir.create(validation_log, recursive = TRUE, showWarnings = FALSE)
  
  train_name <- final_obj$metadata$comparison
  
  summary_df <- c(
    final_obj$metadata,
    final_obj$partition,
    final_obj$training)
  
  summary_df <- as.data.frame(summary_df, stringsAsFactors = FALSE)
  
  write.csv(
    summary_df,
    file.path(training_log, paste0(train_name, ".csv")),
    row.names = FALSE)
  
  write.csv(
    final_obj$validation,
    file.path(validation_log, paste0(train_name, ".csv")),
    row.names = FALSE)
}


# Define main functions --------------------------------------------------------
run_modelling <- function(
    comparison_table,
    file,
    in_dir,
    out_dir) {
  
  cfg <- get_configuration()
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load training object
  message("\nLoading file: ", basename(file))
  train_obj <- readRDS(file)
  train_name <- train_obj$metadata$comparison
  
  # Extract comparison rows for current train cohort
  comparisons <- comparison_table[comparison_table$train_name == train_name,]
  
  if (nrow(comparisons) == 0)
    stop("   > No comparisons found.")
  
  message("   > ", nrow(comparisons), " comparisons found")
  
  # Partition train_obj
  train_obj <- partition_training(train_obj)
  
  # Train and internally validate model
  train_obj <- train_model(
    train_obj = train_obj, 
    cfg = cfg)
  
  # Externally validate model
  final_obj <- validate_model(
    train_obj = train_obj,
    comparisons = comparisons,
    in_dir = in_dir,
    cfg = cfg)
  
  export_logs(
    final_obj = final_obj,
    out_dir = out_dir)
  
  return(final_obj)
}

# Load data and define args ----------------------------------------------------
comparison_table <- read.csv(
  file.path(in_dir, "comparison_table.csv"),
  stringsAsFactors = FALSE)

# For local testing
file <- "results/04_pre_ml/final/FengQ_2015_adenoma_final.rds"

# For Slurm
# args <- commandArgs(trailingOnly = TRUE)
# file <- args[1]
# 
# if (length(args) != 1) {
#   stop("Usage: Rscript 05_ml.R <final.rds>")
# }


# Execute ----------------------------------------------------------------------
modelling_res <- run_modelling(comparison_table, file, in_dir, out_dir)

