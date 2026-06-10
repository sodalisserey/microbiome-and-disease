# Train models on primary cohort datasets
# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
source("src/handler.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/ml_training"


# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary_cohorts.rds") |>
  create_cohort_objects()

test_cohorts <- readRDS("data/test_cohorts.rds") |>
  create_cohort_objects()

no_contrast_cohorts <- readRDS("data/no_contrast_cohorts.rds") |>
  create_cohort_objects()


# Define helper functions ------------------------------------------------------
#' Prepare cohort by creating CrcBiomeScreenObject, splitting/setting taxa, 
#' normalising data, filtering study_condition, running QC, splitting dataset 
#' into training/testing and checking class balance
prep_for_training <- function(cohort,
                              comparison_name,
                              disease,
                              out_dir,
                              partition_ratio = 0.7) {
  
  warnings <- character()
  
  result <- withCallingHandlers(
    {
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
      # TODO output/save normalisation warning
      obj <- NormalizeData(obj, method = "GMPR", level = "Genus")
      
      # Filter for control and disease
      obj <- FilterDataSet(
        obj,
        label = c("control", disease),
        condition_col = "study_condition")
      
      # Perform QC
      # TODO output result?
      obj_qc <- qcByCmdscale(
        obj,
        TaskName = paste0(comparison_name, "_QC"),
        outdir = out_dir,
        normalize_method = "GMPR",
        plot = TRUE)
      
      # Split cohort into training and testing sets
      obj <- SplitDataSet(
        obj_qc,
        label = c("control", disease),
        partition = partition_ratio)
      
      # TODO class balance plot overwrites old one each time
      # TODO why dropping samples (when compared to contingency table)
      
      # Check class balance in training subset
      balance <- checkClassBalance(getModelData(obj)$TrainLabel, 
                                   outdir = out_dir, 
                                   plot = FALSE)
      
      class_ratio <- max(balance$class_counts) / min(balance$class_counts)
      class_imbalance <- balance$is_imbalanced
      class_counts <- as.data.frame(table(getSampleData(obj)$study_condition))
      colnames(class_counts) <- c("condition", "n")
      
      list(
        obj = obj,
        outliers = obj_qc@OutlierSamples,
        class_ratio = class_ratio,
        class_imbalance = class_imbalance,
        class_counts = class_counts)
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      # invokeRestart("muffleWarning"))
    }
  )
  return(list(
    obj = result$obj,
    outliers = result$outliers,
    partition_ratio = partition_ratio,
    class_ratio = result$class_ratio,
    class_imbalance = result$class_imbalance,
    class_counts = class_counts,
    warnings = unique(warnings)
  ))
}


#' Train an ML model based on the training subset of a cohort
train_model <- function(
    obj,
    disease,
    comparison,
    imbalance,
    out_dir) {
  
  safe_comp <- gsub("[/\\\\]", "", comparison)
  
  model_dir <- file.path(out_dir, "models")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Class weights set to true or false according to class imbalance
  class_weights <- if (imbalance == FALSE) {
    message("Class weights = FALSE")
    FALSE
  } else {
    message("Class weights = TRUE")
    TRUE
  }
  
  error_msg <- NULL
  
  suppressMessages(
     obj_rf <- tryCatch({
       obj_rf <- TrainModels(
         obj,
        model_type = "RF",
        TaskName = paste0(safe_comp, "_RF"),
        ClassWeights = class_weights,
        TrueLabel = disease,
        num_cores = 1,
        n_cv = 2
      )
      
      obj_rf <- EvaluateModel(
        obj_rf,
        model_type = "RF",
        outdir = model_dir,
        TaskName = paste0(safe_comp, "_RF_test"),
        TrueLabel = disease,
        PlotAUC = FALSE
      )
      
      file.rename(
        file.path(model_dir, paste0("CrcBiomeScreenObject_", safe_comp, "_RF_test.rds")),
        file.path(model_dir, paste0(safe_comp, "_RF.rds"))
      )
      
      obj_rf

    }, error = function(e) {
      error_msg <<- e$message
      message("RF failed: ", e$message)
      NULL
  })

  # TODO XGBoost model
  )
  
  return(list(
    obj_rf = obj_rf,
    class_weights = class_weights,
    error = error_msg,
    status = if (is.null(error_msg)) "SUCCESS" else "FAILED"))
}

#'
log_qc <- function()

#'
log_training <- function()

# Define main function ---------------------------------------------------------
#'
run_training <- function(cohort_list, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialise storage objects
  results <- list(
    qc_summary = list(),
    training_log = list(),
    models = list()
  )
  
  models <- list()
  
  for (cohort_name in names(cohort_list)) {  
    
    # Process and extract study_conditions from cohort
    cohort <- cohort_list[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check cohort contains both healthy and disease samples
    condition_contrast <- validate_conditions(info, cohort_name)
    
    if (!is.null(condition_contrast)) {
      
      results$training_log[[length(results$training_log) + 1]] <- data.frame(
        comparison = cohort_name,
        status = "SKIPPED",
        reason = condition_contrast
      )
      
      next  
    }
      
    # TRAINING BRANCH 1: healthy x 1 disease
    if (info$n_diseases == 1) {
      
      disease = info$diseases
    
      comparison_name <- paste(cohort_name, disease, sep = "_")
      message("\n", cohort_name, ": healthy vs ", info$diseases)
      
      prep <- prep_for_training(
        cohort = cohort,
        comparison_name = comparison_name,
        disease = disease,
        out_dir = out_dir)
      
      class_counts <- paste(
        paste(
          prep$class_counts$condition,
          prep$class_counts$n,
          sep = "="),
        collapse = "; ")
      
      # TODO export samples x taxa number too
      # TODO make qc_log and training_log functions
      results$qc_summary[[length(results$qc_summary) + 1]] <- data.frame(
        comparison = comparison_name,
        warnings = if (length(prep$warnings) == 0) NA_character_
        else paste(prep$warnings, collapse = " | "),
        class_ratio = prep$class_ratio,
        class_imbalance = prep$class_imbalance,
        class_counts = class_counts,
        outliers = if (length(prep$outliers) == 0) NA_character_ 
        else paste(prep$outliers, collapse = "; "),
        stringsAsFactors = FALSE
      )
      
      mod <- train_model(
        obj = prep$obj,
        disease = disease,
        comparison = comparison_name,
        imbalance = prep$class_imbalance,
        out_dir = out_dir
      )
      
      results$training_log[[length(results$training_log) + 1]] <- data.frame(
        comparison = comparison_name,
        cohort = cohort_name,
        disease = disease,
        status = mod$status,
        error = if (is.null(mod$error)) NA_character_ else mod$error,
        class_weights = mod$class_weights,
        split_ratio = prep$partition_ratio,
        stringsAsFactors = FALSE
      )
      
      # qc_data <- extract_qc(obj)
      # results$qc_summary[[comparison_name]] <- qc_data

      results$models[[comparison_name]] <- mod$obj_rf
      
    # TRAINING BRANCH 2: healthy x >1 disease
    } else {
      message("\n", cohort_name, ": healthy vs ", info$n_diseases, " diseases")
      
      next

      for (disease in info$diseases) {
        message("\n", cohort_name, " subset: healthy vs ", disease)

        # Subset cohort for healthy x 1 disease
        meta <- as.data.frame(colData(cohort))
        keep <- meta$study_condition %in% c("control", disease)
        cohort_subset <- cohort[, keep]

        comparison_name <- paste(cohort_name, disease, sep = "_")
        
      }
    }
  }
  return(results)
}


#' Export training outputs as raw R objects and CSVs
export_models <- function(training, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Bulk save all models
  saveRDS(training$models, file.path(out_dir, "models.rds"))
  
  # Bind data frames
  qc_summary <- dplyr::bind_rows(training$qc_summary)
  training_log <- dplyr::bind_rows(training$training_log)
  
  # Write CSVs
  write_csvs(
    out_dir,
    list(
      qc_summary = qc_summary,
      training_log = training_log
    )
  )
}


# Execute ----------------------------------------------------------------------
training <- run_training(test_cohorts, out_dir)
export_models(training, out_dir)