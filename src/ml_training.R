# Train models on primary cohort datasets
# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
source("src/handler.R")


# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary_cohorts.rds") |>
  create_cohort_objects()

test_cohorts <- readRDS("data/test_cohorts.rds") |>
  create_cohort_objects()


# Define helper functions ------------------------------------------------------
#' Prepare by splitting etc. for modelling
prep_for_training <- function(obj,
                              disease = NULL) {
  
  # Split taxa and normalise
  obj <- SplitTaxas(obj)
  taxa <- getTaxaData(obj)
  colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  setTaxaData(obj) <- taxa
  
  obj <- KeepTaxonomicLevel(obj, level = "Genus")
  
  obj <- NormalizeData(obj, method = "GMPR", level = "Genus")
  
  obj <- FilterDataSet(
    obj,
    label = c("control", disease),
    condition_col = "study_condition"
  )
  # QC
  obj <- qcByCmdscale(
    obj,
    TaskName = paste0(comparison_name, "_QC"),
    normalize_method = "GMPR",
    plot = FALSE
  )
  
  # Split samples into training and testing sets
  obj <- SplitDataSet(
    obj,
    label = c("control", disease),
    partition = 0.7
  )
  
  # TODO class balance plot overwrites old one each time
  # Class balance in training data
  balance <- checkClassBalance(getModelData(obj)$TrainLabel)
  print(table(getSampleData(obj)$study_condition))
  # TODO why dropping samples (when compared to contingency table)
  
  return(list(
    obj = obj,
    balance = balance
  ))
  }

train_model <- function(obj, disease, comparison_name, balance,
                                       model_type = "RF",
                                       n_cv = 2,
                                       num_cores = 1) {
  
  class_weights <- if (balance$suggestion == "balanced") {
    FALSE
    message("Class weights = FALSE")
  } else {
    TRUE
    message("Class weights = TRUE")
  }
  
  suppressMessages(
    rf_result <- tryCatch({
      obj_rf <- TrainModels(
        obj,
        model_type = "RF",
        TaskName = paste0(comparison_name, "_RF"),
        ClassWeights = FALSE,
        TrueLabel = disease,
        num_cores = 1,
        n_cv = 2
      )

      obj_rf <- EvaluateModel(
        obj_rf,
        model_type = "RF",
        TaskName = paste0(comparison_name, "_RF_test"),
        TrueLabel = disease,
        PlotAUC = FALSE
      )
    }, error = function(e) {
      message("RF failed: ", e$message)
      NULL
    })
  )
  return (obj_rf)
}
# Test datasets
# test_cohorts <- read_cohort_list("data/test.txt")
# saveRDS(load_cohorts(test_cohorts), "data/test_cohorts.rds")




# Define main function ---------------------------------------------------------
run_model_training <- function(cohort_list) {
  
  results <- list()
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  out_dir = "results/ml_training"
  setwd(out_dir)
  
  for (cohort_name in names(cohort_list)) {  
    
    # Process and extract study_conditions from cohort
    cohort <- cohort_list[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Filter out cohorts with no control (healthy) samples
    if (!info$healthy_present || info$n_diseases == 0) {
      
      reason <- dplyr::case_when(
        !info$healthy_present ~ "No healthy samples",
        info$n_diseases == 0 ~ "Healthy samples only"
      )
      
      message(reason, " in ", cohort_name, ". Skipping.")
      
      next
    
    # TRAINING BRANCH 1: healthy x 1 disease
    } else if (info$n_diseases == 1) {
      
      disease = info$diseases
    
      comparison_name <- paste(cohort_name, disease, sep = "_")
      message("\n", cohort_name, ": healthy vs ", info$diseases)
      
      obj <- CreateCrcBiomeScreenObject(
        RelativeAbundance = cohort@assays@data@listData$relative_abundance,
        TaxaData = cohort@rowLinks$nodeLab,
        SampleData = cohort@colData
      )
      
      prep <- prep_for_training(obj, disease)
      
      obj <- prep$obj
      balance <- prep$balance
      
      obj_rf <- train_model(
        obj = obj,
        disease = disease,
        comparison_name = comparison_name,
        balance = balance
      )
      
      results[[comparison_name]] <- obj_rf
    
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
      }
    }
  }
  return(results)
  
  setwd("..")
  setwd("..")
}


# Execute ----------------------------------------------------------------------
res <- run_model_training(test_cohorts)
  