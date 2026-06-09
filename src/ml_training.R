# Train models on primary cohort datasets
# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
source("src/handler.R")


# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary_cohorts.rds") |>
  create_cohort_objects()

# test_cohorts <- read_cohort_list("data/test.txt") |>
#   load_cohorts() |>
#   saveRDS("data/test_cohorts.rds")

test_cohorts <- readRDS("data/test_cohorts.rds") |>
  create_cohort_objects()

# Define helper functions ------------------------------------------------------
prepare_crc <- function(obj,
                        label_col = "study_condition",
                        positive_label = "CRC") {
  
  obj <- SplitTaxas(obj)
  obj <- KeepGenusLevel(obj)
  obj <- NormalizeData(obj, method = "TSS")
  
  # standardise labels
  obj$SampleData[[label_col]] <- tolower(obj$SampleData[[label_col]])
  
  obj
}

# Test datasets
# test_cohorts <- read_cohort_list("data/test.txt")
# saveRDS(load_cohorts(test_cohorts), "data/test_cohorts.rds")




# Define main function ---------------------------------------------------------
run_pipeline <- function(cohort_list) {
  
  results <- list()
  
  for (cohort_name in names(cohort_list)) {  
    
    message("\nProcessing: ", cohort_name)
    
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
    
      # Build CrcBiomeScreen object
      obj <- CreateCrcBiomeScreenObject(
        RelativeAbundance = cohort@assays@data@listData$relative_abundance,
        TaxaData = cohort@rowLinks$nodeLab,
        SampleData = cohort@colData
      )

      # Split taxa and normalise
      obj <- SplitTaxas(obj)
      obj <- KeepTaxonomicLevel(obj, level = "Genus")
      obj <- NormalizeData(obj, method = "TSS", level = "Genus")
      
      # Split samples into training and testing sets
      obj <- SplitDataSet(
        obj,
        label = c("control", disease),
        partition = 0.7
      )
      
      # Run QC via cmdscale
      obj <- qcByCmdscale(
        obj,
        TaskName = paste0(comparison_name, "_QC"),
        normalize_method = "TSS",
        plot = FALSE
        )
      
      # TODO class balance plot overwrites old one each time
      checkClassBalance(getModelData(obj)$TrainLabel)
      
      print(table(getSampleData(obj)$study_condition))
      
      # TODO class weights switch on when recommended
      # Train model
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
      
      results[[comparison_name]] <- obj_rf
      
    }
  }
  return(results)
}
# Run --------------------------------------------------------------------------
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
out_dir = "results/ml_training"
setwd(out_dir)

res <- run_pipeline(test_cohorts)

setwd("..")
setwd("..")

ZhuF_2020 <- readRDS("results/ml_training/CrcBiomeScreenObject_ZhuF_2020_schizophrenia_RF_test.rds")
