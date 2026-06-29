# Process ML models for external validation using CrcBiomeScreen package
# 1. Combine individual logs from 05_ml_train/logs/[cohort_name]_model_log.csv

# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
source("R/utils.R")


# Define input/output directories ----------------------------------------------
in_dir <- "results/05_ml_train"
out_dir <- "results/06_ml_validate_prep"


# Define helper functions ------------------------------------------------------
#' 

run_alignment <- function () {
  
}




#'
log_alignment <- function(){}


#'
run_pipeline <- function (
  model_obj,
  cohort_name,
  disease,
  comp_obj,
  out_dir,
  processing_res) {
  
  prep_dir <- file.path(out_dir, "processed")
  dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)
  
  comparison <- paste(cohort_name, disease, sep = "_")
  
  message("Aligning: ", comparison, " with... ")
  
  # Filter for single disease
  suppressWarnings(
    comparison_obj <- FilterDataSet(
      comp_obj$obj,
      label = c("control", disease),
      condition_col = "study_condition"))
  
  
  
  
}
  
# get_within_same_pairs()
# get_within_different_pairs()
# get_cross_same_pairs()
# get_cross_different_pairs()

# Define main functions --------------------------------------------------------
#'

# make for different types: external val
run_validation <- function(
    in_dir,
    out_dir,
    models,
    cohort_list,
    val_mode # = within_diff, cross_same, cross_diff
    ){
  
  # Initialise storage object
  processing_res <- list()
  
  # Compile array output logs
  run_log_compilation(in_dir)
  
  
  # validation_mode =
  # "within_same"
  # "within_different"
  # "cross_same"
  # "cross_different"
}
  

# Load data --------------------------------------------------------------------
models <- readRDS(file.path(paste0(in_dir, "all_models.rds")))


# Execute ----------------------------------------------------------------------
validation_res <- run_evaluation(models, out_dir)
export_validation(validation_res, out_dir)


