# Externally validate performance of ML models using CrcBiomeScreen package


# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
source("R/utils.R")


# Define input/output directories ----------------------------------------------
in_dir <- "results/ml_training"
out_dir <- "results/ml_validation"


# Define helper functions ------------------------------------------------------

# TODO do these need to be separate functions? or move to utils.R
#'
log_qc <- function(){}

#'
log_validation <- function(){}

#'
log_result <- function(){}

#'
evaluate_model <- function(
    train_cohort,
    train_disease,
    val_cohort,
    val_disease,
    external_AUC){
  
  validation_res <- log_validation(
    auc = ...,
    sensitivity = ...,
    specificity = ...,
    n_samples = ...
  )
  
  return(validation_res)
}
  
# get_within_same_pairs()
# get_within_different_pairs()
# get_cross_same_pairs()
# get_cross_different_pairs()

# Define main function ---------------------------------------------------------
# make for different types: external val
run_external_validation <- function(
    models,
    cohort_list,
    val_mode){
  # validation_mode =
  # "within_same"
  # "within_different"
  # "cross_same"
  # "cross_different"
}
  

# Load data --------------------------------------------------------------------
models <- readRDS(file.path(paste0(in_dir, "all_models.rds")))

# Execute ----------------------------------------------------------------------
validation_res <- run_evaluation(models, out_dir, seed)
export_validation(validation_res, out_dir)


