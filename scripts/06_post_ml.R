# Check ML results and compile logs using functions from utils.R
# 1. Define input/output directories and functions
# 2. Execute main function run_ml_check() which runs:
#     a. check_missing_results()
#         * checks for missing results from in_dir/comparisons.txt list
#     b. compile_logs()
#         * compiles all individual training logs from out_dir/training_log and   
#           saves as out_dir/combined_training.csv
#         * compiles all individual validation logs from out_dir/validation_log  
#           and saves as out_dir/combined_validation.csv

# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
library(tools)
source("R/utils.R")

# Define input/output directories ----------------------------------------------
in_dir <- "results/04_pre_ml"
out_dir <- "results/05_ml"


# Define helper function -------------------------------------------------------
check_missing_results <- function(in_dir, out_dir) {
  
  expected <- file_path_sans_ext(
    basename(readLines(file.path(in_dir, "comparisons.txt"))))
  
  expected <- sub("_final$", "", expected)
  
  completed <- sub("_evaluated\\.rds$", "", list.files(file.path(out_dir, "final")))
  
  missing <- setdiff(expected, completed)
  
  if (length(missing) == 0) {
    message("All comparisons evaluated successfully")
  } else {
    message(length(missing), " comparison(s) missing:")
    message("   ", missing)
    message("")
  }
  
  return(missing)
}
  
  
# Compile individual logs, save and return combined csv
compile_logs <- function(in_dir, out_dir, out_file) {
  
  files <- list.files(
    in_dir, 
    pattern = "\\.csv$", 
    full.names = TRUE)
  
  # Exclude output file
  files <- files[basename(files) != out_file]
  
  combined_df <- bind_rows(lapply(files, read_csv, show_col_types = FALSE))
  
  write_csv(combined_df, file.path(out_dir, out_file))
  
  message("Combined logs written as: ", file.path(out_dir, out_file))
}

# Define main function ---------------------------------------------------------
# Check for missing results, compile all logs and save as CSV
run_ml_check <- function(
    in_dir,
    out_dir) {
  
  # Confirm all comparisons have been evaluated
  missing <- check_missing_results(in_dir, out_dir)
  
  if (length(missing) == 0) {

    train_dir <- file.path(out_dir, "training_log")
    
    compile_logs(
      in_dir = train_dir, 
      out_dir = out_dir,
      out_file = "combined_training.csv")
    
    val_dir <- file.path(out_dir, "validation_log")
    
    compile_logs(
      in_dir = val_dir, 
      out_dir = out_dir,
      out_file = "combined_validation.csv")
  }
}


# Execute ----------------------------------------------------------------------
run_ml_check(in_dir, out_dir)



result <- readRDS("results/05_ml/final/FengQ_2015_adenoma_evaluated.rds")


