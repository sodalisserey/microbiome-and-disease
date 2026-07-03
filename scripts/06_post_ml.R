# Check ML results and compile logs using functions from utils.R
# 1. Define input/output directories 
# 2. Execute main function run_ml_check() which:
#     * checks for missing results from expected in_dir/comparisons.txt list
#     * compiles all individual training logs from out_dir/training_log and
#       saves as out_dir/combined_training.csv
#     * compiles all individual validation logs from out_dir/validation_log and
#       saves as out_dir/combined_validation.csv

# Load packages and dependencies -----------------------------------------------
library(tools)
source("R/utils.R")

# Define input/output directories ----------------------------------------------
in_dir <- "results/04_pre_ml"
out_dir <- "results/05_ml"


# Define main function ---------------------------------------------------------
# Check for missing results, compile all logs and save as CSV
run_ml_check <- function(
    in_dir,
    out_dir) {
  
  expected <- file_path_sans_ext(
    basename(readLines(file.path(in_dir, "comparisons.txt"))))
  
  expected <- sub("_final$", "", expected)
  
  completed <- sub("_evaluated\\.rds$", "", list.files(file.path(out_dir, "final")))
  
  missing <- setdiff(expected, completed)
  
  if (length(missing) == 0) {
    message("All comparisons completed successfully")
    
    train_dir <- file.path(out_dir, "training_log")
    compile_logs(
      log_dir = train_dir, 
      file_name = "combined_training.csv")
    
    val_dir <- file.path(out_dir, "validation_log")
    compile_logs(
      val_dir, 
      file_name = "combined_validation.csv")
    
  } else {
    message(length(missing), " comparison(s) missing:")
    print(missing)
  }
}


# Execute ----------------------------------------------------------------------
run_ml_check(in_dir, out_dir)


