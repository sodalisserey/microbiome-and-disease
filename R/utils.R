# Utility functions for reading, writing and processing curatedMetagenomicData cohorts

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(SummarizedExperiment)


# Reading and writing ----------------------------------------------------------
#' Read cohort list from txt file 
read_cohort_list <- function(path) { 
  
  lines <- readLines(path, warn = FALSE)
  
  # Remove comments + trim whitespace
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#", lines)]
  
  unique(lines)
}

#' Load cohorts from list and generate list of cohort objects
load_cohorts <- function(cohorts) {

  out <- lapply(cohorts, function(study) {
    message("\nLoading: ", study)

    tryCatch({

      obj <- curatedMetagenomicData(
        paste0(study, ".relative_abundance"),
        dryrun = FALSE,
        rownames = "short")[[1]]

      meta <- SummarizedExperiment::colData(obj)
      
      if ("disease" %in% colnames(meta)) {
        meta$disease <- factor(clean_condition(as.character(meta$disease)))
        
        SummarizedExperiment::colData(obj) <- meta
      }
    
      obj

    }, error = function(e) {
      message("  ERROR loading ", study, ": ", e$message)
      NULL
    })
  })

  message("Loaded: ", length(cohorts), " cohorts")
  names(out) <- cohorts
  out
}

#' Write multiple data frames to csv files
write_csvs <- function(out_dir, data_list) {
  for (nm in names(data_list)) {
    write.csv(
      data_list[[nm]],
      file.path(out_dir, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }
  message("\nDone. CSVs written to: ", file.path(out_dir))
}


# Processing conditions --------------------------------------------------------
#' Remove / from string to generate clean label
clean_condition <- function(x) {
  gsub("[/-]", "", x)
}

#' Process conditions of a cohort by removing NA or missing values and extracting 
#' column from metadata
process_conditions <- function(
    cohort,
    condition_col = "study_condition",
    healthy_label = "control") {
                           
  # Extract condition column and remove NA samples
  meta <- colData(cohort)
  meta[[condition_col]] <- as.character(meta[[condition_col]])
  
  # Filter valid samples
  keep <- !is.na(meta[[condition_col]]) &
    meta[[condition_col]] != "" &
    meta[[condition_col]] != " "
  
  cohort <- cohort[, keep]
  
  # Extract conditions after filtering
  meta <- colData(cohort)
  
  conditions <- unique(meta[[condition_col]])
  
  diseases <- setdiff(conditions, healthy_label)
  
  list(
    cohort = cohort,
    healthy_present = healthy_label %in% conditions,
    conditions = conditions,
    diseases = diseases,
    n_diseases = length(diseases))
}


#' Check conditions contrast in cohort and print error message if no healthy or
#' disease samples present 
validate_conditions <- function(info, cohort_name) {
  
  if (!info$healthy_present || info$n_diseases == 0) {
    
    reason <- dplyr::case_when(
      !info$healthy_present ~ "No healthy samples",
      info$n_diseases == 0 ~ "Healthy samples only"
    )
    
    message(reason, " in ", cohort_name, ". Skipping.")

    return(reason)
  }
  
  return(NULL)
}


#' Get computing configuration
get_config <- function() {
  
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

#' Extract AUC from CrcBiomeScreenObject
extract_auc <- function(auc_val) {
  auc <- suppressWarnings(as.numeric(sub(".*: ", "", auc_val)))
  if (is.na(auc))
    auc <- NA_real_
  
  message("   AUC = ", auc)
  
  auc
}

#' Compile individual log csvs from array outputs
compile_logs <- function(
    log_dir,
    file_name) {
  
  files <- list.files(
    log_dir, 
    pattern = "\\.csv$", 
    full.names = TRUE)
  
  # Exclude output file
  files <- files[basename(files) != file_name]
  
  combined_df <- bind_rows(lapply(files, read_csv, show_col_types = FALSE))
  
  write_csv(combined_df, file.path(log_dir, file_name))
  
  message("\nCombined logs written to: ", file.path(log_dir))
  
  return(combined_df)
}


read_rds_files <- function(file_dir) {
  
  files <- list.files(
    path = file.path(file_dir),
    pattern = "\\.rds$",
    full.names = TRUE)
  
  names(files) <- file_path_sans_ext(basename(files))
  
  rds_files <- lapply(files, readRDS)
  
  message("\nRDS files read from: ", file_dir)
  return(rds_files)
}

