# Utility functions for reading, writing and processing curatedMetagenomicData cohorts

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(SummarizedExperiment)


# Define functions -------------------------------------------------------------
#' Read cohort list from txt file 
read_cohort_list <- function(path) { 
  
  lines <- readLines(path, warn = FALSE)
  
  # Remove comments + trim whitespace
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#", lines)]
  
  unique(lines)
}

#' Load cohorts from list and generate list of object, abundance and metadata
load_cohorts <- function(cohorts) {
  
  out <- lapply(cohorts, function(study) {
    message("\nLoading: ", study)
    
    tryCatch({
      
      obj <- curatedMetagenomicData(
        paste0(study, ".relative_abundance"),
        dryrun = FALSE,
        rownames = "short")[[1]]
      
      abundance <- assay(obj)   
      metadata  <- colData(obj)
      
      list(
        object = obj,    
        abundance = abundance,
        metadata = metadata   
      )
      
    }, error = function(e) {
      message("  ERROR loading ", study, ": ", e$message)
      NULL
    })
  })
  
  message("Loaded: ", length(cohorts), " cohorts")
  names(out) <- cohorts
  out
}

#' Extract TreeSummarizedExperiment objects from primary cohorts dataset
create_cohort_obj <- function(cohorts) {
  
  lapply(names(cohorts), function(nm) {
    
    x <- cohorts[[nm]]
    
    if (is.list(x) && !is.null(x$object)) {
      x <- x$object
    }
    
    if (is.null(x)) {
      stop("NULL cohort in: ", nm)
    }
    x
  }) |> setNames(names(cohorts))
}

#' Process conditions of a cohort by removing NA or missing values and extracting 
#' column from metadata
process_conditions <- function(cohort,
                            condition_col = "study_condition",
                            healthy_label = "control") {
  # Extract condition column and remove NA samples
  meta <- as.data.frame(colData(cohort))
  meta[[condition_col]] <- as.character(meta[[condition_col]])
  
  keep <- !is.na(meta[[condition_col]]) &
    meta[[condition_col]] != "" &
    meta[[condition_col]] != " "
  
  cohort <- cohort[, keep]
  
  # Extract conditions
  meta <- meta[keep, , drop = FALSE]
  
  conditions <- meta[[condition_col]]
  conditions <- unique(conditions[!is.na(conditions)])
  
  diseases <- setdiff(conditions, healthy_label)
  
  list(
    cohort = cohort,
    healthy_present = healthy_label %in% conditions,
    diseases = diseases,
    n_diseases = length(diseases)
  )
}

#' Remove / from string to generate clean label
clean_label <- function(x) gsub("/", "", x)

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

#' Get computing configuration
get_config <- function() {
  
  is_slurm <- nzchar(Sys.getenv("SLURM_JOB_ID"))
  
  if (is_slurm) {
    cfg <- list(
      num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
      n_cv = 10,
      mode = "slurm")
    message("   Running under Slurm: num_cores = ", cfg$num_cores, ", n_cv = ", cfg$n_cv)
    
  } else {
    cfg <- list(
      num_cores = 1,
      n_cv = 2,
      mode = "local")
    message("   Running locally: num_cores = ", cfg$num_cores, ", n_cv = ", cfg$n_cv)
  }
  # message(paste0("   num_cores = ", cfg$num_cores))
  # message(paste0("   n_cv = ", cfg$n_cv))
  
  return(cfg)
}

#' Extract AUC from CrcBiomeScreenObject
extract_auc <- function(auc_val) {
  tryCatch(
    as.numeric(sub(".*: ", "", auc_val)),
    error = function(e) NA_real_)
}
