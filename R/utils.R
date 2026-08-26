# Utility functions used for the analysis of `curatedMetagenomicData` cohorts


# Functions for reading and writing --------------------------------------------
# Read cohort list from txt file 
read_cohort_list <- function(path) { 
  
  lines <- readLines(path, warn = FALSE)
  
  # Remove comments + trim whitespace
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#", lines)]
  
  unique(lines)
}


# Load cohorts from list and generate list of cohort objects
load_cohorts <- function(cohorts) {

  out <- lapply(cohorts, function(study) {
    message("\nLoading: ", study)

    tryCatch({

      obj <- curatedMetagenomicData(
        paste0(study, ".relative_abundance"),
        dryrun = FALSE,
        rownames = "long")[[1]]

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


# Read a list of RDS files
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


# Write multiple data frames to csv files
write_csvs <- function(data_list, out_dir) {
  for (nm in names(data_list)) {
    write.csv(
      data_list[[nm]],
      file.path(out_dir, paste0(nm, ".csv")),
      row.names = FALSE)
  }
  message("\nDone. CSVs written to: ", file.path(out_dir))
}


# Functions for processing -----------------------------------------------------
# Remove / from string to generate clean label
clean_condition <- function(x) {
  gsub("[/-]", "", x)
}


# Process conditions by removing NA/missing values and extracting metadata
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
  
  # Extract study conditions after filtering
  meta <- colData(cohort)
  
  study_conditions <- unique(meta[[condition_col]])
  
  conditions <- setdiff(study_conditions, healthy_label)
  
  list(
    cohort = cohort,
    case_present = healthy_label %in% study_conditions,
    study_conditions = study_conditions,
    conditions = conditions,
    n_conditions = length(conditions))
}


# Check for presence of both case and controls
validate_conditions <- function(info, cohort_name) {
  
  if (!info$case_present || info$n_conditions == 0) {
    
    reason <- dplyr::case_when(
      !info$case_present ~ "No cases",
      info$n_conditions == 0 ~ "Cases only")
    
    message(reason, " in ", cohort_name, ". Skipping.")
    return(reason)
  }
  
  return(NULL)
}

