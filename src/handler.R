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
    message("\n[handler] Loading: ", study)
    
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
  
  names(out) <- cohorts
  out
}

#' Extract TreeSummarizedExperiment objects from primary cohorts dataset
create_cohort_objects <- function(cohorts) {
  
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


# Read and save cohorts --------------------------------------------------------
primary_cohorts <- read_cohort_list("data/primary.txt")

#saveRDS(load_cohorts(primary_cohorts), "data/primary_cohorts.rds")
