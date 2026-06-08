# Utility functions for reading and writing curatedMetagenomicData cohorts

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(SummarizedExperiment)


# Define functions -------------------------------------------------------------

#'
read_cohort_list <- function(path) {
  
  lines <- readLines(path, warn = FALSE)
  
  # Remove comments + trim whitespace
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#", lines)]
  
  unique(lines)
}

#'
load_cohorts <- function(cohorts) {
  
  cohorts <- read_cohort_list("data/primary.txt")
  
  out <- lapply(cohorts, function(study) {
    message("\n[handler] Loading: ", study)
    
    tryCatch({
      
      obj <- curatedMetagenomicData(
        paste0(study, ".relative_abundance"),
        dryrun = FALSE
      )[[1]]
      
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

# Read and save cohorts --------------------------------------------------------
saveRDS(load_cohorts(primary_cohort_names), "data/primary_cohorts.rds")
