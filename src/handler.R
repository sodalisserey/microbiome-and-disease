# Utility functions for reading and writing curatedMetagenomicData cohorts

# Load packages
library(curatedMetagenomicData)
library(SummarizedExperiment)


#' Load every dataset listed in datasets.R
#' 
#' @param cohorts  Named list as defined in datasets.R 
#' @return Named list of SummarizedExperiment objects (NULL on failure)

load_cohorts <- function(cohorts) {
  
  loaded <- lapply(cohorts, function(study) {
    message("\n[handler] Loading: ", study)
    
    tryCatch(
      curatedMetagenomicData(
        paste0(study, ".relative_abundance"),
        dryrun = FALSE
        )[[1]],
      error = function(e) {
        message("  ERROR loading ", study, ": ", e$message)
        return(NULL)
      }
    )
  })
  
  names(loaded) <- cohorts
  loaded
}


