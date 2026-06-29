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



# ML training and validation ---------------------------------------------------
#'
normalise_obj <- function(
    cohort,
    cohort_name,
    norm_method = "GMPR") {
  
  warnings <- character()
  
  obj <- withCallingHandlers({
    
    obj <- CreateCrcBiomeScreenObject(
      RelativeAbundance = cohort@assays@data@listData$relative_abundance,
      TaxaData = cohort@rowLinks$nodeLab,
      SampleData = cohort@colData)
    
    # Split taxa and keep genus level
    obj <- SplitTaxas(obj)
    taxa <- getTaxaData(obj)
    colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    setTaxaData(obj) <- taxa
    obj <- KeepTaxonomicLevel(obj, level = "Genus")
    
    # Normalise data
    obj <- NormalizeData(obj, method = norm_method, level = "Genus")
    
    # Standardise labels
    obj@SampleData$study_condition <-
      clean_condition(obj@SampleData$study_condition)
    
    obj
    
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  
  return(list(
    obj = obj,
    warnings = warnings,
    norm_method = norm_method))
}

qc_obj <- function(
    cohort,
    cohort_name,
    obj,
    norm_method,
    warnings,
    out_dir) {
  
  qc_dir <- file.path(out_dir, "plots")
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  
  suppressMessages(
    obj <- qcByCmdscale(
      obj,
      TaskName = paste0(cohort_name, "_QC"),
      outdir = qc_dir,
      normalize_method = norm_method,
      plot = TRUE))
  
  file.rename(
    file.path(qc_dir, paste0("cmdscale_", cohort_name, "_QC_", norm_method, ".pdf")),
    file.path(qc_dir, paste0(cohort_name, "_QC_", norm_method, ".pdf")))
  
  return(list(
    outliers = obj@OutlierSamples,
    obj = obj,
    cohort_name = cohort_name,
    norm_method = norm_method,
    warnings = warnings,
    
    n_features = ncol(obj@OriginalNormalizedData),
    n_samples = nrow(obj@OriginalNormalizedData),
    
    n_features_qc = ncol(obj@NormalizedData),
    n_samples_qc = nrow(obj@NormalizedData)))
}

#'
log_qc_obj <- function(
    cohort_name, 
    processed_obj,
    processing_res) {
  
  # initialise containers if missing
  if (is.null(processing_res$qc_summary)) {
    processing_res$qc_summary <- list()
  }
  
  row <- data.frame(
    cohort_name = cohort_name,
    norm_method = processed_obj$norm_method,
    warnings = if (length(processed_obj$warnings) == 0) NA_character_ else
      paste(processed_obj$warnings, collapse = " | "),
    
    # N features (genus level) and samples pre-normalisation/qc
    n_features = processed_obj$n_features,
    n_samples = processed_obj$n_samples,
    
    # Number of outliers identified by QC
    outliers = if (length(processed_obj$outliers) == 0) NA_character_ else 
      paste(processed_obj$outliers, collapse = "; "),
    
    # N features (genus level) and samples post-normalisation/qc
    n_features_qc = processed_obj$n_features_qc,
    n_samples_qc = processed_obj$n_samples_qc,
    
    stringsAsFactors = FALSE)
  
  processing_res$qc_summary[[cohort_name]] <- row
  
  return(processing_res)
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
  return(cfg)
}

#' Extract AUC from CrcBiomeScreenObject
extract_auc <- function(auc_val) {
  tryCatch(
    as.numeric(sub(".*: ", "", auc_val)),
    error = function(e) NA_real_)
}

#' Compile individual log csvs from array outputs
compile_logs <- function(
    in_dir) {
  
  files <- list.files(
    file.path(in_dir, "logs"), 
    pattern = "\\model_log.csv$", 
    full.names = TRUE)
  
  head(read.csv(files[1]))
  
  combined_df <- do.call(
    rbind,
    lapply(files, read.csv, stringsAsFactors = FALSE))
  
  out_file <- file.path(in_dir, "combined_logs.csv")
  write.csv(combined_df, out_file, row.names = FALSE)
  
  message("\nCombined logs written to: ", file.path(in_dir))
}

