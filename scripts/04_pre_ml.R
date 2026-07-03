# Prepare raw cohorts for ML training using CrcBiomeScreen package
# 1. Define output directory and functions
# 2. Load primary.RDS 
# 3. Execute main function run_processing() which runs:
#     a. normalise_cohorts()
#         * normalises list of raw cohorts
#         * saves as out_dir/raw_normalised/[cohort]_raw_normalised.rds
#     b. harmonise_cohorts()
#         * combines intersecting features of all normalised cohorts
#         * sets empty features to 0
#     c. qc_cohorts()
#         * runs QC on normalised and harmonised cohorts
#         * saves as out_dir/norm_processed/[cohort]_norm_processed.rds
#     d. get_comparisons()
#         * loops through cohorts for unique healthy x disease combinations
#         * subsets QC processed cohorts by disease
#         * checks for sufficient sample size (>= 30) and add binary labels
#         * saves as out_dir/final/[comparison]_final.rds
#         * lists file paths of _final.rds objects as out_dir/comparisons.txt
#     e. create_comparison_table()
#         * creates evaluation matrix of unique comparisons
#         * extracts train/val metadata and create comparison_type tags
#         * saves as out_dir/comparison_table.csv


# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
library(dplyr)
source("R/utils.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/04_pre_ml"


# Define helper functions ------------------------------------------------------
# Normalise list of raw cohorts, save as RDS and return objects
normalise_cohorts <- function(
    cohort_list,
    out_dir,
    norm_method = "GMPR") {
  
  message("Normalising raw cohorts:")
  
  # Initialise
  norm_dir <- file.path(out_dir, "raw_normalised")
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  norm_list <- list()
  
  for (cohort_name in names(cohort_list)) {
    
    message("   ", cohort_name)

    cohort <- cohort_list[[cohort_name]]
    
    # Process study conditions
    info <- process_conditions(cohort)
    
    # Skip invalid cohorts
    if (!is.null(validate_conditions(info, cohort_name))) next
    
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
    
    # Write normalised objects as RDS
    saveRDS(
      obj,
      file.path(norm_dir, paste0(cohort_name, "_raw_normalised.rds")))

    # Store in memory
    norm_list[[cohort_name]] <- obj
    
    # TODO Log warnings/errors and norm_method
    
  }
  message(length(norm_list), " cohorts saved to: ", norm_dir)
  
  return(norm_list)
}

# Harmonise features of normalised cohorts and return objects
harmonise_cohorts <- function(
    norm_list,
    global_features) {
  
  # Initialise
  harmonised_list <- list()
  
  for (cohort_name in names(norm_list)) {
    
    obj <- norm_list[[cohort_name]]
    
    # Extract normalised data
    norm_data <- getNormalizedData(obj)
    
    # Find missing genera
    missing_features <- setdiff(global_features, colnames(norm_data))
    
    # Add zero-filled columns
    if (length(missing_features) > 0) {
      
      zero_mat <- matrix(
        0,
        nrow = nrow(norm_data),
        ncol = length(missing_features),
        dimnames = list(NULL, missing_features)
      )
      
      norm_data <- cbind(norm_data, zero_mat)
    }
    
    # Reorder to global feature order
    norm_data <- norm_data[, global_features, drop = FALSE]
    
    # Store back in object
    setNormalizedData(obj) <- norm_data
    
    # Store in memory
    harmonised_list[[cohort_name]] <- obj
  }
  
  message("   All cohorts harmonised")
  return(harmonised_list)
}

# QC list of normalised and harmonised cohorts, save as RDS and return objects
qc_cohorts <- function(
    harmonised_list,
    out_dir,
    norm_method = "GMPR",
    plot = FALSE) {
  
  message("\nRunning QC on normalised cohorts:")
  
  # Initialise
  qc_dir <- file.path(out_dir, "norm_processed")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE) 
  
  if (plot == TRUE) {
    plot_dir <- file.path(qc_dir, "qc_plots") 
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE) 
  }
  qc_list <- list()
  qc_log <- list()
  
  for (cohort_name in names(harmonised_list)) {
    
    message("   ", cohort_name)
    
    harmonised_cohort <- harmonised_list[[cohort_name]]
    
    warnings <- character()
    
    if (plot == TRUE) { 
      
      suppressMessages( 
        obj <- qcByCmdscale( 
          harmonised_cohort, 
          TaskName = paste0(cohort_name, "_QC"), 
          outdir = plot_dir, 
          normalize_method = norm_method, 
          plot = TRUE)) 
      
      file.rename( file.path(
        plot_dir, 
        paste0("cmdscale_", cohort_name, "_QC_", norm_method, ".pdf")), 
        file.path(plot_dir, paste0(cohort_name, "_QC_", norm_method, ".pdf"))) 
      
    } else { 
        suppressMessages( 
          obj <- qcByCmdscale( 
            harmonised_cohort, 
            TaskName = paste0(cohort_name, "_QC"), 
            outdir = plot_dir, 
            normalize_method = norm_method, 
            plot = FALSE)) }
    
    meta <- list(
      cohort_name = cohort_name,
      norm_method = norm_method,
      warnings = if (length(warnings) == 0) NA_character_
      else paste(unique(warnings), collapse = " | "),
      n_features = ncol(obj@OriginalNormalizedData),
      n_samples = nrow(obj@OriginalNormalizedData),
      n_features_qc = ncol(obj@NormalizedData),
      n_samples_qc = nrow(obj@NormalizedData),
      n_outliers = length(obj@OutlierSamples),
      outliers = if (length(obj@OutlierSamples) == 0) NA_character_
      else paste(obj@OutlierSamples, collapse = "; "))
    
    # Store in memory
    qc_list[[cohort_name]]$obj <- obj
    qc_list[[cohort_name]]$metadata <- as.data.frame(meta)
    
    saveRDS(
      list(obj = obj, metadata = meta),
      file.path(
        qc_dir,
        paste0(cohort_name, "_norm_processed.rds")))
  }
  message(length(qc_list), " cohorts saved to: ", qc_dir)
  
  return(qc_list)
}

# Loop through QC-processed cohorts,  generate unique healthy x disease
# combinations, check for sufficient sample size (>=30), add binary validation 
# labels, save as RDS and return objects
get_comparisons <- function(
    qc_list,
    condition_col = "study_condition",
    healthy_label = "control",
    min_n = 30) {
  
  message("\nRetrieving comparisons:")
  
  # Initialise
  comp_dir <- file.path(out_dir, "final")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE) 
  comparison_list <- list()
  
  for (cohort_name in names(qc_list)) {
  
    cohort_obj <- qc_list[[cohort_name]]$obj
    
    # Extract conditions and get unique diseases
    meta <- cohort_obj@SampleData
    conditions <- unique(meta[[condition_col]])
    diseases <- setdiff(conditions, healthy_label)
    
    for (disease in diseases) {
      
      comparison <- paste0(cohort_name, "_", disease)
      
      # Filter sample size
      sub_meta <- meta[meta[[condition_col]] %in% c(healthy_label, disease), ]
      vals <- sub_meta[[condition_col]]
      
      n_control <- sum(vals == healthy_label, na.rm = TRUE)
      n_disease <- sum(vals == disease, na.rm = TRUE)
      
      if (n_control < min_n || n_disease < min_n) {
        message("   Skipping: ", comparison, " (n_control=", n_control,
                ", n_disease=", n_disease, ")")
        next
      }
      
      message("   Filtering: ", comparison)
      
      obj <- cohort_obj
      
      filtered_obj <- FilterDataSet(
        obj,
        label = c(healthy_label, disease),
        condition_col = condition_col)
      
      # Add binary validation label
      filtered_obj@SampleData$validation_condition <-
        ifelse(
          filtered_obj@SampleData[[condition_col]] == disease,
          "disease",
          "control")
      
      # Store result
      comparison_list[[comparison]] <- list(
        obj = filtered_obj,
        metadata = list(
          comparison = comparison,
          cohort = cohort_name,
          disease = disease,
          n_control = n_control,
          n_disease = n_disease))
      
      saveRDS(
        comparison_list[[comparison]],
        file.path(
          comp_dir,
          paste0(comparison, "_final.rds")))
    }
  }
  message(length(comparison_list), " comparisons saved to: ", comp_dir)
  
  # Index file paths as TXT
  files <- list.files(
    path = comp_dir,
    pattern = "_final\\.rds$",
    full.names = TRUE)

  writeLines(
    files, file.path(out_dir,"comparisons.txt"))

  return(comparison_list)
}

# Create evaluation matrix using unique comparisons, extract train/val cohort 
# metadata, add comparison_type tags and save as CSV
create_comparison_table <- function(
    comparison_list,
    out_dir) {
  
  # Build one row per comparison
  comparisons <- do.call(
    rbind,
    lapply(names(comparison_list), function(name) {
      
      obj <- comparison_list[[name]]$obj
      meta <- obj@SampleData
      
      disease <- unique(
        meta$study_condition[meta$validation_condition == "disease"])
      
      cohort <- sub(paste0("_", disease, "$"), "", name)
      
      data.frame(
        comparison = name,
        cohort = cohort,
        disease = disease,
        stringsAsFactors = FALSE)
    }))
  
  # Training cohort
  train_cohorts <- comparisons %>%
    mutate(train_name = paste0(cohort, "_", disease)) %>%
    transmute(
      train_name,
      train_cohort = cohort,
      train_disease = disease)
  
  # Validation datasets
  val_cohorts <- comparisons %>%
    mutate(val_name = paste0(cohort, "_", disease)) %>%
    transmute(
      val_name,
      val_cohort = cohort,
      val_disease = disease)
  
  # Cartesian product
  comparison_table <- merge(train_cohorts, val_cohorts, by = NULL)
  
  comparison_table <- comparison_table %>%
    arrange(train_name, val_name)
  
  # Label comparison type
  comparison_table <- comparison_table |>
    mutate(
      comparison_type = case_when(
        train_cohort == val_cohort &
          train_disease != val_disease ~ "same_cohort_diff_disease",
        
        train_cohort != val_cohort &
          train_disease == val_disease ~ "cross_cohort_same_disease",
        
        train_cohort != val_cohort &
          train_disease != val_disease ~ "cross_cohort_diff_disease",
        
        TRUE ~ "internal"))
  
  write.csv(
    comparison_table,
    file.path(out_dir, "comparison_table.csv"),
    row.names = FALSE)
  
  message("\nComparison table saved to ", out_dir)
  
  return(comparison_table)
}


# Define main functions --------------------------------------------------------
#' Run main ML normalisedaration for cohorts
run_processing <- function(cohort_list, out_dir) {

  # Normalise cohort_list
  norm_list <- normalise_cohorts(
    cohort_list = cohort_list, 
    out_dir = out_dir)
  
  # Extract global features
  global_features <- sort(unique(unlist(
    lapply(norm_list, function(obj)
      colnames(getNormalizedData(obj))))))
  
  message("\nGlobal features detected: ", length(global_features))
  
  # Harmonise norm_list
  harmonised_list <- harmonise_cohorts(
    norm_list = norm_list,
    global_features = global_features)

  # QC harmonised_list
  qc_list <- qc_cohorts(
    harmonised_list = harmonised_list, 
    out_dir = out_dir, 
    plot = TRUE)
  
  # Create unique healthy x disease comparisons
  processing_res <- get_comparisons(qc_list)
  
  # Generate comparison table
  comparison_table <- create_comparison_table(processing_res, out_dir)
  
  return(processing_res)
}


# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary.rds")


# Execute ----------------------------------------------------------------------
processing_res <- run_processing(primary_cohorts, out_dir)

