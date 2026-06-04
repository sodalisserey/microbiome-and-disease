# Pipeline:
# 1. Read dataset registry from datasets.R and load cohorts via handler.R
# 2. For each cohort, extract unique study conditions (disease):
#    * If healthy + 1 disease  -> 
#       - Perform QC:
#           * Print contingency table (age_category x study_condition)
#           * Check class balance (n per group, imbalance ratio) 
#       - Run Lefser and print class/age balance warning message
#    * If healthy + >1 disease -> 
#       - Create pair-wise subsets (healthy + 1 disease)
#       - For each subset, perform QC: 
#           * Print contingency table (age_category x study_condition)
#           * Check class balance (n per group, imbalance ratio) 
#       - Run Lefser and print class/age balance warning message
# 3. Save QC and analysis results as RData and csv
# 4. Save Lefser plots for each comparison

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(patchwork)
library(ggplot2)
source("src/datasets.R")
source("src/handler.R")


# Load analysis cohorts --------------------------------------------------------
primary_cohorts <- load_cohorts(primary_cohort_names)


# Define helper functions ------------------------------------------------------
# TODO: write all documentations for helper functions
#' Process and clean/remove NA
clean_cohort <- function(cohort) {
  
  meta <- as.data.frame(colData(cohort))
  meta$study_condition <- as.character(meta$study_condition)
  
  keep <- !is.na(meta$study_condition) &
    meta$study_condition != "" &
    meta$study_condition != " "
  
  cohort <- cohort[, keep]
  
  return(cohort)
}

#'
# TODO: combine clean and standardise into one function?
standardise_conditions <- function(cohort) {
  
  meta <- as.data.frame(colData(cohort))
  
  meta$study_condition <- factor(
    meta$study_condition,
    levels = c(
      "control",
      sort(setdiff(unique(as.character(meta$study_condition)),
                   "control"))
    )
  )
  
  colData(cohort)$study_condition <- meta$study_condition
  
  cohort
}

#' Extract disease groups and check if healthy + 1 disease or >1 disease
#' @param 
#' @return
get_disease_groups <- function(cohort,
                               condition_col = "study_condition",
                               healthy_label = "control") {
  
  meta <- as.data.frame(colData(cohort))
  
  # Extract condition vector, remove NAs and return unique conditions
  conditions <- meta[[condition_col]]
  conditions <- conditions[!is.na(conditions)]
  conditions <- unique(conditions)
  
  # Disease labels = everything except healthy
  diseases <- setdiff(conditions, healthy_label)
  
  list(
    healthy_present = healthy_label %in% conditions,
    diseases = diseases,
    n_diseases = length(diseases)
  )
}

#' Create contingency tables
#' @param 
#' @return
get_contingency_table <- function(cohort,
                                  age_col = "age_category",
                                  condition_col = "study_condition") 
{
  meta <- as.data.frame(colData(cohort))
  
  if (!all(c(age_col, condition_col) %in% colnames(meta))) {
    stop("Required metadata columns not found")
  }
  
  tbl <- table(
    meta[[age_col]],
    meta[[condition_col]],
    useNA = "ifany"
  )
  
  print(tbl)
  
  return(tbl)
}


#' Check class balance and return recommendation
#' @param 
#' @return
check_class_balance <- function(cohort,
                                condition_col = "study_condition",
                                imbalance_mid = 2,
                                imbalance_high = 600) 
  {
  meta <- as.data.frame(colData(cohort))
  counts <- table(meta[[condition_col]])
  
  n_min <- min(counts)
  n_max <- max(counts)
  
  class_ratio <- n_max / n_min
  
  class_status <- if (class_ratio >= imbalance_high) {
    "HIGH IMBALANCE"
  } else if (class_ratio >= imbalance_mid) {
    "MID IMBALANCE"
  } else {
    "OK"
  }
  
  result <- list(
    class_ratio = class_ratio,
    class_status = class_status
  )
  
  print(result)
  
  return(result)
}


#' Check age balance and return recommendation
check_age_balance <- function(cohort,
                              age_col = "age_category",
                              imbalance = 0.3,
                              condition_col = "study_condition") 
  {
  meta <- as.data.frame(colData(cohort))
  age_tab <- table(meta[[age_col]], meta[[condition_col]])
  age_prop <- prop.table(age_tab, margin = 2)
  group_diffs <- combn(ncol(age_prop), 2, function(cols) {
    sum(abs(age_prop[, cols[1]] - age_prop[, cols[2]]))
  })
  
  age_ratio <- max(group_diffs)
  
  age_status <- if (age_ratio > imbalance) {
    "IMBALANCE"
  } else {
    "OK"
    }
  
  result <- list(
    age_ratio = age_ratio,
    age_status = age_status
  )
  
  print(result)
  
  return(result)
}

#' Run Lefser by setting terminal nodes (age_category disabled)
#' @param 
#' @return
run_lefser <- function(cohort,
                      comparison = NULL,
                      subclassCol = "age_category",
                      seed = 1234) {
  meta <- as.data.frame(colData(cohort))
  
  subclass <- NULL
  
  if ("age_category" %in% colnames(meta)) {
    
    age_tab <- table(meta$age_category, meta$study_condition)
    
    # 1. must have at least 2 age groups
    if (nrow(age_tab) < 2) {
      message("No subclass: only one age group")
    } else {
      
      # 2. check each age group has BOTH classes
      valid_strata <- apply(age_tab, 1, function(x) all(x > 0))
      
      if (all(valid_strata)) {
        subclass <- "age_category"
        message("Subclass enabled: age stratification valid")
      } else {
        message("Subclass disabled: incomplete strata (some age groups lack class balance)")
      }
    }
  }
  
  message(
    "Running LEfSe",
    " | samples = ", ncol(cohort),
    " | features = ", nrow(cohort)
  )
  
  # Terminal node filtering
  tn <- get_terminal_nodes(rownames(cohort))
  cohort <- cohort[tn, , drop = FALSE]
  
  # Relative abundance transform
  cohort <- relativeAb(cohort)
  
  # LEfSe
  set.seed(seed)

  tryCatch({
    suppressWarnings(
      suppressMessages(
        lefser(
          cohort,
          classCol = "study_condition",
          subclassCol = subclass
          )
        )
    )
  }, error = function(e) {
    
    message("ERROR in ", cohort, " | ", comparison, ": ", e$message)
    NULL
  })
}

# Storage objects for results --------------------------------------------------
lefser_results <- list()
contingency_tables <- list()
class_balance_summary <- list()
age_balance_summary   <- list()

# TODO: return sample size number too 
analysis_log <- data.frame(
  cohort = character(),
  comparison = character(),
  status = character(),
  reason = character(),
  stringsAsFactors = FALSE
)

split_log <- data.frame(
  cohort = character(),
  unique_condition = character(),
  combination = character(),
  stringsAsFactors = FALSE
)


# Main analysis ----------------------------------------------------------------
for (cohort_name in names(primary_cohorts)) {  
  cohort <- primary_cohorts[[cohort_name]]
  cohort <- clean_cohort(cohort)
  cohort <- standardise_conditions(cohort)
  
  ## extract study_conditions
  info <- get_disease_groups(cohort)
  
  ## filter out cohorts with no healthy controls
  if (!info$healthy_present) {
    message("No healthy controls in ", cohort_name, ". Skipping.")
    analysis_log <- rbind(
      analysis_log,
      data.frame(
        comparison = cohort_name,
        status = "SKIPPED",
        reason = "No healthy controls"
      )
    )
    next
  }
  ## filter out cohorts with healthy-only samples
  if (info$n_diseases == 0) {
    message("Only healthy samples in ", cohort_name, ". Skipping.")
    analysis_log <- rbind(
      analysis_log,
      data.frame(
        comparison = cohort_name,
        status = "SKIPPED",
        reason = "Healthy samples only"
      )
    )
  ## for cohorts with healthy + 1 disease
  } else if (info$n_diseases == 1) {
    disease <- info$diseases
    message("\n", cohort_name, ": healthy vs ", disease)
    
    ### QC: check and save contingency table, class and age balance
    # TODO: one QC function to avoid repeating code?
    comparison_name <- paste(cohort_name, disease, sep = "_")
    
    contingency_table <- get_contingency_table(cohort)
    class_balance <- check_class_balance(cohort)
    age_balance <- check_age_balance(cohort)
    
    class_balance_summary[[comparison_name]] <- class_balance
    age_balance_summary[[comparison_name]]   <- age_balance
    contingency_tables[[comparison_name]] <- contingency_table
    
    ### Record split log - confirm no of disease conditions = no of comparisons
    # TODO: create function
    meta_sub <- as.data.frame(colData(cohort))
    
    split_log <- rbind(
      split_log,
      data.frame(
        cohort = cohort_name,
        unique_condition = paste(sort(unique(meta_sub$study_condition)), collapse = ";"),
        combination = paste("control", disease, sep = "+")
      )
    )
    
    ### run Lefser if class_status is not "SKIP"
    if (class_balance$class_status != "HIGH IMBALANCE") {
      res <- run_lefser(cohort, comparison = paste("control vs", disease))
      
      if (is.null(res)) {
        analysis_log <- rbind(
          analysis_log,
          data.frame(
            comparison = comparison_name,
            status = "FAILED",
            reason = "LEfSe returned NULL"
          )
        )
        
      } else {
        analysis_log <- rbind(
          analysis_log,
          data.frame(
            comparison = comparison_name,
            status = "SUCCESS",
            reason = paste("n_features =", nrow(res))
          )
        )
        
        lefser_results[[paste(cohort_name, disease, sep = "_")]] <- res
      }
      
    } else {
      message("Skipping Lefser due to class imbalance")
    }
    
    

    ## for cohorts with healthy + >1 disease
  } else {
    message("\n", cohort_name, ": healthy + ", info$n_diseases, " diseases")
    
    for (disease in info$diseases) {
      message("\n", cohort_name, " subset: healthy vs ", disease)
      
      ### create subset of cohort (healthy + 1 disease)
      meta <- as.data.frame(colData(cohort))
      keep <- meta$study_condition %in% c("control", disease)
      cohort_subset <- cohort[, keep]
      cohort_subset <- standardise_conditions(cohort_subset)
      
      meta_sub <- as.data.frame(colData(cohort_subset))
      
      split_log <- rbind(
        split_log,
        data.frame(
          cohort = cohort_name,
          unique_condition = paste(sort(unique(meta_sub$study_condition)), collapse = ";"),
          combination = paste("control", disease, sep = "+")
        )
      )
      
      ### QC: check and save contingency table, class and age balance
      comparison_name <- paste(cohort_name, disease, sep = "_")
      
      contingency_table <- get_contingency_table(cohort_subset)
      class_balance <- check_class_balance(cohort_subset)
      age_balance <- check_age_balance(cohort_subset)
      
      class_balance_summary[[comparison_name]] <- class_balance
      age_balance_summary[[comparison_name]]   <- age_balance
      contingency_tables[[comparison_name]] <- contingency_table
      
      if (class_balance$class_status != "HIGH IMBALANCE") {
        res <- run_lefser(cohort_subset, comparison = paste("control vs", disease))
        
        if (is.null(res)) {
          analysis_log <- rbind(
            analysis_log,
            data.frame(
              comparison = comparison_name,
              status = "FAILED",
              reason = "LEfSe returned NULL"
            )
          )
          
        } else {
          analysis_log <- rbind(
            analysis_log,
            data.frame(
              comparison = comparison_name,
              status = "SUCCESS",
              reason = paste("n_features =", nrow(res))
            )
          )
          
          lefser_results[[paste(cohort_name, disease, sep = "_")]] <- res
        }
        
      } else {
        message("Skipping Lefser due to class imbalance")
        }
    }
  }
}


# Results and visualisation ----------------------------------------------------
# Set output directory
out_dir <- "results/discriminant_taxa"

# Save as RData
save(
  lefser_results,
  contingency_tables,
  class_balance_summary,
  age_balance_summary,
  file = file.path(out_dir, "analysis_results.RData")
)

# Process lefser_results, contingency_tables, class_balance_summary and age_balance_summary 
# TODO: make these into helper functions
lefser_df <- do.call(
  rbind,
  lapply(names(lefser_results), function(nm) {
    
    res <- lefser_results[[nm]]
    
    if (is.null(res) || nrow(res) == 0)
      return(NULL)
    
    df <- as.data.frame(res)
    df$comparison <- nm
    
    df
  })
)

contingency_df <- do.call(
  rbind,
  lapply(names(contingency_tables), function(comp) {
    
    tab <- contingency_tables[[comp]]
    
    df <- as.data.frame(tab)
    
    df$comparison <- comp
    
    df
  })
)

class_balance_df <- do.call(
  rbind,
  lapply(names(class_balance_summary), function(x) {
    data.frame(
      comparison   = x,
      class_ratio  = class_balance_summary[[x]]$class_ratio,
      class_status = class_balance_summary[[x]]$class_status
    )
  })
)

age_balance_df <- do.call(
  rbind,
  lapply(names(age_balance_summary), function(x) {
    data.frame(
      comparison = x,
      age_ratio  = age_balance_summary[[x]]$age_ratio,
      age_status = age_balance_summary[[x]]$age_status
    )
  })
)

# Save results as csv files
# TODO: again, one function
write.csv(
  lefser_df,
  file.path(out_dir, "lefser_results_combined.csv"),
  row.names = FALSE
)

write.csv(
  contingency_df,
  file.path(out_dir, "contingency_tables.csv"),
  row.names = FALSE
)

summary_df <- merge(
  class_balance_df,
  age_balance_df,
  by = "comparison"
)

write.csv(
  summary_df,
  file.path(out_dir, "class_age_summary.csv"),
  row.names = FALSE
)

write.csv(
  split_log,
  file.path(out_dir, "split_log.csv"),
  row.names = FALSE
)

write.csv(
  analysis_log,
  file.path(out_dir, "analysis_log.csv"),
  row.names = FALSE
)


# Create Lefser plots and save
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

# TODO: make plots from csv and not use lefserplots

for (name in names(lefser_results)) {
  
  res <- lefser_results[[name]]
  
  if (is.null(res)) next
  
  if (nrow(res) == 0) {
    message("Skipping empty result: ", name)
    next
  } else if (nrow(res) < 3) {
    message ("Skipping <1 significant feature in", name)
    next
  }
  
  p <- lefserPlot(res) +
    theme(
      plot.margin = margin(10, 40, 10, 60)
    )
  
  safe_name <- gsub("[/\\\\]", "_", name)
  
  ggsave(
    filename = file.path(out_dir, "plots", paste0(safe_name, ".png")),
    plot = p,
    bg = "white",
    width = 12,
    height = 8,
    dpi = 300
  )
}

# TODO: statistifcally analyse features in health/disease
