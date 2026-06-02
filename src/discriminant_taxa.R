# Pipeline:
# 1. Read dataset registry from datasets.R and load cohorts via handler.R
# 2. For each cohort, extract study conditions:
#    * If healthy + 1 disease  -> 
#       - Lefser analysis loop which:
#            * Prints contingency table (age_category x study_condition)
#            * Checks class balance (n per group, imbalance ratio) 
#            * Runs Lefser if class balance OK/CAUTION
#    * If healthy + >1 disease -> 
#       - Create pair-wise subsets (healthy + 1 disease)
#       - For each subset, perform Lefser analysis loop (as above)
# 3. Collect Lefser result tables and plots
# 4. Print class balance summary across all comparisons
# 5. Save and display combined plot via handler.R


# Load packages and dependencies ###############################################
library(lefser)
library(curatedMetagenomicData)
library(patchwork)
library(ggplot2)
source("datasets.R")
source("handler.R")


# Load analysis cohorts ########################################################
primary_cohorts <- load_cohorts(primary_cohort_names)


# Define helper functions ######################################################

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


#' Check class balance and return recommendation
#' @param 
#' @return
check_class_balance <- function(cohort,
                                condition_col = "study_condition",
                                imbalance_warning = 2,
                                imbalance_skip = 5) 
  {
  meta <- as.data.frame(colData(cohort))
  counts <- table(meta[[condition_col]])
  
  n_min <- min(counts)
  n_max <- max(counts)
  
  class_ratio <- n_max / n_min
  
  class_status <- if (class_ratio >= imbalance_skip) {
    "SKIP"
  } else if (class_ratio >= imbalance_warning) {
    "CAUTION"
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
                              condition_col = "study_condition") 
  {
  meta <- as.data.frame(colData(cohort))
  age_tab <- table(meta[[age_col]], meta[[condition_col]])
  age_prop <- prop.table(age_tab, margin = 2)
  group_diffs <- combn(ncol(age_prop), 2, function(cols) {
    sum(abs(age_prop[, cols[1]] - age_prop[, cols[2]]))
  })
  
  age_ratio <- max(group_diffs)
  
  age_status <- if (age_ratio > 0.3) {
    "WARNING"
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
  
  
standardise_conditions <- function(cohort) {
  
  meta <- as.data.frame(colData(cohort))
  
  meta$study_condition <- as.character(meta$study_condition)
  
  meta$study_condition[meta$study_condition == "control"] <- "control"
  meta$study_condition[is.na(meta$study_condition)] <- NA
  
  meta$study_condition <- factor(meta$study_condition,
                                 levels = c("control",
                                            sort(setdiff(unique(meta$study_condition), "control"))))
  
  colData(cohort)$study_condition <- meta$study_condition
  
  cohort
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
      message("Only one age group → no subclass")
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
    "\nRunning LEfSe | ", comparison,
    " | samples = ", ncol(cohort),
    " | features = ", nrow(cohort)
  )
  
  # -----------------------------
  # 1. Terminal node filtering
  # -----------------------------
  tn <- get_terminal_nodes(rownames(cohort))
  cohort <- cohort[tn, , drop = FALSE]
  
  # -----------------------------
  # 2. Relative abundance transform
  # -----------------------------
  cohort <- relativeAb(cohort)
  
  # -----------------------------
  # 4. LEfSe
  # -----------------------------
  set.seed(seed)
  
  tryCatch({
    
    suppressWarnings(
      lefser(
        cohort,
        classCol = "study_condition",
        subclassCol = subclass
      )
    )
    
  }, error = function(e) {
    
    message("ERROR in ", cohort, " | ", comparison, ": ", e$message)
    NULL
  })
}


# Main discriminant taxa analysis loop #########################################
lefser_results <- list()
lefser_plots          <- list()
contingency_tables <- list()
balance_log        <- list()


for (cohort_name in names(primary_cohorts)) {  
  cohort <- primary_cohorts[[cohort_name]]
  cohort <- clean_cohort(cohort)
  cohort <- standardise_conditions(cohort)
  
  info <- get_disease_groups(cohort)
  
  if (!info$healthy_present) {
    message("No healthy controls in ", cohort_name, ". Skipping.")
    next
  }

  if (info$n_diseases == 0) {
    message("Healthy samples only in ", cohort_name, ". Skipping.")
  } 
  else if (info$n_diseases == 1) {
      disease <- info$diseases
      message("\n", cohort_name, ": healthy + ", disease)
      
      # Create and print contingency table
      get_contingency_table(cohort)

      # Check class and age balance
      class_balance <- check_class_balance(cohort)
      check_age_balance(cohort)

      # If class_status = OK, run Lefser
      if (class_balance$class_status != "SKIP") {
        res <- run_lefser(cohort,
                        comparison = paste("control vs", disease)
      )
      
      lefser_results[[paste(cohort_name, disease, sep = "_")]] <- res
      } else {
        message("Skipping Lefser due to class imbalance")
      }

      
  } else {
    message(cohort_name, ": healthy + ", info$n_diseases, " diseases")

    for (disease in info$diseases) {
      message("\n  Comparison: healthy vs ", disease)
      
      # Create subset of cohort
      meta <- as.data.frame(colData(cohort))
      keep <- meta$study_condition %in% c("control", disease)
      cohort_subset <- cohort[, keep]
      cohort_subset <- standardise_conditions(cohort_subset)
      
      # Create and print contingency table
      get_contingency_table(cohort_subset)

      # Check class and age balance
      class_balance <- check_class_balance(cohort_subset)
      check_age_balance(cohort_subset)
      
      # If class_status = OK, run Lefser
      if (class_balance$class_status != "SKIP") {
        res <- run_lefser(cohort_subset,
                          comparison = paste("control vs", disease)
        )
        
        lefser_results[[paste(cohort_name, disease, sep = "_")]] <- res
      } else {
        message("Skipping Lefser due to class imbalance")
      }
    }
  }
} 



# Results and visualisation ####################################################
outdir <- file.path("results/discriminant_taxa")

dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

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
    filename = file.path(outdir, "plots", paste0(safe_name, ".png")),
    plot = p,
    bg = "white",
    width = 12,
    height = 8,
    dpi = 300
  )
}
