# Conduct LEfSe analysis using lefser package
# 1. Define output directory and functions
# 2. Load data using utils.R
# 3. Execute main functions:
#   a. run_analysis()
#     * Clean, extract disease and validate conditions
#     * ANALYSIS BRANCH 1 (control x 1 disease) -> analysis_pipeline() executes:
#           - qc_cohort() 
#           - lefse_analyse() 
#           - log_qc()
#           - log_analysis()
#           - log_result()
#     * ANALYSIS BRANCH 2 (control x >1 disease)  -> 
#           - Create pair-wise subsets
#           - Execute analysis_pipeline()
#   b. export_analysis()
#     * Save lefse_results, contingency_tables, qc_summary and analysis_log  
#       as out_dir/[log].csv
#     * Save bulk lefse_result object as out_dir/lefse_results.rds

# TODO collapse branches via comparisons
# TODO filter by sample size before analysis
# TODO compile contingency table and QC into one output
# TODO update documentation after combining with lefse_plot.R
# TODO conduct at genus level instead of species

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(lefser)
library(ggplot2)
library(dplyr)
library(stringr)
library(rlang)
source("R/utils.R")


# Define output directory ------------------------------------------------------
out_dir <- "results/02_lefse"


# Define analysis functions ----------------------------------------------------
# Run QC on a cohort and return a contingency table, class imbalance and 
# age imbalance ratios and counts
qc_cohort <- function(
    cohort,
    cohort_name,
    disease,
    age_col = "age_category",
    condition_col = "study_condition",
    class_imbalance_mid = 2,
    class_imbalance_high = 5,
    class_imbalance_severe = 10,
    age_imbalance_thresh = 0.25,
    n_status_thresh = 10) {
  
  meta <- as.data.frame(colData(cohort))
  
  # Ensure metadata exist
  if (!all(c(age_col, condition_col) %in% colnames(meta))) {
    stop("Required metadata columns not found")
  }
  
  # Create contingency table data frame
  contingency_df <- as.data.frame(as.table(
    table(
      factor(meta[[age_col]]),
      factor(meta[[condition_col]]),
      useNA = "ifany")))
  
  names(contingency_df) <- c("age_category", "condition", "count")
  
  contingency_df$cohort <- cohort_name

  contingency_df <- contingency_df[, c(
    "cohort", "age_category", "condition", "count"
  )]
  
  # Check class balance
  counts <- table(meta[[condition_col]])
  
  n_min <- min(counts)
  n_max <- max(counts)
  
  class_ratio <- n_max / n_min
  
  class_imbalance <- if (class_ratio > class_imbalance_severe) {
    paste0("severe (>", class_imbalance_severe, ")")
  } else if (class_ratio > class_imbalance_high) {
    paste0("high (>", class_imbalance_high, ")")
  } else if (class_ratio > class_imbalance_mid) {
    paste0("moderate (>", class_imbalance_mid, ")")
  } else {
    "low"
  }
  
  class_balance <- list(
    class_ratio = class_ratio,
    class_imbalance = class_imbalance)
  
  # Check age balance
  age_tab <- table(meta[[age_col]], meta[[condition_col]])
  age_prop <- prop.table(age_tab, margin = 2)
  
  group_diffs <- combn(ncol(age_prop), 2, function(cols) {
    sum(abs(age_prop[, cols[1]] - age_prop[, cols[2]]))
  })
  
  age_ratio <- max(group_diffs)
  
  age_imbalance <- if (age_ratio > age_imbalance_thresh) {
    paste0("high (>", age_imbalance_thresh, ")")
  } else {
    "low"
  }
  
  age_balance <- list(
    age_ratio = age_ratio,
    age_imbalance = age_imbalance)
  
  # Get counts
  counts <- as.list(table(meta[[condition_col]]))
  counts <- as.list(table(meta[[condition_col]]))
  
  n_controls <- counts[["control"]] %||% 0
  n_conditions <- counts[[disease]] %||% 0
  
  n_status <- character()
  
  if (n_controls < n_status_thresh) {
    n_status <- c(n_status, paste0("n controls < ", n_status_thresh))
  }
  
  if (n_conditions < 5) {
    n_status <- c(n_status, paste0("n ", disease, " < ", n_status_thresh))
  }
  
  n_status <- if (length(n_status) == 0) {
    NA_character_
  } else {
    paste(n_status, collapse = " | ")
  }
  
  
  # Return results
  result <- list(
    contingency_table = contingency_df,
    class_balance = class_balance,
    age_balance = age_balance,
    n_controls = n_controls,
    n_conditions = n_conditions,
    n_status = n_status)
  
  return(result)
}

# Analyse a cohort by performing:
# - Terminal node filtering
# - Relative abundance transformation
# - Subclass stratification only if:
#    * age_category exists
#    * at least 2 age groups are present
#    * all age groups contain both control and disease samples
# - LEfSe differential abundance analysis
lefse_analyse <- function(
    cohort,
    comparison = NULL,
    subclassCol = "age_category",
    disease_label = NULL,
    control_label = "control",
    seed = 1234) {
  
  meta <- colData(cohort)
  
  subclass <- NULL
  
  # Enforce ordering of control and disease
  meta$study_condition <- factor(
    meta$study_condition,
    levels = c(control_label, disease_label))
  
  cohort$study_condition <- meta$study_condition
  
  # Enable/disable subclass: age_category when:
  if ("age_category" %in% colnames(meta)) {
    
    age_tab <- table(meta$age_category, meta$study_condition)
    
    ## there is more than one age group
    if (nrow(age_tab) < 2) {
      message("Subclass disabled: only one age group")
      
    } else {
      
      ## and each age group has both control x disease samples
      valid_strata <- apply(age_tab, 1, function(x) all(x > 0))
      
      if (all(valid_strata)) {
        subclass <- "age_category"
        message("Subclass enabled: age stratification valid")
      } else {
        message("Subclass disabled: some age groups lack class balance")
      }
    }
  }
  
  message(
    "Running LEfSe",
    " | samples = ", ncol(cohort),
    " | features = ", nrow(cohort)
  )
  
  # Filter terminal node 
  tn <- get_terminal_nodes(rownames(cohort))
  cohort <- cohort[tn, , drop = FALSE]
  
  # Aggregate species to genus
  abund <- SummarizedExperiment::assay(cohort)
  
  genus <- sub(".*\\|g__", "g__", rownames(abund))
  genus <- sub("\\|.*", "", genus)
  
  keep <- grepl("^g__", genus)
  
  abund <- abund[keep, , drop = FALSE]
  genus <- genus[keep]
  
  if (nrow(abund) == 0) {
    message("ERROR: no genus-level features found")
    return(NULL)
  }
  
  abund_genus <- rowsum(
    abund,
    group = genus)
  
  # Rebuild object
  cohort <- SummarizedExperiment::SummarizedExperiment(
    assays = list(abundance = abund_genus),
    colData = meta)
  
  message(
    "After genus aggregation: features = ",
    nrow(cohort))
  
  # Transform relative abundance
  cohort <- relativeAb(cohort)
  
  # Run LEfSe
  set.seed(seed)
  
  tryCatch({
    suppressWarnings(
      suppressMessages(
        lefser(
          cohort,
          classCol = "study_condition",
          subclassCol = subclass)))
    
  }, error = function(e) {
    message("ERROR in ", comparison, ": ", e$message)
    NULL
  })
}

#' Log QC results from a single comparison including contingency table, class 
#' ratio, class imbalance, age ratio, age imbalance and condition counts, append
#' analysis_res$qc_summary and return analysis_res
log_qc <- function(
    comparison_name, 
    cohort_name, 
    disease, 
    qc, 
    analysis_res) {
  
  # initialise containers if missing
  if (is.null(analysis_res$qc_summary)) {
    analysis_res$qc_summary <- list()
  }
  
  if (is.null(analysis_res$contingency_tables)) {
    analysis_res$contingency_tables <- list()
  }

  row <- data.frame(
    comparison = comparison_name,
    cohort = cohort_name,
    disease = disease,
    
    class_ratio = qc$class_balance$class_ratio,
    class_imbalance = qc$class_balance$class_imbalance,
    
    age_ratio = qc$age_balance$age_ratio,
    age_imbalance = qc$age_balance$age_imbalance,
    
    n_controls = qc$n_controls,
    n_conditions = qc$n_conditions,
    n_status = qc$n_status,
    
    stringsAsFactors = FALSE)
  
  analysis_res$qc_summary[[comparison_name]] <- row
  
  analysis_res$contingency_tables[[comparison_name]] <- 
    qc$contingency_table
  
  return(analysis_res)
}

#' Log LEfSe analysis status for a single comparison, append 
#' analysis_res$analysis_log and return analysis_res
log_analysis <- function(
    comparison_name, 
    cohort_name, 
    disease, 
    status = NULL,
    reason = NULL,
    result, 
    analysis_res) {
  
  # initialise container if missing
  if (is.null(analysis_res$analysis_log)) {
    analysis_res$analysis_log <- list()
  }
  
  create_log <- function(status, reason) {
    data.frame(
      comparison = comparison_name,
      cohort = cohort_name,
      disease = disease,
      status = status,
      reason = reason,
      stringsAsFactors = FALSE)
  }
  
  # Case 1: explicit SKIPPED when condition_invalid == TRUE
  if (!is.null(status) && status == "SKIPPED") {
    analysis_res$analysis_log[[comparison_name]] <-
      create_log("SKIPPED", reason)
    return(analysis_res)
  }
  
  # Case 2: null result
  if (is.null(result)) {
    analysis_res$analysis_log[[comparison_name]] <-
      create_log("FAILED", "LEfSe returned NULL")
    return(analysis_res)
  }
  
  # Case 3: coercion
  result_df <- tryCatch(
    as.data.frame(result),
    error = function(e) NULL)
  
  if (is.null(result_df)) {
    analysis_res$analysis_log[[comparison_name]] <-
      create_log("FAILED", "Could not coerce LEfSe output to data frame")
    return(analysis_res)
  }
  
  # Case 4: success
  analysis_res$analysis_log[[comparison_name]] <-
    create_log("SUCCESS", paste("n_features =", nrow(result_df)))
  
  return(analysis_res)
}

#' Log LEfSe result from a single analysis, append analysis_res$lefse_results
#' and return analysis_res
log_result <- function(
    comparison_name, 
    cohort_name, 
    disease, 
    result, 
    analysis_res) {
  
  # initialise container if missing
  if (is.null(analysis_res$lefse_results)) {
    analysis_res$lefse_results <- list()
  }
  
  # Store metadata
  row <- data.frame(
    comparison = comparison_name,
    cohort = cohort_name,
    disease = disease,
    stringsAsFactors = FALSE)
  
  # No result
  if (is.null(result)) {
    return(analysis_res)
  }
  
  # Safe coercion
  result_df <- tryCatch(
    as.data.frame(result),
    error = function(e) NULL)

  # Failed coercion or empty
  if (is.null(result_df) || nrow(result_df) == 0) {
    return(analysis_res)
  }
  
  # Combine and store data frames
  combined <- cbind(row[rep(1, nrow(result_df)), , drop = FALSE], result_df)
  analysis_res$lefse_results[[comparison_name]] <- combined
  
  return(analysis_res)
}

#' Run analysis pipeline which includes QC checks, LEfSe analysis and logging
analysis_pipeline <- function(
    cohort,
    cohort_name,
    comparison_name,
    disease,
    analysis_res) {
  
  # Run and log QC
  qc <- qc_cohort(
    cohort,
    cohort_name = cohort_name, 
    disease = disease)
  
  analysis_res <- log_qc(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    qc = qc,
    analysis_res = analysis_res)
  
  # Skip LEfSe if either group has too few samples
  if (!is.na(qc$n_status)) {
    
    analysis_res <- log_analysis(
      comparison_name = comparison_name,
      cohort_name = cohort_name,
      disease = disease,
      status = "SKIPPED",
      reason = qc$n_status,
      result = NULL,
      analysis_res = analysis_res)
    
    message("Skipping: sample size < 10")
    
    return(analysis_res)}
  
  # Run and log LEfSe
  result <- lefse_analyse(
    cohort,
    comparison = paste("control vs", disease),
    disease_label = disease)
  
  analysis_res <- log_analysis(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    result = result,
    analysis_res = analysis_res)
  
  # Store results in analysis
  analysis_res <- log_result(
    comparison_name = comparison_name,
    cohort_name = cohort_name,
    disease = disease,
    result = result,
    analysis_res = analysis_res)
  
  return(analysis_res)
}


# Define plotting functions ----------------------------------------------------
# Extract species label from full clade string
extract_taxon <- function(df, rank = c("genus", "species", "family")) {
  
  rank <- match.arg(rank)
  
  prefixes <- c(
    species = "s__",
    genus   = "g__",
    family  = "f__")
  
  prefix <- prefixes[[rank]]
  
  df %>%
    dplyr::mutate(
      taxon = dplyr::case_when(
        stringr::str_detect(features, prefix) ~
          sub(paste0(".*", prefix), prefix, features),
        TRUE ~ features
      ),
      taxon = sub("\\|.*", "", taxon),
      taxon = sub("^[a-z]__", "", taxon),
      taxon = gsub("_", " ", taxon),
      taxon = sub("^[A-Za-z]+:", "", taxon))
}

#' Convert disease label string with underscores to title case
convert_smart_case <- function(x) {
  x <- gsub("_", " ", x)
  
  words <- strsplit(x, " ")[[1]]
  
  words <- ifelse(
    words == toupper(words),
    words,
    stringr::str_to_title(tolower(words)))
  
  paste(words, collapse = " ")
}

#' Prepare LEfSe plot data by filtering scores and assigning control/disease label
process_scores <- function(
    lefse_res_df,
    disease_label,
    control_label,
    max_features) {
  
  lefse_res_df %>%
    dplyr::mutate(scores = as.numeric(scores)) %>%
    dplyr::filter(!is.na(scores)) %>%
    dplyr::slice_max(
      order_by = abs(scores),
      n = max_features,
      with_ties = FALSE) %>%
    
    dplyr::mutate(
      group_label = ifelse(
        scores > 0,
        disease_label,
        control_label),
      
      group_label = factor(
        group_label,
        levels = c(control_label, disease_label)))
}

#' Log plotting status and number of features
log_plot <- function(
    lefse_res_df,
    comparison,
    # TODO remove repeated code
    min_features = 2,
    max_features = 6) {
  
  create_log <- function(
    comp,
    n_features,
    plotted_features = NA,
    status,
    reason = NA,
    df_plot = NULL) {
    
    list(
      log = data.frame(
        comparison = comp,
        n_features = n_features,
        plotted_features = plotted_features,
        status = status,
        reason = reason,
        stringsAsFactors = FALSE
      ),
      df_plot = df_plot
    )
  }
  
  n_features_raw <- nrow(lefse_res_df)
  
  # Case 1: no data
  if (n_features_raw == 0) {
    return(create_log(comparison, 0, NA, "SKIPPED", "No data"))
  }
  
  disease_label <- convert_smart_case(unique(lefse_res_df$disease)[1])
  control_label <- "Control"
  
  # Case 2: too few features
  if (n_features_raw < min_features) {
    return(create_log(
      comparison,
      n_features_raw,
      NA,
      "SKIPPED",
      "Too few features"
    ))
  }
  
  # Process and filter scores for plotting
  df_plot <- process_scores(
    lefse_res_df,
    disease_label,
    control_label,
    max_features
  )
  
  n_features <- nrow(df_plot)
  
  # Case 3: no valid scores after processing
  if (n_features == 0) {
    return(create_log(
      comparison,
      n_features_raw,
      0,
      "SKIPPED",
      "No valid scores"
    ))
  }
  
  # Case 4: set status as "READY" for plotting
  create_log(
    comparison,
    n_features_raw,
    n_features,
    "READY",
    NA,
    df_plot = df_plot
  )
}

# Build a LEfSe bar plot using ggplot2
create_plot <- function(
    lefse_res_df,
    cohort_name,
    disease_label,
    control_label) {
  
  ggplot(
    lefse_res_df,
    aes(
      x = scores,
      y = stats::reorder(taxon, scores, FUN = mean),
      fill = group_label)) +
    
    geom_col(
      width = 0.7,
      colour = "white",
      linewidth = 0.2) +
    
    geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      colour = "grey40") +
    
    scale_fill_manual(
      values = c(
        setNames("#c5cb93ff", control_label),
        setNames("#dfa57eff", disease_label)),
      labels = c("Control", "Case"),
      name = "Enriched in") +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.05))) +
    
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0),
      
      axis.text.y = element_text(
        size = 20,
        face = "italic"),
      
      axis.title = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      
      legend.position = "bottom",
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()) +
    
    labs(
      title = paste0(
        cohort_name,
        " (",
        control_label,
        " vs ",
        disease_label,
        ")"
      ),
      x = "LDA score (log10)",
      y = NULL)
}

# Save a single plot as a pdf
save_plot <- function(
    plot,
    comparison,
    out_dir,
    n_features) {
  
  safe_comp <- gsub("[/\\\\]", "", comparison)
  
  out_path <- file.path(
    out_dir,
    paste0(safe_comp, ".pdf"))
  
  plot_height <- max(4, n_features * 0.45 + 1.5)
  
  ggsave(
    filename = out_path,
    plot = plot,
    width = 7,
    height = plot_height,
    device = "pdf")
  
  invisible(out_path)
}

# Define main functions --------------------------------------------------------
#' Run main LEfSe loop for two analysis branches
run_analysis <- function(primary_cohorts, out_dir) {
  
  # Initialise storage object
  analysis_res <- list()
  
  for (cohort_name in names(primary_cohorts)) {  
    
    # Process and extract study_conditions from cohort
    cohort <- primary_cohorts[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check cohort contains both control and disease samples
    condition_invalid <- validate_conditions(info, cohort_name)
    
    if (!is.null(condition_invalid)) {
      analysis_res <- log_analysis(
        comparison_name,
        cohort_name,
        disease,
        status = "SKIPPED",
        reason = condition_invalid,
        analysis_res = analysis_res)
      next
      
    # ANALYSIS BRANCH 1: control x 1 disease
    } else if (info$n_diseases == 1) {
      
      disease = info$diseases
      
      comparison_name <- paste(cohort_name, disease, sep = "_")
      message("\n", cohort_name, ": control vs ", info$diseases)
      
      analysis_res <- analysis_pipeline(
        cohort = cohort,
        cohort_name = cohort_name,
        disease = disease,
        comparison_name = comparison_name,
        analysis_res = analysis_res)
      
    # ANALYSIS BRANCH 2: control x >1 disease
    } else {
      
      message("\n", cohort_name, ": control vs ", info$n_diseases, " diseases")
      
      for (disease in info$diseases) {
        
        message("\n", cohort_name, " subset: control vs ", disease)
        
        # Subset cohort for control x 1 disease
        meta <- as.data.frame(colData(cohort))
        keep <- meta$study_condition %in% c("control", disease)
        cohort_subset <- cohort[, keep]

        comparison_name <- paste(cohort_name, disease, sep = "_")
        
        analysis_res <- analysis_pipeline(cohort = cohort_subset,
                                 cohort_name = cohort_name,
                                 disease = disease,
                                 comparison_name = comparison_name,
                                 analysis_res = analysis_res)
        }
      }
  }
  sapply(analysis_res$contingency_tables, function(x) {
    if (is.null(x)) return("NULL")
    if (length(x) == 0) return("EMPTY")
    class(x)
  })
  lapply(analysis_res$contingency_tables, dim)
  
  return(analysis_res)
  
  # Export results
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save full object
  saveRDS(analysis_res, file.path(out_dir, "lefse_results.rds"))
  
  # Bind data frames
  analysis_res$analysis_log <- dplyr::bind_rows(analysis_res$analysis_log)
  analysis_res$qc_summary <- dplyr::bind_rows(analysis_res$qc_summary)
  analysis_res$contingency_df <- dplyr::bind_rows(analysis_res$contingency_tables)
  analysis_res$lefse_df <- dplyr::bind_rows(analysis_res$lefse_results)
  
  # Write CSVs
  write_csvs(
    data_list = list(
      lefse_results = analysis_res$lefse_df,
      contingency_tables = analysis_res$contingency_df,
      qc_summary = analysis_res$qc_summary,
      analysis_log = analysis_res$analysis_log),
    out_dir = out_dir)
}

run_plotting <- function(
    analysis_res,
    out_dir,
    min_features = 2,
    max_features = 10) {
  
  plot_dir <- file.path(out_dir, "plots")
  dir.create(file.path(plot_dir), recursive = TRUE,
             showWarnings = FALSE)
  
  # Convert list of comparisons into dataframe
  lefse_res <- dplyr::bind_rows(
    analysis_res$lefse_results)
  
  lefse_res <- extract_taxon(lefse_res, rank = "genus")
  
  plot_log <- lapply(
    unique(lefse_res$comparison),
    function(comp) {
      
      lefse_res_df <- filter(
        lefse_res,
        comparison == comp
      )
      
      res <- log_plot(
        lefse_res_df = lefse_res_df,
        comparison = comp,
        min_features = min_features,
        max_features = max_features
      )
      
      if (res$log$status == "READY") {
        
        cohort_name <- unique(lefse_res_df$cohort)[1]
        disease_label <- convert_smart_case(unique(lefse_res_df$disease)[1])
        
        p <- create_plot(
          res$df_plot,
          cohort_name,
          disease_label,
          "Control")
        
        save_plot(
          p,
          comp,
          plot_dir,
          nrow(res$df_plot))
        
        res$log$status <- "PLOTTED"
      }
      
      res$log
    })
  
  plot_log_df <- bind_rows(plot_log)
  
  write.csv(
    plot_log_df,
    file.path(plot_dir, "_plot_log.csv"),
    row.names = FALSE)
  message("Done. Plots written to: ", plot_dir)
  
  invisible(plot_log_df)
}



# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary.rds")


# Execute ----------------------------------------------------------------------
# analysis_res <- run_analysis(primary_cohorts, out_dir)
run_plotting(analysis_res, out_dir)
