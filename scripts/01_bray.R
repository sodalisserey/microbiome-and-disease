# Perform PERMANOVA analysis using vegan package
# 1. Define output directory and functions
# 2. Load data using utils.R
# 3. Execute main function: run_bray_curtis()
#   * Clean, extract disease and validate conditions
#   * Execute calculate_beta_diversity() and log_beta_diversity()
#   * Export bray_res as out_dir/[cohort_name]_plot.pdf


# Load packages and dependencies -----------------------------------------------
library(vegan)
library(ggplot2)
library(readr)
library(SummarizedExperiment)
source("R/utils.R")

# Define output directory ------------------------------------------------------
out_dir <- "results/01_bray"


# Define helper functions ------------------------------------------------------
# TODO flesh this out if timing permits
# calculate_alpha_diversity <- function(
#     cohort,
#     cohort_name,
#     condition_col = "study_condition") {
#   
#   # Extract abundance matrix
#   abund <- as.matrix(assay(cohort))
#   
#   # Extract metadata
#   meta <- as.data.frame(colData(cohort))
#   
#   # Ensure sample IDs are available
#   if (is.null(rownames(meta))) {
#     rownames(meta) <- colnames(cohort)
#   }
#   
#   # Match samples
#   common_samples <- intersect(colnames(abund), rownames(meta))
#   
#   if (length(common_samples) < 3) {
#     stop("Too few samples")
#   }
#   
#   abund <- abund[, common_samples, drop = FALSE]
#   meta <- meta[common_samples, , drop = FALSE]
#   
#   # Convert to relative abundance
#   abund_rel <- sweep(
#     abund,
#     2,
#     colSums(abund),
#     "/")
#   
#   # Shannon diversity
#   shannon <- diversity(
#     t(abund_rel),
#     index = "shannon")
#   
#   # Return results
#   data.frame(
#     sample = names(shannon),
#     cohort = cohort_name,
#     shannon = shannon,
#     condition = meta[[condition_col]],
#     stringsAsFactors = FALSE)
# }


# Calculate beta-diversity, run PERMANOVA and save plot
calculate_beta_diversity <- function(
    cohort,
    cohort_name,
    condition_col = "study_condition",
    permutations = 9999,
    dispersion_permutations = 999,
    out_dir) {
  
  plot_dir <- file.path(out_dir, "plots")
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  warnings_list <- NA
  error_msg <- NA
  
  result <- tryCatch({
    
    # Extract abundance
    abund <- as.matrix(assay(cohort))
    
    # Extract metadata
    meta <- as.data.frame(colData(cohort))
    
    # Ensure rownames are sample IDs
    if (is.null(rownames(meta))) {
      rownames(meta) <- colnames(cohort)
    }
    
    # Align samples correctly
    common_samples <- intersect(colnames(abund), rownames(meta))
    
    if (length(common_samples) < 3) {
      stop("Too few overlapping samples")
    }
    
    abund <- abund[, common_samples, drop = FALSE]
    meta  <- meta[common_samples, , drop = FALSE]
    
    message(" > Computing Bray-Curtis dissimilarity")
    X <- t(abund)
    bray <- vegdist(X, method = "bray")
    
    message(" > Running PERMANOVA")
    meta[[condition_col]] <- as.factor(meta[[condition_col]])
    
    ad <- adonis2(
      bray ~ meta[[condition_col]],
      permutations = permutations)
    
    message(" > Checking homogeneity of dispersion")
    meta[[condition_col]] <- factor(
      meta[[condition_col]],
      levels = c(
        "control",
        setdiff(unique(meta[[condition_col]]), "control")))
    
    bd <- betadisper(bray, meta[[condition_col]])
    bd_test <- permutest(bd, permutations = dispersion_permutations)
    
    list(
      adonis = ad,
      betadisper = bd,
      dispersion_test = bd_test,
      meta = meta)
    
  }, warning = function(w) {
    warnings_list <<- c(warnings_list, w$message)
    invokeRestart("muffleWarning")
    
  }, error = function(e) {
    error_msg <<- e$message
    return(NULL)
  })
  
  pdf(
    file.path(plot_dir, paste0(cohort_name, "_plot.pdf")),
    width = 6,
    height = 5)

  # Extract eigenvalues
  eig <- bd$eig
  eig_pos <- eig[eig > 0]
  pc_var <- eig_pos / sum(eig_pos) * 100

  condition_levels <- levels(meta[[condition_col]])
  
  my_colours <- setNames(
    colorRampPalette(
      c("#777C3C",
        "#896C74",
        "#DDA303",
        "#419F9B",
        "#B4632D")
    )(length(condition_levels)),
    condition_levels)
  
  plot(
    bd,
    main = cohort_name,
    cex.main = 1.5,
    xlab = paste0("PCoA1 (", round(pc_var[1], 1), "%)"),
    ylab = paste0("PCoA2 (", round(pc_var[2], 1), "%)"),
    col = my_colours[levels(bd$group)],
    seg.col = "grey80",
    hull = FALSE,
    ellipse = TRUE,
    bty = "l",
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex = 1.5)
  
  r2 <- round(ad$R2[1], 3)
  p_value <- ad$`Pr(>F)`[1]
  
  text(
    x = par("usr")[2] - 0.05 * diff(par("usr")[1:2]),
    y = par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
    labels = paste0(
      "R² = ", r2,
      "\np = ", signif(p_value, 3)),
    adj = c(0, 1),
    cex = 1.3)
  
  message(" > Plots saved to: ", plot_dir)
  
  dev.off()
  
  list(
    result = result,
    warnings = warnings_list,
    error = error_msg
  )
}

#' Log PERMANOVA results and append to bray_res()
log_beta_diversity <- function(
    cohort_name,
    conditions,
    n_disease,
    result,
    bray_res) {
  
  # If no results
  if (is.null(result$result)) {
    
    row <- data.frame(
      cohort_name = cohort_name,
      conditions = conditions,
      n_disease = n_disease,
      status = "FAILED",
      r2 = NA_real_,
      p_value = NA_real_,
      dispersion_p = NA_real_,
      n = NA_integer_,
      warnings = if (!is.null(result$warnings)) {
        paste(result$warnings, collapse = " | ")
      } else NA,
      error = if (!is.null(result$error)) {
        paste(result$error, collapse = " | ")
      } else NA,
      stringsAsFactors = FALSE)
    
    bray_res[[cohort_name]] <- row
    return(bray_res)
  }
  
  # If success
  res <- result$result
  ad <- res$adonis
  
  r2 <- if (!is.null(ad$R2)) ad$R2[1] else NA_real_
  p_value <- if (!is.null(ad$`Pr(>F)`)) ad$`Pr(>F)`[1] else NA_real_
  
  bd_test <- res$dispersion_test
  
  disp_p <- NA_real_
  if (!is.null(bd_test) &&
      !is.null(bd_test$tab) &&
      "Pr(>F)" %in% colnames(bd_test$tab)) {
    disp_p <- bd_test$tab$`Pr(>F)`[1]
  }
  
  n_samples <- if (!is.null(res$meta)) nrow(res$meta) else NA_integer_

  row <- data.frame(
    cohort_name = cohort_name,
    conditions = conditions,
    n_disease = n_disease,
    status = "SUCCESS",
    r2 = r2,
    p_value = p_value,
    dispersion_p = disp_p,
    n = n_samples,
    warnings = if (!is.null(result$warnings)) {
      paste(result$warnings, collapse = " | ")
    } else NA,
    error = result$error,
    stringsAsFactors = FALSE
  )
  
  bray_res[[cohort_name]] <- row
  
  return(bray_res)
}

#' Export PERMANOVA results by saving as CSV
export_analysis <- function(bray_res, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Bind log
  permanova_df <- dplyr::bind_rows(bray_res)
  
  # Save CSV
  write_csv(permanova_df, file.path(out_dir, "bray_results.csv"))
}


# Define main function ---------------------------------------------------------
run_bray_curtis <- function(cohorts, out_dir) {
  
  bray_res <- list()
  
  for (cohort_name in names(cohorts)) {
    
    cohort <- cohorts[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check conditions
    condition_invalid <- validate_conditions(info, cohort_name)
    
    if (!is.null(condition_invalid)) {
      bray_res <- log_beta_diversity(
        cohort_name = cohort_name,
        conditions = NA_character_,
        n_disease = NA_character_,
        result = NULL,
        bray_res = bray_res)
      
      next
    }
    
    message("\nCalculating beta-diversity for cohort: ", cohort_name)
    
    # Calculate beta-diversity for filtered cohort
    result <- calculate_beta_diversity(
      cohort = info$cohort, 
      cohort_name = cohort_name,
      out_dir = out_dir)
    
    conditions = paste(unique(info$conditions), collapse = ", ")
    
    # Log results
    bray_res <- log_beta_diversity(
      cohort_name = cohort_name,
      conditions = conditions,
      n_disease = info$n_disease,
      result = result,
      bray_res = bray_res)
  }
  
  message("\nExporting analysis results")
  
  export_analysis(bray_res, out_dir)
  message("Done. Results written to: ", out_dir)
  
  return(bray_res)
}
  
    
# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary.rds") 


# Execute ----------------------------------------------------------------------
bray_res <- run_bray_curtis(primary_cohorts, out_dir)


