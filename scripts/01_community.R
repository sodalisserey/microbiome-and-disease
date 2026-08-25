# Analyse differences in microbial community using `vegan` package

# 1. Define output directory and functions
# 2. Load data using utils.R
# 3. Execute main function run_analysis() which:
#   * Runs calculate_dissimilarity(), which:
#         - quantifies Bray-Curtis dissimilarity
#         - runs PERMANOVA and PERMDISP 
#         - logs results
#   * Runs plot_pcoa()
#   * Exports summary statistics as out_dir/results.csv
#   * Exports plots as out_dir/plots[cohort_name].pdf


# Load packages and dependencies -----------------------------------------------
library(vegan)
library(ggplot2)
library(dplyr)
library(readr)
library(SummarizedExperiment)
source("R/utils.R")

# Define output directory ------------------------------------------------------
out_dir <- "results/01_community"


# Define helper functions ------------------------------------------------------
# Calculate beta-diversity, run PERMANOVA and PERMDISP for each cohort
calculate_dissimilarity <- function(
    cohort,
    cohort_name,
    condition_col = "study_condition",
    permutations = 9999,
    dispersion_permutations = 999,
    n_status_thresh = 5) {
  
  # Initialise
  warnings_list <- character()
  error_msg <- NA
  
  # Analysis
  result <- tryCatch({
    
    ## extract abundance
    abund <- as.matrix(assay(cohort))
    
    ## extract metadata
    meta <- as.data.frame(colData(cohort))
    
    ## ensure sample IDs are rownames
    if (is.null(rownames(meta))) {
      rownames(meta) <- colnames(cohort)
    }
    
    ## align samples 
    common_samples <- intersect(colnames(abund), rownames(meta))
    
    if (length(common_samples) < 3) {
      stop("Too few overlapping samples")
    }
    
    abund <- abund[, common_samples, drop = FALSE]
    meta  <- meta[common_samples, , drop = FALSE]
    
    ## validate condition column
    if (!condition_col %in% colnames(meta)) {
      stop(
        "Condition column not found: ",
        condition_col)
    }
    
    ## ensure factor
    meta[[condition_col]] <- factor(
      meta[[condition_col]])
    
    message(" > Computing Bray-Curtis dissimilarity")
    X <- t(abund)
    bray <- vegdist(X, method = "bray")
    
    message(" > Running PERMANOVA and PERMDISP")
    meta[[condition_col]] <- as.factor(meta[[condition_col]])
    
    ad <- adonis2(
      bray ~ meta[[condition_col]],
      permutations = permutations)
    
    bd <- betadisper(bray, meta[[condition_col]])
    bd_test <- permutest(bd, permutations = dispersion_permutations)
    
    list(
      bray = bray,
      adonis = ad,
      betadisper = bd,
      dispersion_test = bd_test,
      metadata = meta)

  }, warning = function(w) {
    warnings_list <<- c(warnings_list, w$message)
    invokeRestart("muffleWarning")
    
  }, error = function(e) {
    error_msg <<- e$message
    return(NULL)
  })
  
  # If no results
  if (is.null(result)) {
    
    qc <- data.frame(
      cohort_name = cohort_name,
      study_conditions = NA_character_,
      n_conditions = NA_integer_,
      r2 = NA_real_,
      p_value = NA_real_,
      dispersion_p = NA_real_,
      n = NA_integer_,
      n_status = NA_character_,
      warnings = if (length(warnings_list) > 0) {
        paste(warnings_list, collapse = " | ") } else {NA_character_},
      error = error_msg,
      status = "no results",
      stringsAsFactors = FALSE)
    
    return(list(
      dissimilarity = NULL,
      metadata = NULL,
      qc = qc))
  }

  # Extract summary statistics
  ad <- result$adonis
  bd_test <- result$dispersion_test
  meta <- result$metadata
  
  study_conditions <- paste(
    unique(meta[[condition_col]]),
    collapse = ", ")
  
  n_conditions <- nlevels(
    meta[[condition_col]])
  
  r2 <- if (!is.null(ad$R2)) {
    ad$R2[1]
  } else {
    NA_real_
  }
  
  p_value <- if (!is.null(ad$`Pr(>F)`)) {
    ad$`Pr(>F)`[1]
  } else {
    NA_real_
  }
  
  dispersion_p <- NA_real_
  
  if (
    !is.null(bd_test) &&
    !is.null(bd_test$tab) &&
    "Pr(>F)" %in% colnames(bd_test$tab)
  ) {
    dispersion_p <- bd_test$tab$`Pr(>F)`[1]
  }
  
  n_samples <- nrow(meta)
  
  counts <- table(
    meta[[condition_col]])
  
  n_status <- character()
  
  if (length(counts) > 0) {
    
    for (condition in names(counts)) {
      
      if (counts[[condition]] < n_status_thresh) {
        n_status <- c(
          n_status,
          paste0("n ", condition, " < ", n_status_thresh))
      }
    }
  }
  
  if (length(n_status) == 0) {
    n_status <- NA_character_
  } else {
    n_status <- paste(
      n_status,
      collapse = " | ")
  }
  
  # Log successful results
  qc <- data.frame(
    cohort_name = cohort_name,
    study_conditions = study_conditions,
    n_conditions = n_conditions,
    r2 = r2,
    p_value = p_value,
    dispersion_p = dispersion_p,
    n = n_samples,
    n_status = n_status,
    warnings = if (length(warnings_list) > 0) {
      paste(warnings_list, collapse = " | ")} else {NA_character_},
    error = NA_character_,
    status = "SUCCESS",
    stringsAsFactors = FALSE)
  
  return(list(
    dissimilarity = list(
      bray = result$bray,
      adonis = result$adonis,
      betadisper = result$betadisper,
      dispersion_test = result$dispersion_test),
    metadata = result$metadata,
    qc = qc))
}


# Plot PCoA figures for each cohort
plot_pcoa <- function(
    cohort_name,
    processed,
    condition_col = "study_condition",
    out_dir){
  
  # Initialise
  plot_dir <- file.path(out_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  bd <- processed$dissimilarity$betadisper
  meta <- processed$metadata
  
  message(" > Plotting PCoA figure")
  
  # Calculate eigenvalues
  eig <- bd$eig
  eig_pos <- eig[eig > 0]
  
  if (length(eig_pos) < 2) {
    stop("Fewer than two positive eigenvalues available for PCoA plot")
  }
  
  pc_var <- eig_pos / sum(eig_pos) * 100
  
  # Sample coordinates
  scores_sites <- as.data.frame(scores(bd, display = "sites"))
  scores_sites$group <- bd$group
  
  # Centroid coordinates
  scores_centroids <- as.data.frame(scores(bd, display = "centroids"))
  scores_centroids$group <- rownames(scores_centroids)
  
  colnames(scores_centroids)[1:2] <- c("Centroid1", "Centroid2")
  
  # Join sample and centroid coordinates
  plot_df <- scores_sites %>%
    left_join(
      scores_centroids,
      by = "group")
  
  # Set colours
  condition_levels <- levels(meta[[condition_col]])
  
  my_colours <- setNames(
    colorRampPalette(
      c("#777C3C", "#896C74", "#DDA303", "#419F9B", "#B4632D"))
    (length(condition_levels)), condition_levels)

  # Plot figures
  p <- ggplot(plot_df, aes(PCoA1, PCoA2, colour = group)) +
    
    geom_segment(
      aes(
        xend = Centroid1,
        yend = Centroid2),
      colour = "grey80",
      linewidth = 0.4) +
    
    stat_ellipse(
      aes(fill = group),
      geom = "polygon",
      alpha = 0.15,
      colour = NA,
      show.legend = FALSE) +

    scale_colour_manual(values = my_colours) +
    scale_fill_manual(values = my_colours) +
    scale_x_continuous(breaks = scales::breaks_width(0.5)) +
    scale_y_continuous(breaks = scales::breaks_width(0.5)) +
    
    geom_point(size = 1.5) +
    labs(
      title = cohort_name,
      x = paste0("PCoA1 (", round(pc_var[1], 1), "%)"),
      y = paste0("PCoA2 (", round(pc_var[2], 1), "%)"), 
      colour = "Condition",
      fill = "Condition") +
    
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.justification = "left",
      legend.direction = "vertical",
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15)) +
    guides(
      colour = guide_legend(ncol = 1))
  
  n_conditions <- length(condition_levels)
  
  plot_height <- max(4, 2.5 + n_conditions * 0.2)
  
  ggsave(
    file.path(
      plot_dir,
      paste0(cohort_name, ".pdf")),
    plot = p,
    width = 3,
    height = plot_height,
    units = "in")
  
  invisible(p)
}


# Define main function ---------------------------------------------------------
run_analysis <- function(cohorts, out_dir) {
  
  dir.create("results", recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list()
  
  for (cohort_name in names(cohorts)) {
    
    cohort <- cohorts[[cohort_name]]
    info <- process_conditions(cohort)
    
    # Check cohort contains both case and control samples
    validation <- validate_conditions(info, cohort_name)
    if (!is.null(validation)) {next}
    
    message("\nCalculating beta-diversity for cohort: ", cohort_name)
    
    # Calculate dissimilarity
    processed <- calculate_dissimilarity(
      cohort = info$cohort, 
      cohort_name = cohort_name)
    
    results[[cohort_name]] <- processed
    
    # Plot figures
    if (!is.null(processed$dissimilarity)) {
      
      plot_pcoa(
        cohort_name = cohort_name,
        processed = processed,
        out_dir = out_dir)
    }
  }
  
  qc <- dplyr::bind_rows(lapply(results, `[[`, "qc"))
  
  write_csv(qc, file.path(out_dir, "results.csv"))
  
  message(
    "\nResults written to: ", out_dir)
  
  return(qc)
}
  
    
# Load data --------------------------------------------------------------------
primary_cohorts <- readRDS("data/primary.rds") 


# Execute ----------------------------------------------------------------------
analysis_res <- run_analysis(primary_cohorts, out_dir)


