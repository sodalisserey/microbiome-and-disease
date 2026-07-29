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
library(dplyr)
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

  # Extract eigenvalues
  eig <- bd$eig
  eig_pos <- eig[eig > 0]
  pc_var <- eig_pos / sum(eig_pos) * 100
  
  # Plot
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
  
  # Sample coordinates
  scores_sites <- as.data.frame(scores(bd, display = "sites"))
  scores_sites$group <- bd$group
  
  # Group centroids
  scores_centroids <- as.data.frame(scores(bd, display = "centroids"))
  scores_centroids$group <- rownames(scores_centroids)
  
  colnames(scores_centroids)[1:2] <- c("Centroid1", "Centroid2")
  
  # Join centroid coordinates to each sample
  plot_df <- scores_sites %>%
    left_join(
      scores_centroids,
      by = "group"
    ) %>%
    mutate(
      group = recode(
        group,
        TKI_dependent_diarrhoea = "TKI_diarrhoea"
      )
      # For report only
      , group = ifelse(group == "control", "Control", "Case")
    )
  
  # For report only
  my_colours <- c(
    Control = "#777C3C",
    Case = "#B4632D")

  p <- ggplot(plot_df, aes(PCoA1, PCoA2, colour = group)) +
    
    # segments to centroids
    geom_segment(
      aes(
        xend = Centroid1,
        yend = Centroid2
      ),
      colour = "grey80",
      linewidth = 0.4
    ) +
    
    # confidence ellipses
    stat_ellipse(
      aes(fill = group),
      geom = "polygon",
      alpha = 0.15,
      colour = NA,
      show.legend = FALSE
    ) +

    scale_colour_manual(values = my_colours) +
    scale_fill_manual(values = my_colours)+
    
    # samples
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
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    )
  
  ggsave(
    file.path(
      plot_dir,
      paste0(cohort_name, "_pcoa.pdf")),
    plot = p,
    width = 3,
    height = 4,
    units = "in")
  
  message(" > Plots saved to: ", plot_dir)
  
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
  
  n_status_thresh = 0
  
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
      n_status = NA_character_,
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
  
  meta <- res$meta
  n_samples <- if (!is.null(meta)) nrow(meta) else NA_integer_
  
  n_status <- NA_character_
  
  if (!is.null(meta) && "study_condition" %in% colnames(meta)) {
    
    counts <- table(meta[["study_condition"]])
    
    n_status <- character()
    
    if ("control" %in% names(counts) &&
        counts[["control"]] < n_status_thresh) {
      n_status <- c(
        n_status,
        paste0("n control < ", n_status_thresh)
      )
    }
    
    disease_groups <- setdiff(names(counts), "control")
    
    for (d in disease_groups) {
      if (counts[[d]] < n_status_thresh) {
        n_status <- c(
          n_status,
          paste0("n ", d, " < ", n_status_thresh)
        )
      }
    }
    
    if (length(n_status) == 0) {
      n_status <- NA_character_
    } else {
      n_status <- paste(n_status, collapse = " | ")
    }
  }
  
  row <- data.frame(
    cohort_name = cohort_name,
    conditions = conditions,
    n_disease = n_disease,
    status = "SUCCESS",
    r2 = r2,
    p_value = p_value,
    dispersion_p = disp_p,
    n = n_samples,
    n_status = n_status,
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
primary_cohorts <- primary_cohorts[c("GuptaA_2019", "KieserS_2018", "IjazUZ_2017", "IaniroG_2022")]

# Execute ----------------------------------------------------------------------
bray_res <- run_bray_curtis(primary_cohorts, out_dir)


