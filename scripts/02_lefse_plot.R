# Import LEfSe results from csv, generate and save plots
# 1. Define input/output directories and functions
# 2. Load LEfSe data from in_dir/lefse_results.csv
# 3. Execute main function export_plots():
#     * Extract and format species column via extract_species()
#     * Create and log plot via create_plot() and log_plot() 
#     * Save plots as out_dir/[comparison].pdf
#     * Save plot_log into out_dir/_plot_log.csv



# Load packages ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(stringr)


# Define input/output directory ------------------------------------------------
in_dir <- "results/01_lefse_analyse"
out_dir <- "results/02_lefse_plot"


# Define helper functions ------------------------------------------------------
#' Extract species label from full clade string
extract_species <- function(df) {
  df %>%
    dplyr::mutate(
      species = dplyr::case_when(
        stringr::str_detect(features, "s__") ~ sub(".*s__", "s__", features),
        stringr::str_detect(features, "g__") ~ sub(".*g__", "g__", features),
        stringr::str_detect(features, "f__") ~ sub(".*f__", "f__", features),
        TRUE ~ features),
      
      species = sub("\\|.*", "", species),
      species = sub("^[a-z]__", "", species),
      species = gsub("_", " ", species),
      species = sub("^species:", "", species))
}

#' Convert disease label string with underscores to title case
smart_case <- function(x) {
  x <- gsub("_", " ", x)
  
  words <- strsplit(x, " ")[[1]]
  
  words <- ifelse(
    words == toupper(words),
    words,
    stringr::str_to_title(tolower(words)))
  
  paste(words, collapse = " ")
}

#' Prepare LEfSe plot data by filtering scores and assigning healthy/disease label
process_scores <- function(
    lefse_res_df,
    disease_label,
    healthy_label,
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
        healthy_label),
      
      group_label = factor(
        group_label,
        levels = c(healthy_label, disease_label)))
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
  
  disease_label <- smart_case(unique(lefse_res_df$disease)[1])
  healthy_label <- "Healthy"
  
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
    healthy_label,
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

#' Build a LEfSe bar plot using ggplot2
create_plot <- function(
    lefse_res_df,
    cohort_name,
    disease_label,
    healthy_label) {
  
  ggplot2::ggplot(
    lefse_res_df,
    ggplot2::aes(
      x = scores,
      y = stats::reorder(species, scores, FUN = mean),
      fill = group_label)) +
    
    ggplot2::geom_col(
      width = 0.7,
      colour = "white",
      linewidth = 0.2) +
    
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      colour = "grey40") +
    
    ggplot2::scale_fill_manual(
      values = c(
        setNames("#c5cb93ff", healthy_label),
        setNames("#dfa57eff", disease_label)),
      
      name = "Enriched in") +
    
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        hjust = 0.5),
      
      axis.text.y = ggplot2::element_text(
        size = 14,
        face = "italic"),
      
      # TODO legend text size increase
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()) +
    
    ggplot2::labs(
      title = paste0(
        cohort_name,
        " (",
        healthy_label,
        " vs ",
        disease_label,
        ")"
      ),
      x = "LDA score (log10)",
      y = NULL)
}

#' Save a single plot as a pdf
save_plot <- function(
    plot,
    comparison,
    out_dir,
    n_features) {
  
  safe_comp <- gsub("[/\\\\]", "", comparison)
  
  out_path <- file.path(
    out_dir,
    paste0(safe_comp, ".pdf"))
  
  plot_height <- max(4, n_features * 0.35 + 1.5)
  
  ggplot2::ggsave(
    filename = out_path,
    plot = plot,
    width = 9,
    height = plot_height,
    device = "pdf")
  
  invisible(out_path)
}


# Define main function ---------------------------------------------------------
#' Generate, export and log all LEfSe plots
export_plots <- function(
    lefse_res,
    out_dir,
    min_features = 2,
    max_features = 10) {
  
  dir.create(file.path(out_dir), recursive = TRUE,
             showWarnings = FALSE)
  
  lefse_res <- extract_species(lefse_res)
  
  plot_log <- lapply(
    unique(lefse_res$comparison),
    function(comp) {
      
      lefse_res_df <- dplyr::filter(
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
        disease_label <- smart_case(unique(lefse_res_df$disease)[1])
        
        p <- create_plot(
          res$df_plot,
          cohort_name,
          disease_label,
          "Healthy"
        )
        
        save_plot(
          p,
          comp,
          out_dir,
          nrow(res$df_plot)
        )
        
        res$log$status <- "PLOTTED"
      }
      
      res$log
    }
  )
  
  plot_log_df <- dplyr::bind_rows(plot_log)
  
  write.csv(
    plot_log_df,
    file.path(out_dir, "_plot_log.csv"),
    row.names = FALSE)
  message("Done. Plots written to: ", out_dir)
  
  invisible(plot_log_df)
}


# Load data --------------------------------------------------------------------
lefse_res <- read.csv(file.path(in_dir, "lefse_results.csv"), stringsAsFactors = FALSE)


# Execute ----------------------------------------------------------------------
export_plots(lefse_res, out_dir)
