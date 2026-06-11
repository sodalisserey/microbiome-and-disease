# Import LEfSe results from csv, generate and save plots for all comparisons

# Load packages ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(stringr)


# Define input/output directory ------------------------------------------------
in_dir <- "results/lefse_analysis"
out_dir <- "results/lefse_analysis/plots"


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
process_scores <- function(df_sub,
                          disease_label,
                          healthy_label,
                          max_features) {
  
  df_sub %>%
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

#' Build a LEfSe bar plot using ggplot2
create_plot <- function(df_sub,
                               cohort_name,
                               disease_label,
                               healthy_label) {
  ggplot2::ggplot(
    df_sub,
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
        setNames("steelblue", healthy_label),
        setNames("firebrick", disease_label)),
      
      name = "Enriched in") +
    
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12,
        hjust = 0.5),
      
      axis.text.y = ggplot2::element_text(
        size = 8,
        face = "italic"),
      
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
save_plot <- function(plot,
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

#' Log plot features and status
log_plots <- function(df_sub,
                      comparison,
                      out_dir,
                      min_features = 3,
                      max_features = 30) {
  
  create_log <- function(comp, n_features, plotted_features = NA, status, reason = NA) {
    data.frame(
      comparison = comp,
      n_features = n_features,
      plotted_features = plotted_features,
      status = status,
      reason = reason,
      stringsAsFactors = FALSE)
  }
  
  n_features_raw <- nrow(df_sub)
  
  # Case 1: no data
  if (n_features_raw == 0) {
    return(create_log(comparison, 0, NA, "SKIPPED", "No data"))
  }
  
  # Extract cohort and disease info
  cohort_name <- unique(df_sub$cohort)[1]
  disease_label <- smart_case(unique(df_sub$disease)[1])
  healthy_label <- "Healthy"
  
  # Case 2: too few features
  if (n_features_raw < min_features) {
    return(create_log(comparison, n_features_raw, NA, "SKIPPED", "Too few features"))
  }
  
  # Process and filter scores for plotting
  df_plot <- process_scores(df_sub, disease_label, healthy_label, max_features)
  n_features <- nrow(df_plot)
  
  # Case 3: no valid scores after processing
  if (n_features == 0) {
    return(create_log(comparison, n_features_raw, 0, "SKIPPED", "No valid scores"))
  }
  
  # Case 4: create and save plot
  p <- create_plot(df_plot, cohort_name, disease_label, healthy_label)
  save_plot(p, comparison, out_dir, n_features)
  
  create_log(comparison, n_features_raw, n_features, "PLOTTED", NA)
}


# Define main function ---------------------------------------------------------
#' Generate, export and log all LEfSe plots
export_plots <- function(df,
                              out_dir,
                              min_features = 3,
                              max_features = 30) {
  
  dir.create(file.path(out_dir), recursive = TRUE,
             showWarnings = FALSE)
  
  df <- extract_species(df)
  
  plot_log <- lapply(
    unique(df$comparison),
    function(comp) {
      
      df_sub <- dplyr::filter(df, comparison == comp)
      
      log_plots(
        df_sub = df_sub,
        comparison = comp,
        out_dir = out_dir,
        min_features = min_features,
        max_features = max_features)
    })
  
  plot_log_df <- dplyr::bind_rows(plot_log)
  
  write.csv(
    plot_log_df,
    file.path(out_dir, "_plot_log.csv"),
    row.names = FALSE)
  message("Done. Plots written to: ", out_dir)
  
  invisible(plot_log_df)
}


# Load data --------------------------------------------------------------------
df <- read.csv(file.path(in_dir, "lefse_results.csv"), stringsAsFactors = FALSE)


# Execute ----------------------------------------------------------------------
export_plots(df, out_dir)
