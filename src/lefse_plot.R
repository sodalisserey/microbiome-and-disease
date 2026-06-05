# Import LEfSe results from csv, generate and save plots for all comparisons

# Load packages ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(stringr)


# Define input/output directory ------------------------------------------------
in_dir <- "results/lefse_analysis"
out_dir <- "results/lefse_analysis/plots"


# Load data --------------------------------------------------------------------
df <- read.csv(file.path(in_dir, "lefser_results.csv"), stringsAsFactors = FALSE)


# Define helper functions ------------------------------------------------------
#' Extract species label from full clade string
#'
#' @param df A dataframe containing a `features` column with taxonomic strings
#' @return The input dataframe with an additional `species` column containing 
#' simplified species names
extract_species <- function(df) {
  df %>%
    dplyr::mutate(
      species = dplyr::case_when(
        stringr::str_detect(features, "s__") ~ sub(".*s__", "s__", features),
        stringr::str_detect(features, "g__") ~ sub(".*g__", "g__", features),
        stringr::str_detect(features, "f__") ~ sub(".*f__", "f__", features),
        TRUE ~ features
      ),
      species = sub("\\|.*", "", species),
      species = sub("^[a-z]__", "", species),
      species = gsub("_", " ", species)
    )
}

#' Convert disease label to smart title case
#'
#' @param x Character input string, possibly containing underscores
#' 
#' @return Character string with underscores replaced by spaces and title-cased
smart_case <- function(x) {
  x <- gsub("_", " ", x)
  
  words <- strsplit(x, " ")[[1]]
  
  words <- ifelse(
    words == toupper(words),
    words,
    stringr::str_to_title(tolower(words))
  )
  
  paste(words, collapse = " ")
}

#' Prepare LEfSe plot data by filtering scores and assigning healthy/disease label
#'
#' @param df_sub A subset of the main dataframe for a specific comparison
#' @param disease_label Character label for the disease group
#' @param healthy_label Character label for the healthy group
#' @param max_features Maximum number of features to include in the plot
#' 
#' @return A dataframe ready for plotting, with `scores` and `group_label`
get_plot_data <- function(df_sub,
                          disease_label,
                          healthy_label,
                          max_features) {
  
  df_sub %>%
    dplyr::mutate(scores = as.numeric(scores)) %>%
    dplyr::filter(!is.na(scores)) %>%
    dplyr::slice_max(
      order_by = abs(scores),
      n = max_features,
      with_ties = FALSE
    ) %>%
    dplyr::mutate(
      group_label = ifelse(
        scores > 0,
        disease_label,
        healthy_label
      ),
      group_label = factor(
        group_label,
        levels = c(healthy_label, disease_label)
      )
    )
}

#' Build a LEfSe bar plot using ggplot2
#'
#' @param df_sub A dataframe prepared for plotting (output of `get_plot_data()`)
#' @param cohort_name Name of the cohort for the plot title
#' @param disease_label Character label for the disease group
#' @param healthy_label Character label for the healthy group
#' 
#' @return A ggplot object representing the LEfSe bar plot
create_single_plot <- function(df_sub,
                               cohort_name,
                               disease_label,
                               healthy_label) {
  ggplot2::ggplot(
    df_sub,
    ggplot2::aes(
      x = scores,
      y = stats::reorder(species, scores, FUN = mean),
      fill = group_label
    )
  ) +
    ggplot2::geom_col(
      width = 0.7,
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.4,
      colour = "grey40"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        setNames("steelblue", healthy_label),
        setNames("firebrick", disease_label)
      ),
      name = "Enriched in"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        size = 8,
        face = "italic"
      ),
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
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
      y = NULL
    )
}

#' Save a single plot as a PDF file
#'
#' @param plot A ggplot object to save
#' @param comparison Name of the comparison, used for the filename
#' @param out_dir Directory where the plot will be saved
#' @param n_features Number of features in the plot (controls height)
#' 
#' @return Invisible character string of the saved file path
save_single_plot <- function(plot,
                             comparison,
                             out_dir,
                             n_features) {
  
  safe_comp <- gsub("[/\\\\]", "_", comparison)
  
  out_path <- file.path(
    out_dir,
    paste0(safe_comp, ".pdf")
  )
  
  plot_height <- max(4, n_features * 0.35 + 1.5)
  
  ggplot2::ggsave(
    filename = out_path,
    plot = plot,
    width = 9,
    height = plot_height,
    device = "pdf"
  )
  
  invisible(out_path)
}

#' Create a single log entry for plot processing
#'
#' @param comp Name of the comparison
#' @param n_features Number of features in the comparison
#' @param status Status of the plot ("PLOTTED" or "SKIPPED")
#' @param reason Character or NA with optional reason for skipping
#' 
#' @return A one-row dataframe logging the plot status
create_log <- function(comp, n_features, status, reason = NA) {
  data.frame(
    comparison = comp,
    n_features = n_features,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

#' Process comparison, update plot status (skipped/plotted) and reason
#'
#' @param df_sub A subset of the main data.frame for a specific comparison
#' @param comparison Name of the comparison
#' @param out_dir Directory to save plots
#' @param min_features Minimum number of features required to plot
#' @param max_features Maximum number of features to plot
#' 
#' @return A dataframe row recording the comparison, number of features, plotted 
#' features, status and reason
update_log <- function(df_sub,
                            comparison,
                            out_dir,
                            min_features = 3,
                            max_features = 30) {
  
  n_features_raw <- nrow(df_sub)
  
  if (n_features_raw == 0) {
    return(
      create_log(
        comparison,
        0,
        "SKIPPED",
        "No data"
      )
    )
  }
  
  cohort_name <- unique(df_sub$cohort)[1]
  
  disease_label <- smart_case(
    unique(df_sub$disease)[1]
  )
  
  healthy_label <- "Healthy"
  
  if (n_features_raw < min_features) {
    return(
      create_log(
        comparison,
        n_features_raw,
        "SKIPPED",
        "Too few features"
      )
    )
  }
  
  df_plot <- get_plot_data(
    df_sub,
    disease_label,
    healthy_label,
    max_features
  )
  
  n_features <- nrow(df_plot)
  
  if (n_features == 0) {
    return(
      create_log(
        comparison,
        n_features_raw,
        "SKIPPED",
        "No valid scores"
      )
    )
  }
  
  p <- create_single_plot(
    df_plot,
    cohort_name,
    disease_label,
    healthy_label
  )
  
  save_single_plot(
    p,
    comparison,
    out_dir,
    n_features
  )
  
  data.frame(
    comparison = comparison,
    n_features = n_features_raw,
    plotted_features = n_features,
    status = "PLOTTED",
    reason = NA_character_,
    stringsAsFactors = FALSE
  )
}


# Define main function ---------------------------------------------------------
#' Generate and save LEfSe plots for all comparisons
#'
#' @param df A dataframe containing LefSe results, including columns:`comparison`
#' `cohort`, `disease`, `features`, `scores`
#' @param out_dir Directory to save plots and the plot log CSV.
#' @param min_features Minimum number of features to generate plot (default 3)
#' @param max_features Maximum number of features included in a plot (default 30)
#' 
#' @return A dataframe summarizing all comparisons, including:
#'   `comparison`: Comparison name
#'   `n_features`: Number of features in the input data
#'   `plotted_features`: Number of features actually plotted
#'   `status`: "PLOTTED" or "SKIPPED"
#'   `reason`: Reason for skipping, if applicable
create_and_save_plots <- function(df,
                              out_dir,
                              min_features = 3,
                              max_features = 30) {
  
  dir.create(file.path(out_dir), recursive = TRUE,
             showWarnings = FALSE)
  
  df <- extract_species(df)
  
  plot_log <- lapply(
    unique(df$comparison),
    function(comp) {
      
      df_sub <- dplyr::filter(
        df,
        comparison == comp
      )
      
      update_log(
        df_sub = df_sub,
        comparison = comp,
        out_dir = out_dir,
        min_features = min_features,
        max_features = max_features
      )
    }
  )
  
  plot_log_df <- dplyr::bind_rows(plot_log)
  
  write.csv(
    plot_log_df,
    file.path(out_dir, "plot_log.csv"),
    row.names = FALSE
  )
  
  invisible(plot_log_df)
}


# Run --------------------------------------------------------------------------
create_and_save_plots(df, out_dir)
