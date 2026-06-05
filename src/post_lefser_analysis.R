# Perform post-LEfSe analysis:
#  1. Import LEfSe results from csv
#  2. Generate and save LEfSe plots
#  3. Post-LEfSe differential analysis using ANCOM-BC


# Load packages ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(stringr)


# Define input/output directory ------------------------------------------------
in_dir <- "results/lefser_analysis"
out_dir <- "results/post_lefser_analysis"

# Load data --------------------------------------------------------------------
df <- read.csv(file.path(in_dir, "lefser_results.csv"), stringsAsFactors = FALSE)

dir.create(file.path(out_dir), recursive = TRUE, showWarnings = FALSE)


# Define helper functions ------------------------------------------------------
#' Extract species label from full clade string ---------------------------------
df <- df %>%
  mutate(
    species = case_when(
      str_detect(features, "s__") ~ sub(".*s__", "s__", features),
      str_detect(features, "g__") ~ sub(".*g__", "g__", features),
      str_detect(features, "f__") ~ sub(".*f__", "f__", features),
      TRUE ~ features
    ),
    # Remove trailing "|" 
    species = sub("\\|.*", "", species),
    
    species = sub("^[a-z]__", "", species),
    species = gsub("_", " ", species)
  )

# Define helper functions ------------------------------------------------------

log_plot <- function(comp, n_features, status, reason = NA) {
  data.frame(
    comparison = comp,
    n_features = n_features,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

smart_case <- function(x) {
  x <- gsub("_", " ", x)
  
  words <- strsplit(x, " ")[[1]]
  
  words <- ifelse(
    words == toupper(words),
    words,  # keep acronyms like IBD, CRC
    stringr::str_to_title(tolower(words))
  )
  
  paste(words, collapse = " ")
}


# Define main functions --------------------------------------------------------
save_plots <- function(df,
                       out_dir,
                       min_features = 3,
                       max_features = 30) {
  
  dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  
  comparisons <- unique(df$comparison)
  plot_log <- list()
  
  for (comp in comparisons) {
    
    df_sub <- df %>%
      dplyr::filter(comparison == comp)
    
    n_features_raw <- nrow(df_sub)
    
    cohort_name <- unique(df_sub$cohort)[1]
    
    
    # skip empty comparisons early
    if (n_features_raw == 0) {
      plot_log[[length(plot_log) + 1]] <- data.frame(
        comparison = comp,
        n_features = 0,
        status = "SKIPPED",
        reason = "No data"
      )
      next
    }
    
    disease_name <- unique(df_sub$disease)[1]
    disease_label <- smart_case(disease_name)
    healthy_label <- "Healthy"
    
    # ---------------------------------------------------------
    # SKIP CONDITION
    # ---------------------------------------------------------
    if (n_features_raw < min_features) {
      plot_log[[length(plot_log) + 1]] <- data.frame(
        comparison = comp,
        n_features = n_features_raw,
        status = "SKIPPED",
        reason = "Too few features"
      )
      next
    }
    
    # ---------------------------------------------------------
    # LIMIT FEATURES (safer version)
    # ---------------------------------------------------------
    df_sub <- df_sub %>%
      dplyr::mutate(scores = as.numeric(scores)) %>%
      dplyr::filter(!is.na(scores)) %>%
      dplyr::slice_max(order_by = abs(scores),
                       n = max_features,
                       with_ties = FALSE)
    
    n_features <- nrow(df_sub)
    
    if (n_features == 0) next
    
    # ---------------------------------------------------------
    # LABELS
    # ---------------------------------------------------------
    df_sub <- df_sub %>%
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
    
    # ensure fill mapping is correct
    fill_values <- c("steelblue", "firebrick")
    names(fill_values) <- c(healthy_label, disease_label)
    
    # ---------------------------------------------------------
    # PLOT (extra safety)
    # ---------------------------------------------------------
    p <- ggplot2::ggplot(df_sub, ggplot2::aes(
      x = scores,
      y = stats::reorder(species, scores, FUN = mean),
      fill = group_label
    )) +
      ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
      
      ggplot2::scale_fill_manual(
        values = fill_values,
        name = "Enriched in"
      ) +
      
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.05, 0.05))
      ) +
      
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title         = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
        axis.text.y        = ggplot2::element_text(size = 8, face = "italic"),
        legend.position    = "bottom",
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank()
      ) +
      
      ggplot2::labs(
        title = paste0(cohort_name, " (", healthy_label, " vs ", disease_label, ")"),
        x = "LDA score (log10)",
        y = NULL
      )
    
    # ---------------------------------------------------------
    # SAVE
    # ---------------------------------------------------------
    plot_height <- max(4, n_features * 0.35 + 1.5)
    
    safe_comp <- gsub("[/\\\\]", "_", comp) 
    out_path <- file.path(out_dir, "plots", paste0(safe_comp, ".pdf"))
    
    out_path <- file.path(
      out_dir,
      "plots",
      paste0(safe_comp, ".pdf")
    )
    
    ggplot2::ggsave(
      filename = out_path,
      plot = p,
      width = 9,
      height = plot_height,
      device = "pdf"
    )
    
    # optional: also print in interactive sessions
    # print(p)
    
    plot_log[[length(plot_log) + 1]] <- data.frame(
      comparison = comp,
      n_features = n_features_raw,
      plotted_features = n_features,
      status = "PLOTTED",
      reason = NA
    )
  }
  
  plot_log_df <- dplyr::bind_rows(plot_log)
  
  write.csv(
    plot_log_df,
    file.path(out_dir, "plots/plot_log.csv"),
    row.names = FALSE
  )
  
  invisible(plot_log_df)
}


lefse_plots <- save_plots(df, out_dir)


# Run --------------------------------------------------------------------------
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
lefse_plots <- save_plots(df, out_dir)



# TODO: statistifcally analyse features in health/disease

