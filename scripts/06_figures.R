# Load packages and dependencies -----------------------------------------------
library(grid)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(grid)
library(pheatmap)
source("R/utils.R")


# Define input/output directories ----------------------------------------------
in_dir <- "results/05_ml_analysis"
out_dir <- "results/06_figures"


# Define helper functions ------------------------------------------------------
make_barplot <- function(
    stat_results,
    out_dir) {
  
  combo_wide <- stat_results$data$combo_wide
  
  sample_df <- combo_wide %>%
    select(
      train_name,
      n_case,
      n_control
    ) %>%
    mutate(
      train_name = gsub(
        "carcinoma_surgery_history",
        "CSH",
        train_name
      ),
      label = paste0(n_case + n_control)
    ) %>%
    distinct()
  
  sample_long <- sample_df %>%
    pivot_longer(
      cols = c(n_case, n_control),
      names_to = "class",
      values_to = "n"
    ) %>%
    mutate(
      train_name = gsub(
        "carcinoma_surgery_history",
        "CSH",
        train_name
      )
    ) %>%
    mutate(
      class = recode(
        class,
        n_case = "case",
        n_control = "control"
      )
    )
  
  p_samples <- ggplot(
    sample_long,
    aes(
      x = train_name,
      y = n,
      fill = class
    )
  ) +
    geom_col(
      position = "fill",
      width = 0.9
    ) +
    geom_text(
      data = sample_df,
      aes(
        x = factor(train_name, levels = unique(train_name)),
        y = 1.02,
        label = label
      ),
      inherit.aes = FALSE,
      angle = 90,
      size = 5,
      hjust = 0
    ) +
    labs(
      y = "Proportion", x = "",
      fill = "Condition") +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = c(
        case = "#DFA57E",
        control = "#c5cb93ff"
      ),
      labels = c(
        case = "Case",
        control = "Control"
      )
    ) +
    coord_cartesian(
      ylim = c(0, 1.12),
      clip = "off"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, size = 14),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      plot.margin = margin(
        t = 10,   
        r = 10,  
        b = 0,  
        l = 40   
      ))
  
  ggsave(
    file.path(
      out_dir,
      paste0("barplot_class_balance.pdf")),
    plot = p_samples,
    width = 13,
    height = 5,
    units = "in")
  
  
  outlier_df <- combo_wide %>%
    select(
      train_name,
      n_outlier
    ) %>%
    distinct() %>%
    mutate(
      train_name = gsub(
        "carcinoma_surgery_history",
        "CSH",
        train_name
      ),
      label = n_outlier
    )
  
  p_outlier <- ggplot(
    outlier_df,
    aes(
      x = train_name,
      y = n_outlier
    )
  ) +
    geom_col(
      aes(fill = "Outliers"),
      width = 0.9
    ) +
    scale_fill_manual(
      name = NULL, 
      values = c(Outliers = "#FCCD4A")
    ) +
    geom_text(
      data = outlier_df,
      aes(
        x = factor(train_name, levels = unique(train_name)),
        y = n_outlier + max(n_outlier) * 0.03,
        label = label
      ),
      inherit.aes = FALSE,
      angle = 90,
      size = 5,
      hjust = 0
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    coord_cartesian(
      ylim = c(0, 152),
      clip = "off"
    ) +
    labs(
      y = "Number of outliers",
      x = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, size = 14),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      plot.margin = margin(
        t = 10,   
        r = 10,  
        b = 5,  
        l = 40
      )
    )
  
  ggsave(
    file.path(
      out_dir,
      "barplot_outliers.pdf"
    ),
    p_outlier,
    width = 13,
    height = 5,
    units = "in"
  )
  
  
  relationship_df <- combo_wide %>%
    dplyr::count(relationship) %>%
    mutate(
      relationship = factor(
        relationship,
        levels = c(
          "internal",
          "same",
          "continuum",
          "shared",
          "unrelated"
        ),
        labels = c(
          "Internal",
          "Same",
          "Continuum",
          "Shared",
          "Unrelated"
        )
      )
    )
  
  p_relationship <- ggplot(
    relationship_df,
    aes(
      x = relationship,
      y = n
    )
  ) +
    geom_col(
      fill = "#FCCD4A",
      width = 0.9
    ) +
    geom_text(
      aes(
        y = n + max(n) * 0.03,
        label = n
      ),
      size = 6,
      hjust = 0.5,
      vjust = -0.3
    ) +
    labs(
      x = NULL,
      y = "Biological Relationship\n(Counts)"
    ) +
    coord_cartesian(
      ylim = c(0, 826),
      clip = "off"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 20, angle = 50, hjust = 1),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 22),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      plot.margin = margin(
        t = 30,   
        r = 10,  
        b = 5,  
        l = 40
      )
    )
  
  ggsave(
    file.path(out_dir, "barplot_relationship_counts.pdf"),
    p_relationship,
    width = 8.3,
    height = 4
  )
  
  comparison_df <- combo_wide %>%
    dplyr::count(comparison_type) %>%
    mutate(
      comparison_type = factor(
        comparison_type,
        levels = c(
          "same_cohort_same_disease",
          "cross_cohort_same_disease",
          "same_cohort_diff_disease",
          "cross_cohort_diff_disease"
        ),
        labels = c(
          "Same\ncohort\nSame disease",
          "Different\ncohort\nSame disease",
          "Same\ncohort\nDifferent disease",
          "Different\ncohort\nDifferent disease"
        )
      )
    )
  
  p_comparison <-   p_relationship <- ggplot(
    comparison_df,
    aes(
      x = comparison_type,
      y = n
    )
  ) +
    geom_col(
      fill = "#FCCD4A",
      width = 0.9
    ) +
    geom_text(
      aes(
        y = n + max(n) * 0.03,
        label = n
      ),
      size = 6,
      hjust = 0.5,
      vjust = -0.3
    ) +
    labs(
      x = NULL,
      y = "Comparison Type\n(Counts)"
    ) +
    coord_cartesian(
      ylim = c(0, 1000),
      clip = "off"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 22),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      plot.margin = margin(
        t = 30,   
        r = 10,  
        b = 5,  
        l = 40
      )
    )
  
  ggsave(
    file.path(out_dir, "barplot_comparison_counts.pdf"),
    p_comparison,
    width = 7,
    height = 4
  )
  
}

make_scatterplot <- function(
    stat_results,
    metrics,
    out_dir,
    internal) {
  
  combo_wide <- stat_results$data$combo_wide
  cor_res <- stat_results$results$cor_rf_vs_xgb
  
  if (internal) {
    combo_wide <- subset(combo_wide, relationship == "internal")
  }
  
  metric_colours <- c(
    roc_auc = "#419F9B",
    pr_auc  = "#DDA303",
    brier   = "#B4632D")
  
  metric_labels <- c(
    roc_auc = "ROC-AUC",
    pr_auc = "PR-AUC",
    brier = "Brier Score")
  
  for (metric in metrics) {
    
    metric_label <- metric_labels[[metric]]
    
    rf_col <- paste0(metric, "_RF")
    xgb_col <- paste0(metric, "_XGB")
    
    rho <- cor_res$rho[cor_res$metric == metric]
    p_adj <- cor_res$p_adj[cor_res$metric == metric]
    
    p_adj <- ifelse(
      grepl("<", p_adj),
      0,
      as.numeric(p_adj)
    )
    
    p_label <- case_when(
      p_adj < 0.001 ~ "p-adj < 0.001",
      TRUE ~ paste0("p-adj = ", sprintf("%.3f", p_adj))
    )
    
    annotation <- paste0(
      "\nrho = ", sprintf("%.3f", rho),
      "\n", p_label
    )
    
    suffix <- if (internal) "_internal" else ""
    
    p <- ggplot(combo_wide, aes(x = .data[[rf_col]], y = .data[[xgb_col]])) +
      geom_abline(
        intercept = 0,
        colour = "black") + 
      geom_point(
        shape = 16,
        size = 2, 
        colour = metric_colours[[metric]],
        alpha = 0.8) +
      labs(
        x = paste("RF", metric_label),
        y = paste("XGB", metric_label),
        title = paste0(metric, suffix)) +
      theme_classic() +
      scale_x_continuous(limits = c(0,1)) +
      scale_y_continuous(limits = c(0,1)) +
      coord_fixed(ratio = 1) +
      theme(    
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5)) +
      annotate(
        "text",
        x = Inf,
        y = -Inf,
        label = annotation,
        hjust = 1,
        vjust = -0.2,
        size = 4.3)
    
    ggsave(
      file.path(
        out_dir,
        paste0("RF_vs_XGB_", metric, suffix, ".pdf")),
      plot = p,
      width = 3,
      height = 3,
      units = "in")
  }
}

plot_delta_distribution <- function(stat_results, metrics, out_dir) {
  
  combo_wide <- stat_results$data$combo_wide
  train_summary <- stat_results$results$train_summary
  
  metric_labels <- c(
    delta_roc_auc = "Delta ROC-AUC",
    delta_pr_auc = "Delta PR-AUC",
    delta_brier = "Delta Brier Score"
  )
  
  metric_colours <- c(
    delta_roc_auc = "#419F9B",
    delta_pr_auc = "#DDA303",
    delta_brier = "#B4632D"
  )
  
  for (metric in metrics) {
    
    res_name <- metric
    if (metric == "delta_roc_auc") res_name <- "pct_xgb_roc_win"
    if (metric == "delta_pr_auc")  res_name <- "pct_xgb_pr_win"
    if (metric == "delta_brier")  res_name <- "pct_xgb_brier_win"
    
    xgb_win <- train_summary[[res_name]]
    label <- paste0("XGB wins ", sprintf("%.0f%%", xgb_win))
    
    delta <- combo_wide[[paste0(metric)]]
    df <- data.frame(delta = delta)
    
    p <- ggplot(df, aes(x = "", y = delta)) +
      geom_violin(
        fill = NA,
        colour = metric_colours[[metric]],
        linewidth = 1
      ) +
      geom_boxplot(
        width = 0.15,
        outlier.shape = NA,
        fill = "white",
        colour = metric_colours[[metric]],
        linewidth = 0.8
      ) +
      scale_y_continuous(
        limits = c(-0.4,0.4)) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      coord_flip() +
      labs(
        x = NULL,
        y = metric_labels[metric]) +
      theme_classic() +
      theme(    
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14)) +
      annotate(
        "text",
        x = 0.425,
        y = 0.4,
        label = label,
        hjust = 1,
        vjust = -0.3,
        size = 4.3)
    
    ggsave(
      file.path(out_dir, paste0("delta_distribution_", metric, ".pdf")),
      p,
      width = 3,
      height = 3
    )
  }
}

plots <- generate_plots(analysis_res, out_dir)


make_boxplot_ind <- function(
    stat_results,
    metrics,
    out_dir) {
  
  combo_results <- stat_results$data$combo_results
  combo_wide <- stat_results$data$combo_wide
  dunn <- stat_results$results$dunn
  
  metric_labels <- c(
    roc_auc = "ROC-AUC",
    pr_auc = "PR-AUC",
    brier = "Brier Score")
  
  y_limits <- list(
    roc_auc = c(0.4, 1.6),
    pr_auc = c(0, 1.85),
    brier = c(0, 1.35))

  for (metric in metrics) {
    
    plot_data <- combo_results %>%
      filter(!is.na(.data[[metric]]))
    
    plot_data <- combo_results %>%
      filter(!is.na(.data[[metric]])) %>%
      mutate(
        comparison_label = case_when(
          comparison_type == "same_cohort_same_disease" ~ 
            "Internal",
          comparison_type == "cross_cohort_same_disease" ~ 
            "Cohort",
          comparison_type == "same_cohort_diff_disease" ~ 
            "Condition",
          comparison_type == "cross_cohort_diff_disease" ~ 
            "Both"),
        comparison_label = factor(
          comparison_label,
          levels = c(
            "Internal",
            "Cohort",
            "Condition",
            "Both")),
        box_colour = if_else(
          comparison_type == "same_cohort_same_disease",
          "Internal",
          "Other"))
    
    sig_data <- dunn %>%
      filter(
        metric == !!metric,
        P.adj < 0.05) %>%
      separate(
        Comparison,
        into = c("group1", "group2"),
        sep = " - ")
    
    for (model_name in unique(plot_data$mod)) {
      
      model_data <- plot_data %>%
        filter(mod == model_name)
      
      model_sig <- sig_data %>%
        filter(model == model_name)
      
      if (nrow(model_sig) > 0) {
        
        ymax <- max(
          model_data[[metric]],
          na.rm = TRUE
        )
        
        model_sig <- model_sig %>%
          mutate(
            y.position = ymax *
              seq(
                1.05,
                1.05 + 0.1 * (nrow(model_sig) - 1),
                length.out = nrow(model_sig)))
      }
      
      p <- ggplot(
        model_data,
        aes(
          x = comparison_label,
          y = .data[[metric]])) +
        geom_boxplot(
          aes(colour = box_colour),
          size = 0.8,
          alpha = 0.8,
          outlier.size = 1.5
        ) +
        scale_colour_manual(
          values = c(
            Internal = "grey50",
            Other = if (model_name == "RF") "#A7AF5A" else "#896C74"
          ),
          guide = "none"
        ) +
        theme_classic() +
        labs(
          title = paste0(model_name, " ", metric_labels[[metric]]),
          x = "Validation Change/s",
          y = metric_labels[[metric]]) +
        scale_y_continuous(
          limits = y_limits[[metric]],
          breaks = seq(0, 1, by = 0.25)
        ) +
        theme(
          axis.title = element_text(size = 11),
          axis.text = element_text(
            hjust = 0.5,
            size = 11),
          axis.text.x = element_text(size = 11, angle = 50, hjust = 1),
          plot.title = element_text(
            size = 14,
            face = "bold",
            hjust = 0.5),
          plot.margin = margin(
            t = 10,
            r = 10,
            b = 0,
            l = 10)) 
      
      model_sig <- model_sig %>%
        mutate(
          group1 = recode(
            group1,
            "same_cohort_same_disease" = "Internal",
            "cross_cohort_same_disease" = "Cohort",
            "same_cohort_diff_disease" = "Condition",
            "cross_cohort_diff_disease" = "Both"
          ),
          group2 = recode(
            group2,
            "same_cohort_same_disease" = "Internal",
            "cross_cohort_same_disease" = "Cohort",
            "same_cohort_diff_disease" = "Condition",
            "cross_cohort_diff_disease" = "Both"
          )
        )
      
      if (nrow(model_sig) > 0) {
        
        model_sig <- model_sig %>%
          mutate(
            y.position = ymax *
              seq(
                1.15,
                1.15 + 0.2 * (nrow(model_sig) - 1),
                length.out = nrow(model_sig)
              ),
            p_label = case_when(
              P.adj < 0.001 ~ "p < 0.001",
              TRUE ~ paste0("p = ", sprintf("%.3f", P.adj))
            )
          )
        
        p <- p +
          stat_pvalue_manual(
            model_sig,
            label = "p_label",
            tip.length = 0.01,
            size = 3,
            vjust = -0.3,
            bracket.size = 0.5) +
          coord_cartesian(clip = "off")
      }
      
      ggsave(
        file.path(
          out_dir,
          paste0("comp_vs_", metric, "_", model_name, ".pdf")),
        plot = p,
        width = 2.5,
        height = 3.25,
        units = "in")
    }
  }
}

make_boxplot_delta <- function(
    stat_results,
    metrics,
    out_dir) {
  
  combo_results <- stat_results$data$combo_results
  combo_wide <- stat_results$data$combo_wide
  wilcox <- stat_results$results$wilcox
  pairwise <- stat_results$results$pairwise
  
  delta_labels <- c(
    delta_roc_auc = "Delta ROC-AUC",
    delta_pr_auc = "Delta PR-AUC",
    delta_brier = "Delta Brier Score")
  
  metric_colours <- c(
    delta_roc_auc = "#419F9B",
    delta_pr_auc  = "#DDA303",
    delta_brier   = "#B4632D")
  
  combo_wide <- combo_wide %>%
    mutate(
      comparison_label = case_when(
        comparison_type == "same_cohort_same_disease" ~ 
          "Internal",
        comparison_type == "cross_cohort_same_disease" ~ 
          "Cohort",
        comparison_type == "same_cohort_diff_disease" ~ 
          "Condition",
        comparison_type == "cross_cohort_diff_disease" ~ 
          "Both"),
      comparison_label = factor(
        comparison_label,
        levels = c(
          "Internal",
          "Cohort",
          "Condition",
          "Both")),
      box_colour = if_else(
        comparison_type == "same_cohort_same_disease",
        "Internal",
        "Other"))
  
  
  for (metric in metrics) {

    sig_metric <- sub("^delta_", "", metric)
    
    sig_wilcox <- wilcox %>%
      filter(
        .data$metric == .env$sig_metric,
        .data$p_adj < 0.05
      ) %>%
      mutate(
        comparison_label = recode(
          comparison_type,
          "same_cohort_same_disease" = "Internal",
          "cross_cohort_same_disease" = "Cohort",
          "same_cohort_diff_disease" = "Condition",
          "cross_cohort_diff_disease" = "Both"
        ),
        label = case_when(
          p_adj < 0.001 ~ "***",
          p_adj < 0.01 ~ "**",
          TRUE ~ "*"))
    
    sig_pairwise <- pairwise %>%
      dplyr::rename(
        group1 = comparison_type_1,
        group2 = comparison_type_2) %>%
      dplyr::filter(
        metric == .env$sig_metric,
        p_adj < 0.05
      ) %>%
      dplyr::mutate(
        group1 = recode(
          group1,
          "same_cohort_same_disease" = "Internal",
          "cross_cohort_same_disease" = "Cohort",
          "same_cohort_diff_disease" = "Condition",
          "cross_cohort_diff_disease" = "Both"
        ),
        group2 = recode(
          group2,
          "same_cohort_same_disease" = "Internal",
          "cross_cohort_same_disease" = "Cohort",
          "same_cohort_diff_disease" = "Condition",
          "cross_cohort_diff_disease" = "Both"
        )
      )
    
    p <- ggplot(
      combo_wide,
      aes(x = comparison_label, y = .data[[metric]])) +
      geom_hline(
        yintercept = 0,
        colour = "black") +
      geom_boxplot(
        aes(colour = box_colour),
        size = 0.8,
        alpha = 0.8,
        outlier.size = 1.5
      ) +
      scale_colour_manual(
        values = c(
          Internal = "grey50",
          Other = metric_colours[[metric]]
        ),
        guide = "none"
      ) +
      theme_classic() +
      labs(
        title = paste0(delta_labels[[metric]]),
        x = "Validation Change/s",
        y = delta_labels[[metric]]) +
      scale_y_continuous(
        breaks = seq(0, 1, by = 0.25)
      ) +
      theme(
        axis.title = element_text(size = 11),
        axis.text = element_text(
          hjust = 0.5,
          size = 11),
        axis.text.x = element_text(size = 11, angle = 50, hjust = 1),
        plot.title = element_text(
          size = 14,
          face = "bold",
          hjust = 0.5),
        plot.margin = margin(
          t = 10,
          r = 10,
          b = 0,
          l = 10)) 
    
    ymax <- max(combo_wide[[metric]], na.rm = TRUE)
    
    if (nrow(sig_wilcox) > 0) {
      
        sig_wilcox <- sig_wilcox %>%
          mutate(
            y.position = ymax * seq(
              1.15,
              1.15 + 0.1 * (n() - 1),
              length.out = n()
            )
          )
      
      p <- p +
        geom_text(
          data = sig_wilcox,
          aes(
            x = comparison_label,
            y = y.position,
            label = label
          ),
          inherit.aes = FALSE,
          size = 5)
    }
    
    if (nrow(sig_pairwise) > 0) {

      sig_pairwise <- sig_pairwise %>%
        mutate(
          y.position = ymax *      
            seq(
            2.25,
            1.15 + 0.2 * (nrow(sig_pairwise) - 1),
            length.out = nrow(sig_pairwise)
          ),
          label = case_when(
            p_adj < 0.001 ~ "p < 0.001",
            TRUE ~ paste0("p = ", sprintf("%.3f", p_adj))))

      p <- p +
        stat_pvalue_manual(
          sig_pairwise,
          label = "label",
          tip.length = 0.01,
          size = 3,
          vjust = -0.3,
          bracket.size = 0.5) +
        coord_cartesian(clip = "off")
      
    }
    
    ggsave(
      file.path(
        out_dir,
        paste("comp_vs_", metric, ".pdf")),
      plot = p,
      width = 2.5,
      height = 3.25,
      units = "in")
  }
}

generate_plots(analysis_res, out_dir)


make_rel_boxplot_abs <- function(
    stat_results,
    metrics,
    out_dir) {
  
  combo_wide <- stat_results$data$combo_wide
  rel_dunn <- stat_results$results$rel_dunn
  
  metric_labels <- c(
    roc_auc = "ROC-AUC",
    pr_auc = "PR-AUC",
    brier = "Brier Score"
  )
  
  y_limits <- list(
    roc_auc = c(0.4, 1.37),
    pr_auc = c(0, 1.4),
    brier = c(0, 0.9))
  
  models <- c("RF", "XGB")
  
  combo_wide <- combo_wide %>%
    mutate(
      relationship_label = factor(
        relationship,
        levels = c(
          "internal",
          "same",
          "continuum",
          "shared",
          "unrelated"),
        labels = c(
          "Internal",
          "Condition",
          "Continuum",
          "Biology",
          "None")))
  
  for (model in models) {
  
    for (metric in metrics) {
      
      metric_col <- paste0( metric, "_", model)
      
      plot_data <- combo_wide %>%
        select(
          relationship_label,
          value = all_of(metric_col)) %>%
        filter(
          !is.na(value)) %>%
        dplyr::mutate(
          box_colour = if_else(
          relationship_label == "Internal",
          "Internal",
          "Other"
        ))
      
      ns_data <- rel_dunn %>%
        filter(
          model == !!model,
          metric == !!metric,
          P.adj >= 0.05
        ) %>%
        tidyr::separate(
          Comparison,
          into = c("group1", "group2"),
          sep = " - "
        ) %>%
        mutate(
          group1 = recode(
            group1,
            internal = "Internal",
            same = "Condition",
            continuum = "Continuum",
            shared = "Biology",
            unrelated = "None"
          ),
          group2 = recode(
            group2,
            internal = "Internal",
            same = "Condition",
            continuum = "Continuum",
            shared = "Biology",
            unrelated = "None"
          )
        )
      
      p <- ggplot(
        plot_data,
        aes(
          x = relationship_label,
          y = value)) +
        # geom_boxplot(
        #   size = 0.8,
        #   colour = ifelse(
        #     model == "RF",
        #     "#A7AF5A", "#896C74"),
        #   alpha = 0.8,
        #   outlier.size = 1.5) +
        geom_boxplot(
          aes(colour = box_colour),
          size = 0.8,
          alpha = 0.8,
          outlier.size = 1.5
        ) +
        scale_colour_manual(
          values = c(
            Internal = "grey50",
            Other = if (model == "RF") "#A7AF5A" else "#896C74"
          ),
          guide = "none"
        ) +
        labs(
          title = paste(
            model,
            metric_labels[[metric]]),
          x = "Shared Relationship",
          y = metric_labels[[metric]]) +
        theme_classic() +
        scale_y_continuous(
          limits = y_limits[[metric]],
          breaks = seq(0, 1, by = 0.25)
        ) +
        theme(
          axis.title = element_text(size = 11),
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 11),
          axis.text.y = element_text(
            size = 11),
          plot.title = element_text(
            size = 14,
            face = "bold",
            hjust = 0.5),
          plot.margin = margin(
            t = 10,
            r = 10,
            l = 10)) 
      
      if (nrow(ns_data) > 0) {
        
        ymax <- max(plot_data$value, na.rm = TRUE)
        
        ns_data <- ns_data %>%
          mutate(
            y.position = ymax *
              seq(
                1.05,
                1.05 + 0.15 * (n() - 1),
                length.out = n()
              ),
            p_label = "ns"
          )
        
        p <- p +
          ggpubr::stat_pvalue_manual(
            ns_data,
            label = "p_label",
            tip.length = 0.01,
            size = 3,
            bracket.size = 0.5,
            vjust = -0.3)
      }

      ggsave(
        file.path(
          out_dir,
          paste0("rel_vs_", metric, "_", model_name, ".pdf")),
        plot = p,
        width = 2.5,
        height = 3.25,
        units = "in")
      
      
    }
  }
}

generate_plots(analysis_res, out_dir)

# make_rel_boxplot_delta <- function(
#     stat_results,
#     metrics,
#     out_dir) {
#   
#   combo_wide <- stat_results$data$combo_wide
#   
#   delta_labels <- c(
#     delta_roc_auc = "Delta ROC-AUC",
#     delta_pr_auc = "Delta PR-AUC",
#     delta_brier = "Delta Brier Score")
#   
#   metric_colours <- c(
#     delta_roc_auc = "#419F9B",
#     delta_pr_auc  = "#DDA303",
#     delta_brier   = "#B4632D")
#   
#   combo_wide <- combo_wide %>%
#     mutate(
#       relationship_label = factor(
#         relationship,
#         levels = c(
#           "internal",
#           "same",
#           "continuum",
#           "shared",
#           "unrelated"),
#         labels = c(
#           "Internal",
#           "Condition",
#           "Continuum",
#           "Biology",
#           "None")),
#       box_colour = if_else(
#         comparison_type == "same_cohort_same_disease",
#         "Internal",
#         "Other"))
#   
#   for (metric in metrics) {
#     
#     # sig_metric <- sub("^delta_", "", metric)
#     # 
#     # sig_wilcox <- wilcox %>%
#     #   filter(
#     #     .data$metric == .env$sig_metric,
#     #     .data$p_adj < 0.05
#     #   ) %>%
#     #   mutate(
#     #     comparison_label = recode(
#     #       comparison_type,
#     #       "same_cohort_same_disease" = "Internal",
#     #       "cross_cohort_same_disease" = "Cohort",
#     #       "same_cohort_diff_disease" = "Condition",
#     #       "cross_cohort_diff_disease" = "Both"
#     #     ),
#     #     label = case_when(
#     #       p_adj < 0.001 ~ "***",
#     #       p_adj < 0.01 ~ "**",
#     #       TRUE ~ "*"))
#     # 
#     # sig_pairwise <- pairwise %>%
#     #   dplyr::rename(
#     #     group1 = comparison_type_1,
#     #     group2 = comparison_type_2) %>%
#     #   dplyr::filter(
#     #     metric == .env$sig_metric,
#     #     p_adj < 0.05
#     #   ) %>%
#     #   dplyr::mutate(
#     #     group1 = recode(
#     #       group1,
#     #       "same_cohort_same_disease" = "Internal",
#     #       "cross_cohort_same_disease" = "Cohort",
#     #       "same_cohort_diff_disease" = "Condition",
#     #       "cross_cohort_diff_disease" = "Both"
#     #     ),
#     #     group2 = recode(
#     #       group2,
#     #       "same_cohort_same_disease" = "Internal",
#     #       "cross_cohort_same_disease" = "Cohort",
#     #       "same_cohort_diff_disease" = "Condition",
#     #       "cross_cohort_diff_disease" = "Both"
#     #     )
#     #   )
#     
#     delta_col <- paste0("delta_", metric)
#     
#     p <- ggplot(
#       combo_wide,
#       aes(
#         x = relationship_label,
#         y = .data[[delta_col]])) +
#       # geom_boxplot(
#       #   size = 0.8,
#       #   alpha = 0.8,
#       #   outlier.size = 1.5,
#       #   colour = "#6EC4C0") +
#       geom_boxplot(
#         aes(colour = box_colour),
#         size = 0.8,
#         alpha = 0.8,
#         outlier.size = 1.5
#       ) +
#       scale_colour_manual(
#         values = c(
#           Internal = "grey50",
#           Other = metric_colours[[metric]]
#         ),
#         guide = "none"
#       ) +
#       geom_hline(
#         yintercept = 0,
#         linetype = "dashed",
#         colour = "grey50") +
#       labs(
#         x = "Shared Relationship",
#         y = paste0(
#           "Delta ",
#           metric_labels[[metric]],
#           " (RF - XGB)")) +
#       theme_classic() +
#       theme(
#         axis.text.x = element_text(
#           angle = 45,
#           hjust = 0.5))
#     
#     ggsave(
#       file.path(
#         out_dir,
#         paste("rel_vs_", metric, ".pdf")),
#       p,
#       width = 4,
#       height = 3.5,
#       units = "in")
#   }
# }

make_full_heatmap <- function(
    stat_results,
    metrics,
    out_dir) {
  
  combo_wide <- stat_results$data$combo_wide
  
  metric_limits <- lapply(metrics, function(metric) {
    
    vals <- c(
      combo_wide[[paste0(metric, "_RF")]],
      combo_wide[[paste0(metric, "_XGB")]])
    
    range(vals, na.rm = TRUE)
    
  })
  
  names(metric_limits) <- metrics
  
  models <- c("RF", "XGB")
  
  for (model in models) {
    
    for (metric in metrics) {
      
      value_col <- paste0(metric, "_", model)
      
      limits <- metric_limits[[metric]]
      midpoint <- mean(limits)
      breaks <- seq(
        limits[1],
        limits[2],
        length.out = 101
      )
    
      if (metric %in% c("roc_auc", "pr_auc")) {
        
        fill_scale <- colorRampPalette(
          rev(RColorBrewer::brewer.pal(11, "Spectral"))
        )(100)

      } else if (metric == "brier") {
        
        fill_scale <- colorRampPalette(
          RColorBrewer::brewer.pal(11, "Spectral")
        )(100)
        
      }
      
      heat_matrix <- combo_wide %>%
        dplyr::select(
          train_name,
          val_name,
          value = all_of(value_col)
        ) %>%
        tidyr::pivot_wider(
          names_from = val_name,
          values_from = value
        ) %>%
        tibble::column_to_rownames("train_name") %>%
        as.matrix()
      
      rownames(heat_matrix) <- gsub(
        "carcinoma_surgery_history",
        "CSH",
        rownames(heat_matrix)
      )
      
      colnames(heat_matrix) <- gsub(
        "carcinoma_surgery_history",
        "CSH",
        colnames(heat_matrix)
      )
      
      internal_mask <- outer(
        rownames(heat_matrix),
        colnames(heat_matrix),
        FUN = "=="
      )
      
      heat_matrix_plot <- heat_matrix
      heat_matrix_plot[internal_mask] <- NA
      
      cluster_rows <- hclust(
        dist(heat_matrix),
        method = "ward.D2"
      )
      
      cluster_cols <- hclust(
        dist(t(heat_matrix)),
        method = "ward.D2"
      )
      
      p <- pheatmap(
        heat_matrix_plot,
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        color = fill_scale,
        breaks = breaks,
        na_col = "white",
        border_color = NA,
        fontsize_number = 12,
        fontsize = 12,
        fontsize_row = 10,
        fontsize_col = 10,
        angle_col = 90,
        cellwidth = 10,
        cellheight = 10,
        legend = TRUE,
        silent = TRUE,
        main = paste(
          model,
          toupper(metric)
        ))
      
        pdf(
          file.path(out_dir, paste0("clustered_heatmap_", metric, "_", model, ".pdf")),
          width = 9,
          height = 8
        )
        
        grid.draw(p$gtable)
        
        # Move legend to bottom
        p$gtable <- gtable::gtable_add_rows(
          p$gtable,
          heights = grid::unit(1, "cm"),
          pos = nrow(p$gtable)
        )
        
        p$gtable <- gtable::gtable_add_grob(
          p$gtable,
          p$gtable$grobs[[which(p$gtable$layout$name == "legend")]],
          t = nrow(p$gtable),
          l = 2,
          r = ncol(p$gtable) - 1,
          name = "legend_bottom"
        )
    
        grid.text(
          "Validation Combination",
          x = 0.5,
          y = 0.02,
          gp = gpar(fontsize = 12)
        )
          
        grid.text(
          "Training Combination",
          x = 0.02,
          y = 0.5,
          rot = 90,
          gp = gpar(fontsize = 12)
        )
        
        dev.off()
    }
  }
} 

make_condition_heatmap <- function(
    stat_results,
    metrics,
    out_dir) {
  
  combo_results <- stat_results$data$combo_results
  
  models <- c("RF", "XGB")
  
  metric_limits <- setNames(
    lapply(metrics, function(metric) {
      range(combo_results[[metric]], na.rm = TRUE)
    }),
    metrics
  )
  
  names(metric_limits) <- metrics
  
  
  for (model in models) {
    
    for (metric in metrics) {
      
      if (metric %in% c("roc_auc", "pr_auc")) {
        
        fill_scale <- colorRampPalette(
          rev(RColorBrewer::brewer.pal(11, "Spectral"))
        )(100)
        
      } else if (metric == "brier") {
        
        fill_scale <- colorRampPalette(
          RColorBrewer::brewer.pal(11, "Spectral")
        )(100)
        
      }
      
      limits <- metric_limits[[metric]]
      midpoint <- mean(limits)
      breaks <- seq(
        limits[1],
        limits[2],
        length.out = 101
      )

      heatmap_data <- combo_results %>%
        filter(relationship != "internal") %>%
        filter(mod == model) %>%
        group_by(
          train_disease,
          val_disease
        ) %>%
        summarise(
          median_auc = median(.data[[metric]], na.rm = TRUE),
          .groups = "drop")

      heatmap_mat <- heatmap_data %>%
        pivot_wider(
          names_from = val_disease,
          values_from = median_auc
        ) %>%
        tibble::column_to_rownames("train_disease") %>%
        as.matrix()

      rownames(heatmap_mat) <- gsub(
        "carcinoma_surgery_history",
        "CSH",
        rownames(heatmap_mat)
      )

      colnames(heatmap_mat) <- gsub(
        "carcinoma_surgery_history",
        "CSH",
        colnames(heatmap_mat)
      )
      
      rownames(heatmap_mat) <- gsub(
        "MECFS",
        "ME/CFS",
        rownames(heatmap_mat)
      )
      
      colnames(heatmap_mat) <- gsub(
        "MECFS",
        "ME/CFS",
        colnames(heatmap_mat)
      )
      
      rownames(heatmap_mat) <- gsub(
        "prehypertension",
        "pre-hypertension",
        rownames(heatmap_mat)
      )
      
      colnames(heatmap_mat) <- gsub(
        "prehypertension",
        "pre-hypertension",
        colnames(heatmap_mat)
      )

      p2 <- pheatmap(
        heatmap_mat,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        clustering_method = "ward.D2",
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        color = fill_scale,
        breaks = breaks,
        na_col = "white",
        border_color = NA,
        fontsize = 10,
        fontsize_row = 10,
        fontsize_col = 10,
        angle_col = 90,
        cellwidth = 10,
        cellheight = 10,
        main = paste(
          model,
          toupper(metric)
        ))

      pdf(
        file.path(out_dir, paste0("heatmap_condition_", metric, "_", model, ".pdf")),
        width = 6,
        height = 5
      )

      grid.draw(p2$gtable)

      grid.text(
        "Validation Condition",
        x = 0.5,
        y = 0.02,
        gp = gpar(fontsize = 11)
      )

      grid.text(
        "Training Condition",
        x = 0.02,
        y = 0.5,
        rot = 90,
        gp = gpar(fontsize = 11)
      )


      dev.off()
    }
  }
}

make_condition_table <- function(
    stat_results,
    metrics,
    out_dir) {
  
  combo_results <- stat_results$data$combo_results
  
  models <- c("RF", "XGB")
  
  train_condition <- list()
  target_condition <- list()
  
  for (model in models) {
    
    train_condition[[model]] <- list()
    target_condition[[model]] <- list()
    
    model_results <- combo_results %>%
      filter(mod == model)
    
    for (metric in metrics) {      
      
      internal_metric <- model_results %>%
        filter(relationship == "internal") %>%
        group_by(train_disease) %>%
        summarise(
          internal_value = median(
            .data[[metric]],
            na.rm = TRUE
          ),
          .groups = "drop"
        )
      
      # Training disease degradation
      train_metric <- model_results %>%
        filter(relationship != "internal") %>%
        left_join(
          internal_metric,
          by = "train_disease"
        ) %>%
        group_by(train_disease) %>%
        summarise(
          validation_value = median(
            .data[[metric]],
            na.rm = TRUE
          ),
          q1_validation = quantile(
            .data[[metric]],
            0.25,
            na.rm = TRUE
          ),
          q3_validation = quantile(
            .data[[metric]],
            0.75,
            na.rm = TRUE
          ),
          degradation = median(
            internal_value - .data[[metric]],
            na.rm = TRUE
          ),
          q1_degradation = quantile(
            internal_value - .data[[metric]],
            0.25,
            na.rm = TRUE
          ),
          q3_degradation = quantile(
            internal_value - .data[[metric]],
            0.75,
            na.rm = TRUE
          ),
          n_transfers = n()/32,
          .groups = "drop"
        )
      
      # Target disease performance
      target_metric <- model_results %>%
        filter(relationship != "internal") %>%
        group_by(val_disease) %>%
        summarise(
          validation_value = median(
            .data[[metric]],
            na.rm = TRUE
          ),
          q1_validation = quantile(
            .data[[metric]],
            0.25,
            na.rm = TRUE
          ),
          q3_validation = quantile(
            .data[[metric]],
            0.75,
            na.rm = TRUE
          ),
          n_transfers = n()/32,
          .groups = "drop"
        )
    
      train_condition[[model]][[metric]] <- train_metric %>%
        filter(n_transfers >= 2)
      
      target_condition[[model]][[metric]] <- target_metric %>%
        filter(n_transfers >= 2)
    }
  }

    return(list(
      train_condition = train_condition,
      target_condition = target_condition
  ))
}

plots <- make_condition_table(
    stat_results = stat_res,
    metrics = c("roc_auc", "pr_auc", "brier"),
    out_dir = out_dir)

# Define main function ---------------------------------------------------------
generate_plots <- function(analysis_res, out_dir) {
  
  message("Generating plots:")
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  stat_res <- analysis_res$statistic_results
  
  message("  > QC barplots (class balance and outliers")
  make_barplot(
    stat_results = stat_res, 
    out_dir = out_dir)
  
  message("  > RF-XGB correlation (internal-only)")
  make_scatterplot(
    stat_results = stat_res,
    metrics = c("roc_auc", "pr_auc", "brier"),
    out_dir = out_dir,
    internal = TRUE)
  
  message("  > RF-XGB correlation (all comparison types)")
  make_scatterplot(
    stat_results = stat_res,
    metrics = c("roc_auc", "pr_auc", "brier"),
    out_dir = out_dir,
    internal = FALSE)
  
  message("  > Distribution of delta metrics")
  plot_delta_distribution(
    stat_results = stat_res,
    metrics = c("delta_roc_auc", "delta_pr_auc", "delta_brier"),
    out_dir = out_dir)
  
  message("  > RF/XGB boxplots across comparison types")
  make_boxplot_ind( 
    stat_results = stat_res,
    metrics = c("roc_auc", "pr_auc", "brier"), 
    out_dir = out_dir)
  
  message("  > ΔRF-XGB boxplots across comparison types")
  make_boxplot_delta(
    stat_results = stat_res,
    metrics = c("delta_roc_auc", "delta_pr_auc", "delta_brier"),
    out_dir = out_dir)
  
  message("  > Trend line across biological relationship (absolute)")
  make_rel_boxplot_abs(
    stat_results = stat_res,
    metrics = c("roc_auc", "pr_auc", "brier"),
    out_dir = out_dir)
  
  # message("  > Trend line across biological relationship (delta)")
  # make_rel_boxplot_delta(
  #   stat_results = stat_res,
  #   metrics = c("roc_auc", "pr_auc", "brier"),
  #   out_dir = out_dir)
  
  message("  > Heatmaps of all validation results")
  make_full_heatmap(
    stat_results = stat_res,
    metrics = c("roc_auc", "pr_auc", "brier"),
    out_dir = out_dir)
  
  message("  > Heatmaps for conditions")
  make_condition_heatmap(
    stat_results = stat_res,
    metrics = c("roc_auc", "pr_auc", "brier"),
    out_dir)
  
  message("  > Barplots for conditions")
  d <- make_condition_table(
    stat_results = stat_res,
    metrics = c("roc_auc", "pr_auc", "brier"),
    out_dir)
  
  message("\nDone. Results saved to: ", out_dir)
  
  return(d)
  
}


# Load data --------------------------------------------------------------------
analysis_res <- readRDS(file.path(in_dir, "analysis_res.rds"))


# Execute ----------------------------------------------------------------------
plots <- generate_plots(analysis_res, out_dir)

