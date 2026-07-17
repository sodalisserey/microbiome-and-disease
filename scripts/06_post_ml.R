# Check ML results and compile logs using functions from utils.R
# 1. Define input/output directories and functions
# 2. Execute main function run_ml_check() which runs:
#     a. check_results()
#         * checks for missing results from in_dir/comparisons.txt list
#     b. compile_logs()
#         * compiles all individual training logs from out_dir/training_log and   
#           saves as out_dir/combined_training.csv
#         * compiles all individual validation logs from out_dir/validation_log  
#           and saves as out_dir/combined_validation.csv

# Load packages and dependencies -----------------------------------------------
library(CrcBiomeScreen)
library(pROC)
library(PRROC)
library(caret)
library(boot)
library(tools)
library(dplyr)
library(tidyr)
library(readr)
library(grid)
library(ggplot2)
source("R/utils.R")

# Define input/output directories ----------------------------------------------
prep_dir <- "results/04_pre_ml"
in_dir <- "results/05_ml"
out_dir <- "results/06_post_ml"


# Define helper function -------------------------------------------------------
check_results <- function(prep_dir, in_dir) {
  
  expected <- file_path_sans_ext(
    basename(readLines(file.path(prep_dir, "comparisons.txt"))))
  
  expected <- sub("_final$", "", expected)
  
  completed <- sub("_evaluated\\.rds$", "", list.files(file.path(in_dir, "final")))
  
  missing <- setdiff(expected, completed)
  
  if (length(missing) == 0) {
    message("All comparisons evaluated successfully")
  } else {
    message(length(missing), " comparison(s) missing:")
    message("   ", missing)
    message("")
  }
  
  return(missing)
}
  
# Compile individual logs, save and return combined csv
compile_logs <- function(log_dir, out_dir, out_file) {
  
  files <- list.files(
    log_dir, 
    pattern = "\\.csv$", 
    full.names = TRUE)
  
  # Exclude output file
  files <- files[basename(files) != out_file]
  
  combined_df <- bind_rows(lapply(files, read_csv, show_col_types = FALSE))
  
  write_csv(combined_df, file.path(out_dir, out_file))
}

calculate_auc <- function(
    plot = FALSE,
    comparison_key,
    y_true,
    y_prob_rf,
    y_prob_xgb,
    out_dir) {
  
  pr_dir <- file.path(out_dir, "pr")
  roc_dir <- file.path(out_dir, "roc")
  dir.create(pr_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(roc_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate ROC-AUC
  suppressWarnings({
    suppressMessages({
      roc_rf <- roc(y_true, y_prob_rf)
      roc_auc_rf <- as.numeric(auc(roc_rf))
      ci_rf <- ci.auc(roc_rf, method = "delong")
      roc_xgb <- roc(y_true, y_prob_xgb)
      roc_auc_xgb <- as.numeric(auc(roc_xgb))
      ci_xgb <- ci.auc(roc_xgb, method = "delong")
    })})
  
  # Calculate PR-AUC
  scores_pos_rf <- y_prob_rf[y_true == "disease"]
  scores_neg_rf <- y_prob_rf[y_true != "disease"]
  scores_pos_xgb <- y_prob_xgb[y_true == "disease"]
  scores_neg_xgb <- y_prob_xgb[y_true != "disease"]
  
  pr_rf <- pr.curve(scores.class0 = scores_pos_rf, 
                    scores.class1 = scores_neg_rf, curve = TRUE)
  
  pr_auc_rf <- pr_rf$auc.integral
  
  pr_xgb <- pr.curve(scores.class0 = scores_pos_xgb, 
                     scores.class1 = scores_neg_xgb, curve = TRUE)
  
  pr_auc_xgb <- pr_xgb$auc.integral

  # Conditional plotting
  if (plot) {
    
    # For ROC curve
    roc_df <- rbind(
      data.frame(
        FPR = 1 - roc_rf$specificities,
        TPR = roc_rf$sensitivities,
        Model = "RF"),
      data.frame(
        FPR = 1 - roc_xgb$specificities,
        TPR = roc_xgb$sensitivities,
        Model = "XGB"))
    
    roc_labels <- c(
      RF = paste0(
        "RF (AUC = ", sprintf("%.3f", roc_auc_rf),
        ", 95% CI: ", sprintf("%.3f", ci_rf[1]),
        "-", sprintf("%.3f", ci_rf[3]), ")"),
      XGB = paste0(
        "XGB (AUC = ", sprintf("%.3f", roc_auc_xgb),
        ", 95% CI: ", sprintf("%.3f", ci_xgb[1]),
        "-", sprintf("%.3f", ci_xgb[3]), ")"))
    
    roc_plot <- ggplot(roc_df, aes(FPR, TPR, colour = Model)) +
      geom_line(linewidth = 1.2) +
      geom_abline(
        slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      scale_colour_manual(
        values = c("RF" = "#A7AF5A", "XGB" = "#896C74"),
        labels = roc_labels) +
      coord_equal() +
      labs(
        title = paste0(comparison_key, " ROC Curve"),
        x = "False positive rate", y = "True positive rate", colour = NULL) +
      theme_classic() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, margin = margin(r = 15)),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.direction = "vertical",
        legend.text = element_text(size = 16))
    
    # For PR curve
    pr_df <- rbind(
      data.frame(
        Recall = pr_rf$curve[,1],
        Precision = pr_rf$curve[,2],
        Model = "RF"),
      data.frame(
        Recall = pr_xgb$curve[,1],
        Precision = pr_xgb$curve[,2],
        Model = "XGB"))
    
    pr_labels <- c(
      RF = paste0(
        "RF (AUC = ", sprintf("%.3f", pr_auc_rf), ")"),
      XGB = paste0(
        "XGB (AUC = ", sprintf("%.3f", pr_auc_xgb), ")"))
    
    pr_plot <- ggplot(pr_df, aes(Recall, Precision, colour = Model)) +
      geom_line(linewidth = 1.2) +
      scale_colour_manual(
        values = c("RF" = "#A7AF5A", "XGB" = "#896C74"),
        labels = pr_labels) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      coord_equal() +
      labs(
        title = paste0(comparison_key, " PR Curve"),
        x = "Recall", y = "Precision", colour = NULL) +
      theme_classic() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, margin = margin(r = 15)),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.direction = "vertical",
        legend.text = element_text(size = 16))
    
    # Save plots as pdf
    ggsave(
      file.path(roc_dir, paste0(comparison_key, ".pdf")),
      roc_plot,
      width = 6,
      height = 6)
    
    ggsave(
      file.path(pr_dir, paste0(comparison_key, ".pdf")),
      pr_plot,
      width = 6,
      height = 6)
  }
  
  return(list(
    roc_rf = roc_rf,
    roc_auc_rf = roc_auc_rf,
    ci_rf = ci_rf,
    roc_xgb = roc_xgb,
    roc_auc_xgb = roc_auc_xgb,
    ci_xgb = ci_xgb,
    pr_rf = pr_rf,
    pr_auc_rf = pr_auc_rf,
    pr_xgb = pr_xgb,
    pr_auc_xgb = pr_auc_xgb))
}

calculate_brier <- function(y_true, y_prob, positive_class = "disease") {
  
  y_binary <- as.numeric(y_true == positive_class)
  
  mean((y_prob - y_binary)^2)
}

plot_calibration <- function(
    plot = FALSE,
    comparison_key,
    y_true, 
    y_prob_rf,
    y_prob_xgb,
    out_dir) {
  
  # Create df
  calib_data <- data.frame(
    obs = factor(
      y_true,
      levels = c("disease", "control")),
    RF = y_prob_rf,
    XGB = y_prob_xgb)

  # Calculate
  cal <- calibration(
    obs ~ RF + XGB,
    data = calib_data,
    cuts = 10)

  # Save plot as pdf
  if (plot) {
    
    pdf(file.path(out_dir, paste0(comparison_key, ".pdf")),
        width = 6,
        height = 6)
    
    print(xyplot(
      cal,
      main = comparison_key,
      col = c("#A7AF5A", "#896C74"),
      lwd = 2,
      pch = 16,
      panel = function(x, y, ...) {
        
        # 45-degree reference line
        panel.abline(a = 0, b = 1, col = "grey50", lwd = 1.5, lty = "dashed")
        
        panel.xyplot(x, y, ...)
        
        # Bottom axis line
        grid.lines(
          x = unit(c(0, 1), "npc"),
          y = unit(c(0, 0), "npc"),
          gp = gpar(col = "black", lwd = 2))
        
        # Left axis line
        grid.lines(
          x = unit(c(0, 0), "npc"),
          y = unit(c(0, 1), "npc"),
          gp = gpar(col = "black", lwd = 2))
      },
      par.settings = list(
        axis.line = list(col = "transparent"),
        superpose.line = list(col = c("#A7AF5A", "#896C74"), lwd = 2),
        superpose.symbol = list(col = c("#A7AF5A", "#896C74"), pch = 16),
        fontsize = list(
          text = 20)),
      axis=function(side,line.col,...){
        if(side%in%c("left","bottom")){
          axis.default(side=side,ticks="yes",line.col="black",...)
        }
      },
      type = c("b", "p", "l"),
      auto.key = list(lines = TRUE,
                      points = FALSE,
                      columns = 2,
                      space = "bottom")))
    
    
    dev.off()
    
  }

  return(cal)
}

plot_matrix <- function(
    plot = FALSE,
    comparison_key,
    mat_rf,
    mat_xgb,
    out_dir) {
  
  if (plot == FALSE) {
    return(NULL)
  }
  
  # Normalise by row
  norm_rf <- prop.table(mat_rf, margin = 1)
  norm_xgb <- prop.table(mat_xgb, margin = 1)
  
  # Generate data frames
  df_rf <- as.data.frame(as.table(norm_rf))
  names(df_rf) <- c("Reference", "Prediction", "Proportion")
  df_rf$count <- as.vector(mat_rf)
  df_rf$model <- "RF"
  
  df_xgb <- as.data.frame(as.table(norm_xgb))
  names(df_xgb) <- c("Reference", "Prediction", "Proportion")
  df_xgb$count <- as.vector(mat_xgb)
  df_xgb$model <- "XGB"
  
  df <- rbind(df_rf, df_xgb)
  
  plot <- ggplot(df,
                 aes(Prediction, Reference, fill = Proportion)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f%%\n(n = %d)", Proportion * 100, count))) +
    facet_wrap(~ model) +
    labs(title = comparison_key) +
    coord_equal() +
    scale_fill_gradient(low = "grey90", high = "#dfa57eff") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.x = element_text(
        margin = margin(t = 15)),
      axis.title.y = element_text(
        margin = margin(r = 15)),
      strip.text = element_text(
        size = 16))
  
  suppressMessages(ggsave(file.path(out_dir, paste0(comparison_key, ".pdf")), plot))
}

get_relationship <- function(train_disease, val_disease) {
  
  # Same disease
  if (train_disease == val_disease)
    return("same_disease")
  
  # Precursor / disease continuum
  continuum <- list(
    c("adenoma", "CRC"),
    c("IGT", "T2D"),
    c("prehypertension", "hypertension"))
  
  for (pair in continuum) {
    if (all(c(train_disease, val_disease) %in% pair))
      return("continuum")
  }
  
  # Shared biology
  if (all(c(train_disease, val_disease) %in% 
          c("ACVD", "CAD", "HF", "hypertension", "T2D")))
    return("shared_biology")
  
  if (all(c(train_disease, val_disease) %in% 
          c("IBD", "CRC")))
    return("shared_biology")
  
  if (all(c(train_disease, val_disease) %in% 
          c("T1D", "T2D")))
    return("shared_biology")

  return("unrelated")
}

extract_training <- function(mod, train_res_dir) {
  
  message(" > Extracting training results")
  
  meta <- mod$metadata
  obj <- mod$obj
  comparison_key <- meta$comparison
  
  # Initialise directories
  plots_dir <- file.path(train_res_dir, "plots")
  cal_dir <- file.path(plots_dir, "calibration")
  mat_dir <- file.path(plots_dir, "confusion")
  auc_dir <- file.path(plots_dir, "auc")
  
  lapply(
    c(plots_dir, cal_dir, mat_dir, auc_dir),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE)
    
  # Extract predictions data
  y_true = obj@ModelData$TestLabel
  y_prob_rf = obj@EvaluateResult$RF$predictions[, "disease"]
  y_prob_xgb = obj@EvaluateResult$XGBoost$predictions
  
  # Plot ROC-AUC and PR-AUC
  auc <- calculate_auc(
    comparison_key = comparison_key,
    y_true = y_true,
    y_prob_rf = y_prob_rf,
    y_prob_xgb = y_prob_xgb,
    out_dir = auc_dir)
  
  # Calculate Brier score
  brier_rf <- calculate_brier(y_true, y_prob_rf)
  brier_xgb <- calculate_brier(y_true, y_prob_xgb)
  
  # Plot calibration curve
  cal <- plot_calibration(
    comparison_key = comparison_key, 
    y_true = y_true, 
    y_prob_rf = y_prob_rf, 
    y_prob_xgb = y_prob_xgb, 
    out_dir = cal_dir)
  
  # Plot confusion matrix
  mat_rf = obj@EvaluateResult$RF$conf.matrix$table
  mat_xgb = obj@EvaluateResult$XGBoost$ConfusionMatrix$table

  plot_matrix(
    comparison_key = comparison_key,
    mat_rf = mat_rf,
    mat_xgb = mat_xgb,
    out_dir = mat_dir)
    
  # Return training results
  return(list( 
    
    meta = list(
      train_name = meta$comparison,
      train_cohort = meta$cohort,
      train_disease = meta$disease,
      val_name = meta$comparison,
      val_cohort = meta$cohort,
      val_disease = meta$disease,
      comparison_type = "same_cohort_same_disease",
      relationship = "internal"),
    
    results = list(
      rf = list(
        mod = "RF",
        roc_auc = auc$roc_auc_rf,
        ci_lower = as.numeric(auc$ci_rf[1]),
        ci_higher = as.numeric(auc$ci_rf[3]),
        pr_auc = auc$pr_auc_rf,
        brier = brier_rf,
        accuracy = as.numeric(obj@EvaluateResult$RF$conf.matrix$overall[1]),
        balanced_accuracy = as.numeric(obj@EvaluateResult$RF$BalancedAccuracy),
        classification_error = (1 - as.numeric(obj@EvaluateResult$RF$BalancedAccuracy)),
        sensitivity = as.numeric(obj@EvaluateResult$RF$conf.matrix$byClass[1]),
        specificity = as.numeric(obj@EvaluateResult$RF$conf.matrix$byClass[2]),
        precision = as.numeric(obj@EvaluateResult$RF$Precision),
        recall = as.numeric(obj@EvaluateResult$RF$Recall),
        f1 = as.numeric(obj@EvaluateResult$RF$F1)),
      
      xgb = list(
        mod = "XGB",
        roc_auc = auc$roc_auc_xgb,
        ci_lower = as.numeric(auc$ci_xgb[1]),
        ci_higher = as.numeric(auc$ci_xgb[3]),
        pr_auc = auc$pr_auc_xgb,
        brier = brier_xgb,
        accuracy = as.numeric(obj@EvaluateResult$XGBoost$ConfusionMatrix$overall[1]),
        balanced_accuracy = as.numeric(obj@EvaluateResult$XGBoost$BalancedAccuracy),
        classification_error = (1 - as.numeric(obj@EvaluateResult$XGBoost$BalancedAccuracy)),
        sensitivity = as.numeric(obj@EvaluateResult$XGBoost$ConfusionMatrix$byClass[1]),
        specificity = as.numeric(obj@EvaluateResult$XGBoost$ConfusionMatrix$byClass[2]),
        precision = as.numeric(obj@EvaluateResult$XGBoost$Precision),
        recall = as.numeric(obj@EvaluateResult$XGBoost$Recall),
        f1 = as.numeric(obj@EvaluateResult$XGBoost$F1))),
    
    data = list(
      calibration = cal,
      
      matrix = list(
        rf = mat_rf,
        xgb = mat_xgb),
      
      roc = list(
        rf = auc$roc_rf,
        xgb = auc$roc_xgb),
      
      pr = list(
        rf = auc$pr_rf,
        xgb = auc$pr_xgb),
      
      predictions = list(
        y_true = y_true,
        y_prob_rf = y_prob_rf,
        y_prob_xgb = y_prob_xgb),

      objects = list(
        rf = obj@EvaluateResult$RF,
        xgb = obj@EvaluateResult$XGBoost))))
}

extract_validation <- function(mod, train_res, val_res_dir) {
  
  message(" > Extracting validation results")
  
  meta <- mod$metadata
  obj <- mod$obj
  val <- mod$validation
  
  # Initialise directories
  plots_dir <- file.path(val_res_dir, "plots")
  cal_dir <- file.path(plots_dir, "calibration")
  mat_dir <- file.path(plots_dir, "confusion")
  auc_dir <- file.path(plots_dir, "auc")
  
  lapply(
    c(plots_dir, cal_dir, mat_dir, auc_dir),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE)
  
  val_res <- list()

  # Loop through validation results from one training cohort
  for (i in seq_along(val$val_cohort)) {
    
    val_key <- paste(
      meta$comparison, val$val_cohort[i], val$val_disease[i], sep = "_")
    
    val_rf <- obj@PredictResult$RF[[val$val_name[i]]]
    val_xgb <- obj@PredictResult$XGBoost[[val$val_name[i]]]
    
    # Extract predictions data
    y_true = val_rf$predictions$True_label
    y_prob_rf = val_rf$predictions$disease
    y_prob_xgb = val_xgb$predictions$disease
    
    # Plot ROC-AUC and PR-AUC
    auc <- calculate_auc(
      comparison_key = val_key,
      y_true = y_true,
      y_prob_rf = y_prob_rf,
      y_prob_xgb = y_prob_xgb,
      out_dir = auc_dir)

    # Calculate Brier score
    brier_rf <- calculate_brier(y_true, y_prob_rf)
    brier_xgb <- calculate_brier(y_true, y_prob_xgb)

    # Plot calibration curve
    cal <- plot_calibration(
      comparison_key = val_key,
      y_true = y_true,
      y_prob_rf = y_prob_rf,
      y_prob_xgb = y_prob_xgb,
      out_dir = cal_dir)
    
    # Plot confusion matrix
    mat_rf = val_rf$conf.matrix$table
    mat_xgb = val_xgb$conf.matrix$table
    
    plot_matrix(
      comparison_key = val_key,
      mat_rf = mat_rf,
      mat_xgb = mat_xgb,
      out_dir = mat_dir)
    
    val_res[[val_key]] <- list(
      
      meta = list(
        train_name = val$train_name[i],
        train_cohort = val$train_cohort[i],
        train_disease = val$train_disease[i],
        val_name = val$val_name[i],
        val_cohort = val$val_cohort[i],
        val_disease = val$val_disease[i],
        comparison_type = val$comparison_type[i],
        relationship = get_relationship(val$train_disease[i], val$val_disease[i])),
      
      results = list(
        rf = list(
          mod = "RF",
          roc_auc = auc$roc_auc_rf,
          ci_lower = as.numeric(auc$ci_rf[1]),
          ci_higher = as.numeric(auc$ci_rf[3]),
          pr_auc = auc$pr_auc_rf,
          brier = brier_rf,
          accuracy = as.numeric(val_rf$conf.matrix$overall[1]),
          balanced_accuracy = as.numeric(val_rf$BalancedAccuracy),
          classification_error = (1 - as.numeric(val_rf$BalancedAccuracy)),
          sensitivity = as.numeric(val_rf$conf.matrix$byClass[1]),
          specificity = as.numeric(val_rf$conf.matrix$byClass[2]),
          precision = as.numeric(val_rf$Precision),
          recall = as.numeric(val_rf$Recall),
          f1 = as.numeric(val_rf$F1)),
        
        xgb = list(
          mod = "XGB",
          roc_auc = auc$roc_auc_xgb,
          ci_lower = as.numeric(auc$ci_xgb[1]),
          ci_higher = as.numeric(auc$ci_xgb[3]),
          pr_auc = auc$pr_auc_xgb,
          brier = brier_xgb,
          accuracy = as.numeric(val_xgb$conf.matrix$overall[1]),
          balanced_accuracy = as.numeric(val_xgb$BalancedAccuracy),
          classification_error = (1 - as.numeric(val_xgb$BalancedAccuracy)),
          sensitivity = as.numeric(val_xgb$conf.matrix$byClass[1]),
          specificity = as.numeric(val_xgb$conf.matrix$byClass[2]),
          precision = as.numeric(val_xgb$Precision),
          recall = as.numeric(val_xgb$Recall),
          f1 = as.numeric(val_xgb$F1))),
    
      data = list(
        calibration = cal,
        
        matrix = list(
          rf = val_rf$conf.matrix,
          xgb = val_xgb$conf.matrix),
          
        roc = list(
          rf = val_rf$roc.curve,
          xgb = val_xgb$roc.curve),
        
        pr = list(
          rf = auc$pr_rf,
          xgb = auc$pr_xgb),
        
        predictions = list(
          y_true = y_true,
          y_prob_rf = val_rf$predictions$disease,
          y_prob_xgb = val_xgb$predictions$disease),
        
        objects = list(
          rf = val_rf,
          xgb = val_xgb)))
  }

return(val_res)

}

analyse_training <- function(train_res, train_res_dir) {
  
  train_analysis <- list()

  for (comparison_key in names(train_res)) {
    
    training <- train_res[[comparison_key]]
    meta <- training$meta
    results <- training$results
    data <- training$data
    
  # Perform DeLong's test to compare AUCs
  roc_rf <- data$roc$rf
  roc_xgb <- data$roc$xgb
  
  suppressWarnings({
    suppressMessages({
      delong_res <- roc.test(roc_rf, roc_xgb, method = "delong")
  })})
  
  analysis_res <- list(
    delta_roc_auc = results$rf$roc_auc  - results$xgb$roc_auc,
    roc_z_stat = as.numeric(delong_res$statistic),
    roc_ci_lower = as.numeric(delong_res$conf.int[1]),
    roc_ci_upper = as.numeric(delong_res$conf.int[2]),
    roc_pval = as.numeric(delong_res$p.value),
    delta_pr_auc = results$rf$pr_auc  - results$xgb$pr_auc,
    delta_brier = results$xgb$brier  - results$rf$brier)
  
  train_analysis[[comparison_key]] <- c(meta, analysis_res)
    
  }
  
  analysed_df <- do.call(rbind, lapply(train_analysis, as.data.frame))
  
  # Calculate median delta
  median_delta_roc_auc = median(analysed_df$delta_roc_auc, na.rm = TRUE)
  median_delta_pr_auc = median(analysed_df$delta_pr_auc, na.rm = TRUE)
  median_delta_brier = median(analysed_df$delta_brier, na.rm = TRUE)
  
  # Test whether median delta is different from zero
  wilcox_roc <- wilcox.test(
    analysed_df$delta_roc_auc,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE)
  
  wilcox_pr <- wilcox.test(
    analysed_df$delta_pr_auc,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE)
  
  wilcox_brier <- wilcox.test(
    analysed_df$delta_brier,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE)
  
  # FDR correction 
  p_values <- list(
    roc_auc = wilcox_roc$p.value,
    pr_auc = wilcox_pr$p.value,
    brier = wilcox_brier$p.value)
  
  p_adj <- p.adjust(p_values, method = "fdr")
  p_adj <- as.list(p_adj)
  
  # Calculate mean and CIs
  mean(analysed_df$delta_roc_auc)
  sd(analysed_df$delta_roc_auc)
  roc_ttest <- t.test(analysed_df$delta_roc_auc)$conf.int
  
  mean(analysed_df$delta_pr_auc)
  sd(analysed_df$delta_pr_auc)
  pr_ttest <- t.test(analysed_df$delta_pr_auc)$conf.int
  
  mean(analysed_df$delta_brier)
  sd(analysed_df$delta_brier)
  brier_ttest <- t.test(analysed_df$delta_brier)$conf.int
  
  # Calculate percentages of XGB wins
  pct_xgb_roc_win <- mean(analysed_df$delta_roc_auc < 0, na.rm = TRUE) * 100
  pct_xgb_pr_win <- mean(analysed_df$delta_pr_auc < 0, na.rm = TRUE) * 100
  pct_xgb_brier_win <- mean(analysed_df$delta_brier < 0, na.rm = TRUE) * 100
  
  analysed_summary <- list(
    nrow = nrow(analysed_df),
    median_delta_roc_auc = median_delta_roc_auc,
    median_delta_pr_auc = median_delta_pr_auc,
    median_delta_brier = median_delta_brier,
    pct_xgb_roc_win = pct_xgb_roc_win,
    pct_xgb_pr_win = pct_xgb_pr_win,
    pct_xgb_brier_win = pct_xgb_brier_win,
    roc_pval = p_values$roc_auc,
    roc_padj = p_adj$roc_auc,
    pr_pval = p_values$pr_auc,
    pr_padj = p_adj$pr_auc,
    brier_pval = p_values$brier,
    brier_padj = p_adj$brier,
    roc_ttest_lower = roc_ttest[1],
    roc_ttest_higher = roc_ttest[2],
    pr_ttest_lower = pr_ttest[1],
    pr_ttest_higher = pr_ttest[2],
    brier_ttest_lower = brier_ttest[1],
    brier_ttest_higher = brier_ttest[2])
  
  # Save analysis data frame as csv
  write.csv(
    analysed_df,
    file = file.path(train_res_dir, "analysed.csv"),
    row.names = FALSE)
  
  # Save analysis summary as csv
  analysed_summary_df <- as.data.frame(analysed_summary)
  
  transposed_df <- data.frame(
    colnames(analysed_summary_df),
    as.vector(t(analysed_summary_df)),
    stringsAsFactors = FALSE)
  
  write.table(
    transposed_df,
    file = file.path(train_res_dir, "analysed_summary.csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE)
  
  message(" > Training analysis saved to: ", train_res_dir)

  return(list(
    analysed = train_analysis,
    analysed_summary = analysed_summary))
}

analyse_validation <- function(val_res, val_res_dir) {
  
  val_analysis <- list()
  
  for (comparison_key in names(val_res)) {
    
    validation <- val_res[[comparison_key]]
    meta <- validation$meta
    results <- validation$results
    data <- validation$data
    
    # Perform DeLong's test to compare AUCs
    roc_rf <- data$roc$rf
    roc_xgb <- data$roc$xgb
    
    suppressWarnings({
      suppressMessages({
        delong_res <- roc.test(roc_rf, roc_xgb, method = "delong")
      })})
    
    analysis_res <- list(
      delta_roc_auc = results$rf$roc_auc  - results$xgb$roc_auc,
      roc_z_stat = as.numeric(delong_res$statistic),
      roc_ci_lower = as.numeric(delong_res$conf.int[1]),
      roc_ci_upper = as.numeric(delong_res$conf.int[2]),
      roc_pval = as.numeric(delong_res$p.value),
      delta_pr_auc = results$rf$pr_auc  - results$xgb$pr_auc,
      delta_brier = results$xgb$brier  - results$rf$brier)
    
    val_analysis[[comparison_key]] <- c(meta, analysis_res)
  }
  
  analysed_df <- do.call(rbind, lapply(val_analysis, as.data.frame))
  
  # Calculate median delta
  median_delta_roc_auc = median(analysed_df$delta_roc_auc, na.rm = TRUE)
  median_delta_pr_auc = median(analysed_df$delta_pr_auc, na.rm = TRUE)
  median_delta_brier = median(analysed_df$delta_brier, na.rm = TRUE)
  
  # Test whether median delta AUCs are significant overall
  wilcox_roc <- wilcox.test(
    analysed_df$delta_roc_auc,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE)
  
  wilcox_pr <- wilcox.test(
    analysed_df$delta_pr_auc,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE)
  
  wilcox_brier <- wilcox.test(
    analysed_df$delta_brier,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE)
  
  # FDR correction 
  p_values <- list(
    roc_auc = wilcox_roc$p.value,
    pr_auc = wilcox_pr$p.value,
    brier = wilcox_brier$p.value)
  
  p_adj <- p.adjust(p_values, method = "fdr")
  p_adj <- as.list(p_adj)
  
  # Calculate mean and CIs
  mean(analysed_df$delta_roc_auc)
  sd(analysed_df$delta_roc_auc)
  roc_ttest <- t.test(analysed_df$delta_roc_auc)$conf.int
  
  mean(analysed_df$delta_pr_auc)
  sd(analysed_df$delta_pr_auc)
  pr_ttest <- t.test(analysed_df$delta_pr_auc)$conf.int
  
  mean(analysed_df$delta_brier)
  sd(analysed_df$delta_brier)
  brier_ttest <- t.test(analysed_df$delta_brier)$conf.int
  
  # Calculate percentages of XGB wins
  pct_xgb_roc_win <- mean(analysed_df$delta_roc_auc < 0, na.rm = TRUE) * 100
  pct_xgb_pr_win <- mean(analysed_df$delta_pr_auc < 0, na.rm = TRUE) * 100
  pct_xgb_brier_win <- mean(analysed_df$delta_brier < 0, na.rm = TRUE) * 100
  
  analysed_summary <- list(
    nrow = nrow(analysed_df),
    median_delta_roc_auc = median_delta_roc_auc,
    median_delta_pr_auc = median_delta_pr_auc,
    median_delta_brier = median_delta_brier,
    pct_xgb_roc_win = pct_xgb_roc_win,
    pct_xgb_pr_win = pct_xgb_pr_win,
    pct_xgb_brier_win = pct_xgb_brier_win,
    roc_pval = p_values$roc_auc,
    roc_padj = p_adj$roc_auc,
    pr_pval = p_values$pr_auc,
    pr_padj = p_adj$pr_auc,
    brier_pval = p_values$brier,
    brier_padj = p_adj$brier,
    roc_ttest_lower = roc_ttest[1],
    roc_ttest_higher = roc_ttest[2],
    pr_ttest_lower = pr_ttest[1],
    pr_ttest_higher = pr_ttest[2],
    brier_ttest_lower = brier_ttest[1],
    brier_ttest_higher = brier_ttest[2])
  
  # Save analysis data frame as csv
  write.csv(
    analysed_df,
    file = file.path(val_res_dir, "analysed.csv"),
    row.names = FALSE)
  
  # Save analysis summary as csv
  analysed_summary_df <- as.data.frame(analysed_summary)
  
  transposed_df <- data.frame(
    colnames(analysed_summary_df),
    as.vector(t(analysed_summary_df)),
    stringsAsFactors = FALSE)
  
  write.table(
    transposed_df,
    file = file.path(val_res_dir, "analysed_summary.csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE)
  
  message(" > Validation analysis saved to: ", val_res_dir)
  
  return(list(
    analysed = val_analysis,
    analysed_summary = analysed_summary))
  
}

analyse_by_categories <- function(
    train_res, 
    train_analysis, 
    val_res, 
    val_analysis, 
    stat_res_dir) {
  
  message("\nComparing performance across categories")
  
  lapply(
    c(stat_res_dir),
    dir.create, recursive = TRUE, showWarnings = FALSE)
  
  train_analysed <- train_analysis$analysed
  val_analysed <- val_analysis$analysed
  
  # Combine train results and analysis with validation results and analysis
  combo_results <- c(train_res, val_res)
  combo_results_df <- do.call(
    rbind,
    lapply(combo_results, function(x) {
      rbind(
        data.frame(x$meta, x$results$rf, stringsAsFactors = FALSE),
        data.frame(x$meta, x$results$xgb, stringsAsFactors = FALSE))
    }))
  
  combo_analysed <- c(train_analysed, val_analysed)
  combo_analysed_df <- do.call(rbind, lapply(combo_analysed, as.data.frame))
  
  # Save data frames as csv
  write.csv(
    combo_results_df,
    file = file.path(stat_res_dir, "combo_results.csv"),
    row.names = FALSE)
  
  write.csv(
    combo_analysed_df,
    file = file.path(stat_res_dir, "combo_analysed.csv"),
    row.names = FALSE)
  
  message(" > Combined training and validation results saved to: ", stat_res_dir)
  
  # Creating data frame for comparison type analysis
  combo_wide <- combo_results_df %>%
    select(
      train_name, val_name, comparison_type, relationship, 
      mod, roc_auc, pr_auc, brier) %>%
    pivot_wider(
      names_from = mod,
      values_from = c(roc_auc, pr_auc, brier),
      names_sep = "_") %>%
    mutate(
      delta_roc_auc = roc_auc_RF - roc_auc_XGB,
      delta_pr_auc = pr_auc_RF - pr_auc_XGB,
      delta_brier = brier_RF - brier_XGB)
  
  metrics <- c("roc_auc", "pr_auc", "brier")
  
  # Test significant difference in median ΔAUCs within comparison types
  message(" > Testing difference in median ΔAUCs within comparison types")
  
  wilcox_res <- list(
    roc = list(),
    pr = list(),
    brier = list())
  
  for (type in unique(combo_wide$comparison_type)) {
    
    tmp <- subset(
      combo_wide,
      comparison_type == type)
    
    roc_res <- wilcox.test(
      delta_roc_auc ~ 1,
      data = tmp,
      mu = 0,
      conf.int = TRUE)
    
    pr_res <- wilcox.test(
      delta_pr_auc ~ 1,
      data = tmp,
      mu = 0,
      conf.int = TRUE)
    
    brier_res <- wilcox.test(
      delta_brier ~ 1,
      data = tmp,
      mu = 0,
      conf.int = TRUE)
    
    wilcox_res$roc[[type]] <- list(
      p_value = roc_res$p.value,
      ci_lower = roc_res$conf.int[1],
      ci_higher = roc_res$conf.int[2])
    
    wilcox_res$pr[[type]] <- list(
      p_value = pr_res$p.value,
      ci_lower = pr_res$conf.int[1],
      ci_higher = pr_res$conf.int[2])
    
    wilcox_res$brier[[type]] <- list(
      p_value = brier_res$p.value,
      ci_lower = brier_res$conf.int[1],
      ci_higher = brier_res$conf.int[2])
  }
  
  # FDR correction within each metric
  roc_adj <- p.adjust(
    sapply(wilcox_res$roc, `[[`, "p_value"),
    method = "fdr")
  
  pr_adj <- p.adjust(
    sapply(wilcox_res$pr, `[[`, "p_value"),
    method = "fdr")
  
  brier_adj <- p.adjust(
    sapply(wilcox_res$brier, `[[`, "p_value"),
    method = "fdr")
  
  for (type in names(wilcox_res$roc)) {
    
    wilcox_res$roc[[type]]$p_adj <- unname(roc_adj[type])
    wilcox_res$pr[[type]]$p_adj <- unname(pr_adj[type])
    wilcox_res$brier[[type]]$p_adj <- unname(brier_adj[type])
  }
  
  names(wilcox_res) <- metrics
  
  wilcox_df <- do.call(
    rbind,
    lapply(metrics, function(metric) {
      
      do.call(
        rbind,
        lapply(names(wilcox_res[[metric]]), function(type) {
          
          data.frame(
            metric = metric,
            comparison_type = type,
            p_value = wilcox_res[[metric]][[type]]$p_value,
            p_adj = wilcox_res[[metric]][[type]]$p_adj,
            ci_lower = wilcox_res[[metric]][[type]]$ci_lower,
            ci_higher = wilcox_res[[metric]][[type]]$ci_higher,
            stringsAsFactors = FALSE)
        }))
    }))
  
  fn <- "1. wilcox.csv"
  
  write.csv(
    wilcox_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  # Test significance of difference in delta AUC across comparison_types
  message(" > Testing difference in ΔAUC across comparison types")
  
  kruskal_res <- lapply(metrics, function(metric) {
    res <- kruskal.test(
      combo_wide[[paste0("delta_", metric)]] ~ combo_wide$comparison_type)
    
    list(
      statistic = unname(res$statistic),
      df = unname(res$parameter),
      p_value = res$p.value)
  })
  
  names(kruskal_res) <- metrics
  
  pairwise_res <- lapply(metrics, function(metric) {
    res <- pairwise.wilcox.test(
      combo_wide[[paste0("delta_", metric)]],
      combo_wide$comparison_type,
      p.adjust.method = "BH")
    
    res$p.value
  })
  
  names(pairwise_res) <- metrics
  
  kruskal_df <- do.call(
    rbind,
    lapply(names(kruskal_res), function(metric) {
      data.frame(
        metric = metric,
        statistic = kruskal_res[[metric]]$statistic,
        df = kruskal_res[[metric]]$df,
        p_value = kruskal_res[[metric]]$p_value,
        stringsAsFactors = FALSE)
    }))
  
  pairwise_df <- do.call(
    rbind,
    lapply(names(pairwise_res), function(metric) {
      df <- as.data.frame(as.table(pairwise_res[[metric]]))
      colnames(df) <- c(
        "comparison_type_1", "comparison_type_2", "p_adj")
      df$metric <- metric
      df[!is.na(df$p_adj), ]
    }))
  
  fn <- "2. kruskal.csv"
  
  write.csv(
    kruskal_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  fn <- "3. pairwise.csv"
  
  write.csv(
    pairwise_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  # Test correlation between RF vs XGB performance
  message(" > Testing correlations between RF and XGB performance")
  
  cor_rf_vs_xgb <- lapply(metrics, function(metric) {
    res <- cor.test(
      combo_wide[[paste0(metric, "_RF")]],
      combo_wide[[paste0(metric, "_XGB")]],
      method = "spearman",
      exact = FALSE)
    list(
      rho = unname(res$estimate),
      p_value = res$p.value,
      n = sum(complete.cases(
        combo_wide[[paste0(metric, "_RF")]],
        combo_wide[[paste0(metric, "_XGB")]])))
  })
  
  names(cor_rf_vs_xgb) <- metrics
  
  # Adjust p values
  p_adj <- p.adjust(
    sapply(cor_rf_vs_xgb, `[[`, "p_value"),
    method = "fdr")
  
  for (metric in names(cor_rf_vs_xgb)) {
    cor_rf_vs_xgb[[metric]]$p_adj <- unname(p_adj[metric])
  }
  
  names(cor_rf_vs_xgb) <- metrics
  
  cor_rf_vs_xgb_df <- do.call(
    rbind,
    lapply(metrics, function(metric) {
      data.frame(
        metric = metric,
        rho = cor_rf_vs_xgb[[metric]]$rho,
        p_value = format.pval(cor_rf_vs_xgb[[metric]]$p_value, digits = 3),
        p_adj = format.pval(cor_rf_vs_xgb[[metric]]$p_adj, digits = 3),
        n = cor_rf_vs_xgb[[metric]]$n)
    }))
  
  fn <- "4. cor_rf_vs_xgb.csv"
  
  write.csv(
    cor_rf_vs_xgb_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  # Test correlation as biological relationship increases increases
  message(" > Testing correlations between biological relationship and RF/XGB AUCs")
  
  relationship_map <- c(
    "same_disease" = 1,
    "continuum" = 2,
    "shared_biology" = 3,
    "unrelated" = 4)
  
  combo_wide$relationship <- relationship_map[combo_wide$relationship]
  
  # RF models
  cor_rel_rf <- lapply(metrics, function(metric) {
    res <- cor.test(
      combo_wide$relationship,
      combo_wide[[paste0(metric, "_RF")]],
      method = "spearman",
      exact = FALSE)
    list(
      rho = unname(res$estimate),
      p_value = res$p.value,
      n = sum(complete.cases(
        combo_wide$relationship,
        combo_wide[[paste0(metric, "_RF")]])))
  })
  
  names(cor_rel_rf) <- metrics
  
  # Adjust p values
  p_adj <- p.adjust(
    sapply(cor_rel_rf, `[[`, "p_value"),
    method = "fdr")
  
  for (metric in names(cor_rel_rf)) {
    cor_rel_rf[[metric]]$p_adj <- unname(p_adj[metric])
  }
  
  # XGB models
  cor_rel_xgb <- lapply(metrics, function(metric) {
    res <- cor.test(
      combo_wide$relationship,
      combo_wide[[paste0(metric, "_XGB")]],
      method = "spearman",
      exact = FALSE)
    list(
      rho = unname(res$estimate),
      p_value = res$p.value,
      n = sum(complete.cases(
        combo_wide$relationship,
        combo_wide[[paste0(metric, "_XGB")]])))
  })
  
  names(cor_rel_xgb) <- metrics
  
  # Adjust p values
  p_adj <- p.adjust(
    sapply(cor_rel_xgb, `[[`, "p_value"),
    method = "fdr")
  
  for (metric in names(cor_rel_xgb)) {
    cor_rel_xgb[[metric]]$p_adj <- unname(p_adj[metric])
  }
  
  cor_rel_vs_auc_df <- rbind(
    
    do.call(
      rbind,
      lapply(names(cor_rel_rf), function(metric) {
        data.frame(
          model = "RF",
          metric = metric,
          rho = cor_rel_rf[[metric]]$rho,
          p_value = format.pval(cor_rel_rf[[metric]]$p_value, digits = 3),
          p_adj = format.pval(cor_rel_rf[[metric]]$p_adj, digits = 3),
          n = cor_rel_rf[[metric]]$n)
      })
    ),
    
    do.call(
      rbind,
      lapply(names(cor_rel_xgb), function(metric) {
        data.frame(
          model = "XGB",
          metric = metric,
          rho = cor_rel_xgb[[metric]]$rho,
          p_value = format.pval(cor_rel_xgb[[metric]]$p_value, digits = 3),
          p_adj = format.pval(cor_rel_xgb[[metric]]$p_adj, digits = 3),
          n = cor_rel_xgb[[metric]]$n)
      })))
  
  fn <- "5. cor_rel_vs_auc.csv"
  
  write.csv(
    cor_rel_vs_auc_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", file.path(stat_res_dir, fn))
  
  # Test delta as biological relationship increases
  message(" > Testing correlations between biological relationship and RF/XGB ΔAUCs")
  
  cor_delta_res <- lapply(metrics, function(metric) {
    
    res <- cor.test(
      combo_wide$relationship,
      combo_wide[[paste0("delta_", metric)]],
      method = "spearman",
      exact = FALSE)
    
    list(
      rho = unname(res$estimate),
      p_value = res$p.value,
      n = sum(complete.cases(
        combo_wide$relationship,
        combo_wide[[paste0("delta_", metric)]])))
  })
  
  names(cor_delta_res) <- metrics
  
  p_adj <- p.adjust(
    sapply(cor_delta_res, `[[`, "p_value"),
    method = "fdr")
  
  for (metric in names(cor_delta_res)) {
    cor_delta_res[[metric]]$p_adj <- unname(p_adj[metric])
  }
  
  cor_rel_vs_delta_df <- do.call(
    rbind,
    lapply(names(cor_delta_res), function(metric) {
      data.frame(
        metric = metric,
        rho = cor_delta_res[[metric]]$rho,
        p_value = format.pval(cor_delta_res[[metric]]$p_value, digits = 3),
        p_adj = format.pval(cor_delta_res[[metric]]$p_adj, digits = 3),
        n = cor_delta_res[[metric]]$n)
    }))
  
  fn <- "6. cor_rel_vs_delta.csv"
  
  write.csv(
    cor_rel_vs_delta_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  # Summarise effect sizes according to comparison types
  message("\nSummarising statistics across comparison types")
  
  message(" > Calculating effects sizes")
  effect_summary <- combo_wide %>%
    group_by(comparison_type) %>%
    summarise(
      n = n(),
      median_delta_roc_auc = median(delta_roc_auc, na.rm = TRUE),
      IQR_delta_roc_auc = IQR(delta_roc_auc, na.rm = TRUE),
      median_delta_pr_auc = median(delta_pr_auc, na.rm = TRUE),
      IQR_delta_pr_auc = IQR(delta_pr_auc, na.rm = TRUE),
      median_delta_brier = median(delta_brier, na.rm = TRUE),
      IQR_delta_brier = IQR(delta_brier, na.rm = TRUE))
  
  fn <- "7. effect_summary.csv"
  
  write.csv(
    effect_summary,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  # Compare variance across comparison types
  message(" > Calculating variances")
  
  variance <- combo_results_df %>%
    group_by(mod, comparison_type) %>%
    summarise(
      sd_roc = sd(roc_auc, na.rm = TRUE),
      sd_pr = sd(pr_auc, na.rm = TRUE),
      sd_brier = sd(brier, na.rm = TRUE),
      .groups = "drop")
  
  fn <- "8. variance_summary.csv"
  
  write.csv(
    variance,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  return(list(
    data = list(
      combo_results = combo_results_df,
      combo_wide = combo_wide),
    
    results = list(
      wilcox_df = wilcox_df,
      kruskal_df = kruskal_df,
      pairwise_df = pairwise_df,
      cor_rf_vs_xgb_df = cor_rf_vs_xgb_df,
      cor_rel_vs_auc_df = cor_rel_vs_auc_df,
      cor_rel_vs_delta_df = cor_rel_vs_delta_df,
      effect_summary = effect_summary)))
}

make_boxplot_delta <- function(
    combo_wide, 
    metrics,
    plot_dir) {
  
  for (metric in metrics) {
    p <- ggplot(combo_wide,
                aes(x = comparison_type, y = .data[[metric]])) +
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    ggsave(file.path(plot_dir, paste("comp_vs_", metric, ".pdf")), 
           plot = p, width = 6, height = 4, units = "in")
  }
}

make_boxplot_ind <- function(
    combo_results, 
    metrics,
    plot_dir) {
  
  for (metric in metrics) {
    p <- ggplot(combo_results,
                aes(x = comparison_type, y = .data[[metric]], colour=mod)) +
      geom_boxplot() +
      scale_colour_manual(
        values = c(
          "RF" = "#A7AF5A",
          "XGB" = "#896C74")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    
    ggsave(file.path(plot_dir, paste("comp_vs_", metric, ".pdf")), 
           plot = p, width = 6, height = 4, units = "in")
  }
}

make_scatterplot <- function(
    combo_wide, 
    metrics,
    plot_dir) {
  
  for (metric in metrics) {
    
    rf_col <- paste0(metric, "_RF")
    xgb_col <- paste0(metric, "_XGB")
    
    p <- ggplot(combo_wide,
                aes(x = .data[[rf_col]], y = .data[[xgb_col]])) +
      geom_point() +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        colour = "grey50") +
      labs(
        x = paste("RF", metric),
        y = paste("XGB", metric),
        title = paste("RF vs XGB:", metric)) +
      theme_classic()
    
    # TODO add rho/correlation stats
    
    ggsave(
      file.path(plot_dir, paste0("RF_vs_XGB_", metric, ".pdf")),
      plot = p, width = 6, height = 5, units = "in")
  }
}

generate_plots <- function(combo_results, combo_wide, stat_res_dir) {
  
  message("\nGenerating plots")
  
  plot_dir <- file.path(stat_res_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  make_boxplot_delta(
    combo_wide = combo_wide, 
    metrics = c("delta_roc_auc", "delta_pr_auc", "delta_brier"),
    plot_dir = plot_dir)
  
  make_boxplot_ind(
    combo_results = combo_results, 
    metrics = c("roc_auc", "pr_auc", "brier"),
    plot_dir = plot_dir)
  
  make_scatterplot(
    combo_wide = combo_wide, 
    metrics = c("roc_auc", "pr_auc", "brier"),
    plot_dir = plot_dir)

  message("   - results saved to: ", plot_dir)
}

# Define main function ---------------------------------------------------------
# Check for missing results, compile all logs and save as CSV
run_analysis <- function(
    prep_dir,
    in_dir,
    out_dir) {
  
  train_res_dir <- file.path(out_dir, "train_results")
  val_res_dir <- file.path(out_dir, "val_results")
  stat_res_dir <- file.path(out_dir, "stat_results")

  lapply(
    c(out_dir, train_res_dir, val_res_dir, stat_res_dir),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE)
  
  # Confirm all comparisons have been evaluated
  missing <- check_results(prep_dir, in_dir)
  
  # Combine training and validation logs
  if (length(missing) != 0) {

    message("Evaluate all comparisons before combining logs")
    next

  } else {

  train_dir <- file.path(in_dir, "training_log")

  compile_logs(
    log_dir = train_dir,
    out_dir = out_dir,
    out_file = "training_log.csv")

  val_dir <- file.path(in_dir, "validation_log")

  compile_logs(
    log_dir = val_dir,
    out_dir = out_dir,
    out_file = "validation_log.csv")

  message("Combined training and validation logs written to: ", out_dir)
  }
  
  # Read final models
  models <-  list.files(
    file.path(in_dir, "final"),
    pattern = "\\.rds$",
    full.names = TRUE)
  
  # Initialise loop for results extraction
  train_res <- list()
  val_res <- list()
  
  for (file in models) {
    
    mod <- readRDS(file)
    
    # Extract training results
    message("\nAccessing model: ", mod$metadata$comparison)
    
    train <- extract_training(mod, train_res_dir)
    train_res[[train$meta$train_name]] <- train
    
    # Extract validation results
    val <- extract_validation(mod, train_res, val_res_dir)
    val_res <- c(val_res, val)
  }
  
  # Bind data frames and save as csv
  message("\nSaving extracted results")
  
  training_save <- do.call(
    rbind,
    lapply(train_res, function(x) {
      rbind(
        data.frame(x$meta, x$results$rf, stringsAsFactors = FALSE),
        data.frame(x$meta, x$results$xgb, stringsAsFactors = FALSE))
    }))
  
  validation_save <- do.call(
    rbind,
    lapply(val_res, function(x) {
      rbind(
        data.frame(x$meta, x$results$rf, stringsAsFactors = FALSE),
        data.frame(x$meta, x$results$xgb, stringsAsFactors = FALSE))
    }))
  
  write.csv(
    training_save,
    file = file.path(train_res_dir, "results.csv"),
    row.names = FALSE)
  
  write.csv(
    validation_save,
    file = file.path(val_res_dir, "results.csv"),
    row.names = FALSE)
  
  message(" > Training results saved as: ", file.path(train_res_dir, "results.csv"))
  message(" > Plots saved to: ", file.path(train_res_dir, "plots"))
  message(" > Validation results saved as: ", file.path(val_res_dir, "results.csv"))
  message(" > Plots saved to: ", file.path(val_res_dir, "plots"))
  
  # Analyse training and validation results
  message("\nAnalysing extracted results")
  train_analysis <- analyse_training(train_res, train_res_dir)
  val_analysis <- analyse_validation(val_res, val_res_dir)
  
  # Perform statistica tests
  stat_res <- analyse_by_categories(train_res, train_analysis, val_res, val_analysis, stat_res_dir)
  
  # Generate plots
  generate_plots(
    combo_results = stat_res$data$combo_results, 
    combo_wide = stat_res$data$combo_wide, 
    stat_res_dir = stat_res_dir)
    
  return(list(
    train_results = train_res,
    train_analysis = train_analysis,
    validation_results = val_res,
    validation_analysis = val_analysis,
    statistic_results = stat_res))
}


# Execute ----------------------------------------------------------------------
analysis_res <- run_analysis(prep_dir, in_dir, out_dir)

