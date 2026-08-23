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
library(FSA)
library(dplyr)
library(tidyr)
library(readr)
library(grid)
library(ggplot2)
library(ggpubr)
source("R/utils.R")

# Define input/output directories ----------------------------------------------
prep_dir <- "results/03_pre_ml"
in_dir <- "results/04_ml"
out_dir <- "results/05_ml_analysis"


# Define helper functions ------------------------------------------------------
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
    
    roc_labels <- paste0(
        "RF = ", sprintf("%.2f", roc_auc_rf),
        "\nXGB = ", sprintf("%.2f", roc_auc_xgb))
    
    roc_plot <- ggplot(roc_df, aes(FPR, TPR, colour = Model)) +
      geom_line(linewidth = 1.2) +
      geom_abline(
        slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      scale_colour_manual(
        values = c("RF" = "#A7AF5A", "XGB" = "#896C74")) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      coord_equal() +
      labs(
        title = paste0(comparison_key, " ROC Curve"),
        x = "False Positive Rate", y = "True Positive Rate", colour = NULL) +
      theme_classic() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18, margin = margin(r = 40, t = 5)),
        axis.title.y = element_text(size = 18, margin = margin(r = 5)),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        plot.margin = margin(
          r = 10,
          l = 10)) +
      annotate(
        "text",
        x = Inf,
        y = -Inf,
        label = roc_labels,
        hjust = 1,
        vjust = -0.3,
        size = 5.5)
    
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
        "RF (AUC = ", sprintf("%.2f", pr_auc_rf), ")"),
      XGB = paste0(
        "XGB (AUC = ", sprintf("%.2f", pr_auc_xgb), ")"))
    
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
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18, margin = margin(r = 30)),
        axis.title.y = element_text(size = 18),
        plot.title = element_blank(),
        legend.position = "top", 
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        legend.box.just = "left", 
        legend.text = element_text(size = 14))
    
    # Save plots as pdf
    ggsave(
      file.path(roc_dir, paste0(comparison_key, ".pdf")),
      roc_plot,
      width = 4,
      height = 4.5,
      units = "in")
    
    ggsave(
      file.path(pr_dir, paste0(comparison_key, ".pdf")),
      pr_plot,
      width = 4,
      height = 4.5,
      units = "in")
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

calculate_brier <- function(y_true, y_prob) {
  
  y_binary <- as.numeric(y_true == "disease")
  
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
    geom_text(aes(label = sprintf("%.1f%%\n(n = %d)", Proportion * 100, count)),
              size = 4) +
    facet_wrap(~ model) +
    labs(title = comparison_key) +
    coord_equal() +
    scale_fill_gradient(low = "grey90", high = "#dfa57eff",
                        guide = guide_colourbar(direction = "horizontal")) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.title = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.text.y = element_text(angle = 90, hjust = 0.5),
      strip.text = element_text(size = 14))
  
  suppressMessages(
    ggsave(file.path(out_dir, paste0(comparison_key, ".pdf")), 
           plot,
           width = 4,
           height = 2.5))
}

get_relationship <- function(train_disease, val_disease) {
  
  # Same disease
  if (train_disease == val_disease)
    return("same")
  
  # Precursor / disease continuum
  continuum <- list(
    c("adenoma", "CRC"),
    c("IGT", "T2D"),
    c("prehypertension", "hypertension"),
    c("ACVD", "CAD"))
  
  for (pair in continuum) {
    if (all(c(train_disease, val_disease) %in% pair))
      return("continuum")
  }
  
  # Shared biology
  if (all(c(train_disease, val_disease) %in% 
          c("IBD", "CRC")))
    return("shared")
  
  if (all(c(train_disease, val_disease) %in% 
          c("T1D", "T2D")))
    return("shared")
  
  if (all(c(train_disease, val_disease) %in% 
          c("hypertension", "HF")))
    return("shared")
  
  if (all(c(train_disease, val_disease) %in% 
          c("hypertension", "HF", "ACVD", "CAD")))
    return("shared")
  
  if (all(c(train_disease, val_disease) %in% 
          c("IGT", "T2D", "ACVD", "CAD")))
    return("shared")

  return("unrelated")
}

extract_training <- function(mod, train_res_dir, plot) {
  
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
  
  # Extract sample info
  train_labels <- obj@ModelData[["TrainLabel"]]
  
  sample_counts <- as.data.frame(table(train_labels)) %>%
    dplyr::rename(
      class = train_labels,
      n = Freq)
  
  n_case <- sample_counts$n[sample_counts$class == "disease"]
  n_control <- sample_counts$n[sample_counts$class == "control"]
  n_outlier <- length(obj@OutlierSamples)
    
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
    out_dir = auc_dir,
    plot = TRUE)
  
  # Calculate Brier score
  brier_rf <- calculate_brier(y_true, y_prob_rf)
  brier_xgb <- calculate_brier(y_true, y_prob_xgb)
  
  # Plot calibration curve
  cal <- plot_calibration(
    comparison_key = comparison_key, 
    y_true = y_true, 
    y_prob_rf = y_prob_rf, 
    y_prob_xgb = y_prob_xgb, 
    out_dir = cal_dir,
    plot = plot)
  
  # Plot confusion matrix
  mat_rf = obj@EvaluateResult$RF$conf.matrix$table
  mat_xgb = obj@EvaluateResult$XGBoost$ConfusionMatrix$table

  plot_matrix(
    comparison_key = comparison_key,
    mat_rf = mat_rf,
    mat_xgb = mat_xgb,
    out_dir = mat_dir,
    plot = plot)
    
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
      relationship = "internal",
      n_case = n_case,
      n_control = n_control,
      n_outlier = n_outlier),
    
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

extract_validation <- function(mod, train_res, val_res_dir, plot) {
  
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
    
    # Extract sample info
    train_labels <- obj@ModelData[["TrainLabel"]]
    
    sample_counts <- as.data.frame(table(train_labels)) %>%
      dplyr::rename(
        class = train_labels,
        n = Freq)
    
    n_case <- sample_counts$n[sample_counts$class == "disease"]
    n_control <- sample_counts$n[sample_counts$class == "control"]
    n_outlier <- length(obj@OutlierSamples)
    
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
      out_dir = auc_dir,
      plot = plot)

    # Calculate Brier score
    brier_rf <- calculate_brier(y_true, y_prob_rf)
    brier_xgb <- calculate_brier(y_true, y_prob_xgb)

    # Plot calibration curve
    cal <- plot_calibration(
      comparison_key = val_key,
      y_true = y_true,
      y_prob_rf = y_prob_rf,
      y_prob_xgb = y_prob_xgb,
      out_dir = cal_dir,
      plot = plot)
    
    # Plot confusion matrix
    mat_rf = val_rf$conf.matrix$table
    mat_xgb = val_xgb$conf.matrix$table
    
    plot_matrix(
      comparison_key = val_key,
      mat_rf = mat_rf,
      mat_xgb = mat_xgb,
      out_dir = mat_dir,
      plot = plot)
    
    # TODO extract case v control sample size for figures
    
    val_res[[val_key]] <- list(
      
      meta = list(
        train_name = val$train_name[i],
        train_cohort = val$train_cohort[i],
        train_disease = val$train_disease[i],
        val_name = val$val_name[i],
        val_cohort = val$val_cohort[i],
        val_disease = val$val_disease[i],
        comparison_type = val$comparison_type[i],
        relationship = get_relationship(val$train_disease[i], val$val_disease[i]),
        n_case = n_case,
        n_control = n_control,
        n_outlier = n_outlier),
      
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
    delta_brier = results$rf$brier  - results$xgb$brier,
    delta_balanced_accuracy = results$rf$balanced_accuracy  - results$xgb$balanced_accuracy,
    delta_sensitivity = results$rf$sensitivity  - results$xgb$sensitivity,
    delta_precision = results$rf$precision  - results$xgb$precision,
    delta_f1 = results$rf$f1  - results$xgb$f1)
  
  train_analysis[[comparison_key]] <- c(meta, analysis_res)
    
  }
  
  analysed_df <- do.call(rbind, lapply(train_analysis, as.data.frame))
  
  # TODO simplify

  # Calculate delta quartiles
  roc_quartiles <- quantile(
    analysed_df$delta_roc_auc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  pr_quartiles <- quantile(
    analysed_df$delta_pr_auc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  brier_quartiles <- quantile(
    analysed_df$delta_brier, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  balanced_accuracy_quartiles <- quantile(
    analysed_df$delta_balanced_accuracy, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  sensitivity_quartiles <- quantile(
    analysed_df$delta_sensitivity, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  precision_quartiles <- quantile(
    analysed_df$delta_precision, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  f1_quartiles <- quantile(
    analysed_df$delta_f1, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  # Test whether median metrics are different from zero
  wilcox_roc <- wilcox.test(
    analysed_df$delta_roc_auc, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_pr <- wilcox.test(
    analysed_df$delta_pr_auc, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_brier <- wilcox.test(
    analysed_df$delta_brier, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_balanced_accuracy <- wilcox.test(
    analysed_df$delta_balanced_accuracy, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_sensitivity <- wilcox.test(
    analysed_df$delta_sensitivity, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_precision <- wilcox.test(
    analysed_df$delta_precision, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_f1 <- wilcox.test(
    analysed_df$delta_f1, mu = 0, alternative = "two.sided", exact = FALSE)
  
  # FDR correction 
  p_values <- list(
    roc_auc = wilcox_roc$p.value,
    pr_auc = wilcox_pr$p.value,
    brier = wilcox_brier$p.value,
    balanced_accuracy = wilcox_balanced_accuracy$p.value,
    sensitivity = wilcox_sensitivity$p.value,
    precision = wilcox_precision$p.value,
    f1 = wilcox_f1$p.value)
  
  p_adj <- p.adjust(p_values, method = "fdr")
  p_adj <- as.list(p_adj)
  
  # TODO determine if this is needed Calculate mean and CIs
  # mean(analysed_df$delta_roc_auc)
  # sd(analysed_df$delta_roc_auc)
  # roc_ttest <- t.test(analysed_df$delta_roc_auc)$conf.int
  # 
  # mean(analysed_df$delta_pr_auc)
  # sd(analysed_df$delta_pr_auc)
  # pr_ttest <- t.test(analysed_df$delta_pr_auc)$conf.int
  # 
  # mean(analysed_df$delta_brier)
  # sd(analysed_df$delta_brier)
  # brier_ttest <- t.test(analysed_df$delta_brier)$conf.int
  
  # Calculate percentages of XGB wins
  pct_xgb_roc_win <- mean(analysed_df$delta_roc_auc < 0, na.rm = TRUE) * 100
  pct_xgb_pr_win <- mean(analysed_df$delta_pr_auc < 0, na.rm = TRUE) * 100
  pct_xgb_brier_win <- mean(analysed_df$delta_brier > 0, na.rm = TRUE) * 100
  
  analysed_summary <- list(
    nrow = nrow(analysed_df),
    iqr1_delta_roc_auc = roc_quartiles["25%"], 
    median_delta_roc_auc = roc_quartiles["50%"], 
    iqr3_delta_roc_auc = roc_quartiles["75%"], 
    iqr1_delta_pr_auc = pr_quartiles["25%"], 
    median_delta_pr_auc = pr_quartiles["50%"], 
    iqr3_delta_pr_auc = pr_quartiles["75%"], 
    iqr1_delta_brier = brier_quartiles["25%"], 
    median_delta_brier = brier_quartiles["50%"], 
    iqr3_delta_brier = brier_quartiles["75%"], 
    iqr1_delta_balanced_accuracy = balanced_accuracy_quartiles["25%"], 
    median_delta_balanced_accuracy = balanced_accuracy_quartiles["50%"], 
    iqr3_delta_balanced_accuracy = balanced_accuracy_quartiles["75%"], 
    iqr1_delta_sensitivity = sensitivity_quartiles["25%"], 
    median_delta_sensitivity = sensitivity_quartiles["50%"], 
    iqr3_delta_sensitivity = sensitivity_quartiles["75%"], 
    iqr1_delta_precision = precision_quartiles["25%"], 
    median_delta_precision = precision_quartiles["50%"], 
    iqr3_delta_precision = precision_quartiles["75%"], 
    iqr1_delta_f1 = f1_quartiles["25%"], 
    median_delta_f1 = f1_quartiles["50%"], 
    iqr3_delta_f1 = f1_quartiles["75%"], 
    
    roc_pval = p_values$roc_auc,
    roc_padj = p_adj$roc_auc,
    pr_pval = p_values$pr_auc,
    pr_padj = p_adj$pr_auc,
    brier_pval = p_values$brier,
    brier_padj = p_adj$brier,
    balanced_accuracy_pval = p_values$balanced_accuracy,
    balanced_accuracy_padj = p_adj$balanced_accuracy,
    sensitivity_pval = p_values$sensitivity,
    sensitivity_padj = p_adj$sensitivity,
    precision_pval = p_values$precision,
    precision_padj = p_adj$precision,
    f1_pval = p_values$f1,
    f1_padj = p_adj$f1,
    
    pct_xgb_roc_win = pct_xgb_roc_win,
    pct_xgb_pr_win = pct_xgb_pr_win,
    pct_xgb_brier_win = pct_xgb_brier_win
    # roc_ttest_lower = roc_ttest[1],
    # roc_ttest_higher = roc_ttest[2],
    # pr_ttest_lower = pr_ttest[1],
    # pr_ttest_higher = pr_ttest[2],
    # brier_ttest_lower = brier_ttest[1],
    # brier_ttest_higher = brier_ttest[2])
  )
  
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
      delta_brier = results$rf$brier  - results$xgb$brier,
      delta_balanced_accuracy = results$rf$balanced_accuracy  - results$xgb$balanced_accuracy,
      delta_sensitivity = results$rf$sensitivity  - results$xgb$sensitivity,
      delta_precision = results$rf$precision  - results$xgb$precision,
      delta_f1 = results$rf$f1  - results$xgb$f1)
    
    val_analysis[[comparison_key]] <- c(meta, analysis_res)
  }
  
  analysed_df <- do.call(rbind, lapply(val_analysis, as.data.frame))
  
  # Calculate delta quartiles
  roc_quartiles <- quantile(
    analysed_df$delta_roc_auc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  pr_quartiles <- quantile(
    analysed_df$delta_pr_auc, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  brier_quartiles <- quantile(
    analysed_df$delta_brier, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  balanced_accuracy_quartiles <- quantile(
    analysed_df$delta_balanced_accuracy, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  sensitivity_quartiles <- quantile(
    analysed_df$delta_sensitivity, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  precision_quartiles <- quantile(
    analysed_df$delta_precision, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  f1_quartiles <- quantile(
    analysed_df$delta_f1, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  # Test whether median metrics are different from zero
  wilcox_roc <- wilcox.test(
    analysed_df$delta_roc_auc, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_pr <- wilcox.test(
    analysed_df$delta_pr_auc, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_brier <- wilcox.test(
    analysed_df$delta_brier, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_balanced_accuracy <- wilcox.test(
    analysed_df$delta_balanced_accuracy, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_sensitivity <- wilcox.test(
    analysed_df$delta_sensitivity, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_precision <- wilcox.test(
    analysed_df$delta_precision, mu = 0, alternative = "two.sided", exact = FALSE)
  
  wilcox_f1 <- wilcox.test(
    analysed_df$delta_f1, mu = 0, alternative = "two.sided", exact = FALSE)
  
  # FDR correction 
  p_values <- list(
    roc_auc = wilcox_roc$p.value,
    pr_auc = wilcox_pr$p.value,
    brier = wilcox_brier$p.value,
    balanced_accuracy = wilcox_balanced_accuracy$p.value,
    sensitivity = wilcox_sensitivity$p.value,
    precision = wilcox_precision$p.value,
    f1 = wilcox_f1$p.value)
  
  p_adj <- p.adjust(p_values, method = "fdr")
  p_adj <- as.list(p_adj)
  
  # FDR correction 
  p_values <- list(
    roc_auc = wilcox_roc$p.value,
    pr_auc = wilcox_pr$p.value,
    brier = wilcox_brier$p.value,
    balanced_accuracy = wilcox_balanced_accuracy$p.value,
    sensitivity = wilcox_sensitivity$p.value,
    precision = wilcox_precision$p.value,
    f1 = wilcox_f1$p.value)
  
  p_adj <- p.adjust(p_values, method = "fdr")
  p_adj <- as.list(p_adj)
  
  # Calculate mean and CIs
  # mean(analysed_df$delta_roc_auc)
  # sd(analysed_df$delta_roc_auc)
  # roc_ttest <- t.test(analysed_df$delta_roc_auc)$conf.int
  # 
  # mean(analysed_df$delta_pr_auc)
  # sd(analysed_df$delta_pr_auc)
  # pr_ttest <- t.test(analysed_df$delta_pr_auc)$conf.int
  # 
  # mean(analysed_df$delta_brier)
  # sd(analysed_df$delta_brier)
  # brier_ttest <- t.test(analysed_df$delta_brier)$conf.int
  
  # Calculate percentages of XGB wins
  pct_xgb_roc_win <- mean(analysed_df$delta_roc_auc < 0, na.rm = TRUE) * 100
  pct_xgb_pr_win <- mean(analysed_df$delta_pr_auc < 0, na.rm = TRUE) * 100
  pct_xgb_brier_win <- mean(analysed_df$delta_brier > 0, na.rm = TRUE) * 100
  
  analysed_summary <- list(
    nrow = nrow(analysed_df),
    iqr1_delta_roc_auc = roc_quartiles["25%"], 
    median_delta_roc_auc = roc_quartiles["50%"], 
    iqr3_delta_roc_auc = roc_quartiles["75%"], 
    iqr1_delta_pr_auc = pr_quartiles["25%"], 
    median_delta_pr_auc = pr_quartiles["50%"], 
    iqr3_delta_pr_auc = pr_quartiles["75%"], 
    iqr1_delta_brier = brier_quartiles["25%"], 
    median_delta_brier = brier_quartiles["50%"], 
    iqr3_delta_brier = brier_quartiles["75%"], 
    iqr1_delta_balanced_accuracy = balanced_accuracy_quartiles["25%"], 
    median_delta_balanced_accuracy = balanced_accuracy_quartiles["50%"], 
    iqr3_delta_balanced_accuracy = balanced_accuracy_quartiles["75%"], 
    iqr1_delta_sensitivity = sensitivity_quartiles["25%"], 
    median_delta_sensitivity = sensitivity_quartiles["50%"], 
    iqr3_delta_sensitivity = sensitivity_quartiles["75%"], 
    iqr1_delta_precision = precision_quartiles["25%"], 
    median_delta_precision = precision_quartiles["50%"], 
    iqr3_delta_precision = precision_quartiles["75%"], 
    iqr1_delta_f1 = f1_quartiles["25%"], 
    median_delta_f1 = f1_quartiles["50%"], 
    iqr3_delta_f1 = f1_quartiles["75%"], 
    
    roc_pval = p_values$roc_auc,
    roc_padj = p_adj$roc_auc,
    pr_pval = p_values$pr_auc,
    pr_padj = p_adj$pr_auc,
    brier_pval = p_values$brier,
    brier_padj = p_adj$brier,
    balanced_accuracy_pval = p_values$balanced_accuracy,
    balanced_accuracy_padj = p_adj$balanced_accuracy,
    sensitivity_pval = p_values$sensitivity,
    sensitivity_padj = p_adj$sensitivity,
    precision_pval = p_values$precision,
    precision_padj = p_adj$precision,
    f1_pval = p_values$f1,
    f1_padj = p_adj$f1,
    
    pct_xgb_roc_win = pct_xgb_roc_win,
    pct_xgb_pr_win = pct_xgb_pr_win,
    pct_xgb_brier_win = pct_xgb_brier_win
  )
  
  
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
  
  plot_dir <- file.path(stat_res_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("\nComparing performance across categories")
  
  lapply(
    c(stat_res_dir),
    dir.create, recursive = TRUE, showWarnings = FALSE)
  
  train_analysed <- train_analysis$analysed
  val_analysed <- val_analysis$analysed
  
  train_analysed_summary <- train_analysis$analysed_summary
  val_analysed_summary <- val_analysis$analysed_summary
  
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
      n_case, n_control, n_outlier, mod, roc_auc, pr_auc, brier) %>%
    pivot_wider(
      names_from = mod,
      values_from = c(roc_auc, pr_auc, brier),
      names_sep = "_") %>%
    mutate(
      delta_roc_auc = roc_auc_RF - roc_auc_XGB,
      delta_pr_auc = pr_auc_RF - pr_auc_XGB,
      delta_brier = brier_RF - brier_XGB)
  
  metrics <- c("roc_auc", "pr_auc", "brier")
  models <- c("RF", "XGB")
  
  message(" > Testing difference in RF/XGB performance across comparison types")
  
  kruskal_res <- list(
    roc = list(),
    pr = list(),
    brier = list())
  
  for (model in models) {
    for (metric in metrics) {
      formula <- as.formula(
        paste0(metric, "_", model, " ~ comparison_type"))
      
      test <- kruskal.test(
        formula,
        data = combo_wide)
      
      kruskal_res[[metric]][[model]] <- list(
        statistic = unname(test$statistic),
        df = unname(test$parameter),
        p_value = test$p.value)
    }
  }
  
  for (metric in names(kruskal_res)) {
    
    p_adj <- p.adjust(
      sapply(kruskal_res[[metric]], `[[`, "p_value"),
      method = "fdr")
    
    for (model in names(kruskal_res[[metric]])) {
      
      kruskal_res[[metric]][[model]]$p_adj <-
        unname(p_adj[model])
    }
  }
  
  kruskal_df <- do.call(
    rbind,
    lapply(names(kruskal_res), function(metric) {
      
      do.call(
        rbind,
        lapply(names(kruskal_res[[metric]]), function(model) {
          
          data.frame(
            metric = metric,
            model = model,
            statistic = kruskal_res[[metric]][[model]]$statistic,
            df = kruskal_res[[metric]][[model]]$df,
            p_value = kruskal_res[[metric]][[model]]$p_value,
            p_adj = kruskal_res[[metric]][[model]]$p_adj,
            stringsAsFactors = FALSE)
        }))
    }))
  
  dunn_res <- list()
  
  for (metric in metrics) {
    
    dunn_res[[metric]] <- list()
    
    for (model in models) {
      
      formula <- as.formula(
        paste0(metric, "_", model, " ~ comparison_type"))
      
      dunn_test <- dunnTest(
        formula,
        data = combo_wide,
        method = "bh")
      
      dunn_res[[metric]][[model]] <- dunn_test$res
    }
  }
  
  dunn_df <- do.call(
    rbind,
    lapply(names(dunn_res), function(metric) {
      do.call(
        rbind,
        lapply(names(dunn_res[[metric]]), function(model) {
          cbind(
            metric = metric,
            model = model,
            dunn_res[[metric]][[model]])}))}))
  
  fn <- "1a. kruskal.csv"
  
  write.csv(
    kruskal_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)

  fn <- "1b. post_hoc_dunn.csv"
  
  write.csv(
    dunn_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - post-hoc results saved as: ", fn)
  
  message(" > Testing difference in median ΔRF-XGB within comparison types")
  
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
      median_delta = median(tmp$delta_roc_auc, na.rm = TRUE),
      p_value = roc_res$p.value,
      ci_lower = roc_res$conf.int[1],
      ci_higher = roc_res$conf.int[2])
    
    wilcox_res$pr[[type]] <- list(
      median_delta = median(tmp$delta_pr_auc, na.rm = TRUE),
      p_value = pr_res$p.value,
      ci_lower = pr_res$conf.int[1],
      ci_higher = pr_res$conf.int[2])
    
    wilcox_res$brier[[type]] <- list(
      median_delta = median(tmp$delta_brier, na.rm = TRUE),
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
            median_delta = wilcox_res[[metric]][[type]]$median_delta,
            p_value = wilcox_res[[metric]][[type]]$p_value,
            p_adj = wilcox_res[[metric]][[type]]$p_adj,
            ci_lower = wilcox_res[[metric]][[type]]$ci_lower,
            ci_higher = wilcox_res[[metric]][[type]]$ci_higher,
            stringsAsFactors = FALSE)
        }))
    }))
  
  fn <- "2. wilcox.csv"
  
  write.csv(
    wilcox_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  message(" > Testing difference in median ΔRF-XGB across comparison types")
  
  delta_kruskal_res <- lapply(metrics, function(metric) {
    res <- kruskal.test(
      combo_wide[[paste0("delta_", metric)]] ~ combo_wide$comparison_type)
    
    list(
      statistic = unname(res$statistic),
      df = unname(res$parameter),
      p_value = res$p.value)
  })
  
  names(delta_kruskal_res) <- metrics
  
  pairwise_res <- lapply(metrics, function(metric) {
    res <- pairwise.wilcox.test(
      combo_wide[[paste0("delta_", metric)]],
      combo_wide$comparison_type,
      p.adjust.method = "BH")
    
    res$p.value
  })
  
  names(pairwise_res) <- metrics
  
  delta_kruskal_df <- do.call(
    rbind,
    lapply(names(delta_kruskal_res), function(metric) {
      data.frame(
        metric = metric,
        statistic = delta_kruskal_res[[metric]]$statistic,
        df = delta_kruskal_res[[metric]]$df,
        p_value = delta_kruskal_res[[metric]]$p_value,
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
  
  fn <- "3a. kruskal.csv"
  
  write.csv(
    delta_kruskal_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  fn <- "3b. pairwise.csv"
  
  write.csv(
    pairwise_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - post-hoc results saved as: ", fn)
  
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
  
  message(" > Testing correlations between biological relationship and RF/XGB performance")
  
  combo_wide_num <- combo_wide
  
  relationship_map <- c(
    "same" = 1,
    "continuum" = 2,
    "shared" = 3,
    "unrelated" = 4
  )
  
  combo_wide_num$relationship <- unname(
    relationship_map[combo_wide_num$relationship]
  )
  
  cor_rel <- list()
  
  for (model in models) {
    
    cor_rel[[model]] <- lapply(metrics, function(metric) {
      
      res <- cor.test(
        combo_wide_num$relationship,
        combo_wide_num[[paste0(metric, "_", model)]],
        method = "spearman",
        exact = FALSE)
      
      list(
        rho = unname(res$estimate),
        p_value = res$p.value,
        n = sum(complete.cases(
          combo_wide_num$relationship,
          combo_wide_num[[paste0(metric, "_", model)]])))
    })
    
    names(cor_rel[[model]]) <- metrics
    
    p_adj <- p.adjust(
      sapply(cor_rel[[model]], `[[`, "p_value"),
      method = "fdr")
    
    for (metric in metrics) {
      cor_rel[[model]][[metric]]$p_adj <- unname(p_adj[metric])
    }
  } 
  
  cor_rel_df <- do.call(
    rbind,
    lapply(names(cor_rel), function(model) {
      
      do.call(
        rbind,
        lapply(names(cor_rel[[model]]), function(metric) {
          
          data.frame(
            model = model,
            metric = metric,
            rho = cor_rel[[model]][[metric]]$rho,
            p_value = cor_rel[[model]][[metric]]$p_value,
            p_adj = cor_rel[[model]][[metric]]$p_adj,
            n = cor_rel[[model]][[metric]]$n,
            stringsAsFactors = FALSE)
          
        }))
    }))
  
  fn <- "5. cor_rel_vs_auc.csv"
  
  write.csv(
    cor_rel_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", file.path(stat_res_dir, fn))
  
  message(" > Testing difference in RF/XGB performance across biological relationships")
  rel_kruskal_res <- list(
    roc = list(),
    pr = list(),
    brier = list())
  
  for (model in models) {
    for (metric in metrics) {
      formula <- as.formula(
        paste0(metric, "_", model, " ~ relationship"))
      
      test <- kruskal.test(
        formula,
        data = combo_wide)
      
      rel_kruskal_res[[metric]][[model]] <- list(
        statistic = unname(test$statistic),
        df = unname(test$parameter),
        p_value = test$p.value)
    }
  }
  
  for (metric in names(rel_kruskal_res)) {
    
    p_adj <- p.adjust(
      sapply(rel_kruskal_res[[metric]], `[[`, "p_value"),
      method = "fdr")
    
    for (model in names(rel_kruskal_res[[metric]])) {
      
      rel_kruskal_res[[metric]][[model]]$p_adj <-
        unname(p_adj[model])
    }
  }
  
  rel_kruskal_df <- do.call(
    rbind,
    lapply(names(rel_kruskal_res), function(metric) {
      
      do.call(
        rbind,
        lapply(names(rel_kruskal_res[[metric]]), function(model) {
          
          data.frame(
            metric = metric,
            model = model,
            statistic = rel_kruskal_res[[metric]][[model]]$statistic,
            df = rel_kruskal_res[[metric]][[model]]$df,
            p_value = rel_kruskal_res[[metric]][[model]]$p_value,
            p_adj = rel_kruskal_res[[metric]][[model]]$p_adj,
            stringsAsFactors = FALSE)
        }))
    }))
  
  rel_dunn_res <- list()
  
  for (metric in metrics) {
    
    rel_dunn_res[[metric]] <- list()
    
    for (model in models) {
      
      formula <- as.formula(
        paste0(metric, "_", model, " ~ relationship"))
      
      dunn_test <- dunnTest(
        formula,
        data = combo_wide,
        method = "bh")
      
      rel_dunn_res[[metric]][[model]] <- dunn_test$res
    }
  }
  
  rel_dunn_df <- do.call(
    rbind,
    lapply(names(rel_dunn_res), function(metric) {
      do.call(
        rbind,
        lapply(names(rel_dunn_res[[metric]]), function(model) {
          cbind(
            metric = metric,
            model = model,
            rel_dunn_res[[metric]][[model]])}))}))
  
  fn <- "6a. kruskal.csv"
  
  write.csv(
    rel_kruskal_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  fn <- "6b. post_hoc_dunn.csv"
  
  write.csv(
    rel_dunn_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - post-hoc results saved as: ", fn)
  
  message(" > Testing correlations between biological relationship and ΔRF-XGB")
  
  cor_delta_res <- lapply(metrics, function(metric) {
    
    res <- cor.test(
      combo_wide_num$relationship,
      combo_wide_num[[paste0("delta_", metric)]],
      method = "spearman",
      exact = FALSE)
    
    list(
      rho = unname(res$estimate),
      p_value = res$p.value,
      n = sum(complete.cases(
        combo_wide_num$relationship,
        combo_wide_num[[paste0("delta_", metric)]])))
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
  
  fn <- "7. cor_rel_vs_delta.csv"
  
  write.csv(
    cor_rel_vs_delta_df,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  message("\nSummarising statistics across comparison types")
  message(" > Calculating effects sizes")
  
  effect_summary <- combo_wide_num %>%
    group_by(comparison_type) %>%
    summarise(
      n = n(),
      median_delta_roc_auc = median(delta_roc_auc, na.rm = TRUE),
      IQR_delta_roc_auc = IQR(delta_roc_auc, na.rm = TRUE),
      median_delta_pr_auc = median(delta_pr_auc, na.rm = TRUE),
      IQR_delta_pr_auc = IQR(delta_pr_auc, na.rm = TRUE),
      median_delta_brier = median(delta_brier, na.rm = TRUE),
      IQR_delta_brier = IQR(delta_brier, na.rm = TRUE))
  
  fn <- "8. effect_summary.csv"
  
  write.csv(
    effect_summary,
    file = file.path(stat_res_dir, fn),
    row.names = FALSE)
  
  message("   - results saved as: ", fn)
  
  message(" > Calculating variances")
  
  variance <- combo_results_df %>%
    group_by(mod, comparison_type) %>%
    summarise(
      sd_roc = sd(roc_auc, na.rm = TRUE),
      sd_pr = sd(pr_auc, na.rm = TRUE),
      sd_brier = sd(brier, na.rm = TRUE),
      .groups = "drop")
  
  fn <- "9. variance_summary.csv"
  
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
      kruskal = kruskal_df,
      dunn = dunn_df,
      wilcox = wilcox_df,
      delta_kruskal = delta_kruskal_df,
      pairwise = pairwise_df,
      cor_rf_vs_xgb = cor_rf_vs_xgb_df,
      cor_rel = cor_rel_df,
      rel_kruskal = rel_kruskal_df,
      rel_dunn = rel_dunn_df,
      cor_rel_vs_delta = cor_rel_vs_delta_df,
      effect_summary = effect_summary,
      variance = variance,
      train_summary = train_analysed_summary,
      val_summary = val_analysed_summary)))
}


# Define main function ---------------------------------------------------------
# Check for missing results, compile all logs and save as CSV
run_analysis <- function(
    prep_dir,
    in_dir,
    out_dir,
    plot) {
  
  train_res_dir <- file.path(out_dir, "1. train_results")
  val_res_dir <- file.path(out_dir, "2. val_results")
  stat_res_dir <- file.path(out_dir, "3. stat_results")

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
    
    train <- extract_training(mod, train_res_dir, plot)
    train_res[[train$meta$train_name]] <- train
    
    # Extract validation results
    val <- extract_validation(mod, train_res, val_res_dir, plot)
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
  stat_res <- analyse_by_categories(
    train_res, train_analysis, 
    val_res, val_analysis, 
    stat_res_dir)
  
  analysis_res <- list(
    train_results = train_res,
    train_analysis = train_analysis,
    validation_results = val_res,
    validation_analysis = val_analysis,
    statistic_results = stat_res)
  
  saveRDS(analysis_res, file = file.path(out_dir, "analysis_res.RDS"))
  
  return(analysis_res)
}


# Execute ----------------------------------------------------------------------
analysis_res <- run_analysis(prep_dir, in_dir, out_dir, plot = FALSE)
