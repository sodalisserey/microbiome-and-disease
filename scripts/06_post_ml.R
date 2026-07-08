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
library(readr)
library(grid)
library(ggplot2)
source("R/utils.R")

# Define input/output directories ----------------------------------------------
prep_dir <- "results/04_pre_ml"
in_dir <- "results/05_ml"
# in_dir <- "results/hpc/05_ml_update"
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

  return(cal)
}

plot_matrix <- function(
    comparison_key,
    mat_rf,
    mat_xgb,
    out_dir) {

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

extract_training <- function(models, train_res_dir) {
  
  message("\nExtracting model training results")
  
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

  train_res <- list()
  train_dat <- list()

  for (file in models) {
    
    mod <- readRDS(file)
    meta <- mod$metadata
    obj <- mod$obj
    
    comparison_key <- meta$comparison
    
    # Extract predictions data
    y_true = obj@ModelData$TestLabel
    y_prob_rf = obj@EvaluateResult$RF$predictions[, "disease"]
    y_prob_xgb = obj@EvaluateResult$XGBoost$predictions
    
    # Plot ROC-AUC and PR-AUC
    auc <- calculate_auc(
      comparison_key,
      y_true,
      y_prob_rf,
      y_prob_xgb,
      auc_dir)
 
    # Calculate AUC CIs using bootstrapping
    suppressWarnings({
      suppressMessages({
        boot_ci_rf <- ci.auc(
          auc$roc_rf,
          method = "bootstrap",
          boot.n = 2000,
          conf.level = 0.95)
        
        boot_ci_xgb <- ci.auc(
          auc$roc_xgb,
          method = "bootstrap",
          boot.n = 2000,
          conf.level = 0.95)
      })})
    
    # Calculate Brier score
    brier_rf <- calculate_brier(y_true, y_prob_rf)
    brier_xgb <- calculate_brier(y_true, y_prob_xgb)
    
    # Plot calibration curve
    cal <- plot_calibration(
      comparison_key, 
      y_true, 
      y_prob_rf, 
      y_prob_xgb, 
      cal_dir)
    
    # Plot confusion matrix
    mat_rf = obj@EvaluateResult$RF$conf.matrix$table
    mat_xgb = obj@EvaluateResult$XGBoost$ConfusionMatrix$table

    plot_matrix(
      comparison_key,
      mat_rf,
      mat_xgb,
      mat_dir)
    
    # Extract training results
    train_res[[comparison_key]] <- list( 
      
      meta = list(
        comparison = meta$comparison,
        cohort = meta$cohort,
        disease = meta$disease,
        n_control = meta$n_control,
        n_disease = meta$n_disease),
      
      results = list(
        rf = list(
          mod = "RF",
          roc_auc = auc$roc_auc_rf,
          pr_auc = auc$pr_auc_rf,
          brier = brier_rf,
          ci_lower = as.numeric(auc$ci_rf[1]),
          ci_higher = as.numeric(auc$ci_rf[3]),
          boot_ci_lower = as.numeric(boot_ci_rf[1]),
          boot_ci_higher = as.numeric(boot_ci_rf[3]),
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
          pr_auc = auc$pr_auc_xgb,
          brier = brier_xgb,
          ci_lower = as.numeric(auc$ci_xgb[1]),
          ci_higher = as.numeric(auc$ci_xgb[3]),
          boot_ci_lower = as.numeric(boot_ci_xgb[1]),
          boot_ci_higher = as.numeric(boot_ci_xgb[3]),
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
          y_pred_rf = as.numeric(obj@EvaluateResult$RF$predictions),
          y_prob_xgb = y_prob_xgb,
          y_pred_xgb = as.numeric(obj@EvaluateResult$XGBoost$predictions))))
  }
  
  # Bind data frames and save as csv
  training_save <- do.call(
    rbind,
    lapply(train_res, function(x) {
      rbind(
        data.frame(x$meta, x$results$rf, stringsAsFactors = FALSE),
        data.frame(x$meta, x$results$xgb, stringsAsFactors = FALSE))
    }))

  write.csv(
    training_save,
    file = file.path(train_res_dir, "results.csv"),
    row.names = FALSE)

  message(" > results saved as: ", file.path(train_res_dir, "results.csv"))
  message(" > plots saved to: ", file.path(train_res_dir, "plots"))

  return(train_res)
}

analyse_training <- function(train_res, train_res_dir) {
  
  message("\nAnalysing extracted training results")
  
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
    
  train_analysis[[comparison_key]] <- list( 
    comparison = meta$comparison,
    cohort = meta$cohort,
    disease = meta$disease,
    n_control = meta$n_control,
    n_disease = meta$n_disease,
    delta_roc_auc = results$rf$roc_auc - results$xgb$roc_auc,
    roc_z_stat = as.numeric(delong_res$statistic),
    roc_ci_lower = as.numeric(delong_res$conf.int[1]),
    roc_ci_upper = as.numeric(delong_res$conf.int[2]),
    delta_pr_auc = results$rf$pr_auc - results$xgb$pr_auc,
    delta_brier = results$xgb$brier - results$rf$brier,
    roc_pval = as.numeric(delong_res$p.value))
  }
  
  analysed_df <- do.call(rbind, lapply(train_analysis, as.data.frame))
  
  # Calculate median delta
  median_delta_auc = median(analysed_df$delta_roc_auc, na.rm = TRUE)
  median_delta_pr_auc = median(analysed_df$delta_pr_auc, na.rm = TRUE)
  median_delta_brier = median(analysed_df$delta_brier, na.rm = TRUE)
  
  # Calculate percentage where XGB wins
  pct_xgb_auc_win <- mean(analysed_df$delta_roc_auc < 0, na.rm = TRUE) * 100
  pct_xgb_pr_win <- mean(analysed_df$delta_pr_auc < 0, na.rm = TRUE) * 100
  pct_xgb_brier_win <- mean(analysed_df$delta_brier < 0, na.rm = TRUE) * 100
  
  # Test whether median delta is different from zero
  wilcox_auc <- wilcox.test(
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
    roc_auc = wilcox_auc$p.value,
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
  
  analysed_summary <- list(
    nrow = nrow(analysed_df),
    median_delta_auc = median_delta_auc,
    median_delta_pr_auc = median_delta_pr_auc,
    median_delta_brier = median_delta_brier,
    pct_xgb_auc_win = pct_xgb_auc_win,
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
  
  write.csv(
    analysed_summary_df,
    file = file.path(train_res_dir, "analysed_summary.csv"),
    row.names = FALSE)
  
  message(" > results saved to: ", train_res_dir)

  return(list(
    analysed = train_analysis,
    analysed_summary = analysed_summary))
}
  
  
  # TODO: summarise median ΔAUC, Wilcoxon p-value, adjust p_values, %XGB wins
  


extract_validation <- function(models, out_dir) {
  
  # TODO still not working, must fix
  
  rf_val <- list()
  xgb_val <- list()
  
  for (file in models) {
    
    mod <- readRDS(file)
    obj <- mod$obj
    mod_val <- mod$validation
    val_rf <- obj@PredictResult$RF
    val_xgb <- obj@PredictResult$XGBoost
    
    # Extract RF validation results
    for (i in seq_along(val_rf)) {
      
      rf_val[[length(rf_val) + 1]] <- data.frame(
        model = file_path_sans_ext(basename(file)),
        train_cohort = mod_val$train_cohort,
        train_disease = mod_val$train_disease,
        val_cohort = mod_val$val_cohort[i],
        val_disease = mod_val$val_disease[i],
        comparison_type = mod_val$comparison_type[i],
        auroc = as.numeric(val_rf[[i]]$AUC))
    }
    
  }
  
  rf_val <- do.call(rbind, rf_val)
  xgb_val <- do.call(rbind, xgb_val)
  
  return(list(
    rf_val = rf_val,
    xgb_val = xgb_val))
}

# Define main function ---------------------------------------------------------
# Check for missing results, compile all logs and save as CSV
run_analysis <- function(
    prep_dir,
    in_dir,
    out_dir) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  train_res_dir <- file.path(out_dir, "train_results")
  dir.create(train_res_dir, recursive = TRUE, showWarnings = FALSE)
  
  
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
  
  # Extract and analyse training results
  train_res <- extract_training(models, train_res_dir)
  train_analysis <- analyse_training(train_res, train_res_dir)
  

  return(list(
    train_results = train_res,
    train_analysis = train_analysis))
  
  # return(train_res)
}


# Execute ----------------------------------------------------------------------
train_res <- run_analysis(prep_dir, in_dir, out_dir)



# result <- readRDS("results/05_ml/final/QinN_2014_cirrhosis_evaluated.rds")
# obj <- result$obj



# boot_auc <- function(data, idx) {
#   d <- data[idx, ]
#   as.numeric(pROC::auc(d$y_true, d$y_prob))
# }

