# Utility functions for reading and writing curatedMetagenomicData cohorts

library(curatedMetagenomicData)
library(SummarizedExperiment)


#' Load every dataset listed in datasets.R
#' 
#' @param cohorts  Named list as defined in datasets.R 
#' @return Named list of SummarizedExperiment objects (NULL on failure)

load_cohorts <- function(cohorts) {
  
  loaded <- lapply(cohorts, function(cfg) {
    study <- cfg$dataset
    message("\n[handler] Loading: ", study)
    
    tryCatch(
      curatedMetagenomicData(
        paste0(study, ".relative_abundance"),
        dryrun = FALSE
      )[[1]],
      error = function(e) {
        message("  ERROR loading ", study, ": ", e$message)
        return(NULL)
      }
    )
  })
  
  names(loaded) <- names(cohorts)
  loaded
}



#' Assess pair-wise class balance for Lefser analysis
#' 
#' Returns a data.frame with one row per group containing sample counts,
#' the majority:minority imbalance ratio, and three QC flags:
#'
#'   * \code{skip_min_n}   : TRUE if either group is below \code{min_n}
#'                           (tests unreliable)
#'   * \code{warn_ratio}   : TRUE if majority/minority ratio exceeds
#'                           \code{max_ratio} (LDA scores may be inflated)
#'   * \code{warn_small_n} : TRUE if either group is below \code{warn_n}
#'                           (low power, interpret cautiously)
#'
#' The attribute \code{"recommendation"} on the returned data.frame is one of:
#'   "SKIP"    - hard skip (min_n violated)
#'   "CAUTION" - run but interpret carefully (ratio or warn_n violated)
#'   "OK"      - proceed normally
#'
#' @param se         SummarizedExperiment (already subset to the two groups
#'                   of interest and prepared by \code{prepare_se}).
#' @param group_col  colData column carrying class labels.
#' @param min_n      Hard minimum samples per group (default 10).
#' @param warn_n     Soft minimum that triggers a caution warning (default 20).
#' @param max_ratio  Imbalance ratio (larger/smaller) that triggers caution
#'                   (default 3).
#' @return data.frame with columns: group, n, ratio, skip_min_n,
#'         warn_ratio, warn_small_n; attribute "recommendation".
check_class_balance <- function(se,
                                group_col,
                                min_n     = 10,
                                warn_n    = 20,
                                max_ratio = 3) {
  
  groups <- se[[group_col]]
  counts <- sort(table(groups), decreasing = TRUE)   # named integer
  
  if (length(counts) != 2L)
    stop("[handler] check_class_balance expects exactly 2 groups; got ",
         length(counts), ": ", paste(names(counts), collapse = ", "))
  
  n_maj <- counts[[1L]]
  n_min <- counts[[2L]]
  ratio <- n_maj / n_min
  
  df <- data.frame(
    group       = names(counts),
    n           = as.integer(counts),
    ratio       = c(round(ratio, 2), NA_real_),   # ratio on majority row
    skip_min_n  = as.integer(counts) < min_n,
    warn_ratio  = c(ratio > max_ratio, FALSE),
    warn_small_n = as.integer(counts) < warn_n,
    stringsAsFactors = FALSE
  )
  
  recommendation <-
    if (any(df$skip_min_n))       "SKIP"
  else if (any(df$warn_ratio) ||
           any(df$warn_small_n)) "CAUTION"
  else                           "OK"
  
  attr(df, "recommendation") <- recommendation
  attr(df, "min_n")          <- min_n
  attr(df, "warn_n")         <- warn_n
  attr(df, "max_ratio")      <- max_ratio
  df
}


#' Filter samples for body site
#'
#' Reads \code{cfg$body_site_col} and \code{cfg$keep_body_sites} from the
#' cohort config. If \code{keep_body_sites} is empty or the column is absent,
#' the SE is returned unchanged with a message.
#'
#' @param se   SummarizedExperiment
#' @param cfg  One cohort config list from datasets.R
#' @return Filtered SummarizedExperiment
filter_body_site <- function(se, cfg) {
  
  keep_sites <- cfg$keep_body_sites
  site_col   <- cfg$body_site_col
  
  # If keep_sites not defined, disable filtering
  if (is.null(keep_sites) || length(keep_sites) == 0) {
    message("  Body-site filtering disabled for: ", cfg$dataset)
    return(se)
  }
  
  # Warning if body_site_col doesn't exist
  if (!site_col %in% names(colData(se))) {
    message("  WARNING: body_site column '", site_col,
            "' not found in colData for ", cfg$dataset,
            " — skipping body-site filter.")
    return(se)
  }
  
  sites_present <- unique(se[[site_col]])
  keep_mask     <- se[[site_col]] %in% keep_sites
  n_dropped     <- sum(!keep_mask)
  
  if (n_dropped > 0) {
    dropped_sites <- setdiff(sites_present, keep_sites)
    message("  Body-site filter: dropping ", n_dropped, " sample(s) from site(s): ",
            paste(dropped_sites, collapse = ", "))
  } else {
    message("  Body-site filter: all samples already from site(s): ",
            paste(keep_sites, collapse = ", "))
  }
  
  se[, keep_mask]
}


#' Pretty-print the output of \code{check_class_balance}.
#'
#' @param balance_df  data.frame returned by \code{check_class_balance}.
#' @param label       Optional string label (e.g. comparison name) printed
#'                    as a header.
#' @return \code{balance_df} invisibly.
print_balance_report <- function(balance_df, label = NULL) {
  
  rec   <- attr(balance_df, "recommendation")
  min_n <- attr(balance_df, "min_n")
  warn_n <- attr(balance_df, "warn_n")
  max_r <- attr(balance_df, "max_ratio")
  
  # Colour-coded prefix (works in terminals that support ANSI)
  prefix <- switch(rec,
                   OK      = "\033[32m[OK]\033[0m     ",
                   CAUTION = "\033[33m[CAUTION]\033[0m",
                   SKIP    = "\033[31m[SKIP]\033[0m   "
  )
  
  if (!is.null(label))
    cat("\n  Balance check:", label, "\n")
  
  cat(" ", prefix, "\n")
  
  # Print the per-group table
  display <- balance_df
  display$ratio      <- ifelse(is.na(display$ratio), "-",
                               as.character(display$ratio))
  display$skip_min_n  <- ifelse(display$skip_min_n,  "YES (SKIP)",   "no")
  display$warn_ratio  <- ifelse(display$warn_ratio,  "YES (warn)",   "no")
  display$warn_small_n <- ifelse(display$warn_small_n, "YES (warn)", "no")
  print(display, row.names = FALSE)
  
  # Thresholds reminder
  cat(sprintf(
    "  Thresholds: min_n=%d (hard skip) | warn_n=%d | max_ratio=%.0f:1\n",
    min_n, warn_n, max_r
  ))
  
  invisible(balance_df)
}