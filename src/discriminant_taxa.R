# Load lefser package
# Khleborodova, A., Gamboa-Tuz, S.D., Ramos, M., Segata, N., Waldron, L. and Oh, 
# S. 2024. lefser:implementation of metagenomic biomarker discovery tool, LEfSe, 
# in R P. Robinson, ed.Bioinformatics. 40(12), p.btae707.
library(lefser)

# Load curatedMetagenomicData package
# Pasolli, E., Schiffer, L., Manghi, P., Renson, A., Obenchain, V., Truong, D.T., 
# Beghini, F., Malik, F., Ramos, M., Dowd, J.B., Huttenhower, C., Morgan, M., 
# Segata, N. and Waldron, L. 2017. Accessible, curated metagenomic data through 
# ExperimentHub. Nature Methods. 14(11), pp.1023–1024.
library(curatedMetagenomicData)

# Load other packages
library(patchwork)
library(ggplot2)

# Define datasets (CRC-healthy cohorts)
datasets <- c(
  "GuptaA_2019",
  "ThomasAM_2018a",
  "ThomasAM_2018b",
  "ThomasAM_2019_c",
  "VogtmannE_2016",
  "WirbelJ_2018"
)

# Function: load datasets
load_datasets <- function(datasets) {
  
  loaded <- lapply(datasets, function(study) {
    message("Loading: ", study)
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
  
  names(loaded) <- datasets
  return(loaded)
}

all_dat <- load_datasets(datasets)

# Function: get contingency tables
get_contingency_tables <- function(all_dat) {
  lapply(names(all_dat), function(study) {
    dat <- all_dat[[study]]
    if (is.null(dat)) return(NULL)
    cd <- as.data.frame(colData(dat))
    if (!all(c("age_category", "study_condition") %in% names(cd))) {
      message("Missing columns in: ", study)
      return(NULL)
    }
    table(cd$age_category, cd$study_condition)
  }) |> setNames(names(all_dat))
}

# Function: run lefser
run_lefser <- function(all_dat) {
  lapply(names(all_dat), function(study) {
    
    dat <- all_dat[[study]]
    
    if (is.null(dat)) return(NULL)
    message("Running LEfSe: ", study)
    
    # Remove adenoma samples if present
    if ("adenoma" %in% dat$study_condition) {
      message("  Removing adenoma samples from: ", study)
      dat <- dat[, dat$study_condition != "adenoma"]
    }
    
    tn       <- get_terminal_nodes(rownames(dat))
    dattn    <- dat[tn, ]
    dattn_ra <- relativeAb(dattn)
    
    set.seed(1234)
    tryCatch(
      suppressWarnings(
        lefser(dattn_ra, classCol = "study_condition", subclassCol = NULL)
      ),
      error = function(e) {
        message("  ERROR in ", study, ": ", e$message)
        return(NULL)
      }
    )
  }) |> setNames(names(all_dat))
}

# Run functions
contingency_tables <- get_contingency_tables(all_dat)
lefser_results     <- run_lefser(all_dat)

# Contingency tables
for (study in names(contingency_tables)) {
  cat("\n======", study, "======\n")
  if (is.null(contingency_tables[[study]])) {
    cat("No table available\n")
  } else {
    print(contingency_tables[[study]])
  }
}

# Plot lefser results
sapply(lefser_results, function(x) if (is.null(x)) 0 else nrow(x))

plots <- lapply(names(lefser_results), function(study) {
  res <- lefser_results[[study]]
  if (is.null(res) || nrow(res) == 0) return(NULL)
  lefserPlot(res) +
    ggtitle(study) +
    theme(plot.title = element_text(face = "bold"))
})

names(plots) <- names(lefser_results)

# Remove NULLs
plots_valid <- Filter(Negate(is.null), plots)

if (length(plots_valid) == 0) {
  message("No plots to display")
} else {
  wrap_plots(plots_valid, ncol = 1)
}
