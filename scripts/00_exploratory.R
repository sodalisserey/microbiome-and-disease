# Load packages and dependencies -----------------------------------------------
library(SummarizedExperiment)
library(dplyr)
library(grid)
library(sf)
library(giscoR)
library(ggplot2)
source("R/utils.R")

# Define output directory ------------------------------------------------------
out_dir <- "results/00_exploratory"


# Define helper functions ------------------------------------------------------
plot_cohort_on_map <- function(cohorts, out_dir) {
  
  message("\nPlotting cohorts on world map")
  
  # Download country polygons
  world <- gisco_get_countries(
    resolution = "20",
    epsg = 4326)
  
  # Fix geometries
  world <- st_make_valid(world)
  
  # Create points inside countries
  suppressWarnings({
    suppressMessages({
      country_points <- st_point_on_surface(world)
    })})
  
  # Join cohort information
  country_points <- country_points |>
    left_join(
      cohort_locations,
      by = c("ISO3_CODE" = "country")) |>
    filter(!is.na(n_samples))
  
  # Save plot
  p <- ggplot() +
    geom_sf(
      data = world,
      fill = "grey90",
      colour = "white") +
    geom_sf(
      data = country_points,
      aes(size = n_samples),
      colour = "#6EC4C0") +
    theme_minimal() +
    labs(size = "Samples") +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(1.2, "cm"))

  ggsave(
    filename = file.path(out_dir, "cohort_world_map.pdf"),
    plot = p,
    width = 10,
    height = 6)
  
  message(" > Plot saved to: ", out_dir)
}

plot_cohort_samples <- function(cohorts, out_dir) {
  
  message("\nPlotting sample sizes")
  
  case_control <- lapply(names(cohorts), function(x){
    meta <- as.data.frame(colData(cohorts[[x]]))
    meta$cohort <- x
    meta$status <- ifelse(
      meta$study_condition == "control",
      "Control", "Case")
    meta
  })
  
  case_control <- bind_rows(case_control)
  
  cohort_order <- case_control |>
    dplyr::count(cohort) |>
    dplyr::arrange(n) |>
    dplyr::pull(cohort)
  
  case_control$cohort <- factor(
    case_control$cohort,
    levels = cohort_order)
  
  p <- ggplot(case_control,
         aes(cohort, fill = status)) +
    geom_bar(position = "stack") +
    labs(
      y = "Count", x = "",
      fill = "Status") +
    scale_fill_manual(values = c(
      "Control" = "#c5cb93ff", "Case" = "#DFA57E")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16))
  
  ggsave(
    filename = file.path(out_dir, "cohort_sample_size.pdf"),
    plot = p,
    width = 12,
    height = 4)
  
  message(" > Plot saved to: ", out_dir)
}


# Define main function ---------------------------------------------------------
run_exploratory <- function(cohorts, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_cohort_on_map(cohorts, out_dir)
  plot_cohort_samples(cohorts, out_dir)
}


# Load data --------------------------------------------------------------------
cohorts <- readRDS("data/primary.rds") 


# Execute ----------------------------------------------------------------------
run_exploratory(cohorts, out_dir)








