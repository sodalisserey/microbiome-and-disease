# Generate exploratory figures to show sample distribution and composition

# 1. Define output directory and functions
# 2. Load data using utils.R
# 3. Execute main function run_exploratory() which:
#   * Runs plot_distribution
#   * Runs plot_composition
#   * Exports plots as out_dir/[plot_name].pdf


# Load packages and dependencies -----------------------------------------------
library(SummarizedExperiment)
library(countrycode)
library(dplyr)
library(grid)
library(sf)
library(giscoR)
library(ggplot2)
source("R/utils.R")

# Define output directory ------------------------------------------------------
out_dir <- "results/00_exploratory"


# Define helper functions ------------------------------------------------------
plot_distribution <- function(cohorts, out_dir) {
  
  message("Plotting sample distribution")
  
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
  
  # Extract cohort locations
  cohort_locations <- dplyr::bind_rows(
    lapply(cohorts, function(cohort) {
      as.data.frame(colData(cohort))
    })) |>
    dplyr::count(country, name = "n_samples") |>
    mutate(
      continent = countrycode(
        country,
        origin = "iso3c",
        destination = "continent"))
  
  # Join cohort information
  country_points <- country_points |>
    left_join(
      cohort_locations,
      by = c("ISO3_CODE" = "country")) |>
    filter(!is.na(n_samples))
  
  # Summarise continent info
  continent_summary <- cohort_locations %>%
    group_by(continent) %>%
    summarise(
      samples = sum(n_samples),
      .groups = "drop"
    ) %>%
    mutate(
      percentage = samples / sum(samples) * 100)
  
  # Summarise country info
  country_summary <- cohort_locations %>%
    arrange(desc(n_samples))
  
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
    scale_size_continuous(
      breaks = c(500, 1000, 1500, 2000)) +
    theme(
      legend.text = element_text(size = 16),
      legend.title = element_blank(),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(1.2, "cm"),
      legend.position = "bottom",
      legend.direction = "horizontal")
  
  ggsave(
    filename = file.path(out_dir, "cohort_world_map.pdf"),
    plot = p,
    width = 10,
    height = 6)
  
  return(list(
    continents = continent_summary,
    countries = country_summary))
}


plot_composition <- function(cohorts, out_dir) {
  
  message("Plotting sample size composiion")
  
  case_control <- lapply(names(cohorts), function(x){
    meta <- as.data.frame(colData(cohorts[[x]]))
    meta$cohort <- x
    meta$status <- ifelse(
      meta$study_condition == "control",
      "control", "case")
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
  
  cohort_summary <- case_control %>%
    dplyr::count(cohort, status) %>%
    tidyr::pivot_wider(
      names_from = status,
      values_from = n,
      values_fill = 0
    ) %>%
    mutate(
      total = case + control,
      case_prop = case / total,
      control_prop = control / total,
      ratio = paste0(case, ":", control),
      imbalance = case_prop < 0.25 | case_prop > 0.75
    )
  
  star_df <- cohort_summary %>%
    filter(imbalance) %>%
    mutate(
      y = total + 30,
      label = "*"
    )
  
  p <- ggplot(case_control,
         aes(cohort, fill = status)) +
    geom_bar(position = "stack") +
    labs(
      y = "Count", x = "",
      fill = "Condition") +
    scale_fill_manual(
      values = c("control" = "#c5cb93ff", "case" = "#DFA57E"),
      labels = c(control = "Control", case = "Case")) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    geom_text(
      data = star_df,
      aes(
        x = cohort,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 7
    ) + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 18),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal")
  
  ggsave(
    filename = file.path(out_dir, "cohort_sample_size.pdf"),
    plot = p,
    width = 12,
    height = 4.5)
  
  return(cohort_summary)
}


# Define main function ---------------------------------------------------------
run_exploratory <- function(cohorts, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  geoplot <- plot_distribution(cohorts, out_dir)
  barplot <- plot_composition(cohorts, out_dir)
  
  message("\nPlots saved to: ", out_dir)
  
  return(list(
    geoplot = geoplot,
    barplot = barplot))
}


# Load data --------------------------------------------------------------------
cohorts <- readRDS("data/primary.rds") 


# Execute ----------------------------------------------------------------------
exploratory_res <- run_exploratory(cohorts, out_dir)








