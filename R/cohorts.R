# Construct objects from a list of cohort names using functions from utils.R

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(SummarizedExperiment)
source("R/utils.R")

# Construct primary cohorts ----------------------------------------------------
primary <- read_cohort_list("data/primary.txt") |>
  load_cohorts(cohorts = _)

saveRDS(primary, "data/primary.rds")

# Construct secondary cohorts --------------------------------------------------
secondary_case <- read_cohort_list("data/secondary_case.txt") |>
  load_cohorts(cohorts = _)

secondary_control <- read_cohort_list("data/secondary_control.txt") |>
  load_cohorts(cohorts = _)

saveRDS(primary_cohorts, "data/secondary_case.rds")
saveRDS(primary_cohorts, "data/secondary_control.rds")