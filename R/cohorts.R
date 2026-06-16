# Construct cohorts objects from list using functions from utils.R

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(SummarizedExperiment)
source("R/utils.R")

# test_cohorts <- read_cohort_list("data/test.txt") |>
#   load_cohorts(cohorts = _) |>
#   create_cohort_objects()
# 
# saveRDS(test_cohorts, "data/test_cohorts.rds")


# primary_cohorts <- read_cohort_list("data/primary.txt") |>
#   load_cohorts(cohorts = _) |>
#   create_cohort_objects()
# 
# saveRDS(primary_cohorts, "data/primary_cohorts.rds")


# primary_a_cohorts <- read_cohort_list("data/primary_a.txt") |>
#   load_cohorts(cohorts = _) |>
#   create_cohort_objects()
# 
# saveRDS(primary_a_cohorts, "data/primary_a_cohorts.rds")


primary_b_cohorts <- read_cohort_list("data/primary_b.txt") |>
  load_cohorts(cohorts = _) |>
  create_cohort_objects()

saveRDS(primary_b_cohorts, "data/primary_b_cohorts.rds")
