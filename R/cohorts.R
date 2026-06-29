# Construct cohorts objects from list using functions from utils.R

# Load packages and dependencies -----------------------------------------------
library(curatedMetagenomicData)
library(SummarizedExperiment)
source("R/utils.R")

test_cohorts <- read_cohort_list("data/test.txt") |>
  load_cohorts(cohorts = _)

saveRDS(test_cohorts, "data/test.rds")


# primary_cohorts <- read_cohort_list("data/primary.txt") |>
#   load_cohorts(cohorts = _)
#   
# saveRDS(primary_cohorts, "data/primary.rds")


