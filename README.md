# Microbiome-and-Disease
### _Description:_
This pipeline is built to train and cross-validate microbiome-based predictive models across independent stool cohorts from the `curatedMetagenomicData` R/Bioconductor repository. The primary aim is to distinguish disease-specific signatures from shared dysbiosis signals, which have direct implications for the use of microbiome-based classifiers in clinical diagnosis and precision medicine. 

The main objectives are:
  1. Explore sample distribution and composition of primary cohorts
  2. Characterise differences in microbial community structure and genus-level composition
  4. Build machine learning models and internally evaluate performance
  5. Externally validate trained models across cohorts and conditions

### _Workflow_:
The following scripts were made to be run sequentially in R.

#### `R/utils.R`
Defines utility functions used throughout the entire pipeline.

#### `R/cohorts.R`
Builds cohort objects from lists of primary and secondary cohorts in `data/` and saves as RDS.

#### `scripts/00_exploratory.R`
Generates exploratory figures to show sample distribution and composition.

#### `scripts/01_community.R`
Analyses differences in microbial community structure using Bray-Curtis dissimilarity.

#### `scripts/02_taxa.R`
Identifies differential taxa (genera) associated with cases using LEfSe.

#### `scripts/03_processing.R`
Normalises, harmonises and QCs datasets to prepare for training.

#### `scripts/04_ml.R`
Trains and externally validates a single processed model.

This script is executed by a Slurm wrapper found at `ml_array.sh`.

#### `scripts/05_analysis.R`
Analyses training and validation results from machine learning.

#### `scripts/06_figures.R`
Generates figures to help visualise findings.



