# TODO 33 x 33 matrix, then evaluate and label category

run_evaluation <- function() {
  ValidationData <- RunScreening(
    # Input data
    obj = CrcBiomeScreenObject,
    # Model and data splitting
    model_type = "RF",
    partition = 0.7,
    split.requirement = list(
      label = c("control", "CRC"),
      condition_col = "study_condition"
    ),
    ClassWeights = TRUE,
    
    # Cross-validation and parallelization
    n_cv = 10,
    num_cores = 10,
    
    # Task and output naming
    TaskName = "RF_GMPR_toydata",
    
    # External validation
    ValidationData = ValidationData,
    TrueLabel = "CRC"
  )
}