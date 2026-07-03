#!/bin/bash
#SBATCH --job-name=ml_train
#SBATCH --time=06:00:00
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --array=0-33
#SBATCH --error=slurm/jobs/%A_%a.err
#SBATCH --output=slurm/jobs/%A_%a.out

module load miniforge
conda activate r_env

BASE_DIR=/mnt/scratch/qvsg0202/microbiome-and-disease

mkdir -p "$BASE_DIR/slurm/jobs"

# Path to comparisons.txt
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BASE_DIR/results/04_pre_ml/comparisons.txt")

# Path to R script
Rscript "$BASE_DIR/scripts/05_ml.R" "$FILE"


