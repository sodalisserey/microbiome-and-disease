#!/bin/bash
#SBATCH --job-name=ml
#SBATCH --time=10:00:00
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --array=0-33
#SBATCH --error=slurm/jobs/%A_%a.err
#SBATCH --output=slurm/jobs/%A_%a.out

# IN_DIR essentials before Slurm submission: 
#   /final/ (directory with all processed comparison objects)
#   /combinations.txt (paths to all comparison objects in final/)
#   /comparison_table.csv (matrix with all validation tasks)

module load miniforge
conda activate r_env

BASE_DIR=/mnt/scratch/qvsg0202/microbiome-and-disease
IN_DIR="$BASE_DIR/results/03_processing"

mkdir -p "$BASE_DIR/slurm/jobs"

# Path to combinations.txt
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IN_DIR/combinations.txt")

# Path to R script
Rscript "$BASE_DIR/scripts/04_ml.R" "$FILE"


