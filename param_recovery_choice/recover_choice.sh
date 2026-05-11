#!/bin/bash
#SBATCH --job-name=param_recovery_choice
#SBATCH --partition=Calcul
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=./output/choice_%a.log
#SBATCH --error=./error/choice_%a.err
#SBATCH --mail-user=thibault.gajdos-preuss@univ-amu.fr
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-64
#SBATCH --time=08:00:00

echo "============================================"
echo "SLURM Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================"

Rscript ./R/recover_choice.r ${SLURM_ARRAY_TASK_ID}

EXIT_CODE=$?
echo "============================================"
echo "End: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "============================================"
exit ${EXIT_CODE}
