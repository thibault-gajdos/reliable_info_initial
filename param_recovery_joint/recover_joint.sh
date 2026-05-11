#!/bin/bash
#SBATCH --job-name=param_recovery_joint
#SBATCH --partition=Calcul
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=./output/joint_%a.log
#SBATCH --error=./error/joint_%a.err
#SBATCH --mail-user=thibault.gajdos-preuss@univ-amu.fr
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-768
#SBATCH --time=12:00:00

echo "============================================"
echo "SLURM Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================"

Rscript ./R/recover_joint.r ${SLURM_ARRAY_TASK_ID}

EXIT_CODE=$?
echo "============================================"
echo "End: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "============================================"
exit ${EXIT_CODE}

