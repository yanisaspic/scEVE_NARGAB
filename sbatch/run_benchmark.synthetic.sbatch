#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --array=0-1199%10
#SBATCH --time=24:00:00
#SBATCH --exclude=cpu-node-17,cpu-node-91

module load r/4.3.1

DATASETS=({1..1200})

dataset=${DATASETS[$SLURM_ARRAY_TASK_ID]}
mkdir ./results/benchmarks/${dataset}

echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SCRIPT:" run_benchmark.synthetic
echo "DATASET:" $dataset
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

Rscript ./R/benchmark/run_benchmark.base_methods.R $dataset
Rscript ./R/benchmark/run_benchmark.scEVE.R $dataset