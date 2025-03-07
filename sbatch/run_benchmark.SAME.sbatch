#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --array=0-14%1
#SBATCH --time=24:00:00
#SBATCH --exclude=cpu-node-17,cpu-node-91

module load r/4.3.1

DATASETS=("Li_HumCRC_a" \
"Li_HumCRC_b" \
"Baron_HumPan_1" \
"Baron_HumPan_2" \
"Baron_HumPan_3" \
"Baron_HumPan_4" \
"Baron_MouPan_1" \
"Baron_MouPan_2" \
"Darmanis_HumGBM" \
"Tasic_MouBra" \
"JerbyArnon_HumMLM" \
"Gillen_HumEPN" \
"VanGalen_HumAML" \
"Lambrechts_HumNSCLC" \
"Peng_HumPDAC")

dataset=${DATASETS[$SLURM_ARRAY_TASK_ID]}
mkdir ./results/benchmarks/${dataset}

echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SCRIPT:" run_benchmark.SAME
echo "DATASET:" $dataset
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

Rscript ./R/benchmark/run_benchmark.SAME.R $dataset