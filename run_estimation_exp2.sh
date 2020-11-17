#!/bin/bash -l

#$ -P sand               # Set SCC project
#$ -l h_rt=168:00:00   # Specify the hard time limit for the job
# #$ -l h_rt=11:59:00
#$ -j y                  # Merge the error and output streams into a single file
#$ -pe omp 4
#$ -o Simulations/Estimation/EXP2_vary_B

# Set parameters: Seed, N, M, p, q, gamma, alpha, beta, rate_obs, process, process_est, B, H
# Fix N=20, M=50, p=0.7, p=0.3, gamma=2, alpha=0.03, beta=0.01, rate_obs=0.6, process=ER, process_est=ER, H=10
# Vary B = 50,000; 250,000; 500,000; 1000,000

# Give job a name
#$ -N N20M50_BX1
# #$ -N N20M50_BX5
# #$ -N N20M50_BX10
# #$ -N N20M50_BX20

N=20 M=50 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=20 M=50 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=250000 H=10
# N=20 M=50 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=500000 H=10
# N=20 M=50 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=1000000 H=10

echo "==============================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "==============================================="
RANDOM=12111022

module load R/3.5.1
Rscript main_estimation.R $(($RANDOM+$SGE_TASK_ID)) $N $M $p $q $gamma $alpha $beta $rate_obs $process $process_est $B $H $SGE_TASK_ID

echo "===============================================" 
echo "Finished on : $(date)"
echo "==============================================="




