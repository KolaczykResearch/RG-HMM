#!/bin/bash -l

#$ -P sand               # Set SCC project
#$ -l h_rt=168:00:00   # Specify the hard time limit for the job
# #$ -l h_rt=11:59:00
#$ -j y                  # Merge the error and output streams into a single file
#$ -pe omp 4
#$ -o Simulations/Estimation/EXP3_vary_NM

# Set parameters: Seed, N, M, p, q, gamma, alpha, beta, rate_obs, process, process_est, B, H
# Fix p=0.7, p=0.3, gamma=2, alpha=0.03, beta=0.01, rate_obs=0.6, process=ER, process_est=ER, B=50000, H=10
# Vary N = 30, 50, 70; M = 50, 100, 150, 200

# N = 30
# Give job a name
#$ -N N30M50_0.6_ER_BX1
# #$ -N N30M100_0.6_ER_BX1
# #$ -N N30M150_0.6_ER_BX1
# #$ -N N30M200_0.6_ER_BX1
N=30 M=50 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=30 M=100 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=30 M=150 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=30 M=200 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10

# N = 50
# Give job a name
# #$ -N N50M50_0.6_ER_BX1
# #$ -N N50M100_0.6_ER_BX1
# #$ -N N50M150_0.6_ER_BX1
# #$ -N N50M200_0.6_ER_BX1
# N=50 M=50 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=50 M=100 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=50 M=150 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=50 M=200 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10

# N = 70
# Give job a name
# #$ -N N70M50_0.6_ER_BX1
# #$ -N N70M100_0.6_ER_BX1
# #$ -N N70M150_0.6_ER_BX1
# #$ -N N70M200_0.6_ER_BX1
# N=70 M=50 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=70 M=100 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=70 M=150 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10
# N=70 M=200 p=0.7 q=0.3 gamma=2 alpha=0.03 beta=0.01 rate_obs=0.6 process='ER' process_est='ER' B=50000 H=10

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




