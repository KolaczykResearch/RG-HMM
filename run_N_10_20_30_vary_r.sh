#!/bin/bash -l

#$ -P sand               # Set SCC project
# #$ -l h_rt=168:00:00   # Specify the hard time limit for the job
#$ -l h_rt=11:59:00
#$ -j y                  # Merge the error and output streams into a single file
#$ -pe omp 4
#$ -o Simulations/Testing/N_10_20_30_vary_r

# Set parameters: Seed, N, M, p, q, gamma, alpha, beta, rate_obs, process, B, H
# fix scaled_time_obs = 1.34 (M = ...), B = 500,000, varying rate_obs

# N = 10
# Give job a name
# #$ -N N10M3_0.5_ER_BX10
# #$ -N N10M3_0.5_PR_BX10
# #$ -N N10M7_1_ER_BX10
# #$ -N N10M7_1_PR_BX10
# #$ -N N10M10_1.5_ER_BX10
# #$ -N N10M10_1.5_PR_BX10
# rate_obs = 0.5
#N=10 M=3 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=0.5 process='ER' B=500000 H=4
#N=10 M=3 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=0.5 process='PR' B=500000 H=4
# rate_obs = 1
#N=10 M=7 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1 process='ER' B=500000 H=4
#N=10 M=7 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1 process='PR' B=500000 H=4
# rate_obs = 1.5
#N=10 M=10 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=500000 H=4
#N=10 M=10 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=500000 H=4


# N = 20
# Give job a name
# #$ -N N20M7_0.5_ER_BX10
# #$ -N N20M7_0.5_PR_BX10
# #$ -N N20M13_1_ER_BX10
# #$ -N N20M13_1_PR_BX10
# #$ -N N20M20_1.5_ER_BX10
# #$ -N N20M20_1.5_PR_BX10
# rate_obs = 0.5
#N=20 M=7 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=0.5 process='ER' B=500000 H=4
#N=20 M=7 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=0.5 process='PR' B=500000 H=4
# rate_obs = 1
#N=20 M=13 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1 process='ER' B=500000 H=4
#N=20 M=13 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1 process='PR' B=500000 H=4
# rate_obs = 1.5
#N=20 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=500000 H=4
#N=20 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=500000 H=4


# N = 30
# Give job a name
# #$ -N N30M10_0.5_ER_BX10
# #$ -N N30M10_0.5_PR_BX10
# #$ -N N30M20_1_ER_BX10
# #$ -N N30M20_1_PR_BX10
# #$ -N N30M30_1.5_ER_BX10
# #$ -N N30M30_1.5_PR_BX10
# rate_obs = 0.5
#N=30 M=10 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=0.5 process='ER' B=500000 H=4
#N=30 M=10 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=0.5 process='PR' B=500000 H=4
# rate_obs = 1
#N=30 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1 process='ER' B=500000 H=4
#N=30 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1 process='PR' B=500000 H=4
# rate_obs = 1.5
#N=30 M=30 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=500000 H=4
#N=30 M=30 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=500000 H=4

echo "==============================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "==============================================="
RANDOM=12111022

module load R/3.5.1
Rscript main_test.R $(($RANDOM+$SGE_TASK_ID)) $N $M $p $q $gamma $alpha $beta $rate_obs $process $B $H $SGE_TASK_ID

echo "===============================================" 
echo "Finished on : $(date)"
echo "==============================================="




