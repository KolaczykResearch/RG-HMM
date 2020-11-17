#!/bin/bash -l

#$ -P sand               # Set SCC project
#$ -l h_rt=11:59:00      # Specify the hard time limit for the job
#$ -j y                  # Merge the error and output streams into a single file
#$ -pe omp 4
#$ -o Simulations/Testing/N_20

# Set parameters: Seed, N, M, p, q, gamma, alpha, beta, rate_obs, process, B, H

#B=50000
# Give job a name
#$ -N N20M12_1.5_ER_BX1
# #$ -N N20M12_1.5_PR_BX1
# #$ -N N20M20_1.5_ER_BX1
# #$ -N N20M20_1.5_PR_BX1
# #$ -N N20M28_1.5_ER_BX1
# #$ -N N20M28_1.5_PR_BX1
# #$ -N N20M36_1.5_ER_BX1
# #$ -N N20M36_1.5_PR_BX1
N=20 M=12 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=50000 H=4
# N=20 M=12 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=50000 H=4
# N=20 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=50000 H=4
# N=20 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=50000 H=4
# N=20 M=28 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=50000 H=4
# N=20 M=28 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=50000 H=4
# N=20 M=36 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=50000 H=4
# N=20 M=36 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=50000 H=4

#B=250000
# Give job a name
# #$ -N N20M12_1.5_ER_BX5
# #$ -N N20M12_1.5_PR_BX5
# #$ -N N20M20_1.5_ER_BX5
# #$ -N N20M20_1.5_PR_BX5
# #$ -N N20M28_1.5_ER_BX5
# #$ -N N20M28_1.5_PR_BX5
# #$ -N N20M36_1.5_ER_BX5
# #$ -N N20M36_1.5_PR_BX5
# N=20 M=12 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=250000 H=4
# N=20 M=12 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=250000 H=4
# N=20 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=250000 H=4
# N=20 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=250000 H=4
# N=20 M=28 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=250000 H=4
# N=20 M=28 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=250000 H=4
# N=20 M=36 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=250000 H=4
# N=20 M=36 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=250000 H=4


#B=500000
# Give job a name
# #$ -N N20M12_1.5_ER_BX1
# #$ -N N20M12_1.5_PR_BX10
# #$ -N N20M20_1.5_ER_BX10
# #$ -N N20M20_1.5_PR_BX10
# #$ -N N20M28_1.5_ER_BX10
# #$ -N N20M28_1.5_PR_BX10
# #$ -N N20M36_1.5_ER_BX10
# #$ -N N20M36_1.5_PR_BX10
# N=20 M=12 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=500000 H=4
# N=20 M=12 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=500000 H=4
# N=20 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=500000 H=4
# N=20 M=20 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=500000 H=4
# N=20 M=28 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=500000 H=4
# N=20 M=28 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=500000 H=4
# N=20 M=36 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='ER' B=500000 H=4
# N=20 M=36 p=0.9 q=0.1 gamma=2 alpha=0.01 beta=0.01 rate_obs=1.5 process='PR' B=500000 H=4


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




