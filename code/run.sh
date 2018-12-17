#!/bin/bash -l

#$ -P sand            # Set SCC project
# #$ -l h_rt=72:00:00   # Specify the hard time limit for the job
#$ -l h_rt=168:00:00
#$ -N run1           # Give job a name
#$ -j y               # Merge the error and output streams into a single file

# Request a whole processor node of 28 cores with at least 512 GB of RAM
# #$ -pe omp 28
# #$ -l mem_per_core=18G

#$ -pe omp 8

module load R/3.5.0
#R CMD BATCH main.R
Rscript main.R 
