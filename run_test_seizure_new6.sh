#!/bin/bash -l

#$ -P sand               # Set SCC project
#$ -l h_rt=240:00:00     # Specify the hard time limit for the job
# #$ -l h_rt=11:59:00      
#$ -j y                  # Merge the error and output streams into a single file
#$ -pe omp 4
#$ -o ApplicationSeizure/Results

----------------------------------------------------------------------------------------------------
# ----Seizure 1-------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
#### init (0.9, 0.1, 10, 0.01, 0.01)

###### AUTO DETECT: AUTO SEGMENT FINDER
# --------------------------------------
# #$ -N P1S1_BX10_1_0.5_99_123_left
#$ -N P1S1_BX10_1_0.5_178_209_left

# B=500000 H=4 data_path='ApplicationSeizure/Data/nets_EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left.mat' start=99 end=123
B=500000 H=4 data_path='ApplicationSeizure/Data/nets_EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left.mat' start=178 end=209


----------------------------------------------------------------------------------------------------
# ----Seizure 2-------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
###### AUTO DETECT: AUTO SEGMENT FINDER
# --------------------------------------
# #$ -N P1S2_BX10_1_0.5_237_255_left

# B=500000 H=4 data_path='ApplicationSeizure/Data/nets_EP001_clip_2_ecog_ws_1_wstp_0.5_sub_hemi_left.mat' start=237 end=255

# ----------------------------------------------------------------------------------------------------
# -seizure3----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
#### init (0.9, 0.1, 10, 0.01, 0.01)

###### AUTO DETECT: AUTO SEGMENT FINDER
# --------------------------------------
# #$ -N P1S3_BX10_1_0.5_98_164_left
# #$ -N P1S3_BX10_1_0.5_183_193_left

# B=500000 H=4 data_path='ApplicationSeizure/Data/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left.mat' start=98 end=164
# B=500000 H=4 data_path='ApplicationSeizure/Data/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left.mat' start=183 end=193

echo "==============================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "==============================================="
RANDOM=12111022520

module load R/3.5.1
Rscript main_test_seizure.R $(($RANDOM+$SGE_TASK_ID)) $B $H $data_path $start $end $SGE_TASK_ID

echo "===============================================" 
echo "Finished on : $(date)"
echo "==============================================="




