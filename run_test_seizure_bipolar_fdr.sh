#!/bin/bash -l

#$ -P sand               # Set SCC project
#$ -l h_rt=360:00:00     # Specify the hard time limit for the job
# #$ -l h_rt=11:59:00      
#$ -j y                  # Merge the error and output streams into a single file
#$ -pe omp 4
#$ -o ApplicationSeizure/results

# ----------------------------------------------------------------------------------------------------
# ----Seizure 1-------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
#$ -N P1S1_BX10_1_0.5_126_134_left
# #$ -N P1S1_BX10_1_0.5_453_463_left
# #$ -N P1S1_BX10_1_0.5_609_624_left

B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=126 end=134
# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=453 end=463
# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_1_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=609 end=624

# ----------------------------------------------------------------------------------------------------
# ----Seizure 2-------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# #$ -N P1S2_BX10_1_0.5_236_253_left
# #$ -N P1S2_BX10_1_0.5_295_318_left
# #$ -N P1S2_BX10_1_0.5_627_645_left

# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_2_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=236 end=253
# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_2_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=295 end=318
# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_2_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=627 end=645

# ----------------------------------------------------------------------------------------------------
# -Seizure3--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# #$ -N P1S3_BX10_1_0.5_229_238_left
# #$ -N P1S3_BX10_1_0.5_256_269_left
# #$ -N P1S3_BX10_1_0.5_309_316_left
# #$ -N P1S3_BX10_1_0.5_385_394_left

# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=229 end=238
# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=256 end=269
# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=309 end=316
# B=500000 H=4 data_path='DataSeizure/networks_fdr/nets_EP001_clip_3_ecog_ws_1_wstp_0.5_sub_hemi_left_bipolar_fdr_0.5.mat' start=385 end=394


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




