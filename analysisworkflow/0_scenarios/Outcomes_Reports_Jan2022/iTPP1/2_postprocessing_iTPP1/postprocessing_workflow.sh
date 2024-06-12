#!/bin/bash
#
##############################
# Main script for postprocessing OpenMalaria simulations on the cluster. 
# INPUT:
#       SIM_FOLDER = folder with simulation results
#       PARAM_TAB = file with the simulation parameters
#       OUTPUT_FOLDER = folder with postprocessing results
#       FOLLOW_UP = integer representing the survey index to consider for 
#                   evaluating intervention impact
# OUTPUT:
#       The script creates the following folders in the specified OUTPUT_FOLDER:
#               split/ = folder containing scenario parameters for each setting
#               processed/ = folder containing processing results
#
# SYNTHAX: 
#       bash postprocessing_workflow.sh SIM_FOLDER FOLLOW_UP
# 
#
# created 14.09.2019
# monica.golumbeanu@unibas.ch
#############################
ml purge
ml R/3.6.0-foss-2018b

SIM_FOLDER=$1
# FOLLOW_UP=$2
# YEARSBEFINT=$3
# INTYEAR=$4

PARAM_CAT=$SIM_FOLDER"param_ranges_cat.RData"
OM_FOLDER=$SIM_FOLDER"om/"
# ERROR_FOLDER1=$OUTPUT_FOLDER1"err/"

# PARAM_TAB1=$SIM_FOLDER"param_tab.txt"
# OUTPUT_FOLDER1=$SIM_FOLDER"postprocessing_baseline/"
# SPLIT_FOLDER1=$OUTPUT_FOLDER1"split/"
# ERROR_FOLDER1=$OUTPUT_FOLDER1"err/"

# PARAM_TAB2=$SIM_FOLDER"param_tab_prev.txt"
# OUTPUT_FOLDER2=$SIM_FOLDER"postprocessing_outcomes/"
# SPLIT_FOLDER2=$OUTPUT_FOLDER2"split/"
# ERROR_FOLDER2=$OUTPUT_FOLDER2"err/"

PARAM_TAB2=$SIM_FOLDER"param_tab.txt"
OUTPUT_FOLDER2=$SIM_FOLDER"postprocessing_outcomes/"
SPLIT_FOLDER2=$OUTPUT_FOLDER2"split/"
ERROR_FOLDER2=$OUTPUT_FOLDER2"err/"

# create the necessary folders
# mkdir -p $OUTPUT_FOLDER1
# mkdir -p $OUTPUT_FOLDER2
# mkdir -p $SPLIT_FOLDER1
# mkdir -p $SPLIT_FOLDER2
# mkdir -p $ERROR_FOLDER1
mkdir -p $ERROR_FOLDER2

# split the parameter table by setting
# Rscript ../../../analysisworkflow/2_postprocessing/split_param_baseline.R $PARAM_TAB1 $SPLIT_FOLDER1 $PARAM_CAT
# Rscript ../../../analysisworkflow/2_postprocessing/split_param_outcomes.R $PARAM_TAB2 $SPLIT_FOLDER2 $PARAM_CAT

# echo  $SPLIT_FOLDER1
echo  $SPLIT_FOLDER2

# # Submit postprocessing array job
# split_files=(${SPLIT_FOLDER1}*.txt)
# NUM1=${#split_files[@]}
# sbatch -W --array=1-$NUM1 job_postprocessing_baseline.sh $SPLIT_FOLDER1 $OM_FOLDER $OUTPUT_FOLDER1

# Submit postprocessing array job
split_files=(${SPLIT_FOLDER2}*.txt)
NUM2=${#split_files[@]}
sbatch -W --array=1-$NUM2 job_postprocessing_outcomes.sh $SPLIT_FOLDER2 $OM_FOLDER $OUTPUT_FOLDER2
