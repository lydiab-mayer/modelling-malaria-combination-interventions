#!/bin/bash
#SBATCH --mem=1G
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1
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
ml R/4.1.0-foss-2018b

SIM_FOLDER=$1
DATE=$2
YEAR_BASELINE=$3
YEAR_INTERVENTIONA=$4
YEAR_INTERVENTIONB=$5
MIN_INT=$6

PARAM_TAB=$SIM_FOLDER"param_tab.txt"
PARAM_CAT=$SIM_FOLDER"param_ranges_cat.RData"
OM_FOLDER=$SIM_FOLDER"om/"
OUTPUT_FOLDER=$SIM_FOLDER"postprocessing/"
ERROR_FOLDER=$OUTPUT_FOLDER"err/"
SPLIT_FOLDER=$OUTPUT_FOLDER"split/"

# create the necessary folders
mkdir -p $OUTPUT_FOLDER
mkdir -p $SPLIT_FOLDER
mkdir -p $ERROR_FOLDER

# split the parameter table by setting
Rscript ../../../analysisworkflow/2_postprocessing/split_param_table.R $PARAM_TAB $SPLIT_FOLDER $PARAM_CAT

echo $PARAM_TAB
echo $SPLIT_FOLDER
echo $PARAM_CAT

# Submit postprocessing array job
split_files=(${SPLIT_FOLDER}*.txt)
NUM=${#split_files[@]}
sbatch -W --array=1-$NUM job_postprocessing.sh $SPLIT_FOLDER $OM_FOLDER $DATE $YEAR_BASELINE $YEAR_INTERVENTIONA $YEAR_INTERVENTIONB $MIN_INT
