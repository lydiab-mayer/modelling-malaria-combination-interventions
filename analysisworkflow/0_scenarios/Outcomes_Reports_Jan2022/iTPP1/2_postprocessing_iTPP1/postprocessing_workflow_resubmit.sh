#!/bin/bash
#
##############################
# Script for resubmission of postprocessing OpenMalaria simulations on the cluster. 
# created 10.12,2021
# narimane.nekkab@unibas.ch
#############################
ml purge
ml R/3.6.0-foss-2018b

SIM_FOLDER=$1
OM_FOLDER=$SIM_FOLDER"om/"
OUTPUT_FOLDER=$SIM_FOLDER"postprocessing_outcomes/"
SPLIT_FOLDER=$OUTPUT_FOLDER"split_resubmit2/"

# Submit postprocessing array job
split_files=(${SPLIT_FOLDER}*.txt)
NUM=${#split_files[@]}
sbatch -W --array=1-$NUM job_postprocessing_outcomes_resubmit.sh $SPLIT_FOLDER $OM_FOLDER $OUTPUT_FOLDER
