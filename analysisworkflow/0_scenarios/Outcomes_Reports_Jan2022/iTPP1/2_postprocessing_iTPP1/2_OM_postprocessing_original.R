################################
################################
###                          ###
### STEP 2: POST-PROCESSING  ###
###                          ###
################################
################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Main script for running post-processing of OM simulations 
### to aggregate data which will be used to train GP
### 
### Original script:
### Created 12.02.2021
### lydia.burgert@unibas.ch 
###
### Adapted script:
### Saved 01.09.2021
### narimane.nekkab@unibas.ch
###
### R version 3.6.0
###
### -------------------------------------------------------------------------


##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# Library
library(dplyr)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/2_postprocessing/genOMpostprocscripts.R"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp ="..."

# Time
QOS = "1day"


##########################
### RUN POSTPROCESSING ###
##########################

# Run
genOMpostprocscripts(exp, QOS)


