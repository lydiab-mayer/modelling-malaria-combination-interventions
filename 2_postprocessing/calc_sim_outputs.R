################################
### STEP 2: POST-PROCESSING  ###
################################

# -------------------------------------------------------------------------------------------------------------
#
# Support script for running post-processing of OM simulations to aggregate data which will be used to train GP
# 
# Original script:
# Created 15.10.2019
# lydia.braunack-mayer@swisstph.ch 
#
# Adapted from monica.golumbeanu@unibas.ch
#
# R version 3.6.0
#
# -------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------
# SET UP
# -------------------------------------------------------------------------------------------------------------

# Define user
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Source helper functions
source(paste0("/scicore/home/penny/", user, "/M3TPP/analysisworkflow/2_postprocessing/postprocessing_resources.R"))

# Read in command arguments
args <- commandArgs(TRUE)
dir <- args[1]
split_file <- args[2]
date <- args[3]
year_baseline <- as.numeric(args[4])
year_interventionA <- as.numeric(args[5])
year_interventionB <- as.numeric(args[6])
min_int <- as.numeric(args[7])

# # Sample command arguments, retained here for testing
# dir <- "/scicore/home/penny/GROUP/M3TPP/obj6_test/om"
# split_file <- "/scicore/home/penny/GROUP/M3TPP/obj6_test/postprocessing/split/obj6test_seas4mo_Mali_16_5_weibull_0.1227.txt"
# date <- "2030-01-01"
# year_baseline <- 2034
# year_interventionA <- 2039
# year_interventionB <- 2044
# min_int <- 0.25

cat("Command arguments:")
print(paste("dir:", dir))
print(paste("Split file:", split_file))
print(paste("Date of first monitoring:", date))
print(paste("Year for baseline outcomes:", year_baseline))
print(paste("Year for intervention A outcomes:", year_interventionA))
print(paste("Year for intervention B outcomes:", year_interventionB))
print(paste("Minimum intervention age:", min_int))


# -------------------------------------------------------------------------------------------------------------
# PERFORM POST-PROCESSING
# -------------------------------------------------------------------------------------------------------------

# Postprocess the OpenMalaria simulations
cat("Run OM postprocessing function:")
postprocess.om(dir = dir, 
               param.file = split_file,
               date = date,
               year.baseline = year_baseline,
               year.interventionA = year_interventionA,
               year.interventionB = year_interventionB,
               min.int = min_int)
