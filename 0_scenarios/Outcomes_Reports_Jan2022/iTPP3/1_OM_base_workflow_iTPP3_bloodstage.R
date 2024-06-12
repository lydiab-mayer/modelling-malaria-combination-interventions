##############################
# Main script for specifying parameter values and running OpenMalaria simulations on the cluster. 
# 
#
# created 12.02.2021
# lydia.burgert@unibas.ch 
#############################

# Setup
rm(list = ls())
library(tgp)
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/1_OM_basic_workflow/genOMsimscripts.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/generate_param_table.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/create_folders.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/convert_access.R"))

# Insert experiment name here
exp ="iTPP3_bloodstage_OUTCOMESREPORT"

# Create folder in working directory
#create_folders(exp) # <----- run 1st time then comment

# Seasonality type: monthly EIR vs. Fourier transformation (choose 1)
#seasonality_type = "monthly"
seasonality_type = "Fourier"

# !!! Depending on which seasonality is used, choose the correct scaffold.xml to upload
# !!! If using Fourier, make sure number of coefficients correspond (example here I use a 6 coeff model with requires 5 a and 5 b values)

# !!! Copy the scaffold.xml file for your experiment into the ./M3TPP/Experiments/"exp"/OM_jobs folder!!!
# !!! Copy the same scaffold.xml file for your experiment into the .GROUP//M3TPP/Experiments/"exp"/OM_jobs folder!!!
# !!! Copy a seasonality.csv file for your experiment into the ./M3TPP/Experiments/"exp"/ folder!!!

######################################
# Specify the desired parameter values 
######################################

chunk_size = 120000

QOS = "30min"

#############################
# Categorical variables

# Seasonality
if(seasonality_type == "monthly"){
  Seasonality =read.table(paste0("./Experiments/",exp,"/seasonality_monthly.txt"), sep="\t", header = TRUE)
}
if(seasonality_type == "Fourier"){
  Seasonality =read.table(paste0("./Experiments/",exp,"/seasonality_centered_June_Fourier.txt"), header = TRUE)
}

# Biting patterns 
Biting_pattern <- data.frame(Biting_pattern=c("Mali"), indoor=c(0.6), outdoor=c(0.4))

# EIR
EIR = data.frame(EIR = c(2, 4, 8, 16, 32, 64, 128, 256)) #c(1, 2, 4, 8, 16, 32, 64, 128)

# Max age intervention
# age_groups (surveys age group)
# 0.25 = group 1
# 0.5 = group 2
# 0.75 = group 3
# 1 = group 4
# 2 = group 5
# 3 = group 6
# 4 = group 7
# 5 = group 8
# 6 = group 9
# 7 = group 10
# 8 = group 11
# 9 = group 12
# 10 = group 13
# 15 = group 14
# 20 = group 15
# 100 = group 16

MaxAge = data.frame(MaxAge = c(5), maxGroup = c(8)) # age groups specified as per iTPP3 document

# Coverage of healthcare system: probability of healthcare seeking of uncomplicated malaria cases (in %) 
initial_access = data.frame(Access = c(0.1, 0.5)) # 10% = low access, 50% = high access

# Convert access to care to 5 day probabilities for use in XML files
# Sources code from MMC project
Access = round(pmax(convert_access(initial_access * 100), 0.04), 2)

# Timing of cohort recruitment and intervention deployment
Timing = data.frame(Timing = "May",
                    Cohort = "04-25", 
                    Round1 = "05-01", 
                    Round2 = "06-01", 
                    Round3 = "07-01", 
                    Round4 = "08-01")

# Intervention PD parameters
IC50 = data.frame(IC50 = c(0.020831339))

# Combine
param_cat = list(Seasonality = Seasonality,
                 Biting_pattern = Biting_pattern,
                 EIR = EIR,
                 MaxAge = MaxAge,
                 Access = Access,
                 Timing = Timing,
                 IC50 = IC50)

#############################
# Continuous variables and their ranges

# Name of the experiment and parameters
param_ranges_cont = rbind(Coverage1 = c(0.7, 1),
                          Coverage2 = c(0.7, 1),
                          Halflife = c(5, 40),
                          MaxKillingRate = c(2, 30),
                          Slope = c(1, 8))

###############################
# Generate the parameter tables  
###############################

# no. of continuous parameters sampled via lhs
noSamples = 1000 # Increase to 1000 for final run

# no. of OM seeds per sample
noSeeds = 15 # Increase to 15 for final run

# Generate
gen_paramtable(exp, param_ranges_cont, param_cat, noSamples, noSeeds,chunk_size)

########################################################
# Generate submission scripts and run the OM simulations
########################################################\

# Number of outputs to get in OM folder
2*nrow(Seasonality)*nrow(Biting_pattern)*nrow(EIR)*nrow(MaxAge)*nrow(Access)*nrow(IC50)*noSamples*noSeeds

# Run
genOMsimscripts(exp, chunk_size)


