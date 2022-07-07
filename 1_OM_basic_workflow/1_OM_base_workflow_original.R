#####################################
#####################################
###                               ###
### STEP 1: OM SETUP & SIMULATION ###
###                               ###
#####################################
#####################################

### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Main script for specifying parameter values, generating scenario XMLs,
### and running OpenMalaria simulations on the cluster 
### 
### Original script:
### Created 12.02.2021
### lydia.burgert@unibas.ch 
###
### Adapted script:
### Saved 01.09.2021
### Updated 07.07.2022
### narimane.nekkab@unibas.ch
###
### R version 3.6.0
###
### -------------------------------------------------------------------------

closeAllConnections()

##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# Library
library(tgp)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/1_OM_basic_workflow/genOMsimscripts.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/generate_param_table.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/create_folders.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/convert_access.R"))

##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp = "..."

# Create folder in working directory
create_folders(exp) # <----- run 1st time then comment !!!

# Reference folders
GROUP = "/scicore/home/penny/GROUP/M3TPP/"
EXPERIMENT_FOLDER = paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/")

# !!! Copy the scaffold.xml file for your experiment into the ./M3TPP/Experiments/"exp"/OM_jobs folder!!!
# !!! Copy the same scaffold.xml file for your experiment into the .GROUP//M3TPP/Experiments/"exp"/OM_jobs folder!!!
# !!! Copy a seasonality.csv file for your experiment into the ./M3TPP/Experiments/"exp"/ folder!!!

######################################
# Specify the desired parameter values 
######################################

chunk_size = 90000

QOS = "..."

#############################
# Categorical variables

# Seasonality
# !!! Depending on which seasonality is used, choose the correct scaffold.xml to upload
# !!! If using Fourier, make sure number of coefficients correspond (example here I use a 6 coeff model with requires 5 a and 5 b values)
seasonality_file = "seasonality_monthly.txt" 
Seasonality = read.table(paste0(EXPERIMENT_FOLDER,exp,"/",seasonality_file), sep="\t", header = TRUE)

# Biting patterns 
Biting_pattern <- data.frame(Biting_pattern=c("Mali"),indoor=c(0.6),outdoor=c(0.4))

# EIR
EIR= data.frame(EIR=c(10))

# Max age intervention 
# !!! Refer to age group in scaffold.xml 
# e.g. 0.25 = group 1; 2 = group 2; 5 = group 3; 10 = group 4; 15 = group 5 etc.
MaxAge = data.frame(MaxAge=c(10),maxGroup=c(4))

# Intervention decay
LAIdecay <- data.frame(fundecay=c("weibull"),kdecay=c(1 ),LAIdecay=c("exp" ) )

# Coverage of healthcare system: probability of healthcare seeking of uncomplicated malaria cases (in %) 
# Access = data.frame(Access=c(0.1))  # ---> default
# ALTERNATIVE
# Convert access to care to 5 day probabilities for use in XML files
# Sources code from MMC project
initial_access = 0.1
Access = round(pmax(convert_access(initial_access * 100), 0.04), digits = 4)

##############################
### GENERATE PARAMS TABLES ###
##############################

# Combine categorical variables
param_cat = list(Seasonality=Seasonality,
                 Biting_pattern=Biting_pattern,
                 EIR=EIR,
                 MaxAge=MaxAge,
                 LAIdecay=LAIdecay,
                 Access=Access)

# Continuous variables and their ranges
param_ranges_cont = rbind( Coverage = c(0.4, 1),
                           Halflife =c(30,150),
                           Efficacy= c(0.7,1) )

# no. of continuous parameters sampled via lhs
noSamples = 5

# no. of OM seeds per sample
noSeeds=  2

# Generate
gen_paramtable(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size)

###########################################################
### GENERATE SCENARIOS AND RUN OPEN MALARIA SIMULATIONS ###
###########################################################

# Number of outputs to get in OM folder
2*nrow(Seasonality)*nrow(Biting_pattern)*nrow(EIR)*nrow(LAIdecay)*nrow(MaxAge)*nrow(Access)*noSamples*noSeeds

# Run
genOMsimscripts(exp, chunk_size)

##################
### SAVE FILES ###
##################

# Save script
rstudioapi::documentSave()

# Add experiment file, seasonal file, and scaffold to GROUP folder

file.copy(paste0(GROUP,exp,"/1_OM_base_workflow_",exp,".R"), 
          paste0(EXPERIMENT_FOLDER,exp,"/1_OM_base_workflow_",exp,".R"), overwrite=FALSE)

file.copy(paste0(GROUP,exp,"/scaffold_",exp,".R"), 
          paste0(EXPERIMENT_FOLDER,exp,"/OM_JOBS/scaffold.R"), overwrite=FALSE)

file.copy(paste0(GROUP,exp,"/",seasonality_file,"_",exp,".R"), 
          paste0(EXPERIMENT_FOLDER,exp,"/",seasonality_file), overwrite=FALSE)









