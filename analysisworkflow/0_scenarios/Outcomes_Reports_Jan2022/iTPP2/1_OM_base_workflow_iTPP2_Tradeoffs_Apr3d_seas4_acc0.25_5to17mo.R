##############################
# Main script for specifying parameter values and running OpenMalaria simulations on the cluster. 
# 
#
# created 12.02.2021
# lydia.burgert@unibas.ch 
#
# updated 12.15.2021
# josephine.malinga@unibas.ch 
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

# Insert experiment name here
exp ="iTPP2_Tradeoffs_Apr3d_seas4_acc0.10_5to17mo"

# Create folder in working directory
create_folders(exp) # <----- run 1st time then comment

######################################
# Specify the desired parameter values 
######################################

chunk_size = 90000

QOS = "1week"

# Seasonality type: monthly EIR vs. Fourier transformation (choose 1)
Seasonality =read.table(paste0("./Experiments/",exp,"/seasonality_469_months_updatedDec15.txt"), sep="\t", header = TRUE)

#############################
# Categorical variables
#############################

# Biting patterns 
Biting_pattern <- data.frame(Biting_pattern=c("ML"),
                             indoor=c(0.6),outdoor=c(0.4))

# EIR
EIR= data.frame(EIR=c(0.035,0.50,1.25,2.50,3.50,5.50,10.0))

# age minimum and maximum
MaxAge = data.frame(MaxAge=c(1.5833),maxGroup=c(8))

# Intervention decay
vdecay <- data.frame(fundecay=c("weibull"),kdecay=c(0.69),vdecay=c("exp" ) )

# Coverage of healthcare system
# Access = data.frame(Access=c(0.10, 0.25, 0.50))  
Access = data.frame(Access=c(0.04, 0.0972, 0.2415))  

# Combine
param_cat = list(Seasonality=Seasonality,
                 Biting_pattern=Biting_pattern,
                 EIR=EIR,
                 MaxAge=MaxAge,
                 vdecay=vdecay,
                 Access=Access)

#############################
# Continuous variables and their ranges

# Name of the experiment and parameters
param_ranges_cont = rbind( Coverage = c(0.3, 1),
                           Halflife =c(0.5, 1),
                           Efficacy= c(0.5, 1) )

#param_ranges_cont = param_ranges_cont[,-2]
###############################
# Generate the parameter tables  
###############################

# no. of continuous parameters sampled via lhs
noSamples = 1000

# no. of OM seeds per sample
noSeeds= 5

# Generate
#gen_paramtable(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size)
gen_paramtable(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size)

########################################################
# Generate submission scripts and run the OM simulations
########################################################\

# Number of outputs to get in OM folder
2*nrow(Seasonality)*nrow(Biting_pattern)*nrow(EIR)*nrow(vdecay)*nrow(MaxAge)*nrow(Access)*noSamples*noSeeds

# Run
genOMsimscripts(exp, chunk_size)



