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
exp ="iTPP1a_3mo119mo_1round_deploymentTS18"

# Create folder in working directory
# create_folders(exp) # <----- run 1st time then comment

# Optimize file size
chunk_size = 100000

QOS = "30min"

##################
### PARAMETERS ###
##################

# Specify the desired parameter values 
# Note some values are categorical and others continuous
# Values are sampled for continuous variables from defined range

##################
### DEMOGRAPHY ###

# Population size
pop_size = 10000

# Demography dataframe
demography=data.frame(demography="Tanzania",
                      pop_size=pop_size)

##################
### MONITORING ###

# Years of analysis
start_year = 2030
end_year   = 2046

# Number of years burn in
burn_in_years = 20

# Set the burn in relative to start year
burn_in = start_year - burn_in_years

# Monitoring dataframe
monitoring=data.frame(monitoring="15yObs",
                      start_year=start_year,
                      end_year=end_year,
                      burn_in=burn_in)

####################
### INTERVENTION ###

# Efficacy range (in %)
Efficacy= c(0.5, 1)

# Half-life range (in days)
Halflife =c(30, 300)

# Coverage range (in %) (single range, can have multiple)
Coverage = c(0.3, 1)

# Maximum age (refer to reference survey age group) --> default to be chosen soon
# <ageGroup lowerbound="0">
#   <group upperbound="0.165"/> <!-- 2 months -->
#   <group upperbound="0.25"/> <!-- 3 months -->
#   <group upperbound="0.5"/> <!-- 6 months -->
#   <group upperbound="1"/>
#   <group upperbound="2"/>
#   <group upperbound="3"/>
#   <group upperbound="4"/>
#   <group upperbound="5"/>
#   <group upperbound="6"/>
#   <group upperbound="7"/>
#   <group upperbound="8"/>
#   <group upperbound="9"/>
#   <group upperbound="10"/>
#   <group upperbound="15"/>
#   <group upperbound="20"/>
#   <group upperbound="50"/>
#   <group upperbound="100"/>
MaxAge = data.frame(MaxAge=c(10), 
                    maxGroup=c(13))

# Decay: define function, k and, efficacyB
# ---- function: "constant" or "step" or "linear" or "exponential" or "weibull" or "hill" or "smooth-compact"
# ---- k:         shape parameter of distribution. If not specified, default value of 1 is used.
# ---- efficacyB: measure of variation
decay = data.frame(decay="Weibull",
                   fundecay=c("weibull"),
                   effB=c(10))

# Continuous k shape parameter
kdecay = c(0.01,7)

#####################
### HEALTH SYSTEM ###

# # Coverage of healthcare system: probability of healthcare seeking of uncomplicated malaria cases (in %) 
# initial_access = data.frame(access=c(0.1,0.25,0.5)) 
# 
# # Convert access to care to 5 day probabilities for use in XML files
# # Sources code from MMC project
# access = round(pmax(convert_access(initial_access * 100), 0.04), digits = 4)
# 
# # Access dataframe
# Access = data.frame(Access="LowMedHighAccess", 
#                     access=access)

# Coverage of healthcare system: probability of healthcare seeking of uncomplicated malaria cases (in %) 
initial_access = data.frame(access=c(0.1,0.25,0.5)) 

# Convert access to care to 5 day probabilities for use in XML files
# Sources code from MMC project
access = round(pmax(convert_access(initial_access * 100), 0.04), digits = 4)

# Access dataframe
Access = data.frame(Access="LowModHighAccess", 
                    access=access)

##################
### ENTOMOLOGY ###

# Seasonality type: monthly EIR vs. Fourier transformation (by month or with coefficients)
# Check XML structure to match

seasonality = read.table(paste0("./Experiments/",exp,"/seasonality_34569_months_updatedNov26.txt"), sep="\t", header = TRUE)
########################## ONE SEASON ONLY #############################
seasonality = seasonality[which(seasonality$Seasonality == "seas4mo"),]

# Biting patterns 
biting_pattern <- data.frame(biting_pattern="EqualBites",
                             indoor=c(0.5),
                             outdoor=c(0.5))

########################
### UPLOAD EIR TABLE ###

# EIR
EIR_ref = read.table(paste0("./Experiments/",exp,"/REFERENCE_PARAM_TABLE_V11.csv"), sep = ",", header = TRUE)
EIR_ref$access = round(EIR_ref$access, digits = 4)
EIR_ref = EIR_ref[which(EIR_ref$access %in% access$access),]
EIR_ref = EIR_ref[which(EIR_ref$Seasonality %in% seasonality$Seasonality),]
EIR_ref = EIR_ref[which(EIR_ref$EIR >= 0.035),]
EIR_ref = EIR_ref[complete.cases(EIR_ref),]

# Unique EIR
EIR =  data.frame(EIR=unique(EIR_ref$EIR)[order(unique(EIR_ref$EIR))])

###################
### DIAGNOSTICS ###


#############
### MODEL ###


##############################
### GENERATE PARAMS TABLES ###
##############################

# Names of list corresponds to name of variable to include in split-file name

# Categorical variables
param_cat = list(demography=demography,
                 monitoring=monitoring,
                 MaxAge=MaxAge,
                 decay=decay,
                 Access=Access,
                 Seasonality=seasonality,
                 biting_pattern=biting_pattern,
                 EIR=EIR)

#  Continuous variables
param_ranges_cont = rbind(kdecay, Coverage, Halflife, Efficacy)

# Number of continuous parameters to sample via lhs
noSamples = 1000

# Number of OM seeds per sample
noSeeds=  10

# Generate
# gen_paramtable(exp, param_ranges_cont, param_cat, noSamples, noSeeds, chunk_size)

#### New table
gen_paramtable_optimizedCombos(exp, param_ranges_cont, param_cat, noSamples, noSeeds, chunk_size, EIR_ref)



###########################################################
### GENERATE SCENARIOS AND RUN OPEN MALARIA SIMULATIONS ###
###########################################################

# Number of outputs to get in OM folder
param_tab=read.table(paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab.txt"), header = T)
nrow(param_tab)

# Numer of sims
2*nrow(param_tab)

# Run
genOMsimscripts(exp, chunk_size)



