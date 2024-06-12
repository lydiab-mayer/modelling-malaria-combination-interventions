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
exp ="iTPP1a_3mo119mo_1round_deploymentTS18_3cont"

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
#   <group upperbound="0.25"/>
#   <group upperbound="0.5"/>
#   <group upperbound="1"/>
#   <group upperbound="2"/>
#   <group upperbound="3"/>
#   <group upperbound="4"/>
#   <group upperbound="5"/>
#   <group upperbound="10"/>
#   <group upperbound="20"/>
#   <group upperbound="100"/>
MaxAge = data.frame(MaxAge=c(10), 
                    maxGroup=c(8))

# Decay: define function, k and, efficacyB
# ---- function: "constant" or "step" or "linear" or "exponential" or "weibull" or "hill" or "smooth-compact"
# ---- k:         shape parameter of distribution. If not specified, default value of 1 is used.
# ---- efficacyB: measure of variation
decay = data.frame(decay="3decays",
                   fundecay=c("weibull","weibull","weibull"),
                   kdecay=c(1,2,4),
                   effB=c(10,10,10))

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

# Load prevalence table
prevalence_data = read.csv("/scicore/home/penny/GROUP/M3TPP/iTPP1a_3mo59mo_1round_deploymentTS18/prevalence_data.csv")
prev_ref_table = read.table(paste0("/scicore/home/penny/GROUP/M3TPP/iTPP1a_3mo59mo_1round_deploymentTS18/prevalence_table_seed.txt"), header = TRUE)
table(prev_ref_table$Seasonality, prev_ref_table$prev_210_cat)

# Subset, remove prev. > 50%
prev_ref_sub = prev_ref_table %>% filter(!(prev_210_cat == "(0.5,0.6]" | prev_210_cat == "(0.6,0.7]"))

# Remove duplicates
prev_ref_sub = prev_ref_sub %>% filter(!(Seasonality == "seas3mo" & EIR == 0.035 |
                                           Seasonality == "seas3mo" & EIR == 0.075 & access != 0.2412 |
                                           Seasonality == "seas3mo" & EIR == 0.1 & access != 0.0972|
                                           Seasonality == "seas3mo" & EIR == 1.5 & access != 0.0400 |
                                           Seasonality == "seas3mo" & EIR == 2 & access != 0.0972 |
                                           Seasonality == "seas4mo" & EIR == 0.075 & access != 0.2412|
                                           Seasonality == "seas4mo" & EIR == 0.05 & access == 0.2412 |
                                           Seasonality == "seas4mo" & EIR == 0.1 & access == 0.2412 |
                                           Seasonality == "seas4mo" & EIR == 0.25 |
                                           Seasonality == "seas4mo" & EIR == 0.7 |
                                           Seasonality == "seas4mo" & EIR == 1 |
                                           Seasonality == "seas4mo" & EIR == 1.4 |
                                           Seasonality == "seas4mo" & EIR == 2 |
                                           Seasonality == "seas9mo" & EIR == 0.8 & access == 0.2412))

# Because of high number of simulations, need to drop seeds to 5
prev_ref_sub = prev_ref_sub %>%  filter(seed < 6)

# Tables
table(prev_ref_sub$Seasonality, prev_ref_sub$prev_210_cat)
table(prev_ref_sub$Seasonality, prev_ref_sub$prev_210_cat, prev_ref_sub$access)

# Get table of unique combinations
EIR_ref = prev_ref_sub %>% select(Seasonality, access, EIR) %>% distinct()
EIR_ref = EIR_ref[which(EIR_ref$Seasonality %in% seasonality$Seasonality),]
EIR_ref = EIR_ref[complete.cases(EIR_ref),]

# Save
EIR_ref = EIR_ref %>% arrange(EIR, access)
write.table(EIR_ref, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/EIR_reference_table_5seeds.csv"), row.names = F)
write.table(prev_ref_sub, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Prev210_reference_table_5seeds.csv"), row.names = F)

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
param_ranges_cont = rbind(Coverage, Halflife, Efficacy)

# Number of continuous parameters to sample via lhs
noSamples = 1000

# Number of OM seeds per sample
noSeeds=  5

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



