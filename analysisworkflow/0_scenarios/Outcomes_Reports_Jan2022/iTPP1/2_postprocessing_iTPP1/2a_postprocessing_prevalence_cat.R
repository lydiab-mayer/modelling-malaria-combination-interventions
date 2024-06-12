#################################
#################################
###                           ###
### STEP 2a: POST-PROCESSING  ###
###                           ###
#################################
#################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Add reference prevalence range to postprocessing
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

# Group folder
GROUP = "/scicore/home/penny/GROUP/M3TPP/"

# Insert experiment name here
exp = "..."


#################
### LOAD DATA ###
#################

# Get postprocessing aggregate files 
agg_names = list.files(path = paste0(GROUP,exp,"/postprocessing_baseline/"), pattern = "agg", full.names = FALSE)

# Load
agg_data_all <- do.call(rbind, lapply(paste0(GROUP,exp,"/postprocessing_baseline/",agg_names), read.table, header = T))

# Subset prevalence data
prevalence_data = agg_data_all[,c("Seasonality","access","EIR","prev_allInf_210_baseline_yearly")] %>% distinct()

# Check data
min(prevalence_data$prev_allInf_210_baseline_yearly)
max(prevalence_data$prev_allInf_210_baseline_yearly)

# Round and create category
prevalence_data$prev_210_round = round(prevalence_data$prev_allInf_210_baseline_yearly, digits = 2)
prevalence_data$prev_210_cat = cut(prevalence_data$prev_210_round, breaks = seq(0,0.6,0.1))

# Check data
table(prevalence_data$prev_210_cat, prevalence_data$Seasonality)

# Clean up for merge
prevalence_data_merge = prevalence_data[,c("Seasonality","access","EIR","prev_210_cat")] %>% distinct()


########################
### LOAD PARAM TABLE ###
########################

# Load table
param_tab = read.table(paste0(GROUP,exp,"/param_tab.txt"), sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)

# Columns order
param_cols = colnames(param_tab)

# Merge prevalence category
param_tab_prev = merge(param_tab, prevalence_data_merge,
                       by = c("Seasonality","access","EIR"), all.x=T)

# Save new table
write.table(param_tab_prev, paste0(GROUP,exp,"/param_tab_prev.txt"), sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)

