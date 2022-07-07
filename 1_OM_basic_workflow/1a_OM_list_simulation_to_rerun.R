##############################################
##############################################
###                                        ###
### STEP 1a: LIST OF FAILED OM SIMULATIONS ###
###                                        ###
##############################################
##############################################

### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Creates new parameter table for failed simulations
###
### adapted 10.11.2021
### narimane.nekkab@swisstph.ch
###
### -------------------------------------------------------------------------

##############################
### STEP 1; RUN ON CLUSTER ###
##############################

# !!!!! Change experiment name !!!!!

# Run this code (when inside om folder): 
# cd /scicore/home/penny/GROUP/M3TPP/.../om/
# ls >> /scicore/home/penny/GROUP/M3TPP/.../filenames_of_simulations_done.txt

##############
### HEADER ###
##############

rm(list = ls())
library(sjmisc)
library(plyr)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Experiment
exp ="..."

############
### DATA ###
############

# Load list of successful experiments and get names files
f <- read.table(paste0("/scicore/home/penny/GROUP/M3TPP/",exp, "/filenames_of_simulations_done.txt"), header = F, stringsAsFactors = F) 
f_1 <- f[grepl("_out.txt", f$V1), "V1"]
f_2 = gsub("_out.txt","", f_1)

# Load parameter table files (add more a files if more param_tab_X)
a <- read.table(paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab.txt"),header = T, stringsAsFactors = F)

# Name of simulated OM files
simulated <- which(paste0(a$Scenario_Name,"_",a$SEED) %in% f_2)

# Select failed simulations
param_tab_new <- a[-simulated,]

# Table
table(param_tab_new$Seasonality, param_tab_new$EIR)
table(a$Seasonality, a$EIR)

# Percent failures
sim_table = data.frame(Seasonality=a$Seasonality, EIR=a$EIR, value=1) %>% dplyr::group_by(Seasonality, EIR) %>% dplyr::summarise(total_sims = n(value))
fail_table = data.frame(Seasonality=param_tab_new$Seasonality, EIR=param_tab_new$EIR, value=1) %>% dplyr::group_by(Seasonality, EIR) %>% dplyr::summarise(total_fails = n(value))
# Merge
sim_fail_table = merge(sim_table, fail_table, by=c("Seasonality","EIR"), all = T)
# Percent
sim_fail_table$percent = sim_fail_table$total_fails / sim_fail_table$total_sims

# Total
sum(sim_fail_table$total_fails)

# Save new parameter table for new simulations
write.table(param_tab_new,paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/param_tab_resubmit.txt"), sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)


# !!! Follow step 1b !!!


