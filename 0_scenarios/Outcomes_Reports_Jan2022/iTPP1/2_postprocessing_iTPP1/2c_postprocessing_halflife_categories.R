#################################
#################################
###                           ###
### STEP 2c: POST-PROCESSING  ###
###                           ###
#################################
#################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Create halflife categories
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
exp = "iTPP1a_3mo59mo_1round_deploymentTS18"

# Create new directory
destination_folder = paste0(GROUP,exp,"/postprocessing")

# Create destination folder if it doesn't exist
if(!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

#################
### LOAD DATA ###
#################

# Get postprocessing aggregate files 
agg_names = list.files(path = paste0(GROUP,exp,"/postprocessing_kcat/"), pattern = "agg", full.names = FALSE)

# Load
agg_data_all <- do.call(rbind, lapply(paste0(GROUP,exp,"/postprocessing_kcat/",agg_names), read.table, header = T))

# Get postprocessing seed files 
seed_names = list.files(path = paste0(GROUP,exp,"/postprocessing_kcat/"), pattern = "seed", full.names = FALSE)

# Load
seed_data_all <- do.call(rbind, lapply(paste0(GROUP,exp,"/postprocessing_kcat/",seed_names), read.table, header = T))


####################

# Add kdecay categories
agg_data_all$HL_cat = ifelse(agg_data_all$Halflife <= 90, "HL0.90",NA)
agg_data_all$HL_cat = ifelse(agg_data_all$Halflife > 90 & agg_data_all$Halflife <= 120, "HL90.120", agg_data_all$HL_cat)
agg_data_all$HL_cat = ifelse(agg_data_all$Halflife > 120 & agg_data_all$Halflife <= 150, "HL120.150", agg_data_all$HL_cat)
agg_data_all$HL_cat = ifelse(agg_data_all$Halflife > 150 & agg_data_all$Halflife <= 180, "HL150.180", agg_data_all$HL_cat)
agg_data_all$HL_cat = ifelse(agg_data_all$Halflife > 180, "HL180.300", agg_data_all$HL_cat)

# Add kdecay categories
seed_data_all$HL_cat = ifelse(seed_data_all$Halflife <= 90, "HL0.90",NA)
seed_data_all$HL_cat = ifelse(seed_data_all$Halflife > 90 & seed_data_all$Halflife <= 120, "HL90.120", seed_data_all$HL_cat)
seed_data_all$HL_cat = ifelse(seed_data_all$Halflife > 120 & seed_data_all$Halflife <= 150, "HL120.150", seed_data_all$HL_cat)
seed_data_all$HL_cat = ifelse(seed_data_all$Halflife > 150 & seed_data_all$Halflife <= 180, "HL150.180", agg_data_all$HL_cat)
seed_data_all$HL_cat = ifelse(seed_data_all$Halflife > 180, "HL180.300", seed_data_all$HL_cat)


# Categories
split_categories = c("Seasonality","access_cat","prev_210_cat", "kdecay_cat", "HL_cat")

#################### aggregate

# Split
list_tab_splits_agg = split(agg_data_all,  list(agg_data_all[,split_categories[1]],
                                                agg_data_all[,split_categories[2]],
                                                agg_data_all[,split_categories[3]],
                                                agg_data_all[,split_categories[4]],
                                                agg_data_all[,split_categories[5]]), sep="_")


# Write each split to a file named according to the corresponding seasonality and access
for (i in 1:length(list_tab_splits_agg)) {
  
  # File name
  exp_name = str_replace_all(exp, "_", "")
  
  # Full name
  final_table_dest = paste(destination_folder,
                           "/agg_", 
                           exp_name, "_", 
                           names(list_tab_splits_agg)[i], 
                           ".txt", sep="")
  # Save
  write.table(list_tab_splits_agg[[i]], final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
}

#################### seeds

# Split
list_tab_splits_seed = split(seed_data_all,  list(seed_data_all[,split_categories[1]],
                                                  seed_data_all[,split_categories[2]],
                                                  seed_data_all[,split_categories[3]],
                                                  seed_data_all[,split_categories[4]],
                                                  seed_data_all[,split_categories[5]]), sep="_")


# Write each split to a file named according to the corresponding seasonality and access
for (i in 1:length(list_tab_splits_seed)) {
  
  # File name
  exp_name = str_replace_all(exp, "_", "")
  
  # Full name
  final_table_dest = paste(destination_folder,
                           "/seeds_", 
                           exp_name, "_", 
                           names(list_tab_splits_seed)[i], 
                           ".txt", sep="")
  # Save
  write.table(list_tab_splits_seed[[i]], final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
}
