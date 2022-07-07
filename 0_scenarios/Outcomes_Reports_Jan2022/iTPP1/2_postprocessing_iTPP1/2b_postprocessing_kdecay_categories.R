#################################
#################################
###                           ###
### STEP 2b: POST-PROCESSING  ###
###                           ###
#################################
#################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Create k decay categories
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

# Create new directory
destination_folder = paste0(GROUP,exp,"/postprocessing_kcat")

# Create destination folder if it doesn't exist
if(!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

#################
### LOAD DATA ###
#################

# Get postprocessing aggregate files 
agg_names = list.files(path = paste0(GROUP,exp,"/postprocessing_outcomes/"), pattern = "agg", full.names = FALSE)

# Load
agg_data_all <- do.call(rbind, lapply(paste0(GROUP,exp,"/postprocessing_outcomes/",agg_names), read.table, header = T))

# Aggregate
agg_data_all = agg_data_all %>% 
  dplyr::group_by(Scenario_Name, kdecay, prev_210_cat) %>% 
  dplyr::summarise_at(c(names(agg_data_all)[which(names(agg_data_all)=="prev_red1_allInf_allages_yearly"):length(names(agg_data_all))]), median, na.rm=TRUE)

######################################

# Get postprocessing seed files 
seed_names = list.files(path = paste0(GROUP,exp,"/postprocessing_outcomes/"), pattern = "seed", full.names = FALSE)

# Load
seed_data_all <- do.call(rbind, lapply(paste0(GROUP,exp,"/postprocessing_outcomes/",seed_names), read.table, header = T))


####################

# Add kdecay categories
agg_data_all$kdecay_cat = ifelse(agg_data_all$kdecay < 1, "biphasic","sigmoidal")
agg_data_all$kdecay_cat = ifelse(agg_data_all$kdecay > 2.5, "verysigmoidal", agg_data_all$kdecay_cat)
seed_data_all$kdecay_cat = ifelse(seed_data_all$kdecay < 1, "biphasic","sigmoidal")
seed_data_all$kdecay_cat = ifelse(seed_data_all$kdecay > 2.5, "verysigmoidal", seed_data_all$kdecay_cat)

# Categories
split_categories = c("Seasonality","access_cat","prev_210_cat", "kdecay_cat")

#################### aggregate

# Split
list_tab_splits_agg = split(agg_data_all,  list(agg_data_all[,split_categories[1]],
                                                agg_data_all[,split_categories[2]],
                                                agg_data_all[,split_categories[3]],
                                                agg_data_all[,split_categories[4]]), sep="_")


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
                                                seed_data_all[,split_categories[4]]), sep="_")


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
