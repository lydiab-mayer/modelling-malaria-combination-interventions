###########################
# Script to split a parameter table by setting.
# A setting in this case is defined by EIR, seasonality and biting pattern.
# Each split is written to a separate .txt file.
#
# INPUTS: 
#       param_table: full path of the .txt file containing the parameters for each scenario
#       destination_folder: folder where all the split files will be saved
#
# OUTPUTS:
#   For each setting, a .txt file containing the parameters for the scenarios corresponding to that setting
# 
# created 14.05.2019
# monica.golumbeanu@unibas.ch
#
# adapated 24.11.2021
# narimane.nekkab@unibas.ch
########################## 

args = commandArgs(TRUE)
param_table = args[1]
destination_folder = args[2]
param_cat_file = args[3]
# Read the parameter table
param_tab = read.table(param_table, sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)

# Label to be added in front of the names of splits
experiment_name = strsplit(param_tab$Scenario_Name[1],"_")[[1]]
experiment_name = paste0(experiment_name[-length(experiment_name)],collapse="")

# Create destination folder if it doesn't exist
if(!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

##############################################################################
# Split the parameter table before interventions by seasonality and access

# Create new access category
param_tab$access_cat = ifelse(param_tab$access == 0.0400, "lowAccess", param_tab$access)
param_tab$access_cat = ifelse(param_tab$access == 0.0972, "modAccess", param_tab$access_cat)
param_tab$access_cat = ifelse(param_tab$access == 0.2412, "highAccess", param_tab$access_cat)
# param_tab$access_cat = factor(param_tab$access_cat, levels = c("lowAccess","modAccess","highAccess"))

# Categories
split_categories = c("Seasonality","access_cat","prev_210_cat","SEED")

# Split
list_tab_splits = split(param_tab,  list(param_tab[,split_categories[1]],
                                         param_tab[,split_categories[2]],
                                         param_tab[,split_categories[3]],
                                         param_tab[,split_categories[4]]), sep="_" )


# Write each split to a file named according to the corresponding seasonality and access
for (i in 1:length(list_tab_splits)) {
  final_table_dest = paste(destination_folder, experiment_name, "_",
                           names(list_tab_splits)[i], ".txt", sep="")
  write.table(list_tab_splits[[i]], final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
}
