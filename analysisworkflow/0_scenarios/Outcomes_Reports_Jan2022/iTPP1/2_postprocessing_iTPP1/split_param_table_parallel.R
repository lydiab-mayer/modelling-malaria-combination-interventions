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
########################## 

# modified 05.10.2021
# narimane.nekkab@swisstph.ch
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

# Split the parameter table by seasonality, biting pattern and LAI decay shape

# add kdecay category
param_tab$kcat = ifelse(param_tab$kdecay < 2.5, "lowerK", "upperK")

load(param_cat_file)
# split_categories <- names(param_cat)
# split_categories = c("seasonality","kdecay","effB","access")
split_categories = c("Seasonality","EIR","kcat")

list_tab_splits = split(param_tab,  list(param_tab[,split_categories[1]],
                                         param_tab[,split_categories[2]],
                                         param_tab[,split_categories[3]]), sep="_" )

table(list_tab_splits$seas4mo_lowerK_10$kcat)

# list_tab_splits_test = sapply(list_tab_splits, function(x) split(x, ceiling(seq_along(1:nrow(x))/2000)))
# names(list_tab_splits_test) = paste0(names(list_tab_splits))

# Write each split to a file named according to the corresponding EIR and seasonality
for (i in 1:length(list_tab_splits)) {
  final_table_dest = paste(destination_folder,"split/", experiment_name, "_",
                           names(list_tab_splits)[i], ".txt", sep="")
  write.table(list_tab_splits[[i]], final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
}
