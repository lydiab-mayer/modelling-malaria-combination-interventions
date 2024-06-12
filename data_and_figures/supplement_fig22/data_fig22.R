############################################################
#
# Visualises cumulative cases per person year
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())
library(dplyr)
library(tidyr)

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "Obj6_Scen2_PreEryth"

# # !!! Insert your outcome ID here !!!
id <- c("CumCPPY_", "SevCumCPPY_")

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

setting <- sub(".txt.*", "", sub(".*agg_", "", Sys.glob(paste0(GROUP_dr, exp, "/postprocessing/agg*"))))
setting

# !!! Define which setting you want to plot !!! Should be specified as integer (e.g. '1')
setting_id <- 17

# Load unaggregated postprocessing data
df <- read.table(paste0(GROUP_dr, exp, "/postprocessing/seeds_", setting[setting_id], ".txt"), sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)

# !!! Define the parameter set to plot !!! Should be specified as integer/integer vector (e.g. '1')
params <- rownames(param_ranges_cont)
param_index <- c(3, 32)
param_id <- unique(df[, params])[param_index, ]
param_id

dfOut <- vector(mode = "list")
dfOutAggregate <- vector(mode = "list")

for (index in 1:length(param_index)) {
  # Extract desired outcomes
  
  dfOut[[index]] <- df %>%
    select(c(rownames(param_ranges_cont)), "seed", starts_with(c(paste0("Baseline_", id), paste0("SMC_", id), paste0("Layer_", id)))) %>%
    filter(.data[[params[[1]]]] == param_id[[1]][index],
           .data[[params[[2]]]] == param_id[[2]][index],
           .data[[params[[3]]]] == param_id[[3]][index])


  # Change to long format
  dfOut[[index]] <- dfOut[[index]] %>%
    pivot_longer(cols = starts_with(c(paste0("Baseline_", id), paste0("SMC_", id), paste0("Layer_", id))),
                 names_to = c("Intervention", "Outcome", "AgeGroup"),
                 names_sep = "_",
                 values_to = "Value")
  
  # Remove prefixes
  dfOut[[index]]$AgeGroup <- as.numeric(sub("age", "", dfOut[[index]]$AgeGroup))
  
  # Recode variables
  dfOut[[index]] <- dfOut[[index]] %>%
    mutate(seed = as.factor(seed),
           Outcome = case_match(Outcome,
                                "CumCPPY" ~ "Uncomplicated malaria",
                                "SevCumCPPY" ~ "Severe disease"),
           Intervention = case_match(Intervention,
                                     "Baseline" ~ "No intervention",
                                     "SMC" ~ "SMC",
                                     "Layer" ~ "SMC + pre-erythrocytic intervention"),
           Group = paste0(seed, Intervention))
  dfOut[[index]]$Outcome <- factor(dfOut[[index]]$Outcome, levels = c("Uncomplicated malaria", "Severe disease"))
  
  # Calculate median, min and max across seeds
  dfOutAggregate[[index]] <- dfOut[[index]] %>%
    group_by(Outcome,
             Intervention,
             AgeGroup) %>%
    summarise(medianValue = median(Value),
              minValue = min(Value),
              maxValue = max(Value))
}


# ----------------------------------------------------------
# Define plot label
# ----------------------------------------------------------

# Format setting labels
setting_label <- separate(as.data.frame(setting), 
                          col = setting,
                          into = c("Intervention", "Seasonality", "Setting", "EIR", "AgeGroup", "Decay", "Access", "Coverage"), 
                          sep = "_")
setting_label <- setting_label %>%
  mutate(Seasonality = case_match(Seasonality,
                                  "seas4mo" ~ "4 month",
                                  "seas6mo" ~ "6 month"),
         Access = case_match(as.character(Access),
                             "0.04" ~ "10%",
                             "0.2412" ~ "50%"))
setting_label <- setting_label[setting_id, ]


# Format intervention parameter labels and combine with setting labels
param_label <- vector(mode = "list")

for (index in 1:length(param_index)) {
  param_label[[index]] <- paste0("Imperfect seasonal coverage scenario | ",
                                 round(df[param_index[index], ]$AnnualPrev210.2030*100), "% PfPR2-10 | ",
                                 setting_label$Seasonality, " seasonal profile | ",
                                 setting_label$Access, " access to first-line care\n",
                                 round(param_id[index, 1]), " day pre-erythrocytic half-life | ",
                                 round(param_id[index, 2] * 100), "% pre-erythrocytic efficacy | ",
                                 round(param_id[index, 3], 1), " pre-erythrocytic decay shape")
}

# ----------------------------------------------------------
# Write to file
# ----------------------------------------------------------

saveRDS(dfOut, "./data_and_figures/supplement_fig22/data_fig22.rds")
saveRDS(dfOutAggregate, "./data_and_figures/supplement_fig22/data_aggregate_fig22.rds")
saveRDS(param_label, "./data_and_figures/supplement_fig22/label_fig22.rds")

