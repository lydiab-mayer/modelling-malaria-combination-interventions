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

# !!! Insert your experiment name(s) here as a string, e.g. "MyExperiment" !!!
experiments <- c("Obj6_Scen2_PreErythProfiles", "Obj6_Scen2_LayerCounterfactual_PreErythProfiles")

# # !!! Insert your outcome ID here !!!
id <- c("CumCPPY_", "SevCumCPPY_")

GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"

load(paste0(GROUP_dr, experiments[1], "/param_ranges.RData"))
param_ranges_cont


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

# Set up empty dataframes
dfOut <- vector(mode = "list")
dfOutAggregate <- vector(mode = "list")

for (exp in experiments) {

  # Load unaggregated postprocessing data
  df <- read.table(Sys.glob(paste0(GROUP_dr, exp, "/postprocessing/seeds*")), sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
  
  # Extract desired outcomes
  dfOut[[exp]] <- df %>%
    select(c(rownames(param_ranges_cont)), "seed", starts_with(c(paste0("Baseline_", id), paste0("Counterfactual_", id), paste0("Intervention_", id))))

  # Change to long format
  dfOut[[exp]] <- dfOut[[exp]] %>%
    pivot_longer(cols = starts_with(c(paste0("Baseline_", id), paste0("Counterfactual_", id), paste0("Intervention_", id))),
                 names_to = c("Intervention", "Outcome", "AgeGroup"),
                 names_sep = "_",
                 values_to = "Value")
  
  # Remove prefixes
  dfOut[[exp]]$AgeGroup <- as.numeric(sub("age", "", dfOut[[exp]]$AgeGroup))
  
  # Recode variables
  if (grepl("LayerCounterfactual", exp, fixed = TRUE)) {
    counterfactualLabel <- "Pre-liver stage therapeutic"
  } else {
    counterfactualLabel <- "SMC"
  }
  
  dfOut[[exp]] <- dfOut[[exp]] %>%
    mutate(seed = as.factor(seed),
           Outcome = case_match(Outcome,
                                "CumCPPY" ~ "Uncomplicated malaria",
                                "SevCumCPPY" ~ "Severe disease"),
           Intervention = case_match(Intervention,
                                     "Baseline" ~ "No intervention",
                                     "Counterfactual" ~ counterfactualLabel,
                                     "Intervention" ~ "SMC + pre-liver stage therapeutic"),
           TherapeuticProfile = case_match(Halflife,
                                           param_ranges_cont[1, 1] ~ "Short duration, moderate efficacy therapeutic",
                                           param_ranges_cont[1, 2] ~ "Long duration, high efficacy therapeutic"),
           Group = paste0(seed, Intervention))
  dfOut[[exp]]$Outcome <- factor(dfOut[[exp]]$Outcome, levels = c("Uncomplicated malaria", "Severe disease"))
  
  # Add experiment ID
  dfOut[[exp]]$Experiment <- exp
  
  # Calculate median, min and max across seeds
  dfOutAggregate[[exp]] <- dfOut[[exp]] %>%
    group_by(Experiment,
             Outcome,
             Intervention,
             AgeGroup,
             TherapeuticProfile) %>%
    summarise(medianValue = median(Value),
              minValue = min(Value),
              maxValue = max(Value))
}

# Combine into single dataframe
dfOut <- bind_rows(dfOut)
dfOutAggregate <- bind_rows(dfOutAggregate)


# ----------------------------------------------------------
# Define plot label
# ----------------------------------------------------------

# Get setting label
setting <- sub(".txt.*", "", sub(".*seeds_", "", Sys.glob(paste0(GROUP_dr, exp, "/postprocessing/seeds*"))))

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
setting_label


# Format intervention parameter labels and combine with setting labels
param_label <- vector(mode = "list")

for (index in 1:2) {
  param_label[[index]] <- paste0("Imperfect seasonal coverage scenario | ",
                                 round(mean(df$AnnualPrev210.2030*100)), "% PfPR2-10 | ",
                                 setting_label$Seasonality, " seasonal profile | ",
                                 setting_label$Access, " access to first-line care\n",
                                 round(param_ranges_cont[1, index]), " day pre-liver stage half-life | ",
                                 round(param_ranges_cont[2, index] * 100), "% pre-liver stage efficacy | ",
                                 round(param_ranges_cont[3, index], 1), " pre-liver stage decay shape")
}

# ----------------------------------------------------------
# Write to file
# ----------------------------------------------------------

saveRDS(dfOut, "./data_and_figures/supplement_fig22/data_fig22.rds")
saveRDS(dfOutAggregate, "./data_and_figures/supplement_fig22/data_aggregate_fig22.rds")
saveRDS(param_label, "./data_and_figures/supplement_fig22/label_fig22.rds")

