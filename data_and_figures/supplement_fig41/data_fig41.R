# INTRO ----
#
# Visualises sensitivity results
# Written by Lydia Braunack-Mayer


# SET UP ----

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- c("Obj6_Scen3_PreEryth_BloodStage")

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("Reduction_CumCPPY_age5", "Reduction_SevCumCPPY_age5", "Reduction_CumCPPY_age10", "Reduction_SevCumCPPY_age10")

# Load packages
library(dplyr)
library(tidyr)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"

if (!dir.exists(paste0("./Outputs"))) dir.create(paste0("./Outputs"))

load(paste0(GROUP_dr, exp, "/param_ranges_manual.RData"))
param_ranges_cont



# LOAD DATA ----

# Import settings
setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/sensitivity/*"))
setting <- setting[setting != paste0(GROUP_dr, exp, "/gp/trained/sensitivity/err")]
setting <- setting[setting != paste0(GROUP_dr, exp, "/gp/trained/sensitivity/param_ranges_manual.RData")]
setting_id <- sub("_cv_sidx.RData", "", sub(".*agg_", "", sub(".*seeds_", "", setting)))


# Import total effect sizes for each setting
df <- data.frame("S_eff" = c(), "T_eff" = c(), scenario = c())

for (i in 1:length(setting)) {
  
  load(setting[i]) #loads list called sobol_idx_list
  
  sobol_idx_list <- as.data.frame(sobol_idx_list[-3])
  
  sobol_idx_list$S_eff <- sobol_idx_list$S_eff / sum(sobol_idx_list$S_eff) # rescale so total = 1
  sobol_idx_list$T_eff <- sobol_idx_list$T_eff / sum(sobol_idx_list$T_eff) # rescale so total = 1
  
  sobol_idx_list$scenario <- setting_id[i]
  sobol_idx_list$parameter <- rownames(param_ranges_cont)
  
  df <- rbind(df, sobol_idx_list)
  
}

df <- df %>%
  separate(col = scenario, 
           into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Seed", "Temp", "Outcome", "Age"),
           sep = "_",
           remove = FALSE)



# IMPORT BASELINE PREVALENCE ----

# Import from csv
prev <- read.csv(paste0("/scicore/home/penny/brauna0000/M3TPP/Experiments/", exp, "/Outputs/Prevalence_prior_to_intervention.csv"))
prev <- unique(prev[, c("Seasonality", "EIR", "Access", "AnnualPrev210.2030")])

# Merge into data
df <- merge(df, prev, by = c("Seasonality", "EIR", "Access"))



# FORMAT DATA ----

# Subset for predictors of interest
df$OutcomeLabel <- paste0("Reduction_", df$Outcome, "_", df$Age)
df <- df %>%
  filter(OutcomeLabel %in% pred_list)

df <- df %>%
  mutate(Seasonality = case_match(Seasonality,
                                  "seas4mo" ~ "4 month season",
                                  "seas6mo" ~ "6 month season"),
         Access = case_match(Access,
                             "0.04" ~ "10% access to care",
                             "0.2412" ~ "50% access to care"),
         parameter = case_match(parameter,
                                "EfficacyPEV" ~ "Pre-liver stage initial efficacy [30 - 100%]",
                                "HalflifePEV" ~ "Pre-liver stage protection half-life [30 - 500 days]",
                                "KdecayPEV" ~ "Pre-liver stage decay shape [0 - 10]",
                                "EfficacyBSV" ~ "Blood stage initial efficacy [30 - 100%]",
                                "HalflifeBSV" ~ "Blood stage protection half-life [30 - 500 days]",
                                "KdecayBSV" ~ "Blood stage decay shape [0 - 10]"),
         OutcomeLabel = case_match(OutcomeLabel,
                              "Reduction_CumCPPY_age5" ~ "Reduction in age 5\ncumulative uncomp. cases",
                              "Reduction_SevCumCPPY_age5" ~ "Reduction in age 5\ncumulative severe cases",
                              "Reduction_CumCPPY_age10" ~ "Reduction in age 10\ncumulative uncomp. cases",
                              "Reduction_SevCumCPPY_age10" ~ "Reduction in age 10\ncumulative severe cases"),
         label = case_when(T_eff >= 0.08 ~ paste0(round(T_eff*100, 0), "%"),
                           T_eff < 0.08 ~ ""),
         AnnualPrev = paste0(round(AnnualPrev210.2030 * 100), "%")) %>%
  mutate(parameter = factor(parameter, levels = c("Pre-liver stage initial efficacy [30 - 100%]",
                                                  "Pre-liver stage protection half-life [30 - 500 days]",
                                                  "Pre-liver stage decay shape [0 - 10]",
                                                  "Blood stage initial efficacy [30 - 100%]",
                                                  "Blood stage protection half-life [30 - 500 days]",
                                                  "Blood stage decay shape [0 - 10]")),
         OutcomeLabel = factor(OutcomeLabel, levels = c("Reduction in age 5\ncumulative uncomp. cases",
                                                        "Reduction in age 5\ncumulative severe cases",
                                                        "Reduction in age 10\ncumulative uncomp. cases",
                                                        "Reduction in age 10\ncumulative severe cases")),
         AnnualPrev = factor(AnnualPrev, levels = unique(AnnualPrev)[order(as.numeric(sub("%", "", unique(AnnualPrev))))]))




# WRITE TO FILE ----

saveRDS(df, "./data_and_figures/supplement_fig41/data_fig41.rds")
