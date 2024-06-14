# INTRO ----
#
# Visualises sensitivity results
# Written by Lydia Braunack-Mayer


# SET UP ----

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
experiments <- c("Obj6_Scen2_PreEryth", "Obj6_Scen2_BloodStage", "Obj6_Scen2_PreEryth_BloodStage")

# Load packages
library(dplyr)
library(tidyr)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"



# LOAD DATA ----

data <- list()

for (exp in experiments) {
  
  # Load parameter ranges
  load(paste0(GROUP_dr, exp, "/param_ranges_manual.RData"))

  # Import settings
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/sensitivity/*"))
  setting <- setting[setting != paste0(GROUP_dr, exp, "/gp/trained/sensitivity/err")]
  setting <- setting[setting != paste0(GROUP_dr, exp, "/gp/trained/sensitivity/param_ranges_manual.RData")]
  setting_id <- sub("_cv_sidx.RData", "", sub(".*agg_", "", sub(".*seeds_", "", setting)))
  
  # Set up  dataframe to store results
  df <- data.frame("scenario" = setting_id,
                   "maximum" = NA,
                   "minimum" = NA,
                   "quantile0.75" = NA,
                   "median" = NA,
                   "quantile0.25" = NA)
  
  for (set in setting_id) {
    
    # Load Rdata object containing results of the sensitivity analysis, called 'sobol_idx_list'
    load(paste0(GROUP_dr, exp, "/gp/trained/sensitivity/seeds_", set, "_cv_sidx.RData"))
    
    # Store sampled parameter values and corresponding predicted values from sensitivity analysis
    sensData <- sobol_idx_list$SA$X
    sensData$outcome <- sobol_idx_list$SA$y
    names(sensData) <- c(rownames(param_ranges_cont), "outcome")
    
    # Return interquartile range on outcome
    df[df$scenario == set, c("minimum", "quantile0.25", "median", "quantile0.75", "maximum")] <- unname(quantile(sensData[, "outcome"], c(0, 0.25, 0.5, 0.75, 1)))
 
  }

  # Separate out model scenario labels
  df <- df %>%
    separate(col = scenario, 
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Seed", "Temp", "Outcome", "Age"),
             sep = "_",
             remove = FALSE)
  
  # Import baseline prevalence from csv
  prev <- read.csv(paste0("/scicore/home/penny/brauna0000/M3TPP/Experiments/", exp, "/Outputs/Prevalence_prior_to_intervention.csv"))
  prev <- unique(prev[, c("Seasonality", "EIR", "Access", "AnnualPrev210.2030")])
  
  # Merge into data
  df <- merge(df, prev, by = c("Seasonality", "EIR", "Access"))
  
  data[[exp]] <- df
}

# Combine all datasets
data <- do.call("rbind", data)



# FORMAT DATA ----

# Re-label factor variables
data <- data %>%
  mutate(Seasonality = case_match(Seasonality,
                                  "seas4mo" ~ "4 month transmission season",
                                  "seas6mo" ~ "6 month transmission season"),
         Access = case_match(Access,
                             "0.04" ~ "10% access to care",
                             "0.2412" ~ "50% access to care"),
         Experiment = case_match(Experiment,
                                 "Obj6Scen2PreEryth" ~ "Pre-erythrocytic activity alone",
                                 "Obj6Scen2BloodStage" ~ "Blood stage activity alone",
                                 "Obj6Scen2PreErythBloodStage" ~ "Pre-erythrocytic and blood stage activity"),
         Outcome = case_match(Outcome,
                              "CumCPPY" ~ "Reduction in cumulative uncomplicated cases",
                              "SevCumCPPY" ~ "Reduction in cumulative severe cases"),
         Age = case_match(Age,
                          "age5" ~ "5 years old",
                          "age10" ~ "10 years old"),
         AnnualPrev = paste0(round(AnnualPrev210.2030 * 100), "%"))

# Manually override small differences in annual prevalence for same model scenarios
# Note that these differences occur due to stochastic variation in model outputs
data$AnnualPrev[data$AnnualPrev == "14%"] <- "15%"
data$AnnualPrev[data$AnnualPrev == "32%"] <- "33%"
data$AnnualPrev[data$AnnualPrev == "44%"] <- "45%"
data$AnnualPrev[data$AnnualPrev %in% c("63%", "65%")] <- "64%"

# Order factor levels
data <- data %>%
  mutate(Experiment = factor(Experiment, levels = c("Pre-erythrocytic activity alone",
                                                    "Blood stage activity alone",
                                                    "Pre-erythrocytic and blood stage activity")),
         Outcome = factor(Outcome, levels = c("Reduction in cumulative uncomplicated cases", 
                                              "Reduction in cumulative severe cases")),
         Age = factor(Age, levels = c("5 years old", 
                                      "10 years old")),
         AnnualPrev = factor(AnnualPrev, levels = unique(AnnualPrev)[order(as.numeric(sub("%", "", unique(AnnualPrev))))]))




# WRITE TO FILE ----

saveRDS(data, "./data_and_figures/supplement_fig43/data_fig43.rds")
