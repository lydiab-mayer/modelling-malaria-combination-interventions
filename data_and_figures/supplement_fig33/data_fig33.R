# INTRO ----
#
# Visualises parameter relationships
#
# Written by Lydia Braunack-Mayer

# SETUP ----

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "Obj6_Scen4_BloodStage"

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
predictors <- c("Reduction_CumCPPY_age5", "Reduction_SevCumCPPY_age5", "Reduction_CumCPPY_age10", "Reduction_SevCumCPPY_age10")

# Load required packages
library(tidyr)
library(dplyr)

GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont

# Define and load scenarios
setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/", predictors[1], "/*"))
(setting_id <- unique(sub("seeds_", "", sub(paste0("_", predictors[1], ".*"), "", basename(setting)))))
index <- 17
setting_id[index]


# GENERATE DATA ----

# Define number of segments for parameter ranges
N <- 51

# Define empty matrix to store outputs
data <- matrix(NA, nrow = length(predictors) * nrow(param_ranges_cont) * (N - 1), ncol = 8)
colnames(data) <- c("Scenario", "Outcome", "Parameter", "segLower", "segUpper", "quantile0.25", "median", "quantile0.75") 
nRow <- nrow(data)

# Populate matrix with scenario, outcome and parameter names
data[, "Scenario"] <- setting_id[index]
data[, "Outcome"] <- rep(predictors, each = nRow / length(predictors))
data[, "Parameter"] <- rep(rep(rownames(param_ranges_cont), each = N - 1), length(predictors))

# For each parameter, cut the parameter space into segments and store segment boundaries
for (pred in predictors) {
  for (param in rownames(param_ranges_cont)) {
    split <- seq(0, 1, length = N) # Note this assumes that parameters have been scaled to c(0, 1)
    data[data[, "Outcome"] == pred & data[, "Parameter"] == param, "segLower"] <- split[1:(N - 1)]
    data[data[, "Outcome"] == pred & data[, "Parameter"] == param, "segUpper"] <- split[2:N]
  }
}


for (pred in predictors) {
  # Load Rdata object containing results of the sensitivity analysis, called 'sobol_idx_list'
  load(paste0(GROUP_dr, exp, "/gp/trained/sensitivity/seeds_", setting_id[index], "_", pred, "_cv_sidx.RData"))
  
  # Store sampled parameter values and corresponding predicted values from sensitivity analysis
  sensData <- sobol_idx_list$SA$X
  sensData$outcome <- sobol_idx_list$SA$y
  names(sensData) <- c(rownames(param_ranges_cont), pred)
  
  # For each parameter, extract median outcome and quantile range for each segment
  for (i in 1:nRow) {
    if(data[i, "Outcome"] == pred) {
      temp <- sensData[sensData[, data[i, "Parameter"]] >= data[i, "segLower"] & sensData[, data[i, "Parameter"]] < data[i, "segUpper"], ]
      data[i, c("quantile0.25", "median", "quantile0.75")] <- unname(quantile(temp[, pred], c(0.25, 0.5, 0.75)))
    }
  }
}



# FORMAT DATA ----

# Format outputs as dataframe
data <- as.data.frame(data) %>%
  mutate(Outcome = as.factor(Outcome),
         Parameter = as.factor(Parameter),
         segLower = as.numeric(segLower),
         segUpper = as.numeric(segUpper),
         quantile0.25 = as.numeric(quantile0.25),
         median = as.numeric(median),
         quantile0.75 = as.numeric(quantile0.75))

# Reverse scaling of parameter values
for (param in rownames(param_ranges_cont)) {
  data[data$Parameter == param, "segLower"] <- data[data$Parameter == param, "segLower"] * (param_ranges_cont[param, 2] - param_ranges_cont[param, 1]) + param_ranges_cont[param, 1]
  data[data$Parameter == param, "segUpper"] <- data[data$Parameter == param, "segUpper"] * (param_ranges_cont[param, 2] - param_ranges_cont[param, 1]) + param_ranges_cont[param, 1]
}

# Remove missing values, since these represent results outside of the scope of the sensitivity analysis
data <- data[!is.na(data$median), ]

# Format factor variables
data <- data %>%
  mutate(Outcome_age = case_match(Outcome,
                                  c("Reduction_CumCPPY_age5", "Reduction_SevCumCPPY_age5") ~ "5 years",
                                  c("Reduction_CumCPPY_age10", "Reduction_SevCumCPPY_age10") ~ "10 years"),
         Outcome = case_match(Outcome,
                              c("Reduction_CumCPPY_age5", "Reduction_CumCPPY_age10") ~ "Uncomplicated cases",
                              c("Reduction_SevCumCPPY_age5", "Reduction_SevCumCPPY_age10") ~ "Severe cases")) %>%
  mutate(Outcome = factor(Outcome,
                          levels = c("Uncomplicated cases",
                                     "Severe cases")))


## CONSTRUCT ANNOTATION ----
setting_label <- separate(as.data.frame(setting_id), 
                          col = setting_id,
                          into = c("Intervention", "Seasonality", "Setting", "EIR", "Agegroup", "Decay", "Access", "Seed"), 
                          sep = "_")
setting_label <- setting_label %>%
  mutate(Seasonality = case_match(Seasonality,
                                  "seas4mo" ~ "4 month",
                                  "seas6mo" ~ "6 month"),
         Access = case_match(as.character(Access),
                             "0.04" ~ "10%",
                             "0.2412" ~ "50%"))
setting_label <- setting_label[index, ]

tag <- paste0("Random allocation scenario | 33% PfPR2-10 | ",
              setting_label$Seasonality, " seasonal profile | ",
              setting_label$Access, " access to first-line care")



## WRITE TO FILE ----

saveRDS(data, "./data_and_figures/supplement_fig33/data_fig33.rds")
saveRDS(tag, "./data_and_figures/supplement_fig33/label_fig33.rds")
