# INTRO ----
#
# Visualises parameter relationships
#
# Written by Lydia Braunack-Mayer

# SETUP ----

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
experiments <- c("Obj6_v2_PreEryth", "Obj6_Scen2_PreEryth", "Obj6_Scen3_PreEryth", "Obj6_Scen4_PreEryth")

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
predictors <- c("Reduction_CumCPPY_age10", "Reduction_SevCumCPPY_age10")

# Load required packages
library(tidyr)
library(dplyr)

# Define directories
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"

# Define function to generate predictions
predict.grid <- function(param.ranges, grid.ranges, ngrid, model, scale = TRUE) {
  
  require(hetGP)
  
  # Set up
  D <- nrow(param.ranges)
  scale.params <- t(param.ranges)
  scale.grid <- t(grid.ranges)
  
  # Scale grid to c(0, 1)
  if (scale) {
    for (i in 1:D) {
      scale.grid[, i] <- (scale.grid[, i] - scale.params[1, i]) / (scale.params[2, i] - scale.params[1, i])
    }
  }
  
  # Create grid of scenarios
  scenarios <- list()
  
  for (i in 1:D) {
    scenarios[[i]] <- seq(scale.grid[1, i], scale.grid[2, i], length.out = ngrid[i])
  }
  
  scenarios <- expand.grid(scenarios)
  names(scenarios) <- rownames(param.ranges)
  
  
  # Make predictions using emulator
  preds <- predict(x = as.matrix(scenarios), object = model)
  scenarios$mean <- preds$mean
  scenarios$sd2 <- preds$sd2
  scenarios$nugs <- preds$nugs
  
  # Covert parameter values back to original scale
  for (i in rownames(param.ranges)) {
    scenarios[, i] <- scenarios[, i] * (param.ranges[i, 2] - param.ranges[i, 1]) + param.ranges[i, 1]
  }
  
  # Calculate standard error and 95% confidence interval
  scenarios$se <- sqrt(scenarios$sd2 + scenarios$nugs)
  scenarios$cl <- qnorm(0.05, scenarios$mean, scenarios$se)
  scenarios$cu <- qnorm(0.95, scenarios$mean, scenarios$se)
  
  # Calculate target reduction
  scenarios$target <- floor(scenarios$mean / 5) * 5
  
  return(scenarios)
}

# Define wrapper function to generate predictions over multiple outcome variables
wrap.grid <- function(predictors, dir, exp, index, param.ranges, grid.ranges, ngrid) {
  
  # Generate predictions
  df <- data.frame()
  
  for (pred in predictors) {
    
    # Load GP model
    load(paste0(dir, exp, "/gp/trained/", pred, "/seeds_", index, "_", pred, "_cv.RData"))
    
    # Generate model predictions
    temp <- predict.grid(param.ranges = param.ranges, 
                         grid.ranges = grid.ranges, 
                         ngrid = ngrid, 
                         model = cv_result$GP_model)
    
    temp$scenario <- index
    
    # Format predictions
    temp <- temp %>%
      separate(col = scenario,
               into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Seed"),
               sep = "_",
               remove = FALSE)
    temp$pred <- pred
    
    # Store predictions
    df <- rbind(df, temp)
    
  }
  
  # Return outputs
  return(df)
}


# GENERATE DATA ----

# Define empty list
data <- list()

# Define indices for scenarios to plot
indices <- c(7, 8, 17, 18)

for (exp in experiments[1:3]) {
  
  # Load parameter ranges
  load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
  
  # Define and load scenarios
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/", predictors[1], "/*"))
  (setting_id <- unique(sub("seeds_", "", sub(paste0("_", predictors[1], ".*"), "", basename(setting)))))
  
  for (index in indices) {
    
    setting_id[index]
    
    # Generate grid
    df <- wrap.grid(predictors = predictors, 
                    dir = GROUP_dr, 
                    exp = exp, 
                    index = setting_id[index], 
                    param.ranges = param_ranges_cont, 
                    grid.ranges = rbind(Halflife = c(100, 500),
                                        Efficacy = c(0.3, 1.0),
                                        Kdecay = 0.7), # sustained decay
                    ngrid = c(401, 71, 1)) 
    
    # Add column for coverage
    df$Coverage <- NA
  
    data[[setting_id[index]]] <- df
  }
}

for (exp in experiments[4]) {
  
  # Load parameter ranges
  load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
  
  # Define and load scenarios
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/", predictors[1], "/*"))
  (setting_id <- unique(sub("seeds_", "", sub(paste0("_", predictors[1], ".*"), "", basename(setting)))))
  
  for (index in indices) {

    setting_id[index]
    
    # Generate grid
    df <- wrap.grid(predictors = predictors, 
                    dir = GROUP_dr, 
                    exp = exp, 
                    index = setting_id[index], 
                    param.ranges = param_ranges_cont, 
                    grid.ranges = rbind(Coverage = 0.7,
                                        Halflife = c(100, 500),
                                        Efficacy = c(0.3, 1.0),
                                        Kdecay = 0.7), # sustained decay
                    ngrid = c(1, 401, 71, 1)) 
    
    data[[setting_id[index]]] <- df
  }
}

# Combine all datasets
data <- do.call("rbind", data)

# Format data
data <- data %>%
  mutate(pred = case_match(pred,
                           "Reduction_CumCPPY_age10" ~ "Cumulative clinical cases by age 10",
                           "Reduction_SevCumCPPY_age10" ~ "Cumulative severe cases by age 10"),
         targetLabel = paste0(target, "%"),
         Experiment = case_match(Experiment,
                                 "Obj6v2PreEryth" ~ "Scenario 1\n Perfect deployment",
                                 "Obj6Scen2PreEryth" ~ "Scenario 2\nImperfect seasonal coverage",
                                 "Obj6Scen3PreEryth" ~ "Scenario 3\nImperfect deployment",
                                 "Obj6Scen4PreEryth" ~ "Scenario 4\nRandom allocation"),
         Seasonality = case_match(Seasonality,
                                  "seas4mo" ~ "4 month seasonal profile",
                                  "seas6mo" ~ "6 month seasonal profile"),
         Access = case_match(Access,
                             "0.04" ~ "10% access to care",
                             "0.2412" ~ "50% access to care")) %>%
  mutate(targetLabel = factor(targetLabel, levels = unique(targetLabel)[order(unique(target))]))


# WRITE TO FILE ----

saveRDS(data, "./data_and_figures/supplement_fig27/data_fig27.rds")
