############################################################
# import_functions
#
# Helper functions for plotting scripts within M3TPP project workflow
#
# Written by Lydia Burgert
# Adapted by Lydia Braunack-Mayer
############################################################

source(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/2_postprocessing/postprocessing_resources.R"))

# ----------------------------------------------------------
# import_EIRs_cat
# ----------------------------------------------------------

import_EIRs_cat <- function(exp, scenario_id, seeds, timesteps){
  
  # ----------------------------------------------------------
  # This function averages continuous OpenMalaria outputs from multiple random seeds, by setting
  #
  # Inputs
  # exp: string containing the experiment name as it appears on SciCore, e.g. "MyExperiment"
  # scenario_id: character vector containing id(s) for scenario(s) to plot, e.g c("MyExperiment_1", "MyExperiment_2")
  # seeds: integer vector containing selection of seeds to plot, e.g. c(1)
  # timesteps: integer vector containing selection of OpenMalaria timesteps to plot, e.g. 1:73 captures outputs for the first year of simulation
  #
  # Outputs
  # out: list containing 'All', data frame containing all continuous OpenMalaria outputs for the given setting, and
  # 'Average', data frame containing averaged continuous OpenMalaria outputs for the given setting
  #
  # ----------------------------------------------------------
  
  # Load packages
  require(data.table)
  
  # Set up file pathways
  om_results_folder <- paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/om/")
  
  # Set up list to store function outputs
  out_all <- list()
  
  # Identify and concatenate continuous OpenMalaria outputs for the given setting, seed and selection of timesteps
  for (j in scenario_id) {
    for (i in seeds) {
      om_result <- read.table(paste0(om_results_folder, j, "_", i, "_cts.txt"), sep = "\t", header = TRUE)
      om_result$scenario_id <- j; om_result$seed <- i
      out_all[[paste0("scenario_id", j, "_seed", i)]]<- subset(om_result, timestep %in% timesteps)
    }
  }
  
  # Average continuous OpenMalaria outputs across seeds
  out_all <- rbindlist(out_all)
  out_average <- out_all[, lapply(.SD, mean), list(scenario_id, timestep)]
  
  # Return function outputs
  return(list("All" = out_all, "Average" = out_average))
  
}

# ----------------------------------------------------------
# import_monitoring_outcome
# ----------------------------------------------------------

import_monitoring_outcome <- function(exp, scenario_id, measure, seeds, timesteps){
  
  # ----------------------------------------------------------
  # This function averages OpenMalaria monitoring outcomes from multiple random seeds, by setting and age group
  #
  # Inputs
  # exp: string containing the experiment name as it appears on SciCore, e.g. "MyExperiment"
  # scenario_id: character vector containing id(s) for scenario(s) to plot, e.g c("MyExperiment_1", "MyExperiment_2")
  # measure: integer containing the id for the survey measure to plot (see https://github.com/SwissTPH/openmalaria/wiki/MonitoringOptions)
  # seeds: integer vector containing selection of seeds to plot, e.g. c(1)
  # timesteps: integer vector containing selection of OpenMalaria timesteps to plot, e.g. 1:73 captures outputs for the first year of simulation
  #
  # Outputs
  # out: list containing 'All', data frame containing the selected OpenMalaria monitoring outcome for all simulations for a given setting, and
  # 'Average', data frame containing averaged OpenMalaria monitoring outcome for the given setting
  #
  # ----------------------------------------------------------
  
  # Load packages
  require(data.table)
  
  # Set up file pathways
  om_results_folder <- paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/om/")
  
  # Set up list to store function outputs
  out_all <- list()
  
  # Identify and concatenate continuous OpenMalaria outputs for the given setting, seed and selection of timesteps
  for (j in scenario_id) {
    for (i in seeds) {
      om_result <- read.table(paste0(om_results_folder, j, "_", i, "_out.txt"), sep = "\t", header = FALSE)
      names(om_result) <- c("timestep", "agegroup", "measures", "outcome")
      om_result$scenario_id <- j; om_result$seed <- i
      out_all[[paste0("scenario_id", j, "_seed", i)]] <- subset(om_result, timestep %in% timesteps & measures %in% measure)
    }
  }
  
  # Average continuous OpenMalaria outputs across seeds, by age group
  out_all <- rbindlist(out_all)
  out_average <- out_all[, lapply(.SD, mean), list(scenario_id, timestep, agegroup)]
  
  # Return function outputs
  return(list("All" = out_all, "Average" = out_average))
  
}

