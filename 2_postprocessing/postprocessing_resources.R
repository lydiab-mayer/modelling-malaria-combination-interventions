################################
### STEP 2: POST-PROCESSING  ###
################################

# -------------------------------------------------------------------------------------------------------------
#
# Helper functions for running post-processing of OM simulations to aggregate data which will be used to train GP
# 
# Original script:
# Created 14.10.2021
# lydia.braunack-mayer@swisstph.ch 
#
# R version 4.1.0
#
# -------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------
# SET UP
# -------------------------------------------------------------------------------------------------------------

# Import packages
library(dplyr)
library(tidyr)
library(stringr)

options(dplyr.summarise.inform = FALSE)


# -------------------------------------------------------------------------------------------------------------
# DEFINE HANDY HELPER FUNCTIONS REFERRED TO WITHIN KEY CALCULATIONS
# -------------------------------------------------------------------------------------------------------------

extract.agegroups <- function(path, ...) {
  # Function to extract age groups from an OpenMalaria xml
  #
  # Inputs: 
  #   path, a file pathway to an OpenMalaria xml
  #   ..., additional arguments passed to readLines()
  #
  # Outputs: a vector containing the age groups contained within the xml
  
  #require(stringr)
  
  file <- readLines(path, ...)
  
  oline <- grep("<ageGroup lowerbound=\"0\">", file)[2]
  cline <- grep("</ageGroup>", file)[2]
  
  groups <- str_extract_all(file[(oline + 1):(cline - 1)], "[0-9.]+")
  
  out <- c(0)
  
  for(i in 1:length(groups)) {
    out <- c(out, as.numeric(groups[[i]]))
  }
  
  return(out)
}


time.to.date <- function(t, time.step, date) {
  # Function to convert the time step in an OpenMalaria outputs file to a date
  #
  # Inputs: 
  # t: the time steps to convert to date format, a numeric or numeric vector
  # time.step: the time step used in OpenMalaria, a numeric in days. Usually 1 or 5
  # date: the starting date of your OpenMalaria survey period in the format "%Y%m%d"
  #
  # Outputs: a vector containing the dates corresponding to each time step
  
  # Format inputs
  date <- as.Date(date, format = "%Y-%m-%d")
  
  # Create reference sequence of dates assuming no leap years
  ref <- seq(from = date,
             to = date + max(t)*time.step + ceiling(max(t)*time.step/365/4),
             by = "day")
  ref <- ref[strftime(ref, "%m-%d") != "02-29"]
  
  # Return sequence of dates incremented by time.step
  out <- ref[t*time.step]

  # Return outputs
  return(out)
  
}
#time.to.date(t = 1:1095, time.step = 5, date = "2030-01-01")


# -------------------------------------------------------------------------------------------------------------
# DEFINE HELPER FUNCTIONS TO CALCULATE INCIDENCE AND PREVALENCE REDUCTIONS FOR A SINGLE SIMULATION
# -------------------------------------------------------------------------------------------------------------

# # Sample arguments, retained here for testing
# om.result <- read.table("/scicore/home/penny/GROUP/M3TPP/obj6_test/om/obj6_test_1_1_out.txt", header = FALSE)
# measure <- 14
# age.group <- 2:8
# time.step <- 5
# date <- "2030-01-01"
# prevalence <- FALSE

calculate.annual.outcome <- function(om.result, measure, age.group, time.step = 5, date, prevalence = FALSE){
  
  # Function to calculate incidence, prevalence or another outcome by year for a given age group
  #
  # Inputs: 
  # om.result: the outputs of a single OpenMalaria simulation
  # measure: the outcome measure
  # age.group: an integer or integer vector containing the age groups of interest
  # time.step: the time step used in OpenMalaria, a numeric in days. Usually 1 or 5
  # date: the starting date of your OpenMalaria survey period in the format "yyyy-mm-dd"
  # prevalence: TRUE/FALSE indicating if prevalence vs. incidence or another metric should be outputted
  #
  # Outputs: data frame containing the outcome measure divided by the population size of the chosen age group
  # per month
  
  # Load required packages
  #require(dplyr)
  #require(tidyr)
  
  # Format OpenMalaria output file
  colnames(om.result) <- c("time", "age_group", "measure", "value")
  
  # Error monitoring
  if (!(measure %in% om.result$measure)) {
    print(paste0("Data for measure ", measure, " could not be found. Check that this measure has been included in the SurveyOptions section of your xml."))
  }
  
  # Translate from OpenMalaria 5-day time steps to years
  om.result$date <- time.to.date(om.result$time, time.step = time.step, date = date)
  om.result$year <- as.numeric(format(om.result$date, "%Y"))
  
  # Remove first time step from OpenMalaria outputs
  om.result <- om.result[om.result$time != 1, ]
  
  # Remove values for age groups other than those specified
  om.result <- om.result[om.result$age_group %in% age.group, ]
  
  # Remove measures other than that population size and the specified outcome measure
  om.result <- om.result[om.result$measure %in% c(0, measure), ]
  
  # Summarise all measures by summing up over age groups
  om.result <- om.result[, -which(names(om.result) %in% c("age_group"))] %>%
    group_by(measure, time, date, year) %>%
    summarise(value = sum(value))
  
  if (prevalence) {
    
    # Summarise further by averaging over all measures by year
    om.result <- om.result[, -which(names(om.result) %in% c("time", "date"))] %>%
      group_by(measure, year) %>%
      summarise(value = mean(value))
    
  } else {
    
    # Summarise further by summing over outcome measure by year
    om.pop <- om.result[om.result$measure == 0, -which(names(om.result) %in% c("time", "date"))] %>%
      group_by(measure, year) %>%
      summarise(value = mean(value))
    
    om.measure <- om.result[om.result$measure == measure, -which(names(om.result) %in% c("time", "date"))] %>%
      group_by(measure, year) %>%
      summarise(value = sum(value))
    
    om.result <- rbind(om.pop, om.measure)
    
  }
  
  # Transform to long format
  om.result <- pivot_wider(om.result, 
                           id_cols = c(year),
                           names_from = measure,
                           values_from = value,
                           names_prefix = "measure")
  om.result <- as.data.frame(om.result)
  
  # Rename columns
  names(om.result) <- c("year", "npop", "measure")
  
  # # Calculate outcome measure divided by population size
  # om.result$value <- om.result[, paste0("measure", measure)] / om.result[, "measure0"]
  
  # Order resulting data frame
  om.result <- om.result[order(om.result$year), ]
  rownames(om.result) <- NULL
  
  return(om.result)
  
}


# # Sample arguments, retained here for testing
# om.outcome <- calculate.annual.outcome(om.result = read.table("/scicore/home/penny/GROUP/M3TPP/obj6_test/om/obj6_test_1_1_out.txt", header = FALSE),
#                                        measure = 14,
#                                        age.group = 2:8,
#                                        time.step = 5,
#                                        date = "2030-01-01",
#                                        prevalence = FALSE)
# id <- "IncidenceCPPGain"
# year.counterfactual <- 2039
# year.intervention <- 2044
# prevalence <- FALSE

calculate.annual.gain <- function(om.outcome, id, year.counterfactual, year.intervention, prevalence = FALSE) {
  
  # Function to calculate gain in incidence, prevalence or another outcome by year
  #
  # Inputs: 
  # om.outcome: the outputs the function calculate.annual.outcome
  # id: a string that will be used to name columns of the function outputs, e.g. "IncidenceCPPGain"
  # year.counterfactual: the baseline year(s), as integer vector
  # year.intervention: the year(s) in which the intervention occurs , as integer vector
  # prevalence: TRUE/FALSE indicating if prevalence vs. incidence or another metric should be outputted
  #
  # Outputs: data frame containing a single row with the reduction per month
  
  #require(tidyr)
  
  # Set up data
  om.outcome <- om.outcome[om.outcome$year %in% c(year.counterfactual, year.intervention), ]
  om.outcome$year[om.outcome$year %in% year.counterfactual] <- "counterfactual"
  om.outcome$year[om.outcome$year %in% year.intervention] <- "intervention"
  
  # Summarise outcomes
  if (prevalence) {
    
    # Summarise with mean 
    om.outcome <- om.outcome %>%
      group_by(year) %>%
      summarise(value = mean(measure / npop))
    
  } else {
    
    # Summarise with sum
    om.outcome <- om.outcome %>%
      group_by(year) %>%
      summarise(value = sum(measure) / sum(npop))
    
  }
  
  # Calculate gain
  om.outcome <- pivot_wider(om.outcome,
                            names_from = year,
                            values_from = value)
  om.outcome$gain <- om.outcome$counterfactual - om.outcome$intervention
  
  # Format outputs
  om.outcome <- om.outcome[, 3]
  names(om.outcome) <- id
  
  # Return outputs
  return(om.outcome)
  
}

# # Sample arguments, retained here for testing
# dir <- "/scicore/home/penny/GROUP/M3TPP/obj6_test/"
# param.file <- "/scicore/home/penny/GROUP/M3TPP/obj6_test/postprocessing/split/obj6test_seas4mo_Mali_16_5_weibull_0.1227.txt"
# param.table <- read.table(param.file, sep = "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
# scenario.params <- param.table[1, ]
# om.file <- paste(dir, "om/", param.table[1, ]$Scenario_Name, "_", param.table[1, ]$SeedLabel, "_out.txt", sep = "")
# om.result <- read.table(om.file, sep = "\t")
# date <- "2030-01-01"
# year.counterfactual <- 2039
# year.intervention <- 2044
# min.int <- 0.25
#
# report.results(dir, om.result, date, year.counterfactual, year.intervention, min.int, scenario.params)

report.results <- function(dir, om.result, date, year.counterfactual, year.intervention, min.int, scenario.params) {
  
  # Define age groups 
  age.groups <- extract.agegroups(paste0(dir, "scaffold.xml"), warn = FALSE) # All age groups
  age.int <- seq(which(age.groups == min.int), as.numeric(scenario.params["maxGroup"])) # Intervention age group
  # age.int <- seq(which(age.groups == min.int), which(age.groups == 10) - 1) # Children X to 10 years old
  age.210 <- seq(which(age.groups == 2), which(age.groups == 10) - 1) # Children 2 to 10 years old
  #age.05 <- seq(which(age.groups == 0), which(age.groups == 5) - 1) # Children 0 to 5 years old
  
  # Calculate annual prevalence in children 2 to 10 years old
  om.outcome <- calculate.annual.outcome(om.result = om.result, measure = 3, age.group = age.210, time.step = 5, date = date, prevalence = TRUE)
  prev.210 <- om.outcome[om.outcome$year == year.counterfactual, "measure"] / om.outcome[om.outcome$year == year.counterfactual, "npop"]
  names(prev.210) <- paste0("AnnualPrev210.", year.counterfactual)
  rm(om.outcome)

  # Calculate gain in clinical cases per person
  om.outcome <- calculate.annual.outcome(om.result = om.result, measure = 14, age.group = age.int, time.step = 5, date = date)
  inc.gain <- calculate.annual.gain(om.outcome = om.outcome, id = "IncidenceCPPGain", year.counterfactual = year.counterfactual, year.intervention = year.intervention)
  rm(om.outcome)
  
  # Calculate gain in severe cases per person
  om.outcome <- calculate.annual.outcome(om.result = om.result, measure = 78, age.group = age.int, time.step = 5, date = date)
  sev.gain <- calculate.annual.gain(om.outcome = om.outcome, id = "SevereCPPGain", year.counterfactual = year.counterfactual, year.intervention = year.intervention)
  rm(om.outcome) 
  
  # Calculate gain in mortality per person
  om.outcome <- calculate.annual.outcome(om.result = om.result, measure = 74, age.group = age.int, time.step = 5, date = date)
  mor.gain <- calculate.annual.gain(om.outcome = om.outcome, id = "DeathsPPGain", year.counterfactual = year.counterfactual, year.intervention = year.intervention)
  rm(om.outcome)
  
  # Format outputs
  out <- cbind.data.frame(scenario.params$Scenario_Name, scenario.params$SeedLabel, prev.210, inc.gain, sev.gain, mor.gain)
  colnames(out) <- c("Scenario_Name", "seed", names(prev.210), names(inc.gain), names(sev.gain), names(mor.gain))
  rownames(out) <- NULL
  
  # Return outputs
  return(out)
  
}


# -------------------------------------------------------------------------------------------------------------
# DEFINE WRAPPER FUNCTIONS TO CALCULATE INCIDENCE AND PREVALENCE REDUCTIONS FOR MULTIPLE SIMULATIONS
# -------------------------------------------------------------------------------------------------------------

# # Sample arguments, retained here for testing
# dir <- "/scicore/home/penny/GROUP/M3TPP/obj6_test/om"
# param.file <- "/scicore/home/penny/GROUP/M3TPP/obj6_test/postprocessing/split/obj6test_seas4mo_Mali_16_5_weibull_0.1227.txt"
# date <- "2030-01-01"
# year.counterfactual <- 2039
# year.intervention <- 2044
# min.int <- 0.25

postprocess.om <- function(dir, param.file, date, fmonth, months, year.counterfactual, year.intervention, min.int) {
  
  dir <- paste0(dirname(dir), "/")
  dest.agg <- paste0(dir, "postprocessing/agg_", basename(param.file))
  dest.seed <- paste0(dir, "postprocessing/seeds_", basename(param.file))
  
  # Set up function
  param.table <- read.table(param.file, sep = "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  om.outcome <- NULL
  
  # For each row in param.table, process the corresponding OpenMalaria simulation
  for(i in 1:nrow(param.table)) {
    
    skip <- FALSE
    
    print(i)
    om.file <- paste(dir, "om/", param.table[i, ]$Scenario_Name, "_", param.table[i, ]$SeedLabel, "_out.txt", sep = "")
    
    tryCatch(if(file.exists(om.file) & file.info(om.file)$size > 0) {
      
      # Read in file
      om.result <- read.table(om.file, sep = "\t")
      
      out <- report.results(dir = dir, 
                            om.result = om.result,
                            date = date,
                            year.counterfactual =  year.counterfactual,
                            year.intervention = year.intervention,
                            min.int = min.int, 
                            scenario.params = param.table[i, ])

      om.outcome <- data.frame(rbind(om.outcome, out), stringsAsFactors = FALSE)}, 
      error = function(e) {skip <<- TRUE})
    
    if (skip) {next}
  }
  
  # Summarize results over each seed
  om.agg <- om.outcome %>% 
    group_by(Scenario_Name) %>% 
    summarise_at(c(names(om.outcome)[(which(names(om.outcome) == "seed") + 1):length(names(om.outcome))]), median, na.rm = TRUE)
  
  # Prepare results
  no.seed.tab <- unique(param.table[, -c(which(colnames(param.table) %in% c("SEED", "SeedLabel")))])
  seed.tab <- merge(no.seed.tab, om.outcome, by = c("Scenario_Name"))
  tab <- merge(no.seed.tab, om.agg, by = c("Scenario_Name"))
  
  # Write result tables (summarized and with seeds) to files
  write.table(tab, dest.agg, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(seed.tab, dest.seed, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}
