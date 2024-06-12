##########################
# Auxiliary functions for postprocessing OpenMalaria results
#
#######################################
#######################################
###                                 ###
### POSTPROCESSING RESOURCES SCRIPT ###
###                                 ###
#######################################
#######################################

### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Creates functions to calculate outcomes of OM results per simulation
### Produces both timeseries per measure and relative reductions
###
### adapted from:
### original script by: monica.golumbeanu@unibas.ch | 02.10.2019
### modified script by: lydia.burgert@unibas.ch 
###
### created by
### narimane.nekkab@unibas.ch
### 06.12.2021
###
### -------------------------------------------------------------------------

###############
### library ###
###############

require(dplyr)

options(dplyr.summarise.inform = FALSE)

########################
### Helper functions ###
########################

# Created by lydia.branauck-mayer@unibas.ch

############################
### Age groups extration ###

extract.agegroups <- function(path) {
  # Function to extract age groups from an OpenMalaria xml
  # Inputs: path, a file pathway to an OpenMalaria xml
  # Outputs: a vector containing the age groups contained within the xml
  
  require(stringr)
  
  file <- readLines(path)
  
  oline <- grep("<ageGroup lowerbound=\"0\">", file)[2]
  cline <- grep("</ageGroup>", file)[2]
  groups <- str_extract_all(file[(oline + 1):(cline - 1)], "[0-9.]+")
  
  out <- c(0)
  for(i in 1:length(groups)) {
    out <- c(out, as.numeric(groups[[i]]))
  }
  
  return(out)
}

#############################
### Time to date function ###

# time.to.date <- function(t, time.step, date) {
#   # Function to convert the time step in an OpenMalaria outputs file to a date
#   #
#   # Inputs: 
#   # t: the time steps to convert to date format, a numeric or numeric vector
#   # time.step: the time step used in OpenMalaria, a numeric in days. Usually 1 or 5
#   # date: the starting date of your OpenMalaria survey period in the format "%Y%m%d"
#   #
#   # Outputs: a vector containing the dates corresponding to each time step
#   
#   # Format inputs
#   date <- as.Date(date, format = "%Y-%m-%d")
#   
#   # Create reference sequence of dates for a non-leap year (OpenMalaria assumes 365 days in a year)
#   date.ref <- seq(from = as.Date("2021-01-01", format = "%Y-%m-%d"),
#                   to = as.Date("2021-12-31", format = "%Y-%m-%d"),
#                   by = "day")
#   date.ref <- rep(strftime(date.ref, "%m-%d"), 
#                   ceiling(max(t)*time.step/365) + 1)
#   year.ref <- rep(as.numeric(strftime(date, "%Y")):(as.numeric(strftime(date, "%Y")) + ceiling(max(t)*time.step/365)),
#                   each = 365)
#   ref <- paste0(year.ref, "-", date.ref)
#   
#   # Match reference sequence of dates to the start date of your OpenMalaria survey period
#   start <- which(date.ref == strftime(date, "%m-%d"))[1]
#   out <- ref[start - 1 + t*time.step]
#   
#   # Return outputs
#   return(out)
#   
# }

##########################
### OUTCOMES FUNCTIONS ###
##########################

# # Testing
# param_table = read.table(paste0("/scicore/home/penny/GROUP/M3TPP/iTPP1a_postprocessing_test/param_tab.txt"),header = T, stringsAsFactors = F)
# OM_result_file = "/scicore/home/penny/GROUP/M3TPP/iTPP1a_postprocessing_test/om/iTPP1a_postprocessing_test_1_1_out.txt"
# om_result = read.table(OM_result_file, sep="\t")
# scenario_params = param_table[which(param_table$Scenario_Name == "iTPP1a_postprocessing_test_1" & param_table$SEED == 1),]
# results_folder = "/scicore/home/penny/GROUP/M3TPP/iTPP1a_postprocessing_test/om/"
# measure_value = 1
# measure_name = "all"


#############################
### DEMOGRAPHICS FUNCTION ###

getDemographics <- function(om_result, scenario_params){
  
  ###################
  #### AGE GROUP ####
  
  # define age groups 
  age_groups <- c(0.165,0.25,0.5,seq(1,10,1),15,20,50,100)
  # age_groups <- c(0.25,0.5,1,2,3,4,5,10,20,100)
  minIntAge=0.25
  age210 = seq(which(age_groups==2)+1,which(age_groups==10))
  ageint = seq(which(age_groups==minIntAge)+1,as.numeric(scenario_params["maxGroup"]))
  age05 = seq(1,which(age_groups==5))
  
  
  ################
  #### LABELS ####
  
  colnames(om_result) = c("time", "age_group", "measure", "value")
  
  # #########################
  # #### REMOVE MEASURES ####
  # 
  # # # Remove measures without age group
  # measures_to_remove = c(7, 21, 26, 31:36)
  # om_result = om_result[-which(om_result$measure %in% measures_to_remove),]
  
  #################################
  #### REMOVE 1ST OBSERVATIONS ####
  
  # Remove first measurement as it includes all cases until then
  to_remove = which(om_result$time %in% c(1,2))
  om_result = om_result[-to_remove,]
  
  #############################
  #### CREATE NEW TIMESTEP ####
  
  # Timestep 0 starts at startDate so need to remove 2 start dates
  om_result$timestep = om_result$time - 1
  
  # Get date
  # om_result$date = time.to.date(om_result$timestep, 5, "2030-01-01")
  # om_result$date = time.to.date(om_result$timestep, 5, startDate)
  
  ##############################
  #### DEPLOYMENT TIMESTEPS ####
  
  # Get deployment timesteps
  massVaccination = om_result[which(om_result$measure == 22),] %>% 
    # mutate(timestep = factor(time)) %>% 
    group_by(timestep) %>% 
    dplyr::summarise(total = sum(value)) %>% 
    ungroup() %>% 
    filter(total > 0)
  
  # Deployment timesteps
  DeploymentTimesteps = massVaccination$timestep
  
  #############################
  #### REFERENCE TIMESTEPS ####
  
  # List of reference timesteps to evaluate BEFORE: by year
  timesteps_5years_beforeInt_12mo = seq(DeploymentTimesteps[1]-5*73, DeploymentTimesteps[1]-1, 1)
  
  # List of reference timesteps to evaluate BEFORE: by 6 months
  timesteps_5years_beforeInt_6mo = c(seq(DeploymentTimesteps[1]-5*73, DeploymentTimesteps[1]-1-5*73+36, 1),
                                     seq(DeploymentTimesteps[1]-4*73, DeploymentTimesteps[1]-1-4*73+36, 1),
                                     seq(DeploymentTimesteps[1]-3*73, DeploymentTimesteps[1]-1-3*73+36, 1),
                                     seq(DeploymentTimesteps[1]-2*73, DeploymentTimesteps[1]-1-2*73+36, 1),
                                     seq(DeploymentTimesteps[1]-1*73, DeploymentTimesteps[1]-1-1*73+36, 1))
  
  # List of reference timesteps to evaluate AFTER: by year
  timesteps_y1_12mo = seq(DeploymentTimesteps[1], DeploymentTimesteps[2]-1, 1)
  timesteps_y5_12mo = seq(DeploymentTimesteps[5], DeploymentTimesteps[6]-1, 1)
  timesteps_y10_12mo = seq(DeploymentTimesteps[10], DeploymentTimesteps[11]-1, 1)
  
  # List of reference timesteps to evaluate AFTER: by 6 months
  timesteps_y1_6mo = seq(DeploymentTimesteps[1], DeploymentTimesteps[1]-1+36, 1)
  timesteps_y5_6mo = seq(DeploymentTimesteps[5], DeploymentTimesteps[5]-1+36, 1)
  timesteps_y10_6mo = seq(DeploymentTimesteps[10], DeploymentTimesteps[10]-1+36, 1)
  
  ############################
  #### ADD ADJUSTED YEAR  ####
  
  # Reference
  reference = seq(DeploymentTimesteps[1]-5*73,max(om_result$timestep),73)
  
  # Add adjusted year (based on reference timesteps and 1st deployment)
  # om_result$year = c(rep(0,nrow(om_result[which(om_result$timestep < reference[1]),])),
  #                    rep(seq(1,length(reference)-1), each =  length(unique(om_result$age_group))*length(unique(om_result$measure))*73),
  #                    rep(length(reference),nrow(om_result[which(om_result$timestep >= reference[length(reference)]),])))
  om_result$year = 16
  om_result$year = ifelse(om_result$timestep < reference[1], 0, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[1] & om_result$timestep < reference[2], 1, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[2] & om_result$timestep < reference[3], 2, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[3] & om_result$timestep < reference[4], 3, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[4] & om_result$timestep < reference[5], 4, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[5] & om_result$timestep < reference[6], 5, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[6] & om_result$timestep < reference[7], 6, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[7] & om_result$timestep < reference[8], 7, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[8] & om_result$timestep < reference[9], 8, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[9] & om_result$timestep < reference[10], 9, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[10] & om_result$timestep < reference[11], 10, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[11] & om_result$timestep < reference[12], 11, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[12] & om_result$timestep < reference[13], 12, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[13] & om_result$timestep < reference[14], 13, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[14] & om_result$timestep < reference[15], 14, om_result$year)
  om_result$year = ifelse(om_result$timestep >= reference[15] & om_result$timestep < reference[16], 15, om_result$year)
  
  
  # # Add 6 months (based on reference timesteps and 1st deployment)
  om_result = om_result %>%
    dplyr::mutate(sixmonths = NA,
                  sixmonths = ifelse(timestep < reference[1], 0, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-5*73, DeploymentTimesteps[1]-1-5*73+36, 1), 1, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-4*73, DeploymentTimesteps[1]-1-4*73+36, 1), 2, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-3*73, DeploymentTimesteps[1]-1-3*73+36, 1), 3, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-2*73, DeploymentTimesteps[1]-1-2*73+36, 1), 4, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-1*73, DeploymentTimesteps[1]-1-1*73+36, 1), 5, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1], DeploymentTimesteps[1]-1+36, 1), 6, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[2], DeploymentTimesteps[2]-1+36, 1), 7, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[3], DeploymentTimesteps[3]-1+36, 1), 8, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[4], DeploymentTimesteps[4]-1+36, 1), 9, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[5], DeploymentTimesteps[5]-1+36, 1), 10, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[6], DeploymentTimesteps[6]-1+36, 1), 11, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[7], DeploymentTimesteps[7]-1+36, 1), 12, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[8], DeploymentTimesteps[8]-1+36, 1), 13, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[9], DeploymentTimesteps[9]-1+36, 1), 14, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[10], DeploymentTimesteps[10]-1+36, 1), 15, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[11], DeploymentTimesteps[11]-1+36, 1), 16, sixmonths))
  
  
  ##########################
  #### POPULATION: YEAR ####
  
  # population per time step per age group
  total_pop = as.data.frame(om_result[om_result$measure == 0, ])
  
  # calculate total population per age group per year
  total_pop_age_yearly = data.frame(total_pop %>% dplyr::group_by(age_group, year) %>% dplyr::summarise(n = mean(value))) %>% ungroup()
  
  # sum up total population per year
  total_pop_all_yearly = total_pop %>% dplyr::group_by(timestep, year) %>% dplyr::summarise(sum = sum(value)) %>% dplyr::group_by(year) %>% dplyr::summarise(n = mean(sum)) %>% ungroup()
  
  # sum up population intervention age-groups over the years
  pop_int_yearly <- total_pop_age_yearly[total_pop_age_yearly$age_group %in% ageint,] %>% dplyr::group_by(year) %>% dplyr::summarise(n = sum(n)) %>% ungroup()
  pop_210_yearly <- total_pop_age_yearly[total_pop_age_yearly$age_group %in% age210,] %>% dplyr::group_by(year) %>% dplyr::summarise(n = sum(n)) %>% ungroup()
  pop_05_yearly <- total_pop_age_yearly[total_pop_age_yearly$age_group %in% age05,] %>% dplyr::group_by(year) %>% dplyr::summarise(n = sum(n)) %>% ungroup()
  
  # mean intervention population size
  meanpopint_yearly <- mean(pop_int_yearly$n)
  
  ###############################
  #### POPULATION: SIX MONTH ####
  
  # calculate total population per age group per year
  total_pop_age_sixmonthly = total_pop[complete.cases(total_pop),] %>% dplyr::group_by(age_group, sixmonths) %>% dplyr::summarise(n = mean(value)) %>% ungroup()
  
  # sum up total population per year
  total_pop_all_sixmonthly = total_pop[complete.cases(total_pop),] %>% dplyr::group_by(timestep, sixmonths) %>% dplyr::summarise(sum = sum(value)) %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(n = mean(sum)) %>% ungroup()
  
  # sum up population intervention age-groups over the years
  pop_int_sixmonthly <- total_pop_age_sixmonthly[total_pop_age_sixmonthly$age_group %in% ageint,] %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(n = sum(n)) %>% ungroup()
  pop_210_sixmonthly <- total_pop_age_sixmonthly[total_pop_age_sixmonthly$age_group %in% age210,] %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(n = sum(n)) %>% ungroup()
  pop_05_sixmonthly <- total_pop_age_sixmonthly[total_pop_age_sixmonthly$age_group %in% age05,] %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(n = sum(n)) %>% ungroup()
  
  # mean intervention population size
  meanpopint_sixmonthly <- mean(pop_int_sixmonthly$n)
  
  
  ################
  #### TOTALS ####
  ################
  
  # Summarize all measures results by summing up over age groups
  agg_om_result_total_yearly = om_result[,-which(names(om_result)=="age_group")] %>% dplyr::group_by(timestep, year, measure) %>% dplyr::summarise(value = sum(value)) %>% ungroup()
  
  # Summarize all measures results by summing up over age groups
  agg_om_result_total_sixmonthly = om_result[complete.cases(om_result),-which(names(om_result)=="age_group")] %>% dplyr::group_by(timestep, sixmonths, measure) %>% dplyr::summarise(value = sum(value)) %>% ungroup()
  
  ######################
  ### Return results ###
  
  results = list(om_result, 
                 agg_om_result_total_yearly, 
                 agg_om_result_total_sixmonthly, 
                 total_pop_age_yearly, 
                 total_pop_age_sixmonthly,
                 total_pop_all_yearly, 
                 total_pop_all_sixmonthly,
                 pop_int_yearly, 
                 pop_int_sixmonthly,
                 pop_05_yearly, 
                 pop_05_sixmonthly,
                 pop_210_yearly, 
                 pop_210_sixmonthly,
                 ageint, 
                 age05, 
                 age210)
  
  return(results)
  
}

###########################
### PREVALENCE FUNCTION ###

getPrevalenceOutcomes <- function(om_result, 
                                  agg_om_result_total_yearly, agg_om_result_total_sixmonthly, 
                                  total_pop_age_yearly, total_pop_age_sixmonthly,
                                  total_pop_all_yearly, total_pop_all_sixmonthly,
                                  ageint, age05, age210,
                                  measure_value, measure_name, scenario_params){
  
  ####################################
  ############ PREVALENCE ############
  ####################################
  
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  
  # number of infected per age-group per time-step
  n_Inf = as.data.frame(om_result[om_result$measure == measure_value,]) # The number of human hosts with an infection (patent or not) on the reporting timestep
  
  # total number of infected per time-step
  n_Inf_total_yearly = as.data.frame(agg_om_result_total_yearly[agg_om_result_total_yearly$measure == measure_value,])
  n_Inf_total_sixmonthly = as.data.frame(agg_om_result_total_sixmonthly[agg_om_result_total_sixmonthly$measure == measure_value,])
  
  # add the population size in respective age-groups
  n_Inf_yearly = merge(total_pop_age_yearly, n_Inf, by=c("age_group","year"), all = F)
  n_Inf_total_yearly = merge(total_pop_all_yearly, n_Inf_total_yearly, by=c("year"), all = F)
  n_Inf_sixmonthly = merge(total_pop_age_sixmonthly, n_Inf, by=c("age_group","sixmonths"), all = F)
  n_Inf_total_sixmonthly = merge(total_pop_all_sixmonthly, n_Inf_total_sixmonthly, by=c("sixmonths"), all = F)
  
  # divide number of infected by respective age-group
  prev_Inf_agegroups_cont_yearly = n_Inf_yearly %>% mutate(prev = value/n)
  prev_Inf_allages_cont_yearly = n_Inf_total_yearly %>% mutate(prev = value/n)
  prev_Inf_agegroups_cont_sixmonthly = n_Inf_sixmonthly %>% mutate(prev = value/n)
  prev_Inf_allages_cont_sixmonthly = n_Inf_total_sixmonthly %>% mutate(prev = value/n)
  
  # Intervention group
  prev_Inf_int_cont_yearly = prev_Inf_agegroups_cont_yearly %>% 
    dplyr::filter(age_group %in% ageint) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  prev_Inf_int_cont_sixmonthly = prev_Inf_agegroups_cont_sixmonthly %>% 
    dplyr::filter(age_group %in% ageint) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  
  # Children 2-10
  prev_Inf_210_cont_yearly = prev_Inf_agegroups_cont_yearly %>% 
    dplyr::filter(age_group %in% age210) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  prev_Inf_210_cont_sixmonthly = prev_Inf_agegroups_cont_sixmonthly %>% 
    dplyr::filter(age_group %in% age210) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  
  # Children 0-5
  prev_Inf_05_cont_yearly = prev_Inf_agegroups_cont_yearly %>% 
    dplyr::filter(age_group %in% age05) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  prev_Inf_05_cont_sixmonthly = prev_Inf_agegroups_cont_sixmonthly %>% 
    dplyr::filter(age_group %in% age05) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  
  
  ##############################################
  #### PREVALENCE REDUCTION: ALL INFECTIONS ####
  ##############################################
  
  ################
  ### All ages ###
  
  # yearly average prevalence
  prev_Inf_allages_yearly = prev_Inf_allages_cont_yearly %>% dplyr::group_by(year) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  prev_Inf_allages_sixmonthly = prev_Inf_allages_cont_sixmonthly %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  
  # prevalence before and after
  prev_Inf_allages_before_yearly=mean(as.numeric(prev_Inf_allages_yearly$avg[which(prev_Inf_allages_yearly$year %in% c(1:5))]))
  prev_Inf_allages_followup1_yearly=as.numeric(prev_Inf_allages_yearly$avg[which(prev_Inf_allages_yearly$year==6)])
  prev_Inf_allages_followup5_yearly=as.numeric(prev_Inf_allages_yearly$avg[which(prev_Inf_allages_yearly$year==10)])
  prev_Inf_allages_followup10_yearly=as.numeric(prev_Inf_allages_yearly$avg[which(prev_Inf_allages_yearly$year==15)])
  # sixmonths
  prev_Inf_allages_before_sixmonthly=mean(as.numeric(prev_Inf_allages_sixmonthly$avg[which(prev_Inf_allages_sixmonthly$sixmonths %in% c(1:5))]))
  prev_Inf_allages_followup1_sixmonthly=as.numeric(prev_Inf_allages_sixmonthly$avg[which(prev_Inf_allages_sixmonthly$sixmonths==6)])
  prev_Inf_allages_followup5_sixmonthly=as.numeric(prev_Inf_allages_sixmonthly$avg[which(prev_Inf_allages_sixmonthly$sixmonths==10)])
  prev_Inf_allages_followup10_sixmonthly=as.numeric(prev_Inf_allages_sixmonthly$avg[which(prev_Inf_allages_sixmonthly$sixmonths==15)])
  
  # prevalence reduction
  prev_red1_Inf_allages_yearly = ((prev_Inf_allages_before_yearly - prev_Inf_allages_followup1_yearly)/prev_Inf_allages_before_yearly)*100
  prev_red5_Inf_allages_yearly = ((prev_Inf_allages_before_yearly - prev_Inf_allages_followup5_yearly)/prev_Inf_allages_before_yearly)*100
  prev_red10_Inf_allages_yearly = ((prev_Inf_allages_before_yearly - prev_Inf_allages_followup10_yearly)/prev_Inf_allages_before_yearly)*100
  # sixmonths
  prev_red1_Inf_allages_sixmonthly = ((prev_Inf_allages_before_sixmonthly - prev_Inf_allages_followup1_sixmonthly)/prev_Inf_allages_before_sixmonthly)*100
  prev_red5_Inf_allages_sixmonthly = ((prev_Inf_allages_before_sixmonthly - prev_Inf_allages_followup5_sixmonthly)/prev_Inf_allages_before_sixmonthly)*100
  prev_red10_Inf_allages_sixmonthly = ((prev_Inf_allages_before_sixmonthly - prev_Inf_allages_followup10_sixmonthly)/prev_Inf_allages_before_sixmonthly)*100
  
  ##########################
  ### Intervention group ###
  
  # yearly average prevalence
  prev_Inf_int_yearly = prev_Inf_int_cont_yearly %>% dplyr::group_by(year) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  prev_Inf_int_sixmonthly = prev_Inf_int_cont_sixmonthly %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  
  # prevalence before and after
  prev_Inf_int_before_yearly=mean(as.numeric(prev_Inf_int_yearly$avg[which(prev_Inf_int_yearly$year %in% c(1:5))]))
  prev_Inf_int_followup1_yearly=as.numeric(prev_Inf_int_yearly$avg[which(prev_Inf_int_yearly$year==6)])
  prev_Inf_int_followup5_yearly=as.numeric(prev_Inf_int_yearly$avg[which(prev_Inf_int_yearly$year==10)])
  prev_Inf_int_followup10_yearly=as.numeric(prev_Inf_int_yearly$avg[which(prev_Inf_int_yearly$year==15)])
  # sixmonths
  prev_Inf_int_before_sixmonthly=mean(as.numeric(prev_Inf_int_sixmonthly$avg[which(prev_Inf_int_sixmonthly$sixmonths %in% c(1:5))]))
  prev_Inf_int_followup1_sixmonthly=as.numeric(prev_Inf_int_sixmonthly$avg[which(prev_Inf_int_sixmonthly$sixmonths==6)])
  prev_Inf_int_followup5_sixmonthly=as.numeric(prev_Inf_int_sixmonthly$avg[which(prev_Inf_int_sixmonthly$sixmonths==10)])
  prev_Inf_int_followup10_sixmonthly=as.numeric(prev_Inf_int_sixmonthly$avg[which(prev_Inf_int_sixmonthly$sixmonths==15)])
  
  # prevalence reduction
  prev_red1_Inf_int_yearly = ((prev_Inf_int_before_yearly - prev_Inf_int_followup1_yearly)/prev_Inf_int_before_yearly)*100
  prev_red5_Inf_int_yearly = ((prev_Inf_int_before_yearly - prev_Inf_int_followup5_yearly)/prev_Inf_int_before_yearly)*100
  prev_red10_Inf_int_yearly = ((prev_Inf_int_before_yearly - prev_Inf_int_followup10_yearly)/prev_Inf_int_before_yearly)*100
  # sixmonths
  prev_red1_Inf_int_sixmonthly = ((prev_Inf_int_before_sixmonthly - prev_Inf_int_followup1_sixmonthly)/prev_Inf_int_before_sixmonthly)*100
  prev_red5_Inf_int_sixmonthly = ((prev_Inf_int_before_sixmonthly - prev_Inf_int_followup5_sixmonthly)/prev_Inf_int_before_sixmonthly)*100
  prev_red10_Inf_int_sixmonthly = ((prev_Inf_int_before_sixmonthly - prev_Inf_int_followup10_sixmonthly)/prev_Inf_int_before_sixmonthly)*100
  
  
  #################
  ### Ages 2-10 ###
  
  # yearly average prevalence
  prev_Inf_210_yearly = prev_Inf_210_cont_yearly %>% dplyr::group_by(year) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  prev_Inf_210_sixmonthly = prev_Inf_210_cont_sixmonthly %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  
  # prevalence before and after
  prev_Inf_210_before_yearly=mean(as.numeric(prev_Inf_210_yearly$avg[which(prev_Inf_210_yearly$year %in% c(1:5))]))
  prev_Inf_210_followup1_yearly=as.numeric(prev_Inf_210_yearly$avg[which(prev_Inf_210_yearly$year==6)])
  prev_Inf_210_followup5_yearly=as.numeric(prev_Inf_210_yearly$avg[which(prev_Inf_210_yearly$year==10)])
  prev_Inf_210_followup10_yearly=as.numeric(prev_Inf_210_yearly$avg[which(prev_Inf_210_yearly$year==15)])
  # sixmonths
  prev_Inf_210_before_sixmonthly=mean(as.numeric(prev_Inf_210_sixmonthly$avg[which(prev_Inf_210_sixmonthly$sixmonths %in% c(1:5))]))
  prev_Inf_210_followup1_sixmonthly=as.numeric(prev_Inf_210_sixmonthly$avg[which(prev_Inf_210_sixmonthly$sixmonths==6)])
  prev_Inf_210_followup5_sixmonthly=as.numeric(prev_Inf_210_sixmonthly$avg[which(prev_Inf_210_sixmonthly$sixmonths==10)])
  prev_Inf_210_followup10_sixmonthly=as.numeric(prev_Inf_210_sixmonthly$avg[which(prev_Inf_210_sixmonthly$sixmonths==15)])
  
  # prevalence reduction
  prev_red1_Inf_210_yearly = ((prev_Inf_210_before_yearly - prev_Inf_210_followup1_yearly)/prev_Inf_210_before_yearly)*100
  prev_red5_Inf_210_yearly = ((prev_Inf_210_before_yearly - prev_Inf_210_followup5_yearly)/prev_Inf_210_before_yearly)*100
  prev_red10_Inf_210_yearly = ((prev_Inf_210_before_yearly - prev_Inf_210_followup10_yearly)/prev_Inf_210_before_yearly)*100
  # sixmonths
  prev_red1_Inf_210_sixmonthly = ((prev_Inf_210_before_sixmonthly - prev_Inf_210_followup1_sixmonthly)/prev_Inf_210_before_sixmonthly)*100
  prev_red5_Inf_210_sixmonthly = ((prev_Inf_210_before_sixmonthly - prev_Inf_210_followup5_sixmonthly)/prev_Inf_210_before_sixmonthly)*100
  prev_red10_Inf_210_sixmonthly = ((prev_Inf_210_before_sixmonthly - prev_Inf_210_followup10_sixmonthly)/prev_Inf_210_before_sixmonthly)*100
  
  
  #################
  ### Ages 0-5 ###
  
  # yearly average prevalence
  prev_Inf_05_yearly = prev_Inf_05_cont_yearly %>% dplyr::group_by(year) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  prev_Inf_05_sixmonthly = prev_Inf_05_cont_sixmonthly %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  
  # prevalence before and after
  prev_Inf_05_before_yearly=mean(as.numeric(prev_Inf_05_yearly$avg[which(prev_Inf_05_yearly$year %in% c(1:5))]))
  prev_Inf_05_followup1_yearly=as.numeric(prev_Inf_05_yearly$avg[which(prev_Inf_05_yearly$year==6)])
  prev_Inf_05_followup5_yearly=as.numeric(prev_Inf_05_yearly$avg[which(prev_Inf_05_yearly$year==10)])
  prev_Inf_05_followup10_yearly=as.numeric(prev_Inf_05_yearly$avg[which(prev_Inf_05_yearly$year==15)])
  # sixmonths
  prev_Inf_05_before_sixmonthly=mean(as.numeric(prev_Inf_05_sixmonthly$avg[which(prev_Inf_05_sixmonthly$sixmonths %in% c(1:5))]))
  prev_Inf_05_followup1_sixmonthly=as.numeric(prev_Inf_05_sixmonthly$avg[which(prev_Inf_05_sixmonthly$sixmonths==6)])
  prev_Inf_05_followup5_sixmonthly=as.numeric(prev_Inf_05_sixmonthly$avg[which(prev_Inf_05_sixmonthly$sixmonths==10)])
  prev_Inf_05_followup10_sixmonthly=as.numeric(prev_Inf_05_sixmonthly$avg[which(prev_Inf_05_sixmonthly$sixmonths==15)])
  
  # prevalence reduction
  prev_red1_Inf_05_yearly = ((prev_Inf_05_before_yearly - prev_Inf_05_followup1_yearly)/prev_Inf_05_before_yearly)*100
  prev_red5_Inf_05_yearly = ((prev_Inf_05_before_yearly - prev_Inf_05_followup5_yearly)/prev_Inf_05_before_yearly)*100
  prev_red10_Inf_05_yearly = ((prev_Inf_05_before_yearly - prev_Inf_05_followup10_yearly)/prev_Inf_05_before_yearly)*100
  # sixmonths
  prev_red1_Inf_05_sixmonthly = ((prev_Inf_05_before_sixmonthly - prev_Inf_05_followup1_sixmonthly)/prev_Inf_05_before_sixmonthly)*100
  prev_red5_Inf_05_sixmonthly = ((prev_Inf_05_before_sixmonthly - prev_Inf_05_followup5_sixmonthly)/prev_Inf_05_before_sixmonthly)*100
  prev_red10_Inf_05_sixmonthly = ((prev_Inf_05_before_sixmonthly - prev_Inf_05_followup10_sixmonthly)/prev_Inf_05_before_sixmonthly)*100
  
  #####################################################
  ##################### OUTPUT ########################
  #####################################################
  
  ##################
  ### Reductions ###
  
  # Final row with outputs to return
  reductions = cbind.data.frame(scenario_params$Scenario_Name,
                                scenario_params$SEED,
                                prev_red1_Inf_allages_yearly,
                                prev_red5_Inf_allages_yearly,
                                prev_red10_Inf_allages_yearly,
                                prev_red1_Inf_int_yearly,
                                prev_red5_Inf_int_yearly,
                                prev_red10_Inf_int_yearly,
                                prev_red1_Inf_05_yearly,
                                prev_red5_Inf_05_yearly,
                                prev_red10_Inf_05_yearly,
                                prev_red1_Inf_210_yearly,
                                prev_red5_Inf_210_yearly,
                                prev_red10_Inf_210_yearly,
                                prev_red1_Inf_allages_sixmonthly,
                                prev_red5_Inf_allages_sixmonthly,
                                prev_red10_Inf_allages_sixmonthly,
                                prev_red1_Inf_int_sixmonthly,
                                prev_red5_Inf_int_sixmonthly,
                                prev_red10_Inf_int_sixmonthly,
                                prev_red1_Inf_05_sixmonthly,
                                prev_red5_Inf_05_sixmonthly,
                                prev_red10_Inf_05_sixmonthly,
                                prev_red1_Inf_210_sixmonthly,
                                prev_red5_Inf_210_sixmonthly,
                                prev_red10_Inf_210_sixmonthly)
  colnames(reductions) = c("Scenario_Name", "seed",
                           paste0("prev_red1_",measure_name,"Inf_allages_yearly"),
                           paste0("prev_red5_",measure_name,"Inf_allages_yearly"),
                           paste0("prev_red10_",measure_name,"Inf_allages_yearly"),
                           paste0("prev_red1_",measure_name,"Inf_int_yearly"),
                           paste0("prev_red5_",measure_name,"Inf_int_yearly"),
                           paste0("prev_red10_",measure_name,"Inf_int_yearly"),
                           paste0("prev_red1_",measure_name,"Inf_05_yearly"),
                           paste0("prev_red5_",measure_name,"Inf_05_yearly"),
                           paste0("prev_red10_",measure_name,"Inf_05_yearly"),
                           paste0("prev_red1_",measure_name,"Inf_210_yearly"),
                           paste0("prev_red5_",measure_name,"Inf_210_yearly"),
                           paste0("prev_red10_",measure_name,"Inf_210_yearly"),
                           paste0("prev_red1_",measure_name,"Inf_allages_sixmonthly"),
                           paste0("prev_red5_",measure_name,"Inf_allages_sixmonthly"),
                           paste0("prev_red10_",measure_name,"Inf_allages_sixmonthly"),
                           paste0("prev_red1_",measure_name,"Inf_int_sixmonthly"),
                           paste0("prev_red5_",measure_name,"Inf_int_sixmonthly"),
                           paste0("prev_red10_",measure_name,"Inf_int_sixmonthly"),
                           paste0("prev_red1_",measure_name,"Inf_05_sixmonthly"),
                           paste0("prev_red5_",measure_name,"Inf_05_sixmonthly"),
                           paste0("prev_red10_",measure_name,"Inf_05_sixmonthly"),
                           paste0("prev_red1_",measure_name,"Inf_210_sixmonthly"),
                           paste0("prev_red5_",measure_name,"Inf_210_sixmonthly"),
                           paste0("prev_red10_",measure_name,"Inf_210_sixmonthly"))
  
  ##################
  ### Timeseries ###
  
  timeseries = list(prev_Inf_agegroups_cont_yearly,
                    prev_Inf_allages_cont_yearly,
                    prev_Inf_int_cont_yearly,
                    prev_Inf_05_cont_yearly,
                    prev_Inf_210_cont_yearly,
                    prev_Inf_agegroups_cont_sixmonthly,
                    prev_Inf_allages_cont_sixmonthly,
                    prev_Inf_int_cont_sixmonthly,
                    prev_Inf_05_cont_sixmonthly,
                    prev_Inf_210_cont_sixmonthly)
  
  names(timeseries) = c(paste0("prev_",measure_name,"Inf_agegroups_yearly"),
                        paste0("prev_",measure_name,"Inf_allages_yearly"),
                        paste0("prev_",measure_name,"Inf_int_yearly"),
                        paste0("prev_",measure_name,"Inf_05_yearly"),
                        paste0("prev_",measure_name,"Inf_210_yearly"),
                        paste0("prev_",measure_name,"Inf_agegroups_sixmonthly"),
                        paste0("prev_",measure_name,"Inf_allages_sixmonthly"),
                        paste0("prev_",measure_name,"Inf_int_sixmonthly"),
                        paste0("prev_",measure_name,"Inf_05_sixmonthly"),
                        paste0("prev_",measure_name,"Inf_210_sixmonthly"))
  
  
  return(list(reductions, timeseries))
  
}

##########################
### INCIDENCE FUNCTION ###

getIncidenceOutcomes <- function(om_result, 
                                 total_pop_age_yearly, total_pop_age_sixmonthly,
                                 total_pop_all_yearly, total_pop_all_sixmonthly,
                                 pop_int_yearly, pop_int_sixmonthly,
                                 pop_05_yearly, pop_05_sixmonthly,
                                 pop_210_yearly, pop_210_sixmonthly,
                                 ageint, age05, age210,
                                 measure_value, measure_name, scenario_params){
  
  ###################################
  ############ INCIDENCE ############
  ###################################
  
  # Incidence = number of new cases over time
  
  # number of clinical cases
  nInc = as.data.frame(om_result[om_result$measure == measure_value,]) 
  
  # yearly average clinical cases per age-group
  Inc_agegroups_yearly = nInc %>% dplyr::group_by(age_group, year) %>% dplyr::summarise(sum = sum(value)) %>% dplyr::ungroup()
  Inc_agegroups_sixmonthly = nInc %>% dplyr::group_by(age_group, sixmonths) %>% dplyr::summarise(sum = sum(value)) %>% dplyr::ungroup()
  
  ################
  ### All ages ###
  
  # summed yearly clinical cases in intervention age-groups
  Inc_all_yearly = Inc_agegroups_yearly %>% dplyr::group_by(year) %>% dplyr::summarise(inc = sum(sum))
  Inc_all_yearly$cpp = Inc_all_yearly$inc / tail(total_pop_all_yearly$n, 1)
  Inc_all_sixmonthly = Inc_agegroups_sixmonthly %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(inc = sum(sum))
  Inc_all_sixmonthly$cpp = Inc_all_sixmonthly$inc / tail(total_pop_all_sixmonthly$n, 1)
  
  # Clinical incidence before and after
  Inc_all_before_yearly=mean(as.numeric(Inc_all_yearly$cpp[which(Inc_all_yearly$year %in% c(1:5))]))
  Inc_all_followup1_yearly=as.numeric(Inc_all_yearly$cpp[which(Inc_all_yearly$year==6)])
  Inc_all_followup5_yearly=as.numeric(Inc_all_yearly$cpp[which(Inc_all_yearly$year==10)])
  Inc_all_followup10_yearly=as.numeric(Inc_all_yearly$cpp[which(Inc_all_yearly$year==15)])
  # Clinical incidence before and after
  Inc_all_before_sixmonthly=mean(as.numeric(Inc_all_sixmonthly$cpp[which(Inc_all_sixmonthly$sixmonths %in% c(1:5))]))
  Inc_all_followup1_sixmonthly=as.numeric(Inc_all_sixmonthly$cpp[which(Inc_all_sixmonthly$sixmonths==6)])
  Inc_all_followup5_sixmonthly=as.numeric(Inc_all_sixmonthly$cpp[which(Inc_all_sixmonthly$sixmonths==10)])
  Inc_all_followup10_sixmonthly=as.numeric(Inc_all_sixmonthly$cpp[which(Inc_all_sixmonthly$sixmonths==15)])
  
  # Clinical incidence reduction
  Inc_red1_all_yearly = ((Inc_all_before_yearly - Inc_all_followup1_yearly)/Inc_all_before_yearly)*100
  Inc_red5_all_yearly = ((Inc_all_before_yearly - Inc_all_followup5_yearly)/Inc_all_before_yearly)*100
  Inc_red10_all_yearly = ((Inc_all_before_yearly - Inc_all_followup10_yearly)/Inc_all_before_yearly)*100
  # Clinical incidence reduction
  Inc_red1_all_sixmonthly = ((Inc_all_before_sixmonthly - Inc_all_followup1_sixmonthly)/Inc_all_before_sixmonthly)*100
  Inc_red5_all_sixmonthly = ((Inc_all_before_sixmonthly - Inc_all_followup5_sixmonthly)/Inc_all_before_sixmonthly)*100
  Inc_red10_all_sixmonthly = ((Inc_all_before_sixmonthly - Inc_all_followup10_sixmonthly)/Inc_all_before_sixmonthly)*100
  
  # # Clinical incidence averted
  # Inc_avert1_all_yearly = (Inc_all_before_yearly - Inc_all_followup1_yearly)*100000
  # Inc_avert5_all_yearly = (Inc_all_before_yearly - Inc_all_followup5_yearly)*100000
  # Inc_avert10_all_yearly = (Inc_all_before_yearly - Inc_all_followup10_yearly)*100000
  # # Clinical incidence averted
  # Inc_avert1_all_sixmonthly = (Inc_all_before_sixmonthly - Inc_all_followup1_sixmonthly)*100000
  # Inc_avert5_all_sixmonthly = (Inc_all_before_sixmonthly - Inc_all_followup5_sixmonthly)*100000
  # Inc_avert10_all_sixmonthly = (Inc_all_before_sixmonthly - Inc_all_followup10_sixmonthly)*100000
  
  ##########################
  ### Intervention group ###
  
  # summed yearly clinical cases in intervention age-groups
  Inc_int_yearly <- Inc_agegroups_yearly[Inc_agegroups_yearly$age_group %in% ageint,] %>% dplyr::group_by(year) %>% dplyr::summarise(inc = sum(sum)) %>% dplyr::ungroup()
  Inc_int_yearly$cpp <- Inc_int_yearly$inc / tail(pop_int_yearly$n, 1)
  Inc_int_sixmonthly <- Inc_agegroups_sixmonthly[Inc_agegroups_sixmonthly$age_group %in% ageint,] %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(inc = sum(sum)) %>% dplyr::ungroup()
  Inc_int_sixmonthly$cpp <- Inc_int_sixmonthly$inc / tail(pop_int_sixmonthly$n, 1)
  
  # Clinical incidence before and after
  Inc_int_before_yearly=mean(as.numeric(Inc_int_yearly$cpp[which(Inc_int_yearly$year %in% c(1:5))]))
  Inc_int_followup1_yearly=as.numeric(Inc_int_yearly$cpp[which(Inc_int_yearly$year==6)])
  Inc_int_followup5_yearly=as.numeric(Inc_int_yearly$cpp[which(Inc_int_yearly$year==10)])
  Inc_int_followup10_yearly=as.numeric(Inc_int_yearly$cpp[which(Inc_int_yearly$year==15)])
  # Clinical incidence before and after
  Inc_int_before_sixmonthly=mean(as.numeric(Inc_int_sixmonthly$cpp[which(Inc_int_sixmonthly$sixmonths %in% c(1:5))]))
  Inc_int_followup1_sixmonthly=as.numeric(Inc_int_sixmonthly$cpp[which(Inc_int_sixmonthly$sixmonths==6)])
  Inc_int_followup5_sixmonthly=as.numeric(Inc_int_sixmonthly$cpp[which(Inc_int_sixmonthly$sixmonths==10)])
  Inc_int_followup10_sixmonthly=as.numeric(Inc_int_sixmonthly$cpp[which(Inc_int_sixmonthly$sixmonths==15)])
  
  # Clinical incidence reduction
  Inc_red1_int_yearly = ((Inc_int_before_yearly - Inc_int_followup1_yearly)/Inc_int_before_yearly)*100
  Inc_red5_int_yearly = ((Inc_int_before_yearly - Inc_int_followup5_yearly)/Inc_int_before_yearly)*100
  Inc_red10_int_yearly = ((Inc_int_before_yearly - Inc_int_followup10_yearly)/Inc_int_before_yearly)*100
  # Clinical incidence reduction
  Inc_red1_int_sixmonthly = ((Inc_int_before_sixmonthly - Inc_int_followup1_sixmonthly)/Inc_int_before_sixmonthly)*100
  Inc_red5_int_sixmonthly = ((Inc_int_before_sixmonthly - Inc_int_followup5_sixmonthly)/Inc_int_before_sixmonthly)*100
  Inc_red10_int_sixmonthly = ((Inc_int_before_sixmonthly - Inc_int_followup10_sixmonthly)/Inc_int_before_sixmonthly)*100
  
  # Clinical incidence averted
  Inc_avert1_int_yearly = (Inc_int_before_yearly - Inc_int_followup1_yearly)*100000
  Inc_avert5_int_yearly = (Inc_int_before_yearly - Inc_int_followup5_yearly)*100000
  Inc_avert10_int_yearly = (Inc_int_before_yearly - Inc_int_followup10_yearly)*100000
  # Clinical incidence averted
  Inc_avert1_int_sixmonthly = (Inc_int_before_sixmonthly - Inc_int_followup1_sixmonthly)*100000
  Inc_avert5_int_sixmonthly = (Inc_int_before_sixmonthly - Inc_int_followup5_sixmonthly)*100000
  Inc_avert10_int_sixmonthly = (Inc_int_before_sixmonthly - Inc_int_followup10_sixmonthly)*100000
  
  #################
  ### Ages 0-5 ###
  
  # summed yearly clinical cases in 0 - 5s
  Inc_05_yearly <- Inc_agegroups_yearly[Inc_agegroups_yearly$age_group %in% age05,] %>% dplyr::group_by(year) %>% dplyr::summarise(inc = sum(sum)) %>% dplyr::ungroup()
  Inc_05_yearly$cpp <- Inc_05_yearly$inc / tail(pop_05_yearly$n, 1)
  Inc_05_sixmonthly <- Inc_agegroups_sixmonthly[Inc_agegroups_sixmonthly$age_group %in% age05,] %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(inc = sum(sum)) %>% dplyr::ungroup()
  Inc_05_sixmonthly$cpp <- Inc_05_sixmonthly$inc / tail(pop_05_sixmonthly$n, 1)
  
  # Clinical incidence before and after
  Inc_05_before_yearly=mean(as.numeric(Inc_05_yearly$cpp[which(Inc_05_yearly$year %in% c(1:5))]))
  Inc_05_followup1_yearly=as.numeric(Inc_05_yearly$cpp[which(Inc_05_yearly$year==6)])
  Inc_05_followup5_yearly=as.numeric(Inc_05_yearly$cpp[which(Inc_05_yearly$year==10)])
  Inc_05_followup10_yearly=as.numeric(Inc_05_yearly$cpp[which(Inc_05_yearly$year==15)])
  # Clinical incidence before and after
  Inc_05_before_sixmonthly=mean(as.numeric(Inc_05_sixmonthly$cpp[which(Inc_05_sixmonthly$sixmonths %in% c(1:5))]))
  Inc_05_followup1_sixmonthly=as.numeric(Inc_05_sixmonthly$cpp[which(Inc_05_sixmonthly$sixmonths==6)])
  Inc_05_followup5_sixmonthly=as.numeric(Inc_05_sixmonthly$cpp[which(Inc_05_sixmonthly$sixmonths==10)])
  Inc_05_followup10_sixmonthly=as.numeric(Inc_05_sixmonthly$cpp[which(Inc_05_sixmonthly$sixmonths==15)])
  
  # Clinical incidence reduction
  Inc_red1_05_yearly = ((Inc_05_before_yearly - Inc_05_followup1_yearly)/Inc_05_before_yearly)*100
  Inc_red5_05_yearly = ((Inc_05_before_yearly - Inc_05_followup5_yearly)/Inc_05_before_yearly)*100
  Inc_red10_05_yearly = ((Inc_05_before_yearly - Inc_05_followup10_yearly)/Inc_05_before_yearly)*100
  # Clinical incidence reduction
  Inc_red1_05_sixmonthly = ((Inc_05_before_sixmonthly - Inc_05_followup1_sixmonthly)/Inc_05_before_sixmonthly)*100
  Inc_red5_05_sixmonthly = ((Inc_05_before_sixmonthly - Inc_05_followup5_sixmonthly)/Inc_05_before_sixmonthly)*100
  Inc_red10_05_sixmonthly = ((Inc_05_before_sixmonthly - Inc_05_followup10_sixmonthly)/Inc_05_before_sixmonthly)*100
  
  # # Clinical incidence averted
  # Inc_avert1_05_yearly = (Inc_05_before_yearly - Inc_05_followup1_yearly)*100000
  # Inc_avert5_05_yearly = (Inc_05_before_yearly - Inc_05_followup5_yearly)*100000
  # Inc_avert10_05_yearly = (Inc_05_before_yearly - Inc_05_followup10_yearly)*100000
  # # Clinical incidence averted
  # Inc_avert1_05_sixmonthly = (Inc_05_before_sixmonthly - Inc_05_followup1_sixmonthly)*100000
  # Inc_avert5_05_sixmonthly = (Inc_05_before_sixmonthly - Inc_05_followup5_sixmonthly)*100000
  # Inc_avert10_05_sixmonthly = (Inc_05_before_sixmonthly - Inc_05_followup10_sixmonthly)*100000
  
  #################
  ### Ages 2-10 ###
  
  # summed yearly clinical cases in 2 - 10s
  Inc_210_yearly <- Inc_agegroups_yearly[Inc_agegroups_yearly$age_group %in% age210,] %>% dplyr::group_by(year) %>% dplyr::summarise(inc = sum(sum)) %>% dplyr::ungroup()
  Inc_210_yearly$cpp <- Inc_210_yearly$inc / tail(pop_210_yearly$n, 1)
  Inc_210_sixmonthly <- Inc_agegroups_sixmonthly[Inc_agegroups_sixmonthly$age_group %in% age210,] %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(inc = sum(sum)) %>% dplyr::ungroup()
  Inc_210_sixmonthly$cpp <- Inc_210_sixmonthly$inc / tail(pop_210_sixmonthly$n, 1)
  
  # Clinical incidence before and after
  Inc_210_before_yearly=mean(as.numeric(Inc_210_yearly$cpp[which(Inc_210_yearly$year %in% c(1:5))]))
  Inc_210_followup1_yearly=as.numeric(Inc_210_yearly$cpp[which(Inc_210_yearly$year==6)])
  Inc_210_followup5_yearly=as.numeric(Inc_210_yearly$cpp[which(Inc_210_yearly$year==10)])
  Inc_210_followup10_yearly=as.numeric(Inc_210_yearly$cpp[which(Inc_210_yearly$year==15)])
  # Clinical incidence before and after
  Inc_210_before_sixmonthly=mean(as.numeric(Inc_210_sixmonthly$cpp[which(Inc_210_sixmonthly$sixmonths %in% c(1:5))]))
  Inc_210_followup1_sixmonthly=as.numeric(Inc_210_sixmonthly$cpp[which(Inc_210_sixmonthly$sixmonths==6)])
  Inc_210_followup5_sixmonthly=as.numeric(Inc_210_sixmonthly$cpp[which(Inc_210_sixmonthly$sixmonths==10)])
  Inc_210_followup10_sixmonthly=as.numeric(Inc_210_sixmonthly$cpp[which(Inc_210_sixmonthly$sixmonths==15)])
  
  # Clinical incidence reduction
  Inc_red1_210_yearly = ((Inc_210_before_yearly - Inc_210_followup1_yearly)/Inc_210_before_yearly)*100
  Inc_red5_210_yearly = ((Inc_210_before_yearly - Inc_210_followup5_yearly)/Inc_210_before_yearly)*100
  Inc_red10_210_yearly = ((Inc_210_before_yearly - Inc_210_followup10_yearly)/Inc_210_before_yearly)*100
  # Clinical incidence reduction
  Inc_red1_210_sixmonthly = ((Inc_210_before_sixmonthly - Inc_210_followup1_sixmonthly)/Inc_210_before_sixmonthly)*100
  Inc_red5_210_sixmonthly = ((Inc_210_before_sixmonthly - Inc_210_followup5_sixmonthly)/Inc_210_before_sixmonthly)*100
  Inc_red10_210_sixmonthly = ((Inc_210_before_sixmonthly - Inc_210_followup10_sixmonthly)/Inc_210_before_sixmonthly)*100
  
  # # Clinical incidence averted
  # Inc_avert1_210_yearly = (Inc_210_before_yearly - Inc_210_followup1_yearly)*100000
  # Inc_avert5_210_yearly = (Inc_210_before_yearly - Inc_210_followup5_yearly)*100000
  # Inc_avert10_210_yearly = (Inc_210_before_yearly - Inc_210_followup10_yearly)*100000
  # # Clinical incidence averted
  # Inc_avert1_210_sixmonthly = (Inc_210_before_sixmonthly - Inc_210_followup1_sixmonthly)*100000
  # Inc_avert5_210_sixmonthly = (Inc_210_before_sixmonthly - Inc_210_followup5_sixmonthly)*100000
  # Inc_avert10_210_sixmonthly = (Inc_210_before_sixmonthly - Inc_210_followup10_sixmonthly)*100000
  
  ###########################
  ### INCDIENCE OVER TIME ###
  
  # incidence over time 
  nInc_cont_yearly = merge(total_pop_age_yearly, nInc, by.x=c("age_group","year"), by.y=c("age_group","year"), all = F)
  nInc_cont_sixmonthly = merge(total_pop_age_sixmonthly, nInc, by.x=c("age_group","sixmonths"), by.y=c("age_group","sixmonths"), all = F)
  
  # by age group
  Inc_agegroups_cont_yearly = nInc_cont_yearly %>% dplyr::mutate(cpp = value/n) %>% dplyr::select(timestep, age_group, year, n, value)
  Inc_agegroups_cont_sixmonthly = nInc_cont_sixmonthly %>% dplyr::mutate(cpp = value/n) %>% dplyr::select(timestep, age_group, sixmonths, n, value)
  
  # All ages
  Inc_allages_cont_yearly = nInc_cont_yearly %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinc = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(cpp = intinc/intn) %>% 
    dplyr::group_by(year) %>% 
    dplyr::mutate(cppcum = cumsum(cpp)) %>% 
    dplyr::ungroup()
  Inc_allages_cont_sixmonthly = nInc_cont_sixmonthly %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinc = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(cpp = intinc/intn) %>% 
    dplyr::group_by(sixmonths) %>% 
    dplyr::mutate(cppcum = cumsum(cpp)) %>% 
    dplyr::ungroup()
  
  # Intervention age group
  Inc_int_cont_yearly = nInc_cont_yearly %>% 
    dplyr::filter(age_group %in% ageint) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinc = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(cpp = intinc/intn) %>% 
    dplyr::group_by(year) %>% 
    dplyr::mutate(cppcum = cumsum(cpp)) %>% 
    dplyr::ungroup()
  Inc_int_cont_sixmonthly = nInc_cont_sixmonthly %>% 
    dplyr::filter(age_group %in% ageint) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinc = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(cpp = intinc/intn) %>% 
    dplyr::group_by(sixmonths) %>% 
    dplyr::mutate(cppcum = cumsum(cpp)) %>% 
    dplyr::ungroup()
  
  # 0-5 age group
  Inc_05_cont_yearly = nInc_cont_yearly %>% 
    dplyr::filter(age_group %in% age05) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinc = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(cpp = intinc/intn) %>% 
    dplyr::group_by(year) %>% 
    dplyr::mutate(cppcum = cumsum(cpp)) %>% 
    dplyr::ungroup()
  Inc_05_cont_sixmonthly = nInc_cont_sixmonthly %>% 
    dplyr::filter(age_group %in% age05) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinc = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(cpp = intinc/intn) %>% 
    dplyr::group_by(sixmonths) %>% 
    dplyr::mutate(cppcum = cumsum(cpp)) %>% 
    dplyr::ungroup()
  
  # 2-10 age group
  Inc_210_cont_yearly = nInc_cont_yearly %>% 
    dplyr::filter(age_group %in% age210) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinc = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(cpp = intinc/intn) %>% 
    dplyr::group_by(year) %>% 
    dplyr::mutate(cppcum = cumsum(cpp)) %>% 
    dplyr::ungroup()
  Inc_210_cont_sixmonthly = nInc_cont_sixmonthly %>% 
    dplyr::filter(age_group %in% age210) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinc = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(cpp = intinc/intn) %>% 
    dplyr::group_by(sixmonths) %>% 
    dplyr::mutate(cppcum = cumsum(cpp)) %>% 
    dplyr::ungroup()
  
  #####################################################
  ##################### OUTPUT ########################
  #####################################################
  
  ##################
  ### Reductions ###
  
  # Final row with outputs to return
  reductions = cbind.data.frame(scenario_params$Scenario_Name,
                                scenario_params$SEED,
                                Inc_red1_all_yearly,
                                Inc_red5_all_yearly,
                                Inc_red10_all_yearly,
                                # Inc_avert1_all_yearly,
                                # Inc_avert5_all_yearly,
                                # Inc_avert10_all_yearly,
                                Inc_red1_int_yearly,
                                Inc_red5_int_yearly,
                                Inc_red10_int_yearly,
                                Inc_avert1_int_yearly,
                                Inc_avert5_int_yearly,
                                Inc_avert10_int_yearly,
                                Inc_red1_05_yearly,
                                Inc_red5_05_yearly,
                                Inc_red10_05_yearly,
                                # Inc_avert1_05_yearly,
                                # Inc_avert5_05_yearly,
                                # Inc_avert10_05_yearly,
                                Inc_red1_210_yearly,
                                Inc_red5_210_yearly,
                                Inc_red10_210_yearly,
                                # Inc_avert1_210_yearly,
                                # Inc_avert5_210_yearly,
                                # Inc_avert10_210_yearly,
                                Inc_red1_all_sixmonthly,
                                Inc_red5_all_sixmonthly,
                                Inc_red10_all_sixmonthly,
                                # Inc_avert1_all_sixmonthly,
                                # Inc_avert5_all_sixmonthly,
                                # Inc_avert10_all_sixmonthly,
                                Inc_red1_int_sixmonthly,
                                Inc_red5_int_sixmonthly,
                                Inc_red10_int_sixmonthly,
                                Inc_avert1_int_sixmonthly,
                                Inc_avert5_int_sixmonthly,
                                Inc_avert10_int_sixmonthly,
                                Inc_red1_05_sixmonthly,
                                Inc_red5_05_sixmonthly,
                                Inc_red10_05_sixmonthly,
                                # Inc_avert1_05_sixmonthly,
                                # Inc_avert5_05_sixmonthly,
                                # Inc_avert10_05_sixmonthly,
                                Inc_red1_210_sixmonthly,
                                Inc_red5_210_sixmonthly,
                                Inc_red10_210_sixmonthly
                                # Inc_avert1_210_sixmonthly,
                                # Inc_avert5_210_sixmonthly,
                                # Inc_avert10_210_sixmonthly
                                )
  colnames(reductions) = c("Scenario_Name", "seed",
                           paste0(measure_name,"Inc_red1_all_yearly"),
                           paste0(measure_name,"Inc_red5_all_yearly"),
                           paste0(measure_name,"Inc_red10_all_yearly"),
                           paste0(measure_name,"Inc_red1_int_yearly"),
                           paste0(measure_name,"Inc_red5_int_yearly"),
                           paste0(measure_name,"Inc_red10_int_yearly"),
                           paste0(measure_name,"Inc_avert1_int_yearly"),
                           paste0(measure_name,"Inc_avert5_int_yearly"),
                           paste0(measure_name,"Inc_avert10_int_yearly"),
                           paste0(measure_name,"Inc_red1_05_yearly"),
                           paste0(measure_name,"Inc_red5_05_yearly"),
                           paste0(measure_name,"Inc_red10_05_yearly"),
                           paste0(measure_name,"Inc_red1_210_yearly"),
                           paste0(measure_name,"Inc_red5_210_yearly"),
                           paste0(measure_name,"Inc_red10_210_yearly"),
                           paste0(measure_name,"Inc_red1_all_sixmonthly"),
                           paste0(measure_name,"Inc_red5_all_sixmonthly"),
                           paste0(measure_name,"Inc_red10_all_sixmonthly"),
                           paste0(measure_name,"Inc_red1_int_sixmonthly"),
                           paste0(measure_name,"Inc_red5_int_sixmonthly"),
                           paste0(measure_name,"Inc_red10_int_sixmonthly"),
                           paste0(measure_name,"Inc_avert1_int_sixmonthly"),
                           paste0(measure_name,"Inc_avert5_int_sixmonthly"),
                           paste0(measure_name,"Inc_avert10_int_sixmonthly"),
                           paste0(measure_name,"Inc_red1_05_sixmonthly"),
                           paste0(measure_name,"Inc_red5_05_sixmonthly"),
                           paste0(measure_name,"Inc_red10_05_sixmonthly"),
                           paste0(measure_name,"Inc_red1_210_sixmonthly"),
                           paste0(measure_name,"Inc_red5_210_sixmonthly"),
                           paste0(measure_name,"Inc_red10_210_sixmonthly"))
  
  ##################
  ### Timeseries ###
  
  timeseries = list(Inc_allages_cont_yearly,
                    Inc_int_cont_yearly,
                    Inc_05_cont_yearly,
                    Inc_210_cont_yearly,
                    Inc_allages_cont_sixmonthly,
                    Inc_int_cont_sixmonthly,
                    Inc_05_cont_sixmonthly,
                    Inc_210_cont_sixmonthly)
  
  names(timeseries) = c(paste0(measure_name,"Inc_allages_cont_yearly"),
                        paste0(measure_name,"Inc_int_cont_yearly"),
                        paste0(measure_name,"Inc_05_cont_yearly"),
                        paste0(measure_name,"Inc_210_cont_yearly"),
                        paste0(measure_name,"Inc_allages_cont_sixmonthly"),
                        paste0(measure_name,"Inc_int_cont_sixmonthly"),
                        paste0(measure_name,"Inc_05_cont_sixmonthly"),
                        paste0(measure_name,"Inc_210_cont_sixmonthly"))
  
  
  return(list(reductions, timeseries))
}

#########################
### OUTCOMES FUNCTION ###

# Function to calculate prevalence in 2-10 at baseline 
calculate_baseline_prev210 = function(om_result, scenario_params) {
  
  ###################
  #### AGE GROUP ####

  # define age groups
  age_groups <- c(0.165,0.25,0.5,seq(1,10,1),15,20,50,100)
  minIntAge=0.25
  age210 = seq(which(age_groups==2)+1,which(age_groups==10))

  ################
  #### LABELS ####

  colnames(om_result) = c("time", "age_group", "measure", "value")

  #########################
  #### REMOVE MEASURES ####

  # # Remove measures without age group
  measures_to_remove = c(7, 21, 26, 31:36)
  om_result = om_result[-which(om_result$measure %in% measures_to_remove),]

  #################################
  #### REMOVE 1ST OBSERVATIONS ####

  # Remove first measurement as it includes all cases until then
  to_remove = which(om_result$time %in% c(1,2))
  om_result = om_result[-to_remove,]

  #############################
  #### CREATE NEW TIMESTEP ####

  # Timestep 0 starts at startDate so need to remove 2 start dates
  om_result$timestep = om_result$time - 1

  # Get date
  # om_result$date = time.to.date(om_result$timestep, 5, "2030-01-01")
  # om_result$date = time.to.date(om_result$timestep, 5, startDate)

  ##############################
  #### DEPLOYMENT TIMESTEPS ####

  # Get deployment timesteps
  massVaccination = om_result[which(om_result$measure == 22),] %>%
    # mutate(timestep = factor(time)) %>%
    group_by(timestep) %>%
    dplyr::summarise(total = sum(value)) %>%
    ungroup() %>%
    filter(total > 0)

  # Deployment timesteps
  DeploymentTimesteps = massVaccination$timestep

  #############################
  #### REFERENCE TIMESTEPS ####

  # List of reference timesteps to evaluate BEFORE: by year
  timesteps_5years_beforeInt_12mo = seq(DeploymentTimesteps[1]-5*73, DeploymentTimesteps[1]-1, 1)

  # List of reference timesteps to evaluate BEFORE: by 6 months
  timesteps_5years_beforeInt_6mo = c(seq(DeploymentTimesteps[1]-5*73, DeploymentTimesteps[1]-1-5*73+36, 1),
                                     seq(DeploymentTimesteps[1]-4*73, DeploymentTimesteps[1]-1-4*73+36, 1),
                                     seq(DeploymentTimesteps[1]-3*73, DeploymentTimesteps[1]-1-3*73+36, 1),
                                     seq(DeploymentTimesteps[1]-2*73, DeploymentTimesteps[1]-1-2*73+36, 1),
                                     seq(DeploymentTimesteps[1]-1*73, DeploymentTimesteps[1]-1-1*73+36, 1))

  # List of reference timesteps to evaluate AFTER: by year
  timesteps_y1_12mo = seq(DeploymentTimesteps[1], DeploymentTimesteps[2]-1, 1)
  timesteps_y5_12mo = seq(DeploymentTimesteps[5], DeploymentTimesteps[6]-1, 1)
  timesteps_y10_12mo = seq(DeploymentTimesteps[10], DeploymentTimesteps[11]-1, 1)

  # List of reference timesteps to evaluate AFTER: by 6 months
  timesteps_y1_6mo = seq(DeploymentTimesteps[1], DeploymentTimesteps[1]-1+36, 1)
  timesteps_y5_6mo = seq(DeploymentTimesteps[5], DeploymentTimesteps[5]-1+36, 1)
  timesteps_y10_6mo = seq(DeploymentTimesteps[10], DeploymentTimesteps[10]-1+36, 1)

  ############################
  #### ADD ADJUSTED YEAR  ####

  # Reference
  reference = seq(DeploymentTimesteps[1]-5*73,max(om_result$timestep),73)

  # Add adjusted year (based on reference timesteps and 1st deployment)
  om_result$year = c(rep(0,nrow(om_result[which(om_result$timestep < reference[1]),])),
                     rep(seq(1,length(reference)-1), each =  length(unique(om_result$age_group))*length(unique(om_result$measure))*73),
                     rep(length(reference),nrow(om_result[which(om_result$timestep >= reference[length(reference)]),])))

  # # Add 6 months (based on reference timesteps and 1st deployment)
  om_result = om_result %>%
    dplyr::mutate(sixmonths = NA,
                  sixmonths = ifelse(timestep < reference[1], 0, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-5*73, DeploymentTimesteps[1]-1-5*73+36, 1), 1, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-4*73, DeploymentTimesteps[1]-1-4*73+36, 1), 2, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-3*73, DeploymentTimesteps[1]-1-3*73+36, 1), 3, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-2*73, DeploymentTimesteps[1]-1-2*73+36, 1), 4, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1]-1*73, DeploymentTimesteps[1]-1-1*73+36, 1), 5, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[1], DeploymentTimesteps[1]-1+36, 1), 6, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[2], DeploymentTimesteps[2]-1+36, 1), 7, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[3], DeploymentTimesteps[3]-1+36, 1), 8, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[4], DeploymentTimesteps[4]-1+36, 1), 9, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[5], DeploymentTimesteps[5]-1+36, 1), 10, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[6], DeploymentTimesteps[6]-1+36, 1), 11, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[7], DeploymentTimesteps[7]-1+36, 1), 12, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[8], DeploymentTimesteps[8]-1+36, 1), 13, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[9], DeploymentTimesteps[9]-1+36, 1), 14, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[10], DeploymentTimesteps[10]-1+36, 1), 15, sixmonths),
                  sixmonths = ifelse(timestep %in% seq(DeploymentTimesteps[11], DeploymentTimesteps[11]-1+36, 1), 16, sixmonths))


  ##########################
  #### POPULATION: YEAR ####

  # population per time step per age group
  total_pop = as.data.frame(om_result[om_result$measure == 0, ])

  # calculate total population per age group per year
  total_pop_age_yearly = data.frame(total_pop %>% dplyr::group_by(age_group, year) %>% dplyr::summarise(n = mean(value))) %>% ungroup()

  # sum up total population per year
  total_pop_all_yearly = total_pop %>% dplyr::group_by(timestep, year) %>% dplyr::summarise(sum = sum(value)) %>% dplyr::group_by(year) %>% dplyr::summarise(n = mean(sum)) %>% ungroup()

  # sum up population intervention age-groups over the years
  pop_210_yearly <- total_pop_age_yearly[total_pop_age_yearly$age_group %in% age210,] %>% dplyr::group_by(year) %>% dplyr::summarise(n = sum(n)) %>% ungroup()

  ###############################
  #### POPULATION: SIX MONTH ####

  # calculate total population per age group per year
  total_pop_age_sixmonthly = total_pop[complete.cases(total_pop),] %>% dplyr::group_by(age_group, sixmonths) %>% dplyr::summarise(n = mean(value)) %>% ungroup()

  # sum up total population per year
  total_pop_all_sixmonthly = total_pop %>% dplyr::group_by(timestep, sixmonths) %>% dplyr::summarise(sum = sum(value)) %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(n = mean(sum)) %>% ungroup()

  # sum up population intervention age-groups over the years
  pop_210_sixmonthly <- total_pop_age_sixmonthly[total_pop_age_sixmonthly$age_group %in% age210,] %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(n = sum(n)) %>% ungroup()

  ################
  #### TOTALS ####
  ################

  # Summarize all measures results by summing up over age groups
  agg_om_result_total_yearly = om_result[,-which(names(om_result)=="age_group")] %>% dplyr::group_by(timestep, year, measure) %>% dplyr::summarise(value = sum(value)) %>% ungroup()

  # Summarize all measures results by summing up over age groups
  agg_om_result_total_sixmonthly = om_result[complete.cases(om_result),-which(names(om_result)=="age_group")] %>% dplyr::group_by(timestep, sixmonths, measure) %>% dplyr::summarise(value = sum(value)) %>% ungroup()
  
  
  ####################################
  #### PREVALENCE: ALL INFECTIONS ####
  ####################################
  
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  
  # number of infected per age-group per time-step
  n_allInf = as.data.frame(om_result[om_result$measure == 1,]) # The number of human hosts with an infection (patent or not) on the reporting timestep
  
  # total number of infected per time-step
  n_allInf_total_yearly = as.data.frame(agg_om_result_total_yearly[agg_om_result_total_yearly$measure == 1,])
  n_allInf_total_sixmonthly = as.data.frame(agg_om_result_total_sixmonthly[agg_om_result_total_sixmonthly$measure == 1,])
  
  # add the population size in respective age-groups
  n_allInf_yearly <- inner_join(total_pop_age_yearly, n_allInf, by=c("age_group","year"))
  n_allInf_total_yearly <- inner_join(total_pop_all_yearly, n_allInf_total_yearly, by=c("year"))
  n_allInf_sixmonthly <- inner_join(total_pop_age_sixmonthly, n_allInf, by=c("age_group","sixmonths"))
  n_allInf_total_sixmonthly <- inner_join(total_pop_all_sixmonthly, n_allInf_total_sixmonthly, by=c("sixmonths"))
  
  # divide number of infected by respective age-group
  prev_allInf_agegroups_cont_yearly = n_allInf_yearly %>% mutate(prev = value/n)
  prev_allInf_agegroups_cont_sixmonthly = n_allInf_sixmonthly %>% mutate(prev = value/n)
  
  # Children 2-10
  prev_allInf_210_cont_yearly = prev_allInf_agegroups_cont_yearly %>% 
    dplyr::filter(age_group %in% age210) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  prev_allInf_210_cont_sixmonthly = prev_allInf_agegroups_cont_sixmonthly %>% 
    dplyr::filter(age_group %in% age210) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  
  # yearly average prevalence
  prev_allInf_210_yearly = prev_allInf_210_cont_yearly %>% dplyr::group_by(year) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  prev_allInf_210_sixmonthly = prev_allInf_210_cont_sixmonthly %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  
  # prevalence baseline
  prev_allInf_210_baseline_yearly=mean(as.numeric(prev_allInf_210_yearly$avg[which(prev_allInf_210_yearly$year %in% c(1:5))]))
  prev_allInf_210_baseline_sixmonthly=mean(as.numeric(prev_allInf_210_sixmonthly$avg[which(prev_allInf_210_sixmonthly$sixmonths %in% c(1:5))]))
  
  
  #######################################
  #### PREVALENCE: PATENT INFECTIONS ####
  #######################################
  
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  
  # number of infected per age-group per time-step
  n_patentInf = as.data.frame(om_result[om_result$measure == 3,]) # The number of human hosts with an infection (patent or not) on the reporting timestep
  
  # total number of infected per time-step
  n_patentInf_total_yearly = as.data.frame(agg_om_result_total_yearly[agg_om_result_total_yearly$measure == 3,])
  n_patentInf_total_sixmonthly = as.data.frame(agg_om_result_total_sixmonthly[agg_om_result_total_sixmonthly$measure == 3,])
  
  # add the population size in respective age-groups
  n_patentInf_yearly <- inner_join(total_pop_age_yearly, n_patentInf, by=c("age_group","year"))
  n_patentInf_total_yearly <- inner_join(total_pop_all_yearly, n_patentInf_total_yearly, by=c("year"))
  n_patentInf_sixmonthly <- inner_join(total_pop_age_sixmonthly, n_patentInf, by=c("age_group","sixmonths"))
  n_patentInf_total_sixmonthly <- inner_join(total_pop_all_sixmonthly, n_patentInf_total_sixmonthly, by=c("sixmonths"))
  
  # divide number of infected by respective age-group
  prev_patentInf_agegroups_cont_yearly = n_patentInf_yearly %>% mutate(prev = value/n)
  prev_patentInf_agegroups_cont_sixmonthly = n_patentInf_sixmonthly %>% mutate(prev = value/n)
  
  # Children 2-10
  prev_patentInf_210_cont_yearly = prev_patentInf_agegroups_cont_yearly %>% 
    dplyr::filter(age_group %in% age210) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), year=mean(year)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  prev_patentInf_210_cont_sixmonthly = prev_patentInf_agegroups_cont_sixmonthly %>% 
    dplyr::filter(age_group %in% age210) %>% 
    dplyr::group_by(timestep) %>%
    dplyr::summarise(intinf = sum(value), intn=sum(n), sixmonths=mean(sixmonths)) %>% 
    dplyr::mutate(prev = intinf/intn) %>% 
    dplyr::ungroup()
  
  # yearly average prevalence
  prev_patentInf_210_yearly = prev_patentInf_210_cont_yearly %>% dplyr::group_by(year) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  prev_patentInf_210_sixmonthly = prev_patentInf_210_cont_sixmonthly %>% dplyr::group_by(sixmonths) %>% dplyr::summarise(avg = mean(prev)) %>% ungroup()
  
  # prevalence before and after
  prev_patentInf_210_baseline_yearly=mean(as.numeric(prev_patentInf_210_yearly$avg[which(prev_patentInf_210_yearly$year %in% c(1:5))]))
  prev_patentInf_210_baseline_sixmonthly=mean(as.numeric(prev_patentInf_210_sixmonthly$avg[which(prev_patentInf_210_sixmonthly$sixmonths %in% c(1:5))]))
  
  
  #####################################################
  ##################### OUTPUT ########################
  #####################################################
  
  
  # Final row with outputs to return
  baseline_prev210 = cbind.data.frame(scenario_params$Scenario_Name,
                                      scenario_params$SEED,
                                      prev_allInf_210_baseline_yearly,
                                      prev_allInf_210_baseline_sixmonthly,
                                      prev_patentInf_210_baseline_yearly,
                                      prev_patentInf_210_baseline_sixmonthly)
  colnames(baseline_prev210) = c("Scenario_Name", "seed",
                                 "prev_allInf_210_baseline_yearly",
                                 "prev_allInf_210_baseline_sixmonthly",
                                 "prev_patentInf_210_baseline_yearly",
                                 "prev_patentInf_210_baseline_sixmonthly")
  
  return(baseline_prev210)
}

# Function which calculates all outcomes
calculate_outputs <- function(om_result, scenario_params, cont){
  
  #########################
  #### POPULATION DATA ####
  
  # Run function
  demographic_results = getDemographics(om_result, scenario_params)
  
  # Extract results
  om_result = demographic_results[[1]]
  agg_om_result_total_yearly = demographic_results[[2]]
  agg_om_result_total_sixmonthly = demographic_results[[3]] 
  total_pop_age_yearly = demographic_results[[4]]
  total_pop_age_sixmonthly = demographic_results[[5]]
  total_pop_all_yearly = demographic_results[[6]] 
  total_pop_all_sixmonthly = demographic_results[[7]]
  pop_int_yearly = demographic_results[[8]] 
  pop_int_sixmonthly = demographic_results[[9]]
  pop_05_yearly = demographic_results[[10]] 
  pop_05_sixmonthly = demographic_results[[11]]
  pop_210_yearly = demographic_results[[12]] 
  pop_210_sixmonthly = demographic_results[[13]]
  ageint = demographic_results[[14]] 
  age05 = demographic_results[[15]] 
  age210 = demographic_results[[16]]
  

  ####################################
  #### PREVALENCE: ALL INFECTIONS ####
  ####################################
  
  prevalence_allInfections = getPrevalenceOutcomes(om_result = om_result, 
                                                   agg_om_result_total_yearly = agg_om_result_total_yearly, 
                                                   agg_om_result_total_sixmonthly = agg_om_result_total_sixmonthly, 
                                                   total_pop_age_yearly = total_pop_age_yearly, 
                                                   total_pop_age_sixmonthly = total_pop_age_sixmonthly,
                                                   total_pop_all_yearly = total_pop_all_yearly, 
                                                   total_pop_all_sixmonthly = total_pop_all_sixmonthly,
                                                   ageint = ageint, age05 = age05, age210 = age210,
                                                   1, "all", scenario_params)
  prevalence_allInfections_reductions = prevalence_allInfections[[1]]
  prevalence_allInfections_timeseries = prevalence_allInfections[[2]]
  
  #######################################
  #### PREVALENCE: PATENT INFECTIONS ####
  #######################################
  
  prevalence_patentInfections = getPrevalenceOutcomes(om_result = om_result, 
                                                      agg_om_result_total_yearly = agg_om_result_total_yearly, 
                                                      agg_om_result_total_sixmonthly = agg_om_result_total_sixmonthly, 
                                                      total_pop_age_yearly = total_pop_age_yearly, 
                                                      total_pop_age_sixmonthly = total_pop_age_sixmonthly,
                                                      total_pop_all_yearly = total_pop_all_yearly, 
                                                      total_pop_all_sixmonthly = total_pop_all_sixmonthly,
                                                      ageint = ageint, age05 = age05, age210 = age210,
                                                      3, "patent", scenario_params)
  prevalence_patentInfections_reductions = prevalence_patentInfections[[1]]
  prevalence_patentInfections_timeseries = prevalence_patentInfections[[2]]
  
  ###################################
  #### INCIDENCE: ALL INFECTIONS ####
  ###################################
  
  # Extract the clinical case numbers (uncomplicated malaria)
  # An episode of uncomplicated malaria is a period during which an individual has symptoms caused by 
  # malaria parasites present at the time of illness, where the symptoms do not qualifying as severe malaria. 
  
  incidence_allInfections = getIncidenceOutcomes(om_result = om_result, 
                                                 total_pop_age_yearly = total_pop_age_yearly, 
                                                 total_pop_age_sixmonthly = total_pop_age_sixmonthly,
                                                 total_pop_all_yearly = total_pop_all_yearly, 
                                                 total_pop_all_sixmonthly = total_pop_all_sixmonthly,
                                                 pop_int_yearly = pop_int_yearly, 
                                                 pop_int_sixmonthly = pop_int_sixmonthly,
                                                 pop_05_yearly = pop_05_yearly, 
                                                 pop_05_sixmonthly = pop_05_sixmonthly,
                                                 pop_210_yearly = pop_210_yearly, 
                                                 pop_210_sixmonthly = pop_210_sixmonthly,
                                                 ageint = ageint, age05 = age05, age210 = age210,
                                                 1, "all", scenario_params)
  incidence_allInfections_reductions = incidence_allInfections[[1]]
  incidence_allInfections_timeseries = incidence_allInfections[[2]]
  
  ########################################
  #### INCIDENCE: CLINICAL INFECTIONS ####
  ########################################
  
  # Extract the clinical case numbers (uncomplicated malaria)
  # An episode of uncomplicated malaria is a period during which an individual has symptoms caused by 
  # malaria parasites present at the time of illness, where the symptoms do not qualifying as severe malaria. 
  
  incidence_clinicalInfections = getIncidenceOutcomes(om_result = om_result, 
                                                      total_pop_age_yearly = total_pop_age_yearly, 
                                                      total_pop_age_sixmonthly = total_pop_age_sixmonthly,
                                                      total_pop_all_yearly = total_pop_all_yearly, 
                                                      total_pop_all_sixmonthly = total_pop_all_sixmonthly,
                                                      pop_int_yearly = pop_int_yearly, 
                                                      pop_int_sixmonthly = pop_int_sixmonthly,
                                                      pop_05_yearly = pop_05_yearly, 
                                                      pop_05_sixmonthly = pop_05_sixmonthly,
                                                      pop_210_yearly = pop_210_yearly, 
                                                      pop_210_sixmonthly = pop_210_sixmonthly,
                                                      ageint = ageint, age05 = age05, age210 = age210,
                                                      14, "clinical", scenario_params)
  incidence_clinicalInfections_reductions = incidence_clinicalInfections[[1]]
  incidence_clinicalInfections_timeseries = incidence_clinicalInfections[[2]]
  
  ####################################
  #### EXPECTED SEVERE INFECTIONS ####
  ####################################
  
  # Extract the severe episodes of malaria
  # Severe malaria is a potentially life-threatening disease, diagnosable by clinical or laboratory evidence of 
  # vital organ dysfunction, requiring in-patient care. An episode of severe malaria is a period during which an 
  # individual has symptoms, qualifying as severe malaria, caused by malaria parasites present at the time of illness. 
  
  # 78	expectedSevere
  
  incidence_expectedSevereInfections = getIncidenceOutcomes(om_result = om_result, 
                                                            total_pop_age_yearly = total_pop_age_yearly, 
                                                            total_pop_age_sixmonthly = total_pop_age_sixmonthly,
                                                            total_pop_all_yearly = total_pop_all_yearly, 
                                                            total_pop_all_sixmonthly = total_pop_all_sixmonthly,
                                                            pop_int_yearly = pop_int_yearly, 
                                                            pop_int_sixmonthly = pop_int_sixmonthly,
                                                            pop_05_yearly = pop_05_yearly, 
                                                            pop_05_sixmonthly = pop_05_sixmonthly,
                                                            pop_210_yearly = pop_210_yearly, 
                                                            pop_210_sixmonthly = pop_210_sixmonthly,
                                                            ageint = ageint, age05 = age05, age210 = age210,
                                                            78, "expectedSevere", scenario_params)
  incidence_expectedSevereInfections_reductions = incidence_expectedSevereInfections[[1]]
  incidence_expectedSevereInfections_timeseries = incidence_expectedSevereInfections[[2]]
  
  ################################
  #### EXPECTED DIRECT DEATHS ####
  ################################
  
  # 74	expectedDirectDeaths	age group	y	Expected number of direct malaria deaths, from those with severe disease.
  
  incidence_expectedDirDeathsInfections = getIncidenceOutcomes(om_result = om_result, 
                                                               total_pop_age_yearly = total_pop_age_yearly, 
                                                               total_pop_age_sixmonthly = total_pop_age_sixmonthly,
                                                               total_pop_all_yearly = total_pop_all_yearly, 
                                                               total_pop_all_sixmonthly = total_pop_all_sixmonthly,
                                                               pop_int_yearly = pop_int_yearly, 
                                                               pop_int_sixmonthly = pop_int_sixmonthly,
                                                               pop_05_yearly = pop_05_yearly, 
                                                               pop_05_sixmonthly = pop_05_sixmonthly,
                                                               pop_210_yearly = pop_210_yearly, 
                                                               pop_210_sixmonthly = pop_210_sixmonthly,
                                                               ageint = ageint, age05 = age05, age210 = age210,
                                                               74, "expectedDirDeaths", scenario_params)
  incidence_expectedDirDeathsInfections_reductions = incidence_expectedDirDeathsInfections[[1]]
  incidence_expectedDirDeathsInfections_timeseries = incidence_expectedDirDeathsInfections[[2]]
  
  ##################################
  #### EXPECTED HOSPITAL DEATHS ####
  ##################################
  
  # 75	expectedHospitalDeaths	age group	y	Subset of measure 74 which occur in hospital.
  
  incidence_expectedHospDeathsInfections = getIncidenceOutcomes(om_result = om_result, 
                                                                total_pop_age_yearly = total_pop_age_yearly, 
                                                                total_pop_age_sixmonthly = total_pop_age_sixmonthly,
                                                                total_pop_all_yearly = total_pop_all_yearly, 
                                                                total_pop_all_sixmonthly = total_pop_all_sixmonthly,
                                                                pop_int_yearly = pop_int_yearly, 
                                                                pop_int_sixmonthly = pop_int_sixmonthly,
                                                                pop_05_yearly = pop_05_yearly, 
                                                                pop_05_sixmonthly = pop_05_sixmonthly,
                                                                pop_210_yearly = pop_210_yearly, 
                                                                pop_210_sixmonthly = pop_210_sixmonthly,
                                                                ageint = ageint, age05 = age05, age210 = age210,
                                                                75, "expectedHospDeaths", scenario_params)
  incidence_expectedHospDeathsInfections_reductions = incidence_expectedHospDeathsInfections[[1]]
  incidence_expectedHospDeathsInfections_timeseries = incidence_expectedHospDeathsInfections[[2]]
  
  ########################
  #### RETURN RESULTS ####
  ########################
  
  # If continuous == T then time-series returned
  # If continuous == F then reductions returned
  
  if(cont){
    
    # Timeseries
    timeseries = list(prevalence_allInfections_timeseries, 
                      prevalence_patentInfections_timeseries,
                      incidence_allInfections_timeseries,
                      incidence_clinicalInfections_timeseries,
                      incidence_expectedSevereInfections_timeseries,
                      incidence_expectedDirDeathsInfections_timeseries,
                      incidence_expectedHospDeathsInfections_timeseries)
    
    return(timeseries)

  }else{
    
    # Final row with outputs to return
    reductions = cbind.data.frame(prevalence_allInfections_reductions,
                                  prevalence_patentInfections_reductions[,3:ncol(prevalence_patentInfections_reductions)],
                                  incidence_allInfections_reductions[,3:ncol(incidence_allInfections_reductions)],
                                  incidence_clinicalInfections_reductions[,3:ncol(incidence_allInfections_reductions)],
                                  incidence_expectedSevereInfections_reductions[,3:ncol(incidence_expectedSevereInfections_reductions)],
                                  incidence_expectedDirDeathsInfections_reductions[,3:ncol(incidence_expectedDirDeathsInfections_reductions)],
                                  incidence_expectedHospDeathsInfections_reductions[,3:ncol(incidence_expectedHospDeathsInfections_reductions)])
    
    return(reductions)
  }
}

#############################################
### OUTCOMES POSTPROCESSING MAIN FUNCTION ###
#############################################

# Wrapper for looping across all simulation results and gathering postprocessing results in a table
postprocess_OM_baseline = function(results_folder, param_table_file, final_table_dest, final_seed_table_dest) {
  
  #######################
  ### PARAMETER TABLE ###
  
  # Load
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  
  ############################
  ### CALCULATION FUNCTION ###
  
  # Create empty object
  processed_OM_sim = NULL
  
  # Run loop
  for( i in 1:nrow(param_table)) {
    
    print(i)
    
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", 
                           param_table[i,]$SEED, "_out.txt", sep="")
    
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      
      # Read in file
      OM_result = read.table(OM_result_file, sep="\t")
      
      # Identify error to skip
      mtry <- try(calculate_baseline_prev210(OM_result, param_table[i,]),
                  silent = F)
      
      # Skip error or calculate outputs
      if (class(mtry) != "try-error") {
        
        # Get outcomes per simulation
        scenario_row = calculate_baseline_prev210(OM_result, param_table[i,])
        
        # Bind results
        processed_OM_sim = data.frame(rbind(processed_OM_sim, scenario_row), stringsAsFactors = FALSE)
        
      } else {
        
        message("Error")
        
      }
    }
  }
  
  #########################
  ### AGGREGATE RESULTS ###
  
  # Summarize results over seeds and create final results tables --> take median of seed results
  
  aggregated_OM = processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prev_allInf_210_baseline_yearly"):length(names(processed_OM_sim))]), median, na.rm=TRUE)
  
  # Tables
  no_seed_table = param_table[,-c(which(colnames(param_table)=="SEED"))]
  no_seed_table = unique(no_seed_table)
  final_seed_table = merge(no_seed_table, processed_OM_sim, by = c("Scenario_Name"))
  final_table = merge(no_seed_table, aggregated_OM, by = c("Scenario_Name"))
  
  ############
  ### SAVE ###
  
  # Write result tables (summarized and with seeds) to files
  write.table(final_table, final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  write.table(final_seed_table, final_seed_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  
  return(final_table)
}

# Wrapper for looping across all simulation results and gathering postprocessing results in a table
postprocess_OM <- function(results_folder, param_table_file, final_table_dest, final_seed_table_dest){
  
  #######################
  ### PARAMETER TABLE ###
  
  # Load
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  
  ############################
  ### CALCULATION FUNCTION ###
  
  # Create empty object
  processed_OM_sim = NULL
  
  # Run loop
  for( i in 1:nrow(param_table)) {
    
    print(i)
    
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", param_table[i,]$SEED, "_out.txt", sep="")
    
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      
      # Read in file
      OM_result = read.table(OM_result_file, sep="\t")
      
      # Identify error to skip
      mtry <- try(calculate_outputs(OM_result, param_table[i,], cont=F), silent = F)
      
      # Skip error or calculate outputs
      if (class(mtry) != "try-error") {
        
        # Get outcomes per simulation
        scenario_row = calculate_outputs(OM_result, param_table[i,], cont=F)
        
        # Bind results
        processed_OM_sim = data.frame(rbind(processed_OM_sim, scenario_row), stringsAsFactors = FALSE)
        
      } else {
        
        message("Error")
        
      }
    }
  }
  
  #########################
  ### AGGREGATE RESULTS ###
  
  # Summarize results over seeds and create final results tables --> take median of seed results
  aggregated_OM = processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prev_red1_allInf_allages_yearly"):length(names(processed_OM_sim))]), median, na.rm=TRUE)
  
  # Tables
  no_seed_table = param_table[,-c(which(colnames(param_table)=="SEED"))]
  no_seed_table = unique(no_seed_table)
  final_seed_table = merge(no_seed_table, processed_OM_sim, by = c("Scenario_Name"))
  final_table = merge(no_seed_table, aggregated_OM, by = c("Scenario_Name"))
  
  ############
  ### SAVE ###
  
  # Write result tables (summarized and with seeds) to files
  write.table(final_table, final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  write.table(final_seed_table, final_seed_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  
}



########################
### UNUSED FUNCTIONS ###
########################

# Wrapper for looping across all simulation results and gathering postprocessing results in a table
# Version to use during adaptive sampling, returns the postprocessing results instead 
# of writing them to a file
postprocess_OM_as = function(results_folder, param_table, follow_up) {
  processed_OM_sim = NULL
  for( i in 1:nrow(param_table)) {
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", 
                           param_table[i,]$SEED, "_out.txt", sep="")
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      OM_result = read.table(OM_result_file, sep="\t")
      # om_result, scenario_params, total_pop, survey_start, survey_end, int_start, int_end, pulsed_int_start
      scenario_row = calculate_outputs(OM_result, param_table[i,], follow_up)
      processed_OM_sim = rbind(processed_OM_sim, scenario_row)
    }
  }
  
  # # Summarize results over seeds and create final results tables
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prev_red1_allInf_allages"):length(names(processed_OM_sim))]),median,na.rm=TRUE)
  no_seed_table = param_table[,-c(which(colnames(param_table)=="SEED"))]
  no_seed_table = unique(no_seed_table)
  final_seed_table = merge(no_seed_table, processed_OM_sim, by = c("Scenario_Name"))
  # final_table = merge(no_seed_table, aggregated_OM, by = c("Scenario_Name"))
  
  return(final_seed_table)
}

