##########################
# Auxiliary functions for postprocessing OpenMalaria results
#
# created 02.10.2019
# monica.golumbeanu@unibas.ch
# modified by lydia.burgert@unibas.ch 
#
#
#
# modified by josephine.malinga@unibas.ch - December 15, 2021
##########################

# library(rapportools)
# library(survival)
# library(cmprsk)
library(zoo)
require(dplyr)

# Function which calculates prevalence reduction given an OpenMalaria
# simulation result.
calculate_outputs = function(om_result, scenario_params, dosetime1, cont) {
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6
  dosetime1 = 385
  dosetime5 = 677
  dosetime10 = 1042
  
  # define age groups 
  age_groups <- c(0, 0.25, 0.3333, 0.4167, 0.5, 0.5833, 1, 1.5, 1.5833, 2, 
                  2.5833, 3, 3.5833, 4, 4.5833, 5, 5.5833, 6, 6.0833, 6.5833, 7, 
                  7.0833, 8, 9, 10, 11, 15, 20, 100)
  
  #create new smaller age groups; 
  minIntAge=c(0.5833, 4.5833, 6.5833, 10.5833)
  maxAge=c(1.5833,5.5833, 7.5833, 11.5833)
  age05 = seq(which(age_groups==0),which(age_groups==5)-1)
  age10 = seq(which(age_groups==0),which(age_groups==10)-1)
  age210 = seq(which(age_groups==2),which(age_groups==10)-1)
  age610 = seq(which(age_groups==6),which(age_groups==10)-1)
  ageint = seq(which(age_groups==minIntAge[1]),which(age_groups==maxAge[2])-1)
 
  # add age and year
  om_result$time<-as.numeric(om_result$time); om_result$age_group<-as.numeric(om_result$age_group) 
  
  # Remove first measurement as it includes all cases until then - and the first five months to make sure its 10y
  to_remove = which(as.numeric(om_result$time) ==1 )
  om_result = om_result[-to_remove,] 
  
  # Add date, year and month

  # get t=files by outcome, merge and calculate outcomes
  clin<-om_result[om_result$measure==14,]
  names(clin)[4]<-"inc"
  
  severe<-om_result[om_result$measure==78,]
  names(severe)[4]<-"cases"
  
  deaths<-om_result[om_result$measure==74,]
  names(deaths)[4]<-"deaths"
  
  prev<-om_result[om_result$measure==1,]
  names(prev)[4]<-"prev"  
  
  total_pop<-om_result[om_result$measure==0,]
  names(total_pop)[4]<-"n"
  
  prev_all = data.frame(inner_join(total_pop, prev, by=c("time", "age_group")))
  nclin =  data.frame(inner_join(total_pop, clin, by=c("time", "age_group")))
  nsevere =  data.frame(inner_join(total_pop, severe, by=c("time", "age_group")))
  ndeaths =  data.frame(inner_join(total_pop, deaths, by=c("time", "age_group")))
  
  rm(clin, prev, severe, deaths, total_pop)
  
  # sum up population intervention age-groups over the years 
  ######################################
  #calculate output clinical cases
  ######################################   
  
  # Extract the clinical case numbers 
  # set up the before and follow up years
  dosetime1 = 385
  dosetime5 = 677
  dosetime10 = 1042
  
  btime = dosetime1 - (year_to_5day*2)
  base_time = dosetime1 - (year_to_5day*5)
  
  ftime_6before = btime + floor((year_to_5day/2))
  ftime_12before = btime + year_to_5day - 1
  ftime_18before = btime + year_to_5day + floor(year_to_5day/2) - 1
  
  ftime_6after = dosetime1 + floor(year_to_5day/2) 
  ftime_12after = dosetime1 + year_to_5day - 1
  ftime_18after = dosetime1 + year_to_5day + floor(year_to_5day/2) - 1
  
  ftime_6after5 = dosetime5 + floor(year_to_5day/2) 
  ftime_12after5 = dosetime5 + year_to_5day - 1
  ftime_18after5 = dosetime5 + year_to_5day + floor(year_to_5day/2) - 1
  
  ftime_6after10 = dosetime10 + floor(year_to_5day/2) 
  ftime_12after10 = dosetime10 + year_to_5day - 1
  ftime_18after10 = dosetime10 + year_to_5day + floor(year_to_5day/2) - 1
  
  # calculate the incidence numbers for all ages and the average reductions
  incidence_int = nclin[nclin$age_group %in% ageint, -c(3,5)]
  incidence_int <- incidence_int %>% group_by(time) %>% 
    summarize(n = sum(n), inc=sum(inc))
  
  cpp_6before = sum(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_6before,]$inc)/mean(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_6before,]$n)
  cpp_12before = sum(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_12before,]$inc)/mean(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_12before,]$n)
  cpp_18before = sum(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_18before,]$inc)/mean(incidence_int[incidence_int$time>=btime & incidence_int$time<=ftime_18before,]$n)
  
  cpp_6after = sum(incidence_int[incidence_int$time>=dosetime1 & incidence_int$time<=ftime_6after,]$inc)/mean(incidence_int[incidence_int$time>=dosetime1 & incidence_int$time<=ftime_6after,]$n)
  cpp_12after = sum(incidence_int[incidence_int$time>=dosetime1 & incidence_int$time<=ftime_12after,]$inc)/mean(incidence_int[incidence_int$time>=dosetime1 & incidence_int$time<=ftime_12after,]$n)
  cpp_18after = sum(incidence_int[incidence_int$time>=dosetime1 & incidence_int$time<=ftime_18after,]$inc)/mean(incidence_int[incidence_int$time>=dosetime1 & incidence_int$time<=ftime_18after,]$n)
  
  inc_red_int_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  inc_red_int_6m = inc_red_int_6m*(inc_red_int_6m >= 0)
  inc_red_int_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  inc_red_int_12m = inc_red_int_12m*(inc_red_int_12m >= 0)
  inc_red_int_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  inc_red_int_18m = inc_red_int_18m*(inc_red_int_18m >= 0)
  
  cpp_6after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_6after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_12after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_18after5,]$inc)/mean(incidence_int[incidence_int$time>=dosetime5 & incidence_int$time<=ftime_18after5,]$n)
  
  inc_red_int_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  inc_red_int_6m5 = inc_red_int_6m5*(inc_red_int_6m5 >= 0)
  inc_red_int_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  inc_red_int_12m5 = inc_red_int_12m5*(inc_red_int_12m5 >= 0)
  inc_red_int_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  inc_red_int_18m5 = inc_red_int_18m5*(inc_red_int_18m5 >= 0)
  
  cpp_6after10 = sum(incidence_int[incidence_int$time>=dosetime10 & incidence_int$time<=ftime_6after10,]$inc)/mean(incidence_int[incidence_int$time>=dosetime10 & incidence_int$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(incidence_int[incidence_int$time>=dosetime10 & incidence_int$time<=ftime_12after10,]$inc)/mean(incidence_int[incidence_int$time>=dosetime10 & incidence_int$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(incidence_int[incidence_int$time>=dosetime10 & incidence_int$time<=ftime_18after10,]$inc)/mean(incidence_int[incidence_int$time>=dosetime10 & incidence_int$time<=ftime_18after10,]$n)
  
  inc_red_int_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  inc_red_int_6m10 = inc_red_int_6m10*(inc_red_int_6m10 >= 0)
  inc_red_int_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  inc_red_int_12m10 = inc_red_int_12m10*(inc_red_int_12m10 >= 0)
  inc_red_int_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  inc_red_int_18m10 = inc_red_int_18m10*(inc_red_int_18m10 >= 0)
  
  cpp_2nd6before = sum(incidence_int[incidence_int$time>=ftime_6before & incidence_int$time<=ftime_12before,]$inc)/mean(incidence_int[incidence_int$time>=ftime_6before & incidence_int$time<=ftime_12before,]$n)
  cpp_3rd6before = sum(incidence_int[incidence_int$time>=ftime_12before & incidence_int$time<=ftime_18before,]$inc)/mean(incidence_int[incidence_int$time>=ftime_12before & incidence_int$time<=ftime_18before,]$n)
  
  cpp_2nd6after10 = sum(incidence_int[incidence_int$time>=ftime_6after10 & incidence_int$time<=ftime_12after10,]$inc)/mean(incidence_int[incidence_int$time>=ftime_6after10 & incidence_int$time<=ftime_12after10,]$n)
  cpp_3rd6after10 = sum(incidence_int[incidence_int$time>=ftime_12after10 & incidence_int$time<=ftime_18after10,]$inc)/mean(incidence_int[incidence_int$time>=ftime_12after10 & incidence_int$time<=ftime_18after10,]$n)

  inc_red_int_2nd6m10 = 100*(cpp_2nd6before - cpp_2nd6after10)/ cpp_2nd6before
  inc_red_int_2nd6m10 = inc_red_int_2nd6m10*(inc_red_int_2nd6m10 >= 0)
 
  inc_red_int_3rd6m10 = 100*(cpp_3rd6before - cpp_3rd6after10)/ cpp_3rd6before
  inc_red_int_3rd6m10 = inc_red_int_3rd6m10*(inc_red_int_3rd6m10 >= 0)

  # calculate the incidence numbers for all ages and the average reductions
  incidence_10 = nclin[nclin$age_group %in% age10, -c(3,5)]
  incidence_10 <- incidence_10 %>% group_by(time) %>% 
    summarize(n = sum(n), inc=sum(inc))
  
  cpp_6before = sum(incidence_10[incidence_10$time>=btime & incidence_10$time<=ftime_6before,]$inc)/mean(incidence_10[incidence_10$time>=btime & incidence_10$time<=ftime_6before,]$n)
  cpp_12before = sum(incidence_10[incidence_10$time>=btime & incidence_10$time<=ftime_12before,]$inc)/mean(incidence_10[incidence_10$time>=btime & incidence_10$time<=ftime_12before,]$n)
  cpp_18before = sum(incidence_10[incidence_10$time>=btime & incidence_10$time<=ftime_18before,]$inc)/mean(incidence_10[incidence_10$time>=btime & incidence_10$time<=ftime_18before,]$n)
  
  cpp_6after = sum(incidence_10[incidence_10$time>=dosetime1 & incidence_10$time<=ftime_6after,]$inc)/mean(incidence_10[incidence_10$time>=dosetime1 & incidence_10$time<=ftime_6after,]$n)
  cpp_12after = sum(incidence_10[incidence_10$time>=dosetime1 & incidence_10$time<=ftime_12after,]$inc)/mean(incidence_10[incidence_10$time>=dosetime1 & incidence_10$time<=ftime_12after,]$n)
  cpp_18after = sum(incidence_10[incidence_10$time>=dosetime1 & incidence_10$time<=ftime_18after,]$inc)/mean(incidence_10[incidence_10$time>=dosetime1 & incidence_10$time<=ftime_18after,]$n)
  
  inc_red_10_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  inc_red_10_6m = inc_red_10_6m*(inc_red_10_6m >= 0)
  inc_red_10_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  inc_red_10_12m = inc_red_10_12m*(inc_red_10_12m >= 0)
  inc_red_10_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  inc_red_10_18m = inc_red_10_18m*(inc_red_10_18m >= 0)
  
  cpp_6after5 = sum(incidence_10[incidence_10$time>=dosetime5 & incidence_10$time<=ftime_6after5,]$inc)/mean(incidence_10[incidence_10$time>=dosetime5 & incidence_10$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(incidence_10[incidence_10$time>=dosetime5 & incidence_10$time<=ftime_12after5,]$inc)/mean(incidence_10[incidence_10$time>=dosetime5 & incidence_10$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(incidence_10[incidence_10$time>=dosetime5 & incidence_10$time<=ftime_18after5,]$inc)/mean(incidence_10[incidence_10$time>=dosetime5 & incidence_10$time<=ftime_18after5,]$n)
  
  inc_red_10_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  inc_red_10_6m5 = inc_red_10_6m5*(inc_red_10_6m5 >= 0)
  inc_red_10_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  inc_red_10_12m5 = inc_red_10_12m5*(inc_red_10_12m5 >= 0)
  inc_red_10_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  inc_red_10_18m5 = inc_red_10_18m5*(inc_red_10_18m5 >= 0)
  
  cpp_6after10 = sum(incidence_10[incidence_10$time>=dosetime10 & incidence_10$time<=ftime_6after10,]$inc)/mean(incidence_10[incidence_10$time>=dosetime10 & incidence_10$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(incidence_10[incidence_10$time>=dosetime10 & incidence_10$time<=ftime_12after10,]$inc)/mean(incidence_10[incidence_10$time>=dosetime10 & incidence_10$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(incidence_10[incidence_10$time>=dosetime10 & incidence_10$time<=ftime_18after10,]$inc)/mean(incidence_10[incidence_10$time>=dosetime10 & incidence_10$time<=ftime_18after10,]$n)
  
  inc_red_10_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  inc_red_10_6m10 = inc_red_10_6m10*(inc_red_10_6m10 >= 0)
  inc_red_10_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  inc_red_10_12m10 = inc_red_10_12m10*(inc_red_10_12m10 >= 0)
  inc_red_10_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  inc_red_10_18m10 = inc_red_10_18m10*(inc_red_10_18m10 >= 0)
  
  # calculate the incidence numbers for age group under 5 and average reductions
  incidence_u5 = nclin[nclin$age_group %in% age05, -c(3,5)]
  incidence_u5 <- incidence_u5 %>% group_by(time) %>% 
    summarise(n = sum(n), inc=sum(inc))
  
  cpp_6before = sum(incidence_u5[incidence_u5$time>=btime & incidence_u5$time<=ftime_6before,]$inc)/mean(incidence_u5[incidence_u5$time>=btime & incidence_u5$time<=ftime_6before,]$n)
  cpp_12before = sum(incidence_u5[incidence_u5$time>=btime & incidence_u5$time<=ftime_12before,]$inc)/mean(incidence_u5[incidence_u5$time>=btime & incidence_u5$time<=ftime_12before,]$n)
  cpp_18before = sum(incidence_u5[incidence_u5$time>=btime & incidence_u5$time<=ftime_18before,]$inc)/mean(incidence_u5[incidence_u5$time>=btime & incidence_u5$time<=ftime_18before,]$n)
  
  cpp_6after = sum(incidence_u5[incidence_u5$time>=dosetime1 & incidence_u5$time<=ftime_6after,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime1 & incidence_u5$time<=ftime_6after,]$n)
  cpp_12after = sum(incidence_u5[incidence_u5$time>=dosetime1 & incidence_u5$time<=ftime_12after,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime1 & incidence_u5$time<=ftime_12after,]$n)
  cpp_18after = sum(incidence_u5[incidence_u5$time>=dosetime1 & incidence_u5$time<=ftime_18after,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime1 & incidence_u5$time<=ftime_18after,]$n)
  
  inc_red_u5_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  inc_red_u5_6m = inc_red_u5_6m*(inc_red_u5_6m >= 0)
  inc_red_u5_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  inc_red_u5_12m = inc_red_u5_12m*(inc_red_u5_12m >= 0)
  inc_red_u5_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  inc_red_u5_18m = inc_red_u5_18m*(inc_red_u5_18m >= 0)
  
  cpp_6after5 = sum(incidence_u5[incidence_u5$time>=dosetime5 & incidence_u5$time<=ftime_6after5,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime5 & incidence_u5$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(incidence_u5[incidence_u5$time>=dosetime5 & incidence_u5$time<=ftime_12after5,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime5 & incidence_u5$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(incidence_u5[incidence_u5$time>=dosetime5 & incidence_u5$time<=ftime_18after5,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime5 & incidence_u5$time<=ftime_18after5,]$n)
  
  inc_red_u5_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  inc_red_u5_6m5 = inc_red_u5_6m5*(inc_red_u5_6m5 >= 0)
  inc_red_u5_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  inc_red_u5_12m5 = inc_red_u5_12m5*(inc_red_u5_12m5 >= 0)
  inc_red_u5_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  inc_red_u5_18m5 = inc_red_u5_18m5*(inc_red_u5_18m5 >= 0)
  
  cpp_6after10 = sum(incidence_u5[incidence_u5$time>=dosetime10 & incidence_u5$time<=ftime_6after10,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime10 & incidence_u5$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(incidence_u5[incidence_u5$time>=dosetime10 & incidence_u5$time<=ftime_12after10,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime10 & incidence_u5$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(incidence_u5[incidence_u5$time>=dosetime10 & incidence_u5$time<=ftime_18after10,]$inc)/mean(incidence_u5[incidence_u5$time>=dosetime10 & incidence_u5$time<=ftime_18after10,]$n)
  
  inc_red_u5_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  inc_red_u5_6m10 = inc_red_u5_6m10*(inc_red_u5_6m10 >= 0)
  inc_red_u5_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  inc_red_u5_12m10 = inc_red_u5_12m10*(inc_red_u5_12m10 >= 0)
  inc_red_u5_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  inc_red_u5_18m10 = inc_red_u5_18m10*(inc_red_u5_18m10 >= 0)
  
  # calculate the incidence numbers for all ages and the average reductions
  incidence_610 = nclin[nclin$age_group %in% age610, -c(3,5)]
  incidence_610 <- incidence_610 %>% group_by(time) %>% 
    summarise(n = sum(n), inc=sum(inc))
  
  cpp_6before = sum(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_6before,]$inc)/mean(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_6before,]$n)
  cpp_12before = sum(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_12before,]$inc)/mean(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_12before,]$n)
  cpp_18before = sum(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_18before,]$inc)/mean(incidence_610[incidence_610$time>=btime & incidence_610$time<=ftime_18before,]$n)
  
  cpp_6after = sum(incidence_610[incidence_610$time>=dosetime1 & incidence_610$time<=ftime_6after,]$inc)/mean(incidence_610[incidence_610$time>=dosetime1 & incidence_610$time<=ftime_6after,]$n)
  cpp_12after = sum(incidence_610[incidence_610$time>=dosetime1 & incidence_610$time<=ftime_12after,]$inc)/mean(incidence_610[incidence_610$time>=dosetime1 & incidence_610$time<=ftime_12after,]$n)
  cpp_18after = sum(incidence_610[incidence_610$time>=dosetime1 & incidence_610$time<=ftime_18after,]$inc)/mean(incidence_610[incidence_610$time>=dosetime1 & incidence_610$time<=ftime_18after,]$n)
  
  inc_red_610_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  inc_red_610_6m = inc_red_610_6m*(inc_red_610_6m >= 0)
  inc_red_610_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  inc_red_610_12m = inc_red_610_12m*(inc_red_610_12m >= 0)
  inc_red_610_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  inc_red_610_18m = inc_red_610_18m*(inc_red_610_18m >= 0)
  
  cpp_6after5 = sum(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_6after5,]$inc)/mean(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_12after5,]$inc)/mean(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_18after5,]$inc)/mean(incidence_610[incidence_610$time>=dosetime5 & incidence_610$time<=ftime_18after5,]$n)
  
  inc_red_610_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  inc_red_610_6m5 = inc_red_610_6m5*(inc_red_610_6m5 >= 0)
  inc_red_610_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  inc_red_610_12m5 = inc_red_610_12m5*(inc_red_610_12m5 >= 0)
  inc_red_610_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  inc_red_610_18m5 = inc_red_610_18m5*(inc_red_610_18m5 >= 0)
  
  cpp_6after10 = sum(incidence_610[incidence_610$time>=dosetime10 & incidence_610$time<=ftime_6after10,]$inc)/mean(incidence_610[incidence_610$time>=dosetime10 & incidence_610$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(incidence_610[incidence_610$time>=dosetime10 & incidence_610$time<=ftime_12after10,]$inc)/mean(incidence_610[incidence_610$time>=dosetime10 & incidence_610$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(incidence_610[incidence_610$time>=dosetime10 & incidence_610$time<=ftime_18after10,]$inc)/mean(incidence_610[incidence_610$time>=dosetime10 & incidence_610$time<=ftime_18after10,]$n)
  
  inc_red_610_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  inc_red_610_6m10 = inc_red_610_6m10*(inc_red_610_6m10 >= 0)
  inc_red_610_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  inc_red_610_12m10 = inc_red_610_12m10*(inc_red_610_12m10 >= 0)
  inc_red_610_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  inc_red_610_18m10 = inc_red_610_18m10*(inc_red_610_18m10 >= 0)
  
  ######################################
  #calculate output prevalence reduction
  ######################################   
  
  # Calculate the prevalence for all the monitored years 
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  # number of infected per age-group per time-step
  prevalence_10 = prev_all[prev_all$age_group %in% age10, -c(3,5)]
  prevalence_10 <- prevalence_10 %>% group_by(time) %>% 
    summarise(n = sum(n), prev=mean(prev))
  
  cpp_6before = mean(prevalence_10[prevalence_10$time>=btime & prevalence_10$time<=ftime_6before,]$prev/prevalence_10[prevalence_10$time>=btime & prevalence_10$time<=ftime_6before,]$n)
  cpp_12before = mean(prevalence_10[prevalence_10$time>=btime & prevalence_10$time<=ftime_12before,]$prev/prevalence_10[prevalence_10$time>=btime & prevalence_10$time<=ftime_12before,]$n)
  cpp_18before = mean(prevalence_10[prevalence_10$time>=btime & prevalence_10$time<=ftime_18before,]$prev/prevalence_10[prevalence_10$time>=btime & prevalence_10$time<=ftime_18before,]$n)
  
  cpp_6after = mean(prevalence_10[prevalence_10$time>=dosetime1 & prevalence_10$time<=ftime_6after,]$prev/prevalence_10[prevalence_10$time>=dosetime1 & prevalence_10$time<=ftime_6after,]$n)
  cpp_12after = mean(prevalence_10[prevalence_10$time>=dosetime1 & prevalence_10$time<=ftime_12after,]$prev/prevalence_10[prevalence_10$time>=dosetime1 & prevalence_10$time<=ftime_12after,]$n)
  cpp_18after = mean(prevalence_10[prevalence_10$time>=dosetime1 & prevalence_10$time<=ftime_18after,]$prev/prevalence_10[prevalence_10$time>=dosetime1 & prevalence_10$time<=ftime_18after,]$n)
  
  prev_red_10_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  prev_red_10_6m = prev_red_10_6m*(prev_red_10_6m >= 0)
  prev_red_10_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  prev_red_10_12m = prev_red_10_12m*(prev_red_10_12m >= 0)
  prev_red_10_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  prev_red_10_18m = prev_red_10_18m*(prev_red_10_18m >= 0)
  
  cpp_6after5 = mean(prevalence_10[prevalence_10$time>=dosetime5 & prevalence_10$time<=ftime_6after5,]$prev/prevalence_10[prevalence_10$time>=dosetime5 & prevalence_10$time<=ftime_6after5,]$n)
  cpp_12after5 = mean(prevalence_10[prevalence_10$time>=dosetime5 & prevalence_10$time<=ftime_12after5,]$prev/prevalence_10[prevalence_10$time>=dosetime5 & prevalence_10$time<=ftime_12after5,]$n)
  cpp_18after5 = mean(prevalence_10[prevalence_10$time>=dosetime5 & prevalence_10$time<=ftime_18after5,]$prev/prevalence_10[prevalence_10$time>=dosetime5 & prevalence_10$time<=ftime_18after5,]$n)
  
  prev_red_10_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  prev_red_10_6m5 = prev_red_10_6m5*(prev_red_10_6m5 >= 0)
  prev_red_10_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  prev_red_10_12m5 = prev_red_10_12m5*(prev_red_10_12m5 >= 0)
  prev_red_10_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  prev_red_10_18m5 = prev_red_10_18m5*(prev_red_10_18m5 >= 0)
  
  cpp_6after10 = mean(prevalence_10[prevalence_10$time>=dosetime10 & prevalence_10$time<=ftime_6after10,]$prev/prevalence_10[prevalence_10$time>=dosetime10 & prevalence_10$time<=ftime_6after10,]$n)
  cpp_12after10 = mean(prevalence_10[prevalence_10$time>=dosetime10 & prevalence_10$time<=ftime_12after10,]$prev/prevalence_10[prevalence_10$time>=dosetime10 & prevalence_10$time<=ftime_12after10,]$n)
  cpp_18after10 = mean(prevalence_10[prevalence_10$time>=dosetime10 & prevalence_10$time<=ftime_18after10,]$prev/prevalence_10[prevalence_10$time>=dosetime10 & prevalence_10$time<=ftime_18after10,]$n)
  
  prev_red_10_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  prev_red_10_6m10 = prev_red_10_6m10*(prev_red_10_6m10 >= 0)
  prev_red_10_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  prev_red_10_12m10 = prev_red_10_12m10*(prev_red_10_12m10 >= 0)
  prev_red_10_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  prev_red_10_18m10 = prev_red_10_18m10*(prev_red_10_18m10 >= 0)
  
  # calculate the incidence numbers for all ages and the average reductions
  prevalence_u5 = prev_all[prev_all$age_group %in% age05, -c(3,5)]
  prevalence_u5 <- prevalence_u5 %>% group_by(time) %>% 
    summarise(n = sum(n), prev=mean(prev))
  
  cpp_6before = mean(prevalence_u5[prevalence_u5$time>=btime & prevalence_u5$time<=ftime_6before,]$prev/prevalence_u5[prevalence_u5$time>=btime & prevalence_u5$time<=ftime_6before,]$n)
  cpp_12before = mean(prevalence_u5[prevalence_u5$time>=btime & prevalence_u5$time<=ftime_12before,]$prev/prevalence_u5[prevalence_u5$time>=btime & prevalence_u5$time<=ftime_12before,]$n)
  cpp_18before = mean(prevalence_u5[prevalence_u5$time>=btime & prevalence_u5$time<=ftime_18before,]$prev/prevalence_u5[prevalence_u5$time>=btime & prevalence_u5$time<=ftime_18before,]$n)
  
  cpp_6after = mean(prevalence_u5[prevalence_u5$time>=dosetime1 & prevalence_u5$time<=ftime_6after,]$prev/prevalence_u5[prevalence_u5$time>=dosetime1 & prevalence_u5$time<=ftime_6after,]$n)
  cpp_12after = mean(prevalence_u5[prevalence_u5$time>=dosetime1 & prevalence_u5$time<=ftime_12after,]$prev/prevalence_u5[prevalence_u5$time>=dosetime1 & prevalence_u5$time<=ftime_12after,]$n)
  cpp_18after = mean(prevalence_u5[prevalence_u5$time>=dosetime1 & prevalence_u5$time<=ftime_18after,]$prev/prevalence_u5[prevalence_u5$time>=dosetime1 & prevalence_u5$time<=ftime_18after,]$n)
  
  prev_red_u5_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  prev_red_u5_6m = prev_red_u5_6m*(prev_red_u5_6m >= 0)
  prev_red_u5_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  prev_red_u5_12m = prev_red_u5_12m*(prev_red_u5_12m >= 0)
  prev_red_u5_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  prev_red_u5_18m = prev_red_u5_18m*(prev_red_u5_18m >= 0)
  
  cpp_6after5 = mean(prevalence_u5[prevalence_u5$time>=dosetime5 & prevalence_u5$time<=ftime_6after5,]$prev/prevalence_u5[prevalence_u5$time>=dosetime5 & prevalence_u5$time<=ftime_6after5,]$n)
  cpp_12after5 = mean(prevalence_u5[prevalence_u5$time>=dosetime5 & prevalence_u5$time<=ftime_12after5,]$prev/prevalence_u5[prevalence_u5$time>=dosetime5 & prevalence_u5$time<=ftime_12after5,]$n)
  cpp_18after5 = mean(prevalence_u5[prevalence_u5$time>=dosetime5 & prevalence_u5$time<=ftime_18after5,]$prev/prevalence_u5[prevalence_u5$time>=dosetime5 & prevalence_u5$time<=ftime_18after5,]$n)
  
  prev_red_u5_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  prev_red_u5_6m5 = prev_red_u5_6m5*(prev_red_u5_6m5 >= 0)
  prev_red_u5_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  prev_red_u5_12m5 = prev_red_u5_12m5*(prev_red_u5_12m5 >= 0)
  prev_red_u5_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  prev_red_u5_18m5 = prev_red_u5_18m5*(prev_red_u5_18m5 >= 0)
  
  cpp_6after10 = mean(prevalence_u5[prevalence_u5$time>=dosetime10 & prevalence_u5$time<=ftime_6after10,]$prev/prevalence_u5[prevalence_u5$time>=dosetime10 & prevalence_u5$time<=ftime_6after10,]$n)
  cpp_12after10 = mean(prevalence_u5[prevalence_u5$time>=dosetime10 & prevalence_u5$time<=ftime_12after10,]$prev/prevalence_u5[prevalence_u5$time>=dosetime10 & prevalence_u5$time<=ftime_12after10,]$n)
  cpp_18after10 = mean(prevalence_u5[prevalence_u5$time>=dosetime10 & prevalence_u5$time<=ftime_18after10,]$prev/prevalence_u5[prevalence_u5$time>=dosetime10 & prevalence_u5$time<=ftime_18after10,]$n)
  
  prev_red_u5_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  prev_red_u5_6m10 = prev_red_u5_6m10*(prev_red_u5_6m10 >= 0)
  prev_red_u5_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  prev_red_u5_12m10 = prev_red_u5_12m10*(prev_red_u5_12m10 >= 0)
  prev_red_u5_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  prev_red_u5_18m10 = prev_red_u5_18m10*(prev_red_u5_18m10 >= 0)
  
  # calculate the prevalence numbers for all ages and the average reductions
  prevalence_210 = prev_all[prev_all$age_group %in% age210, -c(3,5)]
  prevalence_210  <- prevalence_210  %>% group_by(time) %>% 
    summarise(n = sum(n), prev=sum(prev))
  
  # cases per person
  cpp_6before = mean(prevalence_210[prevalence_210$time>=btime & prevalence_210$time<=ftime_6before,]$prev/prevalence_210[prevalence_210$time>=btime & prevalence_210$time<=ftime_6before,]$n)
  cpp_12before = mean(prevalence_210[prevalence_210$time>=btime & prevalence_210$time<=ftime_12before,]$prev/prevalence_210[prevalence_210$time>=btime & prevalence_210$time<=ftime_12before,]$n)
  cpp_18before = mean(prevalence_210[prevalence_210$time>=btime & prevalence_210$time<=ftime_18before,]$prev/prevalence_210[prevalence_210$time>=btime & prevalence_210$time<=ftime_18before,]$n)
  
  cpp_6after = mean(prevalence_210[prevalence_210$time>=dosetime1 & prevalence_210$time<=ftime_6after,]$prev/prevalence_210[prevalence_210$time>=dosetime1 & prevalence_210$time<=ftime_6after,]$n)
  cpp_12after = mean(prevalence_210[prevalence_210$time>=dosetime1 & prevalence_210$time<=ftime_12after,]$prev/prevalence_210[prevalence_210$time>=dosetime1 & prevalence_210$time<=ftime_12after,]$n)
  cpp_18after = mean(prevalence_210[prevalence_210$time>=dosetime1 & prevalence_210$time<=ftime_18after,]$prev/prevalence_210[prevalence_210$time>=dosetime1 & prevalence_210$time<=ftime_18after,]$n)
  
  prev_red_210_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  prev_red_210_6m = prev_red_210_6m*(prev_red_210_6m >= 0)
  prev_red_210_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  prev_red_210_12m = prev_red_210_12m*(prev_red_210_12m >= 0)
  prev_red_210_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  prev_red_210_18m = prev_red_210_18m*(prev_red_210_18m >= 0)
  
  cpp_6after5 = mean(prevalence_210[prevalence_210$time>=dosetime5 & prevalence_210$time<=ftime_6after5,]$prev/prevalence_210[prevalence_210$time>=dosetime5 & prevalence_210$time<=ftime_6after5,]$n)
  cpp_12after5 = mean(prevalence_210[prevalence_210$time>=dosetime5 & prevalence_210$time<=ftime_12after5,]$prev/prevalence_210[prevalence_210$time>=dosetime5 & prevalence_210$time<=ftime_12after5,]$n)
  cpp_18after5 = mean(prevalence_210[prevalence_210$time>=dosetime5 & prevalence_210$time<=ftime_18after5,]$prev/prevalence_210[prevalence_210$time>=dosetime5 & prevalence_210$time<=ftime_18after5,]$n)
  
  prev_red_210_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  prev_red_210_6m5 = prev_red_210_6m5*(prev_red_210_6m5 >= 0)
  prev_red_210_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  prev_red_210_12m5 = prev_red_210_12m5*(prev_red_210_12m5 >= 0)
  prev_red_210_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  prev_red_210_18m5 = prev_red_210_18m5*(prev_red_210_18m5 >= 0)
  
  cpp_6after10 = mean(prevalence_210[prevalence_210$time>=dosetime10 & prevalence_210$time<=ftime_6after10,]$prev/prevalence_210[prevalence_210$time>=dosetime10 & prevalence_210$time<=ftime_6after10,]$n)
  cpp_12after10 = mean(prevalence_210[prevalence_210$time>=dosetime10 & prevalence_210$time<=ftime_12after10,]$prev/prevalence_210[prevalence_210$time>=dosetime10 & prevalence_210$time<=ftime_12after10,]$n)
  cpp_18after10 = mean(prevalence_210[prevalence_210$time>=dosetime10 & prevalence_210$time<=ftime_18after10,]$prev/prevalence_210[prevalence_210$time>=dosetime10 & prevalence_210$time<=ftime_18after10,]$n)
  
  prev_red_210_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  prev_red_210_6m10 = prev_red_210_6m10*(prev_red_210_6m10 >= 0)
  prev_red_210_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  prev_red_210_12m10 = prev_red_210_12m10*(prev_red_210_12m10 >= 0)
  prev_red_210_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  prev_red_210_18m10 = prev_red_210_18m10*(prev_red_210_18m10 >= 0)
  
  # prevalence at baseline
  prevalence_210_before = mean(prevalence_210[prevalence_210$time>=btime & prevalence_210$time<dosetime1,]$prev/prevalence_210[prevalence_210$time>=btime & prevalence_210$time<dosetime1,]$n)
  
  
  #####################################
  # calculate severe outcomes
  #####################################
  
  # calculate the severe numbers for all ages and the average reductions
  severe_int = nsevere[nsevere$age_group %in% ageint, -c(3,5)]
  severe_int <- severe_int %>% group_by(time) %>% 
    summarise(n = sum(n), cases=sum(cases))
  
  cpp_6before = sum(severe_int[severe_int$time>=btime & severe_int$time<=ftime_6before,]$cases)/mean(severe_int[severe_int$time>=btime & severe_int$time<=ftime_6before,]$n)
  cpp_12before = sum(severe_int[severe_int$time>=btime & severe_int$time<=ftime_12before,]$cases)/mean(severe_int[severe_int$time>=btime & severe_int$time<=ftime_12before,]$n)
  cpp_18before = sum(severe_int[severe_int$time>=btime & severe_int$time<=ftime_18before,]$cases)/mean(severe_int[severe_int$time>=btime & severe_int$time<=ftime_18before,]$n)
  
  cpp_6after = sum(severe_int[severe_int$time>=dosetime1 & severe_int$time<=ftime_6after,]$cases)/mean(severe_int[severe_int$time>=dosetime1 & severe_int$time<=ftime_6after,]$n)
  cpp_12after = sum(severe_int[severe_int$time>=dosetime1 & severe_int$time<=ftime_12after,]$cases)/mean(severe_int[severe_int$time>=dosetime1 & severe_int$time<=ftime_12after,]$n)
  cpp_18after = sum(severe_int[severe_int$time>=dosetime1 & severe_int$time<=ftime_18after,]$cases)/mean(severe_int[severe_int$time>=dosetime1 & severe_int$time<=ftime_18after,]$n)
  
  sev_red_int_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  sev_red_int_6m = sev_red_int_6m*(sev_red_int_6m >= 0)
  sev_red_int_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  sev_red_int_12m = sev_red_int_12m*(sev_red_int_12m >= 0)
  sev_red_int_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  sev_red_int_18m = sev_red_int_18m*(sev_red_int_18m >= 0)
  
  cpp_6after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_6after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_12after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_18after5,]$cases)/mean(severe_int[severe_int$time>=dosetime5 & severe_int$time<=ftime_18after5,]$n)
  
  sev_red_int_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  sev_red_int_6m5 = sev_red_int_6m5*(sev_red_int_6m5 >= 0)
  sev_red_int_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  sev_red_int_12m5 = sev_red_int_12m5*(sev_red_int_12m5 >= 0)
  sev_red_int_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  sev_red_int_18m5 = sev_red_int_18m5*(sev_red_int_18m5 >= 0)
  
  cpp_6after10 = sum(severe_int[severe_int$time>=dosetime10 & severe_int$time<=ftime_6after10,]$cases)/mean(severe_int[severe_int$time>=dosetime10 & severe_int$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(severe_int[severe_int$time>=dosetime10 & severe_int$time<=ftime_12after10,]$cases)/mean(severe_int[severe_int$time>=dosetime10 & severe_int$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(severe_int[severe_int$time>=dosetime10 & severe_int$time<=ftime_18after10,]$cases)/mean(severe_int[severe_int$time>=dosetime10 & severe_int$time<=ftime_18after10,]$n)
  
  sev_red_int_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  sev_red_int_6m10 = sev_red_int_6m10*(sev_red_int_6m10 >= 0)
  sev_red_int_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  sev_red_int_12m10 = sev_red_int_12m10*(sev_red_int_12m10 >= 0)
  sev_red_int_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  sev_red_int_18m10 = sev_red_int_18m10*(sev_red_int_18m10 >= 0)
  
  cpp_2nd6before = sum(severe_int[severe_int$time>=ftime_6before & severe_int$time<=ftime_12before,]$cases)/mean(severe_int[severe_int$time>=ftime_6before & severe_int$time<=ftime_12before,]$n)
  cpp_3rd6before = sum(severe_int[severe_int$time>=ftime_12before & severe_int$time<=ftime_18before,]$cases)/mean(severe_int[severe_int$time>=ftime_12before & severe_int$time<=ftime_18before,]$n)
  
  cpp_2nd6after10 = sum(severe_int[severe_int$time>=ftime_6after10 & severe_int$time<=ftime_12after10,]$cases)/mean(severe_int[severe_int$time>=ftime_6after10 & severe_int$time<=ftime_12after10,]$n)
  cpp_3rd6after10 = sum(severe_int[severe_int$time>=ftime_12after10 & severe_int$time<=ftime_18after10,]$cases)/mean(severe_int[severe_int$time>=ftime_12after10 & severe_int$time<=ftime_18after10,]$n)
  
  sev_red_int_2nd6m10 = 100*(cpp_2nd6before - cpp_2nd6after10)/ cpp_2nd6before
  sev_red_int_2nd6m10 = sev_red_int_2nd6m10*(sev_red_int_2nd6m10 >= 0)
  
  sev_red_int_3rd6m10 = 100*(cpp_3rd6before - cpp_3rd6after10)/ cpp_3rd6before
  sev_red_int_3rd6m10 = sev_red_int_3rd6m10*(sev_red_int_3rd6m10 >= 0)
  
  # calculate the severe numbers for all ages and the average reductions
  severe_10 = nsevere[nsevere$age_group %in% age10, -c(3,5)]
  severe_10 <- severe_10 %>% group_by(time) %>% 
    summarise(n = sum(n), cases=sum(cases))
  
  cpp_6before = sum(severe_10[severe_10$time>=btime & severe_10$time<=ftime_6before,]$cases)/mean(severe_10[severe_10$time>=btime & severe_10$time<=ftime_6before,]$n)
  cpp_12before = sum(severe_10[severe_10$time>=btime & severe_10$time<=ftime_12before,]$cases)/mean(severe_10[severe_10$time>=btime & severe_10$time<=ftime_12before,]$n)
  cpp_18before = sum(severe_10[severe_10$time>=btime & severe_10$time<=ftime_18before,]$cases)/mean(severe_10[severe_10$time>=btime & severe_10$time<=ftime_18before,]$n)
  
  cpp_6after = sum(severe_10[severe_10$time>=dosetime1 & severe_10$time<=ftime_6after,]$cases)/mean(severe_10[severe_10$time>=dosetime1 & severe_10$time<=ftime_6after,]$n)
  cpp_12after = sum(severe_10[severe_10$time>=dosetime1 & severe_10$time<=ftime_12after,]$cases)/mean(severe_10[severe_10$time>=dosetime1 & severe_10$time<=ftime_12after,]$n)
  cpp_18after = sum(severe_10[severe_10$time>=dosetime1 & severe_10$time<=ftime_18after,]$cases)/mean(severe_10[severe_10$time>=dosetime1 & severe_10$time<=ftime_18after,]$n)
  
  sev_red_10_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  sev_red_10_6m = sev_red_10_6m*(sev_red_10_6m >= 0)
  sev_red_10_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  sev_red_10_12m = sev_red_10_12m*(sev_red_10_12m >= 0)
  sev_red_10_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  sev_red_10_18m = sev_red_10_18m*(sev_red_10_18m >= 0)
  
  cpp_6after5 = sum(severe_10[severe_10$time>=dosetime5 & severe_10$time<=ftime_6after5,]$cases)/mean(severe_10[severe_10$time>=dosetime5 & severe_10$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(severe_10[severe_10$time>=dosetime5 & severe_10$time<=ftime_12after5,]$cases)/mean(severe_10[severe_10$time>=dosetime5 & severe_10$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(severe_10[severe_10$time>=dosetime5 & severe_10$time<=ftime_18after5,]$cases)/mean(severe_10[severe_10$time>=dosetime5 & severe_10$time<=ftime_18after5,]$n)
  
  sev_red_10_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  sev_red_10_6m5 = sev_red_10_6m5*(sev_red_10_6m5 >= 0)
  sev_red_10_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  sev_red_10_12m5 = sev_red_10_12m5*(sev_red_10_12m5 >= 0)
  sev_red_10_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  sev_red_10_18m5 = sev_red_10_18m5*(sev_red_10_18m5 >= 0)
  
  cpp_6after10 = sum(severe_10[severe_10$time>=dosetime10 & severe_10$time<=ftime_6after10,]$cases)/mean(severe_10[severe_10$time>=dosetime10 & severe_10$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(severe_10[severe_10$time>=dosetime10 & severe_10$time<=ftime_12after10,]$cases)/mean(severe_10[severe_10$time>=dosetime10 & severe_10$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(severe_10[severe_10$time>=dosetime10 & severe_10$time<=ftime_18after10,]$cases)/mean(severe_10[severe_10$time>=dosetime10 & severe_10$time<=ftime_18after10,]$n)
  
  sev_red_10_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  sev_red_10_6m10 = sev_red_10_6m10*(sev_red_10_6m10 >= 0)
  sev_red_10_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  sev_red_10_12m10 = sev_red_10_12m10*(sev_red_10_12m10 >= 0)
  sev_red_10_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  sev_red_10_18m10 = sev_red_10_18m10*(sev_red_10_18m10 >= 0)
  
  # calculate the incidence numbers for age group under 5 and average reductions
  severe_u5 = nsevere[nsevere$age_group %in% age05, -c(3,5)]
  severe_u5 <- severe_u5 %>% group_by(time) %>% 
    summarise(n = sum(n), cases=sum(cases))
  
  cpp_6before = sum(severe_u5[severe_u5$time>=btime & severe_u5$time<=ftime_6before,]$cases)/mean(severe_u5[severe_u5$time>=btime & severe_u5$time<=ftime_6before,]$n)
  cpp_12before = sum(severe_u5[severe_u5$time>=btime & severe_u5$time<=ftime_12before,]$cases)/mean(severe_u5[severe_u5$time>=btime & severe_u5$time<=ftime_12before,]$n)
  cpp_18before = sum(severe_u5[severe_u5$time>=btime & severe_u5$time<=ftime_18before,]$cases)/mean(severe_u5[severe_u5$time>=btime & severe_u5$time<=ftime_18before,]$n)
  
  cpp_6after = sum(severe_u5[severe_u5$time>=dosetime1 & severe_u5$time<=ftime_6after,]$cases)/mean(severe_u5[severe_u5$time>=dosetime1 & severe_u5$time<=ftime_6after,]$n)
  cpp_12after = sum(severe_u5[severe_u5$time>=dosetime1 & severe_u5$time<=ftime_12after,]$cases)/mean(severe_u5[severe_u5$time>=dosetime1 & severe_u5$time<=ftime_12after,]$n)
  cpp_18after = sum(severe_u5[severe_u5$time>=dosetime1 & severe_u5$time<=ftime_18after,]$cases)/mean(severe_u5[severe_u5$time>=dosetime1 & severe_u5$time<=ftime_18after,]$n)
  
  sev_red_u5_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  sev_red_u5_6m = sev_red_u5_6m*(sev_red_u5_6m >= 0)
  sev_red_u5_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  sev_red_u5_12m = sev_red_u5_12m*(sev_red_u5_12m >= 0)
  sev_red_u5_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  sev_red_u5_18m = sev_red_u5_18m*(sev_red_u5_18m >= 0)
  
  cpp_6after5 = sum(severe_u5[severe_u5$time>=dosetime5 & severe_u5$time<=ftime_6after5,]$cases)/mean(severe_u5[severe_u5$time>=dosetime5 & severe_u5$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(severe_u5[severe_u5$time>=dosetime5 & severe_u5$time<=ftime_12after5,]$cases)/mean(severe_u5[severe_u5$time>=dosetime5 & severe_u5$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(severe_u5[severe_u5$time>=dosetime5 & severe_u5$time<=ftime_18after5,]$cases)/mean(severe_u5[severe_u5$time>=dosetime5 & severe_u5$time<=ftime_18after5,]$n)
  
  sev_red_u5_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  sev_red_u5_6m5 = sev_red_u5_6m5*(sev_red_u5_6m5 >= 0)
  sev_red_u5_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  sev_red_u5_12m5 = sev_red_u5_12m5*(sev_red_u5_12m5 >= 0)
  sev_red_u5_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  sev_red_u5_18m5 = sev_red_u5_18m5*(sev_red_u5_18m5 >= 0)
  
  cpp_6after10 = sum(severe_u5[severe_u5$time>=dosetime10 & severe_u5$time<=ftime_6after10,]$cases)/mean(severe_u5[severe_u5$time>=dosetime10 & severe_u5$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(severe_u5[severe_u5$time>=dosetime10 & severe_u5$time<=ftime_12after10,]$cases)/mean(severe_u5[severe_u5$time>=dosetime10 & severe_u5$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(severe_u5[severe_u5$time>=dosetime10 & severe_u5$time<=ftime_18after10,]$cases)/mean(severe_u5[severe_u5$time>=dosetime10 & severe_u5$time<=ftime_18after10,]$n)
  
  sev_red_u5_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  sev_red_u5_6m10 = sev_red_u5_6m10*(sev_red_u5_6m10 >= 0)
  sev_red_u5_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  sev_red_u5_12m10 = sev_red_u5_12m10*(sev_red_u5_12m10 >= 0)
  sev_red_u5_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  sev_red_u5_18m10 = sev_red_u5_18m10*(sev_red_u5_18m10 >= 0)
  
  # calculate the incidence numbers for all ages and the average reductions
  severe_610 = nsevere[nsevere$age_group %in% age610, -c(3,5)]
  severe_610 <- severe_610 %>% group_by(time) %>% 
    summarise(n = sum(n), cases=sum(cases))
  
  cpp_6before = sum(severe_610[severe_610$time>=btime & severe_610$time<=ftime_6before,]$cases)/mean(severe_610[severe_610$time>=btime & severe_610$time<=ftime_6before,]$n)
  cpp_12before = sum(severe_610[severe_610$time>=btime & severe_610$time<=ftime_12before,]$cases)/mean(severe_610[severe_610$time>=btime & severe_610$time<=ftime_12before,]$n)
  cpp_18before = sum(severe_610[severe_610$time>=btime & severe_610$time<=ftime_18before,]$cases)/mean(severe_610[severe_610$time>=btime & severe_610$time<=ftime_18before,]$n)
  
  cpp_6after = sum(severe_610[severe_610$time>=dosetime1 & severe_610$time<=ftime_6after,]$cases)/mean(severe_610[severe_610$time>=dosetime1 & severe_610$time<=ftime_6after,]$n)
  cpp_12after = sum(severe_610[severe_610$time>=dosetime1 & severe_610$time<=ftime_12after,]$cases)/mean(severe_610[severe_610$time>=dosetime1 & severe_610$time<=ftime_12after,]$n)
  cpp_18after = sum(severe_610[severe_610$time>=dosetime1 & severe_610$time<=ftime_18after,]$cases)/mean(severe_610[severe_610$time>=dosetime1 & severe_610$time<=ftime_18after,]$n)
  
  sev_red_610_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  sev_red_610_6m = sev_red_610_6m*(sev_red_610_6m >= 0)
  sev_red_610_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  sev_red_610_12m = sev_red_610_12m*(sev_red_610_12m >= 0)
  sev_red_610_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  sev_red_610_18m = sev_red_610_18m*(sev_red_610_18m >= 0)
  
  cpp_6after5 = sum(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_6after5,]$cases)/mean(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_12after5,]$cases)/mean(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_18after5,]$cases)/mean(severe_610[severe_610$time>=dosetime5 & severe_610$time<=ftime_18after5,]$n)
  
  sev_red_610_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  sev_red_610_6m5 = sev_red_610_6m5*(sev_red_610_6m5 >= 0)
  sev_red_610_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  sev_red_610_12m5 = sev_red_610_12m5*(sev_red_610_12m5 >= 0)
  sev_red_610_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  sev_red_610_18m5 = sev_red_610_18m5*(sev_red_610_18m5 >= 0)
  
  cpp_6after10 = sum(severe_610[severe_610$time>=dosetime10 & severe_610$time<=ftime_6after10,]$cases)/mean(severe_610[severe_610$time>=dosetime10 & severe_610$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(severe_610[severe_610$time>=dosetime10 & severe_610$time<=ftime_12after10,]$cases)/mean(severe_610[severe_610$time>=dosetime10 & severe_610$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(severe_610[severe_610$time>=dosetime10 & severe_610$time<=ftime_18after10,]$cases)/mean(severe_610[severe_610$time>=dosetime10 & severe_610$time<=ftime_18after10,]$n)
  
  sev_red_610_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  sev_red_610_6m10 = sev_red_610_6m10*(sev_red_610_6m10 >= 0)
  sev_red_610_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  sev_red_610_12m10 = sev_red_610_12m10*(sev_red_610_12m10 >= 0)
  sev_red_610_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  sev_red_610_18m10 = sev_red_610_18m10*(sev_red_610_18m10 >= 0)
  
  
  #####################################
  # calculate mortality/deaths outcomes
  #####################################
  
  # calculate the incidence numbers for all ages and the average reductions
  deaths_10 = ndeaths[ndeaths$age_group %in% age10, -c(3,5)]
  deaths_10 <- deaths_10 %>% group_by(time) %>% 
    summarise(n = sum(n), deaths=sum(deaths))
  
  cpp_6before = sum(deaths_10[deaths_10$time>=btime & deaths_10$time<=ftime_6before,]$deaths)/mean(deaths_10[deaths_10$time>=btime & deaths_10$time<=ftime_6before,]$n)
  cpp_12before = sum(deaths_10[deaths_10$time>=btime & deaths_10$time<=ftime_12before,]$deaths)/mean(deaths_10[deaths_10$time>=btime & deaths_10$time<=ftime_12before,]$n)
  cpp_18before = sum(deaths_10[deaths_10$time>=btime & deaths_10$time<=ftime_18before,]$deaths)/mean(deaths_10[deaths_10$time>=btime & deaths_10$time<=ftime_18before,]$n)
  
  cpp_6after = sum(deaths_10[deaths_10$time>=dosetime1 & deaths_10$time<=ftime_6after,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime1 & deaths_10$time<=ftime_6after,]$n)
  cpp_12after = sum(deaths_10[deaths_10$time>=dosetime1 & deaths_10$time<=ftime_12after,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime1 & deaths_10$time<=ftime_12after,]$n)
  cpp_18after = sum(deaths_10[deaths_10$time>=dosetime1 & deaths_10$time<=ftime_18after,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime1 & deaths_10$time<=ftime_18after,]$n)
  
  dea_red_10_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  dea_red_10_6m = dea_red_10_6m*(dea_red_10_6m >= 0)
  dea_red_10_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  dea_red_10_12m = dea_red_10_12m*(dea_red_10_12m >= 0)
  dea_red_10_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  dea_red_10_18m = dea_red_10_18m*(dea_red_10_18m >= 0)
  
  cpp_6after5 = sum(deaths_10[deaths_10$time>=dosetime5 & deaths_10$time<=ftime_6after5,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime5 & deaths_10$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(deaths_10[deaths_10$time>=dosetime5 & deaths_10$time<=ftime_12after5,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime5 & deaths_10$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(deaths_10[deaths_10$time>=dosetime5 & deaths_10$time<=ftime_18after5,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime5 & deaths_10$time<=ftime_18after5,]$n)
  
  dea_red_10_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  dea_red_10_6m5 = dea_red_10_6m5*(dea_red_10_6m5 >= 0)
  dea_red_10_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  dea_red_10_12m5 = dea_red_10_12m5*(dea_red_10_12m5 >= 0)
  dea_red_10_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  dea_red_10_18m5 = dea_red_10_18m5*(dea_red_10_18m5 >= 0)
  
  cpp_6after10 = sum(deaths_10[deaths_10$time>=dosetime10 & deaths_10$time<=ftime_6after10,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime10 & deaths_10$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(deaths_10[deaths_10$time>=dosetime10 & deaths_10$time<=ftime_12after10,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime10 & deaths_10$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(deaths_10[deaths_10$time>=dosetime10 & deaths_10$time<=ftime_18after10,]$deaths)/mean(deaths_10[deaths_10$time>=dosetime10 & deaths_10$time<=ftime_18after10,]$n)
  
  dea_red_10_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  dea_red_10_6m10 = dea_red_10_6m10*(dea_red_10_6m10 >= 0)
  dea_red_10_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  dea_red_10_12m10 = dea_red_10_12m10*(dea_red_10_12m10 >= 0)
  dea_red_10_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  dea_red_10_18m10 = dea_red_10_18m10*(dea_red_10_18m10 >= 0)
  
  # calculate the incidence numbers for age group under 5 and average reductions
  deaths_u5 = ndeaths[ndeaths$age_group %in% age05, -c(3,5)]
  deaths_u5 <- deaths_u5 %>% group_by(time) %>% 
    summarise(n = sum(n), deaths=sum(deaths))
  
  cpp_6before = sum(deaths_u5[deaths_u5$time>=btime & deaths_u5$time<=ftime_6before,]$deaths)/mean(deaths_u5[deaths_u5$time>=btime & deaths_u5$time<=ftime_6before,]$n)
  cpp_12before = sum(deaths_u5[deaths_u5$time>=btime & deaths_u5$time<=ftime_12before,]$deaths)/mean(deaths_u5[deaths_u5$time>=btime & deaths_u5$time<=ftime_12before,]$n)
  cpp_18before = sum(deaths_u5[deaths_u5$time>=btime & deaths_u5$time<=ftime_18before,]$deaths)/mean(deaths_u5[deaths_u5$time>=btime & deaths_u5$time<=ftime_18before,]$n)
  
  cpp_6after = sum(deaths_u5[deaths_u5$time>=dosetime1 & deaths_u5$time<=ftime_6after,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime1 & deaths_u5$time<=ftime_6after,]$n)
  cpp_12after = sum(deaths_u5[deaths_u5$time>=dosetime1 & deaths_u5$time<=ftime_12after,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime1 & deaths_u5$time<=ftime_12after,]$n)
  cpp_18after = sum(deaths_u5[deaths_u5$time>=dosetime1 & deaths_u5$time<=ftime_18after,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime1 & deaths_u5$time<=ftime_18after,]$n)
  
  dea_red_u5_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  dea_red_u5_6m = dea_red_u5_6m*(dea_red_u5_6m >= 0)
  dea_red_u5_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  dea_red_u5_12m = dea_red_u5_12m*(dea_red_u5_12m >= 0)
  dea_red_u5_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  dea_red_u5_18m = dea_red_u5_18m*(dea_red_u5_18m >= 0)
  
  cpp_6after5 = sum(deaths_u5[deaths_u5$time>=dosetime5 & deaths_u5$time<=ftime_6after5,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime5 & deaths_u5$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(deaths_u5[deaths_u5$time>=dosetime5 & deaths_u5$time<=ftime_12after5,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime5 & deaths_u5$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(deaths_u5[deaths_u5$time>=dosetime5 & deaths_u5$time<=ftime_18after5,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime5 & deaths_u5$time<=ftime_18after5,]$n)
  
  dea_red_u5_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  dea_red_u5_6m5 = dea_red_u5_6m5*(dea_red_u5_6m5 >= 0)
  dea_red_u5_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  dea_red_u5_12m5 = dea_red_u5_12m5*(dea_red_u5_12m5 >= 0)
  dea_red_u5_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  dea_red_u5_18m5 = dea_red_u5_18m5*(dea_red_u5_18m5 >= 0)
  
  cpp_6after10 = sum(deaths_u5[deaths_u5$time>=dosetime10 & deaths_u5$time<=ftime_6after10,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime10 & deaths_u5$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(deaths_u5[deaths_u5$time>=dosetime10 & deaths_u5$time<=ftime_12after10,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime10 & deaths_u5$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(deaths_u5[deaths_u5$time>=dosetime10 & deaths_u5$time<=ftime_18after10,]$deaths)/mean(deaths_u5[deaths_u5$time>=dosetime10 & deaths_u5$time<=ftime_18after10,]$n)
  
  dea_red_u5_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  dea_red_u5_6m10 = dea_red_u5_6m10*(dea_red_u5_6m10 >= 0)
  dea_red_u5_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  dea_red_u5_12m10 = dea_red_u5_12m10*(dea_red_u5_12m10 >= 0)
  dea_red_u5_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  dea_red_u5_18m10 = dea_red_u5_18m10*(dea_red_u5_18m10 >= 0)
  
  # calculate the incidence numbers for all ages and the average reductions
  deaths_610 = ndeaths[ndeaths$age_group %in% age610, -c(3,5)]
  deaths_610 <- deaths_610 %>% group_by(time) %>% 
    summarise(n = sum(n), deaths=sum(deaths))
  
  cpp_6before = sum(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_6before,]$deaths)/mean(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_6before,]$n)
  cpp_12before = sum(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_12before,]$deaths)/mean(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_12before,]$n)
  cpp_18before = sum(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_18before,]$deaths)/mean(deaths_610[deaths_610$time>=btime & deaths_610$time<=ftime_18before,]$n)
  
  cpp_6after = sum(deaths_610[deaths_610$time>=dosetime1 & deaths_610$time<=ftime_6after,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime1 & deaths_610$time<=ftime_6after,]$n)
  cpp_12after = sum(deaths_610[deaths_610$time>=dosetime1 & deaths_610$time<=ftime_12after,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime1 & deaths_610$time<=ftime_12after,]$n)
  cpp_18after = sum(deaths_610[deaths_610$time>=dosetime1 & deaths_610$time<=ftime_18after,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime1 & deaths_610$time<=ftime_18after,]$n)
  
  dea_red_610_6m = 100*(cpp_6before - cpp_6after)/ cpp_6before
  dea_red_610_6m = dea_red_610_6m*(dea_red_610_6m >= 0)
  dea_red_610_12m = 100*( cpp_12before - cpp_12after )/ cpp_12before
  dea_red_610_12m = dea_red_610_12m*(dea_red_610_12m >= 0)
  dea_red_610_18m = 100*( cpp_18before - cpp_18after )/ cpp_18before
  dea_red_610_18m = dea_red_610_18m*(dea_red_610_18m >= 0)
  
  cpp_6after5 = sum(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_6after5,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_6after5,]$n)
  cpp_12after5 = sum(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_12after5,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_12after5,]$n)
  cpp_18after5 = sum(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_18after5,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime5 & deaths_610$time<=ftime_18after5,]$n)
  
  dea_red_610_6m5 = 100*(cpp_6before - cpp_6after5)/ cpp_6before
  dea_red_610_6m5 = dea_red_610_6m5*(dea_red_610_6m5 >= 0)
  dea_red_610_12m5 = 100*( cpp_12before - cpp_12after5 )/ cpp_12before
  dea_red_610_12m5 = dea_red_610_12m5*(dea_red_610_12m5 >= 0)
  dea_red_610_18m5 = 100*( cpp_18before - cpp_18after5 )/ cpp_18before
  dea_red_610_18m5 = dea_red_610_18m5*(dea_red_610_18m5 >= 0)
  
  cpp_6after10 = sum(deaths_610[deaths_610$time>=dosetime10 & deaths_610$time<=ftime_6after10,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime10 & deaths_610$time<=ftime_6after10,]$n)
  cpp_12after10 = sum(deaths_610[deaths_610$time>=dosetime10 & deaths_610$time<=ftime_12after10,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime10 & deaths_610$time<=ftime_12after10,]$n)
  cpp_18after10 = sum(deaths_610[deaths_610$time>=dosetime10 & deaths_610$time<=ftime_18after10,]$deaths)/mean(deaths_610[deaths_610$time>=dosetime10 & deaths_610$time<=ftime_18after10,]$n)
  
  dea_red_610_6m10 = 100*(cpp_6before - cpp_6after10)/ cpp_6before
  dea_red_610_6m10 = dea_red_610_6m10*(dea_red_610_6m10 >= 0)
  dea_red_610_12m10 = 100*( cpp_12before - cpp_12after10 )/ cpp_12before
  dea_red_610_12m10 = dea_red_610_12m10*(dea_red_610_12m10 >= 0)
  dea_red_610_18m10 = 100*( cpp_18before - cpp_18after10 )/ cpp_18before
  dea_red_610_18m10 = dea_red_610_18m10*(dea_red_610_18m10 >= 0)
  
  
  ######################################
  # return results
  ######################################
  
  if(cont==FALSE) {
    # Final row with outputs to return
    return_row = cbind.data.frame(scenario_params$Scenario_Name,scenario_params$SEED,prevalence_210_before,
                                  inc_red_int_6m,inc_red_int_12m,inc_red_int_18m,inc_red_int_6m5,inc_red_int_12m5,inc_red_int_18m5,inc_red_int_6m10,inc_red_int_12m10,inc_red_int_18m10,inc_red_int_2nd6m10,inc_red_int_3rd6m10,
                                  sev_red_int_6m,sev_red_int_12m,sev_red_int_18m,sev_red_int_6m5,sev_red_int_12m5,sev_red_int_18m5,sev_red_int_6m10,sev_red_int_12m10,sev_red_int_18m10,sev_red_int_2nd6m10,sev_red_int_3rd6m10,
                                  prev_red_u5_6m,prev_red_u5_12m,prev_red_u5_18m,prev_red_u5_6m5,prev_red_u5_12m5,prev_red_u5_18m5,prev_red_u5_6m10,prev_red_u5_12m10,prev_red_u5_18m10,
                                  inc_red_u5_6m,inc_red_u5_12m,inc_red_u5_18m,inc_red_u5_6m5,inc_red_u5_12m5,inc_red_u5_18m5,inc_red_u5_6m10,inc_red_u5_12m10,inc_red_u5_18m10,
                                  sev_red_u5_6m,sev_red_u5_12m,sev_red_u5_18m,sev_red_u5_6m5,sev_red_u5_12m5,sev_red_u5_18m5,sev_red_u5_6m10,sev_red_u5_12m10,sev_red_u5_18m10,
                                  dea_red_u5_6m,dea_red_u5_12m,dea_red_u5_18m,dea_red_u5_6m5,dea_red_u5_12m5,dea_red_u5_18m5,dea_red_u5_6m10,dea_red_u5_12m10,dea_red_u5_18m10,
                                  prev_red_10_6m,prev_red_10_12m,prev_red_10_18m,prev_red_10_6m5,prev_red_10_12m5,prev_red_10_18m5,prev_red_10_6m10,prev_red_10_12m10,prev_red_10_18m10,
                                  inc_red_10_6m,inc_red_10_12m,inc_red_10_18m,inc_red_10_6m5,inc_red_10_12m5,inc_red_10_18m5,inc_red_10_6m10,inc_red_10_12m10,inc_red_10_18m10,
                                  sev_red_10_6m,sev_red_10_12m,sev_red_10_18m,sev_red_10_6m5,sev_red_10_12m5,sev_red_10_18m5,sev_red_10_6m10,sev_red_10_12m10,sev_red_10_18m10,
                                  dea_red_10_6m,dea_red_10_12m,dea_red_10_18m,dea_red_10_6m5,dea_red_10_12m5,dea_red_10_18m5,dea_red_10_6m10,dea_red_10_12m10,dea_red_10_18m10,
                                  prev_red_210_6m,prev_red_210_12m,prev_red_210_18m,prev_red_210_6m5,prev_red_210_12m5,prev_red_210_18m5,prev_red_210_6m10,prev_red_210_12m10,prev_red_210_18m10,
                                  inc_red_610_6m,inc_red_610_12m,inc_red_610_18m,inc_red_610_6m5,inc_red_610_12m5,inc_red_610_18m5,inc_red_610_6m10,inc_red_610_12m10,inc_red_610_18m10,
                                  sev_red_610_6m,sev_red_610_12m,sev_red_610_18m,sev_red_610_6m5,sev_red_610_12m5,sev_red_610_18m5,sev_red_610_6m10,sev_red_610_12m10,sev_red_610_18m10,
                                  dea_red_610_6m,dea_red_610_12m,dea_red_610_18m,dea_red_610_6m5,dea_red_610_12m5,dea_red_610_18m5,dea_red_610_6m10,dea_red_610_12m10,dea_red_610_18m10
    )
    
    
    colnames(return_row) = c("Scenario_Name","seed","prevalence_210_before",
                             "inc_red_int_6m","inc_red_int_12m","inc_red_int_18m","inc_red_int_6m5","inc_red_int_12m5","inc_red_int_18m5","inc_red_int_6m10","inc_red_int_12m10","inc_red_int_18m10","inc_red_int_2nd6m10","inc_red_int_3rd6m10",
                             "sev_red_int_6m","sev_red_int_12m","sev_red_int_18m","sev_red_int_6m5","sev_red_int_12m5","sev_red_int_18m5","sev_red_int_6m10","sev_red_int_12m10","sev_red_int_18m10","sev_red_int_2nd6m10","sev_red_int_3rd6m10",
                             "prev_red_u5_6m","prev_red_u5_12m","prev_red_u5_18m","prev_red_u5_6m5","prev_red_u5_12m5","prev_red_u5_18m5","prev_red_u5_6m10","prev_red_u5_12m10","prev_red_u5_18m10",
                             "inc_red_u5_6m","inc_red_u5_12m","inc_red_u5_18m","inc_red_u5_6m5","inc_red_u5_12m5","inc_red_u5_18m5","inc_red_u5_6m10","inc_red_u5_12m10","inc_red_u5_18m10",
                             "sev_red_u5_6m","sev_red_u5_12m","sev_red_u5_18m","sev_red_u5_6m5","sev_red_u5_12m5","sev_red_u5_18m5","sev_red_u5_6m10","sev_red_u5_12m10","sev_red_u5_18m10",
                             "dea_red_u5_6m","dea_red_u5_12m","dea_red_u5_18m","dea_red_u5_6m5","dea_red_u5_12m5","dea_red_u5_18m5","dea_red_u5_6m10","dea_red_u5_12m10","dea_red_u5_18m10",
                             "prev_red_10_6m","prev_red_10_12m","prev_red_10_18m","prev_red_10_6m5","prev_red_10_12m5","prev_red_10_18m5","prev_red_10_6m10","prev_red_10_12m10","prev_red_10_18m10",
                             "inc_red_10_6m","inc_red_10_12m","inc_red_10_18m","inc_red_10_6m5","inc_red_10_12m5","inc_red_10_18m5","inc_red_10_6m10","inc_red_10_12m10","inc_red_10_18m10",
                             "sev_red_10_6m","sev_red_10_12m","sev_red_10_18m","sev_red_10_6m5","sev_red_10_12m5","sev_red_10_18m5","sev_red_10_6m10","sev_red_10_12m10","sev_red_10_18m10",
                             "dea_red_10_6m","dea_red_10_12m","dea_red_10_18m","dea_red_10_6m5","dea_red_10_12m5","dea_red_10_18m5","dea_red_10_6m10","dea_red_10_12m10","dea_red_10_18m10",
                             "prev_red_210_6m","prev_red_210_12m","prev_red_210_18m","prev_red_210_6m5","prev_red_210_12m5","prev_red_210_18m5","prev_red_210_6m10","prev_red_210_12m10","prev_red_210_18m10",
                             "inc_red_610_6m","inc_red_610_12m","inc_red_610_18m","inc_red_610_6m5","inc_red_610_12m5","inc_red_610_18m5","inc_red_610_6m10","inc_red_610_12m10","inc_red_610_18m10",
                             "sev_red_610_6m","sev_red_610_12m","sev_red_610_18m","sev_red_610_6m5","sev_red_610_12m5","sev_red_610_18m5","sev_red_610_6m10","sev_red_610_12m10","sev_red_610_18m10",
                             "dea_red_610_6m","dea_red_610_12m","dea_red_610_18m","dea_red_610_6m5","dea_red_610_12m5","dea_red_610_18m5","dea_red_610_6m10","dea_red_610_12m10","dea_red_610_18m10"
    )
    return(return_row)}
  else{
    out_df= list("prevalence_210"=prev_210,
                 "prevalence_int1"=prev_int1,
                 "prevalence_int2"=prev_int2,
                 "prevalence_int3"=prev_int3,
                 "prevalence_agegroups"=prev_agegroups,
                 "incidence_05"=inc_05,
                 "incidence_int1"=inc_int1,
                 "incidence_int2"=inc_int2,
                 "incidence_int3"=inc_int3,
                 "incidence_agegroups"=inc_agegroups
    )
    return(out_df)
  }
  if(cont==FALSE) { 
    # Final row with outputs to return
    return_row = cbind.data.frame(scenario_params$Scenario_Name,scenario_params$SEED,prevalence_210_before,
                                  inc_red_int_6m,inc_red_int_12m,inc_red_int_18m,inc_red_int_6m5,inc_red_int_12m5,inc_red_int_18m5,inc_red_int_6m10,inc_red_int_12m10,inc_red_int_18m10,inc_red_int_2nd6m10,inc_red_int_3rd6m10,
                                  sev_red_int_6m,sev_red_int_12m,sev_red_int_18m,sev_red_int_6m5,sev_red_int_12m5,sev_red_int_18m5,sev_red_int_6m10,sev_red_int_12m10,sev_red_int_18m10,sev_red_int_2nd6m10,sev_red_int_3rd6m10,
                                  prev_red_u5_6m,prev_red_u5_12m,prev_red_u5_18m,prev_red_u5_6m5,prev_red_u5_12m5,prev_red_u5_18m5,prev_red_u5_6m10,prev_red_u5_12m10,prev_red_u5_18m10,
                                  inc_red_u5_6m,inc_red_u5_12m,inc_red_u5_18m,inc_red_u5_6m5,inc_red_u5_12m5,inc_red_u5_18m5,inc_red_u5_6m10,inc_red_u5_12m10,inc_red_u5_18m10,
                                  sev_red_u5_6m,sev_red_u5_12m,sev_red_u5_18m,sev_red_u5_6m5,sev_red_u5_12m5,sev_red_u5_18m5,sev_red_u5_6m10,sev_red_u5_12m10,sev_red_u5_18m10,
                                  dea_red_u5_6m,dea_red_u5_12m,dea_red_u5_18m,dea_red_u5_6m5,dea_red_u5_12m5,dea_red_u5_18m5,dea_red_u5_6m10,dea_red_u5_12m10,dea_red_u5_18m10,
                                  prev_red_10_6m,prev_red_10_12m,prev_red_10_18m,prev_red_10_6m5,prev_red_10_12m5,prev_red_10_18m5,prev_red_10_6m10,prev_red_10_12m10,prev_red_10_18m10,
                                  inc_red_10_6m,inc_red_10_12m,inc_red_10_18m,inc_red_10_6m5,inc_red_10_12m5,inc_red_10_18m5,inc_red_10_6m10,inc_red_10_12m10,inc_red_10_18m10,
                                  sev_red_10_6m,sev_red_10_12m,sev_red_10_18m,sev_red_10_6m5,sev_red_10_12m5,sev_red_10_18m5,sev_red_10_6m10,sev_red_10_12m10,sev_red_10_18m10,
                                  dea_red_10_6m,dea_red_10_12m,dea_red_10_18m,dea_red_10_6m5,dea_red_10_12m5,dea_red_10_18m5,dea_red_10_6m10,dea_red_10_12m10,dea_red_10_18m10,
                                  prev_red_210_6m,prev_red_210_12m,prev_red_210_18m,prev_red_210_6m5,prev_red_210_12m5,prev_red_210_18m5,prev_red_210_6m10,prev_red_210_12m10,prev_red_210_18m10,
                                  inc_red_610_6m,inc_red_610_12m,inc_red_610_18m,inc_red_610_6m5,inc_red_610_12m5,inc_red_610_18m5,inc_red_610_6m10,inc_red_610_12m10,inc_red_610_18m10,
                                  sev_red_610_6m,sev_red_610_12m,sev_red_610_18m,sev_red_610_6m5,sev_red_610_12m5,sev_red_610_18m5,sev_red_610_6m10,sev_red_610_12m10,sev_red_610_18m10,
                                  dea_red_610_6m,dea_red_610_12m,dea_red_610_18m,dea_red_610_6m5,dea_red_610_12m5,dea_red_610_18m5,dea_red_610_6m10,dea_red_610_12m10,dea_red_610_18m10
    )
    
    colnames(return_row) = c("Scenario_Name","seed","prevalence_210_before",
                             "inc_red_int_6m","inc_red_int_12m","inc_red_int_18m","inc_red_int_6m5","inc_red_int_12m5","inc_red_int_18m5","inc_red_int_6m10","inc_red_int_12m10","inc_red_int_18m10",
                             "sev_red_int_6m","sev_red_int_12m","sev_red_int_18m","sev_red_int_6m5","sev_red_int_12m5","sev_red_int_18m5","sev_red_int_6m10","sev_red_int_12m10","sev_red_int_18m10",
                             "prev_red_u5_6m","prev_red_u5_12m","prev_red_u5_18m","prev_red_u5_6m5","prev_red_u5_12m5","prev_red_u5_18m5","prev_red_u5_6m10","prev_red_u5_12m10","prev_red_u5_18m10",
                             "inc_red_u5_6m","inc_red_u5_12m","inc_red_u5_18m","inc_red_u5_6m5","inc_red_u5_12m5","inc_red_u5_18m5","inc_red_u5_6m10","inc_red_u5_12m10","inc_red_u5_18m10",
                             "sev_red_u5_6m","sev_red_u5_12m","sev_red_u5_18m","sev_red_u5_6m5","sev_red_u5_12m5","sev_red_u5_18m5","sev_red_u5_6m10","sev_red_u5_12m10","sev_red_u5_18m10",
                             "dea_red_u5_6m","dea_red_u5_12m","dea_red_u5_18m","dea_red_u5_6m5","dea_red_u5_12m5","dea_red_u5_18m5","dea_red_u5_6m10","dea_red_u5_12m10","dea_red_u5_18m10",
                             "prev_red_10_6m","prev_red_10_12m","prev_red_10_18m","prev_red_10_6m5","prev_red_10_12m5","prev_red_10_18m5","prev_red_10_6m10","prev_red_10_12m10","prev_red_10_18m10",
                             "inc_red_10_6m","inc_red_10_12m","inc_red_10_18m","inc_red_10_6m5","inc_red_10_12m5","inc_red_10_18m5","inc_red_10_6m10","inc_red_10_12m10","inc_red_10_18m10",
                             "sev_red_10_6m","sev_red_10_12m","sev_red_10_18m","sev_red_10_6m5","sev_red_10_12m5","sev_red_10_18m5","sev_red_10_6m10","sev_red_10_12m10","sev_red_10_18m10",
                             "dea_red_10_6m","dea_red_10_12m","dea_red_10_18m","dea_red_10_6m5","dea_red_10_12m5","dea_red_10_18m5","dea_red_10_6m10","dea_red_10_12m10","dea_red_10_18m10",
                             "prev_red_210_6m","prev_red_210_12m","prev_red_210_18m","prev_red_210_6m5","prev_red_210_12m5","prev_red_210_18m5","prev_red_210_6m10","prev_red_210_12m10","prev_red_210_18m10",
                             "inc_red_610_6m","inc_red_610_12m","inc_red_610_18m","inc_red_610_6m5","inc_red_610_12m5","inc_red_610_18m5","inc_red_610_6m10","inc_red_610_12m10","inc_red_610_18m10",
                             "sev_red_610_6m","sev_red_610_12m","sev_red_610_18m","sev_red_610_6m5","sev_red_610_12m5","sev_red_610_18m5","sev_red_610_6m10","sev_red_610_12m10","sev_red_610_18m10",
                             "dea_red_610_6m","dea_red_610_12m","dea_red_610_18m","dea_red_610_6m5","dea_red_610_12m5","dea_red_610_18m5","dea_red_610_6m10","dea_red_610_12m10","dea_red_610_18m10"
    )
    
    return(return_row)}else{
      out_df= list("prevalence_210"=prev_210,
                   "prevalence_int1"=prev_int1,
                   "prevalence_int2"=prev_int2,
                   "prevalence_int3"=prev_int3,
                   "incidence_int1"=inc_05,
                   "incidence_int1"=inc_int1,
                   "incidence_int2"=inc_int2,
                   "incidence_int3"=inc_int3
      )
      return(out_df)
    }
}


# setwd("/scicore/home/penny/GROUP/M3TPP/iTPP2_Tradeoffs_Apr3d_seas4_acc0.25_5to17mo")
# results_folder <- "/scicore/home/penny/GROUP/M3TPP/iTPP2_Tradeoffs_Apr3d_seas4_acc0.25_5to17mo/om/"
# param_table_file <- "/scicore/home/penny/GROUP/M3TPP/iTPP2_Tradeoffs_Apr3d_seas4_acc0.25_5to17mo/param_tab.txt"

# Wrapper for looping across all simulation results and gathering postprocessing results in a table
postprocess_OM = function(results_folder, param_table_file, final_table_dest, 
                          final_seed_table_dest, dosetime1) {
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  processed_OM_sim = NULL
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
      mtry <- try(calculate_outputs(OM_result, param_table[i,], dosetime1,cont=FALSE),
                  silent = TRUE)
      
      # Skip error or calculate outputs
      if (class(mtry) != "try-error") {
        scenario_row = calculate_outputs(OM_result, param_table[i,], dosetime1,cont=FALSE)
        processed_OM_sim = data.frame(rbind(processed_OM_sim, scenario_row),stringsAsFactors = FALSE)
      } else {
        message("Error")
      }
    }
  }
  
  
  # Summarize results over seeds and create final results tables
  
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prevalence_210_before"):length(names(processed_OM_sim) ) ]),median,na.rm=TRUE)
  
  no_seed_table = param_table[,-c(which(colnames(param_table)=="SEED"))]
  no_seed_table = unique(no_seed_table)
  final_seed_table = merge(no_seed_table, processed_OM_sim, by = c("Scenario_Name"))
  final_table = merge(no_seed_table, aggregated_OM, by = c("Scenario_Name"))
  
  # Write result tables (summarized and with seeds) to files
  write.table(final_table, final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  write.table(final_seed_table, final_seed_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
}
