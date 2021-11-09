########################################################
########################################################
###                                                  ###
### STEP 6: GP-BASED GRID SEARCH OPTIMIZATION PLOTS  ###
###                                                  ###
########################################################
########################################################

##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# Library
library(hetGP)
library(Rsolnp)
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggpubr)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp = "iTPP1_ShapeParameter"

# Group and simulation folder
GROUP = "/scicore/home/penny/GROUP/M3TPP/"
SIM_FOLDER = paste0(GROUP,exp,"/")

# Predictor
pred_list = c("prev_red_allInf_allages","prev_red_allInf_210","prev_red_allInf_int",
              "prev_red_patentInf_allages","prev_red_patentInf_int","prev_red_patentInf_210",
              "clinicalInc_red_all","clinicalInc_red_05","clinicalInc_red_int",
              "nSevere_red_all","nSevere_red_05","nSevere_red_int")
pred_list = pred_list[9]


##################
### LOAD DATA ###
##################

# Get postprocessing seed file name for GP file name
effB10_name = list.files(path = paste0(GROUP,exp,"/gp/GP_grid_optimization_2/",pred_list), pattern = "effB10_", full.names = FALSE)
effB10_optim = read.table(file=paste0(GROUP,exp,"/gp/GP_grid_optimization_2/",pred_list,"/",effB10_name[1]), header = T)
effB10_scenarios = read.table(file=paste0(GROUP,exp,"/gp/GP_grid_optimization_2/",pred_list,"/",effB10_name[2]), header = T)

#################################################################
######################## GGPLOT heatmaps ######################## 

# Get colors
all_viridis_colors = viridis(100)

# Only the points in the kdecay examples
example_data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                              abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1)) &
                                (abs(effB10_scenarios$Halflife-90)==min(abs(effB10_scenarios$Halflife-90)) |
                                   abs(effB10_scenarios$Halflife-180)==min(abs(effB10_scenarios$Halflife-180)) |
                                   abs(effB10_scenarios$Halflife-270)==min(abs(effB10_scenarios$Halflife-270))) &
                                (abs(effB10_scenarios$kdecay-0.25)==min(abs(effB10_scenarios$kdecay-0.25)) |
                                   abs(effB10_scenarios$kdecay-0.5)==min(abs(effB10_scenarios$kdecay-0.5)) |
                                   abs(effB10_scenarios$kdecay-0.75)==min(abs(effB10_scenarios$kdecay-0.75)) |
                                   abs(effB10_scenarios$kdecay-1)==min(abs(effB10_scenarios$kdecay-1)) |
                                   abs(effB10_scenarios$kdecay-2)==min(abs(effB10_scenarios$kdecay-2)) |
                                   abs(effB10_scenarios$kdecay-5)==min(abs(effB10_scenarios$kdecay-5)) |
                                   abs(effB10_scenarios$kdecay-10)==min(abs(effB10_scenarios$kdecay-10)))),]

# Colors
colors = ggthemes_data$gdocs$colors$value[1:length(unique(example_data$kdecay))]

# Plot
plot0a = ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                                             abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),], 
               aes(x = kdecay, y = Halflife)) + 
  geom_tile(aes(fill = mean)) +
  geom_point(data=example_data, aes(x = kdecay, y = Halflife, color = factor(kdecay)), shape = 19, size = 3) +
  scale_fill_viridis(limits = c(0, max(100))) +
  scale_color_manual(values = colors) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & 100% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  guides(color = "none") +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plot0a

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot1_100Eff100Cov_grid_kdecaypoints.png"),
         width = 2000, height = 1777, res = 300)
dev.off()

# Plot
plot0b = ggplot() +
  geom_point(data=example_data, aes(x = kdecay, y = Halflife, color = factor(kdecay)), shape = 19, size = 3) +
  scale_fill_viridis(limits = c(0, max(100))) +
  scale_color_manual(values = colors) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & 100% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30), limits = c(30,300)) +
  guides(color = "none") +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plot0b

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot1_100Eff100Cov_kdecaypoints.png"),
         width = 2100, height = 1777, res = 300)
dev.off()

################################

# Perfect Eff and Cov
plotA = ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                                         abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),], 
               aes(x = kdecay, y = Halflife)) + 
  geom_tile(aes(fill = mean)) + 
  scale_fill_viridis(limits = c(0, max(100))) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & 100% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotA

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot1_100Eff100Cov_grid.png"),
         width = 2000, height = 1777, res = 300)
dev.off()


# Perfect Eff and Cov: contour lines
plotB = ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                                             abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),], 
               aes(x = kdecay, y = Halflife)) + 
  geom_tile(aes(fill = mean)) + 
  geom_contour(aes(z = mean), breaks = 0, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 10, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 20, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 30, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 40, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 50, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 60, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 70, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 80, size = 1, colour = "white") +
  geom_contour(aes(z = mean), breaks = 90, size = 1, colour = "white") +
  scale_fill_viridis(limits = c(0, max(100))) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & 100% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotB

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot2_100Eff100Cov_whitecontour.png"),
         width = 2000, height = 1777, res = 300)
dev.off()


# Perfect Eff and Cov: contour lines by color
plotCa = ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                                             abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),], 
               aes(x = kdecay, y = Halflife)) + 
  geom_contour(aes(z = mean, colour = factor(..level..)), breaks = seq(0,90,10), size = 1) +
  scale_color_manual(values = all_viridis_colors[seq(0,90,10)]) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & 100% coverage", color = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotCa

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot3_100Eff100Cov_linescolor.png"),
         width = 2000, height = 1777, res = 300)
dev.off()


# Perfect Eff and Cov: contour lines by color
plotCb = ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                                             abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),], 
               aes(x = kdecay, y = Halflife)) + 
  geom_contour_filled(aes(z = mean, fill = factor(..level..)), breaks = seq(0,90,10), size = 1) +
  scale_fill_manual(values = all_viridis_colors[seq(0,90,10)], drop = FALSE) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & 100% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotCb

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot4_100Eff100Cov_linesfilled.png"),
         width = 2000, height = 1777, res = 300)
dev.off()

################ uncertainty

# Perfect Eff and Cov: contour lines by color
plotCb2 = ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                                              abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),], 
                aes(x = kdecay, y = Halflife)) + 
  geom_tile(aes(fill = sd2)) + 
  scale_fill_viridis(limits = c(0, max(2.5)), option = "B") +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & 100% coverage", fill = "2 SD") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotCb2


################################

# Imperfect Eff and Cov: contour lines by color
plotD = ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-0.95)==min(abs(effB10_scenarios$Efficacy-0.95)) &
                                             abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),], 
               aes(x = kdecay, y = Halflife)) + 
  geom_contour_filled(aes(z = mean, fill = factor(..level..)), breaks = seq(0,90,10), size = 1) +
  scale_fill_manual(values = all_viridis_colors[seq(0,90,10)], drop = FALSE) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "95% initial efficacy & 100% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotD

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot5_95Eff100Cov_linesfilled.png"),
         width = 2000, height = 1777, res = 300)
dev.off()

# Imperfect Eff and Cov: contour lines by color
plotE = ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-0.9)==min(abs(effB10_scenarios$Efficacy-0.9)) &
                                             abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),], 
               aes(x = kdecay, y = Halflife)) + 
  geom_contour_filled(aes(z = mean, fill = factor(..level..)), breaks = seq(0,90,10), size = 1) +
  scale_fill_manual(values = all_viridis_colors[seq(0,90,10)], drop = FALSE) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "90% initial efficacy & 100% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotE

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot6_90Eff100Cov_linesfilled.png"),
         width = 2000, height = 1777, res = 300)
dev.off()

############################################### Varying coverage

# Subset data
effB10_scenarios_3Cov = effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                                                 (abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1)) |
                                                    abs(effB10_scenarios$Coverage-0.6)==min(abs(effB10_scenarios$Coverage-0.6)) |
                                                    abs(effB10_scenarios$Coverage-0.4)==min(abs(effB10_scenarios$Coverage-0.4)))),]
effB10_scenarios_3Cov$Coverage_cat = paste0(round(effB10_scenarios_3Cov$Coverage*100, digits = -1),"%")
unique(effB10_scenarios_3Cov$Coverage_cat)
effB10_scenarios_3Cov$Coverage_cat = factor(effB10_scenarios_3Cov$Coverage_cat, levels = c("100%","60%","40%"))

plotF = ggplot(data=effB10_scenarios_3Cov, aes(x = kdecay, y = Halflife)) + 
  facet_wrap(.~Coverage_cat) +
  geom_contour_filled(aes(z = mean, fill = factor(..level..)), breaks = seq(0,90,10), size = 1) +
  scale_fill_manual(values = all_viridis_colors[seq(0,90,10)], drop = FALSE) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & varying coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotF

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot7_100Eff_varyingcoverage.png"),
         width = 4500, height = 1600, res = 300)
dev.off()

############################################### Varying efficacy

# Subset data
effB10_scenarios_3Eff = effB10_scenarios[which(abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1)) &
                                                 (abs(effB10_scenarios$Efficacy-0.95)==min(abs(effB10_scenarios$Efficacy-0.95)) |
                                                    abs(effB10_scenarios$Efficacy-0.9)==min(abs(effB10_scenarios$Efficacy-0.9)) |
                                                    abs(effB10_scenarios$Efficacy-0.85)==min(abs(effB10_scenarios$Efficacy-0.85)))),]
effB10_scenarios_3Eff$Coverage_cat = paste0(round(effB10_scenarios_3Cov$Coverage*100, digits = -1),"%")
unique(effB10_scenarios_3Eff$Coverage_cat)
effB10_scenarios_3Eff$Coverage_cat = factor(effB10_scenarios_3Eff$Coverage_cat, levels = c("100%","60%","40%"))

plotF = ggplot(data=effB10_scenarios_3Cov, aes(x = kdecay, y = Halflife)) + 
  facet_wrap(.~Coverage_cat) +
  geom_contour_filled(aes(z = mean, fill = factor(..level..)), breaks = seq(0,90,10), size = 1) +
  scale_fill_manual(values = all_viridis_colors[seq(0,90,10)], drop = FALSE) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & varying coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotF

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot7_100Eff_varyingcoverage.png"),
         width = 4500, height = 1600, res = 300)
dev.off()


##################################### HETEROGENEITY

# Get postprocessing seed file name for GP file name
effB1000_name = list.files(path = paste0(GROUP,exp,"/gp/GP_grid_optimization_2/",pred_list), pattern = "effB1000_", full.names = FALSE)
effB1000_optim = read.table(file=paste0(GROUP,exp,"/gp/GP_grid_optimization_2/",pred_list,"/",effB1000_name[1]), header = T)
effB1000_scenarios = read.table(file=paste0(GROUP,exp,"/gp/GP_grid_optimization_2/",pred_list,"/",effB1000_name[2]), header = T)

# Combine data
data1=effB10_scenarios[which(abs(effB10_scenarios$Efficacy-1)==min(abs(effB10_scenarios$Efficacy-1)) &
                                abs(effB10_scenarios$Coverage-1)==min(abs(effB10_scenarios$Coverage-1))),]
data1$het="effB = 10"
data2=effB1000_scenarios[which(abs(effB1000_scenarios$Efficacy-1)==min(abs(effB1000_scenarios$Efficacy-1)) &
                                 abs(effB1000_scenarios$Coverage-1)==min(abs(effB1000_scenarios$Coverage-1))),]
data2$het="effB = 1000"
data3 = rbind(data1,data2)

# Perfect Eff and Cov
plotA2 = ggplot(data3, 
                  aes(x = kdecay, y = Halflife)) + 
  geom_contour_filled(aes(z = mean, fill = factor(..level..)), breaks = seq(0,90,10), size = 1) +
  facet_wrap(.~het) +
  scale_fill_manual(values = all_viridis_colors[seq(0,90,10)], drop = FALSE) +
  labs(x = "k, shape parameter value", y="Halflife (days)", title = "100% initial efficacy & 100% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(breaks = c(0,0.25,1,2,5,10)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 8)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))
plotA2

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/scenarios_plot8_100Eff100Cov_comparingHet.png"),
         width = 3300, height = 1777, res = 300)
dev.off()


# Find difference
d1 = data1 %>% 
  mutate(mean1=mean) %>% 
  select(-c(mean,sd2,nugs,het))
d2 = data2 %>% 
  mutate(mean2=mean) %>% 
  select(-c(mean,sd2,nugs,het))
d3 = merge(d1,d2,by=c("kdecay","Coverage","Halflife","Efficacy","target_range"), all = T)
d3$diff = d3$mean2 - d3$mean1

# Perfect Eff and Cov
ggplot(d3, aes(x = kdecay, y = Halflife)) + 
  geom_tile(aes(fill = diff)) 


################################################################### WHO mAb example plots

# Perfect Eff and Cov: contour lines by color
ggplot(data=effB10_scenarios[which(abs(effB10_scenarios$Coverage-0.6)==min(abs(effB10_scenarios$Coverage-0.6)) &
                                              abs(effB10_scenarios$kdecay-5)==min(abs(effB10_scenarios$kdecay-5))),], 
                aes(x = Efficacy, y = Halflife)) + 
  geom_contour_filled(aes(z = mean, fill = factor(..level..)), breaks = seq(0,90,10), size = 1) +
  scale_fill_manual(values = all_viridis_colors[seq(0,90,10)], drop = FALSE) +
  labs(x = "Initial efficacy (%)", y="Halflife (days)", title = "60% coverage", fill = "Incidence\nreduction\n(%)") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(breaks = seq(30,300,30)) +
  theme_pubclean() +
  theme(title =  element_text(size = 12)) +
  theme(axis.text = element_text(size = 12), 
        legend.position = "right",
        legend.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey"))


# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/EndPoints_Plot_Optimization_60Cov_linesfilled.png"),
         width = 2000, height = 1777, res = 300)
dev.off()

