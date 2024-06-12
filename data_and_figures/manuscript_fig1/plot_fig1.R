############################################################
#
# Visualises deployment scenarios
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)


# ----------------------------------------------------------
# Construct data
# ----------------------------------------------------------

# Construct baseline seasonality and deployment timing
data <- data.frame(month = 1:12,
                   seas4month = c(0.00001183, 0.0004877,	0.00879777,	0.06945048,	0.2399147,	0.3626748,	0.2399147,	0.06945048,	0.00879777,	0.0004877,	0.00001183,	0.00000013),
                   seas6month = c(0.00310474,	0.01473201,	0.04945658,	0.11746577,	0.19738934,	0.23467193,	0.19738934,	0.11746577,	0.04945658,	0.01473201,	0.00310474,	0.00046293),
                   SMC = c(NA, NA, NA, 0, 0, 0, 0, 0, NA, NA, NA, NA),
                   PEV = c(NA, NA, NA, -0.5, NA, NA, NA, NA, NA, NA, NA, NA),
                   BSV = c(NA, NA, NA, -1, NA, NA, NA, NA, NA, NA, NA, NA))
data <- merge(data, 
              c("Scenario 1.  Perfect deployment", "Scenario 2.  Imperfect seasonal coverage", "Scenario 3.  Imperfect deployment", "Scenario 4.  Random allocation"),
              by = NULL)
data <- data[!is.na(data$month), ]

# Adjust deployment timing for deployment scenarios 2, 3, 4
data[data$y == "Scenario 2.  Imperfect seasonal coverage" & data$month == 8, "SMC"] <- NA
data[data$y == "Scenario 3.  Imperfect deployment" & data$month %in% c(5, 7), "SMC"] <- NA

# Convert to long format
data <- data %>%
  pivot_longer(cols = c(SMC, PEV, BSV),
               names_to = "intervention",
               values_to = "timing") %>%
  pivot_longer(cols = c(seas4month, seas6month),
               names_to = "seasProfile",
               values_to = "seasonality")
data[data$seasProfile == "seas4month", "timing"] <- NA

# Format data
data <- data %>%
  mutate(month = as.factor(month),
         intervention = case_match(intervention,
                                   "SMC" ~ "SP-AQ administration",
                                   "PEV" ~ "Pre-erythrocytic therapeutic administration",
                                   "BSV" ~ "Blood stage therapeutic administration"),
         seasProfile = case_match(seasProfile,
                                  "seas4month" ~ "4 month seasonal profile",
                                  "seas6month" ~ "6 month seasonal profile")) %>%
  mutate(intervention = factor(intervention, levels = c("SP-AQ administration", "Pre-erythrocytic therapeutic administration", "Blood stage therapeutic administration")))


# ----------------------------------------------------------
# Generate plots
# ----------------------------------------------------------

# Define colours
cols <- c("#468AB2", "#ffa74f", "#22223B")

p <- ggplot(data = data)

p <- p + geom_line(aes(x = month, y = seasonality*20, group = seasProfile, linetype = seasProfile), colour = "black")

p <- p + geom_point(aes(x = month, y = timing, colour = intervention), shape = 3, size = 1.2)

p <- p + facet_wrap(. ~ y, scales = "free_x")

p <- p + theme(panel.border = element_rect(colour = "grey45", fill = NA), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "Times", face = "bold", size = 10),
               panel.grid = element_blank(),
               axis.line = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text.x = element_text(family = "Times", colour = "grey45", margin = margin(t = 5)),
               axis.text.y = element_blank(),
               axis.title = element_text(family = "Times", colour = "grey45", size = 10),
               legend.title = element_text(family = "Times"),
               legend.text = element_text(family = "Times"),
               legend.position = "right",
               legend.key = element_blank())

p <- p + scale_colour_manual(values = cols) +
   scale_y_continuous(limits = c(-1, 7.5))

p <- p + labs(x = "Month",
              y = "Entomological innoculation rate",
              colour = "Intervention timing",
              linetype = "Transmission seasonality")


# ----------------------------------------------------------
# Save plots
# ----------------------------------------------------------

p

ggsave(filename = paste0("./data_and_figures/manuscript_fig1/fig1.jpg"),
       plot = last_plot(),
       width = 8,
       height = 4.5,
       dpi = 400)

