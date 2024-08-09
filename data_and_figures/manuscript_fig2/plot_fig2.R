############################################################
#
# Visualises cumulative cases per person year
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

# Read in data
dfOut <- readRDS("./data_and_figures/manuscript_fig2/data_fig2.rds")
dfOutAggregate <- readRDS("./data_and_figures/manuscript_fig2/data_aggregate_fig2.rds")

# Read in plot label
param_label <- readRDS("./data_and_figures/manuscript_fig2/label_fig2.rds")
param_label

# Define colours
cols <- c("#468AB2", "#22223B", "#EB5160", "#ffa74f", "#F4DBD8")

# ----------------------------------------------------------
# Generate plots
# ----------------------------------------------------------

p <- ggplot()

p <- p + geom_rect(data = data.frame("Outcome" = unique(dfOut$Outcome)), aes(xmin = 0.25, xmax = 5, ymin = -Inf, ymax = Inf), fill = cols[5], alpha = 0.5) +
  geom_line(data = dfOutAggregate, aes(x = AgeGroup, y = medianValue, colour = Intervention), linewidth = 0.5) +
  geom_point(data = dfOutAggregate, aes(x = AgeGroup, y = medianValue, colour = Intervention), size = 1) +
  geom_ribbon(data = dfOutAggregate, aes(x = AgeGroup, ymin = minValue, ymax = maxValue, fill = Intervention), alpha = 0.2, linewidth = 0.3)

p <- p + facet_wrap(TherapeuticProfile ~ Outcome, scales = "free_y")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "Times", face = "bold", size = 9),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               axis.text = element_text(family = "Times", colour = "grey45", margin = margin(t = 5)),
               axis.title = element_text(family = "Times", colour = "grey30", face="bold", size = 8), 
               legend.position = "bottom",
               legend.key = element_rect(fill = NA),
               legend.text = element_text(family = "Times"),
               title = element_text(family = "Times", face = "bold", size = 9),
               plot.title.position = "plot")

p <- p + scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, 1))

p <- p + labs(x = "Age (years)",
              y = "Cumulative cases per person",
              colour = "")

# ----------------------------------------------------------
# Save plots
# ----------------------------------------------------------

p

ggsave(filename = paste0("./data_and_figures/manuscript_fig2/fig2.jpg"),
       plot = last_plot(),
       width = 8,
       height = 5.5,
       dpi = 400)

# ----------------------------------------------------------
# Calculate point estimates for manuscript text
# ----------------------------------------------------------

# Exemplar profile A, reductions at 5 years
baseline <- dfOutAggregate %>%
  filter(Intervention == "SMC",
         TherapeuticProfile == "Long duration, high efficacy therapeutic",
         AgeGroup == 5)
intervention <- dfOutAggregate %>%
  filter(Intervention == "SMC + pre-liver stage therapeutic",
         TherapeuticProfile == "Long duration, high efficacy therapeutic",
         AgeGroup == 5)
(baseline$medianValue - intervention$medianValue) / baseline$medianValue

# Exemplar profile A, reductions at 10 years
baseline <- dfOutAggregate %>%
  filter(Intervention == "SMC",
         TherapeuticProfile == "Long duration, high efficacy therapeutic",
         AgeGroup == 10)
intervention <- dfOutAggregate %>%
  filter(Intervention == "SMC + pre-liver stage therapeutic",
         TherapeuticProfile == "Long duration, high efficacy therapeutic",
         AgeGroup == 10)
(baseline$medianValue - intervention$medianValue) / baseline$medianValue

# Exemplar profile B, reductions at 5 years
baseline <- dfOutAggregate %>%
  filter(Intervention == "SMC",
         TherapeuticProfile == "Short duration, moderate efficacy therapeutic",
         AgeGroup == 5)
intervention <- dfOutAggregate %>%
  filter(Intervention == "SMC + pre-liver stage therapeutic",
         TherapeuticProfile == "Short duration, moderate efficacy therapeutic",
         AgeGroup == 5)
(baseline$medianValue - intervention$medianValue) / baseline$medianValue

# Exemplar profile B, years to 0% reduction
baseline <- dfOutAggregate %>%
  filter(Intervention == "SMC",
         TherapeuticProfile == "Short duration, moderate efficacy therapeutic",
         Outcome == "Uncomplicated malaria")
intervention <- dfOutAggregate %>%
  filter(Intervention == "SMC + pre-liver stage therapeutic",
         TherapeuticProfile == "Short duration, moderate efficacy therapeutic",
         Outcome == "Uncomplicated malaria")
reduction <- (baseline$medianValue - intervention$medianValue) / baseline$medianValue
(reduction <- data.frame(AgeGroup = baseline$AgeGroup,
                         reduction = reduction))
