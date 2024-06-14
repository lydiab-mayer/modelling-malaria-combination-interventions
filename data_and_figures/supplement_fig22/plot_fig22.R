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

# Read in data
dfOut <- readRDS("./data_and_figures/supplement_fig22/data_fig22.rds")
dfOutAggregate <- readRDS("./data_and_figures/supplement_fig22/data_aggregate_fig22.rds")

# Read in plot label
param_label <- readRDS("./data_and_figures/supplement_fig22/label_fig22.rds")
param_label

# Define colours
cols <- c("#468AB2", "#22223B", "#EB5160", "#F4DBD8")


# ----------------------------------------------------------
# Generate plots
# ----------------------------------------------------------

p <- ggplot()

p <- p + geom_rect(data = data.frame("Outcome" = unique(dfOut[[1]]$Outcome)), aes(xmin = 0.25, xmax = 5, ymin = -Inf, ymax = Inf), fill = cols[4], alpha = 0.5) +
  geom_line(data = dfOutAggregate[[1]], aes(x = AgeGroup, y = medianValue, colour = Intervention), linewidth = 0.5) +
  geom_point(data = dfOutAggregate[[1]], aes(x = AgeGroup, y = medianValue, colour = Intervention), size = 1) +
  geom_ribbon(data = dfOutAggregate[[1]], aes(x = AgeGroup, ymin = minValue, ymax = maxValue, fill = Intervention), alpha = 0.2, linewidth = 0.3)

p <- p + facet_wrap(. ~ Outcome, scales = "free_y")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "Times", face = "bold", size = 10),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               axis.text = element_text(family = "Times", colour = "grey45", margin = margin(t = 5)),
               axis.title = element_text(family = "Times", colour = "grey30", face="bold", size = 8), 
               legend.text = element_text(family = "Times", size = 10),
               legend.key = element_blank(),
               legend.position = "bottom",
               title = element_text(family = "Times", face = "bold", size = 10),
               plot.title.position = "plot")

p <- p + scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, 1))

p <- p + labs(x = "Age (years)",
              y = "Cumulative cases\nper person",
              colour = "",
              title = "A. Combining a highly efficacious pre-erythrocytic intervention with imperfect seasonal coverage of SMC")


q <- ggplot()

q <- q + geom_rect(data = data.frame("Outcome" = unique(dfOut[[2]]$Outcome)), aes(xmin = 0.25, xmax = 5, ymin = -Inf, ymax = Inf), fill = cols[4], alpha = 0.5) +
  geom_line(data = dfOutAggregate[[2]], aes(x = AgeGroup, y = medianValue, colour = Intervention), linewidth = 0.5) +
  geom_point(data = dfOutAggregate[[2]], aes(x = AgeGroup, y = medianValue, colour = Intervention), size = 1) +
  geom_ribbon(data = dfOutAggregate[[2]], aes(x = AgeGroup, ymin = minValue, ymax = maxValue, fill = Intervention), alpha = 0.2, linewidth = 0.3)

q <- q + facet_wrap(. ~ Outcome, scales = "free_y")

q <- q + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "Times", face = "bold", size = 10),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               axis.text = element_text(family = "Times", colour = "grey45", margin = margin(t = 5)),
               axis.title = element_text(family = "Times", colour = "grey30", face="bold", size = 8), 
               legend.text = element_text(family = "Times", size = 10),
               legend.key = element_blank(),
               legend.position = "bottom",
               title = element_text(family = "Times", face = "bold", size = 10),
               plot.title.position = "plot")

q <- q + scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, 1))

q <- q + labs(x = "Age (years)",
              y = "Cumulative cases\nper person",
              colour = "",
              title = "B. Combining a partially efficacious pre-erythrocytic intervention with imperfect seasonal coverage of SMC")

# ----------------------------------------------------------
# Save plots
# ----------------------------------------------------------

p / q


ggsave(filename = paste0("./data_and_figures/supplement_fig22/fig22.jpg"),
       plot = last_plot(),
       width = 8,
       height = 5.5,
       dpi = 400)

