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
cols <- c("#468AB2", "#22223B", "#EB5160", "#F4DBD8")

# Subset data for exemplar plot
dfOut <- dfOutAggregate[[1]] %>%
  filter(Outcome == "Uncomplicated malaria")

# Extract annotation positions
annotate5 <- dfOut %>%
  filter(AgeGroup == 5, Intervention != "No intervention") %>%
  ungroup() %>%
  select(medianValue) %>%
  pull()
annotate10 <- dfOut %>%
  filter(AgeGroup == 10, Intervention != "No intervention") %>%
  ungroup() %>%
  select(medianValue) %>%
  pull()
  

# ----------------------------------------------------------
# Generate plots
# ----------------------------------------------------------

p <- ggplot()

p <- p + geom_line(data = dfOut, aes(x = AgeGroup, y = medianValue, colour = Intervention), linewidth = 0.5) +
  geom_point(data = dfOut, aes(x = AgeGroup, y = medianValue, colour = Intervention), size = 1) +
  geom_ribbon(data = dfOut, aes(x = AgeGroup, ymin = minValue, ymax = maxValue, fill = Intervention), alpha = 0.2, linewidth = 0.3)

# Add annotations at 5 years old
p <- p + geom_segment(aes(x = 5, y = annotate5[1], xend = 5, yend = annotate5[2]),
                      linewidth = 1) +
  geom_point(aes(x = 5, y = annotate5),
             size = 4) +
  geom_segment(aes(x = 5, y = 9, xend = 5, yend = 3),
               linewidth = 0.25,
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate(geom = "label", 
           x = 5, 
           y = 11,
           label = "Reductions in cumulative incidence are\nfirst evaluated when children reach\nfive years of age",
           family = "Times",
           fill = "white",
           alpha = 0.75,
           label.size = NA)

# Add annotations at 10 years old
p <- p + geom_segment(aes(x = 10, y = annotate10[1], xend = 10, yend = annotate10[2]),
                      linewidth = 1) +
  geom_point(aes(x = 10, y = annotate10),
             size = 4) +
  geom_segment(aes(x = 10, y = 4.5, xend = 10, yend = 8.8),
               linewidth = 0.25,
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate(geom = "label", 
           x = 8.5, 
           y = 2.5,
           label = "Reductions in cumulative incidence are\n evaluated again when children reach\nten years of age",
           family = "Times",
           fill = "white",
           alpha = 0.75,
           label.size = NA)


p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "Times", face = "bold", size = 10),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               axis.text = element_text(family = "Times", colour = "grey45", margin = margin(t = 5)),
               axis.title = element_text(family = "Times", colour = "grey30", face="bold", size = 11), 
               legend.position = "bottom",
               legend.text = element_text(family = "Times", size = 11), 
               plot.title.position = "plot")

p <- p + scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, 10),
                     breaks = seq(0, 10, 1))

p <- p + labs(x = "Age (years)",
              y = "Cumulative uncomplicated cases per person",
              colour = "")


# ----------------------------------------------------------
# Save plot
# ----------------------------------------------------------

p

ggsave(filename = paste0("./data_and_figures/supplement_fig11/fig11.jpg"),
       plot = last_plot(),
       width = 8,
       height = 4,
       dpi = 400)
