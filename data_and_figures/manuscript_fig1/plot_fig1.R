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
# Construct data - panels A and B
# ----------------------------------------------------------

# Construct baseline seasonality data
dataSeason <- data.frame(month = 1:12,
                         seas4month = c(0.00001183, 0.0004877,	0.00879777,	0.06945048,	0.2399147,	0.3626748,	0.2399147,	0.06945048,	0.00879777,	0.0004877,	0.00001183,	0.00000013),
                         seas6month = c(0.00310474,	0.01473201,	0.04945658,	0.11746577,	0.19738934,	0.23467193,	0.19738934,	0.11746577,	0.04945658,	0.01473201,	0.00310474,	0.00046293))
dataSeason <- dataSeason %>%
  pivot_longer(cols = c(seas4month, seas6month),
               names_to = "seasProfile",
               values_to = "seasonality") %>%
  mutate(month = as.factor(month),
         seasProfile = case_match(seasProfile,
                                  "seas4month" ~ "4 month seasonal profile",
                                  "seas6month" ~ "6 month seasonal profile"))


# Construct deployment timing
dataDeploy <- data.frame(month = 1:12,
                         SMC = c(NA, NA, NA, 0, 0, 0, 0, 0, NA, NA, NA, NA),
                         PEV = c(NA, NA, NA, -0.5, NA, NA, NA, NA, NA, NA, NA, NA),
                         BSV = c(NA, NA, NA, -1, NA, NA, NA, NA, NA, NA, NA, NA))
dataDeploy <- merge(dataDeploy, 
                    c("Timing for scenario 1.  Perfect deployment", 
                      "Timing for scenario 2.  Imperfect seasonal coverage", 
                      "Timing for scenario 3.  Imperfect deployment", 
                      "Timing for scenario 4.  Random allocation"),
                    by = NULL)
dataDeploy <- dataDeploy[!is.na(dataDeploy$month), ]

# Adjust deployment timing for deployment scenarios 2, 3, 4
dataDeploy[dataDeploy$y == "Timing for scenario 2.  Imperfect seasonal coverage" & dataDeploy$month == 8, "SMC"] <- NA
dataDeploy[dataDeploy$y == "Timing for scenario 3.  Imperfect deployment" & dataDeploy$month %in% c(5, 7), "SMC"] <- NA

# Convert to long format
dataDeploy <- dataDeploy %>%
  pivot_longer(cols = c(SMC, PEV, BSV),
               names_to = "intervention",
               values_to = "timing")

# Format data
dataDeploy <- dataDeploy %>%
  mutate(month = as.factor(month),
         intervention = case_match(intervention,
                                   "SMC" ~ "SP-AQ",
                                   "PEV" ~ "Pre-erythrocytic therapeutic",
                                   "BSV" ~ "Blood stage therapeutic"),
         coverage = case_match(y,
                            "Timing for scenario 4.  Random allocation" ~ "30 - 100% coverage",
                            .default = "100% coverage")) %>%
  mutate(intervention = factor(intervention, levels = c("SP-AQ", "Pre-erythrocytic therapeutic", "Blood stage therapeutic")))


# ----------------------------------------------------------
# Construct data - panel C
# ----------------------------------------------------------

# Load trace plot data from figure 2
dataTrace <- readRDS("./data_and_figures/manuscript_fig2/data_aggregate_fig2.rds")

# Subset
dataTrace <- dataTrace %>%
  filter(Outcome == "Uncomplicated malaria",
         TherapeuticProfile == "Long duration, high efficacy therapeutic",
         !(Experiment == "Obj6_Scen3_LayerCounterfactual_PreErythProfiles" & Intervention == "SMC + pre-erythrocytic therapeutic"))

# Extract annotation positions
annotate5 <- dataTrace %>%
  filter(AgeGroup == 5,
         Experiment == "Obj6_Scen3_PreErythProfiles",
         Intervention != "No intervention") %>%
  ungroup() %>%
  select(medianValue) %>%
  pull()
annotate10 <- dataTrace %>%
  filter(AgeGroup == 10,
         Experiment == "Obj6_Scen3_PreErythProfiles",
         Intervention != "No intervention") %>%
  ungroup() %>%
  select(medianValue) %>%
  pull()
annotateCounterfactual <- dataTrace %>%
  filter(AgeGroup == 10,
         Intervention %in% c("Pre-erythrocytic therapeutic", "SMC + pre-erythrocytic therapeutic")) %>%
  ungroup() %>%
  select(medianValue) %>%
  pull()


# ----------------------------------------------------------
# Generate plots - panel A
# ----------------------------------------------------------

p <- ggplot(data = dataSeason)

p <- p + geom_line(aes(x = month, y = seasonality, group = seasProfile, linetype = seasProfile), colour = "black")

p <- p + geom_segment(aes(x = 7.1, y = 0.23991470, xend = 8, yend = 0.23991470), linewidth = 0.4) +
  annotate(geom = "label", 
           x = 9.5, 
           y = 0.23991470,
           label = "4 month\nseasonal profile",
           family = "Times",
           fill = "white",
           alpha = 0.75,
           label.size = NA,
           size = 3) +
  geom_segment(aes(x = 8.1, y = 0.11746577, xend = 9, yend = 0.11746577), linewidth = 0.5) +
  annotate(geom = "label", 
           x = 10.5, 
           y = 0.11746577,
           label = "6 month\nseasonal profile",
           family = "Times",
           fill = "white",
           alpha = 0.75,
           label.size = NA,
           size = 3)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "Times", face = "bold", size = 10),
               panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.minor = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(family = "Times", colour = "grey45", margin = margin(t = 5)),
               axis.text.y = element_blank(),
               axis.title = element_text(family = "Times", colour = "grey45", size = 10),
               plot.title = element_text(family = "Times", face = "bold", size = 10),
               legend.position = "none")

p <- p + labs(x = "Month",
              y = "")


# ----------------------------------------------------------
# Generate plots - panel B
# ----------------------------------------------------------

# Define colours
cols <- c("#EB5160", "#22223B", "#85B79D")

q <- ggplot(data = dataDeploy)

q <- q + geom_point(aes(x = month, y = timing, colour = intervention, shape = coverage), size = 1.5)

q <- q + facet_wrap(y ~ ., ncol = 1)

q <- q + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "Times", size = 10, hjust = 0),
               panel.grid.major.x = element_line(colour = "grey80", linetype = "dotted"),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(family = "Times", colour = "grey45", margin = margin(t = 5)),
               axis.text.y = element_blank(),
               axis.title = element_text(family = "Times", colour = "grey45", size = 10),
               legend.title = element_text(family = "Times"),
               legend.text = element_text(family = "Times"),
               legend.position = "bottom",
               legend.box = 'vertical',
               legend.key = element_rect(fill = "white"),
               legend.key.width = unit(0.25, "cm"),
               legend.margin = margin(-7, 0, 0, 0))

q <- q + scale_colour_manual(values = cols) +
  scale_y_continuous(limits = c(-1.2, 0.2)) +
  scale_shape_manual(values = c(16, 1))

q <- q + labs(x = "Month",
              y = "",
              colour = "",
              shape = "")


# ----------------------------------------------------------
# Generate plots - panel C
# ----------------------------------------------------------


# Define colours
cols <- c("#468AB2", "#22223B", "#EB5160", "#ffa74f", "#F4DBD8")

r <- ggplot()

r <- r + geom_line(data = dataTrace, aes(x = AgeGroup, y = medianValue, colour = Intervention), linewidth = 0.5) +
  geom_point(data = dataTrace, aes(x = AgeGroup, y = medianValue, colour = Intervention), size = 1) +
  geom_ribbon(data = dataTrace, aes(x = AgeGroup, ymin = minValue, ymax = maxValue, fill = Intervention), alpha = 0.2, linewidth = 0.3)

# Add annotations at 5 years old
r <- r + #geom_segment(aes(x = 5, y = annotate5[1], xend = 5, yend = annotate5[2]),
  #                     linewidth = 0.8) +
  # geom_point(aes(x = 5, y = annotate5),
  #            size = 5,
  #            shape = 45) +
  geom_segment(aes(x = 5, y = 9, xend = 5, yend = 3),
               linewidth = 0.25,
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate(geom = "label", 
           x = 5, 
           y = 10,
           label = "Reductions in cumulative incidence are\nfirst evaluated when children turn five",
           family = "Times",
           fill = "white",
           alpha = 0.75,
           label.size = NA,
           size = 3)

# Add annotations at 10 years old
r <- r + #geom_segment(aes(x = 10, y = annotate10[1], xend = 10, yend = annotate10[2]),
  #                     linewidth = 0.8) +
  # geom_point(aes(x = 10, y = annotate10),
  #            size = 5,
  #            shape = 45) +
  geom_segment(aes(x = 10, y = 13.5, xend = 10, yend = 11.5),
               linewidth = 0.25,
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate(geom = "label", 
           x = 8.6, 
           y = 14.7,
           label = "Reductions are evaluated\nagain when children\nturn ten",
           family = "Times",
           fill = "white",
           alpha = 0.75,
           label.size = NA,
           size = 3)

# Add annotations for SMC counterfactual
r <- r + #geom_segment(aes(x = 10.2, y = annotateCounterfactual[1], xend = 10.2, yend = annotateCounterfactual[2]),
  #                     linewidth = 0.8) +
  # geom_point(aes(x = 10.2, y = annotateCounterfactual),
  #            size = 5,
  #            shape = 45) +
  geom_segment(aes(x = 10, y = 3.7, xend = 10, yend = 9.5),
               linewidth = 0.25,
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate(geom = "label", 
           x = 8.6, 
           y = 2.5,
           label = "Reductions are evaluated\nboth relative to SMC,\nand to the new therapeutic",
           family = "Times",
           fill = "white",
           alpha = 0.75,
           label.size = NA,
           size = 3)

r <- r + theme(panel.border = element_rect(colour = "grey80", fill = NA), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "Times", face = "bold", size = 10),
               panel.grid.major = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.minor = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_blank(),
               axis.text = element_text(family = "Times", colour = "grey45", margin = margin(t = 5)),
               legend.position = "bottom",
               legend.key = element_rect(fill = NA),
               legend.key.width = unit(0.4, "cm"),
               legend.text = element_text(family = "Times"),
               axis.title = element_text(family = "Times", colour = "grey45", size = 10))

r <- r + scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, 10.2),
                     breaks = seq(0, 10, 1))

r <- r + labs(x = "Age (years)",
              y = "Cumulative cases per person",
              colour = "") +
  guides(colour = guide_legend(nrow = 2))


# ----------------------------------------------------------
# Construct final plot
# ----------------------------------------------------------

layout <- "
AACC
BBCC
"
p + q + r + 
  plot_layout(design = layout,
              heights = c(1, 1.75)) +
  plot_annotation(tag_levels = "A")  & 
  theme(plot.tag = element_text(family = "Times", face = "bold"))

ggsave(filename = paste0("./data_and_figures/manuscript_fig1/fig1.jpg"),
       plot = last_plot(),
       width = 8.2,
       height = 5.5,
       dpi = 400)

