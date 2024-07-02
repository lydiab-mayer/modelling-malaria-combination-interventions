# INTRO ----
#
# Visualises parameter relationships
#
# Written by Lydia Braunack-Mayer

# SETUP ----

rm(list = ls())
require(ggplot2)
require(patchwork)
require(dplyr)

# Load data
dataAB <- readRDS("./data_and_figures/manuscript_fig3/data_fig3_panelsAB.rds")
dataC <- readRDS("./data_and_figures/manuscript_fig3/data_fig3_panelC.rds")

# Load tags
tagAB <- readRDS("./data_and_figures/manuscript_fig3/label_fig3_panelsAB.rds")
tagAB
tagC <- readRDS("./data_and_figures/manuscript_fig3/label_fig3_panelC.rds")
tagC


# GENERATE PLOT FOR 5 YEAR OUTCOMES ----

## Generate plot for halflife ----

# Define colours
col <- "#468AB2"

# Construct plot
p <- ggplot(dataAB[dataAB$Parameter == "Halflife" & dataAB$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               legend.key.width = unit(1.2, "cm"), 
               legend.position = "bottom")

p <- p  + scale_x_continuous(expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 80, 20),
                     limits = c(-10, 80),
                     labels = paste0(seq(0, 80, 20), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

p <- p + labs(x = "Protection half-life (days)", y = "Median red. in age 5\ncum. cases vs SMC")


## Generate plot for efficacy ----

# Define colours
col <- "#EB5160"

# Construct plot
q <- ggplot(dataAB[dataAB$Parameter == "Efficacy" & dataAB$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

q <- q + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               legend.key.width = unit(1.2, "cm"), 
               legend.position = "bottom")

q <- q + scale_x_continuous(breaks = seq(0.2, 1.0, by = 0.2),
                            limits = c(0.2, 1.0),
                            labels = paste0(seq(20, 100, by = 20), "%"),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 80, 20),
                     limits = c(-10, 80),
                     labels = paste0(seq(0, 80, 20), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

q <- q + labs(x = "Initial efficacy", y = "")


## Generate plot for kdecay ----

# Define colours
col <- "#22223B"

# Construct plot
r <- ggplot(dataAB[dataAB$Parameter == "Kdecay" & dataAB$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

r <- r + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               legend.key.width = unit(1.2, "cm"), 
               legend.position = "bottom")

r <- r + scale_x_continuous(breaks = seq(0, 10, by = 1),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 80, 20),
                     limits = c(-10, 80),
                     labels = paste0(seq(0, 80, 20), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

r <- r + labs(x = "Decay shape", y = "")


## Construct panel A

pA <- p + q + r +
  plot_layout(guides = "collect") +
  plot_annotation(title = "A. Parameter relationships with cumulative case outcomes by five years old") &
  theme(legend.position  = "bottom",
        plot.title = element_text(family = "Times", size = 12, face = "bold", vjust = 5))



# GENERATE PLOT FOR 10 YEAR OUTCOMES ----

## Generate plot for halflife ----

# Define colours
col <- "#468AB2"

# Construct plot
p <- ggplot(dataAB[dataAB$Parameter == "Halflife" & dataAB$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               legend.key.width = unit(1.2, "cm"), 
               legend.position = "bottom")

p <- p  + scale_x_continuous(expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 30, 10),
                     limits = c(-10, 30),
                     labels = paste0(seq(0, 30, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

p <- p + labs(x = "Protection half-life (days)", y = "Median red. in age 10\ncum. cases vs SMC")


## Generate plot for efficacy ----

# Define colours
col <- "#EB5160"

# Construct plot
q <- ggplot(dataAB[dataAB$Parameter == "Efficacy" & dataAB$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

q <- q + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               legend.key.width = unit(1.2, "cm"), 
               legend.position = "bottom")

q <- q + scale_x_continuous(breaks = seq(0.2, 1.0, by = 0.2),
                            limits = c(0.2, 1.0),
                            labels = paste0(seq(20, 100, by = 20), "%"),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 30, 10),
                     limits = c(-10, 30),
                     labels = paste0(seq(0, 30, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

q <- q + labs(x = "Initial efficacy", y = "")


## Generate plot for kdecay ----

# Define colours
col <- "#22223B"

# Construct plot
r <- ggplot(dataAB[dataAB$Parameter == "Kdecay" & dataAB$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

r <- r + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 12),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               legend.key.width = unit(1.2, "cm"), 
               legend.position = "bottom")

r <- r + scale_x_continuous(breaks = seq(0, 10, by = 1),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 30, 10),
                     limits = c(-10, 30),
                     labels = paste0(seq(0, 30, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

r <- r + labs(x = "Decay shape", y = "")


## Construct panel B ----

pB <- p + q + r +
  plot_layout(guides = "collect") +
  plot_annotation(title = "B. Parameter relationships with cumulative case outcomes by ten years old") &
  theme(legend.position  = "bottom",
        plot.title = element_text(family = "Times", size = 12, face = "bold", vjust = 5))



# GENERATE PANEL C ----

# Define colours
cols <- c("#d9e8f0", "#ffffff", "#fffaf3", "#fff3e1", "#fee4be", "#FED18C")

# Subset data to desired predictors
dataC <- dataC[dataC$pred == "Cumulative severe cases by age 10", ]

# Construct plot
p <- ggplot(dataC, aes(x = Halflife, y = Efficacy, fill = targetLabel))

p <- p + geom_tile()

p <- p + facet_wrap(. ~ Experiment, ncol = 2)

p <- p + theme(panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(size = 10, face = "bold"),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.text = element_text(size = 9),
               legend.title = element_text(face = "bold", size = 10),
               legend.key.width = unit(0.8, "cm"),
               legend.margin = margin(-3, 0, 10, 0))

p <- p + scale_fill_manual(values = cols) +
  scale_y_continuous(breaks = seq(0.3, 1.0, 0.1),
                     labels = paste0(seq(30, 100, 10), "%"))

p <- p + labs(x = "Protection half-life (days)",
              y = "Initial efficacy")

p <- p + guides(fill = guide_legend(title = "Reduction in cum. severe cases by age 10 vs. SMC", nrow = 1))

## Construct panel C ----

pC <- p +
  plot_annotation(title = "C. Parameter relationships with SMC deployment") &
  theme(legend.position  = "bottom",
        plot.title = element_text(family = "Times", size = 12, face = "bold", vjust = 0))

# CONSTRUCT FINAL FIGURE ----

wrap_elements(pA) / wrap_elements(pB) / wrap_elements(pC) + plot_layout(heights = c(1.15, 1.15, 2))

ggsave(filename = "./data_and_figures/manuscript_fig3/fig3.jpeg",
       plot = last_plot(),
       width = 8.1,
       height = 9,
       dpi = 400)

# GENERATE POINT ESTIMATES FOR MANUSCRIPT TEXT ----

dataAB[dataAB$Parameter == "Kdecay" & dataAB$Outcome_age == "10 years", ] %>%
  filter(segLower == 0.6)

dataAB[dataAB$Parameter == "Kdecay" & dataAB$Outcome_age == "10 years", ] %>%
  filter(segLower == 4)

# Define target
targetOutcome <- 5

# Identify criteria
dataC %>%
  group_by(Experiment) %>%
  filter(target >= targetOutcome,
         pred == "Cumulative severe cases by age 10",
         Efficacy == 0.5) %>%
  summarise(maxHalflife = max(Halflife),
            minHalflife = min(Halflife))

dataC %>%
  group_by(Experiment) %>%
  filter(target >= targetOutcome,
         pred == "Cumulative severe cases by age 10",
         Efficacy == 0.7) %>%
  summarise(maxHalflife = max(Halflife),
            minHalflife = min(Halflife))