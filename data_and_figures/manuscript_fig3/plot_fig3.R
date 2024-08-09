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
dataA <- readRDS("./data_and_figures/manuscript_fig3/data_fig3_panelA.rds")
dataB <- readRDS("./data_and_figures/manuscript_fig3/data_fig3_panelB.rds")

# Load tags
(tagA <- readRDS("./data_and_figures/manuscript_fig3/label_fig3_panelA.rds"))
(tagB <- readRDS("./data_and_figures/manuscript_fig3/label_fig3_panelB.rds"))


# GENERATE PLOT FOR 5 YEAR OUTCOMES ----

plotTheme <- theme(panel.border = element_blank(), 
                   panel.background = element_blank(),
                   panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor = element_blank(),
                   text = element_text(family = "Times", size = 9),
                   strip.background = element_blank(),
                   strip.text = element_text(face = "bold", size = 12),
                   axis.line = element_blank(),
                   axis.text.x = element_text(colour = "grey45"),
                   axis.text.y = element_text(colour = "grey45"),
                   axis.title.x = element_text(colour = "grey30", face = "bold", size = 9, margin = margin(t = 10)),
                   axis.title.y = element_text(colour = "grey30", face = "bold", size = 9, margin = margin(r = 10)),
                   plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.position = "none")

## Generate plot for halflife ----

# Define colours
col <- "#468AB2"

# Construct plot
p <- ggplot(dataA[dataA$Parameter == "Halflife" & dataA$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

p <- p + plotTheme

p <- p  + scale_x_continuous(expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 80, 20),
                     limits = c(-10, 80),
                     labels = paste0(seq(0, 80, 20), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

p <- p + labs(x = "Protection half-life (days)", y = "Median reduction\nin age 5 cumulative\ncases vs SMC")


## Generate plot for efficacy ----

# Define colours
col <- "#EB5160"

# Construct plot
q <- ggplot(dataA[dataA$Parameter == "Efficacy" & dataA$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

q <- q + plotTheme

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
r <- ggplot(dataA[dataA$Parameter == "Kdecay" & dataA$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

r <- r + plotTheme

r <- r + scale_x_continuous(breaks = seq(0, 10, by = 2),
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
  plot_annotation(title = "A. Parameter relationships with cumulative case outcomes") &
  theme(plot.title = element_text(family = "Times", size = 10, face = "bold", vjust = 8))



# GENERATE PLOT FOR 10 YEAR OUTCOMES ----

plotTheme <- theme(panel.border = element_blank(), 
                   panel.background = element_blank(),
                   panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor = element_blank(),
                   text = element_text(family = "Times", size = 9),
                   strip.background = element_blank(),
                   strip.text = element_text(face = "bold", size = 9),
                   axis.line = element_blank(),
                   axis.text.x = element_text(colour = "grey45"),
                   axis.text.y = element_text(colour = "grey45"),
                   axis.title.x = element_text(colour = "grey30", face = "bold", size = 9, margin = margin(t = 10)),
                   axis.title.y = element_text(colour = "grey30", face = "bold", size = 9, margin = margin(r = 10)),
                   plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.text = element_text(size = 9),
                   legend.title = element_blank(),
                   legend.key.height = unit(0.5, "cm"),
                   legend.key.width = unit(1.2, "cm"),
                   legend.key = element_rect(fill="white"),
                   legend.position = "bottom",
                   legend.margin = margin(-3, 0, 0, 0))


## Generate plot for halflife ----

# Define colours
col <- "#468AB2"

# Construct plot
p <- ggplot(dataA[dataA$Parameter == "Halflife" & dataA$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

p <- p + plotTheme

p <- p  + scale_x_continuous(expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 30, 10),
                     limits = c(-10, 30),
                     labels = paste0(seq(0, 30, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

p <- p + labs(x = "Protection half-life (days)", y = "Median reduction\nin age 10 cumulative\ncases vs SMC") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Generate plot for efficacy ----

# Define colours
col <- "#EB5160"

# Construct plot
q <- ggplot(dataA[dataA$Parameter == "Efficacy" & dataA$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

q <- q + plotTheme

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

q <- q + labs(x = "Initial efficacy", y = "") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Generate plot for kdecay ----

# Define colours
col <- "#22223B"

# Construct plot
r <- ggplot(dataA[dataA$Parameter == "Kdecay" & dataA$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

r <- r + plotTheme

r <- r + scale_x_continuous(breaks = seq(0, 10, by = 2),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(0, 30, 10),
                     limits = c(-10, 30),
                     labels = paste0(seq(0, 30, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

r <- r + labs(x = "Decay shape", y = "") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Construct panel B ----

pB <- p + q + r +
  plot_layout(guides = "collect") &
  theme(legend.position  = "bottom")



# GENERATE PANEL C ----

# Define colours
cols <- c("#d9e8f0", "#ffffff", "#fffaf3", "#fff3e1", "#fee4be", "#FED18C")

# Subset data to desired predictors
dataB <- dataB[dataB$pred == "Cumulative severe cases by age 10", ]

# Construct plot
p <- ggplot(dataB, aes(x = Halflife, y = Efficacy, fill = targetLabel))

p <- p + geom_tile()

p <- p + facet_wrap(. ~ Experiment, ncol = 2)

p <- p + theme(panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times", size = 9),
               strip.background = element_blank(),
               strip.text = element_text(size = 9, face = "bold"),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 9, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 9, margin = margin(r = 10)),
               legend.text = element_text(size = 9),
               legend.title = element_text(face = "bold", size = 9),
               legend.key.width = unit(0.8, "cm"),
               legend.margin = margin(-3, 0, 0, 0))

p <- p + scale_fill_manual(values = cols) +
  scale_y_continuous(breaks = seq(0.2, 1.0, 0.2),
                     labels = paste0(seq(20, 100, 20), "%"))

p <- p + labs(x = "Protection half-life (days)",
              y = "Initial efficacy")

p <- p + guides(fill = guide_legend(title = "Reduction in cumulative severe\ncases by age 10 vs SMC", nrow = 1))

## Construct panel C ----

pC <- p +
  plot_annotation(title = "B. Parameter relationships with SMC deployment") &
  theme(legend.position  = "bottom",
        plot.title = element_text(family = "Times", size = 10, face = "bold", vjust = 0))

# CONSTRUCT FINAL FIGURE ----

wrap_elements(pA) / wrap_elements(pB) / wrap_elements(pC) + plot_layout(heights = c(1, 1.2, 2.5))

ggsave(filename = "./data_and_figures/manuscript_fig3/fig3.jpeg",
       plot = last_plot(),
       width = 8.1,
       height = 8.5,
       dpi = 400)

# GENERATE POINT ESTIMATES FOR MANUSCRIPT TEXT ----

dataA[dataA$Parameter == "Kdecay" & dataA$Outcome_age == "10 years", ] %>%
  filter(segLower == 0.6)

dataA[dataA$Parameter == "Kdecay" & dataA$Outcome_age == "10 years", ] %>%
  filter(segLower == 4)

# Define target
targetOutcome <- 5

# Identify criteria
dataB %>%
  group_by(Experiment) %>%
  filter(target >= targetOutcome,
         pred == "Cumulative severe cases by age 10",
         Efficacy == 0.5) %>%
  summarise(maxHalflife = max(Halflife),
            minHalflife = min(Halflife))

dataB %>%
  group_by(Experiment) %>%
  filter(target >= targetOutcome,
         pred == "Cumulative severe cases by age 10",
         Efficacy == 0.7) %>%
  summarise(maxHalflife = max(Halflife),
            minHalflife = min(Halflife))