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
data <- readRDS("./data_and_figures/manuscript_fig5/data_fig5.rds")

# Load tag
tag <- readRDS("./data_and_figures/manuscript_fig5/label_fig5.rds")
tag


# GENERATE PLOT FOR 5 YEAR OUTCOMES ----

plotTheme <- theme(panel.border = element_blank(), 
                   panel.background = element_blank(),
                   panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor = element_blank(),
                   text = element_text(family = "Times", size = 12),
                   strip.background = element_blank(),
                   strip.text = element_text(face = "bold", size = 12),
                   axis.line = element_blank(),
                   axis.text.x = element_text(colour = "grey45"),
                   axis.text.y = element_text(colour = "grey45"),
                   axis.title.x = element_text(colour = "grey30", face = "bold", size = 12, margin = margin(t = 10)),
                   axis.title.y = element_text(colour = "grey30", face = "bold", size = 12, margin = margin(r = 10)),
                   plot.title = element_text(hjust = 0.5, face = "bold"),
                   legend.text = element_text(size = 12),
                   legend.title = element_blank(),
                   legend.key.height = unit(0.5, "cm"),
                   legend.key.width = unit(1.2, "cm"),
                   legend.key = element_rect(fill="white"), 
                   legend.position = "bottom")

## Generate plot for halflife ----

# Define colours
col <- "#468AB2"

# Construct plot
p <- ggplot(data[data$Parameter == "Halflife" & data$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

p <- p + plotTheme

p <- p  + scale_x_continuous(expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(-10, 60, 10),
                     limits = c(-10, 60),
                     labels = paste0(seq(-10, 60, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

p <- p + labs(x = "Protection half-life (days)", y = "Median reduction in age 5\ncumulative cases vs SMC") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Generate plot for efficacy ----

# Define colours
col <- "#EB5160"

# Construct plot
q <- ggplot(data[data$Parameter == "Efficacy" & data$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

q <- q + plotTheme

q <- q + scale_x_continuous(breaks = seq(0.2, 1.0, by = 0.2),
                            limits = c(0.2, 1.0),
                            labels = paste0(seq(20, 100, by = 20), "%"),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(-10, 60, 10),
                     limits = c(-10, 60),
                     labels = paste0(seq(-10, 60, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

q <- q + labs(x = "Initial efficacy", y = "") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Generate plot for kdecay ----

# Define colours
col <- "#22223B"

# Construct plot
r <- ggplot(data[data$Parameter == "Kdecay" & data$Outcome_age == "5 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

r <- r + plotTheme

r <- r + scale_x_continuous(breaks = seq(0, 10, by = 2),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(-10, 60, 10),
                     limits = c(-10, 60),
                     labels = paste0(seq(-10, 60, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

r <- r + labs(x = "Decay shape", y = "") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Construct all panels

p1 <- p + q + r +
  plot_layout(guides = "collect") +
  plot_annotation(title = "A. Blood stage parameter relationships with cumulative case outcomes by five years old") &
  theme(legend.position  = "none",
        plot.title = element_text(family = "Times", size = 10, face = "bold", vjust = 5))



# GENERATE PLOT FOR 10 YEAR OUTCOMES ----

## Generate plot for halflife ----

# Define colours
col <- "#468AB2"

# Construct plot
p <- ggplot(data[data$Parameter == "Halflife" & data$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

p <- p + plotTheme

p <- p  + scale_x_continuous(expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(-10, 30, 10),
                     limits = c(-10, 30),
                     labels = paste0(seq(-10, 30, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

p <- p + labs(x = "Protection half-life (days)", y = "Median reduction in age 10\ncumulative cases vs SMC") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Generate plot for efficacy ----

# Define colours
col <- "#EB5160"

# Construct plot
q <- ggplot(data[data$Parameter == "Efficacy" & data$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

q <- q + plotTheme

q <- q + scale_x_continuous(breaks = seq(0.2, 1.0, by = 0.2),
                            limits = c(0.2, 1.0),
                            labels = paste0(seq(20, 100, by = 20), "%"),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(-10, 30, 10),
                     limits = c(-10, 30),
                     labels = paste0(seq(-10, 30, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

q <- q + labs(x = "Initial efficacy", y = "") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Generate plot for kdecay ----

# Define colours
col <- "#22223B"

# Construct plot
r <- ggplot(data[data$Parameter == "Kdecay" & data$Outcome_age == "10 years", ], aes(x = segLower, y = median, ymin = quantile0.25, ymax = quantile0.75, linetype = Outcome, colour = Outcome, fill = Outcome)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.15, linewidth = 0.1)

r <- r + plotTheme

r <- r + scale_x_continuous(breaks = seq(0, 10, by = 2),
                            expand = expansion(mult = .05, add = 0)) +
  scale_y_continuous(breaks = seq(-10, 30, 10),
                     limits = c(-10, 30),
                     labels = paste0(seq(-10, 30, 10), "%")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_fill_manual(values = rep(col, 2), guide = "none") +
  scale_colour_manual(values = rep(col, 2), guide = "none")

r <- r + labs(x = "Decay shape", y = "") +
  guides(fill = guide_legend(override.aes = list(fill = NA)))


## Construct all panels ----

p2 <- p + q + r +
  plot_layout(guides = "collect") +
  plot_annotation(title = "B. Blood stage parameter relationships with cumulative case outcomes by ten years old") &
  theme(legend.position  = "bottom",
        plot.title = element_text(family = "Times", size = 10, face = "bold", vjust = 7))



# CONSTRUCT FINAL FIGURE ----

wrap_elements(p1) / wrap_elements(p2) + plot_layout(heights = c(0.9, 1))

ggsave(filename = "./data_and_figures/manuscript_fig5/fig5.jpeg",
       plot = last_plot(),
       width = 8,
       height = 6.5,
       dpi = 400)


# GENERATE POINT ESTIMATES FOR MANUSCRIPT TEXT ----

data[data$Parameter == "Efficacy" & data$Outcome_age == "10 years", ] %>%
  filter(median >= 0) %>%
  group_by(Outcome) %>%
  summarise(segLower = min(segLower))

data[data$Parameter == "Efficacy" & data$Outcome_age == "10 years", ] %>%
  filter(median >= 10) %>%
  group_by(Outcome) %>%
  summarise(segLower = min(segLower))
