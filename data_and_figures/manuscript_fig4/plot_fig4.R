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
data <- readRDS("./data_and_figures/manuscript_fig4/data_fig4.rds")

# Load tags
(tag <- readRDS("./data_and_figures/manuscript_fig4/label_fig4.rds"))


# GENERATE PLOT ----

# Define colours
cols <- c("#d9e8f0", "#ffffff", "#fffaf3", "#fff3e1", "#fee4be", "#FED18C")

# Subset data to desired predictors
plot <- data[data$pred == "Cumulative severe cases by age 10", ]

# Construct plot
p <- ggplot(plot, aes(x = Halflife, y = Efficacy, fill = targetLabel))

p <- p + geom_tile()

p <- p + facet_wrap(. ~ Experiment, ncol = 2)

p <- p + theme(panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times", size = 12),
               strip.background = element_blank(),
               strip.text = element_text(size = 12, face = "bold"),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 12, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 12, margin = margin(r = 10)),
               legend.text = element_text(size = 12),
               legend.title = element_text(face = "bold", size = 12),
               legend.key.width = unit(0.8, "cm"),
               legend.margin = margin(-3, 0, 0, 0),
               legend.position = "bottom")

p <- p + scale_fill_manual(values = cols) +
  scale_y_continuous(breaks = seq(0.2, 1.0, 0.2),
                     labels = paste0(seq(20, 100, 20), "%"))

p <- p + labs(x = "Protection half-life (days)",
              y = "Initial efficacy")

p <- p + guides(fill = guide_legend(title = "Reduction in cumulative severe\ncases by age 10 vs SMC", nrow = 1))

p

# WRITE TO FILE ----

ggsave(filename = "./data_and_figures/manuscript_fig4/fig4.pdf",
       plot = last_plot(),
       width = 8.1,
       height = 5,
       dpi = 400)

# GENERATE POINT ESTIMATES FOR MANUSCRIPT TEXT ----

# Define target
targetOutcome <- 5

# Identify criteria
data %>%
  group_by(Experiment) %>%
  filter(target >= targetOutcome,
         pred == "Cumulative severe cases by age 10",
         Efficacy == 0.5) %>%
  summarise(maxHalflife = max(Halflife),
            minHalflife = min(Halflife))

# Match criteria to impact at 5 years old
data %>%
  group_by(Experiment) %>%
  filter(Halflife >= 237,
         pred == "Cumulative severe cases by age 5",
         Efficacy == 0.5) %>%
  summarise(maxImpact = max(mean),
            minImpact = min(mean))

data %>%
  group_by(Experiment) %>%
  filter(target >= targetOutcome,
         pred == "Cumulative severe cases by age 10",
         Efficacy == 0.7) %>%
  summarise(maxHalflife = max(Halflife),
            minHalflife = min(Halflife))