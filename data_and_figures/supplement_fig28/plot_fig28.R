# INTRO ----
#
# Visualises parameter relationships
#
# Written by Lydia Braunack-Mayer

# SETUP ----

rm(list = ls())
require(ggplot2)
require(patchwork)

# Load data
data <- readRDS("./data_and_figures/supplement_fig28/data_fig28.rds")

# Load tag
tag <- readRDS("./data_and_figures/supplement_fig28/label_fig28.rds")
tag


# GENERATE PLOTS ----

# Define colours
cols <- c("#a3c4d8", "#d9e8f0", "#ffffff", "#fffaf3", "#fff3e1", "#fee4be")

# Subset data to desired predictors
data <- data[data$pred == "Cumulative severe cases by age 10", ]

# Construct plot
p <- ggplot(data, aes(x = Halflife, y = Efficacy, fill = targetLabel))

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
               legend.position = "bottom",
               legend.margin = margin(-3, 0, 10, 0))

p <- p + scale_fill_manual(values = cols) +
  scale_y_continuous(breaks = seq(0.3, 1.0, 0.1),
                     labels = paste0(seq(30, 100, 10), "%"))

p <- p + labs(x = "Protection half-life (days)",
              y = "Initial efficacy")

p <- p + guides(fill = guide_legend(title = "Reduction in cum. severe cases by age 10 vs. SMC", nrow = 1))

p



# EXTRACT VALUES FOR MANUSCRIPT TEST ----

# Define target
targetOutcome <- 10

# Identify criteria
data %>%
  group_by(Experiment) %>%
  filter(target >= targetOutcome,
         pred == "Cumulative severe cases by age 10",
         Efficacy == 0.7) %>%
  summarise(maxHalflife = max(Halflife),
            minHalflife = min(Halflife))

data %>%
  group_by(Experiment) %>%
  filter(target >= targetOutcome,
         pred == "Cumulative severe cases by age 10",
         Efficacy == 0.8) %>%
  summarise(maxHalflife = max(Halflife),
            minHalflife = min(Halflife))


# SAVE FIGURE ----

ggsave(filename = "./data_and_figures/supplement_fig28/plot_fig28.jpeg",
       plot = last_plot(),
       width = 8.1,
       height = 4.5,
       dpi = 400)


