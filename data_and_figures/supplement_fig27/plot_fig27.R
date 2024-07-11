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
data <- readRDS("./data_and_figures/supplement_fig27/data_fig27.rds")


# GENERATE PLOTS ----

# Define colours
cols <- c("#d9e8f0", "#ffffff", "#fffaf3", "#fff3e1", "#fee4be", "#FED18C", "#ffb66d", "#ff9021")

# Subset data to desired predictors
data <- data[data$pred == "Cumulative severe cases by age 10", ]

# Construct plot
p <- ggplot(data, aes(x = Halflife, y = Efficacy, fill = targetLabel))

p <- p + geom_tile()

p <- p + facet_grid(Access + Seasonality ~ Experiment)

p <- p + theme(panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(size = 8, face = "bold"),
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
  scale_y_continuous(breaks = seq(0.2, 1.0, 0.2),
                     labels = paste0(seq(20, 100, 20), "%"))

p <- p + labs(x = "Protection half-life (days)",
              y = "Initial efficacy")

p <- p + guides(fill = guide_legend(title = "Reduction in cumulative severe cases by age 10 vs. SMC", nrow = 2))

p



# SAVE FIGURE ----

ggsave(filename = "./data_and_figures/supplement_fig27/fig27.jpeg",
       plot = last_plot(),
       width = 8.1,
       height = 8.5,
       dpi = 400)


