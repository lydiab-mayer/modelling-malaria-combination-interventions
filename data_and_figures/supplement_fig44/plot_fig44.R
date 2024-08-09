# INTRO ----
#
# Visualises sensitivity results
# Written by Lydia Braunack-Mayer

# SETUP ----

rm(list = ls())
require(ggplot2)
require(patchwork)
require(dplyr)

# Load data
data <- readRDS("./data_and_figures/supplement_fig44/data_fig44.rds")


# GENERATE PLOT ----

# Filter data to select model scenarios
data <- data %>%
  filter(Seasonality == "6 month transmission season",
         Access == "10% access to care")

# Define colours
cols <- c("#468AB2", "#ffa74f", "#22223B")


# Construct plot
p <- ggplot(data, aes(x = AnnualPrev, y = median, colour = Experiment))

p <- p + geom_point(position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = quantile0.25, ymax = quantile0.75),
                position = position_dodge(width = 0.25),
                width = 0)

p <- p + facet_wrap(Outcome ~ Age, scales = "free_x")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold", size = 9.5),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.text = element_text(size = 10),
               legend.key = element_blank(),
               legend.title = element_blank(),
               legend.position = "bottom")

p <- p + scale_colour_manual(values = cols) +
  scale_y_continuous(breaks = seq(-10, 60, 10),
                     labels = paste0(seq(-10, 60, 10), "%"))

p <- p + labs(x = expression(bold(paste("Baseline annual ", bolditalic("Pf"), "PR"["2-10"]))),
              y = "Median expected reduction relative to SMC")

p



# SAVE FIGURE ----

ggsave(filename = "./data_and_figures/supplement_fig44/fig44.jpeg",
       plot = last_plot(),
       width = 8,
       height = 5,
       dpi = 400)


