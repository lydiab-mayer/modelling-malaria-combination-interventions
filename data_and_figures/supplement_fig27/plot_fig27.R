# INTRO ----
#
# Visualises sensitivity results
# Written by Lydia Braunack-Mayer

# SETUP ----

rm(list = ls())
require(ggplot2)
require(patchwork)

# Load data
data <- readRDS("/scicore/home/penny/brauna0000/M3TPP/data_and_figures/supplement_fig27/data_fig27.rds")


# GENERATE PLOT ----

# Define colours
cols <- c("#468AB2", "#EB5160", "#22223B", "#ffa74f", "#FED18C", "#4B2840", "#3b597e")
text_cols <- c("black", "black", "white", "black", "black", "white", "white")

# Set font size
fontsize <- 9

# # Subset data
# data <- data %>%
#   filter(Seasonality == "6 month",
#          Access == "10%")

# Construct plot
p <- ggplot(data, aes(x = AnnualPrev, y = T_eff, fill = parameter, label = label))

p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")

p <- p + geom_text(aes(colour = parameter),
                   position = position_stack(vjust = 0.5),
                   family = "Times",
                   size = fontsize*0.28,
                   show.legend = FALSE)

p <- p + facet_wrap(OutcomeLabel ~ Seasonality + Access, scales = "free")

p <- p + theme(panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times"),
               strip.text = element_text(face = "bold", size = 9),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title.x = element_text(colour = "grey30", face="bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face="bold", size = 10, margin = margin(r = 10)),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               plot.title.position = "plot",
               plot.tag.position = c(.5, 0),
               plot.tag = element_text(size = 8, face = "plain"),
               legend.key = element_blank(),
               legend.position = "bottom",
               legend.margin = margin(-3, 0, 10, 0))

p <- p + scale_fill_manual(values = cols) +
  scale_colour_manual(values = text_cols)

p <- p + labs(x = expression(bold(paste("Baseline annual ", bolditalic("Pf"), "PR"["2-10"]))),
              y = "Variation in outcome attributed to parameter",
              fill = "")

p <- p + guides(fill = guide_legend(ncol = 3))

p



# SAVE FIGURE ----

ggsave(filename = "/scicore/home/penny/brauna0000/M3TPP/data_and_figures/supplement_fig27/plot_fig27.jpeg",
       plot = last_plot(),
       width = 8.2,
       height = 9.5,
       dpi = 300)


