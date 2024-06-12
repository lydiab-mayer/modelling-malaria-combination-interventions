# INTRO ----
#
# Visualises sensitivity results
# Written by Lydia Braunack-Mayer

# SETUP ----

rm(list = ls())
require(ggplot2)
require(patchwork)
require(dplyr)
require(tidyr)

# Load data
dataA <- readRDS("./data_and_figures/manuscript_fig6/data_fig6_panelA.rds")
dataB <- readRDS("./data_and_figures/manuscript_fig6/data_fig6_panelB.rds")


# GENERATE PANEL A ----

# Define colours
cols <- c("#468AB2", "#EB5160", "#22223B", "#ffa74f", "#FED18C", "#4B2840", "#3b597e")
text_cols <- c("black", "black", "white", "black", "black", "white", "white")

# Set font size
fontsize <- 9

# Subset data
dataA <- dataA %>%
  filter(Seasonality == "6 month transmission season",
         Access == "10% access to care")

# Construct plot
p <- ggplot(dataA, aes(x = AnnualPrev, y = T_eff, fill = parameter, label = label))

p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")

p <- p + geom_text(aes(colour = parameter),
                   position = position_stack(vjust = 0.5),
                   family = "Times",
                   size = fontsize*0.28,
                   show.legend = FALSE)

p <- p + facet_wrap(OutcomeLabel ~ .)

p <- p + theme(panel.border = element_blank(),
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times"),
               strip.text = element_text(face = "bold", size = 12),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title.x = element_text(colour = "grey30", face="bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face="bold", size = 10, margin = margin(r = 10)),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               plot.title.position = "plot",
               plot.title = element_text(size = 12, face = "bold"),
               plot.tag.position = c(.5, 0),
               plot.tag = element_text(size = 8, face = "plain"),
               legend.key = element_blank(),
               legend.position = "bottom",
               legend.margin = margin(-3, 0, 10, 0))

p <- p + scale_fill_manual(values = cols) +
  scale_colour_manual(values = text_cols)

p <- p + labs(x = expression(bold(paste("Baseline annual ", bolditalic("Pf"), "PR"["2-10"]))),
              y = "Variation in outcome attributed\nto parameter",
              fill = "",
              title = "A")

p <- p + guides(fill = guide_legend(ncol = 3))

pA <- p


# GENERATE PANEL B ----

# Filter data to select model scenarios
dataB <- dataB %>%
  filter(Seasonality == "6 month transmission season",
         Access == "10% access to care")

# Define colours
cols <- c("#468AB2", "#ffa74f", "#22223B")


# Construct plot
p <- ggplot(dataB, aes(x = AnnualPrev, y = median, colour = Experiment))

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
               strip.text = element_text(face = "bold", size = 12),
               axis.line = element_blank(),
               axis.text.x = element_text(colour = "grey45"),
               axis.text.y = element_text(colour = "grey45"),
               axis.title.x = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(t = 10)),
               axis.title.y = element_text(colour = "grey30", face = "bold", size = 10, margin = margin(r = 10)),
               plot.title.position = "plot",
               plot.title = element_text(face = "bold", size = 12),
               legend.text = element_text(size = 11),
               legend.key = element_blank(),
               legend.title = element_blank(),
               legend.position = "bottom")

p <- p + scale_colour_manual(values = cols) +
  scale_y_continuous(breaks = seq(-20, 80, 20),
                     labels = paste0(seq(-20, 80, 20), "%"))

p <- p + labs(x = expression(bold(paste("Baseline annual ", bolditalic("Pf"), "PR"["2-10"]))),
              y = "Median expected reduction relative to SMC",
              title = "B")

pB <- p



# CONSTRUCT FINAL FIGURE ----

pA / pB

ggsave(filename = "./data_and_figures/manuscript_fig6/fig6.jpeg",
       plot = last_plot(),
       width = 8.5,
       height = 9.5,
       dpi = 400)


# GENERATE VALUES FOR MANUSCRIPT TEXT ----

data_text <- dataB %>%
  select(Seasonality:Access, Experiment, Outcome:Age, median) %>%
  filter(Experiment %in% c("Pre-erythrocytic activity alone", "Pre-erythrocytic and blood stage activity")) %>%
  pivot_wider(names_from = Experiment,
              values_from = median)

data_text <- data_text %>%
  mutate(diff = `Pre-erythrocytic and blood stage activity` - `Pre-erythrocytic activity alone`) %>%
  group_by(Age) %>%
  summarise(maxdiff = max(diff),
            mindiff = min(diff))
data_text

