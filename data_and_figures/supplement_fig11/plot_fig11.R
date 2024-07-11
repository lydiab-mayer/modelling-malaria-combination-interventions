############################################################
#
# Visualises sample product profiles
#
# Written by Lydia Braunack-Mayer
############################################################

# Setup ------------------------------------------------------------------------

rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Define colours
cols <- c("#468AB2", "#22223B", "#EB5160", "#ffa74f")

# Define weibull decay function
weibull <- function(t, L, k) exp(-(t/L)^k * log(2))



# Sample profile 1 -------------------------------------------------------------

## Define data to plot ---------------------------------------------------------

# Set up data
df <- data.frame(t = 1:730,
                 "SPAQ" = NA,
                 "PE" = NA,
                 "BS" = NA)

# Construct intervention decay profiles
SPAQ <- 1 * weibull(df$t, L = 31.1, k = 5.4)
df$SPAQ <- SPAQ
df$SPAQ[31:730] <- SPAQ[1:700]
df$SPAQ[62:730] <- SPAQ[1:669]
df$SPAQ[92:730] <- SPAQ[1:639]
df$SPAQ[123:730] <- SPAQ[1:608]
df$PE <- 0.7 * weibull(df$t, L = 150, k = 1)

# Format data
df <- df %>%
  pivot_longer(cols = SPAQ:BS,
               names_to = "intervention")
df <- df %>%
  mutate(intervention = case_match(intervention,
                                   "SPAQ" ~ "SP-AQ",
                                   "PE" ~ "Pre-erythrocytic therapeutic",
                                   "BS" ~ "Blood stage therapeutic")) %>%
  mutate(intervention = factor(intervention, levels = c("SP-AQ", "Pre-erythrocytic therapeutic", "Blood stage therapeutic")))


## Construct plot --------------------------------------------------------------

p <- ggplot()

p <- p + geom_line(data = df, aes(x = t, y = value, colour = intervention), linewidth = 1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times"),
               axis.line = element_blank(),
               axis.text = element_text(colour = "grey45", margin = margin(t = 5)),
               axis.title = element_text(colour = "grey30", face="bold", size = 10), 
               legend.text = element_text(size = 10),
               legend.key = element_blank(),
               legend.position = "bottom",
               legend.margin = margin(0, 0, 10, 0))

p <- p + scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = seq(0, 700, 100)) +
  scale_y_continuous(labels = paste0(seq(0, 100, 25), "%"))

p <- p + labs(x = "Day",
              y = "Protective efficacy",
              colour = "")

p


# Sample profile 2 -------------------------------------------------------------

## Define data to plot ---------------------------------------------------------

# Set up data
df <- data.frame(t = 1:730,
                 "SPAQ" = NA,
                 "PE" = NA,
                 "BS" = NA)

# Construct intervention decay profiles
SPAQ <- 1 * weibull(df$t, L = 31.1, k = 5.4)
df$SPAQ <- SPAQ
df$SPAQ[31:730] <- SPAQ[1:700]
df$SPAQ[62:730] <- SPAQ[1:669]
df$SPAQ[92:730] <- SPAQ[1:639]
df$SPAQ[123:730] <- SPAQ[1:608]
df$PE <- 0.5 * weibull(df$t, L = 150, k = 1)
df$BS <- 0.9 * weibull(df$t, L = 200, k = 4)

# Format data
df <- df %>%
  pivot_longer(cols = SPAQ:BS,
               names_to = "intervention")
df <- df %>%
  mutate(intervention = case_match(intervention,
                                   "SPAQ" ~ "SP-AQ",
                                   "PE" ~ "Pre-erythrocytic therapeutic",
                                   "BS" ~ "Blood stage therapeutic")) %>%
  mutate(intervention = factor(intervention, levels = c("SP-AQ", "Pre-erythrocytic therapeutic", "Blood stage therapeutic")))


## Construct plot --------------------------------------------------------------

q <- ggplot()

q <- q + geom_line(data = df, aes(x = t, y = value, colour = intervention), linewidth = 1)

q <- q + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               strip.background = element_blank(),
               panel.grid.major.y = element_line(colour = "grey80", linetype = "dotted"),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times"),
               axis.line = element_blank(),
               axis.text = element_text(colour = "grey45", margin = margin(t = 5)),
               axis.title = element_text(colour = "grey30", face="bold", size = 10), 
               legend.text = element_text(size = 10),
               legend.key = element_blank(),
               legend.position = "bottom",
               legend.margin = margin(0, 0, 10, 0))

q <- q + scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = seq(0, 700, 100)) +
  scale_y_continuous(labels = paste0(seq(0, 100, 25), "%"))

q <- q + labs(x = "Day",
              y = "Protective efficacy",
              colour = "")

q



# Construct final plot ---------------------------------------------------------

p + q +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        plot.tag = element_text(family = "Times", face = "bold"))


# Save plot --------------------------------------------------------------------

ggsave(filename = "data_and_figures/supplement_fig11/fig11.jpg",
       plot = last_plot(),
       width = 8,
       height = 3,
       dpi = 400)

