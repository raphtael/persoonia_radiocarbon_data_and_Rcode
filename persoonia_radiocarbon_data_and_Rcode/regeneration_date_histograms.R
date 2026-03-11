################################################################################
# Plot histograms of reconstructed regeneration date
# Reconstruct age from size using DBH allometry where available,
# and height allometry otherwise
# Input: trees_recruitment_figure_data.csv
################################################################################

# Load libraries
library(data.table)
library(ggplot2)
library(cowplot)

################################################################################
# Read data
dt <- fread("data_tree_regeneration_date.csv")

################################################################################
# Reconstruct tree age from size (Weibull models)
dt[, age := 4 + 58.07545 * (1 - exp(-0.02031 * dbh_cm^1.54540))]
dt[is.na(dbh_cm) & !is.na(height_m),
   age := 60.125140 * (1 - exp(-0.001396 * height_m^3.272210))]

################################################################################
# Compute regeneration date
dt[, regeneration_date := year - age]

################################################################################
# Fire panel
dt_fire <- dt[disturbance_category == "fire"]
data_text <- data.table(
  fire_disturbance_type = "Fire Severity 1 - Crown Burn",
  x = c(1979, 2005),
  y = c(1.8, 1.8),
  label = c(1983, 2009)
)
fig_fire <- ggplot(dt_fire, aes(regeneration_date)) +
  geom_histogram(binwidth = 2) +
  facet_wrap(~fire_disturbance_type, scales = "free_y", ncol = 1, drop = TRUE) +
  geom_vline(xintercept = c(1983, 2009), colour = "red") +
  geom_text(data = data_text, aes(x = x, y = y, label = label), colour = "red") +
  scale_x_continuous(breaks = seq(1950, 2030, by = 10)) +
  coord_cartesian(xlim = c(1955, 2021)) +
  labs(x = "Regeneration date", y = "Count") +
  theme_bw()
fig_fire

################################################################################
# Logging panel
dt_logging <- dt[disturbance_category == "logging"]
fig_logging <- ggplot(dt_logging, aes(regeneration_date)) +
  geom_histogram(binwidth = 2) +
  facet_wrap(~logging_decade, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = seq(1950, 2030, by = 10)) +
  coord_cartesian(xlim = c(1955, 2021)) +
  labs(x = "Regeneration date", y = "Count") +
  theme_bw()
fig_logging

################################################################################
# Combined figure
fig_regeneration_date_histograms <- plot_grid(
  fig_fire + ggtitle("a."),
  fig_logging + ggtitle("b."),
  ncol = 2,
  align = "v"
)
fig_regeneration_date_histograms


