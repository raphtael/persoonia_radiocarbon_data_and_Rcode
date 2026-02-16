################################################################################
# Persoonia arborea phenology
# Flowering and fruiting GAM analyses
# Input: data_phenology.csv
################################################################################

# Load libraries
library(mgcv)
library(ggplot2)
library(viridis)

# Utility function
inv_logit <- function(x) 1 / (1 + exp(-x))

# Plot settings
month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Load data
phenology <- read.csv("data_phenology.csv")

################################################################################
# Flowering versus day of year

data_flowering <- subset(phenology, 
                         !is.na(flower_binary) & year >= 2020 & dbh_cm >= 1)

gam_flowering <- gam(
  flower_binary ~ s(day_of_the_year, bs = "cc", k = 6),
  knots = list(day_of_the_year = c(1, 365)),
  family = binomial,
  method = "REML",
  data = data_flowering
)

pred_flowering <- data.frame(day_of_the_year = 1:365)
fits <- predict(gam_flowering, newdata = pred_flowering, 
                type = "link", se.fit = TRUE)

pred_flowering$mean <- inv_logit(fits$fit)
pred_flowering$ymin <- inv_logit(fits$fit - 1.96 * fits$se.fit)
pred_flowering$ymax <- inv_logit(fits$fit + 1.96 * fits$se.fit)

fig_flowering <- ggplot(data_flowering, aes(x = day_of_the_year)) +
  geom_ribbon(data = pred_flowering,
              aes(ymin = ymin, ymax = ymax),
              alpha = 0.3,
              fill = viridis(option = "cividis", n = 6)[2]) +
  geom_line(data = pred_flowering,
            aes(y = mean),
            linewidth = 0.8) +
  geom_rug(data = subset(data_flowering, flower_binary == 1),
           sides = "t", colour = "grey70", linewidth = 0.3) +
  geom_rug(data = subset(data_flowering, flower_binary == 0),
           sides = "b", colour = "grey70", linewidth = 0.3) +
  scale_x_continuous(breaks = month_breaks,
                     labels = month_labels,
                     limits = c(1, 365)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "", y = "Probability of flowering") +
  theme_bw()

fig_flowering

################################################################################
# Fruiting versus day of year (by tree size)

data_fruiting <- subset(phenology,
                        !is.na(fruit_binary) & year >= 2020 & dbh_cm > 0)

gam_fruiting <- gam(
  fruit_binary ~ s(day_of_the_year, bs = "cc", k = 6) + log(dbh_cm),
  knots = list(day_of_the_year = c(1, 365)),
  family = binomial,
  method = "REML",
  data = data_fruiting
)

pred_fruiting_by_size <- expand.grid(
  day_of_the_year = 1:365,
  dbh_cm = c(1, 5, 10, 20, 30, 40)
)

fits <- predict(gam_fruiting, newdata = pred_fruiting_by_size,
                type = "link", se.fit = TRUE)

pred_fruiting_by_size$mean <- inv_logit(fits$fit)
pred_fruiting_by_size$ymin <- inv_logit(fits$fit - 1.96 * fits$se.fit)
pred_fruiting_by_size$ymax <- inv_logit(fits$fit + 1.96 * fits$se.fit)

fig_fruiting_by_size <- ggplot(data_fruiting, aes(x = day_of_the_year)) +
  geom_ribbon(data = pred_fruiting_by_size,
              aes(ymin = ymin, ymax = ymax,
                  fill = factor(dbh_cm),
                  group = factor(dbh_cm)),
              alpha = 0.25) +
  geom_line(data = pred_fruiting_by_size,
            aes(y = mean,
                colour = factor(dbh_cm),
                group = factor(dbh_cm)),
            linewidth = 0.8) +
  geom_rug(data = subset(data_fruiting, fruit_binary == 1),
           sides = "t", colour = "grey70", linewidth = 0.3) +
  geom_rug(data = subset(data_fruiting, fruit_binary == 0),
           sides = "b", colour = "grey70", linewidth = 0.3) +
  scale_x_continuous(breaks = month_breaks,
                     labels = month_labels,
                     limits = c(1, 365)) +
  scale_colour_viridis_d(option = "cividis", name = "DBH (cm)") +
  scale_fill_viridis_d(option = "cividis", name = "DBH (cm)") +
  labs(x = "", y = "Probability of fruiting") +
  theme_bw()

fig_fruiting_by_size

################################################################################
# Fruiting versus tree size (by season)

pred_fruiting_by_dbh <- expand.grid(
  day_of_the_year = month_breaks[c(4, 9)],  # April and September
  dbh_cm = seq(0.1, 44, by = 0.1)
)

fits <- predict(gam_fruiting, newdata = pred_fruiting_by_dbh, se.fit = TRUE)

pred_fruiting_by_dbh$mean <- inv_logit(fits$fit)
pred_fruiting_by_dbh$ymin <- inv_logit(fits$fit - 1.96 * fits$se.fit)
pred_fruiting_by_dbh$ymax <- inv_logit(fits$fit + 1.96 * fits$se.fit)
pred_fruiting_by_dbh$month <- factor(pred_fruiting_by_dbh$day_of_the_year,
                                      levels = month_breaks[c(4, 9)],
                                      labels = c("April", "September"))

fig_fruiting_by_dbh <- ggplot(data_fruiting, aes(x = dbh_cm)) +
  geom_ribbon(data = pred_fruiting_by_dbh,
              aes(ymin = ymin, ymax = ymax, fill = month),
              alpha = 0.3) +
  geom_line(data = pred_fruiting_by_dbh,
            aes(y = mean, colour = month)) +
  geom_rug(data = subset(data_fruiting, fruit_binary == 1),
           sides = "t", colour = "grey70", linewidth = 0.3) +
  geom_rug(data = subset(data_fruiting, fruit_binary == 0),
           sides = "b", colour = "grey70", linewidth = 0.3) +
  scale_colour_manual(values = c("lightblue2", "red3"), name = "Month") +
  scale_fill_manual(values = c("lightblue2", "red3"), name = "Month") +
  labs(x = "DBH (cm)", y = "Probability of fruiting") +
  theme_bw()

fig_fruiting_by_dbh