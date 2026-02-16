################################################################################
# Age-size allometry calibration
# Fit multiple growth models and compare goodness of fit
# Inputs: data_age_dbh_allometry.csv, data_age_height_allometry.csv
################################################################################

# Load libraries
library(ggplot2)
library(minpack.lm)
library(gsl)
library(data.table)

################################################################################
# Utility functions

# Compute goodness of fit metrics
compute_gof <- function(model, response) {
  res <- residuals(model)
  pred <- fitted(model)
  obs <- response
  r2 <- cor(pred, obs)^2
  rmse <- sd(res)
  cv <- rmse / mean(obs)
  aic <- AIC(model)
  return(data.frame(aic = round(aic, 1), r2 = round(r2, 2), 
                    rmse = round(rmse, 2), cv = round(cv, 2)))
}

# Hypergeometric function for Bontemps-Duplat model
# From: https://stats.stackexchange.com/questions/33451
gauss_2f1 <- function(a, b, c, x) {
  hyperg_2F1(c - a, b, c, 1 - 1 / (1 - x)) / (1 - x)^b
}

################################################################################
# Age versus DBH analysis

age_dbh <- read.csv("data_age_dbh_allometry.csv")

# Fit nonlinear models with inverse variance weights
weights_dbh <- 1 / age_dbh$se_age_yr^2

# Power function
nls_power <- nlsLM(
  age_yr ~ 4 + a * dbh_cm^b,
  start = c(a = 0.5, b = 0.7),
  data = age_dbh,
  weights = weights_dbh
)

# Saturating functions
nls_michaelis_menten <- nlsLM(
  age_yr ~ 4 + a * dbh_cm / (b + dbh_cm),
  start = c(a = 100, b = 10),
  data = age_dbh,
  weights = weights_dbh
)

nls_mitscherlich <- nlsLM(
  age_yr ~ 4 + a * (1 - exp(-log(2) / b * dbh_cm)),
  start = c(a = 100, b = 10),
  data = age_dbh,
  weights = weights_dbh
)

# Sigmoid saturating functions
nls_chapman_richards <- nlsLM(
  age_yr ~ 4 + a * (1 - exp(-log(2) / b * dbh_cm))^pow,
  start = c(a = 100, b = 10, pow = 1),
  data = age_dbh,
  weights = weights_dbh
)

nls_gompertz <- nlsLM(
  age_yr ~ 4 + a * exp(-b * exp(-d * dbh_cm)),
  start = c(a = 10, b = 1, d = 1),
  data = age_dbh,
  weights = weights_dbh
)

nls_weibull <- nlsLM(
  age_yr ~ 4 + a * (1 - exp(-b * dbh_cm^pow)),
  start = c(a = 100, b = 1, pow = 1),
  data = age_dbh,
  weights = weights_dbh
)

# Bontemps-Duplat models (sigmoid converging to power function)
nls_bontemps_3p <- nlsLM(
  age_yr ~ 4 + (exp(K) * exp(R) * log(1 + (-1 + exp(m2)) * (dbh_cm / exp(K))^exp(m2))) / 
    (exp(m2) - 1),
  start = c(K = log(4), R = log(4), m2 = log(2)),
  data = age_dbh,
  weights = weights_dbh
)

nls_bontemps_4p <- nlsLM(
  age_yr ~ 4 + -R * dbh_cm^(1 + m1 * m2) * 
    gauss_2f1(1, m1 + m2^(-1), 1 + m1 + m2^(-1), 
              m1 / ((-1 + m1) * (K / dbh_cm)^m2)) / 
    (K^(m1 * m2) * (-1 + m1) * (1 + m1 * m2)),
  start = c(K = 2, R = 6, m1 = 0.5, m2 = 2),
  lower = c(0, 0, 0, 0),
  upper = c(Inf, Inf, 1, Inf),
  data = age_dbh,
  weights = weights_dbh
)

# Compare model fits
compute_gof(nls_power, age_dbh$age_yr)
compute_gof(nls_michaelis_menten, age_dbh$age_yr)
compute_gof(nls_mitscherlich, age_dbh$age_yr)
compute_gof(nls_chapman_richards, age_dbh$age_yr)
compute_gof(nls_gompertz, age_dbh$age_yr)
compute_gof(nls_weibull, age_dbh$age_yr)
compute_gof(nls_bontemps_3p, age_dbh$age_yr)
compute_gof(nls_bontemps_4p, age_dbh$age_yr)

# Best model summary
summary(nls_weibull)

# Generate predictions
pred_dbh <- data.table(dbh_cm = seq(0, 36.5, length = 200))
pred_dbh[, power := predict(nls_power, newdata = pred_dbh)]
pred_dbh[, michaelis_menten := predict(nls_michaelis_menten, newdata = pred_dbh)]
pred_dbh[, mitscherlich := predict(nls_mitscherlich, newdata = pred_dbh)]
pred_dbh[, chapman_richards := predict(nls_chapman_richards, newdata = pred_dbh)]
pred_dbh[, gompertz := predict(nls_gompertz, newdata = pred_dbh)]
pred_dbh[, weibull := predict(nls_weibull, newdata = pred_dbh)]
pred_dbh[, bontemps_3p := predict(nls_bontemps_3p, newdata = pred_dbh)]
pred_dbh[, bontemps_4p := predict(nls_bontemps_4p, newdata = pred_dbh)]

pred_dbh_long <- melt(pred_dbh, id.vars = "dbh_cm", 
                      variable.name = "model", value.name = "age_yr")

# Figure: best model (Weibull)
fig_age_dbh <- ggplot(age_dbh, aes(x = dbh_cm, y = age_yr, col = method)) +
  geom_errorbar(aes(ymin = age_yr - 1.96 * se_age_yr, 
                    ymax = age_yr + 1.96 * se_age_yr),
                width = 0) +
  geom_point() +
  geom_line(data = pred_dbh_long[model == "weibull"], 
            aes(y = age_yr), col = 'black') +
  scale_y_continuous(breaks = seq(0, 80, by = 20)) +
  labs(x = "DBH (cm)", y = "Estimated tree age (years)", colour = "Method") +
  theme_bw()

fig_age_dbh

# Figure: all models comparison
fig_age_dbh_all <- ggplot(age_dbh, aes(x = dbh_cm, y = age_yr)) +
  geom_errorbar(aes(ymin = age_yr - 1.96 * se_age_yr, 
                    ymax = age_yr + 1.96 * se_age_yr),
                width = 0) +
  geom_point(aes(shape = method)) +
  geom_line(data = pred_dbh_long, aes(y = age_yr, colour = model)) +
  scale_y_continuous(breaks = seq(0, 80, by = 20)) +
  facet_wrap(~model) +
  labs(x = "DBH (cm)", y = "Estimated tree age (years)", 
       shape = "Method", colour = "Model") +
  theme_bw()

fig_age_dbh_all

################################################################################
# Age versus height analysis

age_height <- read.csv("data_age_height_allometry.csv")

# Fit nonlinear models with inverse variance weights
weights_height <- 1 / age_height$se_age_yr^2

# Power function
nls_power <- nlsLM(
  age_yr ~ a * height_m^b,
  start = c(a = 0.5, b = 0.7),
  data = age_height,
  weights = weights_height
)

# Saturating functions
nls_michaelis_menten <- nlsLM(
  age_yr ~ a * height_m / (b + height_m),
  start = c(a = 100, b = 10),
  data = age_height,
  weights = weights_height
)

nls_mitscherlich <- nlsLM(
  age_yr ~ a * (1 - exp(-log(2) / b * height_m)),
  start = c(a = 100, b = 10),
  data = age_height,
  weights = weights_height
)

# Sigmoid saturating functions
nls_chapman_richards <- nlsLM(
  age_yr ~ a * (1 - exp(-log(2) / b * height_m))^pow,
  start = c(a = 100, b = 10, pow = 1),
  data = age_height,
  weights = weights_height
)

nls_gompertz <- nlsLM(
  age_yr ~ a * exp(-b * exp(-d * height_m)),
  start = c(a = 10, b = 1, d = 1),
  data = age_height,
  weights = weights_height
)

nls_weibull <- nlsLM(
  age_yr ~ a * (1 - exp(-b * height_m^pow)),
  start = c(a = 100, b = 1, pow = 1),
  data = age_height,
  weights = weights_height
)

# Bontemps-Duplat models
nls_bontemps_3p <- nlsLM(
  age_yr ~ (exp(K) * exp(R) * log(1 + (-1 + exp(m2)) * (height_m / exp(K))^exp(m2))) / 
    (exp(m2) - 1),
  start = c(K = log(4), R = log(5), m2 = log(2)),
  data = age_height,
  weights = weights_height
)

nls_bontemps_4p <- nlsLM(
  age_yr ~ -(R * height_m^(1 + exp(m1) * m2) * 
    gauss_2f1(1, exp(m1) + m2^(-1), 1 + exp(m1) + m2^(-1),
              exp(m1) / ((-1 + exp(m1)) * (K / height_m)^m2))) / 
    (K^(exp(m1) * m2) * (-1 + exp(m1)) * (1 + exp(m1) * m2)),
  start = c(K = 2, R = 6, m1 = -1, m2 = 2),
  lower = c(0, 0, -10, 0),
  upper = c(Inf, Inf, 0, 20),
  data = age_height,
  weights = weights_height
)

# Compare model fits
compute_gof(nls_power, age_height$age_yr)
compute_gof(nls_michaelis_menten, age_height$age_yr)
compute_gof(nls_mitscherlich, age_height$age_yr)
compute_gof(nls_chapman_richards, age_height$age_yr)
compute_gof(nls_gompertz, age_height$age_yr)
compute_gof(nls_weibull, age_height$age_yr)
compute_gof(nls_bontemps_3p, age_height$age_yr)
compute_gof(nls_bontemps_4p, age_height$age_yr)

# Best model summary
summary(nls_weibull)

# Generate predictions
pred_height <- data.table(height_m = seq(0, 18.6, length = 200))
pred_height[, power := predict(nls_power, newdata = pred_height)]
pred_height[, michaelis_menten := predict(nls_michaelis_menten, newdata = pred_height)]
pred_height[, mitscherlich := predict(nls_mitscherlich, newdata = pred_height)]
pred_height[, chapman_richards := predict(nls_chapman_richards, newdata = pred_height)]
pred_height[, gompertz := predict(nls_gompertz, newdata = pred_height)]
pred_height[, weibull := predict(nls_weibull, newdata = pred_height)]
pred_height[, bontemps_3p := predict(nls_bontemps_3p, newdata = pred_height)]
pred_height[, bontemps_4p := predict(nls_bontemps_4p, newdata = pred_height)]

pred_height_long <- melt(pred_height, id.vars = "height_m",
                         variable.name = "model", value.name = "age_yr")

# Figure: best model (Weibull)
fig_age_height <- ggplot(age_height, aes(x = height_m, y = age_yr, col = method)) +
  geom_errorbar(aes(ymin = age_yr - 1.96 * se_age_yr,
                    ymax = age_yr + 1.96 * se_age_yr),
                width = 0) +
  geom_point() +
  geom_line(data = pred_height_long[model == "weibull"],
            aes(y = age_yr), colour = "black") +
  scale_y_continuous(breaks = seq(0, 80, by = 20)) +
  labs(x = "Height (m)", y = "Estimated tree age (years)", colour = "Method") +
  theme_bw()

fig_age_height

# Figure: all models comparison
fig_age_height_all <- ggplot(age_height, aes(x = height_m, y = age_yr)) +
  geom_errorbar(aes(ymin = age_yr - 1.96 * se_age_yr,
                    ymax = age_yr + 1.96 * se_age_yr),
                width = 0) +
  geom_point(aes(shape = method)) +
  geom_line(data = pred_height_long, aes(y = age_yr, colour = model)) +
  scale_y_continuous(breaks = seq(0, 80, by = 20)) +
  facet_wrap(~model) +
  labs(x = "Height (m)", y = "Estimated tree age (years)",
       shape = "Method", colour = "Model") +
  theme_bw()

fig_age_height_all