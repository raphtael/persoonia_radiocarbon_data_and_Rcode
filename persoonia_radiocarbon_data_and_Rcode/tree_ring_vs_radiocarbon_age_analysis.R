################################################################################
# Tree ring versus radiocarbon age figure
# Bayesian log-log linear regression with error in variable
# Input: data_tree_ring_vs_radiocarbon.csv
################################################################################

# Load libraries
library(ggplot2)
library(brms)

################################################################################
# Load data
data_tree_ring_vs_radiocarbon <- read.csv('data_tree_ring_vs_radiocarbon.csv')

# Set priors for Bayesian measurement error model
priors <- c( 
  set_prior("normal(0, 0.5)", class = "Intercept"),
  set_prior("normal(1, 0.5)", class = "b")
)   

# Fit model with measurement error in radiocarbon age
model <- brm(
  log(tree_ring_age_yr) ~ me(log(f14_age_yr), se_log_f14_age_yr),
  family = gaussian,
  data = data_tree_ring_vs_radiocarbon,
  chains = 2,
  iter = 2000,
  prior = priors,
  backend = 'cmdstanr'
) 
summary(model)

# Generate predictions
newdata <- expand.grid(
  f14_age_yr = 28:76, 
  se_log_f14_age_yr = mean(data_tree_ring_vs_radiocarbon$se_log_f14_age_yr) # 0.048
) 

fits <- fitted(model, newdata = newdata)
newdata$log_age_mean <- fits[, 'Estimate']
newdata$log_age_se <- fits[, 'Est.Error']

# Back-transform to original scale
newdata$tree_ring_age_yr <- exp(newdata$log_age_mean + newdata$log_age_se^2 / 2)
newdata$ymin <- exp(newdata$log_age_mean - 2 * newdata$log_age_se)
newdata$ymax <- exp(newdata$log_age_mean + 2 * newdata$log_age_se)

# Create figure
ggplot(data = data_tree_ring_vs_radiocarbon, aes(f14_age_yr, tree_ring_age_yr)) + 
  geom_point() +
  geom_abline(lty = 2) + 
  geom_line(data = newdata, aes(f14_age_yr, tree_ring_age_yr)) + 
  geom_ribbon(data = newdata, aes(x = f14_age_yr, ymin = ymin, ymax = ymax), 
              alpha = 0.3) +
  geom_segment(aes(x = min_f14_age_yr, xend = max_f14_age_yr, 
                   y = tree_ring_age_yr, yend = tree_ring_age_yr)) +
  scale_x_log10() + 
  scale_y_log10() + 
  coord_cartesian(ylim = c(24, 80), xlim = c(24, 80)) +
  labs(x = 'Estimated tree age using F14C (years)',
       y = 'Estimated tree ring age (years)') +
  theme_bw()
  
  
  
  
              