################################################################################
# Age required to reach 1.3 m height
# Logistic regression analysis
# Input: data_regen.csv
################################################################################

# Load library
library(ggplot2)

# Utility function
inv_logit <- function(x) 1 / (1 + exp(-x))

################################################################################
# Load data
regen <- read.csv("data_age_to_1.3m_height.csv")

# Fit logistic regression
glm_height <- glm(persoonia_regen_height_m > 1.3 ~ max_tree_age_yrs, 
                  family = binomial, 
                  data = regen)
summary(glm_height)

# Calculate age at which 50% of trees reach 1.3 m
# At 50% probability, logit = 0, so age = -intercept / slope
age_at_50pct <- -coef(glm_height)[1] / coef(glm_height)[2]
round(age_at_50pct, 1)  # 4.7 years

# Accounting for 1 year lag between harvest and regeneration
round(age_at_50pct - 1, 1)  # 3.7 years

# Generate predictions
pred_height <- data.frame(max_tree_age_yrs = seq(2, 6, length = 50))

fits <- predict(glm_height, newdata = pred_height, se.fit = TRUE)
pred_height$mean <- inv_logit(fits$fit)
pred_height$ymin <- inv_logit(fits$fit - 1.96 * fits$se.fit)
pred_height$ymax <- inv_logit(fits$fit + 1.96 * fits$se.fit)

# Create figure
fig_height_age <- ggplot(pred_height, aes(x = max_tree_age_yrs)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.3) +
  geom_line(aes(y = mean)) +
  labs(x = "Maximum tree age (years)",
       y = "Probability of reaching 1.3 m height") +
  theme_bw()
fig_height_age