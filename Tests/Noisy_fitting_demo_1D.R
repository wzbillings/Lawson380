###
# Testing Noisy vs. Non-noisy fitting in 1D for manuscript.
# @author Zane Billings
# @date 2019-12-13

# Source in needed scripts from elsewhere in the Project.
require(here)
source(here::here("Scripts", "Dimensionless_exploration.R"))
source(here::here("Scripts", "Helpers.R"))

# Some changing some settings here
set.seed(300)
par(mfrow = c(3,2))

# Generate data with no noise.
cat("No noise\n")
test_analytic <- generate_dimensionless_logistic_data(.1, 0.1, 50, 0.1, TRUE)
prepped_analytic <- prep_data(test_analytic)
model_2 <- model_logistic_data_dimensionless(prepped_analytic, TRUE)

# Generate data with almost no noise.
cat("\n0.04% noise\n")
test_low_noise <- generate_noisy_dimensionless_logistic_data(0.1, 0.1, 0.0004, 50, 0.1, TRUE)
prepped_low_noise <- prep_data(test_low_noise)
model_low_noise <- model_logistic_data_dimensionless(prepped_low_noise, TRUE)

# Generate data with noise.
cat("\n4% noise\n")
test_high_noise <- generate_noisy_dimensionless_logistic_data(0.1, 0.1, 0.04, 50, 0.1, TRUE)
prepped_high_noise <- prep_data(test_high_noise)
model_high_noise <- model_logistic_data_dimensionless(prepped_high_noise, TRUE)
