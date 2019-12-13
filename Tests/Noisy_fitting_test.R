###
# Testing Noisy vs. Non-noisy fitting.
# @author Zane Billings
# @date 2019-10-21; 2019-12-13

# Source in needed scripts from elsewhere in the Project.
require(here)
source(here::here("Scripts","Sample_data_generation.R"))
source(here::here("Scripts", "Least_squares_methods.R"))
source(here::here("Scripts", "Helpers.R"))

# Some changing some settings here
set.seed(300)
par(mfrow = c(3,2))

# Generate data with no noise.
cat("No noise\n")
test_analytic <- generate_analytic_logistic_data(100,1000,0.2,25,0.1, TRUE)
prepped_analytic <- prep_data(test_analytic)
model_2 <- model_logistic_data(prepped_analytic, TRUE)

# Generate data with almost no noise.
cat("\n0.04% noise\n")
test_low_noise <- generate_noisy_analytic_logistic_data(100,1000,0.2,25,0.1,0.0004, TRUE)
prepped_low_noise <- prep_data(test_low_noise)
model_low_noise <- model_logistic_data(prepped_low_noise, TRUE)

# Generate data with noise.
cat("\n4% noise\n")
test_high_noise <- generate_noisy_analytic_logistic_data(100,1000,1,25,0.1,0.04, TRUE)
prepped_high_noise <- prep_data(test_high_noise)
model_high_noise <- model_logistic_data(prepped_high_noise, TRUE)

# Generate data with ridiculous noise.
# cat("\n10% noise\n")
# test_bigg_noise <- generate_noisy_analytic_logistic_data(100,1000,1,25,0.1,0.1, TRUE)
# prepped_bigg_noise <- prep_data(test_bigg_noise)
# model_bigg_noise <- model_logistic_data(prepped_bigg_noise, TRUE)


# If noise parameter is 0, this works.
# If noise parameter is 0.00004, works perfectly. BUT not for 0.04!
# Does stronger initial growth give better answers?
# If we have a smaller gap between P0 and K, do we get better answers?