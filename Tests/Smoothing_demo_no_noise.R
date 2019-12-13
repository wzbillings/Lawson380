###
# Plots for Demonstrating Smoothing in Manuscript without Noise
# @author Zane Billings
# @data 2019-12-13
# Note that this code is for the most part copied straight from the file
# 'Tests/Smoothing_test.Rmd'.
###

require(here)
source(here::here("Scripts", "Least_squares_methods.R"))
source(here::here("Scripts", "Sample_data_generation.R"))
source(here::here("Scripts", "Helpers.R"))
set.seed(300)

# Generate data without noise.
test_2p_nn <- generate_analytic_logistic_data(100, 1000, 0.2, 25, 0.1, TRUE)
prep_2p_nn <- prep_data(test_2p_nn)

# Test fits with different amounts of smoothing.
cat("Regular; no noise.\n")
plain_2p_nn <- model_logistic_data(prep_2p_nn)
cat("Smoothing = 0.001; no noise\n")
ridge1_2p_nn <- model_logistic_data_smoothing(prep_2p_nn, 0.001)
cat("Smoothing = 0.1; no noise\n")
ridge2_2p_nn <- model_logistic_data_smoothing(prep_2p_nn, 0.1)
cat("Smoothing = 10; no noise\n")
ridge2_2p_nn <- model_logistic_data_smoothing(prep_2p_nn, 10)
