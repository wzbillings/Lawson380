###
# Plots for Demonstrating Smoothing in Manuscript WITH Noise
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

# Generate data with 4% gaussian noise
test_2p_wn <- generate_noisy_analytic_logistic_data(100, 1000, 0.2, 25, 0.1, 0.04, TRUE)
prep_2p_wn <- prep_data(test_2p_wn)
cat("Regular; 4% noise.\n")
plain_2p_wn <- model_logistic_data(prep_2p_wn)
cat("Smoothing = 0.001; 4% noise\n")
ridge1_2p_wn <- model_logistic_data_smoothing(prep_2p_wn, 0.001, print_res = TRUE)
cat("Smoothing = 0.1; 4% noise\n")
ridge2_2p_wn <- model_logistic_data_smoothing(prep_2p_wn, 0.1, print_res = TRUE)
cat("Smoothing = 1; 4% noise\n")
ridge2_2p_wn <- model_logistic_data_smoothing(prep_2p_wn, 1, print_res = TRUE)
cat("Smoothing = 10; 4% noise\n")
ridge2_2p_wn <- model_logistic_data_smoothing(prep_2p_wn, 10, print_res = TRUE)
cat("Smoothing = 15; 4% noise\n")
ridge2_2p_wn <- model_logistic_data_smoothing(prep_2p_wn, 15, print_res = TRUE)
