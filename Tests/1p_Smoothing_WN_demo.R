###
# Smoothing Demo 1D, With Noise for Manuscript
# @author Zane Billings
# @date 2019-12-13
###

require(here)
source(here::here("Scripts", "Dimensionless_exploration.R"))
source(here::here("Scripts", "Helpers.R"))
set.seed(300)

test_1p_wn <- generate_noisy_dimensionless_logistic_data(0.1, 0.1, 0.05, 50, 0.1, TRUE)

prep_1p_wn <- prep_data(test_1p_wn)
cat("Regular; no noise.\n")
plain_1p_wn <- model_logistic_data_dimensionless(prep_1p_wn)
cat("Smoothing = 0.001; 4% noise\n")
ridge_1p_wn <- model_logistic_data_dimensionless_smoothing(prep_1p_wn, 0.001, print_res = TRUE)
cat("Smoothing = 0.01; 4% noise\n")
ridge_1p_wn <- model_logistic_data_dimensionless_smoothing(prep_1p_wn, 0.01, print_res = TRUE)
cat("Smoothing = 0.1; 4% noise\n")
ridge_1p_wn <- model_logistic_data_dimensionless_smoothing(prep_1p_wn, 0.1, print_res = TRUE)
cat("Smoothing = 1; 4% noise\n")
ridge_1p_wn <- model_logistic_data_dimensionless_smoothing(prep_1p_wn, 1, print_res = TRUE)
cat("Smoothing = 10 ; 4% noise\n")
ridge_1p_wn <- model_logistic_data_dimensionless_smoothing(prep_1p_wn, 10, print_res = TRUE)
cat("Smoothing = 20 ; 4% noise\n")
ridge_1p_wn <- model_logistic_data_dimensionless_smoothing(prep_1p_wn, 20, print_res = TRUE)
cat("Smoothing = 100 ; 4% noise\n")
ridge_1p_wn <- model_logistic_data_dimensionless_smoothing(prep_1p_wn, 100, print_res = TRUE)


