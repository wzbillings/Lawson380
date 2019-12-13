###
# Smoothing Demo 1D, No Noise for Manuscript
# @author Zane Billings
# @date 2019-12-13
###

require(here)
source(here::here("Scripts", "Dimensionless_exploration.R"))
source(here::here("Scripts", "Helpers.R"))
set.seed(300)

test_1p_nn <- generate_dimensionless_logistic_data(.1, 0.1, 50, 0.1, TRUE)
prep_1p_nn <- prep_data(test_1p_nn)
cat("Regular; no noise.\n")
plain_1p_nn <- model_logistic_data_dimensionless(prep_1p_nn)
cat("Smoothing = 0.001; no noise\n")
ridge_1p_nn <- model_logistic_data_dimensionless_smoothing(prep_1p_nn, 0.001)
cat("Smoothing = 0.01; no noise\n")
ridge_1p_nn <- model_logistic_data_dimensionless_smoothing(prep_1p_nn, 0.01)
cat("Smoothing = 0.1; no noise\n")
ridge_1p_nn <- model_logistic_data_dimensionless_smoothing(prep_1p_nn, 0.1)
cat("Smoothing = 1; no noise\n")
ridge_1p_nn <- model_logistic_data_dimensionless_smoothing(prep_1p_nn, 1)
cat("Smoothing = 10 ; no noise\n")
ridge_1p_nn <- model_logistic_data_dimensionless_smoothing(prep_1p_nn, 10)


