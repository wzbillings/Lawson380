###
# Testing Dimensionless Generation
# @author Zane Billings
# @date 2019-11-24
# Contains all the messing around for dimensionaless exploration so that script
#  can be sourced correctly.
###

par(mfrow = c(2,2))
ex1 <- generate_dimensionless_logistic_data(.25, 0.1, 50, 0.1, TRUE)
ex2 <- generate_dimensionless_logistic_data(.25, 2, 50, 0.1, TRUE)
ex3 <- generate_dimensionless_logistic_data(.50, 0.1, 50, 0.1, TRUE)
ex4 <- generate_dimensionless_logistic_data(.50, 2, 50, 0.1, TRUE)
par(mfrow = c(1,1))
ex5 <- generate_dimensionless_logistic_data(0.01, 2, 50, 0.1, TRUE)

par(mfrow = c(2,1))
prep_ex1 <- prep_data(ex1)
fit_ex1 <- model_logistic_data(prep_ex1)

prep_ex2 <- prep_data(ex2)
fit_ex2 <- model_logistic_data(prep_ex2)

prep_ex3 <- prep_data(ex3)
fit_ex3 <- model_logistic_data(prep_ex3)

prep_ex4 <- prep_data(ex4)
fit_ex4 <- model_logistic_data(prep_ex4)

prep_ex5 <- prep_data(ex5)
fix_ex5 <- model_logistic_data(prep_ex5, TRUE)


# Test w/ and w/o noise.
set.seed(100)
par(mfrow = c(2,2))
no_noise <- generate_dimensionless_logistic_data(0.1, .2, 25, 0.1, TRUE)
prep_nn <- prep_data(no_noise)
fit_nn <- model_logistic_data(prep_nn, TRUE)

noisy_4p <- generate_noisy_dimensionless_logistic_data(0.1, .2, 0.004, 25, 0.1, TRUE)
prep_wn <- prep_data(noisy_4p)
fit_wn <- model_logistic_data(prep_wn, TRUE)

