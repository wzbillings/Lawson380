source('~/Stuff/Research/Lawson_380/Lawson_380_R_code/Sample_data_generation.R')
source('~/Stuff/Research/Lawson_380/Lawson_380_R_code/Least_squares_methods.R')
set.seed(300)
test <- generate_noisy_analytic_logistic_data(100,1000,0.2,25,0.1,0.00004, TRUE)
this <- prep_data(test)
model <- model_logistic_data(this)

test_analytic <- generate_analytic_logistic_data(100,1000,0.2,25,0.1, TRUE)
prepped_analytic <- prep_data(test_analytic)
model_2 <- model_logistic_data(prepped_analytic)

table <- data.frame(test$t, test_analytic$P, test$P)
colnames(table) <- c("times", "regular", "noisy")

test3 <- generate_noisy_analytic_logistic_data(100,1000,1,25,0.1,0.004, TRUE)
this3 <- prep_data(test3)
model3 <- model_logistic_data(this3)

# How much noise can you add before getting negative issues?
#  Test: 0.04 and increase noise parameter to about 0.16.

# If noise parameter is 0, this works.
# If noise parameter is 0.00004, works perfectly. BUT not for 0.0004!