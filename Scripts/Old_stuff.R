###
# Old Stuff for Lawson 380 Project
# @author Zane Billings
# @date 2019-11-24
# Old functions that do not serve a purpose in the current "pipeline".
###

model_every_other_odd <- function(input_data, make_plot = FALSE) {
  # Estimate time series parameters from data points 1, 3, 5, ... last odd.
  
  # Configure input correctly
  #list[P_forward, P_present, time_step] <- prep_data(input_data)
  P_forward <- input_data[[1]]
  P_present <- input_data[[2]]
  time_step <- input_data[[3]]
  
  # Only take odd indices.
  odd_indices <- seq(from = 1, to = length(P_forward), by = 2)
  P_forward <- P_forward[odd_indices]
  P_present <- P_present[odd_indices]
  
  # Run the OLS model.
  b <- (1/time_step)*log(P_forward/P_present)
  A <- cbind(rep(1, length(P_present)), -P_present)
  params <- MASS::ginv(t(A) %*% A) %*% t(A) %*% b
  
  # Now knowing that m = r/K, we solve for K.
  r_hat <- params[1]
  m_hat <- params[2]
  K_hat <- (r_hat/m_hat)
  cat("The estimated growth rate is:", r_hat,
      "\nThe estimated carrying capacity is:", K_hat, "\n")
  
  # Optionally make a linear plot showing the data the regression equation
  #  is obtained from.
  if (make_plot == TRUE) {
    plot(P_present, b)
  }
  
  # Collect outputs into a vector and return.
  output_params <- as.numeric(r_hat, K_hat)
  return(output_params)
  
}

model_every_other_even <- function(input_data, make_plot = FALSE) {
  # Estimate time series parameters from data points 2, 4, 6, ... last even.
  
  # Configure input correctly
  #list[P_forward, P_present, time_step] <- prep_data(input_data)
  P_forward <- input_data[[1]]
  P_present <- input_data[[2]]
  time_step <- input_data[[3]]
  
  # Only take even indices.
  even_indices <- seq(from = 2, to = length(P_forward), by = 2)
  P_forward <- P_forward[even_indices]
  P_present <- P_present[even_indices]
  
  # Run the OLS model.
  b <- (1/time_step)*log(P_forward/P_present)
  A <- cbind(rep(1, length(P_present)), -P_present)
  params <- MASS::ginv(t(A) %*% A) %*% t(A) %*% b
  
  # Now knowing that m = r/K, we solve for K.
  r_hat <- params[1]
  m_hat <- params[2]
  K_hat <- (r_hat/m_hat)
  cat("The estimated growth rate is:", r_hat,
      "\nThe estimated carrying capacity is:", K_hat, "\n")
  
  # Optionally make a linear plot showing the data the regression equation
  #  is obtained from.
  if (make_plot == TRUE) {
    plot(P_present, b)
  }
  
  # Collect outputs into a vector and return.
  output_params <- as.numeric(r_hat, K_hat)
  return(output_params)
}

calculate_SSR <- function(model) {
  # Calculate the sum of squared residuals for a model.
  model_residuals <- model$residuals
  SSR <- sum(model_residuals ^ 2)
  return(SSR)
}

model_random_sample <- function(input_data) {
  # Estimate time series parameters from a random subset of half the
  #  data points.
  num_samples <- floor(0.5 * nrow(input_data))
  sample_data <- input_data[sample(nrow(input_data), num_samples), ]
  model <- model_logistic_data(sample_data)
  return(model)
}