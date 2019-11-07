###
# Zane Billings
# Logistic data fitting
# Created 2019-10-06
#  Fitting a least-squares model to logistic time-series data in several 
#  different ways to see if we can get known parameters back.
###

require(MASS) # ginv() function source for generalize inverse computation.

prep_data <- function(output_data) {
  # This function takes in a list of the time points and population size values,
  #  and puts it into the form used by the regression methods.
  # This saves a little time on calculations when trying multiple methods.

  # Extract the two vectors from the dataframe to make life easier.
  t <- output_data$t
  P <- output_data$P
  
  # Calculate the maximum time value and the time step from the list of times,
  #  assuming evenly spaced time points.
  max_t <- max(t)
  time_step <- t[[2]] - t[[1]]
  
  # We need to use P_n and P_{n+1} for the least squares method, so calculate
  #  those here.
  P_forward <- P[2:max_t]
  P_present <- P[1:(max_t - 1)]
  
  processed_data <- list(
    "forward" <- P_forward,
    "present" <- P_present,
    "step" <- time_step
  )
  
  # return(c(P_forward,P_present,time_step))
  return(processed_data)
}

model_logistic_data <- 
  function(input_data, make_plot = FALSE){
    # This function uses a derived least-squares method with an Euler
    #  discretization in order to solve for the parameters of the logistic
    #  equation modeling the inputted data.
    # The input should be a data frame returned from one of the generators.
    
    # Configure input correctly
    #list[P_forward, P_present, time_step] <- prep_data(input_data)
    P_forward <- input_data[[1]]
    P_present <- input_data[[2]]
    time_step <- input_data[[3]]
    
    # The least squares algorithm for this model can be shown to estimate r and
    #   m = r/K when b = 1/t((P+1)/P) and A = (1 - P). The solution will be
    #   of the form x = [a; m].
    # These three lines calculate the regression model.
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

model_beginning_only <- function(input_data, make_plot = FALSE) {
  # Estimate time series parameters from data before the inflection point of
  #  the logistic curve.
  
  # How to find the inflection point? 
  # Source: https://socratic.org/questions/how-do-you-find-the-inflection-point-of-a-logistic-function
  # The inflection point is at (ln(A)/K, K/2), so the inflection point is at
  #  the spot where the population is half the carrying capacity.
  # BUT we are estimating the carrying capacity!
  
  # Configure input correctly
  #list[P_forward, P_present, time_step] <- prep_data(input_data)
  P_forward <- input_data[[1]]
  P_present <- input_data[[2]]
  time_step <- input_data[[3]]
  
  # Assume the model goes to K eventually.
  inflection_point <- max(P_forward) / 2
  
  # Only take values less than the inflection point.
  P_forward <- P_forward[P_forward <= inflection_point]
  P_present <- P_present[1:length(P_forward)]
  
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

model_end_only <- function(input_data, make_plot = FALSE) {
  # Estimate time series parameters from data after the inflection point of
  #  the logistic curve.
  
  # Configure input correctly
  #list[P_forward, P_present, time_step] <- prep_data(input_data)
  P_forward <- input_data[[1]]
  P_present <- input_data[[2]]
  time_step <- input_data[[3]]
  
  # Assume the model goes to K eventually.
  inflection_point <- max(P_forward) / 2
  
  # Only take values less than the inflection point.
  P_forward <- P_forward[P_forward >= inflection_point]
  P_present <- P_present[1:length(P_forward)]
  
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

# model_random_sample <- function(input_data) {
#   # Estimate time series parameters from a random subset of half the 
#   #  data points.
#   num_samples <- floor(0.5 * nrow(input_data))
#   sample_data <- input_data[sample(nrow(input_data), num_samples), ]
#   model <- model_logistic_data(sample_data)
#   return(model)
# }

calculate_SSR <- function(model) {
  # Calculate the sum of squared residuals for a model.
  model_residuals <- model$residuals
  SSR <- sum(model_residuals ^ 2)
  return(SSR)
}