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
  
  return(processed_data)
}

model_logistic_data <- 
  function(list, make_plot = FALSE){
    # This function uses a derived least-squares method with an Euler
    #  discretization in order to solve for the parameters of the logistic
    #  equation modeling the inputted data.
    # The input should be a list where the first entry is the time values and the
    #  second entry is the population size values at each time.
    
    P_forward <- list[[1]]
    P_present <- list[[2]]
    time_step <- list[[3]]
    
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
    K_hat <- (1/m_hat)*r_hat
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

model_beginning_only <- function() {
  # Estimate time series parameters from data before the inflection point of
  #  the logistic curve.
}

model_end_only <- function() {
  # Estimate time series parameters from data after the inflection point of
  #  the logistic curve.
}

model_every_other_odd <- function() {
  # Estimate time series parameters from data points 1, 3, 5, ... last odd.
}

model_every_other_even <- function() {
  # Estimate time series parameters from data points 2, 4, 6, ... last even.
}

model_random_sample <- function() {
  # Estimate time series parameters from a random subset of data points.
}

calculate_SSR <- function() {
  # Calculate the sum of squared residuals for a model.
}