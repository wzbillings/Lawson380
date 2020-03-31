###
# Helper Functions For Analysis
# @author Zane Billings
# @date 2019-11-24
# Provides helper functions so they aren't in multiple places at once.
###

prep_data <- function(output_data, time_step = NULL) {
  # This function takes in a list of the time points and population size values,
  #  and puts it into the form used by the regression methods.
  # This saves a little time on calculations when trying multiple methods.
  
  # Extract the two vectors from the dataframe to make life easier.
  t <- output_data$t
  P <- output_data$P
  
  # Calculate the maximum time value and the time step from the list of times,
  #  assuming evenly spaced time points.
  max_t <- max(t)
  if (is.null(time_step)) {
    time_step <- t[[2]] - t[[1]]
  }
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

validate_inputs <- function(...) {
  check_these <- list(...)
  for (param in check_these) {
    if (param <= 0 || mode(param) != "numeric") {
      stop("Invalid parameter. \nAll parameters should be positive real numbers.")
    }
  }
}

random_partition <- function(data, k) {
  # data is a data frame of observations
  # k is the number of partitions to form
  # Returns a list of k random partitions of the data
  num_obs <- nrow(data)
  part_size <- floor(num_obs/k)
  rem <- seq(1, num_obs, 1)
  parts <- list()
  for (i in 1:(k - 1)) { 
    choices <- sample(rem, part_size, replace = F)
    parts[[i]] <- data[choices, ]
    rem <- setdiff(rem, choices)
  }
  parts[[k]] <- data[rem, ]
  return(parts)
}

calculate_SSE <- function(true, fits) {
  # Calculate the sum of squared errors for a model.
  SSE <- sum((true - fits)^2)
  return(SSE)
}