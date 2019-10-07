###
# Zane Billings
# Sample data generation 
# Created 2019-09-13
#  This script will generate time series data using several different methods.
#  The data replicates 1D logistic population growth.
###

require(deSolve) # Provides methods for solving ODE problems.

solve_logistic_equation <- 
  function(P0, K, r, t) {
  # This function solves the logistic equation with the given parameters at
  #  returns the value of P at time t, using the known analytic solution.
  
  # The analytic solution to the logistic growth equation is P = K/(1+Ae^(-rt)),
  #  where A = (K - P0)/P0, with P0 the initial population size.
  A = (K - P0)/P0
  P = K/(1 + A*exp(-r*t))
  return(P)
  }

calc_logistic_derivative <-
  function(P, t, params) {
    # This function calculates the derivative of the logistic different equation
    # with supplied parameters and value P at time t.
    with(as.list(c(params)), { 
    
    dn <- r*P*(1-P/K)
    return(list(dn))
    
    })
  }

# The logistic growth equation is: dP/dt = rP(1-P/k), where P is the population
#  size, t is the time, r is the population intrinsic growth rate, and k is the
#  carrying capacity. The initial population size, P0 is used in the
#  analytic solution to the equation.

generate_analytic_logistic_data <- 
  function(P_naught, K, r, max_t, time_step, make_plot = FALSE) {
  # This function generates sample time-series data with given parameters
  #  by using the analytic solution to the logistic equation.
  
  # Generate a list of time steps to solve for.
  t <- seq(from = 0, to = max_t, by = time_step)
  # Use the analytic solution to solve for P at each time.
  P <- solve_logistic_equation(P_naught, K, r, t)
  
  # Optionally generate a plot of the logistic curve.
  if (make_plot) {
    plot(x = t, y = P, ylab = "P(t)", type = "b")
  }
  # Output a list of results to be returned.
  output <- data.frame(t,P)
  return(output)
}

generate_noisy_analytic_logistic_data <- 
  function(P_naught, K, r, max_t, time_step, noise, make_plot = FALSE) {
  # This function generates sample logistic data using the analytic solution
  #  but with added Gaussian noise with SD equal to noise% of K.
  
  # Generate the time steps to be used.
  t <- seq(from = 0, to = max_t, by = time_step)
  
  # Calculate the analytic solution at each time step.
  P <- solve_logistic_equation(P_naught, K, r, t)
  P[0] <- P_naught
  
  # Generate noise
  noise_mean <- 0
  noise_sd <- noise * K
  noise_vector <- rnorm(length(P), noise_mean, noise_sd)
  
  # Add noise to P values.
  P <- P + noise_vector
  
  # Optionally, plot the noisy logistic data.
  if (make_plot == TRUE) {
    plot(x = t, y = P, ylab = "P(t)")
  }
  
  # Make the two data series to be returned into a list.
  output <- data.frame(t,P)
  return(output)
}

generate_euler_logistic_data <- 
  function(P_naught, K, r, max_t, time_step, make_plot = FALSE) {
  # This function uses the Euler discretization of the logistic growth DE in 
  #  order to generate time series data.

  # Create a vector of times to use
  t <- seq(from = 0, to = max_t, by = time_step)
  # Intialize P_naught
  P <- numeric(length(t))
  P[1] <- P_naught
  
  # Calculate all of the values of P-naught
  for (j in 2:length(t)) {
    P[j] <- P[j - 1] + r*P[j - 1] * (1 - P[j - 1]/K) * time_step
  }
  
  # Optionally, plot the data.
  if (make_plot == TRUE) {
    plot(x = t, y = P, ylab = "P(t)")
  }
  
  # Create output dataframe
  output <- data.frame(t,P)
  return(output)
}

generate_RK4_logistic_data <- 
  function(P_naught, K, r, max_t, time_step, make_plot = FALSE) {
  # This function uses an ODE solver method in order to generate logistic time
  #  series data rather than the analytic solution.
  solver_output <- deSolve::ode(y = c("P" = P_naught), 
                    times = seq(from = 0, to = max_t, by = time_step), 
                    func = calc_logistic_derivative, 
                    parms = c("K" = K, "r" = r))
  output <- as.data.frame(solver_output)
  
  # Optionally, generate a plot.
  if (make_plot == TRUE) {
    plot(output, ylab = "P(t)")
  }
  
  # Return output data frame
  return(output)
}


# Explore sensitivities
# use ODE solver to generate data instead! Check if this breaks everything!
# Use euler method to generate data--if we still get error back from this, we
#  can tell how much the least squares machinery is breaking.

# What techniques work best with no noise?
# Testing data series before or after inflection point -- where does the error
#  show up at?
# What if we sample every other data point? How does this affect error?
#  Residuals of this curve with the better curve--residuals of missing points?
#  Metric: sum of squared residuals.
# Sum of residuals is net error squared, this is our objective fcn.
# Then we can sample other subsets, e.g. just the front, just the back, and
#  random sampling. Use multiple fits to fit data and one fit to test.
# Weighting least squares with smoothing--kind of like Lagrange multipliers.

# Need to add input testing to each function: 
# Probably want to implement this with tryCatch() or stopifnot().
# #Check if any parameters are not strictly positive.
# not_zero_condition <- any(P_naught <= 0, K <= 0, r <= 0, max_t <= 0)
# If any params are not strictly positive or not numeric, stop.
# if (not_zero_condition == TRUE || is.na(not_zero_condition)) {
#   stop("Invalid parameter. \nAll parameters should be strictly positive reals.")
# }

# Ideally all of the data generation functions should be combined into one 
#  function with a "method" argument".

# Re-test noisy logistic fit to get back parameter.
# Try Euler method with deSolve!
# Compare Euler w/ Analytic
# Figure out ODE function with something easier.
