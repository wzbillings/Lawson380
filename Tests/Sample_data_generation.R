###
# Sample data generation function
# Author: Zane Billings
# Date created: 2019-09-13
# Last modified: 2019-10-19
# This script will provides several different methods for simulating 1D
#  logistic population growth data.
###



# The logistic growth equation is: dP/dt = rP(1-P/k), where P is the population
#  size, t is the time, r is the population intrinsic growth rate, and k is the
#  carrying capacity. The initial population size, P0 is used in the
#  analytic solution to the equation.

generate_analytic_logistic_data <-
  function(P_naught, K, r, max_t, time_step, make_plot = FALSE) {
  # This function generates sample time-series data with given parameters
  #  by using the analytic solution to the logistic equation.

  # Make sure that all inputs are strictly positive integers.
  validate_inputs(P_naught, K, r, max_t, time_step)


  # Use the analytic solution to solve for P at each time.


  # Optionally generate a plot of the logistic curve.

  # Output a dataframe of results to be returned.

}

generate_noisy_analytic_logistic_data <-
  function(P_naught, K, r, max_t, time_step, noise, make_plot = FALSE) {
  # This function generates sample logistic data using the analytic solution
  #  but with added Gaussian noise with SD equal to noise% of K.

  # Make sure that all inputs are strictly positive integers.
  validate_inputs(P_naught, K, r, max_t, time_step)

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

  # Make the two data series to be returned into a data frame.
  output <- data.frame(t,P)
  return(output)
}

generate_euler_logistic_data <-
  function(P_naught, K, r, max_t, time_step, make_plot = FALSE) {
  # This function uses the Euler discretization of the logistic growth DE in
  #  order to generate time series data.

  # Make sure that all inputs are strictly positive integers.
  validate_inputs(P_naught, K, r, max_t, time_step)

  # Create a vector of times to use.
  t <- seq(from = 0, to = max_t, by = time_step)
  # Intialize P_naught.


  # Calculate all of the values of P.


  # Optionally, plot the data.
  if (make_plot == TRUE) {
    plot(x = t, y = P, ylab = "P(t)")
  }

  # Create output dataframe
  output <- data.frame(t,P)
  return(output)
}

generate_noisy_euler_logistic_data <-
  function(P_naught, K, r, max_t, time_step, noise, make_plot = FALSE) {
  # This function uses the Euler discretization of the logistic growth DE in
  #  order to generate time series data. Noise is added to the data.

  # Make sure that all inputs are strictly positive integers.
  validate_inputs(P_naught, K, r, max_t, time_step)

  # Create a vector of times to use.
  t <- seq(from = 0, to = max_t, by = time_step)
  # Intialize P_naught.
  P <- numeric(length(t))
  P[1] <- P_naught

  # Calculate all of the values of P.
  for (j in 2:length(t)) {
    P[j] <- P[j - 1] + r*P[j - 1] * (1 - P[j - 1]/K) * time_step
  }

  # Generate noise
  noise_mean <- 0
  noise_sd <- noise * K
  noise_vector <- rnorm(length(P), noise_mean, noise_sd)

  # Add noise to P values.
  P <- P + noise_vector

  # Optionally, plot the data.
  if (make_plot == TRUE) {
    plot(x = t, y = P, ylab = "P(t)")
  }

  # Create output dataframe
  output <- data.frame(t,P)
  return(output)

}

generate_solver_logistic_data <-
  function(P_naught, K, r, max_t, time_step, make_plot = FALSE) {
  # This function uses an ODE solver method in order to generate logistic time
  #  series data rather than the analytic solution.

  # Make sure that all inputs are strictly positive integers.
  validate_inputs(P_naught, K, r, max_t, time_step)

  # Use the default ode solver method to generate time series data.


  # The output is a "deSolve" object which can be coerced into a data frame.
  output <- as.data.frame(solver_output)

  # Optionally, generate a plot.
  if (make_plot == TRUE) {
    plot(output, ylab = "P(t)")
  }

  # Return output data frame
  return(output)
}

generate_noisy_solver_logistic_data <-
  function(P_naught, K, r, max_t, time_step, noise, make_plot = FALSE) {
    # This function uses an ODE solver method in order to generate logistic time
    #  series data rather than the analytic solution. Adds noise.

    # Make sure that all inputs are strictly positive integers.
    validate_inputs(P_naught, K, r, max_t, time_step)

    # Use the default ode solver method to generate time series data.
    solver_output <- deSolve::ode(
      y = c(P = P_naught),
      times = seq(from = 0, to = max_t, by = time_step),
      func = calculate_logistic_derivative,
      parms = c(r = r, K = K)
    )

    # The output is a "deSolve" object which can be coerced into a data frame.
    output <- as.data.frame(solver_output)

    # Generate noise
    noise_mean <- 0
    noise_sd <- noise * K
    noise_vector <- rnorm(length(output$P), noise_mean, noise_sd)

    # Add noise to P values.
    output$P <- output$P + noise_vector

    # Optionally, generate a plot.
    if (make_plot == TRUE) {
      plot(output, ylab = "P(t)")
    }

    # Return output data frame
    return(output)
}

# The dimensionless population growth model in 1D is given as
# dx/dt = rx(1 - x), where x = P/k. The solution is
# x(t) = (e^(rt)*x0) / (1-x0 + e^(et)*x0)

generate_dimensionless_logistic_data <-
  function(P0, r, max_t, time_step, make_plot = FALSE) {
    # This function generates sample time-series data with given parameters
    #  by using the analytic solution to the logistic equation.

    # Make sure that all inputs are strictly positive integers.
    validate_inputs(P0, r, max_t, time_step)

    # Generate a list of time steps to solve for.
    t <- seq(from = 0, to = max_t, by = time_step)
    # Use the analytic solution to solve for x at each time.
    P <- (exp(r*t)*P0) / (1-P0 + exp(r*t)*P0)

    # Optionally generate a plot of the logistic curve.
    if (make_plot) {
      plot(x = t, y = P, ylab = "x(t)")
    }
    # Output a dataframe of results to be returned.
    output <- data.frame(t,P)
    return(output)
  }


generate_noisy_dimensionless_logistic_data <-
  function(P0, r, noise, max_t, time_step, make_plot = FALSE) {
    # This function generates sample time-series data with given parameters
    #  by using the analytic solution to the logistic equation.

    # Make sure that all inputs are strictly positive integers.


    # Generate a list of time steps to solve for.
    t <- seq(from = 0, to = max_t, by = time_step)
    # Use the analytic solution to solve for x at each time.
    P <- (exp(r*t)*P0) / (1-P0 + exp(r*t)*P0)

    # Add noise
    noise_vec <- rnorm(length(P), mean = 0, sd = noise)
    P <- P + noise_vec

    # Optionally generate a plot of the logistic curve.
    if (make_plot) {
      plot(x = t, y = P, ylab = "x(t)")
    }
    # Output a dataframe of results to be returned.
    output <- data.frame(t,P)
    return(output)
  }

