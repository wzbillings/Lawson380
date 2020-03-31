###
# Dimensionless Model Exploration
# @author Zane Billings
# @date 2019-11-17, 2019-11-18
# Exploring data-generation and least squares fitting for the dimensionless
# predator-prey model, which is reduced to have only one parameter.
###


### COPIED CODE FIX THIS ###
require(MASS) # ginv() function source for generalize inverse computation.

model_logistic_data_dimensionless <- 
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
    A <- 1 - P_present
    params <- MASS::ginv(t(A) %*% A) %*% t(A) %*% b
    
    r_hat <- params[1]
    cat("The estimated growth rate is:", r_hat,"\n")
    
    # Optionally make a linear plot showing the data the regression equation
    #  is obtained from.
    if (make_plot == TRUE) {
      plot(1 - P_present, b)
    }
    
    # Collect outputs into a vector and return.
    output_params <- as.numeric(r_hat)
    return(output_params)
  }

model_logistic_data_dimensionless_smoothing <- 
  function(input_data, smoothing_val, make_plot = FALSE, print_res = FALSE){
    # This function uses a derived least-squares method with an Euler
    #  discretization in order to solve for the parameters of the logistic
    #  equation modeling the inputted data.
    # The input should be a data frame returned from one of the generators.
    
    validate_inputs(smoothing_val)
    
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
    A <- 1 - P_present
    ATA <- t(A) %*% A
    ID <- diag(nrow = nrow(ATA))
    # Would QR make a difference here? Probably not...
    # But when do you smooth vs. QR/SVD?
    # SVD and throwing away extreme eigenvalues.
    params <- MASS::ginv(ATA + smoothing_val * ID) %*% t(A) %*% b
    
    r_hat <- params[1]
    if (print_res) {
      cat("The estimated growth rate is:", r_hat,"\n")
    }
    # Optionally make a linear plot showing the data the regression equation
    #  is obtained from.
    if (make_plot == TRUE) {
      plot(1 - P_present, b)
    }
    
    # Collect outputs into a vector and return.
    output_params <- as.numeric(r_hat)
    return(output_params)
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
    validate_inputs(P0, r, max_t, time_step)
    
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

