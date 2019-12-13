###
# Zane Billings
# Logistic data fitting
# Created 2019-10-06
#  Fitting a least-squares model to logistic time-series data in several 
#  different ways to see if we can get known parameters back.
###

require(MASS) # ginv() function source for generalize inverse computation.

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
      plot(P_present, 1/b)
    }
    
    # Collect outputs into a vector and return.
    output_params <- unlist(
      list("r" = as.numeric(r_hat), 
           "K" = as.numeric(K_hat))
    )
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
    plot(P_present, 1/b)
  }
  
  # Collect outputs into a vector and return.
  output_params <- unlist(
    list("r" = as.numeric(r_hat), 
         "K" = as.numeric(K_hat))
  )
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
    plot(P_present, 1/b)
  }
  
  # Collect outputs into a vector and return.
  output_params <- unlist(
    list("r" = as.numeric(r_hat), 
         "K" = as.numeric(K_hat))
  )
  return(output_params)
}

model_logistic_data_smoothing <- 
  function(input_data, smoothing_value, make_plot = FALSE){
    # This function uses a derived least-squares method with an Euler
    #  discretization in order to solve for the parameters of the logistic
    #  equation modeling the inputted data.
    # This method implements Tihkonov Regularization (l2 smoothing) on the
    #  least-squares parameters.
    # The input should be a list returned from prep_data().
    
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
    ATA <- t(A) %*% A
    ID <- diag(nrow(ATA))
    params <- MASS::ginv(ATA + smoothing_value * ID) %*% t(A) %*% b
    
    # Now knowing that m = r/K, we solve for K.
    r_hat <- params[1]
    m_hat <- params[2]
    K_hat <- (r_hat/m_hat)
    cat("The estimated growth rate is:", r_hat,
        "\nThe estimated carrying capacity is:", K_hat, "\n")
    
    # Optionally make a linear plot showing the data the regression equation
    #  is obtained from.
    if (make_plot == TRUE) {
      plot(P_present, 1/b)
    }
    
    # Collect outputs into a list and return.
    output_params <- unlist(
      list("r" = as.numeric(r_hat), 
           "K" = as.numeric(K_hat))
      )
    return(output_params)
  }