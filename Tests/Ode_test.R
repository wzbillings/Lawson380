###
# Testing deSolve for logistic equation
# Author: Zane Billings
# Date created: 2019-10-19
# Last modified: 2019-10-19
# In order to debug problems with the set of functions being used for my 
#  research with Dr. Lawson, I am testing out the ode() function independently
#  to see if I can (A) recreate results from examples in the vignette and (B)
#  get the method to work so I can implement it elsewhere.
###

library(deSolve)

model_parameters <- c("r" = 0.2, "K" = 1000)
initial_state <- c("P" = 100)
time_steps <- seq(from = 0, to = 25, by = 0.1)
log_der <- function(t, P, parms) {
  with(
    data = as.list(c(P, parms)),
    expr = {
      # Calculate rate of pop growth.
      dP = r * P * (1 - P/K)
      # Return rate as a list.
      return(list(c(dP)))
    } # End expr argument.
  ) # End with.
} # End function.

output <- ode(y = initial_state,
              times = time_steps,
              func = log_der,
              parms = model_parameters)