source("Sample_data_generation.R")

P_naught = 100
K = 1000
r = 0.2
max_t = 25
time_step = 0.1

generate_RK4_logistic_data(P_naught, K, r, max_t, time_step)
calc_logistic_derivative(P_naught,K,r,0)


# Solve P' = P with ODE
ex_func <- function(t,P) {
  return(P)
}

deSolve::ode(y = c("P" = 10),
             times = 1:25,
             func = ex_func,
             parms = c())






