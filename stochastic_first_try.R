library(hmer)
library(lhs)
library(deSolve)
library(ggplot2)
library(reshape2)
library(purrr)
library(tidyverse)
set.seed(123)
aggregate_points <- function(data, input_names, func = 'mean') {
  unique_points <- unique(data[,input_names])
  uids <- apply(unique_points, 1, rlang::hash)
  data_by_point <- purrr::map(uids, function(h) {
    data[apply(data[,input_names], 1, rlang::hash) == h,]
  })
  aggregate_func <- get(func)
  aggregate_values <- do.call('rbind.data.frame', purrr::map(data_by_point, function(x) {
    output_data <- x[,!(names(x) %in% input_names)]
    apply(output_data, 2, aggregate_func)
  }))
  return(setNames(cbind.data.frame(unique_points, aggregate_values), names(data)))
}

############################# HELPER FUNCTIONS #############################
# Gillespie algorithm implementation
# N is a list containing N$M, the initial number of people in each compartment (comprising deaths), N$Post, N$Pre, N$h, T is the
# time up to which we want to simulate (400th day by default). The while loop is there to "record" only the wanted times: for 
# example, if dt=1, then we want to record only the state of the system at t=0,1,2,3,4,5,... and if the rexp(1,h0) is not enough
# to send us to the next integer time, then we will do more than one cycle and record the state only when we reach the next 
# integer time.
gillespied=function (N, T=400, dt=1, ...)
{
  tt=0
  n=T%/%dt
  x=N$M
  S=t(N$Post-N$Pre)
  u=nrow(S)
  v=ncol(S)
  xmat=matrix(0,ncol=u,nrow=n)
  i=1
  target=0
  repeat {
    h=N$h(x, tt, ...)
    h0=sum(h)
    if (h0<1e-10)
      tt=1e99
    else
      tt=tt+rexp(1,h0)
    while (tt>=target) {
      xmat[i,]=x
      i=i+1
      target=target+dt
      if (i>n)
        #			        return(ts(xmat,start=0,deltat=dt))
        return(xmat)
    }
    j=sample(v,1,prob=h)
    x=x+S[,j]
  }
}

# Initial setup: transition matrices, hazard function, and initial conditions
Num <- 1000
N=list()
N$M=c(900,100,0,0,0)
# N$Pre=matrix(c(1,0,1,1,0,1),ncol=2,byrow=TRUE)
# N$Post=matrix(c(2,0,0,2,0,0),ncol=2,byrow=TRUE)
# Columns of N$Pre are 5, for S, E, I, R, and Deaths. Rows are for transitions: births, S to E, S to D, E to I, E to D, I to D,
# I to R, R to S, R to D. For example: the first row is all of zeros, since we need no person in S, nor in the other compartments 
# for a new birth to happen. 
# Columns of N$Post are 5 as for N$Pre. Rows are for transitions as for N$Pre. For example: the element in (1,1) is 1 since a new 
# birth goes straight into the S compartment. 
N$Pre = matrix(c(0,0,0,0,0,
                 1,0,0,0,0,
                 1,0,0,0,0,
                 0,1,0,0,0,
                 0,1,0,0,0,
                 0,0,1,0,0,
                 0,0,1,0,0,
                 0,0,0,1,0,
                 0,0,0,1,0), ncol = 5, byrow = TRUE)
N$Post = matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,0,0,1,
                  0,0,1,0,0,
                  0,0,0,0,1,
                  0,0,0,0,1,
                  0,0,0,1,0,
                  1,0,0,0,0,
                  0,0,0,0,1), ncol = 5, byrow = TRUE)
# N$Pre=matrix(c(1,0,0,0,1,0,1,0,0),ncol=3,byrow=TRUE)
# N$Post=matrix(c(0,1,0,0,0,1,0,0,1),ncol=3,byrow=TRUE)
# Here x is the vector with the number of people in each of the compartments (deaths excluded), t is the time and th is a vector
# with the parameters. N$h return a vector with all the transition rates at time t.
N$h=function(x,t,th=rep(1,9))
{
  Num = x[1]+x[2]+x[3]+x[4]
  if (t > 270) tns <- th[5]
  else if (t > 180) tns <- (th[5]-th[4])*t/90+3*th[4]-2*th[5]
  else if (t > 100) tns <- th[4]
  else tns <- (th[4]-th[3])*t/100+th[3]
  return(
    c(th[1]*Num,
      tns*x[3]*x[1]/Num,
      th[2]*x[1],
      th[6]*x[2],
      th[2]*x[2],
      (th[7]+th[2])*x[3],
      th[8]*x[3],
      th[9]*x[4],
      th[2]*x[4])
  )
}

get_results <- function(params, nreps = 100, outs, times, raw = FALSE) {
  tseq <- 0:max(times)
  arra <- array(0, dim = c(max(tseq)+1, 5, nreps))
  for(i in 1:nreps) arra[,,i] <- gillespied(N,T=max(times) + 1 + 0.001,dt=1,th=params)
  if(raw) return(arra)
  collected <- list()
  for (i in 1:nreps) {
    relev <- c(arra[times+1, which(c("S", "E", "I", "R", "D") %in% outs), i])
    names <- unlist(purrr::map(outs, ~paste0(., times, sep = "")))
    relev <- setNames(relev, names)
    collected[[i]] <- relev
  }
  input_dat <- setNames(data.frame(matrix(rep(params, nreps), ncol = length(params), byrow = TRUE)), names(params))
  return(cbind(input_dat, do.call('rbind', collected)))
}

chosen_params <- c(
  b = 1/(76*365),
  mu = 1/(76*365),
  beta1 = 0.214, beta2 = 0.107, beta3 = 0.428,
  sigma = 1/7,
  alpha = 1/50,
  gamma = 1/14,
  omega = 1/365
)

solution <- get_results(chosen_params, outs = c("I", "R"), 
                        times = c(25, 40, 100, 200), raw = TRUE)

plot(0:200, ylim=c(0,700), ty="n", xlab = "Time", ylab = "Number")
for(j in 3:4) for(i in 1:100) lines(0:200, solution[,j,i], col=(3:4)[j-2], lwd=0.3)
legend('topleft', legend = c('Infected', "Recovered"), lty = 1, 
       col = c(3,4), inset = c(0.05, 0.05))

plot(0:200, ylim=c(0,1000), ty="n", xlab = "Time", ylab = "Number", main = "Susceptibles")
for(i in 1:100) lines(0:200, solution[,1,i], col='black', lwd=0.3, 
                      xlab = "Time", ylab = "Number", main = "Susceptibles")

#this is for a task, not important (till line 185)
higher_foi_params <- c(
  b = 1/(76*365),
  mu = 1/(76*365),
  beta1 = 0.3, beta2 = 0.1, beta3 = 0.5,
  sigma = 1/7,
  alpha = 1/50,
  gamma = 1/14,
  omega = 1/365
)
higher_foi_solution <- get_results(higher_foi_params, outs = c("I", "R"), 
                                   times = c(25, 40, 100, 200), raw = TRUE)
plot(0:200, ylim=c(0,700), ty="n", xlab = "Time", ylab = "Number")
for(j in 3:4) for(i in 1:100) lines(0:200, higher_foi_solution[,j,i], col=(3:4)[j-2], lwd=0.3)
legend('topleft', legend = c('Infected', "Recovered"), lty = 1, 
       col = c(3,4), inset = c(0.05, 0.05))

smaller_gamma_params <- c(
  b = 1/(60*365),
  mu = 1/(76*365),
  beta1 = 0.214, beta2 = 0.107, beta3 = 0.428,
  sigma = 1/7,
  alpha = 1/50,
  gamma = 1/20,
  omega = 1/365
)
smaller_gamma_solution <- get_results(smaller_gamma_params, outs = c("I", "R"), 
                                      times = c(25, 40, 100, 200), raw = TRUE)
plot(0:200, ylim=c(0,700), ty="n", xlab = "Time", ylab = "Number")
for(j in 3:4) for(i in 1:100) lines(0:200, smaller_gamma_solution[,j,i], 
                                    col=(3:4)[j-2], lwd=0.3)
legend('topleft', legend = c('Infected', "Recovered"), lty = 1, 
       col = c(3,4), inset = c(0.05, 0.05))

# same ranges as deterministic workshop
ranges = list(
  b = c(1e-5, 1e-4), # birth rate
  mu = c(1e-5, 1e-4), # rate of death from other causes
  beta1 = c(0.2, 0.3), # infection rate at time t=0
  beta2 = c(0.1, 0.2), # infection rates at time t=100
  beta3 = c(0.3, 0.5), # infection rates at time t=270
  sigma = c(0.07, 0.21), # rate of becoming infectious after infection
  alpha = c(0.01, 0.025), # rate of death from the disease
  gamma = c(0.05, 0.08), # recovery rate
  omega = c(0.002, 0.004) # rate at which immunity is lost following recovery
)

#only considering targets up to t=200, to eliminate bimodality
targets <- list(
  I25 = list(val = 115.88, sigma = 5.79),
  I40 = list(val = 137.84, sigma = 6.89),
  I100 = list(val = 26.34, sigma = 1.317),
  I200 = list(val = 0.68, sigma = 0.034),
  R25 = list(val = 125.12, sigma = 6.26),
  R40 = list(val = 256.80, sigma = 12.84),
  R100 = list(val = 538.99, sigma = 26.95),
  R200 = list(val = 444.23, sigma = 22.21)
)

# targets <- list(
#   I25 = list(val = 115.88, sigma = 5.79/2),
#   I40 = list(val = 137.84, sigma = 6.89/2),
#   I100 = list(val = 26.34, sigma = 1.317/2),
#   I200 = list(val = 0.68, sigma = 0.034/2),
#   R25 = list(val = 125.12, sigma = 6.26/2),
#   R40 = list(val = 256.80, sigma = 12.84/2),
#   R100 = list(val = 538.99, sigma = 26.95/2),
#   R200 = list(val = 444.23, sigma = 22.21/2)
# )
initial_lhs <- maximinLHS(100, 9)
validation_lhs <- lhs::randomLHS(50, 9)

initial_points <- setNames(data.frame(t(apply(initial_lhs, 1, 
                                              function(x) x * purrr::map_dbl(ranges, diff) + 
                                              purrr::map_dbl(ranges, ~.[[1]])))), names(ranges))

validation_points <- setNames(data.frame(t(apply(validation_lhs, 1, 
                                                 function(x) x * purrr::map_dbl(ranges, diff) + 
                                                 purrr::map_dbl(ranges, ~.[[1]])))), names(ranges))

output_list <- list()
for (i in 1:nrow(initial_points)) {
  model_out <- get_results(unlist(initial_points[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  output_list[[i]] <- model_out
}
validation_list <- list()
for (i in 1:nrow(validation_points)) {
  model_out <- get_results(unlist(validation_points[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  validation_list[[i]] <- model_out
}


all_output <- data.frame(do.call('rbind', output_list))
all_valid <- data.frame(do.call('rbind', validation_list))
output_names <- c("I25", "I40", "I100", "I200",
                  "R25", "R40", "R100", "R200")

stoch_emulators <- variance_emulator_from_data(all_output, output_names, ranges)

plot_actives(stoch_emulators$variance)
plot_actives(stoch_emulators$expectation)

emulator_plot(stoch_emulators$expectation$I40, params = c('sigma', 'alpha'),
     fixed_vals = all_output[1, names(ranges)[-c(6,7)]], plot_type = 'var') +
   geom_point(data = all_output[1,], aes(x = sigma, y = alpha))

chosen_df <- data.frame(b = 1/(76*365),
                        mu = 1/(76*365),
                        beta1 = 0.214, beta2 = 0.107, beta3 = 0.428,
                        sigma = 1/7,
                        alpha = 1/50,
                        gamma = 1/14,
                        omega = 1/365
)

emulator_plot(subset_emulators(stoch_emulators, output_names), plot_type = 'nimp', 
              targets = targets, params = c('alpha', 'sigma'))

emulator_plot(subset_emulators(stoch_emulators, output_names), plot_type = 'nimp', targets = targets,
              params = c('alpha', 'sigma'),
              fixed_vals = chosen_params[!names(chosen_params) %in% c('alpha', 'sigma')], 
              cb=T) +geom_point(aes(x=1/7, y=1/50), size=3)

vd <- validation_diagnostics(subset_emulators(stoch_emulators, output_names), 
                       targets, all_valid, plt=TRUE)

new_points <- generate_new_runs(stoch_emulators, 200, targets)

plot_wrap(new_points, ranges)

min_val <- list()
max_val <- list()
new_ranges <- list()
for (i in 1:length(ranges)) {
    par <- names(ranges)[[i]]
    min_val[[par]] <- max(min(new_points[,par])-0.05*diff(range(new_points[,par])), 
                      ranges[[par]][1])
    max_val[[par]] <- min(max(new_points[,par])+0.05*diff(range(new_points[,par])),
                      ranges[[par]][2])
    new_ranges[[par]] <- c(min_val[[par]], max_val[[par]])
}

t_sample <- sample(1:nrow(new_points), 100)
new_training <- new_points[t_sample,]
new_validation <- new_points[-t_sample,]

new_output_list <- list()
for (i in 1:nrow(new_training)) {
  new_model_out <- get_results(unlist(new_training[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  new_output_list[[i]] <- new_model_out
}
new_validation_list <- list()
for (i in 1:nrow(new_validation)) {
  new_model_out <- get_results(unlist(new_validation[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  new_validation_list[[i]] <- new_model_out
}
new_all_output <-  data.frame(do.call('rbind',new_output_list))
new_all_valid <- data.frame(do.call('rbind', new_validation_list))

new_stoch_emulators <- variance_emulator_from_data(new_all_output, output_names, new_ranges)

vd <- validation_diagnostics(new_stoch_emulators, targets, new_all_valid, plt=TRUE)

#modifying 2 of the emulators to make them more robust
# new_stoch_emulators$expectation$R40 <- new_stoch_emulators$expectation$R40$mult_sigma(2)
# new_stoch_emulators$expectation$I100 <- new_stoch_emulators$expectation$I100$mult_sigma(2)
# new_stoch_emulators$expectation$I200 <- new_stoch_emulators$expectation$I200$mult_sigma(2)
# 
# new_stoch_emulators$expectation$R100 <- new_stoch_emulators$expectation$R100$mult_sigma(6)
# new_stoch_emulators$expectation$R200 <- new_stoch_emulators$expectation$R200$mult_sigma(6)


#vd <- validation_diagnostics(new_stoch_emulators, targets, all_valid, plt=TRUE)

new_new_points <- generate_new_runs(c(new_stoch_emulators,stoch_emulators), 200, targets)

plot_wrap(new_new_points, ranges)

#new new points are kind of the same as new points unfortunately
wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges))

new_new_output_list <- list()
for (i in 1:nrow(new_new_points)) {
  new_new_model_out <- get_results(unlist(new_new_points[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  new_new_output_list[[i]] <- new_new_model_out
}
new_new_all_output <-  data.frame(do.call('rbind',new_new_output_list))
all_points <- list(all_output, new_all_output, new_new_all_output)

#aggregating results with your function
x1 <- aggregate_points(all_output, names(ranges))
x2 <- aggregate_points(new_all_output, names(ranges))
x3 <- aggregate_points(new_new_all_output, names(ranges))

#plotting the aggregated results. Unfortunately again wave 2 seems not to improve wave 1.
simulator_plot(list(x1,x2,x3),targets, barcol = "white")

