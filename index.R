## ----setup, include=F-----------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache =T)
chosen_params <- list(b = 1/(76*365), mu = 1/(76*365), beta1 = 0.214, beta2 = 0.107, beta3 = 0.428, lambda = 1/7, alpha = 1/50, gamma = 1/14, omega = 1/365)
library(hmer)
library(lhs)
library(deSolve)
library(ggplot2)
library(reshape2)
library(purrr)
library(tidyverse)
library(progressr)
handlers("progress")
set.seed(123)
############################# HELPER FUNCTIONS #############################
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
  epsilon = 1/7,
  alpha = 1/50,
  gamma = 1/14,
  omega = 1/365
)


solution <- get_results(chosen_params, outs = c("I", "R"),
                        times = c(25, 40, 100, 200), raw = TRUE)



## 
## we get a 3-dimensional array `solution` such that `solution[t,j,i]` contains the number of individuals in the j-th compartment at time $t$ for the i-th run of the model at `chosen_params`. In particular, $t$ can take values $1,2,...,201$, $j$ can take values $1,2,3,4,5$ corresponding to $S, E, I, R, D$ (where $D$ stands for the cumulative number of deaths occurred), and $i$ can be $1,2,3,...,100$.


## ---- fig.height=7, fig.width=7-------------------------------------------------------------------------------------------------
plot(0:200, ylim=c(0,700), ty="n", xlab = "Time", ylab = "Number")
for(j in 3:4) for(i in 1:100) lines(0:200, solution[,j,i], col=(3:4)[j-2], lwd=0.3)
legend('topleft', legend = c('Infected', "Recovered"), lty = 1, 
       col = c(3,4), inset = c(0.05, 0.05))


## ---- fig.height=7, fig.width=7-------------------------------------------------------------------------------------------------
plot(0:200, ylim=c(0,1000), ty="n", xlab = "Time", ylab = "Number", main = "Susceptibles")
for(i in 1:100) lines(0:200, solution[,1,i], col='black', lwd=0.3, 
                      xlab = "Time", ylab = "Number", main = "Susceptibles")



## -------------------------------------------------------------------------------------------------------------------------------
ranges = list(
  b = c(1e-5, 1e-4), # birth rate
  mu = c(1e-5, 1e-4), # rate of death from other causes
  beta1 = c(0.05, 0.3), # infection rate at time t=0
  beta2 = c(0.1, 0.2), # infection rates at time t=100
  beta3 = c(0.3, 0.5), # infection rates at time t=270
  epsilon = c(0.01, 0.21), # rate of becoming infectious after infection
  alpha = c(0.01, 0.025), # rate of death from the disease
  gamma = c(0.01, 0.08), # recovery rate
  omega = c(0.002, 0.004) # rate at which immunity is lost following recovery
)


## -------------------------------------------------------------------------------------------------------------------------------
targets <- list(
  I25 = c(98.51, 133.25),
  I40 = c(117.17, 158.51),
  I100 = c(22.39, 30.29),
  I200 = c(0.578, 0.782),
  R25 = c(106.34, 143.9),
  R40 = c(218.28, 295.32),
  R100 = c(458.14, 619.84),
  R200 = c(377.6, 510.86)
)


## -------------------------------------------------------------------------------------------------------------------------------
initial_LHS_training <- maximinLHS(100, 9)
initial_LHS_validation <- maximinLHS(50, 9)
initial_LHS <- rbind(initial_LHS_training, initial_LHS_validation)


## -------------------------------------------------------------------------------------------------------------------------------
initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                                        function(x) x * purrr::map_dbl(ranges, diff) + 
                                        purrr::map_dbl(ranges, ~.[[1]])))), names(ranges))


## -------------------------------------------------------------------------------------------------------------------------------
initial_results <- list()
with_progress({
  p <- progressor(nrow(initial_points))
for (i in 1:nrow(initial_points)) {
  model_out <- get_results(unlist(initial_points[i,]), nreps = 15, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  initial_results[[i]] <- model_out
  p(message = sprintf("Run %g", i))
}
})


## -------------------------------------------------------------------------------------------------------------------------------
wave0 <- data.frame(do.call('rbind', initial_results))
all_training <- wave0[1:1500,]
all_valid <- wave0[1501:2250,]
output_names <- c("I25", "I40", "I100", "I200", "R25", "R40", "R100", "R200")


stoch_emulators <- variance_emulator_from_data(all_training, output_names, ranges)

plot_actives(stoch_emulators$variance)
plot_actives(stoch_emulators$expectation)


emulator_plot(stoch_emulators$expectation$I100, params = c('epsilon', 'alpha'),

     fixed_vals = all_training[1, names(ranges)[-c(6,7)]], plot_type = 'var') +

      geom_point(data = all_training[1,], aes(x = epsilon, y = alpha))

emulator_plot(stoch_emulators, plot_type = 'nimp', 
              targets = targets, params = c('epsilon', 'alpha'))

emulator_plot(stoch_emulators, plot_type = 'nimp', targets = targets,

              params = c('epsilon', 'alpha'),

               fixed_vals = chosen_params[!names(chosen_params) %in% c('epsilon', 'alpha')],

               cb=T) +geom_point(aes(x=1/7, y=1/50), size=3)

vd <- validation_diagnostics(stoch_emulators$expectation, targets, all_valid, plt=TRUE, row=2)

new_points <- generate_new_runs(stoch_emulators, 150, targets)

plot_wrap(new_points, ranges)

new_results <- list()
for (i in 1:nrow(new_points)) {
  model_out <- get_results(unlist(new_points[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  new_results[[i]] <- model_out
}
wave1 <- data.frame(do.call('rbind', new_results))
new_all_training <- wave1[1:5000,]
new_all_valid <- wave1[5001:7500,]

new_stoch_emulators <- variance_emulator_from_data(new_all_training, output_names, ranges, 
                                                   check.ranges=TRUE)


## ---- fig.height=6, fig.width=8-------------------------------------------------------------------------------------------------
vd <- validation_diagnostics(new_stoch_emulators, targets, new_all_valid, plt=TRUE, row=2)


## -------------------------------------------------------------------------------------------------------------------------------
new_new_points <- generate_new_runs(c(new_stoch_emulators, stoch_emulators), 150, targets)


## -------------------------------------------------------------------------------------------------------------------------------
new_new_results <- list()
for (i in 1:nrow(new_new_points)) {
  model_out <- get_results(unlist(new_new_points[i,]), nreps = 100, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  new_new_results[[i]] <- model_out
}
wave2 <- data.frame(do.call('rbind', new_new_results))
new_new_all_training <- wave2[1:10000,]
new_new_all_valid <- wave2[10001:15000,]
new_new_stoch_emulators <- variance_emulator_from_data(new_new_all_training, output_names, ranges, 
                                                   check.ranges=TRUE)
new_new_new_points <- generate_new_runs(c(new_new_stoch_emulators, new_stoch_emulators, stoch_emulators), 150, targets)

new_new_new_results <- list()
for (i in 1:nrow(new_new_new_points)) {
  model_out <- get_results(unlist(new_new_new_points[i,]), nreps = 150, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200))
  new_new_new_results[[i]] <- model_out
}
wave3 <- data.frame(do.call('rbind', new_new_new_results))
new_new_new_all_training <- wave3[1:15000,]
new_new_new_all_valid <- wave3[15001:22500,]
new_new_new_stoch_emulators <- variance_emulator_from_data(new_new_new_all_training, output_names, ranges, 
                                                   check.ranges=TRUE)
new_new_new_new_points <- generate_new_runs(c(new_new_new_stoch_emulators, new_new_stoch_emulators, new_stoch_emulators, stoch_emulators), 150, targets)



## ----  fig.height=8, fig.width=8------------------------------------------------------------------------------------------------
wave_points(list(initial_points, new_points, new_new_points, new_new_new_points, new_new_new_new_points), 
            input_names = names(ranges)[4:8]) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6))


## ----  fig.height=8, fig.width=8------------------------------------------------------------------------------------------------
new_new_training_list <- list()
with_progress({
  p <- progressor(nrow(new_new_points))
for (i in 1:nrow(new_new_points)) {
  model_out <- get_results(unlist(new_new_points[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200, 300, 400, 450))
  new_new_training_list[[i]] <- model_out
  p(message = sprintf("Run %g", i))
}
})
new_new_all_training <-  data.frame(do.call('rbind',new_new_training_list))
new_new_new_training_list <- list()
with_progress({
  p <- progressor(nrow(new_new_new_points))
for (i in 1:nrow(new_new_new_points)) {
  model_out <- get_results(unlist(new_new_new_points[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200, 300, 400, 450))
  new_new_new_training_list[[i]] <- model_out
  p(message = sprintf("Run %g", i))
}
})
new_new_new_all_training <-  data.frame(do.call('rbind',new_new_new_training_list))
new_new_new_new_training_list <- list()
with_progress({
  p <- progressor(nrow(new_new_new_new_points))
for (i in 1:nrow(new_new_new_new_points)) {
  model_out <- get_results(unlist(new_new_new_new_points[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200, 300, 400, 450))
  new_new_new_new_training_list[[i]] <- model_out
  p(message = sprintf("Run %g", i))
}
})
new_new_new_new_all_training <-  data.frame(do.call('rbind',new_new_new_new_training_list))



all_training_aggregated <- aggregate_points(all_training, names(ranges))

new_all_training_aggregated <- aggregate_points(new_all_training, names(ranges))

new_new_all_training_aggregated <- aggregate_points(new_new_all_training, names(ranges))

 new_new_new_all_training_aggregated <- aggregate_points(new_new_new_all_training, names(ranges))

new_new_new_new_all_training_aggregated <- aggregate_points(new_new_new_new_all_training, names(ranges))

 all_aggregated <- list(all_training_aggregated, new_all_training_aggregated, new_new_all_training_aggregated,

                       new_new_new_all_training_aggregated)

## ``
 all_aggregated <- list(all_training_aggregated, new_all_training_aggregated, new_new_all_training_aggregated)
                        
## 
## We are now ready to compare the performance of parameter sets at the end of each wave, passing the list `all_aggregated` and the list `targets` to the function `simulator_plot`:

## 
## ``{r,  fig.height=8, fig.width=8}
all_aggregated <- list(all_training_aggregated, new_all_training_aggregated, new_new_new_new_all_training_aggregated)

 simulator_plot(all_aggregated, targets, barcol = "grey")

## ``

## ----  fig.height=8, fig.width=8------------------------------------------------------------------------------------------------
wave_values(all_aggregated, targets, l_wid=1) 


## Since this is an iterative process, at the end of each wave we need to decide whether to perform a new wave or to stop. One possible stopping criterion consists of comparing the emulator uncertainty and the target uncertainty. If the former is larger, another wave can be performed, since new, more confident emulators can potentially help further reduce the non-implausible space. If the uncertainty of emulators is smaller than the uncertainty in the targets, improving the performance of emulators would not make a substantial difference, and additional waves would not be beneficial. We may also choose to stop the iterations when we get emulators that provide us with full fitting points at a sufficiently high rate. In such a case, rather than spending time training new emulators, we can generate new points with the function `generate_new_runs` using the current emulators until we find enough full fitting ones. Finally, we might end up with all the input space deemed implausible at the end of a wave. In this situation, we would deduce that there are no parameter sets that give an acceptable match with the data: in particular, this would raise doubts about the adequacy of the chosen model, or input and/or output ranges.


## ---- fig.height=7, fig.width=7-------------------------------------------------------------------------------------------------
solution <- get_results(chosen_params, outs = c("I", "R"), 
                        times = c(25, 40, 100, 200, 300, 400, 500), raw = TRUE)
plot(0:500, ylim=c(0,700), ty="n", xlab = "Time", ylab = "Number")
for(j in 3:4) for(i in 1:100) lines(0:500, solution[,j,i], col=(3:4)[j-2], lwd=0.3)
legend('topleft', legend = c('Infected', "Recovered"), lty = 1, 
       col = c(3,4), inset = c(0.05, 0.05))


## Investigate how different parameter sets give rise to different levels of bimodality.

## 
## ``{info, title = "R tip", collapsible = TRUE, ref.label = "rtip1"}

## ``

## 

## Let us see what happens when we simulate a higher rate of infection between each infectious and susceptible:

## 
## ``{r, fig.height=8, fig.width=8}

## higher_beta_params <- c(

##   b = 1/(76*365),

##   mu = 1/(76*365),

##   beta1 = 0.3, beta2 = 0.1, beta3 = 0.5,

##   sigma = 1/7,

##   alpha = 1/50,

##   gamma = 1/14,

##   omega = 1/365

## )

## higher_beta_solution <- get_results(higher_beta_params, outs = c("I", "R"),

##                         times = c(25, 40, 100, 200, 300, 400, 500), raw = TRUE)

## plot(0:500, ylim=c(0,700), ty="n", xlab = "Time", ylab = "Number")

## for(j in 3:4) for(i in 1:100) lines(0:500, solution[,j,i], col=(3:4)[j-2], lwd=0.3)

## for(j in 3:4) for(i in 1:100) lines(0:500, higher_beta_solution[,j,i],

##                                     col=(c("#E69F00","#CC79A7"))[j-2], lwd=0.3)

## legend('topleft', legend = c('Infected', "Recovered", "Infected - higher beta",

##                              "Recovered - higher beta"), lty = 1,

##        col = c(3,4,"#E69F00","#CC79A7"), inset = c(0.01, 0.01))

## ``

## 
## We see that the peaks of recovered and infectious individuals during the first wave have now increased, as expected. Also note how increasing the transmission rate reduced the role of bimodality (very few trajectories in orange take off after the first wave): the fact that more individuals get now infected during the first wave (around $t=50$), makes the chances of a second "wave" around time $t=300$ lower.


## -------------------------------------------------------------------------------------------------------------------------------
bimodal_targets <- list(
  I25 = list(val = 115.88, sigma = 5.79),
  I40 = list(val = 137.84, sigma = 6.89),
  I100 = list(val = 26.34, sigma = 1.317),
  I200 = list(val = 0.68, sigma = 0.034),
  I300 = list(val = 29.55, sigma = 1.48),
  I400 = list(val = 15.67, sigma = 5.36),
  R25 = list(val = 125.12, sigma = 6.26),
  R40 = list(val = 256.80, sigma = 12.84),
  R100 = list(val = 538.99, sigma = 26.95),
  R200 = list(val = 444.23, sigma = 22.21),
  R300 = list(val = 371.08, sigma = 15.85),
  R400 = list(val = 569.39, sigma = 26.52)
)


## -------------------------------------------------------------------------------------------------------------------------------
bimodal_initial_results <- list()
with_progress({
  p <- progressor(nrow(initial_points))
for (i in 1:nrow(initial_points)) {
  model_out <- get_results(unlist(initial_points[i,]), nreps = 50, outs = c("I", "R"), 
                           times = c(25, 40, 100, 200, 300, 400, 450))
  bimodal_initial_results[[i]] <- model_out
  p(message = sprintf("Run %g", i))
}
})
bimodal_wave0 <- data.frame(do.call('rbind', bimodal_initial_results))
bimodal_all_training <- bimodal_wave0[1:5000,]
bimodal_all_valid <- bimodal_wave0[5001:7500,]
bimodal_output_names <- c("I25", "I40", "I100", "I200", "I300", "I400", "I450", 
                          "R25", "R40", "R100", "R200", "R300", "R400", "R450")


## To train bimodal emulators we use the function `bimodal_emulator_from_data`, which requires the training data, the names of the outputs to emulate, and the ranges of the parameters:

## 
## ``{r, fig.height=7, fig.width=7}

## bimodal_emulators <- bimodal_emulator_from_data(bimodal_all_training,

##                                                 bimodal_output_names, ranges)

## ``

## 
## Behind the scenes, this function does the following:

## 
## - First it looks at the provided training data and identifies which of the outputs are bimodal.

## 
## - For the outputs where bimodality is found, the repetitions at each parameter set are clustered into two subsets, based on the mode they belong to.

## 
## - For outputs without bimodality, stochastic emulators are trained as before. For outputs with bimodality, variance and mean emulators are trained separately for each of the two modes. To access the mean emulators for the first mode, we type `bimodal_emulators$mode1$expectation`, and to access the variance emulators for the second mode we type `bimodal_emulators$mode2$variance`.

## 
## - Finally, an emulator for the proportion of points in each mode is also trained (this is a single emulator, as in the deterministic case). This emulator can be accessed by typing `bimodal_emulators$prop` and will be referred to as proportion emulator.

## 
## Let us now plot the expectation of the mean emulator of $R400$ and similarly for each mode:

## 
## ``{r, fig.height=7, fig.width=9}

## emulator_plot(bimodal_emulators$mode1$expectation$R400, params = c('alpha', 'epsilon'))

## emulator_plot(bimodal_emulators$mode2$expectation$R400, params = c('alpha', 'epsilon'))

## ``

## 
## We immediately notice that there is large difference between the two plots, with mode 2 containing higher values than mode 1. This indicates that mode 1 is the one where the infection dies out.

## 
## It is also instructive to plot the expectation of the proportion emulator. If we use `plot_actives`:

## 
## ``{r, fig.height=2, fig.width=9}

## plot_actives(bimodal_emulators$prop)

## ``

## 
## we see that $\beta_2$ and $\omega$ are active variables for the proportion emulator. So let us plot such emulator on those two parameters:

## 

## ``{r, fig.height=7, fig.width=9}

## emulator_plot(bimodal_emulators$prop, params = c('beta2', 'omega'))

## ``

## 
## The plot shows that at points where epsilon (rate of becoming infectious) and alpha (rate of death from the disease) are high, the first mode, where the infection dies out, dominates. However, from the plot it is also evident that there is always a

## substantial risk of the disease dying out, since the proportion of points in the first mode never drops below $0.5$.

## 
## Let us conclude this section with a plot showing the predictions and relative uncertainties for mode 1 and mode 2 for the $I$ and $R$ outputs. We first create a dataframe with the parameter set `chosen_params` and extract the emulators' predictions (with `get_exp`) and uncertainties (with `get_cov`) at such parameter set, for each mode:

## 

## ``{r, fig.height=7, fig.width=7}

## chosen_df <- data.frame(b = 1/(76*365),

##                          mu = 1/(76*365),

##                          beta1 = 0.214, beta2 = 0.107, beta3 = 0.428,

##                          epsilon = 1/7,

##                          alpha = 1/50,

##                          gamma = 1/14,

##                          omega = 1/365

## )

## em_preds1 <- purrr::map_dbl(bimodal_emulators$mode1$expectation, ~.$get_exp(chosen_df))

## em_vars1 <- purrr::map_dbl(bimodal_emulators$mode1$expectation, ~.$get_cov(chosen_df))

## em_preds2 <- purrr::map_dbl(bimodal_emulators$mode2$expectation, ~.$get_exp(chosen_df))

## em_vars2 <- purrr::map_dbl(bimodal_emulators$mode2$expectation, ~.$get_cov(chosen_df))

## ``

## We then plot:

## 
## - the solutions for "I" and "R" at `chosen_params`,

## 
## - unbroken red and black lines representing the mode 1 and mode 2 expectations for the $I$ and $R$ outputs,

## 
## - dotted red and black lines representing $3$-sigma upper and lower bounds for those expectations.

## 
## 
## ``{r, fig.height=7, fig.width=7}

## plot(0:500, ylim=c(0,700), ty="n", xlab = "Time", ylab = "Number")

## for(j in 3:4) for(i in 1:100) lines(0:500, solution[,j,i], col=(3:4)[j-2], lwd=0.3)

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds1[1:7])

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds1[8:14])

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds1[1:7]-3*sqrt(em_vars1[1:7]), lty = 2)

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds1[1:7]+3*sqrt(em_vars1[1:7]), lty = 2)

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds1[8:14]-3*sqrt(em_vars1[8:14]), lty = 2)

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds1[8:14]+3*sqrt(em_vars1[8:14]), lty = 2)

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds2[1:7], col = 'red')

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds2[8:14], col = 'red')

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds2[1:7]-3*sqrt(em_vars1[1:7]),

##       lty = 2, col = 'red')

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds2[1:7]+3*sqrt(em_vars1[1:7]),

##       lty = 2, col = 'red')

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds2[8:14]-3*sqrt(em_vars1[8:14]),

##       lty = 2, col = 'red')

## lines(c(25, 40, 100, 200, 300, 400, 450), em_preds2[8:14]+3*sqrt(em_vars1[8:14]),

##       lty = 2, col = 'red')

## legend('topleft', legend = c("Infected", "Recovered"), lty = 1, col = c(3,4),

##        inset = c(0.05, 0.05))

## ``

## 
## We clearly see how `bimodal_emulators_from_data` captures the bimodal structure of the outputs.


## Confirm that mode 1 is where the infection dies out by comparing the plots of the expectation of the mean emulator for $I300$ for both modes.


## We use `emulator_plot` to produce the suggested plots:

## 
## ``{r, fig.height=8, fig.width=8}

## emulator_plot(bimodal_emulators$mode1$expectation$I300, params = c('alpha', 'epsilon'))

## emulator_plot(bimodal_emulators$mode2$expectation$I300, params = c('alpha', 'epsilon'))

## ``

## 
## The values in the plot for mode 1 are clearly lower than those for mode 2: this reflects the fact that mode 2 experiences a second wave around $t=300$, while in mode 1 the infection dies out after the first wave (around $t=40$).


## 
## Since we are now in a bimodal setting, we want to regard a point as non-implausible for a given target if it is valid with respect to either mode. This is the default behaviour of the package, when dealing with bimodal emulators.

## 
## For a given output and a given point, the implausibility for mode 1 and the implausibility for mode 2 are calculated, and the minimum of them is taken. The maximum (or second-maximum, third-maximum etc) of this collection of minimised implausibilities is then selected, depending on the user choice. For example, to plot the maximum of these minimised implausibilities, we set `plot_type` to `nimp`:

## 
## ``{r, fig.height=7, fig.width=9}

## emulator_plot(subset_emulators(bimodal_emulators, bimodal_output_names[-c(7,14)]),

##               plot_type = 'nimp', targets = bimodal_targets, params = c('alpha', 'epsilon'))

## ``

## 
## Here we used the function `subset_emulators` to remove the emulators for the outputs at $t=450$ (the seventh and fourteenth emulators), for which we did not set targets (to keep the symmetry with the deterministic case).


## Set the argument `plot_type` to `imp` to produce implausibility plots for each output, for mode 1 and mode 2. Which implausibility plots are the same for each mode, and which are different? Why?


## Setting `plot_type='imp'` in `emulator_plot` and passing it either `subset_emulators(bimodal_emulators, output_names[-c(7,14)])$mode1` or `subset_emulators(bimodal_emulators, output_names[-c(7,14)])$mode2`, we get:

## 
## ``{r fig.width = 7, fig.height = 6}

## emulator_plot(subset_emulators(bimodal_emulators, bimodal_output_names[-c(7,14)])$mode1, plot_type = 'imp',

##               targets = bimodal_targets, params = c('alpha', 'epsilon'))

## emulator_plot(subset_emulators(bimodal_emulators, bimodal_output_names[-c(7,14)])$mode2, plot_type = 'imp',

##               targets = bimodal_targets, params = c('alpha', 'epsilon'))

## ``

## The implausibility plots are the same between modes for early times ($t=25,40,100$) and different between modes for later outputs. This makes sense, since bimodality, which enters the picture after the first wave, does not play a role in earlier times.


## ----  fig.height=8, fig.width=8------------------------------------------------------------------------------------------------
vd <- validation_diagnostics(subset_emulators(bimodal_emulators, bimodal_output_names[-c(7,14)]), 
                       bimodal_targets, bimodal_all_valid, plt=TRUE)


## 
## You may have noticed that in the first column we have more points plotted than we have validation points. This is because we do diagnostics for each mode for each output, i.e. we compare the predictions of the mean emulators for mode 1 and mode 2 with the model output values, which are also clustered in two subsets, due to bimodality. In fact, for $I400$ you can see the mode separation in the left plot (some runs are concentrated around zero and the rest are all above 15). Similarly, for $I200$ you can distinguish the two modes by the step change in uncertainty (points on the left, near zero, have smaller intervals than points more on the right).


## -------------------------------------------------------------------------------------------------------------------------------
bimodal_emulators$mode1$expectation$I200 <- bimodal_emulators$mode1$expectation$I200$mult_sigma(2)
bimodal_emulators$mode2$expectation$I200 <- bimodal_emulators$mode2$expectation$I200$mult_sigma(2)


## ----cache=F--------------------------------------------------------------------------------------------------------------------
new_points <- generate_new_runs(subset_emulators(bimodal_emulators, 
                                bimodal_output_names[-c(7,14)]), 200, bimodal_targets, nth=1)

