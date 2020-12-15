library(deSolve)
library(tidyverse)
# Number of stochastic load profiles to generate
# nrep <- 100000
nrep <- 1000
# Load parameters
# gamma is dead-to-live load ratio 
gamma = 0.25
# R_0 is the characteristic strength
R_0 = 3000

alpha_d = 1.25; alpha_l = 1.5

## tau(t):  this will be stochastically generated for live load applications
zip <- function(a, b){
  idx <- order(c(seq_along(a), seq_along(b)))
  unlist(c(a,b)[idx])
}

# A = 0.0977, B = 5.1023 in Vancouver
# A = 0.3222, B = 17.0689 in Quebec City
A <- 0.3222
B <- 17.0689


# roof factor, convert ground snow load to roof snow load
r_flat_mean <- 0.6
r_flat_COV <- 0.45


sigmalog <- function(mu,  sigma){
  result <- sqrt(log(1 + exp(log(sigma*sigma) - 2 * log(mu))));
  result
}

Transfactor <- function(N,Bias,CoV){
  mu <- Bias
  sigma <- Bias * CoV
  input1 <- log(mu) - sigmalog(mu, sigma)*sigmalog(mu, sigma) / 2
  input2 <- sigmalog(mu, sigma)
  rlnorm(N,meanlog = input1, sdlog = input2)
}

# NS is the number of segment in each year
NS <- 10

# snow generation function
snow_generating<- function(N,NS,A,B){
  
  # p_0 is the probability of no-snow in one year (10 segment)
  p0 <- exp(-exp(A*B))
  
  # p_e is the probability of snow in a segment
  p_e <- 1.00 - exp(-1.00/NS*(exp(A*B)))
  
  # parameters for distribution g
  Bstr <- A*B/(A*B+3.9019)
  Astr <- A*B+3.9019
  
  # generate N random number to decide whether it snows in one segment
  r_n <- runif(N) 
  # generated snow laod
  g_s <- rep(0,N)
  # id that it snows in that segment. with r_n < p_e
  # it does not snow in that winter segment if r_n > p_e
  id_snow <- which (r_n <= p_e)
  # p is used to generate snow load
  p <- runif(length(id_snow))
  # roof factor that convert the ground snow load to the roof snow load 
  r <- Transfactor(length(p),Bias = r_flat_mean, CoV = r_flat_COV)
  # Final generation
  g_s[id_snow] <- r*(Bstr + 1.0 / Astr * (-log(-NS * log(1.0 - p_e + p_e * p))))
  g_s
}




# D_d is the standardized dead load assumed to follow N(1.05,0.1)
D_d <- list()
T_s <- list()
load_s <- list()

# the snow load or the wind load is named as D_e
# load_p is the segment snow/wind load
# T_e is time between the occurrence
# T_p is duration of the load
T_e <- list()
T_p <- list()
load_p <- list()
# D_e is the snow/wind load
D_e <- list()
N <- 50

for (i in 1:nrep) {
  # D_d is the standardized dead load assumed to follow N(1.05,0.1)
  D_d[[i]] <- rnorm(1, 1.05, 0.1)
  # generate some periods
  
  # the size of the load in each period
  
  # The occurrence between snow loads is 7 months exp(1/0.5833)
  T_e[[i]] <- rexp(N, 1/0.5833)
  # The duration time of each snow load segment is 2 weeks
  T_p[[i]] <- vector(mode = 'list',length = N) %>%
    lapply(function(x) x = rexp(NS, 1/0.03835))
  load_p[[i]] <- vector(mode = 'list',length = N) %>%
    lapply(function(x) x = snow_generating(10,10,A,B))
  
  D_e[[i]] <- stepfun(cumsum(zip(T_e[[i]], T_p[[i]])), zip(rep(0, N+1), load_p[[i]]))
  
}
save(D_d, D_e, gamma, R_0, alpha_d, alpha_l, file="rel_loads_snow_Queb.rda")

