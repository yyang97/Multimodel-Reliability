library(deSolve)

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


G_cov <- 0.150
G_mean <- (1+3.050*G_cov)/(1+2.592*G_cov)
A <- 1.282/(G_cov*G_mean)
B <- G_mean - (0.577/A) 
NS <- 12



wind_generating <- function(N,NS,A,B){
  
# p_0 is the probability of no-snow in one year (10 segment)
  p0 <- exp(-exp(A*B))
  
  # p_e is the probability of snow in a segment
  p_e <- 1.00 - exp(-1.00/NS*(exp(A*B)))
  
  # parameters for distribution W_50
  Wind_velocity_50 <- B + 1.0 / A * (-log(-NS * log(1.0 - p_e + p_e * 49/50)))
  
  # generate N random number to decide whether it snows in one segment
  r_n <- runif(N) 
  # generated snow laod
  g_s <- rep(0,N)
  # id that it snows in that segment. with r_n < p_e
  # it does not snow in that winter segment if r_n > p_e
  id_snow <- which (r_n <= p_e)
  # p is used to generate snow load
  p <- runif(length(id_snow))
  # Final generation
  Wind_velocity <- (B + 1.0 / A * (-log(-NS * log(1.0 - p_e + p_e * p))))
  g_s[id_snow] <- Wind_velocity^2*Transfactor(length(p),0.68,0.22)/
    ((Wind_velocity_50)^2*0.68)
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
N <- 100

for (i in 1:nrep) {
  # D_d is the standardized dead load assumed to follow N(1.05,0.1)
  D_d[[i]] <- rnorm(1, 1.05, 0.1)
  # generate some periods
  
  # the size of the load in each period
  
  # The occurrence between wind loads is 29 days exp(1/0.07945)
  T_e[[i]] <- rexp(N, 1/0.07945)
  # The duration time of each wind load segment is 1 day exp(1/0.00274)
  T_p[[i]] <- rexp(N, 1/0.00274)
  load_p[[i]] <- wind_generating(N,NS,A,B)
  
  D_e[[i]] <- stepfun(cumsum(zip(T_e[[i]], T_p[[i]])), zip(rep(0, N+1), load_p[[i]]))
  
}
save(D_d, D_e, gamma, R_0, alpha_d, alpha_l, file="rel_loads_wind.rda")







##----------------------generate snow load-------------
N <- 100
# zip modified so that it can be applied to two lists
zip_modified <- function(a,b){
  idx <- order(c(seq_along(a),seq_along(b)))
  unlist(c(a,b)[idx])
}

# roof factor, convert ground snow load to roof snow load
# this is the parameter for vancouver 
r_flat_mean <- 0.6
r_flat_COV <- 0.45


sigmalog <- function(mu,  sigma){
  result <- sqrt(log(1 + exp(log(sigma*sigma) - 2 * log(mu))));
  result
}

rooffactor <- function(N,r_flat_mean,r_flat_COV){
  mu <- r_flat_mean
  sigma <- r_flat_COV * r_flat_mean
  input1 <- log(mu) - sigmalog(mu, sigma)*sigmalog(mu, sigma) / 2
  input2 <- sigmalog(mu, sigma)
  rlnorm(N,meanlog = input1, sdlog = input2)
}

# NS is the number of segment in each year
NS <- 10
# snow generating function
snow_generating<- function(N,NS){
  # A = 4.3915  B = 0.1115 in vancouver 
  A <- 0.3222
  B <- 17.0689
  
  # p_0 is the probability of no-snow in one year (10 segment)
  # p_e = 0.15 for vancouver
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
  r <- rooffactor(length(p),r_flat_mean = r_flat_mean, r_flat_COV = r_flat_COV)
  # Final generation
  g_s[id_snow] <- r*(Bstr + 1.0 / Astr * (-log(-NS * log(1.0 - p_e + p_e * p))))
  g_s
}



# D_snow is the snow load
# D_snow has exponential duration with mean 2 weeks
# the time between occurrences of D_e is exponentially distributed with mean 7 months
# 2 weeks mean 0.03835 year

# generate the period  between each 10 segments (e.g. from April to Nov.)

T_snow_between <- rexp(N, 1/(7/12))

# generate the duration time (two weeks per segment) 
# generate the snow load of each segment in the winter
require(tidyverse)
T_snow_duration <- vector(mode = 'list',length = N) %>%
  lapply(function(x) x = rexp(NS, 1/0.03835))

load_snow <- vector(mode = 'list',length = N) %>%
  lapply(function(x) x = snow_generating(NS))


# generate the final snow laod
D_snow <- stepfun(cumsum(zip_modified(T_snow_between, T_snow_duration)),
                  zip_modified(rep(0, N+1), load_snow))

# check the result
t <- seq(from = 0, to = sum(zip_modified(T_snow_between, T_snow_duration)),
         by = 0.01)
plot(t, D_snow(t),type = 'l',main = 'snow load in vancouver')




##----------generate wind load---------------
# wind generating function
# G_mean = 0.33 G_cov = 0.39
wind_generating <- function(N,G_mean,G_cov){
  A <- 1.282/(G_cov*G_mean)
  B <- G_mean - (0.577/A) 
  p <- runif(N)
  wind <- B + 1/A*(-log(-log(p)))
  wind
}


snow_generating<- function(N,NS,G_mean,G_cov){
  A <- 1.282/(G_cov*G_mean)
  B <- G_mean - (0.577/A) 
  
  # p_0 is the probability of no-snow in one year (10 segment)

  p0 <- exp(-exp(A*B))
  
  # p_e is the probability of snow in a segment
  p_e <- 1.00 - exp(-1.00/NS*(exp(A*B)))
  
  # parameters for distribution W_50
  Wind_50 <- B + 1.0 / A * (-log(-NS * log(1.0 - p_e + p_e * 49/50)))
  
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
  r <- rooffactor(length(p),r_flat_mean = r_flat_mean, r_flat_COV = r_flat_COV)
  # Final generation
  g_s[id_snow] <- r*(Bstr + 1.0 / Astr * (-log(-NS * log(1.0 - p_e + p_e * p))))
  g_s
}

# D_wind is the wind load
# D_e has exponential duration with mean 1 day
# the time between occurrences of D_e is exponentially distributed with mean one year
# 2 weeks mean 0.03835 year
T_wind_e <- rexp(N, 1)
T_wind_p <- rexp(N, 1/(1/365))
load_wind <- wind_generating(N,G_mean = 0.33, G_cov = 0.39)
D_wind <- stepfun(cumsum(zip(T_wind_e, T_wind_p)), zip(rep(0, N+1), load_wind))


# check the results
t <- seq(from = 0, to = sum(T_wind_e)+sum(T_wind_p), by = 0.0001)
plot(t, D_wind(t),type = 'l',main = "wind load")
D_wind


N <- 100

NS <- 4
G_cov <- 0.150
G_mean <- (1+3.050*G_cov)/(1+2.592*G_cov)
A <- 1.282/(G_cov*G_mean)
B <- G_mean - (0.577/A) 

# p_0 is the probability of no-snow in one year (10 segment)

p0 <- exp(-exp(A*B))

# p_e is the probability of snow in a segment
p_e <- 1.00 - exp(-1.00/NS*(exp(A*B)))

# parameters for distribution W_50
Wind_velocity_50 <- B + 1.0 / A * (-log(-NS * log(1.0 - p_e + p_e * 49/50)))

# generate N random number to decide whether it snows in one segment
r_n <- runif(N+1) 
# generated snow laod
g_s <- rep(0,N+1)
# id that it snows in that segment. with r_n < p_e
# it does not snow in that winter segment if r_n > p_e
id_snow <- which (r_n <= p_e)
# p is used to generate snow load
p <- runif(length(id_snow))
# Final generation
g_s[id_snow] <- (B + 1.0 / A * (-log(-NS * log(1.0 - p_e + p_e * p))))^2*rooffactor(100,0.68,0.22)/
  ((Wind_velocity_50)^2*rooffactor(100,0.68,0.22))
g_s
mean(g_s)

rooffactor(100,0.68,0.22)/rooffactor(100,0.68,0.22)



T_wind <- rexp(N, 1/0.03835)


length(T_wind)
length(g_s)
D_wind <- stepfun(cumsum(T_wind),g_s)


t <- seq(from = 0, to = max(cumsum(T_wind)),
         by = 0.01)
plot(t, D_wind(t),type = 'l',main = 'wind_load')

