library(deSolve)

# Number of stochastic load profiles to generate
nrep <- 100000
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
NS <- 1



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
#  g_s[id_snow] <- Wind_velocity*Transfactor(length(p),0.68,0.22)/
#    ((Wind_velocity_50))
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
  
  # The occurrence between wind loads is 29 days exp(1/0.07945)
  T_e[[i]] <- rexp(N, 1/1)
  # The duration time of each wind load segment is 1 day exp(1/0.00274)
  T_p[[i]] <- rexp(N, 1/0.01918)
  load_p[[i]] <- wind_generating(N,NS,A,B)
  
  D_e[[i]] <- stepfun(cumsum(zip(T_e[[i]], T_p[[i]])), zip(rep(0, N+1), load_p[[i]]))
  
}
save(D_d, D_e, gamma, R_0, alpha_d, alpha_l, file="rel_loads_wind.rda")


