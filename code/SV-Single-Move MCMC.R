# # =======================================================================
# # Univariate stationary AR(1) stochastic volatility
# # A single-move sampler for the volatility based on Kim et al. (1998) algorithm
# #
# # y_t = exp(h_t/2)*epsilon_t,                    epsilon_t ~ N(0,1)),
# # h_t = mu + phih(h_{t-1}-mu) + sig_nu*nu_t,     nu_t ~ N(0,1),
# # h_1 = ~ N(mu, sig_nu2/(1-phi^2) ), 
# # 
# # See Kim, Shephard, and Chib (1998) 
# # (c) 2021, Hamza Bouguerra . Email: hamza93bouguerra@gmail.com
# # =======================================================================



rm(list=ls())
library(truncnorm)  # To draw proposals from truncated normal distribution
library(invgamma)   # To draw from inverse-gammal distribution
library(stochvol)   # Efficient Bayesian Inference for Stochastic Volatility (SV) by Kastner G (2016a)
require(svMisc)     # Fancier text progress
library(readxl)     # To read excel files
library(dygraphs)   # Good package to represent time series of volatility
source("SV-Single-Move MCMC-Functions.R")
options(scipen = 999,digits = 5)  


# I use the same data as the reference paper of Kim, Shephard, and Chib (1998) 

# The daily observations of weekday close exchange rates for
# the U.K. Sterling £/U.S. Dollar $ exchange rate from 1/10/81 to 28/6/85 with 
# sample size = 946

xrates = read_excel("xrates.xlsx")
USXUK = xrates$USXUK
obs = length(USXUK)

y = 100 * diff(log(USXUK)) - mean(diff(log(USXUK)))


t = length(y)
draws = 15000
burnin = 0.2*draws

#Initializing parameters
mu = 0
phi = 0.95
tho_2 = 0.02

#Initializing states
h = rep(0,t)

#To store outputs
h_store = matrix(0, draws-burnin, t)
tho_2_store = rep(NULL, draws-burnin)
phi_store = rep(NULL, draws-burnin)
mu_store = rep(NULL, draws-burnin)

# Gibbs sampling algorithm for SV model
for (i in 1:draws) {
  progress(i*(100/draws))
  for (j in 1:t) {
    #Sample the volatility
    h[j] = Sample_ht(mu, phi, h, tho_2, y, j)
  }
  
  #Sample the parameters
  tho_2 = Sample_tho_2( phi, mu, h)
  phi = Sample_phi(mu, phi, tho_2, h, t)
  mu = Sample_mu(phi, tho_2, h, t)
  
  if(i > burnin){
    i = i-burnin
    #Store the parameters
    h_store[i,] = h
    tho_2_store[i] = tho_2
    phi_store[i] = phi
    mu_store[i] = mu
  }
  
  
  if (i == draws)
    cat("Done!\n")
  
}



#Display trace plots for the parameters with their densities
par(mfrow=c(3,3), mar=c(2,1,2,1)) # 3 rows and 1 columns
plot(mu_store, type = 'l', main = "Trace of mu (thin = 1)" , ylab = "")
plot(density(mu_store),  main = "Density of mu", ylab = "")
acf(mu_store,lag.max = 100)
plot(phi_store, type = 'l', main = "Trace of phi (thin = 1)", ylab = "")
plot(density(phi_store),  main = "Density of phi", ylab = "")
acf(phi_store,lag.max = 100)
plot(tho_2_store, type = 'l', main = "Trace of tho_2 (thin = 1)", ylab = "")
plot(density(tho_2_store),  main = "Density of tho_2", ylab = "")
acf(tho_2_store,lag.max = 100)



# Comparaison with stochvol results
# stochol = svsample(y)

# stochol_sv = colMeans(stochol$latent)

load("C:/Users/Lenovo/Desktop/SV Test/svdata.RData")

dygraph(ts(data.frame( colMeans(h_store),stochol_sv ) ))
cor(data.frame( colMeans(h_store),stochol_sv ) )



c(mean( mu_store),
  mean(phi_store),
  mean(sqrt(tho_2_store)))

colMeans(stochol$para)


