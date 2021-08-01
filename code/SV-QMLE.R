# Using Kalman Filter to estimate Stochastic volatility model
# Load the data ( times series of th observations )

rm(list = ls())


xrates = read_excel("xrates.xlsx")
USXUK = xrates$USXUK
obs = length(USXUK)

y = 100 * diff(log(USXUK)) - mean(diff(log(USXUK)))

r = log(y^2) 


kalmanF = function(par,r){
  #Initial values
  alpha = par[1]
  beta = par[2]
  sigma = par[3]
  
  h_predict = p_predict = h_predict = p_predict = vt = f = h_update = p_update = K = c()
  
  h_predict[1] = alpha/(1-beta)
  p_predict[1] = sigma^2/(1-beta^2)
  
  for (i in 1:length(r)) {
    
    vt[i] = r[i] - h_predict[i] + 1.27 #Prediction error
    f[i] = p_predict[i] + pi^2/2 #Kalman error variance
    
    K[i] = p_predict[i]/f[i] #Kalman Gain
    
    h_update[i] = h_predict[i] + K[i]*vt[i] #Filter update (posteriori)
    p_update[i] = p_predict[i] *(1-K[i]) #Variance of update error (posteriori)
    
    h_predict[i+1] = alpha + beta*h_update[i] #Filter prediction (priori)
    p_predict[i+1] = beta^2*p_update[i] + sigma^2 #Variance of prediction error (priori)
    
  }
  output = list("h_update" = h_update, "p_update" = p_update, "h_predict" = h_predict,
                "p_predict" = p_predict, "vt" = vt, "f" = f, "r" = r, "K" =K)
  return(output)
}

loglik <- function(par,r){
  klfil = kalmanF(par, r)
  beta <- par[1:2]
  f <- klfil$f
  y <- r
  eps <- klfil$vt
  ll = -0.5*length(y)*log(2*pi)-0.5*sum(log(f))-0.5*sum((eps)^2/f)
  return(-ll)
}

par = c(alpha = 0.1, beta = 0.9, sigma = 0.15)
loglik(par,r)

p0 = c(mu = -9, phi = 0.90 , sigmaeta = 0.3) #Initial
m <- optim(par = p0,fn = loglik,method="L-BFGS-B",
           lower =c(-Inf , -0.9999 , 0.001 ),
           upper =c(Inf , 0.9999 , Inf ),
           hessian = TRUE,
           r=r)

m$par

KL = kalmanF(m$par, r)

library(dygraphs)
mcmc = svsample(y)
mcmc_svol = colMeans(mcmc$latent)
dygraph(ts(data.frame( KL$h_predict[2:length(KL$h_predict)],KL$h_update,   mcmc_svol )))
dygraph(ts(data.frame( "Kalman filter estimate"=KL$h_update
                       ,"Stochvol MCMC estimation" =  mcmc_svol  )))
cor(ts(data.frame( "Kalman filter estimate"=KL$h_update,"Stochvol MCMC estimation" =  mcmc_svol  )))

colMeans(mcmc$para)
c(m$par[1]/(1-m$par[2]), m$par[2],m$par[3])

#Smooth
h = c()
p = c = c()
h[length(r)] = KL$h_update[length(r)]
p[length(r)] = KL$p_update[length(r)]

for(i in (length(r)-1):1){
  c = m$par[2]/KL$p_update[i+1]
  h[i] = KL$h_update[i] + KL$p_update[i]* c*(h[i+1] - KL$h_predict[i+1])
  p[i] = KL$p_update[i] + KL$p_update[i]* c^2*(p[i+1] - KL$p_predict[i+1])
}

dygraph(ts(data.frame( "Kalman filter estimate" = KL$h_update,
                       "Stochvol MCMC estimation" =  mcmc_svol,
                       "Kalman smoother estimate" = h  )))%>%
  dyLegend(width = 350)


