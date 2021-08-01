# # =======================================================================
# # univariate stationary AR(1) stochastic volatility
# # See Kim, Shephard, and Chib 1998) 
# # (c) 2021, Bouguerra Hamza . Email: hamza93bouguerra@gmail.com
# # =======================================================================


#Full conditional distribution of h_t:

log_f_star = function(proposal_h,y){
  out = -0.5*proposal_h - 0.5*y^2 * exp(-proposal_h)
  return(out)
}

log_g_star = function(proposal_h,hstar,y){
  out = -0.5*proposal_h-0.5*y^2*(exp(-hstar)*(1+hstar)-proposal_h*exp(-hstar))
  return(out)
}

Sample_ht = function(mu, phi, h, tho_2, y, k){
  if (k==1){
    h0 = rnorm(1,mu, sqrt(tho_2/(1-phi^2)))
    hstar = mu + phi*( (h0-mu) + (h[k+1]-mu) )/(1+phi^2)
    sig_h_2 = tho_2
  }else if (k >= 2 && k <= (t-1)){
    hstar = mu + phi*( (h[k-1]-mu) + (h[k+1]-mu) )/(1+phi^2)
    sig_h_2 = tho_2/(1+phi^2)
  } else {
    hstar = mu + phi*(h[k-1]-mu)
    sig_h_2 = tho_2
  }
  
  mu_t = hstar + (sig_h_2/2)*(y[k]^2*exp(-hstar)-1)
  
  #Accept-reject procedure to sample h_t
  repeat {
    proposal_h = rnorm(1,mu_t,sqrt(sig_h_2))
    u = runif(1)
    alpha = exp(log_f_star(proposal_h,y[k]))/exp(log_g_star(proposal_h,hstar,y[k]))
    if(alpha>u) {
      h = proposal_h
      break
    }
  }

  return(h)
  
}


#Full conditional distribution of tho_2
Sample_tho_2 = function(phi, mu, h){
  sigma_r = 5
  S_sigma = 0.01*sigma_r

  ## sample sig2
  sum1 = (t+sigma_r)/2
  tmp1 = (h[1] - mu)^2 * (1 - phi^2)
  tmp = sum( c( (h[2:length(h)]-mu)-phi*(h[1:(length(h)-1)]-mu) )^2  )

  delta1 = (S_sigma + tmp1 + tmp)/2
  tho_2 = 1/rgamma(1 , shape = sum1, rate = delta1);
  #tho_2 = rinvgamma(1 , shape = sum1, rate = delta1);
  
  return(tho_2)
  
}


#Full conditional distribution of phi

#Distribution of prior phi
p_phi = function(phi){
  a = 20
  b = 1.5
  out = ( ( (1+phi)/2  )^(a-1) ) * ( ( (1-phi)/2  )^(b-1) )
  return(out)
  
}

g_phi = function(mu, phi, tho_2, h){
  out = log(p_phi(phi)) - (h[1]-mu)^2 * (1-phi^2)/(2*tho_2) + 0.5*log(1-phi^2)
  return(out)
}


Sample_phi = function(mu, phi, tho_2, h, t){
  phi_hat = sum( (h[2:t]-mu)*(h[1:(t-1)]-mu) )/ sum( (h[1:(t-1)]-mu)^2 )
  V_phi = tho_2/sum( (h[1:(t-1)]-mu)^2 )
  
  #Metropolis-Hastings algorithm: 
  
  #Sample phi_star: the proposal
  phi_star = rtruncnorm(1,a = -0.9999 , b = 0.9999 , mean = phi_hat, sd = sqrt(V_phi))
  
  #Acceptance probability
  alp = exp(g_phi(mu, phi_star, tho_2, h) - g_phi(mu, phi, tho_2, h) );
  if (alp>runif(1)) phi = phi_star;
  
  out = phi
  return(out)
}


#Full conditional distribution of mu
Sample_mu = function(phi, tho_2, h, t){
  sig_mu_2 = tho_2*(  (t-1)*(1-phi)^2 + (1-phi^2) )^(-1)
  mu_hat = sig_mu_2* ((1-phi^2)/tho_2 * h[1] + (1-phi)/tho_2 * sum(h[2:t]-phi*h[1:(t-1)]))
  
  mu = rnorm(1,mu_hat ,sqrt(sig_mu_2) )
  out = mu
  return(out)
}



