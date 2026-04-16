#BIC approximation of BARS in surface fitting with tensor product spline

library(splines)
library(EBARS)

# sample (k,xi) through reversible jump mcmc
Birth <- function(k,c,lam){
  # compute b_k
  # k: the number of splines' knots
  return(c*min(1,dpois(k+1,lam)/dpois(k,lam)))
}
Death <- function(k,c,lam){
  # compute d_k
  return(c*min(1,dpois(k-1,lam)/dpois(k,lam)))
}
Relocate <- function(k,c,lam){
  # compute r_k
  return(1-Birth(k,c,lam)-Death(k,c,lam))
}

# predict z values in (x,y) through SBARS estimation
binbars_fit <- function(x,y,conf){
  # conf: return value of SBARS function
  # value: \hat{f}(x,y)
  n = length(x)
  t_x = (x-conf$range$range_1[1])/(conf$range$range_1[2]-conf$range$range_1[1])
  t_y = (y-conf$range$range_2[1])/(conf$range$range_2[2]-conf$range$range_2[1])
  B = tensor_spline(cbind(t_x,t_y),conf$knots$xi_1,conf$knots$xi_2,intercept_1=conf$intercepts[1],intercept_2=conf$intercepts[2])
  z = B%*%conf$beta
  return(z)
}

# Bayes adaptive regression spline to surface fitting
BinBARS <- function(x,y,z,c=0.3,mu=50,burns=100,effes=100,lams=c(1,1),intercepts=c(TRUE,TRUE)){
  # x,y,z:data set (x_i,y_i;z_i)_{1<=i<=n}
  # value: conf, a list of n objects, each object is a list(nums,knots,beta_mle)
  # and loop times, accept prob, marginal like ratio and proposal ratio

  z = matrix(z,ncol=1)
  nsamp = length(x)
  x.min = min(x); x.max = max(x)
  y.min = min(y); y.max = max(y)
  t_x = (x-x.min)/(x.max-x.min)
  t_y = (y-y.min)/(y.max-y.min)

  # initialize
  k_x = 1;xi_x = runif(k_x)
  k_y = 1;xi_y = runif(k_y)
  n=0 # v records loops times
  B_hat = tensor_spline(cbind(t_x,t_y),xi_x,xi_y,intercept_1=intercepts[1],intercept_2=intercepts[2])
  BTB = t(B_hat)%*%B_hat
  beta_hat = chol2inv(chol(BTB+0.01*sum(diag(BTB))/nrow(BTB)
                           *diag(nrow(BTB))))%*%t(B_hat)%*%z
  rss_sd = norm(z-B_hat%*%beta_hat,type='f')

  while(n<burns+effes){
    J = runif(1,0,1) #J decides to renew which direction
    if(J<0.5) {
      k = k_x; xi = xi_x; lam = lams[1]
    } else {
      k = k_y; xi = xi_y; lam = lams[2]
    }

    sche = runif(1,0,1)
    birth = Birth(k,c,lam)
    death = Death(k,c,lam)
    i = sample(1:k,1)

    if(k==1) {death = 0}

    # birth scheme
    if(sche<birth){
      xi_can = (rbeta(1,xi[i]*mu,(1-xi[i])*mu)-0.5)*(1-1e-3)+0.5
      xi_next = sort(c(xi,xi_can))
      k_next = k + 1
      # compute proposal density ratio
      mix = 0
      for(j in 1:k){
        mix = mix + dbeta(xi_can, xi[j]*mu, (1-xi[j])*mu)
      }
      q_ratio = (Death(k_next,c,lam)/(k+1))/(Birth(k,c,lam)/k*mix)
    }

    # death scheme
    else if(sche<birth+death){
      xi_del = xi[i]
      xi_next = xi[-i]
      k_next = k - 1
      # compute proposal density ratio
      mix = 0
      for(j in 1:k_next){
        mix = mix + dbeta(xi_del,xi_next[j]*mu,(1-xi_next[j])*mu)
      }
      q_ratio = (Birth(k_next,c,lam)/k_next*mix)/(Death(k,c,lam)/k)
    }

    # relocation scheme
    else {
      xi_chan = xi[i]
      xi_can = (rbeta(1,xi_chan*mu,(1-xi_chan)*mu)-0.5)*(1-1e-3)+0.5
      xi_next = sort(c(xi[-i],xi_can))
      k_next = k
      # compute proposal density ratio
      q_ratio = dbeta(xi_chan,xi_can*mu,(1-xi_can)*mu)/
        dbeta(xi_can,xi_chan*mu,(1-xi_chan)*mu)
    }

    #compute accept probability to decide whether to move
    if(J<0.5) {B_hat_next = tensor_spline(cbind(t_x,t_y),xi_next,xi_y,intercept_1=intercepts[1],intercept_2=intercepts[2])}
    else {B_hat_next = tensor_spline(cbind(t_x,t_y),xi_x,xi_next,intercept_1=intercepts[1],intercept_2=intercepts[2])}

    BTB_next = t(B_hat_next)%*%B_hat_next
    beta_hat_next = chol2inv(chol(BTB_next+0.01*sum(diag(BTB_next))/nrow(BTB_next)
                                  *diag(nrow(BTB_next))))%*%t(B_hat_next)%*%z

    #margin_like_ratio = n^((k-k')/2)*(RSS_k/RSS_k')^(n/2)
    rss_sd_next = norm(z-B_hat_next%*%beta_hat_next,type='f')
    like_ratio = nsamp^((ncol(B_hat)-ncol(B_hat_next))/2)*((rss_sd/
                                            (rss_sd_next+1e-7))^nsamp)
    #acc_ratio = min(1,A) with A = q_ratio*margin_like_ratio*lambda^(k'-k)
    acc_ratio = q_ratio*like_ratio*lam^(k_next-k)
    acc = runif(1,0,1)

    #decide whether to move or not
    if(acc<acc_ratio){
      n = n + 1
      if(J<0.5) {k_x = k; xi_x = xi_next}
      else {k_y = k; xi_y = xi_next}
      B_hat = B_hat_next;beta_hat = beta_hat_next;rss_sd = rss_sd_next
    }
  }

  results = list(knots=list(xi_1=xi_x,xi_2=xi_y),beta=beta_hat,
                 range=list(range_1=c(x.min,x.max),range_2=c(y.min,y.max)),intercepts=intercepts)
  return(results)
}
