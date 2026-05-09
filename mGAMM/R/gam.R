##################################################
## Bayesian estimate of the generalized additive model.
##################################################
library(rstan)
library(loo)
library(splines2)

##################################################

setwd("~/Documents/librerias-R/mGAMM/")
load("data/data_gam.rda")
attach(data)
pairs(data)

##################################################
## Design matrice for B-splines basis function
k1 = 6
cuantiles1 = seq(0,1,length.out=k1)[-c(1,k1)]
knots1 = quantile(x1, cuantiles1)
knots2 = quantile(x2, cuantiles1)
XI1 <- ibs(x1, knots = knots1, degree = 1, intercept = TRUE)
XI2 <- ibs(x2, knots = knots2, degree = 1, intercept = TRUE)

## Penalized matrix
diffMatrix = function(k, d = 2){
  if( (d<1) || (d %% 1 != 0) )stop("d must be a positive integer value");
  if( (k<1) || (k %% 1 != 0) )stop("k must be a positive integer value");
  if(d >= k)stop("d must be lower than k");
  out = diag(k);
  for(i in 1:d){
    out = diff(out);
  }
  return(out)
}
(D1 = diffMatrix(k=k1, d=2))
(S1 = t(D1)%*%D1 + diag(1,k1)*10e-4)


##################################################
## generalized additive models
data.gam.non <- list(Y=y, n=length(y) ,
                     X=cbind(rep(1,length(y))), p=1,
                     Xspl=list(XI1,XI2), k1=k1, S1=S1, padd=2,
                     zero = rep(0,1+k1))

para.gam.non <- c("betas","tau","gamas","lambda")

fit.gam.non <- stan("stan/gam_non.stan",
                   data=data.gam.non,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

##################################################
## summary, plots and criteria

print(fit.gam.non, pars=para.gam.non)
stan_trace(fit.gam.non, pars=para.gam.non)
stan_dens(fit.gam.non, pars=para.gam.non)

mu = get_posterior_mean(fit.gam.non,"mu")
estima.gam.non = get_posterior_mean(fit.gam.non, para.gam.non)[,"mean-all chains"]
plot(y,mu[,"mean-all chains"])


ERRORS <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}
ERRORS(y,mu[,"mean-all chains"])

loo_gam_non = fit.gam.non
log_lik_gam_non = extract_log_lik(loo_gam_non, merge_chains = F)
r_eff_gam_non = relative_eff(log_lik_gam_non)
loo(log_lik_gam_non, r_eff=r_eff_gam_non)
waic(log_lik_gam_non, r_eff=r_eff_gam_non)

##################################################

##################################################
## generalized additive models with increasing contrains

data.gam.incr <- list(Y=y, n=length(y) ,
                      X=cbind(rep(0,length(y))), p=1,
                      XI=cbind(rep(0,length(y))), pI=1,
                      Xspl=list(XI1), k1=k1, S1=S1, padd=1,
                      XIspl=list(XI2), k1I=k1, S1I=S1, pIadd=1,
                      zero = rep(0,1+k1+k1))
para.gam.incr <- c("betas","betasI","tau","gamas","lambda","gamasI","lambdaI")

fit.gam.incr <- stan("stan/gam_incr.stan",
                    data=data.gam.incr,
                    chains=3,warmup=1000,iter=2000,thin=2,cores=4 )

##################################################
## summary, plots and criteria

print(fit.gam.incr, pars=para.gam.incr)
stan_trace(fit.gam.incr, pars=para.gam.incr)
stan_dens(fit.gam.incr, pars=para.gam.incr)

mu = get_posterior_mean(fit.gam.incr,"mu")
estima.gam.incr = get_posterior_mean(fit.gam.incr, para.gam.incr)[,"mean-all chains"]
plot(y,mu[,"mean-all chains"])

ERRORS <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}
ERRORS(y,mu[,"mean-all chains"])

loo_gam_incr = fit.gam.incr
log_lik_gam_incr = extract_log_lik(loo_gam_incr, merge_chains = F)
r_eff_gam_incr = relative_eff(log_lik_gam_incr)
loo(log_lik_gam_incr, r_eff=r_eff_gam_incr)
waic(log_lik_gam_incr, r_eff=r_eff_gam_incr)

##################################################
