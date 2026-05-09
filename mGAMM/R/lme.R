##################################################
## Bayesian estimate of the linear regression model parameters.
##################################################
library(rstan)
library(loo)
library(tidyverse)

##################################################

setwd("~/Documents/librerias-R/mGAMM/")
load("data/data_lme.rda")
attach(data)
pairs(data)

##################################################
## linear regression model with repeated measures without constraints:

data.lme.non <- list( Y=y , n=length(y),
                      X=cbind(1,x1,x2), p=3,
                      Z=cbind(1,time), q=2,
                      N=n_distinct(id) , id=id )
para.lme.non <- c("betas","tau","sig")

fit.lme.non <- stan("stan/lme_non.stan",
                    data=data.lme.non,
                    chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

##################################################
## summary, plots and criteria

print(fit.lme.non, pars=para.lme.non)
stan_trace(fit.lme.non, pars=para.lme.non)
stan_dens(fit.lme.non, pars=para.lme.non)

estima.lme.non = get_posterior_mean(fit.lme.non, para.lme.non)[,"mean-all chains"]
plot(time,mu[,"mean-all chains"])
mu = get_posterior_mean(fit.lme.non,"mu")

ERRORS <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}
ERRORS(y,mu[,"mean-all chains"])

loo_lme_non = fit.lme.non
log_lik_lme_non = extract_log_lik(loo_lme_non, merge_chains = F)
r_eff_lme_non = relative_eff(log_lik_lme_non)
loo(log_lik_lme_non, r_eff=r_eff_lme_non)
waic(log_lik_lme_non, r_eff=r_eff_lme_non)


##################################################
## linear regression model with repeated measures with increasing constraints

data.lme.incr <- list(Y=y, n=length(y),
                      X=cbind(1,x2), p=2,
                      XI=cbind(x1), pI=1,
                      Z=cbind(rep(1,length(y))), q=1,
                      ZI=cbind(time), qI=1,
                      N=n_distinct(id), id=id )
para.lme.incr <- c("betas","betasI","tau","sig","sigI")

fit.lme.incr <- stan("stan/lme_incr.stan",
                     data=data.lme.incr,
                     chains=3,warmup=1000,iter=2000,thin=2,cores=4 )

##################################################
## summary, plots and criteria

print(fit.lme.incr, pars=para.lme.non)
stan_trace(fit.lme.incr, pars=para.lme.non)
stan_dens(fit.lme.incr, pars=para.lme.non)

mu = get_posterior_mean(fit.lme.incr,"mu")
estima.lme.incr = get_posterior_mean(fit.lme.incr, para.lme.incr)[,"mean-all chains"]
plot(time,mu[,"mean-all chains"])

ERRORS <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}
ERRORS(y,mu[,"mean-all chains"])

loo_lme_incr = fit.lme.incr
log_lik_lme_incr = extract_log_lik(loo_lme_incr, merge_chains = F)
r_eff_lme_incr = relative_eff(log_lik_lme_incr)
loo(log_lik_lme_incr, r_eff=r_eff_lme_incr)
waic(log_lik_lme_incr, r_eff=r_eff_lme_incr)

##################################################
