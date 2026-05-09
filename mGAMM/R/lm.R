##################################################
## Bayesian estimate of the linear regression model parameters.
##################################################
library(rstan)
library(loo)
##################################################

setwd("~/Documents/librerias-R/mGAMM/")
load("data/data_lm.rda")
attach(data)
pairs(data)

##################################################
## linear regression model

data.lm.non <- list( Y=y , n=length(y),
                     X=cbind(1,x1,x2) , p=3 )

para.lm.non <- c("betas","tau")

fit.lm.non <- stan("stan/lm_non.stan",
                   data=data.lm.non,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

##################################################
## summary, plots and criteria

print(fit.lm.non, pars=para.lm.non)
stan_trace(fit.lm.non, pars=para.lm.non)
stan_dens(fit.lm.non, pars=para.lm.non)

mu = get_posterior_mean(fit.lm.non,"mu")
estima.lm.non = get_posterior_mean(fit.lm.non, para.lm.non)[,"mean-all chains"]
plot(eta,mu[,"mean-all chains"])
plot(Y,mu[,"mean-all chains"])


ERRORS <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}
ERRORS(y,mu[,"mean-all chains"])

loo_lm_non = fit.lm.non
log_lik_lm_non = extract_log_lik(loo_lm_non, merge_chains = F)
r_eff_lm_non = relative_eff(log_lik_lm_non)
loo(log_lik_lm_non, r_eff=r_eff_lm_non)
waic(log_lik_lm_non, r_eff=r_eff_lm_non)

##################################################

##################################################
## Linea regression model witg increasing  constrains

data.lm.incr <- list(Y=Yerror, n=length(Ytrue),
                     X=cbind(1,x2), p=2,
                     XI=cbind(x1), pI=1)
para.lm.incr <- c("betas","betasI","tau")

fit.lm.incr <- stan("stan/lm_incr.stan",
                    data=data.lm.incr,
                    chains=3,warmup=1000,iter=2000,thin=2,cores=4 )

##################################################
## summary, plots and criteria

print(fit.lm.incr, pars=para.lm.incr)
stan_trace(fit.lm.incr, pars=para.lm.incr)
stan_dens(fit.lm.incr, pars=para.lm.incr)

mu = get_posterior_mean(fit.lm.incr,"mu")
estima.lm.incr = get_posterior_mean(fit.lm.incr, para.lm.incr)[,"mean-all chains"]
plot(eta,mu[,"mean-all chains"])
plot(Ytrue,mu[,"mean-all chains"])

ERRORS <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}
ERRORS(Ytrue,mu[,"mean-all chains"])

loo_lm_incr = fit.lm.incr
log_lik_lm_incr = extract_log_lik(loo_lm_incr, merge_chains = F)
r_eff_lm_incr = relative_eff(log_lik_lm_incr)
loo(log_lik_lm_incr, r_eff=r_eff_lm_incr)
waic(log_lik_lm_incr, r_eff=r_eff_lm_incr)

##################################################
