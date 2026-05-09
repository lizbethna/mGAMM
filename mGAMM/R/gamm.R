##################################################
## Bayesian estimate of the generalized additive model.
##################################################
library(rstan)
library(loo)
library(splines2)

##################################################

setwd("~/Documents/librerias-R/mGAMM/")
load("data/data_gamm.rda")
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
data.gamm.non <- list(Y=y, n=length(y) ,
                      X=cbind(rep(1,length(y))), p=1,
                      Z=cbind(1,time), q=2,
                      Xspl=list(XI1,XI2), k1=k1, S1=S1, padd=2,
                      zero = rep(0,1+k1),
                      id=id, N=n_distinct(id))

para.gamm.non <- c("betas","gamas","tau","sig","rho","lambda")

fit.gamm.non <- stan("stan/gamm_non.stan",
                   data=data.gamm.non,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

##################################################
## summary, plots and criteria

print(fit.gamm.non, pars=para.gamm.non)
stan_trace(fit.gamm.non, pars=para.gamm.non)
stan_dens(fit.gamm.non, pars=para.gamm.non)

mu = get_posterior_mean(fit.gamm.non,"mu")
estima.gamm.non = get_posterior_mean(fit.gamm.non, para.gamm.non)[,"mean-all chains"]
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

loo_gamm_non = fit.gamm.non
log_lik_gamm_non = extract_log_lik(loo_gamm_non, merge_chains = F)
r_eff_gamm_non = relative_eff(log_lik_gamm_non)
loo(log_lik_gamm_non, r_eff=r_eff_gamm_non)
waic(log_lik_gamm_non, r_eff=r_eff_gamm_non)

##################################################

##################################################
## generalized additive models with increasing contrains

data.gamm.incr <- list(Y=y, n=length(y) ,
                       X=cbind(rep(1,length(y))), p=1,
                       XI=cbind(rep(0,length(y))), pI=1,
                       Z=cbind(1,time), q=2,
                       ZI=cbind(rep(0,length(y))), qI=1,
                       Xspl=list(XI2), k1=k1, S1=S1, padd=1,
                       XIspl=list(XI1), kI1=k1, S1I=S1, pIadd=1,
                       zero = rep(0,1+k1+k1),
                       id=id, N=n_distinct(id))
para.gamm.incr <- c("betas","betasI","gamas","gamasI","tau","sig","rho","rhoI","lambda","lambdaI")

fit.gamm.incr <- stan("stan/gamm_incr.stan",
                    data=data.gamm.incr,
                    chains=3,warmup=1000,iter=2000,thin=2,cores=4 )

##################################################
## summary, plots and criteria

print(fit.gamm.incr, pars=para.gamm.incr)
stan_trace(fit.gamm.incr, pars=para.gamm.incr)
stan_dens(fit.gamm.incr, pars=para.gamm.incr)

mu = get_posterior_mean(fit.gamm.incr,"mu")
estima.gamm.incr = get_posterior_mean(fit.gamm.incr, para.gamm.incr)[,"mean-all chains"]
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

loo_gamm_incr = fit.gamm.incr
log_lik_gamm_incr = extract_log_lik(loo_gamm_incr, merge_chains = F)
r_eff_gamm_incr = relative_eff(log_lik_gamm_incr)
loo(log_lik_gamm_incr, r_eff=r_eff_gamm_incr)
waic(log_lik_gamm_incr, r_eff=r_eff_gamm_incr)

##################################################
