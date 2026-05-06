##################################################
## Bayesian estimate of the linear regression model parameters.
##################################################
library(rstan)
library(loo)
##################################################

N <- 100
set.seed(123)

## Parameters
tau <- 1.0   # sd
beta   <- c(1,2,3)   # regression coefficients

## explanatory variables
x1 <- rnorm(N)
x2 <- rbinom(N,size=1,prob=0.5)
X = cbind(1,x1,x2)   # design matrix

## Random error
epsilon = rnorm(N,mean=0,sd=tau)
eta = X%*%beta
Y = eta + epsilon
Y = as.vector(Y)

##################################################

data.lm.non <- list( Y=Y , n=length(Y),
                     X=X , p=ncol(X) )
  
para.lm.non <- c("betas","tau")
  
fit.lm.non <- stan("R/lm_non.stan",
                   data=data.lm.non,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

ERRORS <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}
##################################################
## summary, plotss and criteria

print(fit.lm.non, pars=para.lm.non)
stan_trace(fit.lm.non, pars=para.lm.non)
stan_dens(fit.lm.non, pars=para.lm.non)

mu = get_posterior_mean(fit.lm.non,"mu")
ERRORS(Y,mu[,"mean-all chains"])

estima.lm.non = get_posterior_mean(fit.lm.non, para.lm.non)[,"mean-all chains"]

plot(eta,mu[,"mean-all chains"])
plot(Y,mu[,"mean-all chains"])

loo_lm_non = fit.lm.non
log_lik_lm_non = extract_log_lik(loo_lm_non, merge_chains = F)
r_eff_lm_non = relative_eff(log_lik_lm_non)
loo(log_lik_lm_non, r_eff=r_eff_lm_non)
waic(log_lik_lm_non, r_eff=r_eff_lm_non)

##################################################