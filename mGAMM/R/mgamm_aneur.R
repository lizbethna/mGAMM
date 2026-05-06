#' @title Abdominal aortic aneurysm progression data
#' @description We consider a data set that corresponds to longitudinal
#' measurements of grades of abdominal aortic aneurysms, measured by
#' ultrasound examination of the diameter of the aorta.
#' @param
#' @return
#' @export
#' @examples
#'

# libraries
library(VGAM)
##### library(mgcv)
library(rstan)
##### library(rjags)
##### load.module("glm")
##### library(plot3D)
library(dplyr)
library(ggplot2)
##### library(corrplot)
library(splines)
library(splines2)
##### library(cgam)
##### library(MASS)
library(msm)
##### library(nlme)
#setwd("/Users/lizbethna/Documents/librerias-R/mgam/R/")
#dir <- "/Users/lizbethna/Documents/librerias-R/mgam/R/"
library(loo)


# Data
data(aneur)
help(aneur)
attach(aneur)

N = n_distinct(aneur$ptnum) # number of subjects
K = max(table(aneur$ptnum)) # number of times
table(table(aneur$ptnum))

J = n_distinct(aneur$state)   # categories for ordinal response variable
Y_diam = array(NA,dim=c(N,K))
Y_state = array(NA,dim=c(N,K))
X_age = array(NA,dim=c(N,K))
Ki = table(aneur$ptnum) # number of observations per subject

Ni = c(0,cumsum(Ki))
for(i in 1:N){
  aneur_i = aneur[aneur$ptnum==i,]
  for(k in 1:Ki[i]){
    Y_diam[i,k] = aneur_i$diam[k]
    Y_state[i,k] = aneur_i$state[k]
    X_age[i,k] = aneur_i$age[k]
  }
}

# Considering only data having more than one screen (diam!=29, or diam<29 & dim>29)
idx3 = c()
for(i in 1:N){
  if( min(Y_diam[i,1:Ki[i]])!=max(Y_diam[i,1:Ki[i]])){
    idx3 = c(idx3,i)
  }
}
Y3_diam = Y_diam[idx3,]
Y3_state = Y_state[idx3,]
X3_age = X_age[idx3,]
N3 = length(idx3)
Ki3 = Ki[idx3]
aneur3 = aneur%>%filter(aneur$ptnum%in%idx3)


ggplot(data=aneur3,
       mapping=aes(x=age,y=diam,group=ptnum)) +
  geom_line(color=aneur3$ptnum) +
  theme_bw() +
  xlab("Age at examination in years") + ylab("Aortic diameter in mm") +
  ggtitle("Profiles aortic diameter by patient")

ggplot(data=aneur3,
       mapping=aes(x=age,y=state,group=ptnum)) +
  geom_line(color=aneur3$ptnum) +
  theme_bw() +
  xlab("Age at examination in years") + ylab("States of aneurysm severity") +
  ggtitle("Profiles states of aneurysm severity by patient")


# First analysis: continuous response variable

# Centred variables

y = aneur3$diam -29
x1 = aneur3$age -60
x2 = aneur3$age -60
id = as.numeric(as.factor(aneur3$ptnum))

n = length(y)
N = n_distinct(id)
Ni = c(0,cumsum(table(id)))


# The number of knots K is considered a priori
k1 = 5 #knots
k2 = 5 #knots
knots1 = quantile(x1, seq(0,1,length.out=k1)[-c(1,k1)])
knots2 = quantile(x2, seq(0,1,length.out=k2)[-c(1,k2)])


# The number of knots is chosen large enough to avoid over-smoothing,
# but small enough to avoid excessive computational cost.


# Generate the design matrix X for the B-splines

# Generate a basis matrix for Natural Cubic Splines
X2 <- ns(x = x2, knots = knots2, intercept = TRUE)
matplot(x2, X2)


# Generate the design matrix XI1 for the I-splines

# ibs: integrated basis splines
# degree = 3 cubic splines
XI1 <- ibs(x1, knots = knots1, degree = 1, intercept = TRUE)
matplot(x1, XI1)
abline(v = knots1, h = knots1, lty = 2, col = "gray")

# 4. Define the penalization S1 and S2

#Este es el código que produce la matriz de diferenciación.
#No es el óptimo, pero funciona.
#“k” es el número de b-splines y
#“d” el orden de la diferenciación.
#Adjunto el artículo donde discutimos esto (página 7).

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
(D2 = diffMatrix(k=k2, d=2))
(S1 = t(D1)%*%D1 + diag(1,k1)*10e-4)
(S2 = t(D2)%*%D2 + diag(1,k2)*10e-4)



# Linear regression model
# Non constraints
#data.lm.non <- list( Y=aneur$diam, n=length(aneur$age),
#                     X=cbind(1,aneur$age), p=2)
data.lm.non <- list( Y=y, n=length(y),
                     X=cbind(1,x1), p=2)
param.lm = c("betas","tau")
fit.lm <- stan("R/lm_non.stan",
               data=data.lm.non,
               chains=3,warmup=500,iter=1000,thin=2,cores=6 )
print(fit.lm, pars=param.lm, digits=5)
stan_trace(fit.lm,pars=param.lm)
stan_dens(fit.lm,pars=param.lm)
pairs(fit.lm, pars = c("betas"), las = 1)
mu = get_posterior_mean(fit.lm,"mu")
plot(x1,mu[,"mean-all chains"])


# Linear regression model
# Increasing constraints
#data.lm.incr <- list( Y=aneur$diam, n=length(aneur$diam),
#                      X=cbind(rep(1,length(aneur$diam))), p=1,
#                      XI=cbind(aneur$age), pI=1)
data.lm.incr <- list( Y=y, n=length(y),
                      X=cbind(rep(1,length(y))), p=1,
                      XI=cbind(x1), pI=1)
param.lm.incr = c("betas","betasI","tau")
fit.lm.incr <- stan("R/lm_incr.stan",
               data=data.lm.incr,
               chains=3,warmup=500,iter=1000,thin=1,cores=6 )

print(fit.lm.incr, pars=param.lm.incr, digits=5)
stan_trace(fit.lm.incr,pars=param.lm.incr)
stan_dens(fit.lm.incr,pars=param.lm.incr)
pairs(fit.lm.incr, pars = c("betas","betasI"), las = 1)
mu = get_posterior_mean(fit.lm.incr,"mu")
plot(x1,mu[,"mean-all chains"])



# Linear Mixed Effects
# Non constraints
# data.lme.non <- list( Y=aneur$diam, n=length(aneur$diam),
#                       X=cbind(1,aneur$age), p=2,
#                       Z=cbind(aneur$age), q=1,
#                       Ni=c(0,cumsum(table(aneur$ptnum))), N=n_distinct(aneur$ptnum))
data.lme.non <- list( Y=y, n=length(y),
                      X=cbind(1,x1), p=2,
                      Z=cbind(x1), q=1,
                      id=id, N=N)
fit.lme.non <- stan("R/lme_non.stan",
                    data=data.lme.non,
                    chains=3,warmup=500,iter=1000,thin=1,cores=6 )
param.lme.non = c("betas","sig","tau")
print(fit.lme.non, pars=param.lme.non, digits=5)
stan_trace(fit.lme.non,pars=param.lme.non)
stan_dens(fit.lme.non,pars=param.lme.non)
pairs(fit.lme.non, pars = c("betas"), las = 1)
mu = get_posterior_mean(fit.lme.non,"mu")
plot(x1,mu[,"mean-all chains"])

# Linear Mixed Effects
# Increasing constraints
# data.lme.incr <- list( Y=aneur$diam, n=length(aneur$diam),
#                        X=cbind(rep(1,length(aneur$diam))), p=1,
#                        XI=cbind(aneur$age), pI=1,
#                        Z=cbind(rep(0,length(aneur$diam))), q=1,
#                        ZI=cbind(aneur$age), qI=1,
#                        Ni=c(0,cumsum(table(aneur$ptnum))), N=n_distinct(aneur$ptnum))
data.lme.incr <- list( Y=y, n=length(y),
                       X=cbind(rep(1,length(y))), p=1,
                       XI=cbind(x1), pI=1,
                       Z=cbind(rep(0,length(y))), q=1,
                       ZI=cbind(x1), qI=1,
                       id=id, N=N)
fit.lme.incr <- stan("R/lme_incr.stan",
                     data=data.lme.incr,
            chains=3,warmup=1000,iter=2000,thin=2,cores=6 )
param.lme.incr = c("betas","betasI","sig","tau")
print(fit.lme.incr, pars=param.lme.incr, digits=5)
stan_trace(fit.lme.incr,pars=param.lme.incr)
stan_dens(fit.lme.incr,pars=param.lme.incr)
pairs(fit.lme.incr, pars = c("betas","betasI"), las = 1)
mu = get_posterior_mean(fit.lme.incr,"mu")
plot(x1,mu[,"mean-all chains"])



# Generalized Additive Models
# Spline NO constraints
#knots1 = quantile(aneur$age, seq(0,1,length.out=k1)[-c(1,k1)])
#XI1 <- ibs(aneur$age, knots = knots1, degree = 1, intercept = TRUE)
# data.gam.non <- list( Y=aneur$diam, n=length(aneur$diam),
#                       X=cbind(rep(1,length(aneur$diam))), p=1,
#                       Xspl=ibs(aneur$age, knots = quantile(aneur$age, seq(0,1,length.out=k1)[-c(1,k1)]), degree = 1, intercept = TRUE),
#                       k1=k1, S1=S1,
#                       zero=rep(0,1+k1))
data.gam.non <- list(Y=y, n=length(y),
                     X=cbind(rep(1,length(y))), p=1,
                     Xspl=XI1,k1=k1, S1=S1,
                     zero=rep(0,1+k1))
param.gam.non = c("betas", "invtau2","tau", "lambda","rho")
fit.gam.non <- stan("R/gam_non.stan",
                    data=data.gam.non,
                    chains=3,warmup=1000,iter=2000,thin=2,cores=6)
print(fit.gam.non, pars=param.gam.non, digits=5)
stan_trace(fit.gam.non, pars=param.gam.non)
stan_plot(fit.gam.non, pars=c("betas"), point_est = "mean", show_density = TRUE)
stan_plot(fit.gam.non, pars=c("gamas"), point_est = "mean", show_density = TRUE)
stan_plot(fit.gam.non, pars=c("invtau2","tau","lambda","rho"), point_est = "mean", show_density = TRUE)
stan_dens(fit.gam.non, pars=c("betas"))
stan_dens(fit.gam.non, pars=c("gamas"))
stan_dens(fit.gam.non, pars=c("invtau2","tau", "lambda","rho"))
pairs(fit.gam.non, pars = c("betas","gamas"), las = 1)
mu1 = get_posterior_mean(fit.gam.non,"mu1")
plot(x1,mu1[,"mean-all chains"])



# Generalized Additive Models
# Increasing constraints
#knots1 = quantile(aneur$age, seq(0,1,length.out=k1)[-c(1,k1)])
#XI1 <- ibs(aneur$age, knots = knots1, degree = 1, intercept = TRUE)
# data.gam.incr <- list(Y=aneur$diam, n=length(aneur$diam),
#                       X=cbind(rep(1,length(aneur$diam))), p=1,
#                       XI=cbind(rep(0,length(aneur$diam))), pI=1,
#                       Xspl=cbind(rep(0,length(aneur$diam))), k1=1, S1=as.matrix(1,1,1),
#                       XIspl=ibs(aneur$age, knots = quantile(aneur$age, seq(0,1,length.out=k1)[-c(1,k1)]), degree = 1, intercept = TRUE),
#                       kI1=k1, S1I=S1,
#                       zero=rep(0,1+k1+1))
data.gam.incr <- list(Y=y, n=length(y),
                      X=cbind(rep(1,length(y))), p=1,
                      XI=cbind(rep(0,length(x1)), pI=1,
                      Xspl=cbind(rep(0,length(y))), k1=1, S1=as.matrix(1,1,1),
                      XIspl=XI1, kI1=k1, S1I=S1,
                      zero=rep(0,1+k1+1))
param.gam.incr = c("betas", "invtau2","tau", "lambda","rho")
fit.gam.incr <- stan("R/gam_incr.stan",
                    data=data.gam.incr,
                    chains=3,warmup=1000,iter=2000,thin=2,cores=6)
print(fit.gam.incr, pars=param.gam.incr, digits=5)
stan_trace(fit.gam.incr, pars=param.gam.incr)
stan_plot(fit.gam.incr, pars=c("betas"), point_est = "mean", show_density = TRUE)
stan_plot(fit.gam.incr, pars=c("gamas"), point_est = "mean", show_density = TRUE)
stan_plot(fit.gam.incr, pars=c("invtau2","tau","lambda","rho"), point_est = "mean", show_density = TRUE)
stan_dens(fit.gam.incr, pars=c("betas"))
stan_dens(fit.gam.incr, pars=c("gamas"))
stan_dens(fit.gam.incr, pars=c("invtau2","tau", "lambda","rho"))
pairs(fit.gam.incr, pars = c("betas","gamas"), las = 1)
mu1 = get_posterior_mean(fit.gam.incr,"mu1")
plot(x1,mu1[,"mean-all chains"])


# Generalized Additive Models random effects
# Non constraints
# data.gamm.non <- list(Y=aneur$diam, n=length(aneur$diam),
#                       X=cbind(rep(1,length(aneur$diam))), p=1,
#                       Z=cbind(aneur$age), q=1,
#                       Xspl=ibs(aneur$age, knots = quantile(aneur$age, seq(0,1,length.out=k1)[-c(1,k1)]), degree = 1, intercept = TRUE),
#                       k1=k1, S1=S1,
#                       zero=rep(0,1+k1),
#                       Ni=c(0,cumsum(table(aneur$ptnum))), N=n_distinct(aneur$ptnum))
data.gamm.non <- list(Y=y, n=length(y) ,
                 X=cbind(rep(1,length(y))), p=1,
                 Z=cbind(x1), q=1,
                 Xspl=XI1, k1=k1, S1=S1,
                 zero = rep(0,1+k1),
                 id=id, N=N)
param.gamm.non = c("betas", "invtau2","tau", "lambda","rho", "invsig2","sig")
fit.gamm.non <- stan("R/gamm_non.stan",
            data=data.gamm.non,
            chains=3,warmup=500,iter=1000,thin=2,cores=6)
print(fit.gamm.non, pars=param.gamm.non, digits=5)
stan_trace(fit.gamm.non, pars=param.gamm.non)
stan_trace(fit.gamm.non, pars="gamas")
stan_plot(fit.gamm.non, pars=c("betas"), point_est = "mean", show_density = TRUE)
stan_plot(fit.gamm.non, pars=c("invtau2","tau", "invsig2","sig",  "lambda","rho"), point_est = "mean", show_density = TRUE)
stan_dens(fit.gamm.non, pars=c("betas"))
stan_dens(fit.gamm.non, pars=c("invtau2","tau", "invsig2","sig", "lambda","rho"))
pairs(fit.gamm.non, pars = c(c("betas","gamas")), las = 1)
mu1 = get_posterior_mean(fit.gamm.non,"mu1")
plot(x1,mu1[,"mean-all chains"])
mu2 = get_posterior_mean(fit.gamm.non,"mu2")
plot(x1,mu2[,"mean-all chains"])



# 8. Spline con restricciones creciente
## 8.2. LME: Spline con restricciones creciente
data.gamm.incr <- list(Y=y, n=length(y) ,
                       X=cbind(rep(1,length(y))), p=1,
                       XI=cbind(rep(0,length(y))), pI=1,
                       Z=cbind(rep(0,length(y))), q=1,
                       ZI=cbind(x1), qI=1,
                       Xspl=cbind(rep(0,length(y))), k1=1, S1=as.matrix(1,1,1),
                       XIspl=XI1, kI1=k1, S1I=S1,
                       zero = rep(0,1+k1+1),
                       id=id, N=N)
fit.gamm.incr <- stan("R/gamm_incr.stan",
            data=data.gamm.incr,
            chains=3,warmup=500,iter=1000,thin=2,cores=6)
param.gamm.incr = c("betas", "invtau2","tau", "lambda","rho", "invsig2","sig")
print(fit.gamm.incr, pars=param.gamm.incr, digits=5)
stan_trace(fit.gamm.incr, pars=param.gamm.incr)
stan_plot(fit.gamm.incr, pars=c("betas","gamas"), point_est = "mean", show_density = TRUE)
stan_plot(fit.gamm.incr, pars=c("invtau2","tau", "invsig2","sig",  "lambda","rho"), point_est = "mean", show_density = TRUE)
stan_dens(fit.gamm.incr, pars=c("betas","gamas"))
stan_dens(fit.gamm.incr, pars=c("invtau2","tau", "invsig2","sig", "lambda","rho"))
pairs(fit.gamm.incr, pars = c("betas","gamas"), las = 1)
mu1 = get_posterior_mean(fit.gamm.incr,"mu1")
plot(x1,mu1[,"mean-all chains"])
mu2 = get_posterior_mean(fit.gamm.incr,"mu2")
plot(x1,mu2[,"mean-all chains"])




# Generalized Additive Models random effects
# Spline Increasing constraints
data.gamm.add.non <- list(Y=y, n=length(y) ,
                      X=cbind(rep(1,length(y))), p=1,
                      Z=cbind(x1), q=1,
                      Xspl=XI1, k1=k1, S1=S1,
                      Zspl=XI1, k2=k1, S2=S1,
                      zero = rep(0,1+k1+k1),
                      id=id, N=N)
param.gamm.add.non = c("betas", "invtau2","tau", "invsig2","sig",
                       "lambda1","rho1", "lambda2","rho2")
fit.gamm.add.non <- stan("R/gamm_add_non.stan",
                     data=data.gamm.add.non,
                     chains=3,warmup=500,iter=1000,thin=2,cores=6)
print(fit.gamm.add.non, pars=param.gamm.add.non, digits=5)
stan_trace(fit.gamm.add.non, pars=param.gamm.add.non)
stan_plot(fit.gamm.add.non, pars=c(c("betas","gamas")), point_est = "mean", show_density = TRUE)
stan_plot(fit.gamm.add.non, pars=c("invtau2","tau", "invsig2","sig","lambda1","rho1","lambda2","rho2"), point_est = "mean", show_density = TRUE)
stan_dens(fit.gamm.add.non, pars=c(c("betas","gamas")))
stan_dens(fit.gamm.add.non, pars=c("invtau2","tau", "invsig2","sig","lambda1","rho1","lambda2","rho2"))
pairs(fit.gamm.add.non, pars = c(c("betas","gamas")), las = 1)
mu1 = get_posterior_mean(fit.gamm.add.non,"mu1")
plot(x1,mu1[,"mean-all chains"])
mu2 = get_posterior_mean(fit.gamm.add.non,"mu2")
plot(x1,mu2[,"mean-all chains"])
mu3 = get_posterior_mean(fit.gamm.add.non,"mu3")
plot(x1,mu3[,"mean-all chains"])
mu4 = get_posterior_mean(fit.gamm.add.non,"mu4")
plot(x1,mu4[,"mean-all chains"])



# 8. Spline con restricciones creciente
## 8.2. LME: Spline con restricciones creciente
data.gamm.add.incr <- list(Y=y, n=length(y) ,
                       X=cbind(rep(1,length(y))), p=1,
                       XI=cbind(rep(0,length(y))), pI=1,
                       Z=cbind(rep(0,length(y))), q=1,
                       ZI=cbind(x1), qI=1,
                       Xspl=cbind(rep(0,length(y))), k1=1, S1=as.matrix(1,1,1),
                       XIspl=XI1, kI1=k1, S1I=S1,
                       Zspl=cbind(rep(0,length(y))), k2=1, S2=as.matrix(1,1,1),
                       ZIspl=XI1, kI2=k1, S2I=S1,
                       zero = rep(0,1+k1+1+1+k1),
                       id=id, N=N)
fit.gamm.add.incr <- stan("R/gamm_add_incr.stan",
                      data=data.gamm.add.incr,
                      chains=3,warmup=500,iter=1000,thin=2,cores=6)
param.gamm.add.incr = c("betas", "invtau2","tau", "invsig2","sig",
                        "lambda1","rho1", "lambda2","rho2")
print(fit.gamm.add.incr, pars=param.gamm.add.incr, digits=5)
stan_trace(fit.gamm.add.incr, pars=param.gamm.add.incr)
stan_plot(fit.gamm.add.incr, pars=c("betas","gamas"), point_est = "mean", show_density = TRUE)
stan_plot(fit.gamm.add.incr, pars=c("invtau2","tau", "invsig2","sig",  "lambda1","rho1","lambda2","rho2"), point_est = "mean", show_density = TRUE)
stan_dens(fit.gamm.add.incr, pars=c("betas","gamas"))
stan_dens(fit.gamm.add.incr, pars=c("invtau2","tau", "invsig2","sig", "lambda1","rho1","lambda2","rho2"))
pairs(fit.gamm.add.incr, pars = c("betas","gamas"), las = 1)
mu1 = get_posterior_mean(fit.gamm.add.incr,"mu1")
plot(x1,mu1[,"mean-all chains"])
mu2 = get_posterior_mean(fit.gamm.add.incr,"mu2")
plot(x1,mu2[,"mean-all chains"],col=as.numeric(id))
mu3 = get_posterior_mean(fit.gamm.add.incr,"mu3")
plot(x1,mu3[,"mean-all chains"])
mu4 = get_posterior_mean(fit.gamm.add.incr,"mu4")
plot(x1,mu4[,"mean-all chains"])




# Comparar


loo3_sample = fit.lme.non.reslope
loo6_sample = fit.lme.incr.reslope
loo9_sample = fit.lme.add.non.reslope
loo12_sample = fit.lme.add.incr.reslope

### we have to extract those log-likelihood terms that we so carefully had Stan calculate for us:
log_lik_3 =extract_log_lik(loo3_sample, merge_chains = F)
log_lik_6 =extract_log_lik(loo6_sample, merge_chains = F)
log_lik_9 =extract_log_lik(loo9_sample, merge_chains = F)
log_lik_12 =extract_log_lik(loo12_sample, merge_chains = F)

r_eff_3 =relative_eff(log_lik_3)
r_eff_6 =relative_eff(log_lik_6)
r_eff_9 =relative_eff(log_lik_9)
r_eff_12 =relative_eff(log_lik_12)


###  look at the results for each model, first the one with mu estimated:
(loo_3 <- loo(log_lik_3, r_eff=r_eff_3))
(loo_6 <- loo(log_lik_6, r_eff=r_eff_6))
(loo_9 <- loo(log_lik_9, r_eff=r_eff_9))
(loo_12 <- loo(log_lik_12, r_eff=r_eff_12))

### The second model fits better than the first one, since its looic is smaller.


###  look at the results for each model, first the one with mu estimated:
(waic_3 <- waic(log_lik_3, r_eff=r_eff_3))
(waic_6 <- waic(log_lik_6, r_eff=r_eff_6))
(waic_9 <- waic(log_lik_9, r_eff=r_eff_9))
(waic_12 <- waic(log_lik_12, r_eff=r_eff_12))

### The second model fits better than the first one, since its looic is smaller.

