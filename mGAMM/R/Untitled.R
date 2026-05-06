##################################################
### jagam_9_simulations_incr_crossval_N20_re1_casoA
##################################################
# library(VGAM)
# library(mgcv)
library(rstan)
# #library(rjags)
# load.module("glm")
# library(plot3D)
# library(dplyr)
library(ggplot2)
# library(corrplot)
library(splines)
library(splines2)
# library(cgam)
# library(MASS)
# library(msm)
# library(nlme)
#
# library(dlookr)
# library(gtools)
# library(gdata)
# library(ggeffects)
# library(sjmisc)
# library(mnormt)
library(loo)
### http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/

#dir <- "~/Documents/Sabatico/Examples/Simula_f1increasing 231028/Resultados N20 re1 casoA/"

### Specifications to simulate data
set.seed(12345)
TT <- 8   # number of time for repetitions
N <- 20   # number of subjects

k1 = 6
k2 = 6
cuantiles1 = c(0.2,0.4, 0.6,0.8)
cuantiles2 =  c(0.2,0.4, 0.6,0.8)

### Subjects and time
id <- rep(1:N, each=TT)   # identify subjects
time <- rep(1:TT, times=N)  # time

f1 <- function(x) 10 + 0.1*(((x-5))^3)    ### Si es monotona
f1error <- function(x) 10 + 0.1*(((x-5))^3) + rnorm(length(x),0,1)   ### Si es monotona, mas un error
#f1 <- function(x) +.5*x - .1*((x-1)^2) + 0.1*((x-4)^3) - 0.08*(I(x-6>0)*(x-6)^4)   ### No es monotona, se busca que sea monotona
f2 <- function(x) x*((x-1)^2)*((x^2)*(10^2) - (3^4)*((x-1)^2))   ### cualquier otra forma


### Fixed effects
beta0 <- 6

### Random effects
sigma2 <- c(0.7,0.15)

### Variance
tau2 <- 1.7   # variance of the latent classes

### Funciones
### Right Normal Truncada I[Z < trb]
rnormright <- function(trb,mu,sig){
  rp <- pnorm(trb, mean=mu, sd=sig)
  u <- rp*runif(1)
  q <- qnorm(u, mean=mu, sd=sig)
  if(!is.finite(q)){ q = trb }
  return(q)
}
### Left Normal Truncada I[tra < Z]
rnormleft <- function(tra,mu,sig){
  rp <- pnorm(tra, mean=mu, sd=sig)
  u <- rp + (1-rp)*runif(1)
  q <- qnorm(u, mean=mu, sd=sig)
  if(!is.finite(q)){ q = tra }
  return(q)
}

ERRORS <- function(obs,est){
  error = obs-est
  error2 = error^2
  mae = mean(abs(error))
  mse = mean(error2)
  rmse = sqrt(mse)
  return(c("mae"=mae, "mse"=mse, "rmse"=rmse))
}

### Subjects and time
time1 <- time ###+ jitter(rep(0,N*TT),factor=10)
time2 <- time^2
time3 <- time^3
time4 <- time^4

### Indicator for the beginning of observations for each subject for long tables format
offset <- c(1,cumsum(table(id))+1)
n = N*TT
Ni = c(0,cumsum(table(id)))+1 ### offset

# 4. Definir la penalización $S1$ y $S2$
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


### http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/




### http://ritsokiguess.site/docs/2019/06/25/going-to-the-loo-using-stan-for-model-comparison/


### SIMULAR DATOS

x1 <- rnorm(N*TT, time, 0.2)
x2 <- runif(N*TT, 0, 1)

### Random effects
gama0 <- rnorm(N,0,sqrt(sigma2[1]))  # random effects
gama1 <- abs(rnorm(N,0,sqrt(sigma2[2])))  # random effects

### Linear predictor and response

f1_obs = f1(x1)
f2_obs = f2(x2)
f1_con_error = f1error(x1)
f1diff = f1_con_error-f1(x1)


### Linear predictor and response
eta <- rep(0,N*TT)   # linear predictor
Ytrue <- rep(NA,N*TT)   # response variable
etaerror = Yerror = rep(0,N*TT)   # error
for(i in 1:N){
  for(t in offset[i]){
    eta[t] <- beta0 + f1(x1[t]) + f2(x2[t]) + gama0[id[t]] + gama1[id[t]]*time[t]   ### random slope
    Ytrue[t] <- rnorm(1, eta[t], sqrt(tau2))
    etaerror[t] <- eta[t] + f1diff[t]    ### random intercept
    Yerror[t] <- Ytrue[t] + f1diff[t]
  }
  for(t in (offset[i]+1):(offset[i+1]-1)){
    eta[t] <- beta0 + f1(x1[t]) + f2(x2[t]) + gama0[id[t]] + gama1[id[t]]*time[t]   ### random slope
    Ytrue[t] <- rnorm(1, eta[t], sqrt(tau2))
    etaerror[t] <- eta[t] + f1diff[t]    ### random intercept
    Yerror[t] <- Ytrue[t] + f1diff[t]
  }
}


# 2. Generar la matriz diseño $X$ para los B-splines
# Generate a basis matrix for Natural Cubic Splines
knots1 = quantile(x1, cuantiles1)
knots2 = quantile(x2, cuantiles2)
X1 <- ns(x = x1, knots = knots1, intercept = TRUE)
X2 <- ns(x = x2, knots = knots2, intercept = TRUE)



## 1.2 lm: Lineal fit without constraints:
data.lm.non <- list( Y=Yerror , n=length(Ytrue),
                     X=cbind(1,x1,x2) , p=3 )
para.lm.non <- c("betas","tau")

fit.lm.non <- stan("R/lm_non.stan",
                   data=data.lm.non,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

print(fit.lm.non, pars=para.lm.non)
stan_trace(fit.lm.non, pars=para.lm.non)
stan_dens(fit.lm.non, pars=para.lm.non)

mu = get_posterior_mean(fit.lm.non,"mu")
ERRORS(Ytrue,mu[,"mean-all chains"])

estima.lm.non = get_posterior_mean(fit.lm.non, para.lm.non)[,"mean-all chains"]

plot(eta,mu[,"mean-all chains"])
plot(Ytrue,mu[,"mean-all chains"])

loo_lm_non = fit.lm.non
log_lik_lm_non = extract_log_lik(loo_lm_non, merge_chains = F)
r_eff_lm_non = relative_eff(log_lik_lm_non)
loo(log_lik_lm_non, r_eff=r_eff_lm_non)
waic(log_lik_lm_non, r_eff=r_eff_lm_non)


## 2.2 lm: Lineal con restricciones creciente  :

data.lm.incr <- list(Y=Yerror, n=length(Ytrue),
                     X=cbind(1,x2), p=2,
                     XI=cbind(x1), pI=1)
para.lm.incr <- c("betas","betasI","tau")

fit.lm.incr <- stan("R/lm_incr.stan",
                    data=data.lm.incr,
                    chains=3,warmup=1000,iter=2000,thin=2,cores=4 )

print(fit.lm.incr, pars=para.lm.incr)
stan_trace(fit.lm.incr, pars=para.lm.incr)
stan_dens(fit.lm.incr, pars=para.lm.incr)

mu = get_posterior_mean(fit.lm.incr,"mu")
ERRORS(Ytrue,mu[,"mean-all chains"])

estima.lm.incr = get_posterior_mean(fit.lm.incr, para.lm.incr)[,"mean-all chains"]
plot(eta,mu[,"mean-all chains"])
plot(Ytrue,mu[,"mean-all chains"])

loo_lm_incr = fit.lm.incr
log_lik_lm_incr = extract_log_lik(loo_lm_incr, merge_chains = F)
r_eff_lm_incr = relative_eff(log_lik_lm_incr)
loo(log_lik_lm_incr, r_eff=r_eff_lm_incr)
waic(log_lik_lm_incr, r_eff=r_eff_lm_incr)



# 3. Generar la matriz diseño $XI1$ para los I-splines
### ibs: integrated basis splines
### degree = 3 cubic splines
XI1 <- ibs(x1, knots = knots1, degree = 1, intercept = TRUE)
XI2 <- ibs(x2, knots = knots2, degree = 1, intercept = TRUE)


## 3.2 GAM: For a spline-based fit without constraints:
data.gam.non <- list(Y=Yerror, n=length(Ytrue) ,
                     X=cbind(rep(1,length(Ytrue))), p=1,
                     Xspl=list(XI1,XI2), k1=k1, S1=S1, padd=2,
                     zero = rep(0,1+k1))
para.gam.non <- c("betas","tau","gamas","lambda")

fit.gam.non <- stan("R/gam_non.stan",
                    data=data.gam.non,
                    chains=3,warmup=1000,iter=2000,thin=2,cores=4)

print(fit.gam.non, pars=para.gam.non)
stan_trace(fit.gam.non, pars=para.gam.non)
stan_dens(fit.gam.non, pars=para.gam.non)

mu1 = get_posterior_mean(fit.gam.non,"mu")
ERRORS(Ytrue,mu1[,"mean-all chains"])

estima.gam.non = get_posterior_mean(fit.gam.non, para.gam.non)[,"mean-all chains"]

loo_gam_non = fit.gam.non
log_lik_gam_non = extract_log_lik(loo_gam_non, merge_chains = F)
r_eff_gam_non = relative_eff(log_lik_gam_non)
loo(log_lik_gam_non, r_eff=r_eff_gam_non)
waic(log_lik_gam_non, r_eff=r_eff_gam_non)

plot(eta,mu1[,"mean-all chains"])
plot(Ytrue,mu1[,"mean-all chains"])

## 3.2. GAM: Spline con restricciones creciente

data.gam.incr <- list(Y=Yerror, n=length(Ytrue) ,
                      X=cbind(rep(0,length(Ytrue))), p=1,
                      XI=cbind(rep(0,length(Ytrue))), pI=1,
                      Xspl=list(XI1), k1=k1, S1=S1, padd=1,
                      XIspl=list(XI2), k1I=k2, S1I=S1, pIadd=1,
                      zero = rep(0,1+k1+k2))
para.gam.incr <- c("betas","betasI","tau","gamas","lambda","gamasI","lambdaI")

fit.gam.incr <- stan("R/gam_incr.stan",
                     data=data.gam.incr,
                     chains=3,warmup=1000,iter=2000,thin=2,cores=4)

mu = get_posterior_mean(fit.gam.incr,"mu")
ERRORS(Ytrue,mu[,"mean-all chains"])

estima.gam.incr = get_posterior_mean(fit.gam.incr, para.gam.incr)[,"mean-all chains"]

loo_gam_incr = fit.gam.incr
log_lik_gam_incr = extract_log_lik(loo_gam_incr, merge_chains = F)
r_eff_gam_incr = relative_eff(log_lik_gam_incr)
loo(log_lik_gam_incr, r_eff=r_eff_gam_incr)
waic(log_lik_gam_incr, r_eff=r_eff_gam_incr)



### ESTIMAR

## 5.2 LME: Lineal fit without constraints:
data.lme.non <- list( Y=Yerror , n=length(Ytrue),
                      X=cbind(1,x1,x2), p=3,
                      Z=cbind(1,time), q=2,
                      N=N , id=id )
para.lme.non <- c("betas","tau","sig")

fit.lme.non <- stan("R/lme_non.stan",
                    data=data.lme.non,
                    chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

print(fit.lme.non, pars=para.lme.non)
stan_trace(fit.lme.non, pars=para.lme.non)
stan_dens(fit.lme.non, pars=para.lme.non)

mu = get_posterior_mean(fit.lme.non,"mu")
ERRORS(Ytrue,mu[,"mean-all chains"])

estima.lme.non = get_posterior_mean(fit.lme.non, para.lme.non)[,"mean-all chains"]

plot(time,mu[,"mean-all chains"])

loo_lme_non = fit.lme.non
log_lik_lme_non = extract_log_lik(loo_lme_non, merge_chains = F)
r_eff_lme_non = relative_eff(log_lik_lme_non)
loo(log_lik_lme_non, r_eff=r_eff_lme_non)
waic(log_lik_lme_non, r_eff=r_eff_lme_non)


## 6.2 LME: Lineal con restricciones creciente  :

data.lme.incr <- list(Y=Yerror, n=length(Ytrue),
                      X=cbind(1,x2), p=2,
                      XI=cbind(x1), pI=1,
                      Z=cbind(rep(1,length(Ytrue))), q=1,
                      ZI=cbind(time), qI=1,
                      N=N, id=id )
para.lme.incr <- c("betas","betasI","tau","sig","sigI")

fit.lme.incr <- stan("R/lme_incr.stan",
                     data=data.lme.incr,
                     chains=3,warmup=1000,iter=2000,thin=2,cores=4 )

print(fit.lme.incr, pars=para.lme.non)
stan_trace(fit.lme.incr, pars=para.lme.non)
stan_dens(fit.lme.incr, pars=para.lme.non)

mu = get_posterior_mean(fit.lme.incr,"mu")
ERRORS(Ytrue,mu[,"mean-all chains"])

estima.lme.incr = get_posterior_mean(fit.lme.incr, para.lme.incr)[,"mean-all chains"]
plot(time,mu[,"mean-all chains"])

loo_lme_incr = fit.lme.incr
log_lik_lme_incr = extract_log_lik(loo_lme_incr, merge_chains = F)
r_eff_lme_incr = relative_eff(log_lik_lme_incr)
loo(log_lik_lme_incr, r_eff=r_eff_lme_incr)
waic(log_lik_lme_incr, r_eff=r_eff_lme_incr)


## 7.2 LME: For a spline-based fit without constraints:
data.gamm.non <- list(Y=Yerror, n=length(Ytrue) ,
                      X=cbind(rep(1,length(Ytrue))), p=1,
                      Z=cbind(1,time), q=2,
                      Xspl=list(XI1,XI2), k1=k1, S1=S1, padd=2,
                      zero = rep(0,1+k1),
                      id=id, N=N)
para.gamm.non <- c("betas","gamas","tau","sig","rho","lambda")

fit.gamm.non <- stan("R/gamm_non.stan",
                     data=data.gamm.non,
                     chains=3,warmup=1000,iter=2000,thin=2,cores=4)

print(fit.gamm.non, pars=para.gamm.non)
stan_trace(fit.gamm.non, pars=para.gamm.non)
stan_dens(fit.gamm.non, pars=para.gamm.non)

mu1 = get_posterior_mean(fit.gamm.non,"MUspl")
mu2 = get_posterior_mean(fit.gamm.non,"MUre")
mu3 = get_posterior_mean(fit.gamm.non,"mu")
ERRORS(Ytrue,mu1[,"mean-all chains"])
ERRORS(Ytrue,mu3[,"mean-all chains"])

loo_gamm_non = fit.gamm.non
log_lik_gamm_non = extract_log_lik(loo_gamm_non, merge_chains = F)
r_eff_gamm_non = relative_eff(log_lik_gamm_non)
loo(log_lik_gamm_non, r_eff=r_eff_gamm_non)
waic(log_lik_gamm_non, r_eff=r_eff_gamm_non)


## 8.2. LME: Spline con restricciones creciente

data.gamm.incr <- list(Y=Yerror, n=length(Ytrue) ,
                       X=cbind(rep(1,length(Ytrue))), p=1,
                       XI=cbind(rep(0,length(Ytrue))), pI=1,
                       Z=cbind(1,time), q=2,
                       ZI=cbind(rep(0,length(Ytrue))), qI=1,
                       Xspl=list(XI2), k1=k2, S1=S2, padd=1,
                       XIspl=list(XI1), kI1=k1, S1I=S1, pIadd=1,
                       zero = rep(0,1+k1+k2),
                       id=id, N=N)
para.gamm.incr <- c("betas","betasI","gamas","gamasI","tau","sig","rho","rhoI","lambda","lambdaI")

fit.gamm.incr <- stan("R/gamm_incr.stan",
                      data=data.gamm.incr,
                      chains=3,warmup=1000,iter=2000,thin=2,cores=4)

print(fit.gamm.incr, pars=para.gamm.incr)
stan_trace(fit.gamm.incr, pars=para.gamm.incr)
stan_dens(fit.gamm.incr, pars=para.gamm.incr)

mu = get_posterior_mean(fit.gamm.incr,"mu")
ERRORS(Ytrue,mu[,"mean-all chains"])

estima.gamm.incr = get_posterior_mean(fit.gamm.incr, para.gamm.incr)[,"mean-all chains"]
plot(time,mu[,"mean-all chains"])
plot(Ytrue,mu[,"mean-all chains"])

loo_gamm_incr = fit.gamm.incr
log_lik_gamm_incr = extract_log_lik(loo_gamm_incr, merge_chains = F)
r_eff_gamm_incr = relative_eff(log_lik_gamm_incr)
loo(log_lik_gamm_incr, r_eff=r_eff_gamm_incr)
waic(log_lik_gamm_incr, r_eff=r_eff_gamm_incr)



### dev, portrait, 5.4 x 3.7
### Graphics
data1 <- data.frame("id"=id, "Ytrue"=Ytrue, "eta"=eta,
                    "time"=time, "x1"=x1, "x2"=x2,
                    "f1_con_error"=f1_con_error, "f1_obs"=f1_obs,
                    "Yerror"=Yerror)
str(data1)

p0 <- ggplot(data1, aes(x=time, y=Ytrue, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=time, y=Ytrue, group=1),
              method="loess", se=TRUE, level=0.95,
              lwd=1, alpha=0.5, color="red") +
  ggtitle("Scenario 2") +
  theme(plot.title = element_text(size=17),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text=element_text(size=13)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=13)) +
  xlab(expression(paste("Time ",t))) +
  ylab(expression(paste("Response ", y)))
p0


p1 <- ggplot(data1, aes(x=x1, y=Ytrue, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x1, y=Ytrue, group=1),
              method="loess", se=TRUE, level=0.95,
              lwd=1, alpha=0.5, color="red") +
  ggtitle("Scenario 2") +
  theme(plot.title = element_text(size=17),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text=element_text(size=13)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=13)) +
  xlab(expression(paste("Variable ",x[1]))) +
  ylab(expression(paste("Response ", y)))
p1

p2 <- ggplot(data1, aes(x=x2, y=Ytrue, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x2, y=Ytrue, group=1),
              method="loess", se=TRUE, level=0.95,
              lwd=1, alpha=0.5, color="red") +
  ggtitle("Scenario 2") +
  theme(plot.title = element_text(size=17),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text=element_text(size=13)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=13)) +
  xlab(expression(paste("Variable ",x[2]))) +
  ylab(expression(paste("Response ", y)))
p2

p3 <- ggplot(data1, aes(x=time, y=eta, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=time, y=eta, group=1),
              method="loess", se=TRUE, level=0.95,
              lwd=1, alpha=0.5, color="red") +
  ggtitle("Scenario 2") +
  theme(plot.title = element_text(size=17),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text=element_text(size=13)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=13)) +
  xlab("time t") +
  ylab(expression(paste("Linear predictor ", eta)))
p3

p4 <- ggplot(data1, aes(x=x1, y=f1_con_error, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x1, y=f1_obs, group=1),
              method="loess", se=TRUE, level=0.100,
              lwd=1, alpha=0.5, color="red") +
  ggtitle("Scenario 2") +
  theme(plot.title = element_text(size=17),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text=element_text(size=13)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=13)) +
  xlab(expression(paste("Variable ",x[1]))) +
  ylab(expression(paste(f[omega](x[1]))))
p4


p5 <- ggplot(data1, aes(x=x1, y=Yerror, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x1, y=beta0+f1_obs+f2_obs, group=1),
              method="loess", se=TRUE, level=0.95,
              lwd=1, alpha=0.5, color="red") +
  ggtitle("Scenario 2") +
  theme(plot.title = element_text(size=17),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text=element_text(size=13)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=13)) +
  xlab(expression(paste("Variable ",x[1]))) +
  ylab(expression(paste("Response ", y)))
p5


p6 <- ggplot(data1, aes(x=x2, y=Yerror, group=id)) +
  theme_bw()+
  geom_point(cex=0.8,color=id) +
  geom_line(lwd=0.1,color=id) +
  stat_smooth(data=data1, aes(x=x2, y=beta0+f1_obs+f2_obs, group=1),
              method="loess", se=TRUE, level=0.95,
              lwd=1, alpha=0.5, color="red") +
  ggtitle("Scenario 2") +
  theme(plot.title = element_text(size=17),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.text=element_text(size=13)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=13)) +
  xlab(expression(paste("Variable ",x[2]))) +
  ylab(expression(paste("Response ", y)))
p6

