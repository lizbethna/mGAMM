\name{gamm_decr}
\alias{gamm_decr}
\title{Generalized additive model wit repeated measures with decreasing constrains}
\description{
Bayesian estimate of the generalized additive model.
}
\usage{
      gamm_decr(Y=Y, X=X,XD=XD, Z=Z,ZD=ZD, Xspl=Xspl,XDspl=XDspl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  Vector of response variables.
}
  \item{Z}{
  Design matrix of the explanatory variables with repeated measures.
}
\item{X}{
  Design matrix of the explanatory variables.
}
  \item{Xspl}{
  Design matrix for the B-spline basis function.
}
}
\details{
Bayesian estimate of the parameters generalized additive model.

The vector of response \eqn{Y} is

\deqn{ Y = X\beta + Z\gamma + g(Xspl) + XD\beta{D} + ZD\gamma{D} + g(XDspl)+  \epsilon }
\deqn{ \gamma \sim N(0,\sigma^2)}
\deqn{ \gamma{D} \sim N(0,\sigma{D}^2)}
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{\gamma} is the regression parameters vector for repeated measures, \eqn{g(Xspl)} spline function, \eqn{\sigma^2} and \eqn{\tau^2} are the variance parameters.
Design matrices \eqn{XD}, \eqn{ZD} and \eqn{XDspl}, and parameters \eqn{\beta{D}}, \eqn{\gamma{D}}, \eqn{\sigma{D}^2} for decreasing constrains.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
load("data/data_gamm.rda")
attach(data)

## Design matrice for B-splines basis function
k1 = 6
cuantiles1 = seq(0,1,length.out=k1)[-c(1,k1)]
knots1 = quantile(x1, cuantiles1)
knots2 = quantile(x2, cuantiles1)
XD1 <- ibs(x1, knots = knots1, degree = 1, intercept = TRUE)
XD2 <- ibs(x2, knots = knots2, degree = 1, intercept = TRUE)

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


## generalized additive models
data.gamm.decr <- list(Y=y, n=length(y) ,
                       X=cbind(rep(1,length(y))), p=1,
                       XD=cbind(rep(0,length(y))), pI=1,
                       Z=cbind(1,time), q=2,
                       ZD=cbind(rep(0,length(y))), qI=1,
                       Xspl=list(XD2), k1=k1, S1=S1, padd=1,
                       XDspl=list(XD1), kD1=k1, S1D=S1, pDadd=1,
                       zero = rep(0,1+k1+k1),
                       id=id, N=n_distinct(id))

# Stan
para.gamm.decr <- c("betas","betasD","gamas","gamasD", "tau","sig","rho","rhoD","lambda","lambdaD")

fit.gamm.decr <- stan("stan/gamm_decr.stan",
                   data=data.gamm.decr,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

## summary, plots and criteria

print(fit.gamm.decr, pars=para.gamm.decr)
stan_trace(fit.gamm.decr, pars=para.gamm.decr)
stan_dens(fit.gamm.decr, pars=para.gamm.decr)

mu = get_posterior_mean(fit.gamm.decr,"mu")
estima.gamm.decr = get_posterior_mean(fit.gamm.decr, para.gamm.decr)[,"mean-all chains"]
plot(y,mu[,"mean-all chains"])

# criteria
loo_gamm_decr = fit.gamm.decr
log_lik_gamm_decr = extract_log_lik(loo_gamm_decr, merge_chains = F)
r_eff_gamm_decr = relative_eff(log_lik_gamm_decr)
loo(log_lik_gamm_decr, r_eff=r_eff_gamm_decr)
waic(log_lik_gamm_decr, r_eff=r_eff_gamm_decr)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
