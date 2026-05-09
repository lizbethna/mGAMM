\name{gamm_incr}
\alias{gamm_incr}
\title{Generalized additive model wit repeated measures with increasing constrains}
\description{
Bayesian estimate of the generalized additive model.
}
\usage{
      gamm_incr(Y=Y, X=X,XI=XI, Z=Z,ZI=ZI, Xspl=Xspl,XIspl=XIspl)
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

\deqn{ Y = X\beta + Z\gamma + g(Xspl) + XI\beta{I} + ZI\gamma{I} + g(XIspl)+  \epsilon }
\deqn{ \gamma \sim N(0,\sigma^2)}
\deqn{ \gamma{I} \sim N(0,\sigma{I}^2)}
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{\gamma} is the regression parameters vector for repeated measures, \eqn{g(Xspl)} spline function, \eqn{\sigma^2} and \eqn{\tau^2} are the variance parameters.
Design matrices \eqn{XI}, \eqn{ZI} and \eqn{XIspl}, and parameters \eqn{\beta{I}}, \eqn{\gamma{I}}, \eqn{\sigma{I}^2} for increasing constrains.
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


## generalized additive models
data.gamm.incr <- list(Y=y, n=length(y) ,
                       X=cbind(rep(1,length(y))), p=1,
                       XI=cbind(rep(0,length(y))), pI=1,
                       Z=cbind(1,time), q=2,
                       ZI=cbind(rep(0,length(y))), qI=1,
                       Xspl=list(XI2), k1=k1, S1=S1, padd=1,
                       XIspl=list(XI1), kI1=k1, S1I=S1, pIadd=1,
                       zero = rep(0,1+k1+k1),
                       id=id, N=n_distinct(id))

# Stan
para.gamm.incr <- c("betas","betasI","gamas","gamasI", "tau","sig","rho","rhoI","lambda","lambdaI")

fit.gamm.incr <- stan("stan/gamm_incr.stan",
                   data=data.gamm.incr,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

## summary, plots and criteria

print(fit.gamm.incr, pars=para.gamm.incr)
stan_trace(fit.gamm.incr, pars=para.gamm.incr)
stan_dens(fit.gamm.incr, pars=para.gamm.incr)

mu = get_posterior_mean(fit.gamm.incr,"mu")
estima.gamm.incr = get_posterior_mean(fit.gamm.incr, para.gamm.incr)[,"mean-all chains"]
plot(y,mu[,"mean-all chains"])

# criteria
loo_gamm_incr = fit.gamm.incr
log_lik_gamm_incr = extract_log_lik(loo_gamm_incr, merge_chains = F)
r_eff_gamm_incr = relative_eff(log_lik_gamm_incr)
loo(log_lik_gamm_incr, r_eff=r_eff_gamm_incr)
waic(log_lik_gamm_incr, r_eff=r_eff_gamm_incr)

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
