\name{lme_decr}
\alias{lme_decr}
\title{Liner Model with repeated measures with decreasing constrains}
\description{
Bayesian estimate of the linear regression model parameters.
}
\usage{
      lme_decr(Y=Y, X=X, XD, Z, ZD)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  Vector of response variables.
}
  \item{X}{
  Design matrix of the explanatory variables.
}
  \item{Z}{
  Design matrix of the random effects or repeated measures.
}
  \item{XD}{
  Design matrix of the explanatory variables with incrreasing constrains.
}
  \item{ZD}{
  Design matrix of the random effects or repeated measures with decreasing constrains.
}
}
\details{
Bayesian estimate of the regression parameters vector and variance parameter of the linear regression model under normal distribution.

The vector of response \eqn{Y} is

\deqn{ Y = X\beta + XD\beta{D} + Z\gamma + ZD\gamma{D} + \epsilon }
\deqn{ \gamma \sim N(0,\sigma^2)}
\deqn{ \gammaI \sim N(0,\sigma{D}^2)}
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{Z} is the random effects design matrix, \eqn{\gamma} is the random effects regression parameters vector, \eqn{\tau^2} is the variance parameter, and \eqn{\sigma^2} is the variance parameter of the random effects. The design matrices \eqn{XD}, design matrix \eqn{ZD} and varianza parameters \eqn{\sigma{D}^2} are related for the decreasing constrains.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
load("data/data_gamm.rda")
attach(data)
data.lme_decr <- list(Y=y, n=length(y),
                      X=cbind(1,x2), p=2,
                      XD=cbind(x1), pI=1,
                      Z=cbind(rep(1,length(y))), q=1,
                      ZD=cbind(time), qI=1,
                      N=n_distinct(id), id=id )
para.lme.decr <- c("betas","betasI","tau","sig","sigI")
# Stan
fit.lme.decr <- stan("stan/lme_decr.stan",
                   data=data.lme.decr)

# summary, plots
print(fit.lme.decr, pars=para.lme.decr)
stan_trace(fit.lme.decr, pars=para.lme.decr)
stan_dens(fit.lme.decr, pars=para.lme.decr)

# criteria
loo_lme_decr = fit.lme.decr
log_lik_lme_decr = extract_log_lik(loo_lme_decr, merge_chains = F)
r_eff_lme_decr = relative_eff(log_lik_lme_decr)
loo(log_lik_lme_decr, r_eff=r_eff_lme_decr)
waic(log_lik_lme_decr, r_eff=r_eff_lme_decr)

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
