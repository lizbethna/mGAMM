\name{lme_incr}
\alias{lme_incr}
\title{Liner Model with repeated measures with increasing constrains}
\description{
Bayesian estimate of the linear regression model parameters.
}
\usage{
      lme_incr(Y=Y, X=X, XI, Z, ZI)
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
  \item{XI}{
  Design matrix of the explanatory variables with incrreasing constrains.
}
  \item{ZI}{
  Design matrix of the random effects or repeated measures with increasing constrains.
}
}
\details{
Bayesian estimate of the regression parameters vector and variance parameter of the linear regression model under normal distribution.

The vector of response \eqn{Y} is

\deqn{ Y = X\beta + XI\beta{I} + Z\gamma + ZI\gamma{I} + \epsilon }
\deqn{ \gamma \sim N(0,\sigma^2)}
\deqn{ \gammaI \sim N(0,\sigma{I}^2)}
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{Z} is the random effects design matrix, \eqn{\gamma} is the random effects regression parameters vector, \eqn{\tau^2} is the variance parameter, and \eqn{\sigma^2} is the variance parameter of the random effects. The design matrices \eqn{XI}, design matrix \eqn{ZI} and varianza parameters \eqn{\sigma{I}^2} are related for the increasing constrains.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
load("data/data_gamm.rda")
attach(data)
data.lme.incr <- list(Y=y, n=length(y),
                      X=cbind(1,x2), p=2,
                      XI=cbind(x1), pI=1,
                      Z=cbind(rep(1,length(y))), q=1,
                      ZI=cbind(time), qI=1,
                      N=n_distinct(id), id=id )
para.lme.incr <- c("betas","betasI","tau","sig","sigI")
# Stan
fit.lme.incr <- stan("stan/lme_incr.stan",
                   data=data.lme.incr)

# summary, plots
print(fit.lme.incr, pars=para.lme.incr)
stan_trace(fit.lme.incr, pars=para.lme.incr)
stan_dens(fit.lme.incr, pars=para.lme.incr)

# criteria
loo_lme_incr = fit.lme.incr
log_lik_lme_incr = extract_log_lik(loo_lme_incr, merge_chains = F)
r_eff_lme_incr = relative_eff(log_lik_lme_incr)
loo(log_lik_lme_incr, r_eff=r_eff_lme_incr)
waic(log_lik_lme_incr, r_eff=r_eff_lme_incr)

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
