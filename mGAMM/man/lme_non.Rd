\name{lme_non}
\alias{lme_non}
\title{Liner Model with repeated measures}
\description{
Bayesian estimate of the linear regression model parameters.
}
\usage{
      lme_non(Y=Y, X=X, Z=Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  Vector of response variables.
}
  \item{X}{
  Design matrix of the explanatory variables.
}
}
\details{
Bayesian estimate of the regression parameters vector and variance parameter of the linear regression model under normal distribution.

The vector of response \eqn{Y} is

\deqn{ Y = X\beta + Z\gamma+ \epsilon }
\deqn{ \gamma \sim N(0,\sigma^2)}
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{Z} is the random effects design matrix, \eqn{\gamma} is the random effects regression parameters vector, \eqn{\tau^2} is the variance parameter, and \eqn{\sigma^2} is the variance parameter of the random effects.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
load("data/data_gamm.rda")
attach(data)
data.lme.non <- list( Y=y , n=length(y),
                      X=cbind(1,x1,x2), p=3,
                      Z=cbind(1,time), q=2,
                      N=n_distinct(id) , id=id )
para.lme.non <- c("betas","tau","sig")

# Stan
fit.lme.non <- stan("stan/lme_non.stan",
                   data=data.lme.non)

# summary, plots
print(fit.lme.non, pars=para.lme.non)
stan_trace(fit.lme.non, pars=para.lme.non)
stan_dens(fit.lme.non, pars=para.lme.non)

# criteria
loo_lme_non = fit.lme.non
log_lik_lme_non = extract_log_lik(loo_lme_non, merge_chains = F)
r_eff_lme_non = relative_eff(log_lik_lme_non)
loo(log_lik_lme_non, r_eff=r_eff_lme_non)
waic(log_lik_lme_non, r_eff=r_eff_lme_non)

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
