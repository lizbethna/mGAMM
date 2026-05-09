\name{lm_incr}
\alias{lm_incr}
\title{Liner Model with increasing constrains}
\description{
Bayesian estimate of the linear regression model parameters.
}
\usage{
      lm_incr(Y=Y, X=X, XI=XI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  Vector of response variables.
}
  \item{X}{
  Design matrix of the explanatory variables.
}
  \item{XI}{
  Design matrix of the explanatory variables for increasing constrains.
}
}
\details{
Bayesian estimate of the regression parameters vector and variance parameter of the linear regression model under normal distribution.

The vector of response \eqn{Y} is

\deqn{ Y = X\beta + X\beta + XI\beta{I} + \epsilon }
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{XI} is the design matrix for increasing constrains, \eqn{\beta{I}>0} is the regression parameters vector with increasing contrains, and \eqn{\tau^2} is the variance parameter.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
load("data/data_gam.rda")
attach(data)
data.lm.incr <- list( Y=y , n=length(y),
                     X=cbind(1,x2) , p=2,
                     XI=cbind(x1) , pI=1)
para.lm.incr <- c("betas","betasI","tau")

# Stan
fit.lm.incr <- stan("stan/lm_incr.stan",
                   data=data.lm.incr)

# summary, plots
print(fit.lm.incr, pars=para.lm.incr)
stan_trace(fit.lm.incr, pars=para.lm.incr)
stan_dens(fit.lm.incr, pars=para.lm.incr)

# criteria
loo_lm_incr = fit.lm.incr
log_lik_lm_incr = extract_log_lik(loo_lm_incr, merge_chains = F)
r_eff_lm_incr = relative_eff(log_lik_lm_incr)
loo(log_lik_lm_incr, r_eff=r_eff_lm_incr)
waic(log_lik_lm_incr, r_eff=r_eff_lm_incr)

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
