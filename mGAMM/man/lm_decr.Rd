\name{lm_decr}
\alias{lm_decr}
\title{Liner Model with decreasing constrains}
\description{
Bayesian estimate of the linear regression model parameters.
}
\usage{
      lm_decr(Y=Y, X=X, XD=XD)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  Vector of response variables.
}
  \item{X}{
  Design matrix of the explanatory variables.
}
  \item{XD}{
  Design matrix of the explanatory variables with decreasing constrains.
}
}
\details{
Bayesian estimate of the regression parameters vector and variance parameter of the linear regression model under normal distribution.

The vector of response \eqn{Y} is

\deqn{ Y = X\beta + XD\beta{D} + \epsilon }
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{X} is the design matrix for decreasing constrains, \eqn{\beta{D}<0} is the regression parameters vector with decreasing constrains,  and \eqn{\tau^2} is the variance parameter.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
load("data/data_gam.rda")
attach(data)
data.lm.decr <- list( Y=y , n=length(y),
                     X=cbind(1,x2) , p=2,
                     XD=cbind(x1), pD=1)
para.lm.decr <- c("betas","betasD","tau")

# Stan
fit.lm.decr <- stan("stan/lm_decr.stan",
                   data=data.lm.decr)

# summary, plots
print(fit.lm.decr, pars=para.lm.decr)
stan_trace(fit.lm.decr, pars=para.lm.decr)
stan_dens(fit.lm.decr, pars=para.lm.decr)

# criteria
loo_lm_decr = fit.lm.decr
log_lik_lm_decr = extract_log_lik(loo_lm_decr, merge_chains = F)
r_eff_lm_decr = relative_eff(log_lik_lm_decr)
loo(log_lik_lm_decr, r_eff=r_eff_lm_decr)
waic(log_lik_lm_decr, r_eff=r_eff_lm_decr)

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
