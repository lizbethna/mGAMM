\name{lm_non}
\alias{lm_non}
\title{Liner Model}
\description{
Bayesian estimate of the linear regression model parameters.
}
\usage{
      lm_non(Y=Y, X=X)
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

\deqn{ Y = X\beta + \epsilon }
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, and \eqn{\tau^2} is the variance parameter.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
load("data/data_gam.rda")
attach(data)
data.lm.non <- list( Y=y , n=length(y),
                     X=cbind(1,x1,x2) , p=3 )
para.lm.non <- c("betas","tau")

# Stan
fit.lm.non <- stan("stan/lm_non.stan",
                   data=data.lm.non)

# summary, plots
print(fit.lm.non, pars=para.lm.non)
stan_trace(fit.lm.non, pars=para.lm.non)
stan_dens(fit.lm.non, pars=para.lm.non)

# criteria
loo_lm_non = fit.lm.non
log_lik_lm_non = extract_log_lik(loo_lm_non, merge_chains = F)
r_eff_lm_non = relative_eff(log_lik_lm_non)
loo(log_lik_lm_non, r_eff=r_eff_lm_non)
waic(log_lik_lm_non, r_eff=r_eff_lm_non)

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
