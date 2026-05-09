\name{gam_non}
\alias{gam_non}
\title{Generalized additive model}
\description{
Bayesian estimate of the generalized additive model.
}
\usage{
      gam_non(Y=Y, X=X, Xspl=Xspl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  Vector of response variables.
}
  \item{X}{
  Design matrix of the explanatory variables.
}
  \item{Xspl}{
  Design matrix for the B-spline bassis function.
}
}
\details{
Bayesian estimate of the parameters generalized additive model.

The vector of response \eqn{Y} is

\deqn{ Y = X\beta + g(Xspl)+ \epsilon }
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{g(Xspl)} spline function, and \eqn{\tau^2} is the variance parameter.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
load("data/data_gam.rda")
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
data.gam.non <- list(Y=y, n=length(y) ,
                     X=cbind(rep(1,length(y))), p=1,
                     Xspl=list(XI1,XI2), k1=k1, S1=S1, padd=2,
                     zero = rep(0,1+k1))

# Stan
para.gam.non <- c("betas","tau","gamas","lambda")

fit.gam.non <- stan("stan/gam_non.stan",
                   data=data.gam.non,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

## summary, plots and criteria

print(fit.gam.non, pars=para.gam.non)
stan_trace(fit.gam.non, pars=para.gam.non)
stan_dens(fit.gam.non, pars=para.gam.non)

mu = get_posterior_mean(fit.gam.non,"mu")
estima.gam.non = get_posterior_mean(fit.gam.non, para.gam.non)[,"mean-all chains"]
plot(y,mu[,"mean-all chains"])

# criteria
loo_gam_non = fit.gam.non
log_lik_gam_non = extract_log_lik(loo_gam_non, merge_chains = F)
r_eff_gam_non = relative_eff(log_lik_gam_non)
loo(log_lik_gam_non, r_eff=r_eff_gam_non)
waic(log_lik_gam_non, r_eff=r_eff_gam_non)

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
