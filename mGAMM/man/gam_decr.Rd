\name{gam_decr}
\alias{gam_decr}
\title{Generalized additive model with decreasing function}
\description{
Bayesian estimate of the generalized additive model.
}
\usage{
      gam_decr(Y=Y, X=X, XD=XD, Xspl=Xspl, XDspl=XDspl)
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
  \item{Xspl}{
  Design matrix of the B-splines basis function.
}
  \item{XDspl}{
  Design matrix of the B-splines basis function with decreasing constrains.
}
}
\details{
Bayesian estimate of the parameters generalized additive model with decreasing constrains.

The vector of response \eqn{Y} is

\deqn{ Y = X\beta + XD\beta{D} + g(Xspl) + g(XDspl) + \epsilon }
\deqn{ \epsilon \sim N(0,\tau^2)}

where \eqn{X} is the design matrix, \eqn{\beta} is the regression parameters vector, \eqn{g(Xspl)} spline function, and \eqn{\tau^2} is the variance parameter. Design matrix \eqn{XD}, vector parameters \eqn{\beta{D}}, and \eqn{g(XDspl)} B-splines for decreasing constrains.
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
data.gam.decr <- list(Y=y, n=length(y) ,
                      X=cbind(rep(0,length(y))), p=1,
                      XD=cbind(rep(0,length(y))), pI=1,
                      Xspl=list(XD1), k1=k1, S1=S1, padd=1,
                      XDspl=list(XD2), k1I=k1, S1I=S1, pIadd=1,
                      zero = rep(0,1+k1))

# Stan
para.gam.decr <- c("betas","betasI","tau","gamas","lambda","gamasI","lambdaI")

fit.gam.decr <- stan("stan/gam_decr.stan",
                   data=data.gam.decr,
                   chains=3, warmup=1000, iter=2000, thin=2, cores=4 )

## summary, plots and criteria

print(fit.gam.decr, pars=para.gam.decr)
stan_trace(fit.gam.decr, pars=para.gam.decr)
stan_dens(fit.gam.decr, pars=para.gam.decr)

mu = get_posterior_mean(fit.gam.decr,"mu")
estima.gam.decr = get_posterior_mean(fit.gam.decr, para.gam.decr)[,"mean-all chains"]
plot(y,mu[,"mean-all chains"])

# criteria
loo_gam_decr = fit.gam.decr
log_lik_gam_decr = extract_log_lik(loo_gam_decr, merge_chains = F)
r_eff_gam_decr = relative_eff(log_lik_gam_decr)
loo(log_lik_gam_decr, r_eff=r_eff_gam_decr)
waic(log_lik_gam_decr, r_eff=r_eff_gam_decr)

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
