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
n <- 100
set.seed(123)

## Parameters
tau <- 1.0   # sd
beta   <- c(1,2,3)   # regression coefficients

## explanatory variables
x1 <- rnorm(n)
x2 <- rbinom(n,size=1,prob=0.5)
X = cbind(1,x1,x2)   # design matrix

## Random error
epsilon = rnorm(n,mean=0,sd=tau)
Y = X%*%beta + epsilon

##################################################
# Fitting
lm_non(Y,X)

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
