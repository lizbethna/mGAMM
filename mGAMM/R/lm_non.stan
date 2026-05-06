/*
Bayesian
Linear regression model
Non constraints
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
}

parameters{
  vector[p] betas;
  real<lower=0> invtau2;
}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau;
  mu = X*betas ; // expected response
  tau = pow(invtau2,-0.5);
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  betas ~ normal(0,9.9e+06);
  Y ~ normal(mu, tau);   // response
}

generated quantities {
  vector[n] log_lik;
  for(i in 1:n){
    log_lik[i] = normal_lpdf(Y[i] | mu[i], tau);
  }
}

