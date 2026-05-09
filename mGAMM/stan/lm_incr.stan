/*
Bayesian
Linear regression model
Increasing constraints
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> pI;
  matrix[n,pI] XI; // fixed effects with increasing constraints
}

parameters{
  vector[p] betas;
  vector<lower=0>[pI] betasI;
  real<lower=0> invtau2;
}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau;
  mu = X*betas + XI*betasI ; // expected response
  tau = pow(invtau2,-0.5);
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  betas ~ normal(0,9.9e+06);
  betasI ~ normal(0,9.9e+06);
  Y ~ normal(mu, tau);   // response
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n){
    log_lik[i] = normal_lpdf(Y[i] | mu[i], tau);
  }
}
