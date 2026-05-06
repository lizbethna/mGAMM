/*
Bayesian
Linear mixed effects
Non constraints
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> q;
  matrix[n,q] Z; // random effects
  int<lower=1> N; // subjects
  int<lower=1> id[n]; // id of subjects
}

parameters{
  vector[p] betas;
  vector[q] bre[N]; // random effects
  real<lower=0> invtau2;
  real<lower=0> invsig2;
}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau;
  real<lower=0> sig;
  for(t in 1:n){
    mu[t] = X[t,]*betas + Z[t,]*bre[id[t]]  ; // expected response
  }
  tau = pow(invtau2,-0.5);
  sig = pow(invsig2,-0.5);
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  invsig2 ~ gamma(.05,.005); // precision parameter prior
  betas ~ normal(0,9.9e+06);
  for(i in 1:N){
    bre[i] ~ normal(0,sig);
  }
  Y ~ normal(mu, tau);   // response
}


generated quantities {
  vector[n] log_lik;
  for(t in 1:n){
    log_lik[t] = normal_lpdf(Y[t] | mu[t], tau);
  }
}
