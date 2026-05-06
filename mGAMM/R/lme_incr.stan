/*
Bayesian
Linear mixed effects
Increasing constraints
*/

data{
  int<lower=1> n;
  vector[n] Y;// observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> q;
  matrix[n,q] Z; // random effects
  int<lower=1> pI;
  matrix[n,pI] XI; // fixed effects with increasing constraints
  int<lower=1> qI;
  matrix[n,qI] ZI; // random effects with increasing constraints
  int<lower=1> N;
  int<lower=1> id[n]; // id of subjects
}

parameters{
  vector[p] betas;
  vector<lower=0>[pI] betasI;
  vector[q] bre[N];   // random effects
  vector<lower=0>[qI] breI[N];   // random effects
  real<lower=0> invtau2;
  real<lower=0> invsig2;
  real<lower=0> invsigI2;
}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau;
  real<lower=0> sig;
  real<lower=0> sigI;
  for(t in 1:n){
    mu[t] = X[t,]*betas + XI[t,]*betasI + Z[t,]*bre[id[t]] + ZI[t,]*breI[id[t]] ; // expected response
  }
  tau = pow(invtau2,-0.5);
  sig = pow(invsig2,-0.5);
  sigI = pow(invsigI2,-0.5);
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  invsig2 ~ gamma(.05,.005); // precision parameter prior
  invsigI2 ~ gamma(.05,.005); // precision parameter prior
  betas ~ normal(0,9.9e+06);
  betasI ~ normal(0,9.9e+06);
  for(i in 1:N){
    bre[i] ~ normal(0,sig);
    breI[i] ~ normal(0,sigI);
  }
  Y ~ normal(mu, tau);   // response
}


generated quantities {
  vector[n] log_lik;
  for(t in 1:n){
    log_lik[t] = normal_lpdf(Y[t] | mu[t], tau);
  }
}



