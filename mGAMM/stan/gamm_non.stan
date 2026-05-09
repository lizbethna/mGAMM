/*
Bayesian
Generalized Additive Mixed Models
Non constraint
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> q;
  matrix[n,q] Z; // random effects
  int<lower=1> k1;
  int<lower=1> padd;
  matrix[n,k1] Xspl[padd]; // additive effects
  vector[1+k1] zero;
  matrix[k1,k1] S1; //
  int<lower=1> N; // subjects
  int<lower=1> id[n]; // id of subjects
}

parameters{
  vector[p] betas;
  matrix[padd,k1] gamas;
  vector[q] bre[N];   // random effects
  real<lower=0> invtau2;
  real<lower=0> invsig2;
  vector<lower=0>[padd] lambda;
}

transformed parameters{
  vector[n] mu;
  vector[n] MUspl;
  vector[n] MUre;
  real<lower=0> tau;
  real<lower=0> sig;
  vector[padd] rho;
  matrix[k1,k1] K1[padd];
  MUspl = X*betas;
  for(j in 1:padd){
    MUspl = MUspl + Xspl[j]*to_vector(gamas[j,]);
  }
  for(t in 1:n){
    MUre[t] = Z[t,]*bre[id[t]] ;
  }
  // expected response
  mu = MUspl + MUre;
  tau = pow(invtau2,-0.5);
  sig = pow(invsig2,-0.5);
  for(j in 1:padd){
    K1[j] = S1 * lambda[j] ;
    rho[j] = log(lambda[j]);
  }
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  invsig2 ~ gamma(.5,.05); // precision parameter prior
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate!
  // prior for s(x1) ...
  betas ~ normal(0,9.9e+06);
  for(j in 1:padd){
    gamas[j,] ~ multi_normal(zero[(1+1):(1+k1)],K1[j]);
    lambda[j] ~ gamma(.05,.005);
  }
  // smoothing parameter priors CHECK...
  for(i in 1:N){
    bre[i] ~ normal(0, sig);
  }
  Y ~ normal(mu, tau);   // response
}

generated quantities {
  vector[n] log_lik;
  for(t in 1:n){
    log_lik[t] = normal_lpdf(Y[t] | mu[t], tau);
  }
}


