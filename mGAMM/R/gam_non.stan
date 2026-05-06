/*
Bayesian
Generalized additive models
Non constraints
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> k1;
  int<lower=1> padd;
  matrix[n,k1] Xspl[padd]; // additive effects
  vector[1+k1] zero;
  matrix[k1,k1] S1;
}

parameters{
  vector[p] betas;
  matrix[padd,k1] gamas;
  real<lower=0> invtau2;
  vector<lower=0>[padd] lambda;
}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau;
  vector[padd] rho;
  matrix[k1,k1] K1[padd];
  mu = X*betas;
  for(j in 1:padd){
    mu = mu + Xspl[j]*to_vector(gamas[j,]);
  }
  tau = pow(invtau2,-0.5);
  for(j in 1:padd){
    K1[j] = S1 * lambda[j] ;
    rho[j] = log(lambda[j]);
  }
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate!
  // prior for s(x1) and s(x2)
  betas ~ normal(0,9.9e+06);
  for(j in 1:padd){
    gamas[j,] ~ multi_normal(zero[(1+1):(1+k1)],K1[j]);
    lambda[j] ~ gamma(.05,.005);
  }
  Y ~ normal(mu, tau);   // response
}

generated quantities {
  vector[n] log_lik;
  for(i in 1:n){
    log_lik[i] = normal_lpdf(Y[i] | mu[i], tau);
  }
}


