/*
Bayesian
Generalized additive models
Increasing constraints
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> pI;
  matrix[n,pI] XI; // fixed effects with increasing constrains
  int<lower=1> k1;
  int<lower=1> padd;
  matrix[n,k1] Xspl[padd]; // additive effects
  int<lower=1> k1I;
  int<lower=1> pIadd;
  matrix[n,k1I] XIspl[pIadd]; // additive effects with increasing constrains
  vector[1+k1+k1I] zero;
  matrix[k1,k1] S1;
  matrix[k1I,k1I] S1I;
}

parameters{
  vector[p] betas;
  vector<lower=0>[pI] betasI;
  matrix[padd,k1] gamas;
  matrix<lower=0>[pIadd,k1I] gamasI;
  real<lower=0> invtau2;
  vector<lower=0>[padd] lambda;
  vector<lower=0>[pIadd] lambdaI;
}

transformed parameters{
  vector[n] mu;
  vector[n] MUspl;
  vector[n] MUIspl;
  real<lower=0> tau;
  vector[padd] rho;
  vector[pIadd] rhoI;
  matrix[k1,k1] K1[padd];
  matrix[k1I,k1I] K1I[pIadd];
  MUspl = X*betas ;
  for(j in 1:padd){
    MUspl = MUspl + Xspl[j]*to_vector(gamas[j,]) ;
  }
  MUIspl = XI*betasI ;
  for(j in 1:pIadd){
    MUIspl = MUIspl + XIspl[j]*to_vector(gamasI[j,]) ;
  }
  mu =  MUspl + MUIspl ; // expected response
  tau = pow(invtau2,-0.5);
  for(j in 1:padd){
    K1[j] = S1 * lambda[j] ;
    rho[j] = log(lambda[j]);
  }
  for(j in 1:pIadd){
    K1I[j] = S1I * lambdaI[j] ;
    rhoI[j] = log(lambdaI[j]);
  }
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate!
  // prior for s(x1) and s(x2)...
  betas ~ normal(0,9.9e+06);
  betasI ~ normal(0,9.9e+06);
  for(j in 1:padd){
    gamas[j,] ~ multi_normal(zero[(1+1):(1+k1)],K1[j]);
    lambda[j] ~ gamma(.05,.005);
  }
  for(j in 1:pIadd){
    gamasI[j,] ~ multi_normal(zero[(1+1):(1+k1I)],K1I[j]);
    lambdaI[j] ~ gamma(.05,.005);
  }
  Y ~ normal(mu, tau);   // response

}

generated quantities {
  vector[n] log_lik;
  for(i in 1:n){
    log_lik[i] = normal_lpdf(Y[i] | mu[i], tau);
  }
}


