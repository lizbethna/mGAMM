/*
Bayesian
Generalized additive models
Decreasing constraints
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> pD;
  matrix[n,pD] XD; // fixed effects with decreasing constrains
  int<lower=1> k1;
  int<lower=1> padd;
  matrix[n,k1] Xspl[padd]; // additive effects
  int<lower=1> k1D;
  int<lower=1> pDadd;
  matrix[n,k1D] XDspl[pDadd]; // additive effects with decreasing constrains
  vector[1+k1+k1D] zero;
  matrix[k1,k1] S1;
  matrix[k1D,k1D] S1D;
}

parameters{
  vector[p] betas;
  vector<upper=0>[pD] betasD;
  matrix[padd,k1] gamas;
  matrix<upper=0>[pDadd,k1D] gamasD;
  real<lower=0> invtau2;
  vector<lower=0>[padd] lambda;
  vector<lower=0>[pDadd] lambdaD;
}

transformed parameters{
  vector[n] mu;
  vector[n] MUspl;
  vector[n] MUDspl;
  real<lower=0> tau;
  vector[padd] rho;
  vector[pDadd] rhoD;
  matrix[k1,k1] K1[padd];
  matrix[k1D,k1D] K1D[pDadd];
  MUspl = X*betas ;
  for(j in 1:padd){
    MUspl = MUspl + Xspl[j]*to_vector(gamas[j,]) ;
  }
  MUDspl = XD*betasD ;
  for(j in 1:pDadd){
    MUDspl = MUDspl + XDspl[j]*to_vector(gamasD[j,]) ;
  }
  mu =  MUspl + MUDspl ; // expected response
  tau = pow(invtau2,-0.5);
  for(j in 1:padd){
    K1[j] = S1 * lambda[j] ;
    rho[j] = log(lambda[j]);
  }
  for(j in 1:pDadd){
    K1D[j] = S1D * lambdaD[j] ;
    rhoD[j] = log(lambdaD[j]);
  }
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate!
  // prior for s(x1) and s(x2)...
  betas ~ normal(0,9.9e+06);
  betasD ~ normal(0,9.9e+06);
  for(j in 1:padd){
    gamas[j,] ~ multi_normal(zero[(1+1):(1+k1)],K1[j]);
    lambda[j] ~ gamma(.05,.005);
  }
  for(j in 1:pDadd){
    gamasD[j,] ~ multi_normal(zero[(1+1):(1+k1)],K1D[j]);
    lambdaD[j] ~ gamma(.05,.005);
  }
  Y ~ normal(mu, tau);   // response

}

generated quantities {
  vector[n] log_lik;
  for(i in 1:n){
    log_lik[i] = normal_lpdf(Y[i] | mu[i], tau);
  }
}


