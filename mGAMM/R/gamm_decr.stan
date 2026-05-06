/*
Bayesian
Generalized Additive Mixed Models
Decreasing constraint
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> pD;
  matrix[n,pD] XD; // fixed effects with decreasing constrains
  int<lower=1> q;
  matrix[n,q] Z; // random effects
  int<lower=1> qD;
  matrix[n,qD] ZD; //random effects with decreasing constrains
  int<lower=1> k1;
  int<lower=1> padd;
  matrix[n,k1] Xspl[padd]; // additive effects
  matrix[k1,k1] S1;
  int<lower=1> kD1;
  int<lower=1> pDadd;
  matrix[n,kD1] XDspl[pDadd]; // additive effects with decreasing constrains
  matrix[kD1,kD1] S1D;
  vector[1+k1+kD1] zero;
  int<lower=1> N;
  int<lower=1> id[n];
}

parameters{
  vector[p] betas;
  vector<upper=0>[pD] betasD;
  matrix[padd,k1] gamas;
  matrix<upper=0>[padd,kD1] gamasD;
  vector[q] bre[N];   // random effects
  vector<upper=0>[qD] breD[N];   // random effects
  real<lower=0> invtau2;
  real<lower=0> invsig2;
  vector<lower=0>[padd] lambda;
  vector<lower=0>[pDadd] lambdaD;
}

transformed parameters{
  vector[n] mu;
  vector[n] MUspl;
  vector[n] MUDspl;
  vector[n] MUre;
  real<lower=0> tau;
  real<lower=0> sig;
  vector[padd] rho;
  vector[pDadd] rhoD;
  matrix[k1,k1] K1[padd];
  matrix[kD1,kD1] KD1[pDadd];

  MUspl = X*betas ;
  for(j in 1:padd){
    MUspl = MUspl + Xspl[j]*to_vector(gamas[j,]) ;
  }
  MUDspl = XD*betasD ;
  for(j in 1:pDadd){
    MUDspl = MUDspl + XDspl[j]*to_vector(gamasD[j,]) ;
  }
  for(t in 1:n){
    MUre[t] = Z[t,]*bre[id[t]] + ZD[t,]*breD[id[t]] ;
  }
  mu =  MUspl + MUDspl + MUre; // expected response

  tau = pow(invtau2,-0.5);
  sig = pow(invsig2,-0.5);
  for(j in 1:padd){
    K1[j] = S1 * lambda[j] ;
    rho[j] = log(lambda[j]);
  }
  for(j in 1:pDadd){
    KD1[j] = S1D * lambdaD[j] ;
    rhoD[j] = log(lambdaD[j]);
  }
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  invsig2 ~ gamma(.5,.05); // precision parameter prior
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate!
  // prior for s(x1) ...
  betas ~ normal(0,9.9e+06);
  betasD ~ normal(0,9.9e+06);
  for(j in 1:padd){
    gamas[j,] ~ multi_normal(zero[(1+1):(1+k1)],K1[j]);
    lambda[j] ~ gamma(.05,.005);
  }
  for(j in 1:pDadd){
    gamasD[j,] ~ multi_normal(zero[(1+1):(1+kD1)],KD1[j]);
    lambdaD[j] ~ gamma(.05,.005);
  }
  for(i in 1:N){
    bre[i] ~ normal(0, sig);
    breD[i] ~ normal(0, sig);
  }

  Y ~ normal(mu, tau);   // response
}

generated quantities {
  vector[n] log_lik;
  for(t in 1:n){
    log_lik[t] = normal_lpdf(Y[t] | mu[t], tau);
  }

}


