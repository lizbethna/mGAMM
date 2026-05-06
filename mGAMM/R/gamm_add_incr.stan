/*
Bayesian
*/

data{
  int<lower=1> n;
  vector[n] Y;
  int<lower=1> p;
  matrix[n,p] X;
  int<lower=1> pI;
  matrix[n,pI] XI;
  int<lower=1> q;
  matrix[n,q] Z; //random effects
  int<lower=1> qI;
  matrix[n,qI] ZI; //random effects
  int<lower=1> k1;
  matrix[n,k1] Xspl;
  matrix[k1,k1] S1;
  int<lower=1> kI1;
  matrix[n,kI1] XIspl;
  matrix[kI1,kI1] S1I;
  int<lower=1> k2;
  matrix[n,k2] Zspl;
  matrix[k2,k2] S2;
  int<lower=1> kI2;
  matrix[n,kI2] ZIspl;
  matrix[kI2,kI2] S2I;
  vector[1+k1+k2+kI1+kI2] zero;
  int<lower=1> N;
  int<lower=0> Ni[N+1];
}

parameters{
  vector[p] betas;
  vector[k1] gamas;
  vector[q] bre[N]; // random effects
  vector[k2] gre[N]; // random effects
  real<lower=0> invtau2;
  real<lower=0> invsig2;
  real<lower=0> lambda1;
  real<lower=0> lambda2;
  vector<lower=0>[pI] betasI;
  vector<lower=0>[kI1] gamasI;
  vector<lower=0>[qI] breI[N]; // random effects
  vector<lower=0>[kI2] greI[N]; // random effects
  real<lower=0> lambda1I;
  real<lower=0> lambda2I;

}

transformed parameters{
  vector[n] mu;
  real<lower=0> tau;
  real<lower=0> sig;
  real rho1;
  real rho2;
  matrix[k1,k1] K1;
  matrix[k2,k2] K2;
  real rho1I;
  real rho2I;
  matrix[kI1,kI1] KI1;
  matrix[kI2,kI2] KI2;
  for(t in 1:n){
    mu[t] = X[t,]*betas + XI[t,]*betasI + Z[t,]*bre[id[t]] + ZI[t,]*breI[id[t]] + Xspl[t,]*gamas + XIspl[t,]*gamasI + dot_product(Zspl[t,],gre[id[t]]) + dot_product(ZIspl[t,],greI[id[t]]) ; // expected response
  }
  tau = pow(invtau2,-0.5);
  sig = pow(invsig2,-0.5);
  K1 = S1 * lambda1 ;
  K2 = S2 * lambda2 ;
  rho1 = log(lambda1);
  rho2 = log(lambda2);
  KI1 = S1I * lambda1I ;
  KI2 = S2I * lambda2I ;
  rho1I = log(lambda1I);
  rho2I = log(lambda2I);
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  invsig2 ~ gamma(.5,.05); // precision parameter prior
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate!
  // prior for s(x1) ...
  betas ~ normal(0,9.9e+06);
  betasI ~ normal(0,9.9e+06);
  gamas ~ multi_normal(zero[(1+1):(1+k1)],K1);
  gamasI ~ multi_normal(zero[(1+1):(1+kI1)],KI1);
  // smoothing parameter priors CHECK...
  lambda1 ~ gamma(.05,.005);
  lambda1I ~ gamma(.05,.005);
  lambda2 ~ gamma(.05,.005);
  lambda2I ~ gamma(.05,.005);
  for(i in 1:N){
    bre[i] ~ normal(0, sig);
    breI[i] ~ normal(0, sig);
    gre[i] ~ multi_normal(zero[(1+1):(1+k2)],K2);
    greI[i] ~ multi_normal(zero[(1+1):(1+kI2)],KI2);
  }
  Y ~ normal(mu, tau);   // response

}

generated quantities {
  vector[n] log_lik;
  vector[n] mu1;
  vector[n] mu2;
  vector[n] mu3;
  vector[n] mu4;
  for(t in 1:n){
    log_lik[t] = normal_lpdf(Y[t] | mu[t], tau);
    mu1[t] = X[t,]*betas + Z[t,]*bre[id[t]] ; // expected response
    mu2[t] = X[t,]*betas + XI[t,]*betasI + Z[t,]*bre[id[t]] + ZI[t,]*breI[id[t]] ;
    mu3[t] = X[t,]*betas + XI[t,]*betasI + Z[t,]*bre[id[t]] + ZI[t,]*breI[id[t]] + Xspl[t,]*gamas + dot_product(Zspl[t,],gre[id[t]]) ; // expected response
    mu4[t] = X[t,]*betas + XI[t,]*betasI + Z[t,]*bre[id[t]] + ZI[t,]*breI[i] + Xspl[t,]*gamas + XIspl[t,]*gamasI + dot_product(Zspl[t,],gre[id[t]]) + dot_product(ZIspl[t,],greI[id[t]]) ;
  }

}

