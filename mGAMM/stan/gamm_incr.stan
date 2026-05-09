/*
Bayesian
Generalized Additive Mixed Models
Increasing constraint
*/

data{
  int<lower=1> n;
  vector[n] Y; // observations
  int<lower=1> p;
  matrix[n,p] X; // fixed effects
  int<lower=1> pI;
  matrix[n,pI] XI; // fixed effects with increasing constrains
  int<lower=1> q;
  matrix[n,q] Z; //random effects
  int<lower=1> qI;
  matrix[n,qI] ZI; // random effects with increasing constrains
  int<lower=1> k1;
  int<lower=1> padd;
  matrix[n,k1] Xspl[padd]; // additive effects
  matrix[k1,k1] S1;
  int<lower=1> kI1;
  int<lower=1> pIadd;
  matrix[n,kI1] XIspl[pIadd]; // additive effects with increasing constrains
  matrix[kI1,kI1] S1I;
  vector[1+k1+kI1] zero;
  int<lower=1> N;
  int<lower=1> id[n];
}

parameters{
  vector[p] betas;
  vector<lower=0>[pI] betasI;
  matrix[padd,k1] gamas;
  matrix<lower=0>[padd,kI1] gamasI;
  vector[q] bre[N];   // random effects
  vector<lower=0>[qI] breI[N];   // random effects
  real<lower=0> invtau2;
  real<lower=0> invsig2;
  vector<lower=0>[padd] lambda;
  vector<lower=0>[pIadd] lambdaI;
}

transformed parameters{
  vector[n] mu;
  vector[n] MUspl;
  vector[n] MUIspl;
  vector[n] MUre;
  real<lower=0> tau;
  real<lower=0> sig;
  vector[padd] rho;
  vector[pIadd] rhoI;
  matrix[k1,k1] K1[padd];
  matrix[kI1,kI1] KI1[pIadd];

  MUspl = X*betas ;
  for(j in 1:padd){
    MUspl = MUspl + Xspl[j]*to_vector(gamas[j,]) ;
  }
  MUIspl = XI*betasI ;
  for(j in 1:pIadd){
    MUIspl = MUIspl + XIspl[j]*to_vector(gamasI[j,]) ;
  }
  for(t in 1:n){
    MUre[t] = Z[t,]*bre[id[t]] + ZI[t,]*breI[id[t]] ;
  }
  mu =  MUspl + MUIspl + MUre; // expected response

  tau = pow(invtau2,-0.5);
  sig = pow(invsig2,-0.5);
  for(j in 1:padd){
    K1[j] = S1 * lambda[j] ;
    rho[j] = log(lambda[j]);
  }
  for(j in 1:pIadd){
    KI1[j] = S1I * lambdaI[j] ;
    rhoI[j] = log(lambdaI[j]);
  }
}

model {
  invtau2 ~ gamma(.05,.005); // precision parameter prior
  invsig2 ~ gamma(.5,.05); // precision parameter prior
  // Parametric effect priors CHECK invtau2=1/56^2 is appropriate!
  // prior for s(x1) ...
  betas ~ normal(0,9.9e+06);
  betasI ~ normal(0,9.9e+06);
  for(j in 1:padd){
    gamas[j,] ~ multi_normal(zero[(1+1):(1+k1)],K1[j]);
    lambda[j] ~ gamma(.05,.005);
  }
  for(j in 1:pIadd){
    gamasI[j,] ~ multi_normal(zero[(1+1):(1+kI1)],KI1[j]);
    lambdaI[j] ~ gamma(.05,.005);
  }
  for(i in 1:N){
    bre[i] ~ normal(0, sig);
    breI[i] ~ normal(0, sig);
  }
  Y ~ normal(mu, tau);   // response
}

generated quantities {
  vector[n] log_lik;
  for(t in 1:n){
    log_lik[t] = normal_lpdf(Y[t] | mu[t], tau);
  }

}


