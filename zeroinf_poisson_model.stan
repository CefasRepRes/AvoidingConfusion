functions {
real lse(vector x){
   real mx = max(x);
   vector[num_elements(x)] xexp = x - mx;
   return(mx + log(sum(exp(xexp))));
}
  
real zip_lpmf(array[] int w, vector q, vector P, vector mu){
  real ret;
  vector[4] ll;
  ll[1] = log(1-q[1]) + log(1-q[2]) + poisson_lpmf(w[1]|mu[1]*P[1] + mu[2]*(1-P[2])) + poisson_lpmf(w[2]|mu[1]*(1-P[1]) + mu[2]*P[2]);
  ll[2] = log(1-q[1]) + log(q[2]) + poisson_lpmf(w[1]|mu[1]*P[1]) + poisson_lpmf(w[2]|mu[1]*(1-P[1]));
  ll[3] = log(q[1]) + log(1-q[2]) + poisson_lpmf(w[1]|mu[2]*(1-P[2])) + poisson_lpmf(w[2]|mu[2]*P[2]);
  if (w[1] == 0 && w[2] == 0){
    ll[4] = log(q[1]) + log(q[2]); 
    ret = lse(ll);
  }
  else {
    ret = lse(ll[1:3]);
  }
  return(ret);
}
  
}


data {
  int<lower=0> L;
  
  array[L, 2] int<lower=0> w;
  array[2, 2] int<lower=0> C;
  matrix<lower=0> [2, 2] alpha_P;
  matrix<lower=0> [2, 2] alpha_q;
  vector[L] X;
}

transformed data {
 array[2] int sumC;
 for (m in 1:2){
   sumC[m] = sum(C[m,]);
 }
}

parameters {
  vector <lower=0,upper=1>[2] P;
  vector <lower=0,upper=1>[2] q;
  vector[2] intercept;
  vector[2] betas;
}

transformed parameters {

  array[L] vector[2] mu;
  for (n in 1:L){
    for (m in 1:2){
      mu[n][m] = exp(intercept[m] + betas[m]*X[n]);
    }
  }
}

model {
  
  intercept ~ normal(0, 100);
  betas ~ normal(0, 10);
  
  for (m in 1:2) {
    P[m] ~ beta(alpha_P[1,m], alpha_P[2,m]); 
  }
  
  for (m in 1:2) {
    q[m] ~ beta(alpha_q[1,m], alpha_q[2,m]); 
  }
  
  for (m in 1:2) {
       C[m,m] ~ binomial(sumC[m], P[m]);
  }
  
  for (n in 1:L) {
    w[n] ~ zip(q,P,mu[n]);
  }
  
}

