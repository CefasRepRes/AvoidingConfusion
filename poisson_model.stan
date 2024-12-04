data {
  int<lower=0> L;
  int<lower=0> M;
  array[L, M] int<lower=0> w;
  array[M, M] int<lower=0> C;
  matrix<lower=0> [M, M] alpha_P;
}

parameters {
  vector <lower=0>[M] mu;
  array[M] simplex[M] P;
}

transformed parameters {
  vector[M] theta = rep_vector(0.0,M);
  for (m in 1:M){
    for (n in 1:M){
      theta[m] += mu[n] * P[n,m];
    }
  }
}

model {
  
  mu ~ gamma(1,0.001);
  
  for (m in 1:M) {
   P[m] ~ dirichlet(alpha_P[m,]'); 
  }
  
  for (m in 1:M) {
    C[m,] ~ multinomial(P[m]);
  }
  
  for (m in 1:M) {
    for (n in 1:L) {
      w[n,m] ~ poisson(theta[m]);
    }
  }
  
}
