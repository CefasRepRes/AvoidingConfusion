data{
  int<lower=0> L;
  int<lower=0> M;
  array[L, M] int<lower=0> w;
  array[M, M] int<lower=0> C;
}
transformed data{
  array[L] int y;
  for (i in 1:L){
    y[i] = sum(w[i,]);
  }
}
parameters{
  real log_mu;
  real <lower=0> k_nb;
  real <lower=0> b_par;
  simplex [M] alpha_b;
  array[M] simplex[M] P;
  array[L] simplex[M] q;
}
transformed parameters{
  real mu = exp(log_mu);
  vector [M] alphas = b_par * alpha_b;
  vector [M] mus = mu * alpha_b;
  vector [M] intravar = mu * (b_par/(b_par + 1.0)) * alpha_b .* (rep_vector(1.0, M) - alpha_b);
  vector [M] intervar = ((mu/(b_par + 1.0)) * (mu + 1.0 + (mu/k_nb)) * alpha_b .* (rep_vector(1.0, M) - alpha_b)) + (mu * (1.0 + (mu/k_nb)) * alpha_b .* alpha_b);
  vector [M] totalvar = intravar + intervar;
  array [M, M] real covar;
  array [M, M] real corr;
  
  for (i in 1:M){
    for (j in 1:M){
      covar[i,j] = (mu^2) * alpha_b[i] * alpha_b[j] * (b_par - k_nb)/(k_nb * (b_par + 1));
      corr[i,j] = covar[i,j]/sqrt(totalvar[i] * totalvar[j]);
    }
  }
  
  array[L] vector[M] PQ;
  
  for (j in 1:L){
    for (k in 1:M){
      PQ[j][k] = 0.0;
      for(i in 1:M){
        PQ[j][k] += P[i,k] * q[j][i];
      }
    }
  }
}
model{
  // priors
  for (m in 1:M) {
    P[m] ~ dirichlet(rep_vector(1.0, M)); 
  }
  log_mu ~ normal(0,1000);
  k_nb ~ normal(0,1000);
  b_par ~ normal(0,1000);
  for (i in 1:L) {
    q[i] ~ dirichlet(alphas); 
    w[i,] ~ multinomial(PQ[i]);
    y[i] ~ neg_binomial_2(mu, k_nb);
  }
  for (m in 1:M) {
    C[m,] ~ multinomial(P[m]);
  }
}
