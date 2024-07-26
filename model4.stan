data {
  int nsites1;
  int nsites2;
  int nsites3;
  int C1[nsites1];
  int C2[nsites2];
  int y[nsites3];
  vector[nsites1] selev1;
  vector[nsites2] selev2;
  vector[nsites3] selev3;
}
parameters {
  real alpha;
  real beta;
}
model {
  vector[nsites1] lambda1;
  vector[nsites2] lambda2;
  vector[nsites3] psi;
  
  // Priors
  alpha ~ uniform(-10, 10);
  beta ~ normal(0, 100);

  // Likelihood
  for (i in 1:nsites1){
    lambda1[i] = exp(alpha + beta * selev1[i]);
    C1[i] ~ poisson(lambda1[i]);
  }
  for (j in 1:nsites2){
    lambda2[j] = exp(alpha + beta * selev2[j]);
    C2[j] ~ poisson(lambda2[j]) T[1, ];
  }
  for (k in 1:nsites3){
    psi[k] = inv_cloglog(alpha + beta * selev3[k]);
    y[k] ~ bernoulli(psi[k]);
  }
}

generated quantities {
  real mean_lam = exp(alpha);
}
