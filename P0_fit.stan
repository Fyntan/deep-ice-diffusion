//Linear regression stan model

data {
  int<lower=0> N;
  vector[N] freq;
  vector[N] spec;
  
  real alpha_mu;
  real<lower = 0> alpha_sigma;
  real beta_mu;
  real<lower = 0> beta_sigma;
  real phi_scale;
}

parameters {
  real<lower = 0> alpha;
  real beta;
  real<lower = 0> phi;
}

transformed parameters {
  vector[N] y_hat;
  for(n in 1:N){
    y_hat[n] = alpha*freq[n]^-beta;
  }
}

model {
  
  alpha ~ normal(alpha_mu, alpha_sigma);
  beta ~ normal(beta_mu, beta_sigma);
  
  phi ~ normal(0, phi_scale);
  
  spec ~ gamma(phi, phi ./ y_hat);
}
