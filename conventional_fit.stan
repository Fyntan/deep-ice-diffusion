//Conventional diffusion model (beta = 0)

data {
  int<lower = 0> N;
  vector[N] freq;
  vector[N] spec;
  real<lower = 0> dz;
  
  real alpha_mu;
  real<lower = 0> alpha_sigma;
  real noise_mu;
  real<lower = 0> noise_sigma;
  real sigma_mu;
  real<lower = 0> sigma_sigma;
  real phi_scale;
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> noise;
  real<lower = 0> sigma;
  real<lower = 0> phi;
}

transformed parameters {
  vector[N] y_hat;
  real epsilon = noise^2*dz;
  vector[N] k = 2*pi()*freq;

  for(n in 1:N){
      vector[N] tf;
    tf[n] = exp(-k[n]^2*sigma^2);
    y_hat[n] = alpha*tf[n] + epsilon;
  }
}

model {
  
  alpha ~ normal(alpha_mu, alpha_sigma);
  noise ~ normal(noise_mu, noise_sigma);
  sigma ~ normal(sigma_mu, sigma_sigma);
  
  phi ~ normal(0, phi_scale);
  
  spec ~ gamma(phi, phi ./ y_hat);
  
}
