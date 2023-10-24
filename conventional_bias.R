#Instead of simulating a diffused red noisy time series, do it all in the spectral domain

n <- 1000
dz <- 0.05
alpha_true <- 1
beta_true <- 1.5
sigma_true <- 0.1
noise_true <- 0.07

freq <- (0:(n - 1))/(n*dz)
freq[which(freq > 0.5/dz)] <- freq[which(freq > 0.5/dz)] - 1/dz

k <- 2*pi*freq
tf2 <- exp(-abs(k)^2 * sigma_true^2)

Po <- alpha_true*abs(freq)^-beta_true

P_sim <- c(0, Po[2:n])*tf2 + noise_true^2*dz

#New method

P_mod <- cmdstan_model("~/Desktop/Deep_Diffusion_Length/Stan Modelling/diffusion_length_fit.stan")

alpha_mu <- 1
alpha_sigma <- 0.1
beta_mu <- 1
beta_sigma <- 2
sigma_mu <- 0.4
sigma_sigma <- 0.4
noise_mu <- 0.08
noise_sigma <- 0.1

new_data <- list(N = length(freq[2:501]), freq = freq[2:501], spec = P_sim[2:501], dz = dz,
                  alpha_mu = alpha_mu, alpha_sigma = alpha_sigma, beta_mu = beta_mu, beta_sigma = beta_sigma,
                  noise_mu = noise_mu, noise_sigma = noise_sigma, sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, phi_scale = 1)

P_fit <- P_mod$sample(data = new_data, chains = 4, parallel_chains = 4, refresh = 500)

alpha_est <- mean(P_fit$draws('alpha'))
beta_est <- mean(P_fit$draws('beta'))
sigma_est <- mean(P_fit$draws('sigma'))
noise_est <- mean(P_fit$draws('noise'))

P_est <- alpha_est*freq[2:501]^-beta_est * exp(-(2*pi*freq[2:501]*sigma_est)^2) + noise_est^2*dz

#Convetional method

alpha_mu <- 10
alpha_sigma <- 100
sigma_mu <- 0.4
sigma_sigma <- 0.4
noise_mu <- 0.08
noise_sigma <- 0.1

conv_mod <- cmdstan_model("conventional_fit.stan")

conv_data <- list(N = length(freq[2:501]), freq = freq[2:501], spec = P_sim[2:501], dz = dz,
                  alpha_mu = alpha_mu, alpha_sigma = alpha_sigma,
                  noise_mu = noise_mu, noise_sigma = noise_sigma, sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, phi_scale = 1)

conv_fit <- conv_mod$sample(data = conv_data, chains = 4, parallel_chains = 4, refresh = 500)

alpha_conv_est <- mean(conv_fit$draws('alpha'))
sigma_conv_est <- mean(conv_fit$draws('sigma'))
noise_conv_est <- mean(conv_fit$draws('noise'))

P_conv_est <- alpha_conv_est * exp(-(2*pi*freq[2:501]*sigma_conv_est)^2) + noise_conv_est^2*dz

#par(mfrow = c(1, 2))
#
# plot(freq[2:501], P_sim[2:501], 'l', log = 'xy', lwd = 3, xlab = 'Frequency (arbitrary units)', ylab = 'PSD (arbitrary units)')
# lines(freq[2:501], P_est, col = 'purple', lwd = 3, lty = 2)
# lines(freq[2:501], P_conv_est, col = 'blue', lwd = 3)
# 
# legend(0.4, 8e2, c('theoretical spectrum', 'conventional best fit'), col = c('black', 'blue'), lwd = 3, bty = 'n', seg.len = 1, x.intersp = 0.6, y.intersp = 1)
# 
# plot(freq[2:501]^2, P_sim[2:501], 'l', log = 'y', xlim = c(0, 25), lwd = 3, xlab = 'Frequency^2 (arbitrary units)', ylab = 'PSD (arbitrary units)')
# lines(freq[2:501]^2, P_est, col = 'purple', lwd = 3, lty = 2)
# lines(freq[2:501]^2, P_conv_est, col = 'blue', lwd = 3)
# legend(12, 8e2, c('theoretical spectrum', 'conventional best fit'), col = c('black', 'blue'), lwd = 3, bty = 'n', seg.len = 1, x.intersp = 0.6, y.intersp = 1)

results <- tibble(freq = freq[2:501],
                  P_sim = P_sim[2:501],
                  P_new_est = P_est,
                  P_conv_est = P_conv_est)

write.csv(results, "bias_data.csv")




