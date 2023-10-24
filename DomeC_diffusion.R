#Bayesian MIS 19 diffusion length estimate

library(tidyverse)
library(cmdstanr)
#library(dplyr)

DomeC <- read_tsv("/Users/fshaw/Desktop/Diffusion_Length/Vasileios_etal-2021.tab")

DomeC <- DomeC[-c(1:37),] %>%
  separate(names(DomeC), into = c("depth_top", "depth_bot", "depth", "age", "d18O"), sep = "\t") %>%
  transmute(depth_top = as.numeric(depth_top),
            depth_bot = as.numeric(depth_bot),
            depth = as.numeric(depth),
            age = as.numeric(age),
            d18O = as.numeric(d18O)
  )

P_conv_mod <- cmdstan_model("conventional_fit.stan")

#General function for finding spectra from potentially unevenly spaced data
find_spec <- function(x, d18O, smooth = NULL){
  
  reg_x <- seq(x[1], x[length(x)], length.out = length(x))
  
  reg_d18O <- approx(x, d18O, reg_x)$y
  
  dx <- reg_x[2] - reg_x[1]
  ts <- ts(reg_d18O, start = reg_x[1], deltat = dx)
  sp <- spectrum(ts, plot = F, spans = smooth, taper = 0.1, pad = 0, fast = T, detrend = T)
  
  freq <- sp$freq
  spec <- sp$spec
  
  return(list(freq = freq, spec = spec))
}

alpha_mu <- 0.1
alpha_sigma <- 0.5
sigma_mu <- 0.4
sigma_sigma <- 0.4
noise_mu <- 0.07
noise_sigma <- 0.02

dz <- mean(diff(DomeC$depth))
N <- length(DomeC$depth)

d18O <- approx(DomeC$d18O, xout = 1:length(DomeC$d18O))$y

DomeC_ts <- ts(data = d18O, start = DomeC$depth[1], deltat = dz)

L <- 500 #No. of data points per section
#n <- N - L #For highest resolution
n <- floor(N/L) #For lowest resolution that uses all data points given L
section <- list()
depth_est <- vector(length = n)
spec <- list()
DomeC_results <- list()
sigma_mean <- vector(length = n)
sigma_sd <- vector(length = n)
DomeC_P_fit <- list()

for (i in 1:n){
  
  #section[[i]] <- window(DomeC_ts, start = DomeC$depth[i], end = DomeC$depth[i + L - 1]) #Highest resolution
  section[[i]] <- window(DomeC_ts, start = DomeC$depth[N - i * L + 1], end = DomeC$depth[N - (i - 1) * L]) #Lowest resolution given L, starting from deepest section
  
  depth_est[i] <- time(section[[i]])[L/2]
  
  #section_NAs <- sum(is.na(DomeC$d18O[i:(i + L - 1)])) #Highest resolution
  section_NAs <- sum(is.na(DomeC$d18O[(N - i * L + 1):(N - (i - 1) * L)])) #Lowest resolution given L, starting from deepest section
  
  if (section_NAs/L < 0.25){
    
    spec[[i]] <- find_spec(x = time(section[[i]]), d18O = section[[i]])
    
    data <- list(N = length(spec[[i]]$freq), freq = spec[[i]]$freq, spec = spec[[i]]$spec, dz = dz,
                 alpha_mu = alpha_mu, alpha_sigma = alpha_sigma, sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, noise_mu = noise_mu, noise_sigma = noise_sigma, phi_scale = 1)
    
    DomeC_fit <- P_conv_mod$sample(data = data, chains = 4, parallel_chains = 4, refresh = 500)
    
    DomeC_results[[i]] <- tibble(alpha = as.vector(DomeC_fit$draws('alpha')),
                                 sigma = as.vector(DomeC_fit$draws('sigma')))
    
    sigma_mean[i] <- mean(DomeC_results[[i]]$sigma)
    sigma_sd[i] <- sd(DomeC_results[[i]]$sigma)
    
    DomeC_P_fit[[i]] <- mean(DomeC_results[[i]]$alpha) * exp(-(2*pi*spec[[i]]$freq*mean(DomeC_results[[i]]$sigma))^2) + mean(DomeC_results[[i]]$noise)^2*dz
    
  } else {
    
    sigma_mean[i] <- NA
    sigma_sd[i] <- NA
    
  }
  
  if (i %% 5 == 0) print(i)
  
}

# plot(depth_est, sigma_mean, 'l')
# lines(depth_est, sigma_mean + 2*sigma_sd, col = 'red')
# lines(depth_est, sigma_mean - 2*sigma_sd, col = 'red')

conv_data <- tibble(depth = depth_est, sigma_mean = sigma_mean, sigma_min = sigma_mean - 2*sigma_sd, sigma_max = sigma_mean + 2*sigma_sd)
conv_data[conv_data < 0] <- 0

write.csv(conv_data, "conv_sigmas.csv")

