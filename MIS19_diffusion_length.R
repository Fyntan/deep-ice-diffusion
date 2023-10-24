#Bayesian MIS 19 diffusion length estimate

library(tidyverse)
library(cmdstanr)
library(dplyr)
library(gridExtra)
library(RColorBrewer)

# Data --------------------------------------------------------------------

#Load Dome C data as a tibble
DomeC <- read_tsv("/Users/fshaw/Desktop/Diffusion_Length/Vasileios_etal-2021.tab")

#Transform the data
DomeC <- DomeC[-c(1:37),] %>%
  separate(names(DomeC), into = c("depth_top", "depth_bot", "depth", "age", "d18O"), sep = "\t") %>%
  transmute(depth_top = as.numeric(depth_top),
            depth_bot = as.numeric(depth_bot),
            depth = as.numeric(depth),
            age = as.numeric(age),
            d18O = as.numeric(d18O)
  )

#Define time-series for analysis
MIS1 <- DomeC[1:3150, ]
MIS5 <- DomeC[13926:15701, ] 
MIS9 <- DomeC[22900:23500, ]
MIS19 <- DomeC[28650:28890, ]

#Gaps
#MIS1: Depths 119.955 - 120.615, Times 3.14848 - 3.17031
#          Depths 234.355 - 234.795, Times 7.26902 - 7.28199
#          Depths 288.805 - 290.345, Times 9.32127 - 9.37726
#MIS5: Depths 1583.505 - 1594.395, Times 120.2833 - 121.1313
#      Depths 1638.505 - 1649.395, Times 124.5008 - 125.3630
#      Depths 1693.505 - 1704.395, Times 128.5310 - 129.2680


# Run First! --------------------------------------------------------------

#Data from another script "DomeC_diffusion.R", finding conventional diffusion length estimates for the whole Dome C ice core
conv_sigmas <- read.csv("conv_sigmas.csv")

#Data from another script "conventional_bias.R", visualising the bias of conventional estimator on simulated data
bias_data <- read_csv('bias_data.csv')


# Stan Models -------------------------------------------------------------

#Model for power law fit for P0
P0_mod <- cmdstan_model("P0_fit.stan")

#Model for diffusion length fit
P_mod <- cmdstan_model("diffusion_length_fit.stan")


# Functions ---------------------------------------------------------------

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

#Function that calculates the frequency limits for our fitting procedures (Not used in manuscript results)
define_freq_range <- function(min_amp, max_amp, sigma_min, sigma_max){
  
  min_freq <- sqrt(-log(min_amp))/(2*pi*sigma_max)
  max_freq <- sqrt(-log(max_amp))/(2*pi*sigma_min)
  
  return(list(min = min_freq, max = max_freq))
}

#Applies Stan power law fit
find_P0_fit <- function(freq, spec, alpha_mu = 0.1, alpha_sigma = 1, beta_mu = 1.5, beta_sigma = 1, phi_scale = 1){
  
  stan_data <- list(N = length(freq), freq = freq, spec = spec, 
                    alpha_mu = alpha_mu, alpha_sigma = alpha_sigma, beta_mu = beta_mu, beta_sigma = beta_sigma, phi_scale = phi_scale)
  
  fit <- P0_mod$sample(data = stan_data, chains = 4, parallel_chains = 4, refresh = 500)
  
  return(fit)
}

#Applies Stan diffusion length fit
find_P_fit <- function(freq, spec, dz, alpha_mu = 0.1, alpha_sigma = 1, beta_mu = 1.5, beta_sigma = 1, noise_mu = 0.08, noise_sigma = 0.1, sigma_mu = 0.4, sigma_sigma = 0.4, phi_scale = 1){
  
  stan_data <- list(N = length(freq), freq = freq, spec = spec, dz = dz,
                    alpha_mu = alpha_mu, alpha_sigma = alpha_sigma, beta_mu = beta_mu, beta_sigma = beta_sigma,
                    noise_mu = noise_mu, noise_sigma = noise_sigma, sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, phi_scale = phi_scale)
  
  fit <- P_mod$sample(data = stan_data, chains = 4, parallel_chains = 4, refresh = 500)
  
  return(fit)
}

#Encompassing function for estimating diffusion length
new_diffusion_length_estimate <- function(P0_age, P0_d18O, P_age = MIS19$age, P_depth = MIS19$depth, P_d18O = MIS19$d18O,
                                          range = freq_range){
  
  #Mean annual layer thickness * 1000 over the diffused record. Used for converting between time and depth domains
  bdot <- diff(range(P_depth))/diff(range(P_age))
  
  #Find power law fit of P0
  P0 <- find_spec(x = P0_age, d18O = P0_d18O)
  P0_freq_inds <- which(P0$freq > range$min*bdot & P0$freq < range$max*bdot)
  P0_fit <- find_P0_fit(freq = P0$freq[P0_freq_inds], spec = P0$spec[P0_freq_inds])
  
  alpha_mu_time <- mean(P0_fit$draws("alpha"))
  alpha_sigma_time <- sd(P0_fit$draws("alpha"))
  beta_mu <- mean(P0_fit$draws("beta"))
  beta_sigma <- sd(P0_fit$draws("beta"))
  
  P0_est <- alpha_mu_time*P0$freq[P0_freq_inds]^-beta_mu
  
  #Find diffusion length fit
  P <- find_spec(x = P_depth, d18O = P_d18O)
  P_freq_inds <- which(P$freq > range$min & P$freq < range$max)
  
  alpha_mu <- alpha_mu_time*bdot #Converting prior from time to depth domain
  alpha_sigma <- alpha_sigma_time*bdot #Converting prior from time to depth domain
  
  dz <- P_depth[2] - P_depth[1]
  
  P_fit <- find_P_fit(freq = P$freq[P_freq_inds], spec = P$spec[P_freq_inds], dz = dz, alpha_mu = alpha_mu, alpha_sigma = alpha_sigma, beta_mu = beta_mu, beta_sigma = beta_sigma,
                      sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, noise_mu = noise_mu, noise_sigma = noise_sigma)
  
  alpha_est <- P_fit$draws("alpha")
  beta_est <- P_fit$draws("beta")
  sigma_est <- P_fit$draws("sigma")
  noise_est <- P_fit$draws("noise")
  
  P_est <- mean(alpha_est)*P$freq^-mean(beta_est)*exp(-(2*pi*P$freq*mean(sigma_est))^2) + mean(noise_est)^2*dz
  
  params <- tibble(alpha_ests = as.vector(alpha_est), 
                   beta_ests = as.vector(beta_est),
                   sigma_ests = as.vector(sigma_est),
                   noise_ests = as.vector(noise_est),
                   alpha_mu = alpha_mu,
                   alpha_sigma = alpha_sigma,
                   beta_mu = beta_mu,
                   beta_sigma = beta_sigma,
                   sigma_mu = sigma_mu,
                   sigma_sigma = sigma_sigma,
                   noise_mu = noise_mu,
                   noise_sigma = noise_sigma
  )
  
  P0_spec <- tibble(freq = P0$freq[P0_freq_inds],
                    spec = P0$spec[P0_freq_inds],
                    fit = P0_est)
  
  P0_spec_full <- tibble(freq_full = P0$freq,
                         spec_full = P0$spec)
  
  
  MIS19_fit <- tibble(freq = P$freq,
                      spec = P$spec,
                      fit = P_est)
  
  return(list(params = params, P0_spec = P0_spec, MIS19_fit = MIS19_fit, P0_spec_full = P0_spec_full))
}


# Finding diffusion lengths for individual interglaicials -----------------

#Selecting freq range
MIS19_bdot <- diff(range(MIS19$depth))/diff(range(MIS19$age))
freq_range <- list(min = 0.25/MIS19_bdot, max = 2.5/MIS19_bdot) #Manually chosen and converted to m^-1 units

#Calculating appropriate freq range, gives roughly the same result as manual method
# P <- find_spec(x = MIS19$depth, d18O = MIS19$d18O)
# noise_level <- mean(P$spec[30:121]) #Manually selected from highest frequencies
# P0_level <- mean(P$spec[1:10]) #Manually selected from lowest frequencies
# freq_range <- define_freq_range(min_amp = exp(-1), max_amp = noise_level/P0_level, sigma_min = 0.15, sigma_max = 0.6) #Computed frequency range (closely matches manual method)


#Defining diffusion length and noise priors
sigma_mu <- 0.4
sigma_sigma <- 0.4
noise_mu <- 0.07
noise_sigma <- 0.02

#Running function for each MIS
MIS5_results <- new_diffusion_length_estimate(P0_age = MIS5$age, P0_d18O = MIS5$d18O, range = freq_range)
MIS9_results <- new_diffusion_length_estimate(P0_age = MIS9$age, P0_d18O = MIS9$d18O, range = freq_range)
MIS1_results <- new_diffusion_length_estimate(P0_age = MIS1$age, P0_d18O = MIS1$d18O, range = freq_range)


# Finding diffusion length using mean P0 estimate -------------------------

#Mean alpha prior
mean_alpha_mu <- mean(c(MIS1_results$params$alpha_mu[1], MIS5_results$params$alpha_mu[1], MIS9_results$params$alpha_mu[1]))
mean_alpha_sigma <- sd(c(MIS1_results$params$alpha_mu[1], MIS5_results$params$alpha_mu[1], MIS9_results$params$alpha_mu[1]))

#Mean beta prior
mean_beta_mu <- mean(c(MIS1_results$params$beta_mu[1], MIS5_results$params$beta_mu[1], MIS9_results$params$beta_mu[1]))
mean_beta_sigma <- sd(c(MIS1_results$params$beta_mu[1], MIS5_results$params$beta_mu[1], MIS9_results$params$beta_mu[1]))

P <- find_spec(x = MIS19$depth, d18O = MIS19$d18O)
dz <- MIS19$depth[2] - MIS19$depth[1]
P_freq_inds <- which(P$freq > freq_range$min & P$freq < freq_range$max)

mean_P_fit <- find_P_fit(freq = P$freq[P_freq_inds], spec = P$spec[P_freq_inds], dz = dz, 
                         alpha_mu = mean_alpha_mu, alpha_sigma = mean_alpha_sigma, beta_mu = mean_beta_mu, beta_sigma = mean_beta_sigma,
                         sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, noise_mu = noise_mu, noise_sigma = noise_sigma)

alpha_est <- mean_P_fit$draws("alpha")
beta_est <- mean_P_fit$draws("beta")
sigma_est <- mean_P_fit$draws("sigma")
noise_est <- mean_P_fit$draws("noise")

P_est <- mean(alpha_est)*P$freq^-mean(beta_est)*exp(-(2*pi*P$freq*mean(sigma_est))^2) + mean(noise_est)^2*dz

mean_params <- tibble(MIS = "mean", alpha_ests = as.vector(alpha_est), beta_ests = as.vector(beta_est), sigma_ests = as.vector(sigma_est), noise_ests = as.vector(noise_est),
                      alpha_mu = mean_alpha_mu, alpha_sigma = mean_alpha_sigma, beta_mu = mean_beta_mu, beta_sigma = mean_beta_sigma, sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, noise_mu = noise_mu, noise_sigma = noise_sigma)

mean_P0_spec <- tibble(freq = P$freq[P_freq_inds]*MIS19_bdot,
                       spec = (mean_alpha_mu/MIS19_bdot)*(P$freq[P_freq_inds]*MIS19_bdot)^-mean_beta_mu,
                       fit = (mean_alpha_mu/MIS19_bdot)*(P$freq[P_freq_inds]*MIS19_bdot)^-mean_beta_mu)

mean_MIS19_fit <- tibble(freq = P$freq,
                         spec = P$spec,
                         fit = P_est)


# Data manipulation for plotting ------------------------------------------

#Combining all parameter results
MIS_params <- bind_rows(MIS1 = MIS1_results$params, MIS5 = MIS5_results$params, MIS9 = MIS9_results$params, mean = mean_params, .id = "MIS")
MIS_params$MIS <- factor(MIS_params$MIS, levels = c('MIS1', 'MIS5', 'MIS9', 'mean'))
MIS_params_long <- MIS_params %>% 
  select(alpha_ests, beta_ests, sigma_ests, noise_ests, MIS) %>% 
  pivot_longer(., cols = c(alpha_ests, beta_ests, sigma_ests, noise_ests), names_to = "parameter")

#Table of average parameter values
MIS_params_average <- aggregate(MIS_params, list(MIS_params$MIS), mean) %>%
  select(Group.1, alpha_ests, beta_ests, sigma_ests, noise_ests)

#Alpha estimates and prior data frame
prior_df_alpha <- MIS_params %>% 
  group_by(MIS) %>% 
  select(alpha_mu, alpha_sigma, MIS) %>% 
  distinct() %>% 
  summarise(
    parameter = "alpha_ests",
    x = seq(0, 0.05, length.out = 1000),
    d = dnorm(x, alpha_mu, alpha_sigma)
  )

#Beta estimates and prior data frame
prior_df_beta <- MIS_params %>% 
  group_by(MIS) %>% 
  select(beta_mu, beta_sigma, MIS) %>% 
  distinct() %>% 
  summarise(
    parameter = "beta_ests",
    x = seq(0.5, 3, length.out = 1000),
    d = dnorm(x, beta_mu, beta_sigma)
  )

#Sigma estimates and prior data frame
prior_df_sigma <- MIS_params %>% 
  group_by(MIS) %>% 
  select(sigma_mu, sigma_sigma, MIS) %>% 
  distinct() %>% 
  summarise(
    parameter = "sigma_ests",
    x = seq(0.1, 0.5, length.out = 1000),
    d = dnorm(x, sigma_mu, sigma_sigma)
  )

#Noise estimates and prior data frame
prior_df_noise <- MIS_params %>% 
  group_by(MIS) %>% 
  select(noise_mu, noise_sigma, MIS) %>% 
  distinct() %>% 
  summarise(
    parameter = "noise_ests",
    x = seq(0.055, 0.09, length.out = 1000),
    d = dnorm(x, noise_mu, noise_sigma)
  )

#Combining parameter estimates and prior data frames
prior_df <- bind_rows(prior_df_alpha, prior_df_beta, prior_df_sigma, prior_df_noise)

#Windowing and combining P0 fits in a data frame
P0_fits <- bind_rows(MIS1 = MIS1_results$P0_spec, MIS5 = MIS5_results$P0_spec, MIS9 = MIS9_results$P0_spec, mean = mean_P0_spec, .id = "MIS")
P0_fits$MIS <- factor(P0_fits$MIS, levels = c('MIS1', 'MIS5', 'MIS9', 'mean'))

#Combining MIS 19 fits in a data frame
MIS19_fits <- bind_rows(MIS5 = MIS5_results$MIS19_fit, MIS9 = MIS9_results$MIS19_fit, MIS1 = MIS1_results$MIS19_fit, mean = mean_MIS19_fit, .id = "MIS")
MIS19_fits$MIS <- factor(MIS19_fits$MIS, levels = c('MIS1', 'MIS5', 'MIS9', 'mean'))

#Small data frame for new MIS 19 estimate in Dome C sigmas plot
new_est_with_errors <- tibble(depth = range(MIS19$depth), mean = rep(mean(mean_params$sigma_ests), 2),
                              min = rep(mean(mean_params$sigma_ests) - 2*sd(mean_params$sigma_ests), 2), 
                              max = rep(mean(mean_params$sigma_ests) + 2*sd(mean_params$sigma_ests), 2))


# MIS 19 diffusion length estimate using the conventional method -----------------

P_conv_mod <- cmdstan_model("conventional_fit.stan")

data <- list(N = length(P$freq), freq = P$freq, spec = P$spec, dz = 0.11,
             alpha_mu = 1, alpha_sigma = 10, noise_mu = 0.07, noise_sigma = 0.02, sigma_mu = 0.4, sigma_sigma = 0.4, phi_scale = 1)

conventional_fit <- P_conv_mod$sample(data = data, chains = 4, parallel_chains = 4, refresh = 500)

conventional_results <- tibble(alpha = as.vector(conventional_fit$draws('alpha')),
                               sigma = as.vector(conventional_fit$draws('sigma')),
                               noise = as.vector(conventional_fit$draws('noise')))

conventional_P_fit <- mean(conventional_results$alpha) * exp(-(2*pi*P$freq*mean(conventional_results$sigma))^2) + mean(conventional_results$noise)^2*0.1






# Paper figures -----------------------------------------------------------

#For figure 2 (P0 fits)
P0_x_min <- min(MIS1_results$P0_spec_full$freq_full, MIS5_results$P0_spec_full$freq_full, MIS9_results$P0_spec_full$freq_full)
P0_x_max <- max(MIS1_results$P0_spec_full$freq_full, MIS5_results$P0_spec_full$freq_full, MIS9_results$P0_spec_full$freq_full)
P0_y_min <- min(MIS1_results$P0_spec_full$spec_full, MIS5_results$P0_spec_full$spec_full, MIS9_results$P0_spec_full$spec_full)
P0_y_max <- max(MIS1_results$P0_spec_full$spec_full, MIS5_results$P0_spec_full$spec_full, MIS9_results$P0_spec_full$spec_full)

#Need average 1000 year layer thickness to create extra x axis on spectra plots
MIS1_bdot <- diff(range(MIS1$depth))/diff(range(MIS1$age))
MIS5_bdot <- diff(range(MIS5$depth))/diff(range(MIS5$age))
MIS9_bdot <- diff(range(MIS9$depth))/diff(range(MIS9$age))
MIS19_bdot <- diff(range(MIS19$depth))/diff(range(MIS19$age))

#Colour palette for figures
colfun <- brewer.pal(n = 6, name = 'Dark2') #Dark2 or Set2?

#Figure 1: Dome C d18O timeseries with highlighted interglacials
pdf(file = "figures/fig_1.pdf",
    width = 10,
    height = 4)

plot(DomeC$age, DomeC$d18O, 'l', lwd = 1, xlim = c(28, 775), main = expression(paste(bold("Dome C "), bold(delta^{18}), bold("O record"))), xlab = "age (ka)", ylab = expression(paste(delta^{18}, "O")),
     panel.first = rect(xleft = c(MIS1$age[1], MIS5$age[1], MIS9$age[1], MIS19$age[1]),
                        ybottom = -70,
                        xright = c(MIS1$age[length(MIS1$age)], MIS5$age[length(MIS5$age)], MIS9$age[length(MIS9$age)], MIS19$age[length(MIS19$age)]),
                        ytop = -40,
                        col = paste(colfun[1:4], "70", sep = ""),
                        border = NA
     ))
legend(660, -43.2, c("MIS 1", "MIS 5", "MIS 9", "MIS 19"), col = colfun[1:4], fill = paste(colfun[1:4], "70", sep = ""), cex = 1,  bty = "n", x.intersp = 0.3, y.intersp = 1, seg.len = 0.7)

dev.off()

#Figure 2: P0 fits on a) MIS 1, b) MIS 5 and c) MIS 9 spectra

fig_2a <- ggplot(data = MIS1_results$P0_spec) +
  geom_rect(aes(xmin = freq[1], xmax = freq[length(freq)], ymin = 0, ymax = Inf, fill = 'Fitted region'), alpha = 1) +
  geom_line(data = MIS1_results$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "Full spectrum")) +
  geom_line(aes(x = freq, y = fit, colour = "Fit"), lwd = 1) +
  scale_x_continuous(trans = 'log', limits = c(P0_x_min, P0_x_max), breaks = 10^(-1:2), sec.axis = sec_axis(~ . / MIS1_bdot, name = 'Frequency (1/m)', breaks = 10^(-2:1))) +
  scale_y_continuous(trans = 'log', limits = c(P0_y_min, P0_y_max), breaks = 10^(-8:2)) +
  scale_colour_manual(name = '', values = c('black', colfun[1]), breaks = c('Full spectrum', 'Fit')) +
  scale_fill_manual(name = '', values = 'grey92') +
  labs(title = 'MIS 1', tag = 'a)') +
  xlab("Frequency (1/kyr)") +
  ylab("PSD") +
  theme_bw() +
  theme(legend.position = "none",
        plot.tag = element_text())

fig_2b <- ggplot(data = MIS5_results$P0_spec) +
  geom_rect(aes(xmin = freq[1], xmax = freq[length(freq)], ymin = 0, ymax = Inf, fill = 'Fitted region'), alpha = 1) +
  geom_line(data = MIS5_results$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "Full spectrum")) +
  geom_line(aes(x = freq, y = fit, colour = "Fit"), lwd = 1) +
  scale_x_continuous(trans = 'log', limits = c(P0_x_min, P0_x_max), breaks = 10^(-1:2), sec.axis = sec_axis(~ . / MIS5_bdot, name = 'Frequency (1/m)', breaks = 10^(-2:1))) +
  scale_y_continuous(trans = 'log', limits = c(P0_y_min, P0_y_max), breaks = 10^(-8:2)) +
  scale_colour_manual(name = '', values = c('black', colfun[2]), breaks = c('Full spectrum', 'Fit')) +
  scale_fill_manual(name = '', values = 'grey92') +
  labs(title = 'MIS 5', tag = 'b)') +
  xlab("Frequency (1/kyr)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.tag = element_text())

#Fake points outside the axes limits, created so the data from the first two plots appears in the third plot's legend
fake1 <- tibble(x = 1000, y = 1000)
fake2 <- tibble(x = 1000, y = 1000)

fig_2c <- ggplot(data = MIS9_results$P0_spec) +
  geom_rect(aes(xmin = freq[1], xmax = freq[length(freq)], ymin = 0, ymax = Inf, fill = 'Fitted region'), alpha = 1) +
  geom_line(data = MIS9_results$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "Full spectrum")) +
  geom_line(aes(x = freq, y = fit, colour = "MIS 9 fit"), lwd = 1) +
  geom_line(data = fake1, aes(x, y, colour = "MIS 1 fit")) +
  geom_line(data = fake1, aes(x, y, colour = "MIS 5 fit")) +
  scale_x_continuous(trans = 'log', limits = c(P0_x_min, P0_x_max), breaks = 10^(-1:2), sec.axis = sec_axis(~ . / MIS9_bdot, name = 'Frequency (1/m)', breaks = 10^(-2:1))) +
  scale_y_continuous(trans = 'log', limits = c(P0_y_min, P0_y_max), breaks = 10^(-8:2)) +
  scale_colour_manual(name = '', values = c('black', colfun[1:3]), breaks = c('Full spectrum', 'MIS 1 fit', 'MIS 5 fit', 'MIS 9 fit')) +
  scale_fill_manual(name = '', values = 'grey92') +
  labs(title = 'MIS 9', tag = 'c)') +
  xlab("Frequency (1/kyr)") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.tag = element_text())

pdf(file = "figures/fig_2.pdf",
    width = 12,
    height = 4)

grid::grid.draw(cbind(ggplotGrob(fig_2a), ggplotGrob(fig_2b), ggplotGrob(fig_2c), size = "last"))

dev.off()

#Figure 3: a) MIS 19 fits and b) diffusion length histograms

fig_3a <- ggplot() + 
  geom_rect(data = MIS19_fits[which(MIS19_fits$freq > freq_range$min & MIS19_fits$freq < freq_range$max), ], aes(xmin = freq[1], xmax = freq[length(freq)], ymin = 0, ymax = Inf, fill = 'Fitted region'), alpha = 1) +
  geom_line(data = MIS19_fits, aes(x = freq, y = spec, linetype = 'MIS 19 spectrum')) +
  geom_line(data = MIS19_fits, aes(x = freq, y = fit, col = MIS)) +
  scale_x_continuous(trans = 'log', breaks = 10^(-1:1), sec.axis = sec_axis(~ . / MIS19_bdot, name = 'Frequency (1/m)', breaks = 10^(-2:1))) +
  scale_y_continuous(trans = 'log', breaks = 10^(-8:2)) +
  scale_linetype_manual(name = '', values = 1) +
  scale_colour_manual(name = 'MIS', values = colfun[1:4]) +
  scale_fill_manual(name = '', values = 'grey92') +
  guides(linetype = guide_legend(order = 1), col = guide_legend(order = 2), fill = guide_legend(order = 3)) +
  labs(tag = 'a)') +
  xlab("Frequency (1/kyr)") +
  ylab("PSD") +
  theme_bw()

fig_3b <- ggplot(MIS_params) +
  geom_histogram(alpha = 0.5, bins = 50, aes(x = sigma_ests, y = ..density.., fill = MIS), position = 'identity') +
  geom_line(data = prior_df_sigma, aes(x = x, y = d, size = 'Prior')) +
  scale_fill_manual(name = 'MIS', values = colfun[1:4]) +
  scale_size_manual(name = '', values = 0.5) +
  scale_x_continuous(limits = c(0.18, 0.42)) +
  labs(tag = 'b)') +
  xlab("Diffusion length (m)") +
  theme_bw()

pdf(file = "figures/fig_3.pdf",
    width = 12,
    height = 4)

grid.arrange(fig_3a, fig_3b, nrow = 1)

dev.off()

#Figure 4a & 4b: Bias of conventional estimator on diffused power law spectra
#NOTE!!!: This figure requires set "bias_data" from another script: conventional_bias.R

bias_data <- read_csv('bias_data.csv')

fig_4a <- ggplot(data = bias_data, aes(x = freq)) +
  geom_line(aes(y = P_sim, col = 'Simulated spectrum'), col = 'black', lwd = 1) +
  geom_line(aes(y = P_conv_est, col = 'Conventional fit'), col = colfun[5], lwd = 1) +
  geom_line(aes(y = P_new_est, col = 'New fit'), linetype = 2, col = colfun[4], lwd = 1) +
  scale_x_continuous(trans = 'log', breaks = 10^(-1:1)) +
  scale_y_continuous(trans = 'log', breaks = 10^(-3:2)) +
  labs(tag = 'a)') +
  xlab("Frequency (arbitrary units)") +
  ylab("PSD (arbitrary units)") +
  theme_bw()

fig_4b <- ggplot(data = bias_data, aes(x = freq^2)) +
  geom_line(aes(y = P_sim, col = 'Simulated spectrum'), lwd = 1) +
  geom_line(aes(y = P_conv_est, col = 'Conventional fit'), lwd = 1) +
  geom_line(aes(y = P_new_est, col = 'New fit'), linetype = 2, lwd = 1) +
  scale_x_continuous(breaks = (0:3)*10, limits = c(0, 25)) +
  scale_y_continuous(trans = 'log', breaks = 10^(-3:2)) +
  scale_colour_manual(name = '', values = c('black', colfun[5], colfun[4]), breaks = c('Simulated spectrum', 'Conventional fit', 'New fit')) +
  labs(tag = 'b)') +
  xlab("Frequency^2 (arbitrary units)") +
  ylab("PSD (arbitrary units)") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

pdf(file = "figures/fig_4.pdf",
    width = 10,
    height = 4)

grid::grid.draw(cbind(ggplotGrob(fig_4a), ggplotGrob(fig_4b), size = "last"))

dev.off()

#Figure 5: a) New vs conventional MIS 19 fit b) Conventional diffusion length estimate depth profile for EDC and new estimate for comparison
#NOTE!!!: This figure requires data set "conv_sigmas" from another script: "DomeC_diffusion.R"

conv_sigmas <- read_csv('conv_sigmas.csv')

fig_5a <- ggplot(subset(MIS19_fits, MIS == 'mean')) +
  geom_line(aes(x = freq, y = spec, col = 'MIS 19')) +
  geom_line(aes(x = freq, y = fit, col = 'New fit'), size = 0.7) +
  geom_line(aes(x = freq, y = conventional_P_fit, col = 'Conventional fit'), size = 0.7) +
  scale_colour_manual(values = c('black', colfun[5], colfun[4]), breaks = c('MIS 19', 'Conventional fit', 'New fit')) +
  scale_x_continuous(trans = 'log', breaks = 10^(-1:1)) +
  scale_y_continuous(trans = 'log', breaks = 10^(-8:2)) +
  labs(tag = 'a)') +
  xlab("Frequency (1/m)") +
  ylab("PSD") +
  labs(col = "") +
  theme_bw()
# ^+ theme(legend.position = c(0.85, 0.9))

fig_5b_main <- ggplot(data = conv_sigmas) +
  geom_line(aes(x = depth, y = sigma_mean, colour = 'Conventional')) +
  geom_ribbon(aes(x = depth, ymin = sigma_min, ymax = sigma_max), colour = NA, fill = colfun[5], alpha = 0.3) +
  geom_point(aes(x = mean(MIS19$depth), y = mean(mean_params$sigma_ests), colour = 'New')) +
  geom_errorbar(aes(x = mean(MIS19$depth), ymin = mean(mean_params$sigma_ests) - 2*sd(mean_params$sigma_ests), ymax = mean(mean_params$sigma_ests) + 2*sd(mean_params$sigma_ests), colour = 'New'), width = 50) +
  geom_point(aes(x = mean(3147, 3190), y = 0.5, colour = 'Pol')) +
  geom_errorbar(aes(x = mean(3147, 3190), ymin = 0.4, ymax = 0.6, colour = 'Pol'), width = 50) +
  scale_colour_manual(name = 'Estimate', values = c(colfun[5], colfun[4], colfun[6])) +
  labs(tag = 'b)') +
  xlab("Depth (m)") +
  ylab("Diffusion length (m)") +
  theme_bw()

fig_5b_inset <- ggplot(data = conv_sigmas) +
  geom_line(aes(x = depth, y = sigma_mean), colour = colfun[5]) +
  geom_ribbon(aes(x = depth, ymin = sigma_min, ymax = sigma_max), colour = NA, fill = colfun[5], linetype = 1, alpha = 0.3) +
  geom_point(aes(x = mean(MIS19$depth), y = mean(mean_params$sigma_ests)), colour = colfun[4]) +
  geom_errorbar(aes(x = mean(MIS19$depth), ymin = mean(mean_params$sigma_ests) - 2*sd(mean_params$sigma_ests), ymax = mean(mean_params$sigma_ests) + 2*sd(mean_params$sigma_ests)), colour = colfun[4], width = 30) +
  geom_point(aes(x = mean(3147, 3190), y = 0.5),  colour = colfun[6]) +
  geom_errorbar(aes(x = mean(3147, 3190), ymin = 0.4, ymax = 0.6), colour = colfun[6], width = 30) +
  xlim(2800, 3200) +
  ylim(0.15, 0.8) +
  theme_bw() +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank())

fig_5b <- fig_5b_main +
  annotation_custom(ggplotGrob(fig_5b_inset), xmin = 0, xmax = 2500, ymin = 0.3, ymax = 0.75)

pdf(file = "figures/fig_5.pdf",
    width = 12,
    height = 4)

grid.arrange(fig_5a, fig_5b, nrow = 1)

dev.off()




