#Simulating a power law "climate" spectrum and fitting with conventional method and new method

library(cmdstanr)

setwd('/Users/fshaw/Desktop/Diffusion_Length_Estimate_Manuscript/Final Submission/Final Code')

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


#Applies Stan diffusion length fit
find_P_fit <- function(freq, spec, dz, alpha_mu = 0.1, alpha_sigma = 1, beta_mu = 1.5, beta_sigma = 1, noise_mu = 0.08, noise_sigma = 0.1, sigma_mu = 0.4, sigma_sigma = 0.4, phi_scale = 1){
  
  stan_data <- list(N = length(freq), freq = freq, spec = spec, dz = dz,
                    alpha_mu = alpha_mu, alpha_sigma = alpha_sigma, beta_mu = beta_mu, beta_sigma = beta_sigma,
                    noise_mu = noise_mu, noise_sigma = noise_sigma, sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, phi_scale = phi_scale)
  
  fit <- P_mod$sample(data = stan_data, chains = 4, parallel_chains = 4, refresh = 500, save_warmup = TRUE)
  
  return(fit)
}


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

P_mod <- cmdstan_model("./diffusion_length_fit.stan")

alpha_mu <- 1
alpha_sigma <- 0.1
beta_mu <- 1
beta_sigma <- 2
sigma_mu <- 0.4
sigma_sigma <- 0.4
noise_mu <- 0.08
noise_sigma <- 0.1

P_fit <- find_P_fit(freq = freq[2:501], spec = P_sim[2:501], dz = dz, alpha_mu = alpha_mu, alpha_sigma = alpha_sigma, beta_mu = beta_mu, beta_sigma = beta_sigma,
                    sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, noise_mu = noise_mu, noise_sigma = noise_sigma)

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

conv_mod <- cmdstan_model("./conventional_fit.stan")

conv_data <- list(N = length(freq[2:501]), freq = freq[2:501], spec = P_sim[2:501], dz = dz,
                  alpha_mu = alpha_mu, alpha_sigma = alpha_sigma,
                  noise_mu = noise_mu, noise_sigma = noise_sigma, sigma_mu = sigma_mu, sigma_sigma = sigma_sigma, phi_scale = 1)

conv_fit <- conv_mod$sample(data = conv_data, chains = 4, parallel_chains = 4, refresh = 500)

alpha_conv_est <- mean(conv_fit$draws('alpha'))
sigma_conv_est <- mean(conv_fit$draws('sigma'))
noise_conv_est <- mean(conv_fit$draws('noise'))

P_conv_est <- alpha_conv_est * exp(-(2*pi*freq[2:501]*sigma_conv_est)^2) + noise_conv_est^2*dz

results <- tibble(freq = freq[2:501],
                  P_sim = P_sim[2:501],
                  P_new_est = P_est,
                  P_conv_est = P_conv_est)

write.csv(results, "./conventional_bias.csv")




