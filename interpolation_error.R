#Testing effect of equidistant age and gap interpolation on spectra

#Must run MIS19_diffusion_length.R first!

library(PaleoSpec)

setwd('/Users/fshaw/Desktop/Diffusion_Length_Estimate_Manuscript/Final Submission/Final Code')

#Function for simulating a timeseries with the same gaps spectral parameters, gaps and non-equidistant age as input data
sim_real <- function(data, alpha, beta){
  
  dt <- mean(diff(data$age))
  
  data_sim <- SimProxySeries(a = alpha, b = beta, nt = (length(data$age) - 1)*10, var.noise = 0.07^2)
  age_sim <- seq(data$age[1], data$age[length(data$age)] - dt/10, by = dt/10)
  
  breaks_even <- seq(data$age[1], data$age[length(data$age)] - dt, by = dt)
  breaks_uneven <- data$age[-1] - diff(data$age)/2
  
  data_even <- AvgToBin(age_sim, data_sim, breaks = breaks_even)
  data_uneven <- AvgToBin(age_sim, data_sim, breaks = breaks_uneven)
  
  gaps <- is.na(data$d18O)
  data_gappy <- data_uneven
  data_gappy$avg[gaps[-((length(gaps) - 1):length(gaps))]] <- NA
  
  data_interp <- approx(data_gappy$centers, data_gappy$avg, data_even$centers, rule = 2)
  
  return(list(even = list(age = data_even$centers, d18O = data_even$avg),
              gappy = list(age = data_gappy$centers, d18O = data_gappy$avg),
              interp = list(age = data_interp$x, d18O = data_interp$y)))
  
}

#Function to apply power law fit to input data, where alpha and beta are the true values 
spec_fit <- function(age, d18O){
  
  dt <- mean(diff(age))
  
  ts <- ts(d18O, start = age[1], deltat = mean(diff(age)))
  sp <- spectrum(ts, plot = F, spans = NULL, taper = 0.1, pad = 0, fast = T, detrend = T)
  
  bdot <- diff(range(MIS19$depth))/diff(range(MIS19$age))
  freq_range <- list(min = 0.25/MIS19_bdot, max = 2.5/MIS19_bdot)
  
  freq_inds <- which(sp$freq > freq_range$min*bdot & sp$freq < freq_range$max*bdot)
  
  fit <- find_P0_fit(freq = sp$freq[freq_inds], sp$spec[freq_inds], alpha_mu = 0.1, alpha_sigma = 1, beta_mu = 1.5, beta_sigma = 1)
  
  P0_fit <- mean(fit$draws('alpha'))*sp$freq[freq_inds]^-mean(fit$draws('beta'))
  
  P0_fit_min <- fit$summary()$q5[-(1:4)]
  P0_fit_max <- fit$summary()$q95[-(1:4)]
  
  P0_fit_plot <- tibble(freq = sp$freq[freq_inds],
                        spec = sp$spec[freq_inds], 
                        fit = P0_fit,
                        fit_min = P0_fit_min,
                        fit_max = P0_fit_max)
  
  P0_spec_full <- tibble(freq_full = sp$freq,
                         spec_full = sp$spec)
  
  return(list(P0_fit_plot = P0_fit_plot, P0_spec_full = P0_spec_full, fit_summary = fit$summary()))
  
}

#Function encompassing previous two functions, and outputting both the "true" equidistant result and the interpolated result
spec_error <- function(data, alpha, beta){
  
  sim <- sim_real(data, alpha, beta)
  
  spec_even <- spec_fit(age = sim$even$age, d18O = sim$even$d18O)
  spec_interp <- spec_fit(age = sim$interp$age, d18O = sim$interp$d18O)
  
  return(list(spec_even = spec_even, spec_interp = spec_interp, spec_sim = sim$spec_sim))
}

#Finding MIS 1/5/9 alpha and beta values for simulated data
#This requires the results from the main script "MIS19_diffusion_length.R"
MIS1_beta <- mean(MIS1_results$P0_params$beta_ests)
MIS1_alpha <- mean(MIS1_results$P0_params$alpha_ests)/(mean(diff(MIS1$age))^(1-MIS1_beta))

MIS5_beta <- mean(MIS5_results$P0_params$beta_ests)
MIS5_alpha <- mean(MIS5_results$P0_params$alpha_ests)/(mean(diff(MIS5$age))^(1-MIS5_beta))

MIS9_beta <- mean(MIS9_results$P0_params$beta_ests)
MIS9_alpha <- mean(MIS9_results$P0_params$alpha_ests)/(mean(diff(MIS9$age))^(1-MIS9_beta))

#Applied to each interglacial
MIS1_error <- spec_error(data = MIS1, alpha = MIS1_alpha, beta = MIS1_beta)
MIS5_error <- spec_error(data = MIS5, alpha = MIS5_alpha , beta = MIS5_beta)
MIS9_error <- spec_error(data = MIS9, alpha = MIS9_alpha , beta = MIS9_beta)

#Plot results
fig_A1a <- ggplot(data = MIS1_error$spec_even$P0_fit_plot) +
  geom_rect(aes(xmin = freq[1], xmax = freq[length(freq)], ymin = 0, ymax = Inf, fill = 'Fitted region'), alpha = 1) +
  geom_line(data = MIS1_error$spec_even$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "True spectrum")) +
  geom_line(aes(x = freq, y = fit, colour = "True Fit"), lwd = 1) +
  geom_line(data = MIS1_error$spec_interp$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "Interpolated spectrum")) +
  geom_line(data = MIS1_error$spec_interp$P0_fit_plot, aes(x = freq, y = fit, colour = "Interpolated Fit"), lwd = 1) +
  geom_ribbon(aes(x = freq, ymin = fit_min, ymax = fit_max), colour = NA, fill = 'grey50', alpha = 0.4) +
  scale_x_continuous(trans = 'log', limits = c(P0_x_min, P0_x_max), breaks = 10^(-1:2)) +#, sec.axis = sec_axis(~ . / MIS1_bdot, name = 'Frequency (1/m)', breaks = 10^(-2:1))) +
  scale_y_continuous(trans = 'log', limits = c(P0_y_min, P0_y_max), breaks = 10^(-8:2)) +
  scale_colour_manual(name = '', values = c('black', 'grey50', colfun[1], 'green2'), breaks = c('True spectrum', 'True Fit', 'Interpolated spectrum', 'Interpolated Fit')) +
  scale_fill_manual(name = '', values = 'grey92') +
  labs(title = 'a) MIS 1 parameters') +
  xlab("Frequency (1/kyr)") +
  ylab("PSD") +
  theme_bw() +
  theme(legend.position = "none",
        plot.tag = element_text())

fig_A1b <- ggplot(data = MIS5_error$spec_even$P0_fit_plot) +
  geom_rect(aes(xmin = freq[1], xmax = freq[length(freq)], ymin = 0, ymax = Inf, fill = 'Fitted region'), alpha = 1) +
  geom_line(data = MIS5_error$spec_even$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "True spectrum")) +
  geom_line(aes(x = freq, y = fit, colour = "True Fit"), lwd = 1) +
  geom_line(data = MIS5_error$spec_interp$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "Interpolated spectrum")) +
  geom_line(data = MIS5_error$spec_interp$P0_fit_plot, aes(x = freq, y = fit, colour = "Interpolated Fit"), lwd = 1) +
  geom_ribbon(aes(x = freq, ymin = fit_min, ymax = fit_max), colour = NA, fill = 'grey50', alpha = 0.4) +
  scale_x_continuous(trans = 'log', limits = c(P0_x_min, P0_x_max), breaks = 10^(-1:2)) +#, sec.axis = sec_axis(~ . / MIS5_bdot, name = 'Frequency (1/m)', breaks = 10^(-2:1))) +
  scale_y_continuous(trans = 'log', limits = c(P0_y_min, P0_y_max), breaks = 10^(-8:2)) +
  scale_colour_manual(name = '', values = c('black', 'grey50', colfun[2], 'orange2'), breaks = c('True spectrum', 'True Fit', 'Interpolated spectrum', 'Interpolated Fit')) +
  scale_fill_manual(name = '', values = 'grey92') +
  labs(title = 'b) MIS 5 parameters') +
  xlab("Frequency (1/kyr)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.tag = element_text())

#Fake points outside the axes limits, created so the data from the first two plots appears in the third plot's legend
fake <- tibble(x = 1000, y = 1000)

fig_A1c <- ggplot(data = MIS9_error$spec_even$P0_fit_plot) +
  geom_rect(aes(xmin = freq[1], xmax = freq[length(freq)], ymin = 0, ymax = Inf, fill = 'Fitted region'), alpha = 1) +
  geom_line(data = MIS9_error$spec_even$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "True spectrum")) +
  geom_line(aes(x = freq, y = fit, colour = "True Fit"), lwd = 1) +
  geom_line(data = MIS9_error$spec_interp$P0_spec_full, aes(x = freq_full, y = spec_full, colour = "MIS 9 interpolated spectrum")) +
  geom_line(data = MIS9_error$spec_interp$P0_fit_plot, aes(x = freq, y = fit, colour = "MIS 9 interpolated fit"), lwd = 1) +
  geom_line(data = fake, aes(x, y, colour = "MIS 1 interpolated spectrum")) +
  geom_line(data = fake, aes(x, y, colour = "MIS 1 interpolated fit"), lwd = 1) +
  geom_line(data = fake, aes(x, y, colour = "MIS 5 interpolated spectrum")) +
  geom_line(data = fake, aes(x, y, colour = "MIS 5 interpolated fit"), lwd = 1) +
  geom_ribbon(aes(x = freq, ymin = fit_min, ymax = fit_max), colour = NA, fill = 'grey50', alpha = 0.4) +
  scale_x_continuous(trans = 'log', limits = c(P0_x_min, P0_x_max), breaks = 10^(-1:2)) +#, sec.axis = sec_axis(~ . / MIS1_bdot, name = 'Frequency (1/m)', breaks = 10^(-2:1))) +
  scale_y_continuous(trans = 'log', limits = c(P0_y_min, P0_y_max), breaks = 10^(-8:2)) +
  scale_colour_manual(name = '', values = c('black', 'grey50', colfun[1:3], 'green2', 'orange2', 'purple2'), breaks = c('True spectrum', 'True Fit', 'MIS 1 interpolated spectrum', 'MIS 5 interpolated spectrum', 'MIS 9 interpolated spectrum', 'MIS 1 interpolated fit', 'MIS 5 interpolated fit', 'MIS 9 interpolated fit')) +
  scale_fill_manual(name = '', values = 'grey92') +
  labs(title = 'c) MIS 9 parameters') +
  xlab("Frequency (1/kyr)") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.tag = element_text())

pdf(file = "./fig_A1.pdf",
    width = 12,
    height = 4)

grid::grid.draw(cbind(ggplotGrob(fig_A1a), ggplotGrob(fig_A1b), ggplotGrob(fig_A1c), size = "last"))

dev.off()
