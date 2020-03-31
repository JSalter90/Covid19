#### Predicting future number of cases ####
library(ukcovidtools)
library(ggplot2)
library(EpiEstim)
library(tidyverse)

# Need:
# Number of cases I_0, ..., I_t up to current time t
# R(t) at current time (likely a distribution)
# w, distribution for incubation period/serial interval (likely have uncertainty on both parameters of this distribution)

# Data
ts = ukcovidtools::getUKCovidTimeseries() # contains count data by region, UA, against time
load("R0timeseries.RData") # posterior for R(t) by UA

# Distribution for w
serialIntervals = tibble(
  mean_si_estimate = c(3.96, 6.3, 4.22, 4.56, 3.95, 5.21, 4.7, 7.5,6.6),
  mean_si_estimate_low_ci = c(3.53, 5.2, 3.43, 2.69,-4.47, -3.35, 3.7, 5.3, 0.7),
  mean_si_estimate_high_ci = c(4.39, 7.6, 5.01, 6.42, 12.51,13.94, 6.0, 19.0, 19.0),
  std_si_estimate = c(4.75,4.2, 0.4, 0.95, 4.24, 4.32, 2.3, 3.4, NA),
  std_si_estimate_low_ci = c(4.46, 3.1, NA, NA, 4.03, 4.06, 1.6, NA, NA),
  std_si_estimate_high_ci = c(5.07, 5.3, NA, NA, 4.95, 5.58, 3.5, NA, NA),
  sample_size = c(468,48,135,93,45,54,28,16,90),
  population = c("China", "Shenzhen","Taijin","Singapore","Taijin","Singapore", "SE Asia", "Wuhan","Italy"),
  source = c(
    "Zhanwei Du et al. Serial Interval of COVID-19 among Publicly Reported Confirmed Cases. Emerging Infectious Disease journal 26, (2020)",
    "Bi, Q. et al. Epidemiology and Transmission of COVID-19 in Shenzhen China: Analysis of 391 cases and 1,286 of their close contacts. Infectious Diseases (except HIV/AIDS) (2020) doi:10.1101/2020.03.03.20028423",
    "Tindale, L. et al. Transmission interval estimates suggest pre-symptomatic spread of COVID-19. Epidemiology (2020) doi:10.1101/2020.03.03.20029983",
    "Tindale, L. et al. Transmission interval estimates suggest pre-symptomatic spread of COVID-19. Epidemiology (2020) doi:10.1101/2020.03.03.20029983",
    "Ganyani, T. et al. Estimating the generation interval for COVID-19 based on symptom onset data. Infectious Diseases (except HIV/AIDS) (2020) doi:10.1101/2020.03.05.20031815",
    "Ganyani, T. et al. Estimating the generation interval for COVID-19 based on symptom onset data. Infectious Diseases (except HIV/AIDS) (2020) doi:10.1101/2020.03.05.20031815",
    "Nishiura, H., Linton, N. M. & Akhmetzhanov, A. R. Serial interval of novel coronavirus (COVID-19) infections. Int. J. Infect. Dis. (2020) doi:10.1016/j.ijid.2020.02.060",
    "Li, Q. et al. Early Transmission Dynamics in Wuhan, China, of Novel Coronavirus-Infected Pneumonia. N. Engl. J. Med. (2020) doi:10.1056/NEJMoa2001316",
    "Cereda, D. et al. The early phase of the COVID-19 outbreak in Lombardy, Italy. arXiv [q-bio.PE] (2020)")
)

wtSIs = serialIntervals %>% summarise(
  mean_si = weighted.mean(mean_si_estimate,sample_size,na.rm = TRUE),
  min_mean_si = weighted.mean(mean_si_estimate_low_ci,sample_size,na.rm = TRUE),
  max_mean_si = weighted.mean(mean_si_estimate_high_ci,sample_size,na.rm = TRUE),
  std_si  = weighted.mean(ifelse(is.na(std_si_estimate_low_ci),NA,1)*std_si_estimate,sample_size,na.rm = TRUE),
  min_std_si  = weighted.mean(std_si_estimate_low_ci,sample_size,na.rm = TRUE),
  max_std_si  = weighted.mean(std_si_estimate_high_ci,sample_size,na.rm = TRUE)
  #total = sum(sample_size)
) %>% mutate(
  std_mean_si = (max_mean_si - min_mean_si) / 3.92, # TODO: fit gamma
  std_std_si = (max_std_si - min_std_si) / 3.92
)

#### Now just need to simulate from everything ####
# Select a single location
UA <- 'Devon'
UA_data <- subset(R0timeseries, GSS_NM == UA)
I_ts <- UA_data$incidence
plot(UA_data$date, I_ts, pch = 19, xlab = 'Date', ylab = 'Total cases')

# Find R(t) at current time t
Rt <- UA_data$`Median(R)`[dim(UA_data)[1]]
# For now just fix at the median

# Serial interval
w_mean <- wtSIs$mean_si
w_sd <- wtSIs$std_si
#alpha <- w_mean^2 / w_sd^2
#beta <- w_mean / w_sd^2
w_sim <- discr_si(length(I_ts):1, w_mean, w_sd) # discrete distn of the serial interval


# Predictions
tt_pred <- c(I_ts, rep(NA, 10))
s <- sum(is.na(tt_pred))
for (i in 1:s){
  current_t <- length(tt_pred) - sum(is.na(tt_pred)) + 1
  tt_pred[current_t] <- Rt * sum(discr_si((current_t-1):1, w_mean, w_sd)*tt_pred[1:(current_t-1)])
}
tt_pred
plot(tt_pred, pch = 19) # seems reasonable

#### Function for generating forecasts ####
#' @param UA unitary authority or region
#' @param R0estimates output from EpiEstim, time series of R(t) by region, date
#' @param Rt_sim fix R(t) at its median value (FALSE), or sample from its posterior (TRUE)
#' @param Nsim number of simulations
#' @param Nt number of days to predict for
#' @param w_mean mean of the distribution for the serial interval
#' @param w_sd standard deviation of the distribution for the serial interval
#'
GenerateForecast <- function(UA, R0estimates, Rt_sim = FALSE, Nsim = 1000, Nt = 14, w_mean = 4.559007, w_sd = 4.530451){
  
  UA_data <- subset(R0timeseries, GSS_NM == UA) # subsets by chosen region
  I_ts <- UA_data$incidence # finds past time series of cases
  t <- length(I_ts) # current time t
  start_cases <- UA_data$cumulative_cases[1] # number of cases at startof time period
  
  if (Rt_sim == TRUE){
    # Draw Nsim samples from distribution for R(t)
    mean_r <- UA_data$`Mean(R)`[dim(UA_data)[1]]
    sd_r <- UA_data$`Std(R)`[dim(UA_data)[1]]
    alpha <- mean_r^2 / sd_r^2
    beta <- mean_r / sd_r^2
    Rt <- rgamma(Nsim, alpha, beta)
  }
  
  else {
    Rt <- rep(UA_data$`Median(R)`[dim(UA_data)[1]], Nsim) # fix at median
  }

  Preds <- matrix(0, Nsim, t + Nt)
  
  for (i in 1:Nsim){
    Preds[i,1:t] <- I_ts
    tmpPred <- rep(NA, Nt) # fill in as get new data so can append non-NA to past data
    for (j in 1:Nt){
      SIweight <- discr_si((t+j-1):1, w_mean, w_sd) # finds discrete distn of the serial interval
      PastData <- c(I_ts, tmpPred[which(!is.na(tmpPred))])
      PredMean <- Rt[i] * sum(SIweight*PastData) # mean of Poisson
      Preds[i,j+t] <- tmpPred[j] <- rpois(1, PredMean)
    }
    Preds[i,] <- cumsum(Preds[i,]) + start_cases
  }
  return(list(Preds = Preds, Rt = Rt))
}

#### Loop over all UAs ####
region_list <- RegionNames$GSS_NM
Nsim <- 1000
Nt <- 14
AllUAPreds <- array(0, dim = c(Nsim, 22+Nt, length(region_list)))
for (k in 1:length(region_list)){
  tmpForecast <- GenerateForecast('Devon', R0timeseries, Rt_sim = TRUE)
  AllUAPreds[,,k] <- tmpForecast$Preds
}

EnglandPreds <- apply(AllUAPreds, c(1,2), sum)
plot(1:dim(EnglandPreds)[2], (EnglandPreds[1,]), type = 'l', xlab = 'Day', ylab = 'Incidences', ylim = c(0,100000))
abline(v = t, col = 'red', lwd = 2)
for (i in 2:Nsim){lines(1:dim(EnglandPreds)[2], (EnglandPreds[i,]))}
qqs <- apply((EnglandPreds), 2, quantile, probs = c(0.025, 0.5, 0.975))
lines(1:dim(EnglandPreds)[2], qqs[2,], col = 'blue', lwd = 3, lty = 1)
lines(1:dim(EnglandPreds)[2], qqs[1,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(EnglandPreds)[2], qqs[3,], col = 'blue', lwd = 3, lty = 2)

ind <- which(region_list == 'Devon')
par(mfrow=c(1,2), mar=c(4,2,2,2))
plot(1:dim(AllUAPreds)[2], log(AllUAPreds[1,]), type = 'l', xlab = 'Day', ylab = 'Incidences', ylim = c(0,8))
abline(v = t, col = 'red', lwd = 2)
for (i in 2:Nsim){lines(1:dim(AllUAPreds)[2], log(AllUAPreds[i,]))}
qqs <- apply(log(AllUAPreds), 2, quantile, probs = c(0.025, 0.5, 0.975))
lines(1:dim(AllUAPreds)[2], qqs[2,], col = 'blue', lwd = 3, lty = 1)
lines(1:dim(AllUAPreds)[2], qqs[1,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(AllUAPreds)[2], qqs[3,], col = 'blue', lwd = 3, lty = 2)




FixedRt <- GenerateForecast('Devon', R0timeseries, Rt_sim = FALSE)
SampledRt <- GenerateForecast('Devon', R0timeseries, Rt_sim = TRUE)

png('Plots/Devon_EpiEstim.png', width = 800, height = 500)
par(mfrow=c(1,2), mar=c(4,2,2,2))
plot(1:dim(FixedRt$Preds)[2], log(FixedRt$Preds[1,]), type = 'l', xlab = 'Day', ylab = 'Incidences', ylim = c(0,8))
abline(v = t, col = 'red', lwd = 2)
for (i in 2:Nsim){lines(1:dim(FixedRt$Preds)[2], log(FixedRt$Preds[i,]))}
qqs <- apply(log(FixedRt$Preds), 2, quantile, probs = c(0.025, 0.5, 0.975))
lines(1:dim(FixedRt$Preds)[2], qqs[2,], col = 'blue', lwd = 3, lty = 1)
lines(1:dim(FixedRt$Preds)[2], qqs[1,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(FixedRt$Preds)[2], qqs[3,], col = 'blue', lwd = 3, lty = 2)

plot(1:dim(SampledRt$Preds)[2], log(SampledRt$Preds[1,]), type = 'l', xlab = 'Day', ylab = 'Incidences', ylim=c(0,8))
abline(v = t, col = 'red', lwd = 2)
for (i in 2:Nsim){lines(1:dim(SampledRt$Preds)[2], log(SampledRt$Preds[i,]))}
qqs2 <- apply(log(SampledRt$Preds), 2, quantile, probs = c(0.025, 0.5, 0.975))
lines(1:dim(SampledRt$Preds)[2], qqs2[2,], col = 'blue', lwd = 3, lty = 1)
lines(1:dim(SampledRt$Preds)[2], qqs2[1,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(SampledRt$Preds)[2], qqs2[3,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(SampledRt$Preds)[2], qqs[2,], col = 'green', lwd = 3, lty = 1)
dev.off()


par(mfrow=c(1,2), mar=c(4,4,2,2))
plot(1:dim(FixedRt$Preds)[2], (FixedRt$Preds[1,]), type = 'l', xlab = 'Day', ylab = 'Incidences', ylim = c(0,600), main = 'Devon, fixed R(t)')
abline(v = t, col = 'red', lwd = 2)
for (i in 2:Nsim){lines(1:dim(FixedRt$Preds)[2], (FixedRt$Preds[i,]))}
qqs <- apply((FixedRt$Preds), 2, quantile, probs = c(0.025, 0.5, 0.975))
lines(1:dim(FixedRt$Preds)[2], qqs[2,], col = 'blue', lwd = 3, lty = 1)
lines(1:dim(FixedRt$Preds)[2], qqs[1,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(FixedRt$Preds)[2], qqs[3,], col = 'blue', lwd = 3, lty = 2)

plot(1:dim(SampledRt$Preds)[2], (SampledRt$Preds[1,]), type = 'l', xlab = 'Day', ylab = 'Incidences', ylim=c(0,600), main = 'Devon, sampled R(t)')
abline(v = t, col = 'red', lwd = 2)
for (i in 2:Nsim){lines(1:dim(SampledRt$Preds)[2], (SampledRt$Preds[i,]))}
qqs2 <- apply((SampledRt$Preds), 2, quantile, probs = c(0.025, 0.5, 0.975))
lines(1:dim(SampledRt$Preds)[2], qqs2[2,], col = 'blue', lwd = 3, lty = 1)
lines(1:dim(SampledRt$Preds)[2], qqs2[1,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(SampledRt$Preds)[2], qqs2[3,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(SampledRt$Preds)[2], qqs[2,], col = 'green', lwd = 3, lty = 1)





plot(1:dim(NewPred$Preds)[2], NewPred$Preds[1,], type = 'l', xlab = 'Day', ylab = 'Incidences')
abline(v = t, col = 'red', lwd = 2)
for (i in 2:Nsim){lines(1:dim(NewPred$Preds)[2], NewPred$Preds[i,])}
qqs <- apply(NewPred$Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))
lines(1:dim(NewPred$Preds)[2], qqs[2,], col = 'blue', lwd = 3, lty = 1)
lines(1:dim(NewPred$Preds)[2], qqs[1,], col = 'blue', lwd = 3, lty = 2)
lines(1:dim(NewPred$Preds)[2], qqs[3,], col = 'blue', lwd = 3, lty = 2)



