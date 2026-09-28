# Clean R environment
rm(list = ls())

setwd("C:/Users/Utilisateur/Documents/R Script/Financial Econometrics/Final Assignement")

library(tidyverse)    # to organize data
library(tidyquant)    # to download from Yahoo Finance
library(forecast)     # for ACF Plot within ggplot (could also be done using acf())
library(ggpubr)       # for arrangement of the grids
library(rugarch)
library(xtable)

# download stock from Yahoo Finance 


# which ticker:
which.ticker <- "ORCL"

prices <- tq_get(which.ticker,
                 get = "stock.prices",
                 from = "2000-01-01",
                 to = "2026-08-09")

prices <- prices |> drop_na(adjusted)

returns <- prices |>
  arrange(date) |>
  mutate(ret = adjusted/lag(adjusted) - 1) |>
  select(symbol, date, ret)
returns

# remove missing value

returns <- returns |>
  drop_na(ret)


# compute log returns

log.returns <- prices |>
  arrange(date) |>
  mutate(logret = log(adjusted) - log(lag(adjusted))) |>
  select(symbol, date, logret)
log.returns

# remove missing value

log.returns <- log.returns |>
  drop_na(logret)

# collect all data DAILY for ORCL


date     <- prices$date[-1]
price    <- prices$adjusted[-1]
ret      <- returns$ret
grossret <- ret + 1
logret   <- log.returns$logret
sqret    <- logret^2
absret   <- abs(logret)

stock.df <- data.frame(date)
stock.df <- cbind(stock.df, price, ret, grossret, logret, sqret, absret)

# height and weight of plot that is saved
plot_width <- 6
plot_height <- 4

# plot the stock

a <- prices |>
  ggplot(aes(x = date, y = adjusted)) +
  geom_line() + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  )
a

# save pdf

pdf( paste0(c("Plots/",which.ticker,"_Prices_Daily.pdf"),collapse = ""), width=plot_width, height=plot_height)
print(a)
dev.off()


# plot return
b <- returns |>
  ggplot(aes(x = date, y = ret)) +
  geom_line() + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  )
b

# save pdf

pdf( paste0(c("Plots/",which.ticker,"_Returns_Daily.pdf"),collapse = ""), width=plot_width, height=plot_height)
print(b)
dev.off()

# plot log.return
c <- log.returns |>
  ggplot(aes(x = date, y = logret)) +
  geom_line() + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  )
c

# save pdf

pdf( paste0(c("Plots/",which.ticker,"_Log_Returns_Daily.pdf"),collapse = ""), width=plot_width, height=plot_height)
print(c)
dev.off()

# plot squared log.return
d <- stock.df |>
  ggplot(aes(x = date, y = sqret)) +
  geom_line() + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  )
d

# save pdf

pdf( paste0(c("Plots/",which.ticker,"_Squared_Log_Returns_Daily.pdf"),collapse = ""), width=plot_width, height=plot_height)
print(d)
dev.off()

# Plot ACF

aa <-ggAcf(ret, lag = 50, ci = 0.95) +
  geom_line() + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  )
aa

# save the pdf
pdf( paste0(c("Plots/",which.ticker,"_ACF_alt_Returns_Daily.pdf"),collapse = ""), width=plot_width, height=plot_height)
print(aa)
dev.off()

# squared returns ACF

bb <-ggAcf(sqret, lag = 50, ci = 0.95) +
  geom_line() + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  )
bb

# save the pdf
pdf( paste0(c("Plots/",which.ticker,"_ACF_alt_Squared_Returns_Daily.pdf"),collapse = ""), width=plot_width, height=plot_height)
print(bb)
dev.off()


# Estimate all models

# estimate the a GARCH(1,1) model with normal innovations using the rugarch 
# package

spec.uGARCH <- rugarch::ugarchspec(
  variance.model = list( model = "sGARCH",
                         garchOrder = c(1,1)),
  mean.model    = list( armaOrder = c(0,0),
                        include.mean = TRUE),
  distribution.model = "norm")

# estimate the model
fit.N.ORCL <-  ugarchfit(data = log.returns$logret, spec = spec.uGARCH)

# with t-distributed innovations

spec.uGARCH.t <- rugarch::ugarchspec(
  variance.model = list( model = "sGARCH",
                         garchOrder = c(1,1)),
  mean.model    = list( armaOrder = c(0,0),
                        include.mean = TRUE),
  distribution.model = "std")

# estimate the model
fit.t.ORCL <-  ugarchfit(data = log.returns$logret, spec = spec.uGARCH.t)

# EGARCH - N

spec.uGARCH.EGARCH <- rugarch::ugarchspec(
  variance.model = list( model = "eGARCH",
                         garchOrder = c(1,1)),
  mean.model    = list( armaOrder = c(0,0),
                        include.mean = TRUE),
  distribution.model = "norm")

# estimate the model
fit.EGARCH.ORCL <-  ugarchfit(data = log.returns$logret, spec = spec.uGARCH.EGARCH)

# GJR-GARCH - N

# with t-distributed innovations

spec.uGARCH.GJR <- rugarch::ugarchspec(
  variance.model = list( model = "gjrGARCH",
                         garchOrder = c(1,1)),
  mean.model    = list( armaOrder = c(0,0),
                        include.mean = TRUE),
  distribution.model = "norm")

# estimate the model
fit.GJR.ORCL <-  ugarchfit(data = log.returns$logret, spec = spec.uGARCH.GJR)


# build a tables for ORCL


tab.ORCL <- matrix(NA,nrow = 4, ncol = 8)

# GARCH-N
tab.ORCL[1,1:4]  <- coef(fit.N.ORCL)

# GARCH-t
tab.ORCL[2,1:5]  <- coef(fit.t.ORCL)

# EGARCH
tmp <- coef(fit.EGARCH.ORCL)
tab.ORCL[3,1:4]  <- tmp[1:4]
tab.ORCL[3,8]    <- tmp[5]

# GJR-GARCH
tmp <- coef(fit.GJR.ORCL)
tab.ORCL[4,1:4]  <- tmp[1:4]
tab.ORCL[4,8]    <- tmp[5]

colnames(tab.ORCL) <- c("$\\mu$", "$\\omega$", "$\\alpha$", "$\\beta$", 
                       "$\\nu$", "$\\phi$", "$\\theta$",
                       "$\\gamma$") #, "$\\delta$")

rownames(tab.ORCL) <- c( "GARCH-N", "GARCH-t",
                        "EGARCH", "GJR-GARCH" )

print(xtable(caption = "Parameter estimates of GARCH models", tab.ORCL,  
             label = "ORCL", digits=c(0,5,4,3,3,1,2,2,2)), 
      include.rownames = TRUE, include.colnames = TRUE, 
      sanitize.text.function = I, hline.after = c(-1,-1,0,4,4), 
      type="latex", file = "Data/ORCL.tex")


# organize the results of sigma.t

T.obs <- length(log.returns$logret)

sig.ORCL <- matrix(NA, nrow = T.obs, ncol = 4)

colnames(sig.ORCL) <- c( "GARCH_N", "GARCH_t",
                        "EGARCH", "GJR_GARCH")


sig.ORCL[,1] <- sigma(fit.N.ORCL)
sig.ORCL[,2] <- sigma(fit.t.ORCL)
sig.ORCL[,3] <- sigma(fit.EGARCH.ORCL)
sig.ORCL[,4] <- sigma(fit.GJR.ORCL)


df.sig.ORCL <- as.data.frame(sig.ORCL)
df.sig.ORCL <- cbind(log.returns$date,df.sig.ORCL)

k <- autoplot(ts(sig.ORCL)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  )

pdf( paste0(c("Plots/Sigma_ORCL.pdf"),collapse = ""), width=12, height=plot_height)
print(k)
dev.off()

Estimate.VaR.ES.models.IS <- function(ret, alpha, date) {
  
  # Prepare data
  
  ret <- 100 * ret
  ret <- as.vector(unlist(ret))
  
  n_days <- length(ret)
  
  if (length(date) != n_days) {
    stop("date and ret must have the same length.")
  }
  
  VaR.ES.df <- data.frame(
    date = date,
    ret  = ret
  )
  
  pars <- list()
  
  
  # Help :
  # extract in-sample forecasts from Optimize_models()
  # added an in_sample = TRUE argument
  # check modified_recursive_models.R for more details

  
  get_IS_forecasts <- function(x) {
    
    if (!is.null(x$ISforecasts)) {
      return(x$ISforecasts)
    }
    
    if (!is.null(x$OOSforecasts)) {
      return(x$OOSforecasts)
    }
    
    stop("Optimize_models() did not return forecasts.")
  }
  
  # 1. Historical Simulation - 250 day window
  
  HistSim_QES <- function(x, alpha) {
    
    x <- x[is.finite(x)]
    
    x_Q <- as.numeric(
      quantile(x, probs = alpha, type = 7, na.rm = TRUE)
    )
    
    x_ES <- mean(
      x[x <= x_Q],
      na.rm = TRUE
    )
    
    return(
      list(
        Q = x_Q,
        ES = x_ES
      )
    )
  }
  
  n_HSwindow <- 250
  
  VaR.ES.df$q_HistSim <- NA_real_
  VaR.ES.df$e_HistSim <- NA_real_
  
  if (n_days > n_HSwindow) {
    
    for (t in (n_HSwindow + 1):n_days) {
      
      HistSim_QES_help <- HistSim_QES(
        ret[(t - n_HSwindow):(t - 1)],
        alpha
      )
      
      VaR.ES.df$q_HistSim[t] <- HistSim_QES_help$Q
      VaR.ES.df$e_HistSim[t] <- HistSim_QES_help$ES
    }
  }
  
  VaR.ES.df$hit_HistSim <- ret < VaR.ES.df$q_HistSim
  
  # 2. GARCH(1,1)-N
  
  GARCH_N.spec <- ugarchspec(
    mean.model = list(
      include.mean = FALSE,
      armaOrder = c(0, 0)
    ),
    variance.model = list(
      model = "sGARCH",
      garchOrder = c(1, 1)
    ),
    distribution.model = "norm"
  )
  
  GARCH_N.fit <- ugarchfit(
    spec = GARCH_N.spec,
    data = ret
  )
  
  sigma_GARCH_N <- as.numeric(
    rugarch::sigma(GARCH_N.fit)
  )
  
  VaR.ES.df$q_GARCH_N <-
    qnorm(alpha) * sigma_GARCH_N
  
  VaR.ES.df$e_GARCH_N <-
    -dnorm(qnorm(alpha)) / alpha *
    sigma_GARCH_N
  
  VaR.ES.df$hit_GARCH_N <- ret < VaR.ES.df$q_GARCH_N
  
  pars$GARCH_N <- coef(GARCH_N.fit)
  
  # 3. GJR-GARCH(1,1)-t
  
  GJRGARCH_t.spec <- ugarchspec(
    mean.model = list(
      include.mean = FALSE,
      armaOrder = c(0, 0)
    ),
    variance.model = list(
      model = "gjrGARCH",
      garchOrder = c(1, 1)
    ),
    distribution.model = "std"
  )
  
  GJRGARCH_t.fit <- ugarchfit(
    spec = GJRGARCH_t.spec,
    data = ret
  )
  
  sigma_GJRGARCH_t <- as.numeric(
    rugarch::sigma(GJRGARCH_t.fit)
  )
  
  # Estimated degrees of freedom
  nu <- as.numeric(
    coef(GJRGARCH_t.fit)["shape"]
  )
  
  # Student-t quantile standardized to variance = 1
  q_t <- rugarch::qdist(
    distribution = "std",
    p = alpha,
    mu = 0,
    sigma = 1,
    shape = nu
  )
  
  # VaR
  VaR.ES.df$q_GJRGARCH_t <- q_t * sigma_GJRGARCH_t
  
  # Expected Shortfall of standardized Student-t
  ES_t <- -(
    dt(q_t, df = nu) *
      (nu + q_t^2) /
      ((nu - 1) * alpha)
  ) *
    sqrt((nu - 2) / nu)
  
  VaR.ES.df$e_GJRGARCH_t <- ES_t * sigma_GJRGARCH_t
  
  VaR.ES.df$hit_GJRGARCH_t <- ret < VaR.ES.df$q_GJRGARCH_t
  
  pars$GJRGARCH_t <- coef(GJRGARCH_t.fit)
  
  # 4. AS-CAViaR / ES-CAViaR
  
  CAREAS.IS <- Optimize_models(
    r = ret,
    alpha = alpha,
    model = "CAREAS",
    in_sample = TRUE
  )
  
  CAREAS.forecasts <- get_IS_forecasts(CAREAS.IS)
  
  if (nrow(CAREAS.forecasts) != n_days) {
    stop(
      paste0(
        "CAREAS returned ", nrow(CAREAS.forecasts),
        " forecasts, but there are ", n_days, " returns."
      )
    )
  }
  
  VaR.ES.df$q_CAREAS <- CAREAS.forecasts[, 1]
  
  VaR.ES.df$e_CAREAS <- CAREAS.forecasts[, 2]
  
  VaR.ES.df$hit_CAREAS <- ret < VaR.ES.df$q_CAREAS
  
  pars$CAREAS <- CAREAS.IS$parameters

  # 5. SAV-CAViaR / ES-CAViaR
  
  CARESAV.IS <- Optimize_models(
    r = ret,
    alpha = alpha,
    model = "CARESAV",
    in_sample = TRUE
  )
  
  CARESAV.forecasts <- get_IS_forecasts(CARESAV.IS)
  
  if (nrow(CARESAV.forecasts) != n_days) {
    stop(
      paste0(
        "CARESAV returned ", nrow(CARESAV.forecasts),
        " forecasts, but there are ", n_days, " returns."
      )
    )
  }
  
  VaR.ES.df$q_CARESAV <- CARESAV.forecasts[, 1]
  
  VaR.ES.df$e_CARESAV <- CARESAV.forecasts[, 2]
  
  VaR.ES.df$hit_CARESAV <- ret < VaR.ES.df$q_CARESAV
  
  pars$CARESAV <- CARESAV.IS$parameters
  
  # 6. One-factor GAS model
  
  GAS1F.IS <- Optimize_models(
    r = ret,
    alpha = alpha,
    model = "GAS1F",
    in_sample = TRUE
  )
  
  GAS1F.forecasts <- get_IS_forecasts(GAS1F.IS)
  
  if (nrow(GAS1F.forecasts) != n_days) {
    stop(
      paste0("GAS1F returned ",nrow(GAS1F.forecasts),
        " forecasts, but there are ",n_days," returns."
      )
    )
  }
  
  VaR.ES.df$q_GAS1F <- GAS1F.forecasts[, 1]
  
  VaR.ES.df$e_GAS1F <- GAS1F.forecasts[, 2]
  
  VaR.ES.df$hit_GAS1F <- ret < VaR.ES.df$q_GAS1F
  
  pars$GAS1F <- GAS1F.IS$parameters
  
  # Return results
  
  return(list(VaR.ES.df = VaR.ES.df,pars = pars))
}
source("modified_recursive_models.R")

results <- Estimate.VaR.ES.models.IS(ret = logret, alpha = 0.025, date = date)

# plot the different VaR models 

#[n.start:n.days]

results$VaR.ES.df.Q <- results$VaR.ES.df |>
  ggplot(aes(x = date, y = ret)) +
  geom_line()+ 
  geom_line(aes(y=q_GARCH_N, colour = "GARCH-N"),size=0.1) + 
  geom_line(aes(y=q_GJRGARCH_t, colour = "GJR-t"),size=0.1) +
  geom_line(aes(y=q_GAS1F, colour = "GAS1F"),size=0.1)+
  geom_line(aes(y=q_HistSim, colour = "HistSim"),size=0.1)+  
  geom_line(aes(y=q_CAREAS, colour = "AS-CAViaR"),size=0.1)+ 
  geom_line(aes(y=q_CARESAV, colour = "SAV-CAViaR"),size=0.1)+    
  theme_bw() +
  theme(axis.line = element_line(colour = c("black")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_color_manual(name = "VaR models", 
                     values = c("GARCH-N" = "darkred", 
                                "GJR-t"   = "orange",
                                "GAS1F"   = "darkorange",
                                "HistSim"     = "brown",
                                "AS-CAViaR"   = "darkblue",
                                "SAV-CAViaR"  = "blue",
                                "return"      = "black"))
results$VaR.ES.df.Q
pdf( paste0(c("Plots/ORCL_Q.pdf"),collapse = ""), width=12, height=plot_height)
print(results$VaR.ES.df.Q)
dev.off()

# plot the different ES models 

#[n.start:n.days]

results$VaR.ES.df.ES <- results$VaR.ES.df |>
  ggplot(aes(x = date, y = ret)) +
  geom_line()+ 
  geom_line(aes(y=e_GARCH_N, colour = "GARCH-N"),size=0.1) + 
  geom_line(aes(y=e_GJRGARCH_t, colour = "GJR-t"),size=0.1) +
  geom_line(aes(y=e_GAS1F, colour = "GAS1F"),size=0.1)+
  geom_line(aes(y=e_HistSim, colour = "HistSim"),size=0.1)+  
  geom_line(aes(y=e_CAREAS, colour = "AS-CAViaR"),size=0.1)+ 
  geom_line(aes(y=e_CARESAV, colour = "SAV-CAViaR"),size=0.1)+    
  theme_bw() +
  theme(axis.line = element_line(colour = c("black")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_color_manual(name = "ES models", 
                     values = c("GARCH-N" = "darkred", 
                                "GJR-t"   = "orange",
                                "GAS1F"   = "darkorange",
                                "HistSim"     = "brown",
                                "AS-CAViaR"   = "darkblue",
                                "SAV-CAViaR"  = "blue",
                                "return"      = "black"))
results$VaR.ES.df.ES
pdf( paste0(c("Plots/ORCL_ES.pdf"),collapse = ""), width=12, height=plot_height)
print(results$VaR.ES.df.ES)
dev.off()

results$VaR.ES.df.hits <- cbind(results$VaR.ES.df$hit_GARCH_N, results$VaR.ES.df$hit_GJRGARCH_t,
                         results$VaR.ES.df$hit_HistSim, results$VaR.ES.df$hit_CARESAV,
                         results$VaR.ES.df$hit_CAREAS,results$VaR.ES.df$hit_GAS1F)

hits.vola <- rbind(colSums(results$VaR.ES.df.hits[,1:6], na.rm = TRUE)/6688)

rownames(hits.vola) <- c("ORCL")
colnames(hits.vola) <- c("GARCH-N", "GJR-GARCH-t", "HistSim", 
                         "SAV-CAViaR","AS-CAViaR","GAS1F")

print(xtable(caption = "Hit ratio for VaR 2.5%", hits.vola,  
             label = "hits_vola", digits=6), 
      include.rownames = TRUE, include.colnames = TRUE, 
      sanitize.text.function = I, hline.after = c(-1,0,1), 
      type="latex", file = "Data/hits_vola.tex")


## FORECASTING EXERCISE ##

# use obs until 31/12/2019 as in-sample estimation period

T.IS_alt  <- 5029
T.OOS_alt <- T.obs - T.IS

#returns.IS <- returns[1:T.IS]
#returns.OOS <- returns[(T.IS+1):(T.IS+T.OOS)]

# Here we only want to forecast on the last 200

T.IS <- 6488
T.OOS <- T.obs - T.IS

# Fixed-window function

forecast_fixed <- function(spec, returns, N_IS, N_OOS) {
  
  forecasts <- numeric(N_OOS)
  
  # Estimate model only once using information only up to t
  fit <- ugarchfit(spec = spec, data = returns[1:N_IS],solver = "hybrid")
 
  # Fixed parameter estimates
  coef_fixed <- coef(fit)
  
  # Rebuild specification using the FIXED coefficients
  spec_fixed <- getspec(fit)
  setfixed(spec_fixed) <- as.list(coef_fixed)
  
  for (i in 1:N_OOS) {
    # End of estimation sample
    t <- N_IS + i - 1
    
    # Fixed-length estimation window
    estimation_data <- returns[1:t]
    
    # One-step-ahead forecast for t+1
    fc <- ugarchforecast(spec_fixed, data = estimation_data, n.ahead = 1)
    
    forecasts[i] <- as.numeric(sigma(fc)^2)
  }
  
  forecasts
}

# run fixed-window forecasting for our three models
fixed_garch_n <- forecast_fixed(spec.uGARCH, returns$ret, T.IS, T.OOS)
fixed_garch_t <- forecast_fixed(spec.uGARCH.t, returns$ret, T.IS, T.OOS)
fixed_gjr_n <- forecast_fixed(spec.uGARCH.GJR, returns$ret, T.IS, T.OOS)

# generate a data frame that contains the actual data and the forecasts

date.OOS <- c(1:T.OOS)

df.simfix <- data.frame(cbind(date.OOS,returns$ret[(T.IS+1):(T.IS+T.OOS)],
                              fixed_garch_n,fixed_garch_t,fixed_gjr_n))

# plot the realized and forecasted stock returns:
p <- df.simfix |>
  ggplot(aes(x = date.OOS, y = returns$ret[(T.IS+1):(T.IS+T.OOS)]^2)) +
  geom_line(alpha = 0.4)+ 
  geom_line(aes(y = returns$ret[(T.IS+1):(T.IS+T.OOS)]^2, colour = "Realized Returns (Squared)"),
            linewidth = 0.6, alpha = 0.4) +
  geom_line(aes(y=fixed_garch_n, colour = "GARCH-N"),size=0.1) + 
  geom_line(aes(y=fixed_garch_t, colour = "GARCH-t"),size=0.1) +   
  geom_line(aes(y=fixed_gjr_n, colour = "GJR-GARCH-N"),size=0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = c("black")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_color_manual(name = "Realized and Forecasted Series\n(Fixed Estimation Window)", 
                     values = c("Realized Returns (Squared)"="black",
                                "GARCH-N" = "green", 
                                "GARCH-t" = "blue",
                                "GJR-GARCH-N" = "red"
                                ))
p

pdf( paste0(c("Plots/Simfix.pdf"),collapse = ""), width=10, height=plot_height)
print(p)
dev.off()

# Expanding-window function

forecasts_exp <- function(spec, returns, N_IS, N_OOS) {
  
  forecasts <- numeric(N_OOS)
  
  for (i in 1:N_OOS) {
    
    t <- N_IS + i - 1
    
    # Expanding information set: 1,...,t
    data_t <- returns[1:t]
    
    # Re-estimate model
    fit <- ugarchfit(spec = spec, data = data_t,solver = "hybrid")
    
    # One-step-ahead forecast
    fc <- ugarchforecast(fit,n.ahead = 1)
    
    forecasts[i] <- as.numeric(sigma(fc)^2)
  }
  
  forecasts
}

# run fixed-window forecasting for our three models
exp_garch_n <- forecasts_exp(spec.uGARCH, returns$ret, T.IS, T.OOS)
exp_garch_t <- forecasts_exp(spec.uGARCH.t, returns$ret, T.IS, T.OOS)
exp_gjr_n <- forecasts_exp(spec.uGARCH.GJR, returns$ret, T.IS, T.OOS)

# generate a data frame that contains the actual data and the forecasts

df.simexp <- data.frame(cbind(date.OOS,returns$ret[(T.IS+1):(T.IS+T.OOS)],
                              exp_garch_n,exp_garch_t,exp_gjr_n))

# plot the realized and forecasted stock returns:
pp <- df.simexp |>
  ggplot(aes(x = date.OOS, y = returns$ret[(T.IS+1):(T.IS+T.OOS)]^2)) +
  geom_line(alpha = 0.4)+ 
  geom_line(aes(y = returns$ret[(T.IS+1):(T.IS+T.OOS)]^2, colour = "Realized Returns (Squared)"),
            linewidth = 0.6, alpha =0.4) +
  geom_line(aes(y=exp_garch_n, colour = "GARCH-N"),size=0.1) + 
  geom_line(aes(y=exp_garch_t, colour = "GARCH-t"),size=0.1) +   
  geom_line(aes(y=exp_gjr_n, colour = "GJR-GARCH-N"),size=0.1) +
  theme_bw() +
  theme(axis.line = element_line(colour = c("black")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_color_manual(name = "Realized and Forecasted Series\n(Expanding Estimation Window)", 
                     values = c("Realized Returns (Squared)"="black",
                                "GARCH-N" = "green", 
                                "GARCH-t" = "blue",
                                "GJR-GARCH-N" = "red"
                     ))
pp

pdf( paste0(c("Plots/Simexp.pdf"),collapse = ""), width=10, height=plot_height)
print(pp)
dev.off()

# FORECAST EVALUATION:

# Compute OOS statistics: Loss functions

# 1. MSE

mse.garch_n_fixed    = mean((returns$ret[(T.IS+1):(T.IS+T.OOS)]^2 - fixed_garch_n)^2)
mse.garch_t_fixed    = mean((returns$ret[(T.IS+1):(T.IS+T.OOS)]^2 - fixed_garch_t)^2)
mse.gjr_n_fixed      = mean((returns$ret[(T.IS+1):(T.IS+T.OOS)]^2 - fixed_gjr_n)^2)

mse.garch_n_exp    = mean((returns$ret[(T.IS+1):(T.IS+T.OOS)]^2 - exp_garch_n)^2)
mse.garch_t_exp    = mean((returns$ret[(T.IS+1):(T.IS+T.OOS)]^2 - exp_garch_t)^2)
mse.gjr_n_exp      = mean((returns$ret[(T.IS+1):(T.IS+T.OOS)]^2 - exp_gjr_n)^2)

df.mse <- data.frame(cbind(mse.garch_n_fixed,mse.garch_n_exp, mse.garch_t_fixed,
                           mse.garch_t_exp, mse.gjr_n_fixed,mse.gjr_n_exp))

# 2 QLIKE (from Patton [2011] Appendix which is the QLIKE loss function up to additive and multiplicative constants)

qlike <- function(realized,forecast){
  loss <- realized/forecast - log(realized/forecast) - 1
  loss <- mean(loss)
}

qlike.garch_n_fixed    = qlike(returns$ret[(T.IS+1):(T.IS+T.OOS)]^2, fixed_garch_n)
qlike.garch_t_fixed    = qlike(returns$ret[(T.IS+1):(T.IS+T.OOS)]^2, fixed_garch_t)
qlike.gjr_n_fixed      = qlike(returns$ret[(T.IS+1):(T.IS+T.OOS)]^2, fixed_gjr_n)

qlike.garch_n_exp    = qlike(returns$ret[(T.IS+1):(T.IS+T.OOS)]^2, exp_garch_n)
qlike.garch_t_exp    = qlike(returns$ret[(T.IS+1):(T.IS+T.OOS)]^2, exp_garch_t)
qlike.gjr_n_exp      = qlike(returns$ret[(T.IS+1):(T.IS+T.OOS)]^2, exp_gjr_n)

df.qlike <- data.frame(cbind(qlike.garch_n_fixed,qlike.garch_n_exp, qlike.garch_t_fixed,
                           qlike.garch_t_exp, qlike.gjr_n_fixed,qlike.gjr_n_exp))

# MAYBE PLOT EACH FORECAST AGAINST ITS COUNTERPART ON SEPARATE PLOTS AND ADD IN APPENDIX FOR CLARITY

# Diebold-Mariano tests DM:

# 1.F MSE DM - Fixed Window

# These are pair-wise test: 

dm.mat_fix <- matrix(NA,nrow = 3,ncol = 3)


dm12_fix <- dm.test( fixed_garch_n, fixed_garch_t, alternative = "two.sided", h = 1,
                 power = 2, varestimator = "acf" )

dm.mat_fix[1,2] <- dm12_fix$p.value


dm13_fix <- dm.test( fixed_garch_n, fixed_gjr_n, alternative = "two.sided", h = 1,
                 power = 2, varestimator = "acf" )

dm.mat_fix[1,3] <- dm13_fix$p.value


dm23_fix <- dm.test( fixed_garch_t, fixed_gjr_n, alternative = "two.sided", h = 1,
                 power = 2, varestimator = "acf" )

dm.mat_fix[2,3] <- dm23_fix$p.value

print("DM tests for equal predictive ability with MSE - p-values (fixed-window)")
print(dm.mat_fix)

# 1.E MSE DM - Expanding Window

dm.mat_exp <- matrix(NA,nrow = 3,ncol = 3)


dm12_exp <- dm.test( exp_garch_n, exp_garch_t, alternative = "two.sided", h = 1,
                     power = 2, varestimator = "acf" )

dm.mat_exp[1,2] <- dm12_exp$p.value


dm13_exp <- dm.test( exp_garch_n, exp_gjr_n, alternative = "two.sided", h = 1,
                     power = 2, varestimator = "acf" )

dm.mat_exp[1,3] <- dm13_exp$p.value


dm23_exp <- dm.test( exp_garch_t, exp_gjr_n, alternative = "two.sided", h = 1,
                     power = 2, varestimator = "acf" )

dm.mat_exp[2,3] <- dm23_exp$p.value

print("DM tests for equal predictive ability with MSE - p-values (expanding-window)")
print(dm.mat_exp)

# 2. DM test with QLIKE

# Although a package exist (Forecomp), a stable solution can be coded from scratch

dm_qlike <-function(L1,L2){
  
  # Compute loss differential
  d <- L1 - L2
  
  # DM statistic
  dbar <- mean(d)
  T <- length(d)
  
  # For one-step-ahead forecasts:
  gamma0 <- var(d)
  
  DM <- dbar / sqrt(gamma0 / T)
  
  # Two-sided p-value
  p_value <- 2 * (1 - pnorm(abs(DM)))
  
  c(p_value = p_value)
}

# 2.F QLIKE DM - fixed window

dm.mat_fix_qlike <- matrix(NA,nrow = 3,ncol = 3)


dm12_fix_qlike <- dm_qlike( fixed_garch_n, fixed_garch_t)

dm.mat_fix_qlike[1,2] <- dm12_fix_qlike


dm13_fix_qlike <- dm_qlike( fixed_garch_n, fixed_gjr_n)

dm.mat_fix_qlike[1,3] <- dm13_fix_qlike


dm23_fix_qlike <- dm_qlike( fixed_garch_t, fixed_gjr_n)

dm.mat_fix_qlike[2,3] <- dm23_fix_qlike

print("DM tests for equal predictive ability with QLIKE - p-values (fixed-window)")
print(dm.mat_fix_qlike)

# 2.E QLIKE DM - expanding window

dm.mat_exp_qlike <- matrix(NA,nrow = 3,ncol = 3)


dm12_exp_qlike <- dm_qlike( exp_garch_n, exp_garch_t)

dm.mat_exp_qlike[1,2] <- dm12_exp_qlike


dm13_exp_qlike <- dm_qlike( exp_garch_n, exp_gjr_n)

dm.mat_exp_qlike[1,3] <- dm13_exp_qlike


dm23_exp_qlike <- dm_qlike( exp_garch_t, exp_gjr_n)

dm.mat_exp_qlike[2,3] <- dm23_exp_qlike

print("DM tests for equal predictive ability with QLIKE - p-values (expanding-window)")
print(dm.mat_exp_qlike)

# 10 days with highest volatility forecasts
# These are the same across models
print("10 days with highest volatility forecasts index")
order(exp_garch_t, decreasing=TRUE)[1:10]

# plot the highest volatility period
a_june <- prices[6642:6663,] |>
  ggplot(aes(x = date, y = adjusted)) +
  geom_line() + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  )

pdf( paste0(c("Plots/10_day_high_vol.pdf"),collapse = ""), width=10, height=plot_height)
print(a_june)
dev.off()

# 21-step ahead volatility forecast

# Fixed-window function for 21 steps

forecast_fixed_21 <- function(spec, returns, N_IS, N_OOS, h) {
  
  forecasts <- numeric(N_OOS-h+1)
  
  # Estimate model only once using information only up to t
  fit <- ugarchfit(spec = spec, data = returns[1:N_IS],solver = "hybrid")
  
  # Fixed parameter estimates
  coef_fixed <- coef(fit)
  
  # Rebuild specification using the FIXED coefficients
  spec_fixed <- getspec(fit)
  setfixed(spec_fixed) <- as.list(coef_fixed)
  
  for (i in 1:N_OOS-h+1) {
    # End of estimation sample
    t <- N_IS + i - 1
    
    # Fixed-length estimation window
    estimation_data <- returns[1:t]
    
    # One-step-ahead forecast for t+1
    fc <- ugarchforecast(spec_fixed, data = estimation_data, n.ahead = h)
    
    forecasts[i] <- as.numeric(sigma(fc)[h]^2)
  }
  
  forecasts
}

# Expanding-window function for 21 steps

forecasts_exp_21 <- function(spec, returns, N_IS, N_OOS, h) {
  
  forecasts <- numeric(N_OOS-h+1)
  
  for (i in 1:N_OOS-h+1) {
    
    t <- N_IS + i - 1
    
    # Expanding information set: 1,...,t
    data_t <- returns[1:t]
    
    # Re-estimate model
    fit <- ugarchfit(spec = spec, data = data_t,solver = "hybrid")
    
    # One-step-ahead forecast
    fc <- ugarchforecast(fit,n.ahead = h)
    
    forecasts[i] <- as.numeric(sigma(fc)[h]^2)
  }
  
  forecasts
}

# run forecasting for our 2 estimation schemes
fixed_garch_n_21 <- forecast_fixed_21(spec.uGARCH, returns$ret, T.IS, T.OOS, 21) 
exp_garch_n_21 <- forecasts_exp_21(spec.uGARCH, returns$ret, T.IS, T.OOS, 21)

# generate a data frame that contains the actual data and the forecasts
date.OOS_21 <- c(21:T.OOS)

df.sim_21 <- data.frame(
  date = date.OOS_21,
  realized_sq = returns$ret[(T.IS+21):(T.IS+T.OOS)]^2,
  fixed_garch_n = fixed_garch_n[21:200],
  fixed_garch_n_21 = fixed_garch_n_21,
  exp_garch_n = exp_garch_n[21:200],
  exp_garch_n_21 = exp_garch_n_21
)

# plot the realized and the 2 forecasted stock returns:
p_21 <- df.sim_21 |>
  ggplot(aes(x = date, y = realized_sq)) +
  geom_line(alpha = 0.4)+ 
  geom_line(aes(y = realized_sq, colour = "Realized Returns (Squared)"),
            linewidth = 0.6, alpha =0.4) +
  geom_line(aes(y=fixed_garch_n, colour = "GARCH-N"),size=0.1) + 
  geom_line(aes(y=fixed_garch_n_21, colour = "GARCH-N (21-step)"),size=0.1) + 
  theme_bw() +
  theme(axis.line = element_line(colour = c("black")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_color_manual(name = "Iterated Single and Multi-step Forecasts\n(Fixed Estimation Window)", 
                     values = c("Realized Returns (Squared)"="black",
                                "GARCH-N" = "green", 
                                "GARCH-N (21-step)" = "blue"
                     ))
p_21
pdf( paste0(c("Plots/Simfix_21.pdf"),collapse = ""), width=10, height=plot_height)
print(p_21)
dev.off()

pp_21 <- df.sim_21 |>
  ggplot(aes(x = date, y = realized_sq)) +
  geom_line(alpha = 0.4)+ 
  geom_line(aes(y = realized_sq, colour = "Realized Returns (Squared)"),
            linewidth = 0.6, alpha =0.4) +
  geom_line(aes(y=exp_garch_n, colour = "GARCH-N"),size=0.1) + 
  geom_line(aes(y=exp_garch_n_21, colour = "GARCH-N (21-step)"),size=0.1) + 
  theme_bw() +
  theme(axis.line = element_line(colour = c("black")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(
    x = NULL,
    y = NULL
  ) +
  scale_color_manual(name = "Iterated Single and Multi-step Forecasts\n(Expanding Estimation Window)", 
                     values = c("Realized Returns (Squared)"="black",
                                "GARCH-N" = "green", 
                                "GARCH-N (21-step)" = "blue"
                     ))
pp_21
pdf( paste0(c("Plots/Simexp_21.pdf"),collapse = ""), width=10, height=plot_height)
print(pp_21)
dev.off()

# Compute OOS statistics: Loss functions

# For 21-step

mse.garch_n_fixed_21    = mean((df.sim_21$realized_sq - fixed_garch_n_21)^2)
mse.garch_n_exp_21    = mean((df.sim_21$realized_sq - exp_garch_n_21)^2)

qlike.garch_n_fixed_21    = qlike(df.sim_21$realized_sq, fixed_garch_n_21)
qlike.garch_n_exp_21    = qlike(df.sim_21$realized_sq, exp_garch_n_21)

# For One-step but Realigned

mse.garch_n_fixed_alt    = mean((df.sim_21$realized_sq - df.sim_21$fixed_garch_n)^2)
mse.garch_n_exp_alt    = mean((df.sim_21$realized_sq - df.sim_21$exp_garch_n)^2)

qlike.garch_n_fixed_alt    = qlike(df.sim_21$realized_sq, df.sim_21$fixed_garch_n)
qlike.garch_n_exp_alt    = qlike(df.sim_21$realized_sq, df.sim_21$exp_garch_n)

df.mse_21 <- data.frame(cbind(mse.garch_n_fixed_21,mse.garch_n_fixed_alt,
                           mse.garch_n_exp_21, mse.garch_n_exp_alt))

df.qlike_21 <- data.frame(cbind(qlike.garch_n_fixed_21,qlike.garch_n_fixed_alt,
                           qlike.garch_n_exp_21, qlike.garch_n_exp_alt))

# Diebold-Mariano tests DM:

# MSE DM - 21-day-ahead

dm_fix_21 <- dm.test( df.sim_21$fixed_garch_n, fixed_garch_n_21, alternative = "two.sided", h = 1,
                     power = 2, varestimator = "acf" )
dm_exp_21 <- dm.test( df.sim_21$exp_garch_n, exp_garch_n_21, alternative = "two.sided", h = 1,
                      power = 2, varestimator = "acf" )

dm.mat_fix_21 <- dm_fix_21$p.value
dm.mat_exp_21 <- dm_exp_21$p.value

# QLIKE DM - 21-day-ahead

dm_fix_qlike_21 <- dm_qlike( fixed_garch_n_21, df.sim_21$fixed_garch_n)
dm_exp_qlike_21 <- dm_qlike( exp_garch_n_21, df.sim_21$exp_garch_n)