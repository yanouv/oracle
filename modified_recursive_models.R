#' @title One-factor GAS model
#' @description This function calculates the average FZ loss of the one factor GAS model similar to PSC (2019)'s MATLAB function GAS_onefactor_LL3.m
#' @import stats
#'
#' @param theta Vector of parameters. theta = [beta,gamma,b,c]. Set beta = pnorm(beta), gamma = exp(gamma), b=-exp(b), c = pnorm(c). omega = 0 for identification
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return one factor GAS model quantile and ES forecasts
#' @export
GAS_onefactor <- function( theta, r, alpha, tau = -1 ){

  # read out the parameters:
  beta  = pnorm(theta[1])
  gamma = exp(theta[2])
  b     = - exp(theta[3])
  c     = pnorm(theta[4])
  a     = c*b

  omega = 0
  kalpha = 1         # normalizing positive constant for Hessian. If numerical optimization fails, use a value smaller than unity

  T <- length(r)

  kappa_t   = array(NaN,dim=c(T,1))
  forecasts = array(NaN,dim=c(T,2))       # quantile and ES

  loss      = array(0,dim=c(T,1))

  # initial values:
  kappa_t[1,1]  = omega/(1-beta)     # unconditional average of kappa
  forecasts[1,] = cbind(a * exp(kappa_t[1,1]), b * exp(kappa_t[1,1]))

  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1]<= forecasts[1,1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-forecasts[1,1])))
    # smooth approximation to the indicator function
  }

  loss[1,1] <- -1/alpha/forecasts[1,2]*hits*(forecasts[1,1]-r[1]) - 1/forecasts[1,2]*(forecasts[1,2]-forecasts[1,1]) + log(-forecasts[1,2])

  for (tt in 2:T){

    hitsL <- hits   # lagged hit variable

    # forcing variable
    kappa_t[tt,1] <- omega + beta * kappa_t[tt-1,1] + gamma/kalpha/forecasts[tt-1,2]*(1/alpha* hitsL *r[tt-1] - forecasts[tt-1,2])

    # forecasts
    forecasts[tt,] = cbind(a * exp(kappa_t[tt,1]), b * exp(kappa_t[tt,1]))

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= forecasts[tt,1])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-forecasts[tt,1])))
      # smooth approximation to the indicator function
    }

    loss[tt,1] <- -1/alpha/forecasts[tt,2]*hits*(forecasts[tt,1]-r[tt]) - 1/forecasts[tt,2]*(forecasts[tt,2]-forecasts[tt,1]) + log(-forecasts[tt,2])

  }

  # average loss
  Eloss = mean(loss)



  # checks
  if (sum(forecasts[,1] < forecasts[,2], na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(forecasts, na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= forecasts[,1], na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }

  # also save the forcing variable
  forecasts <- cbind(forecasts,kappa_t)


  # read out parameters
  parameters <- cbind(pnorm(theta[1]),exp(theta[2]), -exp(theta[3])*pnorm(theta[4]), - exp(theta[3]) )

  return(list(Eloss = Eloss, loss = loss, forecasts = forecasts, parameters = parameters))
}


#' @title One-factor GAS model function handle
#' @description This function calculates the average FZ loss of the one factor GAS model similar to PSC (2019)'s MATLAB function GAS_onefactor_LL3.m. It is used as a function handle in the optimization.
#' @import stats
#'
#' @param theta Vector of parameters. theta = [beta,gamma,b,c]. Set beta = pnorm(beta), gamma = exp(gamma), b=-exp(b), c = pnorm(c). omega = 0 for identification
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return one factor GAS model function handle
GAS_onefactor_fh <- function( theta, r, alpha, tau = -1 ){

  # read out the parameters:
  beta  = pnorm(theta[1])
  gamma = exp(theta[2])
  b     = - exp(theta[3])
  c     = pnorm(theta[4])
  a     = c*b

  omega = 0
  kalpha = 1         # normalizing positive constant for Hessian. If numerical optimization fails, use a value smaller than unity

  T <- length(r)

  kappa_t   = array(NaN,dim=c(T,1))
  forecasts = array(NaN,dim=c(T,2))       # quantile and ES

  loss      = array(0,dim=c(T,1))

  # initial values:
  kappa_t[1,1]  = omega/(1-beta)     # unconditional average of kappa
  forecasts[1,] = cbind(a * exp(kappa_t[1,1]), b * exp(kappa_t[1,1]))

  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1]<= forecasts[1,1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-forecasts[1,1])))
    # smooth approximation to the indicator function
  }

  loss[1,1] <- -1/alpha/forecasts[1,2]*hits*(forecasts[1,1]-r[1]) - 1/forecasts[1,2]*(forecasts[1,2]-forecasts[1,1]) + log(-forecasts[1,2])

  for (tt in 2:T){

    hitsL <- hits   # lagged hit variable

    # forcing variable
    kappa_t[tt,1] <- omega + beta * kappa_t[tt-1,1] + gamma/kalpha/forecasts[tt-1,2]*(1/alpha* hitsL *r[tt-1] - forecasts[tt-1,2])

    # forecasts
    forecasts[tt,] = cbind(a * exp(kappa_t[tt,1]), b * exp(kappa_t[tt,1]))

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= forecasts[tt,1])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-forecasts[tt,1])))
      # smooth approximation to the indicator function
    }

    loss[tt,1] <- -1/alpha/forecasts[tt,2]*hits*(forecasts[tt,1]-r[tt]) - 1/forecasts[tt,2]*(forecasts[tt,2]-forecasts[tt,1]) + log(-forecasts[tt,2])

  }

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(forecasts[,1] < forecasts[,2], na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(forecasts, na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= forecasts[,1], na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }


  return(Eloss)
}


#' @title Two-factor GAS model
#' @description This function calculates the average FZ loss of the one factor GAS model similar to PSC (2019)'s MATLAB function GAS_twofactor_LL3.m
#' @import stats
#'
#' @param theta Vector of parameters. theta = [w1,w2,b1,b2,a11,a12,a21,a22]
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return two factor GAS model quantile and ES forecasts
#' @export
GAS_twofactor <- function( theta, r, alpha, tau = -1 ){

  # read out the parameters:
  w1  = theta[1]
  w2  = theta[2]
  b1  = pnorm(theta[3])     # make sure this is in (0,1)
  b2  = pnorm(theta[4])
  a11 = theta[5]
  a12 = theta[6]
  a21 = theta[7]
  a22 = theta[8]

  # assemble into vectors and matrices
  b   = diag(as.vector(rbind(b1,b2)))
  A   = rbind(cbind(a11,a12),cbind(a21,a22))
  w   = rbind(w1,w2)

  T <- length(r)

  # Initial value of forecasts. Use unconditional sample VaR and ES:
  FC0 <- sample_forecasts( r = r, alpha = alpha )
  FC0 <- as.matrix(unlist(FC0))

  lamE = array(NaN,dim=c(T,1))
  lamV = array(NaN,dim=c(T,1))
  loss = array(0,dim=c(T,1))
  forecasts = array(NaN,dim=c(T,2))

  # initial values:
  lamE[1,1] = 1
  lamV[1,1] = 1
  forecasts[1,] = t(FC0)

  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1] <= forecasts[1,1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-forecasts[1,1])))
    # smooth approximation to the indicator function
  }

  loss[1,1] <- -1/alpha/forecasts[1,2]*hits*(forecasts[1,1]-r[1]) - 1/forecasts[1,2]*(forecasts[1,2]-forecasts[1,1]) + log(-forecasts[1,2])

  for (tt in 2:T){

    hitsL <- hits   # lagged hit variable

    # forcing variables
    lamV[tt,1] <- -forecasts[tt-1,1]*(hitsL - alpha)
    lamE[tt,1] <- 1/alpha*hitsL*r[tt-1] - forecasts[tt-1,2]

    forecasts[tt,] <- t(w + b %*% rbind(forecasts[tt-1,1],forecasts[tt-1,2]) + A %*% rbind(lamV[tt,1],lamE[tt,1]))

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= forecasts[tt,1])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-forecasts[tt,1])))
      # smooth approximation to the indicator function
    }

    loss[tt,1] <- -1/alpha/forecasts[tt,2]*hits*(forecasts[tt,1]-r[tt]) - 1/forecasts[tt,2]*(forecasts[tt,2]-forecasts[tt,1]) + log(-forecasts[tt,2])

  }

  # forcing variable
  lam = cbind(lamV,lamE)

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(forecasts[,1] < forecasts[,2], na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(forecasts, na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= forecasts[,1], na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }

  # also save the forcing variable
  forecasts <- cbind(forecasts,lam)


  # read out parameters
  parameters <- rbind( w1, w2, b1, b2, a11, a12, a21, a22 )

  return(list(Eloss = Eloss, loss = loss, forecasts = forecasts, parameters = parameters))
}


#' @title Two-factor GAS model function handle
#' @description This function calculates the average FZ loss of the one factor GAS model similar to PSC (2019)'s MATLAB function GAS_twofactor_LL3.m. This function is used as a function handle.
#' @import stats
#'
#' @param theta Vector of parameters. theta = [w1,w2,b1,b2,a11,a12,a21,a22]
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return two factor GAS model function handle
GAS_twofactor_fh <- function( theta, r, alpha, tau = -1 ){

  # read out the parameters:
  w1  = theta[1]
  w2  = theta[2]
  b1  = pnorm(theta[3])     # make sure this is in (0,1)
  b2  = pnorm(theta[4])
  a11 = theta[5]
  a12 = theta[6]
  a21 = theta[7]
  a22 = theta[8]

  # assemble into vectors and matrices
  b   = diag(as.vector(rbind(b1,b2)))
  A   = rbind(cbind(a11,a12),cbind(a21,a22))
  w   = rbind(w1,w2)

  T <- length(r)

  # Initial value of forecasts. Use unconditional sample VaR and ES:
  FC0 <- sample_forecasts( r = r, alpha = alpha )
  FC0 <- as.matrix(unlist(FC0))

  lamE = array(NaN,dim=c(T,1))
  lamV = array(NaN,dim=c(T,1))
  loss = array(0,dim=c(T,1))
  forecasts = array(NaN,dim=c(T,2))

  # initial values:
  lamE[1,1] = 1
  lamV[1,1] = 1
  forecasts[1,] = t(FC0)

  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1] <= forecasts[1,1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-forecasts[1,1])))
    # smooth approximation to the indicator function
  }

  loss[1,1] <- -1/alpha/forecasts[1,2]*hits*(forecasts[1,1]-r[1]) - 1/forecasts[1,2]*(forecasts[1,2]-forecasts[1,1]) + log(-forecasts[1,2])

  for (tt in 2:T){

    hitsL <- hits   # lagged hit variable

    # forcing variables
    lamV[tt,1] <- -forecasts[tt-1,1]*(hitsL - alpha)
    lamE[tt,1] <- 1/alpha*hitsL*r[tt-1] - forecasts[tt-1,2]

    forecasts[tt,] <- t(w + b %*% rbind(forecasts[tt-1,1],forecasts[tt-1,2]) + A %*% rbind(lamV[tt,1],lamE[tt,1]))

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= forecasts[tt,1])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-forecasts[tt,1])))
      # smooth approximation to the indicator function
    }

    loss[tt,1] <- -1/alpha/forecasts[tt,2]*hits*(forecasts[tt,1]-r[tt]) - 1/forecasts[tt,2]*(forecasts[tt,2]-forecasts[tt,1]) + log(-forecasts[tt,2])

  }

  # forcing variable
  lam = cbind(lamV,lamE)

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(forecasts[,1] < forecasts[,2], na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(forecasts, na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= forecasts[,1], na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }


  return(Eloss)
}


#' @title GARCH model with FZ loss function
#' @description This function calculates the average FZ loss of the GARCH model similar to PSC (2019)'s MATLAB function garch_FZ_LL.m
#' @import stats
#'
#' @param theta Vector of parameters. theta = [beta,gamma,b,c]. Set beta = pnorm(beta), gamma = exp(gamma), b lies in (0,1), a is positive. omega = 1 for identification
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param omega this has to be fixed to a constanz. The default is omega = 1.
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return GARCH model quantile and ES forecasts
#' @export
GARCH_FZ <- function( theta, r, alpha, omega = 1, tau = -1 ){

  # read out the parameters:
  beta  = pnorm(theta[1])
  gamma = exp(theta[2])
  b     = - exp(theta[3])
  c     = pnorm(theta[4])
  a     = c*b


  # unconditional sample variance of return residuals
  h0    = var(r)

  T <- length(r)


  h_t       = array(NaN,dim=c(T,1))       # forcing variable
  forecasts = array(NaN,dim=c(T,2))       # quantile and ES

  loss      = array(0,dim=c(T,1))

  # initial values:
  h_t[1,1]  = h0
  forecasts[1,] = cbind(a * exp(h_t[1,1]), b * exp(h_t[1,1]))

  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1] <= forecasts[1,1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-forecasts[1,1])))
    # smooth approximation to the indicator function
  }

  loss[1,1] <- -1/alpha/forecasts[1,2]*hits*(forecasts[1,1]-r[1]) - 1/forecasts[1,2]*(forecasts[1,2]-forecasts[1,1]) + log(-forecasts[1,2])

  for (tt in 2:T){

    hitsL <- hits   # lagged hit variable

    # forcing variable
    h_t[tt,1] <- omega + beta * h_t[tt-1,1] + gamma * (r[tt-1]^2)

    # forecasts
    forecasts[tt,] = cbind(a * exp(h_t[tt,1]), b * exp(h_t[tt,1]))

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= forecasts[tt,1])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-forecasts[tt,1])))
      # smooth approximation to the indicator function
    }

    loss[tt,1] <- -1/alpha/forecasts[tt,2]*hits*(forecasts[tt,1]-r[tt]) - 1/forecasts[tt,2]*(forecasts[tt,2]-forecasts[tt,1]) + log(-forecasts[tt,2])

  }

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(forecasts[,1] < forecasts[,2], na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(forecasts, na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= forecasts[,1], na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }

  # also save the forcing variable
  forecasts <- cbind(forecasts,h_t)


  # read out parameters
  parameters <- cbind(pnorm(theta[1]),exp(theta[2]), -exp(theta[3])*pnorm(theta[4]), - exp(theta[3]) )

  return(list(Eloss = Eloss, loss = loss, forecasts = forecasts, parameters = parameters))
}


#' @title GARCH model with FZ loss function handle
#' @description This function calculates the average FZ loss of the GARCH model similar to PSC (2019)'s MATLAB function garch_FZ_LL.m. It is used as the function handle.
#' @import stats
#'
#' @param theta Vector of parameters. theta = [beta,gamma,b,c]. Set beta = pnorm(beta), gamma = exp(gamma), b lies in (0,1), a is positive. omega = 1 for identification
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param omega this has to be fixed to a constanz. The default is omega = 1.
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return GARCH model function handle
GARCH_FZ_fh <- function( theta, r, alpha, omega = 1, tau = -1 ){

  # read out the parameters:
  beta  = pnorm(theta[1])
  gamma = exp(theta[2])
  b     = - exp(theta[3])
  c     = pnorm(theta[4])
  a     = c*b


  # unconditional sample variance of return residuals
  h0    = var(r)

  T <- length(r)


  h_t       = array(NaN,dim=c(T,1))       # forcing variable
  forecasts = array(NaN,dim=c(T,2))       # quantile and ES

  loss      = array(0,dim=c(T,1))

  # initial values:
  h_t[1,1]  = h0
  forecasts[1,] = cbind(a * exp(h_t[1,1]), b * exp(h_t[1,1]))

  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1] <= forecasts[1,1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-forecasts[1,1])))
    # smooth approximation to the indicator function
  }

  loss[1,1] <- -1/alpha/forecasts[1,2]*hits*(forecasts[1,1]-r[1]) - 1/forecasts[1,2]*(forecasts[1,2]-forecasts[1,1]) + log(-forecasts[1,2])

  for (tt in 2:T){

    hitsL <- hits   # lagged hit variable

    # forcing variable
    h_t[tt,1] <- omega + beta * h_t[tt-1,1] + gamma * (r[tt-1]^2)

    # forecasts
    forecasts[tt,] = cbind(a * exp(h_t[tt,1]), b * exp(h_t[tt,1]))

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= forecasts[tt,1])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-forecasts[tt,1])))
      # smooth approximation to the indicator function
    }

    loss[tt,1] <- -1/alpha/forecasts[tt,2]*hits*(forecasts[tt,1]-r[tt]) - 1/forecasts[tt,2]*(forecasts[tt,2]-forecasts[tt,1]) + log(-forecasts[tt,2])

  }

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(forecasts[,1] < forecasts[,2], na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(forecasts, na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= forecasts[,1], na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }

  return(Eloss)
}



#' @title AS CARE model
#' @description This function calculates the asymmetric Laplace loss function of the asymmetric slope (AS) CARE model. Compare Taylor (2017).
#'
#' @param theta Vector of parameters. theta = [beta,gamma]. Set beta = [beta1,beta2,beta3,beta4] in the quantile equation and gamma = [gamma1, gamma2, gamma3] in ES equation.
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return AS-CARE model quantile and ES forecasts
#' @export
CARE_AS <- function( theta, r, alpha, tau = -1 ){

  # read out the parameters:
  beta  = theta[1:4]
  gamma = theta[5:7]

  T = length(r)

  Q  <- rep(0,T)
  x  <- rep(0,T)
  ES <- rep(0,T)

  loss <- rep(0,T)

  # initial values: use unconditional sample average
  FC0  <- sample_forecasts( r = r, alpha = alpha )
  FC0  <- as.matrix(unlist(FC0))

  Q[1]  <- FC0[1,1]
  ES[1] <- FC0[2,1]
  x[1]  <- FC0[1,1] - FC0[2,1]


  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1] <= Q[1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-Q[1])))
    # smooth approximation to the indicator function
  }

  loss[1] <- -log((alpha-1)/ES[1]) - ((r[1]-Q[1])*(alpha - hits))/(alpha*ES[1])


  for (tt in 2:T) {
    Q[tt]  <- beta[1] + beta[2]*abs(r[tt-1])*(r[tt-1] > 0) + beta[3]*abs(r[tt-1])*(r[tt-1] <= 0) + beta[4]*Q[tt-1]
    x[tt]  <- ( Q[tt-1] >= r[tt-1]) * ( gamma[1] + gamma[2]* (Q[tt-1] - r[tt-1]) + gamma[2] * x[tt-1]) + ( Q[tt-1] < r[tt-1]) * x[tt-1]
    ES[tt] <- Q[tt] - x[tt]

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= Q[tt])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-Q[tt])))
      # smooth approximation to the indicator function
    }

    loss[tt] <- -log((alpha-1)/ES[tt]) - ((r[tt]-Q[tt])*(alpha - hits))/(alpha*ES[tt])
  }

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(Q < ES, na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(cbind(Q,ES), na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= Q, na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }

  # also save the forcing variable
  forecasts <- cbind(Q,ES)


  # read out parameters
  parameters <- cbind( beta[1], beta[2], beta[3], beta[4], gamma[1], gamma[2], gamma[3] )

  return(list(Eloss = Eloss, loss = loss, forecasts = forecasts, parameters = parameters))
}


#' @title AS CARE model function handle
#' @description This function calculates the asymmetric Laplace loss function of the asymmetric slope (AS) CARE model. Compare Taylor (2017). It is used as the function handle.
#'
#' @param theta Vector of parameters. theta = [beta,gamma]. Set beta = [beta1,beta2,beta3,beta4] in the quantile equation and gamma = [gamma1, gamma2, gamma3] in ES equation.
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return AS-CARE model function handle
CARE_AS_fh <- function( theta, r, alpha, tau = -1 ){

  # read out the parameters:
  beta  = theta[1:4]
  gamma = theta[5:7]

  T = length(r)

  Q  <- rep(0,T)
  x  <- rep(0,T)
  ES <- rep(0,T)

  loss <- rep(0,T)

  # initial values: use unconditional sample average
  FC0  <- sample_forecasts( r = r, alpha = alpha )
  FC0  <- as.matrix(unlist(FC0))

  Q[1]  <- FC0[1,1]
  ES[1] <- FC0[2,1]
  x[1]  <- FC0[1,1] - FC0[2,1]


  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1] <= Q[1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-Q[1])))
    # smooth approximation to the indicator function
  }

  loss[1] <- -log((alpha-1)/ES[1]) - ((r[1]-Q[1])*(alpha - hits))/(alpha*ES[1])


  for (tt in 2:T) {
    Q[tt]  <- beta[1] + beta[2]*abs(r[tt-1])*(r[tt-1] > 0) + beta[3]*abs(r[tt-1])*(r[tt-1] <= 0) + beta[4]*Q[tt-1]
    x[tt]  <- ( Q[tt-1] >= r[tt-1]) * ( gamma[1] + gamma[2]* (Q[tt-1] - r[tt-1]) + gamma[2] * x[tt-1]) + ( Q[tt-1] < r[tt-1]) * x[tt-1]
    ES[tt] <- Q[tt] - x[tt]

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= Q[tt])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-Q[tt])))
      # smooth approximation to the indicator function
    }

    loss[tt] <- -log((alpha-1)/ES[tt]) - ((r[tt]-Q[tt])*(alpha - hits))/(alpha*ES[tt])
  }

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(Q < ES, na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(cbind(Q,ES), na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= Q, na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }


  return(Eloss)
}


#' @title SAV CARE model
#' @description This function calculates the asymmetric Laplace loss function of the symmetric absolute value  (SAV) CARE model. Compare Taylor (2017).
#'
#' @param theta Vector of parameters. theta = [beta,gamma]. Set beta = [beta1,beta2,beta3] in the quantile equation and gamma = [gamma1, gamma2, gamma3] for ES equation.
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return SAV-CARE model quantile and ES forecasts
#' @export
CARE_SAV <- function( theta, r, alpha, tau = -1 ){

  # read out the parameters:
  beta  = theta[1:3]
  gamma = theta[4:6]

  T = length(r)

  Q  <- rep(0,T)
  x  <- rep(0,T)
  ES <- rep(0,T)

  loss <- rep(0,T)

  # initial values: use unconditional sample average
  FC0 <- sample_forecasts( r = r, alpha = alpha )
  FC0 <- as.matrix(unlist(FC0))

  Q[1]  <- FC0[1,1]
  ES[1] <- FC0[2,1]
  x[1]  <- FC0[1,1] - FC0[2,1]


  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1] <= Q[1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-Q[1])))
    # smooth approximation to the indicator function
  }

  loss[1] <- -log((alpha-1)/ES[1]) - ((r[1]-Q[1])*(alpha - hits))/(alpha*ES[1])


  for (tt in 2:T) {
    Q[tt]  <- beta[1] + beta[2]*abs(r[tt-1]) + beta[3]*Q[tt-1]
    x[tt]  <- ( Q[tt-1] >= r[tt-1]) * ( gamma[1] + gamma[2]* (Q[tt-1] - r[tt-1]) + gamma[2] * x[tt-1]) + ( Q[tt-1] < r[tt-1]) * x[tt-1]
    ES[tt] <- Q[tt] - x[tt]

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= Q[tt])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-Q[tt])))
      # smooth approximation to the indicator function
    }

    loss[tt] <- -log((alpha-1)/ES[tt]) - ((r[tt]-Q[tt])*(alpha - hits))/(alpha*ES[tt])
  }

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(Q < ES, na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(cbind(Q,ES), na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= Q, na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }

  # also save the forcing variable
  forecasts <- cbind(Q,ES)


  # read out parameters
  parameters <- cbind( beta[1], beta[2], beta[3], gamma[1], gamma[2], gamma[3] )

  return(list(Eloss = Eloss, loss = loss, forecasts = forecasts, parameters = parameters))
}


#' @title SAV CARE model function handle
#' @description This function calculates the asymmetric Laplace loss function of the symmetric absolute value  (SAV) CARE model. Compare Taylor (2017). It is used as a function handle.
#'
#' @param theta Vector of parameters. theta = [beta,gamma]. Set beta = [beta1,beta2,beta3] in the quantile equation and gamma = [gamma1, gamma2, gamma3] for ES equation.
#' @param r vector of returns
#' @param alpha VaR and ES level
#' @param tau smoothing parameter for indicator function/ hit function. tau = 5, 10, 20, 50 from high to low and tau = -1 (default) no smoothing
#'
#' @return SAV-CARE model function handle
CARE_SAV_fh <- function( theta, r, alpha, tau = -1 ){

  # read out the parameters:
  beta  = theta[1:3]
  gamma = theta[4:6]

  T = length(r)

  Q  <- rep(0,T)
  x  <- rep(0,T)
  ES <- rep(0,T)

  loss <- rep(0,T)

  # initial values: use unconditional sample average
  FC0 <- sample_forecasts( r = r, alpha = alpha )
  FC0 <- as.matrix(unlist(FC0))

  Q[1]  <- FC0[1,1]
  ES[1] <- FC0[2,1]
  x[1]  <- FC0[1,1] - FC0[2,1]


  # hit variable with or without smoothing
  if (tau == -1){
    hits = (r[1] <= Q[1])    # indicator function
  } else {
    hits = 1/(1+exp(tau*(r[1]-Q[1])))
    # smooth approximation to the indicator function
  }

  loss[1] <- -log((alpha-1)/ES[1]) - ((r[1]-Q[1])*(alpha - hits))/(alpha*ES[1])


  for (tt in 2:T) {
    Q[tt]  <- beta[1] + beta[2]*abs(r[tt-1]) + beta[3]*Q[tt-1]
    x[tt]  <- ( Q[tt-1] >= r[tt-1]) * ( gamma[1] + gamma[2]* (Q[tt-1] - r[tt-1]) + gamma[2] * x[tt-1]) + ( Q[tt-1] < r[tt-1]) * x[tt-1]
    ES[tt] <- Q[tt] - x[tt]

    # hit variable with or without smoothing
    if (tau == -1){
      hits = (r[tt] <= Q[tt])    # indicator function
    } else {
      hits = 1/(1+exp(tau*(r[tt]-Q[tt])))
      # smooth approximation to the indicator function
    }

    loss[tt] <- -log((alpha-1)/ES[tt]) - ((r[tt]-Q[tt])*(alpha - hits))/(alpha*ES[tt])
  }

  # average loss
  Eloss = mean(loss)


  # checks
  if (sum(Q < ES, na.rm = TRUE) > 0 ){
    Eloss <- 1e6
    # ES < VaR for homogenous zero degree oss function
  }
  if (max(max(cbind(Q,ES), na.rm = TRUE)) >= 0){
    Eloss <- 1e6
    # both forecasts have to be negative
  }
  if (sum(r <= Q, na.rm = TRUE) == 0){
    Eloss <- 1e6
    # there should be at least some hits
  }

  return(Eloss)
}


#' @title Sample VaR and ES forecast
#' @description This function calculates the sample quantile and ES from the data. This can be used as a starting value. Similar to PZC(2019) sample_VE.m MATLAB function.
#'
#' @param r return time series
#' @param alpha VaR and ES probability level
#'
#' @return sample quantile and ES
#' @export
sample_forecasts <- function( r, alpha ){

  T <- length(r)

#  N <- ncol(r)
#  p <- length(alpha)

#  forecasts = array(NaN, dim = c(N,2,p))

#  for (ii in 1:N) {
#    data = r[,ii]
#    for (aa in 1:p){
#      forecasts[ii,1,aa] = quantile(data,alpha[aa])
#      if (sum( ( data <= forecasts[ii,1,aa])) > 0){
#        forecasts[ii,2,aa] = mean(data[data<=forecasts[ii,1,aa]])
#      }
#    }
#  }

  forecasts = array(NaN,dim = c(2,1))

  forecasts[1,1] <- quantile(r,alpha, na.rm = TRUE)
  forecasts[2,1] <- mean(r[r<=forecasts[1,1]], na.rm = TRUE)

  return(list(forecasts = forecasts))

}


#' @title Optimize models
#' @description Estimate the model and generate in-sample fits and OOS forecasts
#' @import rugarch stats
#'
#' @param r return time series
#' @param alpha quantile and ES probability level
#' @param model Specifies the model that is estimated and OOS forecasts are obtained. Choices are "GARCHGZ", "GAS1F", "CAREAS" and "CARESAV"
#' @param splitdate cut off between IS and OOS
#'
#' @return VEall
#' @export
Optimize_models <- function( r, alpha, model, splitdate = round(1/3*length(r),0),
                             in_sample = FALSE ){

  T <- length(r)

  # separate data in an in- (IS) and out-of-sample (OOS) period
  if (in_sample) {
    data_IS  <- r
    data_OOS <- r
    m <- 0
    n <- T
  } else {
    data_IS <- r[1:splitdate]
    
    m <- length(data_IS)
    n <- T - m
    
    data_OOS <- r[-(1:m)]
  }
  # initialize forecast vector
  VEall <- array(NaN,dim = c(n,2))

  if (model == "GARCHFZ"){

  # get starting values: estimate a GARCH model on the residuals

  residual = data_IS - mean(data_IS)

  # fit GARCH model
  GARCH_spec <- ugarchspec(mean.model = list(include.mean = FALSE,armaOrder=c(0,0)) )
  GARCH_fit  <- ugarchfit(spec=GARCH_spec, residual)

  # simulate the forcing variable from the GARCH process
  GARCHparam = as.vector(GARCH_fit@fit$coef)


  # estimate GARCH model with FZ loss function:

  GARCHFZparam <- array(NaN, dim = c(4,1))

  FCbar <- sample_forecasts( r = residual/rugarch::sigma(GARCH_fit), alpha = alpha )
  FCbar <- as.matrix(unlist(FCbar))

  # initial parameter values:
  theta0 <- cbind(qnorm(GARCHparam[3]), log(GARCHparam[2]), log(-FCbar[2,1]),qnorm(FCbar[1,1]/FCbar[2,1]))

  # Optimize the GARCH via the FZ loss - create a function handle and optimize model.
  # Decrease the amount of smoothing in every step.
  fun1 <- function(theta) suppressWarnings(GARCH_FZ_fh( theta = theta, r = residual, alpha = alpha, omega = GARCHparam[1], tau = 5 ))
  fit1 <- try(stats::optim(par = theta0, fn = fun1, method = "Nelder-Mead"), silent = TRUE)

  fun2 <- function(theta) suppressWarnings(GARCH_FZ_fh( theta = theta, r = residual, alpha = alpha, omega = GARCHparam[1], tau = 20 ))
  fit2 <- try(stats::optim(par = fit1$par, fn = fun2, method = "Nelder-Mead"), silent = TRUE)

  fun3 <- function(theta) suppressWarnings(GARCH_FZ_fh( theta = theta, r = residual, alpha = alpha, omega = GARCHparam[1] ))
  fit3 <- try(stats::optim(par = fit2$par, fn = fun3, method = "Nelder-Mead"), silent = TRUE)

  GARCHFZparam = fit3$par

  # get OOS forecasts: forecasts
  GARCH_FZ_results <- GARCH_FZ( theta = GARCHFZparam, r = r - mean(data_OOS), alpha = alpha, omega = GARCHparam[1])
  mu_OOS = mean(data_OOS)

  # get OOS VaR and ES forecasts from GARCH_FZ model
  VEall = mu_OOS + GARCH_FZ_results$forecasts[-(1:m),1:2]
  # VEall = GARCH_FZ_results$forecasts[-(1:m),1:2]


  parameters = GARCH_FZ_results$parameters

  } else if (model == "GAS1F"){

  # estimate one-factor GAS model with FZ loss function:

  GAS1FZparam <- array(NaN, dim = c(4,1))

  FCbar <- sample_forecasts( r = data_OOS, alpha = alpha )
  FCbar <- as.matrix(unlist(FCbar))

  # initial parameter values:
  theta0 = cbind(qnorm(0.95),log(0.005),log(-FCbar[2,1]),qnorm(FCbar[1,1]/FCbar[2,1]))

  # Optimize the GARCH via the FZ loss - create a function handle and optimize model.
  # Decrease the amount of smoothing in every step.
  fun1 <- function(theta) suppressWarnings(GAS_onefactor_fh( theta = theta, r = data_OOS, alpha = alpha, tau = 5 ))
  fit1 <- try(stats::optim(par = theta0, fn = fun1, method = "Nelder-Mead"), silent = TRUE)

  fun2 <- function(theta) suppressWarnings(GAS_onefactor_fh( theta = theta, r = data_OOS, alpha = alpha, tau = 20 ))
  fit2 <- try(stats::optim(par = fit1$par, fn = fun2, method = "Nelder-Mead"), silent = TRUE)

  fun3 <- function(theta) suppressWarnings(GAS_onefactor_fh( theta = theta, r = data_OOS, alpha = alpha ))
  fit3 <- try(stats::optim(par = fit2$par, fn = fun3, method = "Nelder-Mead"), silent = TRUE)

  GAS1FZparam = fit3$par

  # get OOS forecasts: forecasts
  GAS1FZ_results <- GAS_onefactor( theta = GAS1FZparam, r = r , alpha = alpha )

  # get OOS VaR and ES forecasts from GAS1FZ model
  if (in_sample) {
    VEall = GAS1FZ_results$forecasts[,1:2]
  } else {
    VEall = GAS1FZ_results$forecasts[-(1:m),1:2]
  }

  parameters = GAS1FZ_results$parameters

  } else if (model == "GAS2F"){

    # estimate one-factor GAS model with FZ loss function:

    GAS2FZparam <- array(NaN, dim = c(8,1))

    FCbar <- sample_forecasts( r = data_OOS, alpha = alpha )
    FCbar <- as.matrix(unlist(FCbar))

    # try three different initial parameter values:
    theta0_vec = rbind(cbind((1-0.97)*FCbar[1,1],(1-0.97)*FCbar[2,1],qnorm(0.97),qnorm(0.97),0,0,0,0),cbind((1-0.98)*FCbar[1,1],(1-0.98)*FCbar[2,1],qnorm(0.98),qnorm(0.98),0,0,0,0),cbind((1-0.99)*FCbar[1,1],(1-0.99)*FCbar[2,1],qnorm(0.99),qnorm(0.99),0,0,0,0))

    # Optimize the GARCH via the FZ loss - create a function handle and optimize model.
    # Decrease the amount of smoothing in every step.
    param_out <- array(NaN, dim = c(nrow(theta0_vec),ncol(theta0_vec)+1))


    for (values in 1:nrow(theta0_vec)){

      theta0 <- as.vector(theta0_vec[values,])

      # try this starting value
      tau_init = GAS_twofactor(theta = theta0, r = data_OOS, alpha = alpha, tau = 5)
      if (tau_init$Eloss < 1e6){

        fun1 <- function(theta) suppressWarnings(GAS_twofactor_fh( theta = theta, r = data_OOS, alpha = alpha, tau = 5 ))
        fit1 <- try(stats::optim(par = theta0, fn = fun1, method = "Nelder-Mead"), silent = TRUE)
        theta1 <- fit1$par

      } else {
        theta1 <- theta0
      }

      tau_5 = GAS_twofactor(theta = theta1, r = data_OOS, alpha = alpha, tau = 20)
      if (tau_5$Eloss < 1e6){

        fun2 <- function(theta) suppressWarnings(GAS_twofactor_fh( theta = theta, r = data_OOS, alpha = alpha, tau = 20 ))
        fit2 <- try(stats::optim(par = theta1, fn = fun2, method = "Nelder-Mead"), silent = TRUE)
        theta2 <- fit2$par

      } else {
        theta2 <- theta0
      }

      tau_20 = GAS_twofactor(theta = theta2, r = data_OOS, alpha = alpha)
      if (tau_20$Eloss < 1e6){

        fun3 <- function(theta) suppressWarnings(GAS_twofactor_fh( theta = theta, r = data_OOS, alpha = alpha ))
        fit3 <- try(stats::optim(par = theta2, fn = fun3, method = "Nelder-Mead"), silent = TRUE)
        theta3 <- fit3$par

      } else {
        theta3 <- theta0
      }

      tau_0 = GAS_twofactor(theta = theta3, r = data_OOS, alpha = alpha)
      GAS2param = as.vector(theta3)

      out <- c(GAS2param,tau_0$Eloss)
      param_out[values,] <- out

    }

    idx <- which.min(param_out[,ncol(param_out)])

    GAS2FZparam <- as.vector(param_out[idx,-ncol(param_out)])


    # get OOS forecasts: forecasts
    GAS2FZ_results <- GAS_twofactor( theta = GAS2FZparam, r = r , alpha = alpha )

    # get OOS VaR and ES forecasts from GAS1FZ model
    if (in_sample) {
      VEall = GAS2FZ_results$forecasts[,1:2]
    } else {
      VEall = GAS2FZ_results$forecasts[-(1:m),1:2]
    }

    parameters = GAS2FZ_results$parameters

  } else if (model == "CAREAS"){

  # estimate AS-CARE model with asymmetric Laplace loss function:

  ASCAREparam <- array(NaN, dim = c(7,1))

  # initial parameter values:
  beta_AS  =  cbind(-0.0003,-0.05,-0.15,0.8)
  gamma_AS = cbind(0.00017, 0.125,0.84)
  theta0   = cbind(beta_AS,gamma_AS)

  # Optimize the GARCH via the FZ loss - create a function handle and optimize model.
  # Decrease the amount of smoothing in every step.
  fun1 <- function(theta) suppressWarnings(CARE_AS_fh( theta = theta, r = r, alpha = alpha, tau = 5 ))
  fit1 <- try(stats::optim(par = theta0, fn = fun1, method = "Nelder-Mead"), silent = TRUE)

  fun2 <- function(theta) suppressWarnings(CARE_AS_fh( theta = theta, r = r, alpha = alpha, tau = 20 ))
  fit2 <- try(stats::optim(par = fit1$par, fn = fun2, method = "Nelder-Mead"), silent = TRUE)

  fun3 <- function(theta) suppressWarnings(CARE_AS_fh( theta = theta, r = r, alpha = alpha ))
  fit3 <- try(stats::optim(par = fit2$par, fn = fun3, method = "Nelder-Mead"), silent = TRUE)

  ASCAREparam = fit3$par

  # get OOS forecasts: forecasts
  ASCARE_results <- CARE_AS( theta = ASCAREparam, r = r , alpha = alpha )

  # get OOS VaR and ES forecasts from AS-CARE model
  if (in_sample) {
    VEall = ASCARE_results$forecasts[,1:2]
  } else {
    VEall = ASCARE_results$forecasts[-(1:m),1:2]
  }

  parameters = ASCARE_results$parameters


  } else if (model == "CARESAV"){

  # estimate SAV-CARE model with asymmetric Laplace loss function:

  SAVCAREparam <- array(NaN, dim = c(6,1))


  # initial parameter values:
  beta_SAV  =  cbind(-0.0003,-0.1,0.8)
  gamma_SAV = cbind(0.00017, 0.125,0.84)
  theta0    = cbind(beta_SAV,gamma_SAV)

  # Optimize the GARCH via the FZ loss - create a function handle and optimize model.
  # Decrease the amount of smoothing in every step.
  fun1 <- function(theta) suppressWarnings(CARE_SAV_fh( theta = theta, r = data_OOS, alpha = alpha, tau = 5 ))
  fit1 <- try(stats::optim(par = theta0, fn = fun1, method = "Nelder-Mead"), silent = TRUE)

  fun2 <- function(theta) suppressWarnings(CARE_SAV_fh( theta = theta, r = data_OOS, alpha = alpha, tau = 20 ))
  fit2 <- try(stats::optim(par = fit1$par, fn = fun2, method = "Nelder-Mead"), silent = TRUE)

  fun3 <- function(theta) suppressWarnings(CARE_SAV_fh( theta = theta, r = data_OOS, alpha = alpha ))
  fit3 <- try(stats::optim(par = fit2$par, fn = fun3, method = "Nelder-Mead"), silent = TRUE)

  SAVCAREparam = fit3$par

  # get OOS forecasts: forecasts
  SAVCARE_results <- CARE_SAV( theta = SAVCAREparam, r = r , alpha = alpha )

  # get OOS VaR and ES forecasts from SAV-CARE model
  if (in_sample) {
    VEall = SAVCARE_results$forecasts[,1:2]
  } else {
    VEall = SAVCARE_results$forecasts[-(1:m),1:2]
  }

  parameters = SAVCARE_results$parameters

  }

 return(list(OOSforecasts = VEall, parameters = parameters))
}
