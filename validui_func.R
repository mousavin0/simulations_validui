#* FUNCTIONS*####

simulation_ui <- function(FUN = generatingdata_ui, n,
                          rep = 500, seed = 1127, rho, vecy, vect,
                          param = "ACE", trueValue = 0, subset = "refit", linear = FALSE, sigma_correction = TRUE) {
  # now set the seed before finding the true value of the parameter and the bias in the main code
  # set.seed(seed)
  ui <- list()
  pb <- tkProgressBar(min = 0, max = rep, initial = 0)
  starttime <- Sys.time()
  cov <- 0
  for (i in 1:rep) {

    # FUN=generatingdata_ui
    data <- FUN(n = n, rho = rho, vecy = vecy, vect = vect, linear)
    ui[[i]] <- ui.causal(
      Y = data$Y, T = data$T, X = as.data.frame(data$X), rho1 = rho, rho0 = rho, param = param,
      subset = subset, regularization_alpha = 1, sigma_correction = sigma_correction
    )
    if (trueValue >= ui[[i]]$DR$ci[1] & trueValue <= ui[[i]]$DR$ci[2]) {
      cov <- cov + 1
    }
    print(cov / i)
    flush.console()
    setTkProgressBar(pb, i, label = paste(
      "n=", n, " ", "rho=", rho, " \n", i,
      "of", rep, "done", " started at",
      starttime
    ))
  }
  close(pb)

  return(list(ui = ui, cov = cov))
}

# #lowdimensional data
# generatingdata_ui2= function(n,rho,vecy,vect,xdim,ycoef,linear=FALSE,sigma0=1,sigma1=1,xvar=1){# Generating data
#  l=xdim
#   X=rmvn(n,rep(0,l),diag(xvar,l))
#   X=cbind(X,X[,1]^2,X[,1]^3,X[,1]^4,
#                   X[,1]^5,X[,1]^6,X[,1]^7,X[,1]^8,X[,1]^9,X[,1]^10)
#   error<-rmvn(n,c(0,0,0),matrix(c(1,rho*sigma0,rho*sigma1,rho*sigma0,sigma0^2,rho*rho*sigma0*sigma1,rho*sigma1,rho*rho*sigma0*sigma1,sigma1^2),ncol=3))
#   tstar=X[,1]* vect+error[,1]
#   T=as.numeric(tstar>=0)
#
#   if(linear){
#     Y0=1+X[,1]* vecy+X[,(l+1):(l+9)]%*%ycoef +error[,2]
#     Y1=2+X[,1]* vecy+X[,(l+1):(l+9)]%*%ycoef+error[,3]
#   }else{
#     ###E(Y(t)|X,T=t) which is equal to E(Y(t)|X)+rho*sigma*lambda must be appoximate linear
#     Y0=1+X[,1]* vecy+X[,(l+1):(l+9)]%*%ycoef+rho*sigma0*lambda0(X[,1]* vect)+error[,2]
#     # (2^-2*X[,1]^2+3^-3*X[,1]^3+4^-4*X[,1]^4+
#     #                          5^-5*X[,1]^5+6^-6*X[,1]^6+7^-7*X[,1]^7+
#     #                          8^-8*X[,1]^8+9^-9*X[,1]^9+10^-10*X[,1]^10)+
#
#     Y1=2+X[,1]* vecy+X[,(l+1):(l+9)]%*%ycoef-rho*sigma1*lambda1(X[,1]* vect)+error[,3]
# }
#
#   Y <- ifelse(T == 1, Y1, Y0)
#   return(list(X=X,Y=Y,T=T,Y1=Y1,Y0=Y0))
# }
# lowdimensional data
generatingdata_ui <- function(n, rho, vecy, vect, linear = FALSE, sigma0 = 1, sigma1 = 1, xvar = 1) { # Generating data
  # vec=seq(0.2,0.01,by=-0.01)
  l <- max(length(vect), length(vecy))
  # X=rmvn(n,rep(0,l),diag(xvar,l))
  X <- matrix(rnorm(n * l, mean = 0, sd = xvar), n, l)
  error <- rmvn(n, c(0, 0, 0), matrix(c(1, rho * sigma0, rho * sigma1, rho * sigma0, sigma0^2, rho * rho * sigma0 * sigma1, rho * sigma1, rho * rho * sigma0 * sigma1, sigma1^2), ncol = 3))
  gX <- X %*% vect
  mX <- X %*% vecy
  tstar <- gX + error[, 1]
  T <- as.numeric(tstar >= 0)

  if (linear) {
    Y0 <- 1 + mX + error[, 2]
    Y1 <- 2 + mX + error[, 3]
  } else {
    ### E(Y(t)|X,T=t) which is equal to E(Y(t)|X)+rho*sigma*lambda must be appoximate linear
    Y0 <- 1 + mX + rho * sigma0 * lambda0(gX) + error[, 2]
    Y1 <- 2 + mX - rho * sigma1 * lambda1(gX) + error[, 3]
  }





  Y <- ifelse(T == 1, Y1, Y0)
  # cov.sel.high.lasso(Y[T==1],X[T==1,])

  return(list(X = X, Y = Y, T = T, Y1 = Y1, Y0 = Y0))
}

generatingdata_ui_high <- function(n, rho, vecy, vect = 0.3, linear = FALSE, sigma0 = 1, sigma1 = 1, xvar = 1) { # Generating data
  p <- n
  # X=rmvn(n,rep(0,p),diag(xvar,p))
  X <- matrix(rnorm(n * p, mean = 0, sd = xvar), n, p)
  error <- rmvn(n, c(0, 0, 0), matrix(c(1, rho * sigma0, rho * sigma1, rho * sigma0, sigma0^2, rho * rho * sigma0 * sigma1, rho * sigma1, rho * rho * sigma0 * sigma1, sigma1^2), ncol = 3))

  gamma <- vect * c(
    1, 1 / 2, 1 / 3, 1 / 4, 1 / 5,
    rep(1, 5),
    rep(0, p - 10)
  )
  b <- c(
    1, 1 / 2, 1 / 3, 1 / 4, 1 / 5,
    1, 1 / 2, 1 / 3, 1 / 4, 1 / 5,
    rep(0, p - 10)
  )

  gX <- X %*% gamma
  mX <- X %*% b * vecy
  tstar <- gX + error[, 1]
  T <- as.numeric(tstar >= 0)


  if (linear) {
    Y0 <- 1 + mX + error[, 2]
    Y1 <- 2 + mX + error[, 3]
  } else {
    Y0 <- 1 + mX + rho * sigma0 * lambda0(gX) + error[, 2]
    Y1 <- 2 + mX - rho * sigma1 * lambda1(gX) + error[, 3]
  }

  Y <- ifelse(T == 1, Y1, Y0)
  #

  return(list(X = X, Y = Y, T = T, Y1 = Y1, Y0 = Y0))
}









generatingdata_ui_logistic <- function(n, rho, vecy, vect, linear = FALSE, sigma0 = 1, sigma1 = 1, xvar = 1) { # Generating data

  l <- max(length(vect), length(vecy))
  X <- matrix(rnorm(n * l, mean = 0, sd = xvar), n, l)
  
  eta <- rlogis(n,0,1)
  xi <- rho*sigma1*eta + rnorm(n,0,sigma1)
  
  gX <- X %*% vect
  mX <- X %*% vecy
  tstar <- gX + eta
  T <- as.numeric(tstar >= 0)
  Y1 <- 2 + mX - rho * sigma1 * lambda1(gX) + xi
  #Y0 is set to zero, doesnt matter since we are interested in y1 only
  Y <- ifelse(T == 1, Y1,0)
  
  return(list(X = X, Y = Y, T = T, Y1 = Y1, Y0 = Y0))
}

generatingdata_ui_high_logistic <- function(n, rho, vecy, vect = 0.3, linear = FALSE, sigma0 = 1, sigma1 = 1, xvar = 1) { # Generating data
  p <- n
  X <- matrix(rnorm(n * p, mean = 0, sd = xvar), n, p)
  
  eta <- rlogis(n,0,1)
  xi <- rho*sigma1*eta + rnorm(n,0,sigma1)
  
  gamma <- vect * c(
    1, 1 / 2, 1 / 3, 1 / 4, 1 / 5,
    rep(1, 5),
    rep(0, p - 10)
  )
  b <- c(
    1, 1 / 2, 1 / 3, 1 / 4, 1 / 5,
    1, 1 / 2, 1 / 3, 1 / 4, 1 / 5,
    rep(0, p - 10)
  )
  
  gX <- X %*% gamma
  mX <- X %*% b * vecy
  tstar <- gX + eta
  T <- as.numeric(tstar >= 0)
  
  Y1 <- 2 + mX - rho * sigma1 * lambda1(gX) + xi
  #Y0 is set to zero, doesnt matter since we are interested in y1 only
  Y <- ifelse(T == 1, Y1, 0)

  return(list(X = X, Y = Y, T = T, Y1 = Y1, Y0 = Y0))
}









results_ui <- function(uilist, trueValue, true_b = 0, var_detect_percent = c(3, 4)) {
  len <- length(uilist)
  q <- qnorm(0.975)
  ui_lower <- c()
  ui_upper <- c()
  confounding_bias_dr <- c()
  confounding_bias_or <- c()
  estimated_bias_dr <- c()
  estimated_bias_or <- c()
  naivDR <- c()
  est_se <- c()
  XYhat <- matrix(nrow = len, ncol = length(var_detect_percent))
  XThat <- matrix(nrow = len, ncol = length(var_detect_percent))
  # est_b=c()
  for (i in 1:len) {
    r <- uilist[[i]]
    ui_lower[i] <- r$DR$ui[1]
    ui_upper[i] <- r$DR$ui[2]

    naivDR[i] <- r$DR$naiv
    confounding_bias_dr[i] <- r$DR$naiv - trueValue
    confounding_bias_or[i] <- r$OR$naiv - trueValue

    estimated_bias_dr[i] <- r$DR$conf.bias[1]
    estimated_bias_or[i] <- r$OR$conf.bias[1]


    est_se[i] <- r$DR$se

    XYhat[i, ] <- is.element(var_detect_percent, r$XYhat)
    XThat[i, ] <- is.element(var_detect_percent, r$XThat)

    # est_b[i]=r$conf.biascor
  }


  return(list(
    XYhat = colSums(XYhat) / len,
    XThat = colSums(XThat) / len,
    confounding_bias_or = confounding_bias_or,
    confounding_bias_dr = confounding_bias_dr,
    estimated_bias_or = estimated_bias_or,
    estimated_bias_dr = estimated_bias_dr,
    cov = c(
      sum((ui_lower) <= trueValue & (ui_upper) >= trueValue) / len,
      sum((naivDR - true_b - q * est_se) <= trueValue & (naivDR - true_b + q * est_se) >= trueValue) / len
    ),
    
    #CHANGES AFTER REVISION
    tauhat = (ui_lower + ui_upper) / 2,
    tauhat_trueb = naivDR - true_b,
    
    # confbias=b,
    # biasbhat=mean(est_b)-b,
    # sdbhat=sd(est_b),
    sdtauhat = sd((ui_lower + ui_upper) / 2),
    sdhattauhat = mean(est_se)
  ))
}
