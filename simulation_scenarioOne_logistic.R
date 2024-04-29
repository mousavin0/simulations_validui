# load("scenarioOne.Rdata")
# install R packages for ubuntu: https://cran.r-project.org/bin/linux/ubuntu/
install.packages('devtools')
library(devtools)
# install("C:\\Users\\Nimo0033\\OneDrive\\PhD\\validui\\code\\packages\\ui.causal")
# document("C:\\Users\\Nimo0033\\OneDrive\\PhD\\validui\\code\\packages\\ui.causal")
# install("C:\\Users\\Nimo0033\\OneDrive\\PhD\\validui\\code\\packages\\ui")
# document("C:\\Users\\Nimo0033\\OneDrive\\PhD\\validui\\code\\packages\\ui")

install_github('https://github.com/stat4reg/hdim.ui.git')
library(hdim.ui)


#CHANGES AFTER REVISION
#setwd('/home/sajjad/Neli/Rfiles_papertwo/')
#source("validui_func.R")
source("/Users/niloofarmoosavi/Library/CloudStorage/OneDrive-Personal/PhD/validui/code/final/Paper Two-final/validui_func.R")



library(glmnet)
library(MASS)
library(tcltk)
library(hdm)
library(mgcv)





# First Scenario #####
simlow <- list()
nvec <- c(1000, 5000, 20000)

#CHANGED AFTER REVISION
rhovec <- seq(0.8, 0.3, by = -0.1)
#rhovec=0.2
alphavec <- c(-0.07, -0.035, -0.015)
vect <- c(0, 1)

rep <- 500




## lowdim#########################




#####
# install.packages("ramify")
# library(ramify)
# cov_scen1_doubleselection_sigCor_notheory=resize(cov_scen1_doubleselection_sigCor_notheory,3,6,across="columns")
# cov_scen1_refit_trueb_farrellMinusConstant=resize(cov_scen1_refit_trueb_farrellMinusConstant,3,6,across="columns")
# cov_scen1_refit_sigNotCor_ours=resize(cov_scen1_refit_sigNotCor_ours,3,6,across="columns")

cov_scen1_refit_trueb_farrellMinusConstant <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))


#
# cov_scen1_noselection_sigNotCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
# cov_scen1_noselection_sigCor_minna=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
#
cov_scen1_refit_sigNotCor_ours <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))

# cov_scen1_refit_sigCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
#
# cov_scen1_doubleselection_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
# cov_scen1_doubleselection_sigNotCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
cov_scen1_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))


#CHANGES AFTER REVISION
bias_scen1_refit_trueb_farrellMinusConstant <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
bias_scen1_refit_sigNotCor_ours <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
bias_scen1_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))



sigma1 <- 1


library('tcltk')

set.seed(1127)
# Finding true values of the parameter####
pa <- c()
el <- c()
n <- 10000
param <- matrix(NA, nrow = length(alphavec), ncol = length(rhovec))
elambda <- 0
for (i in 1:length(alphavec)) {
  for (j in 1:length(rhovec)) {
    for (k in 1:rep) {
      rho <- rhovec[j]
      vecy <- c(-0.5, alphavec[i])
      vect <- vect
      sigma1 <- 1
      xvar <- 1
      l <- max(length(vect), length(vecy))
      X <- matrix(rnorm(n * l, mean = 0, sd = xvar), n, l)
      gX <- X %*% vect
      lam <- lambda1(gX)
      
      
      # eta <- rlogis(n,0,1)
      # xi <- rho*sigma1*eta + rnorm(n,0,sigma1)
      # 
      # gX <- X %*% vect
      # mX <- X %*% vecy
      # tstar <- gX + eta
      # T <- as.numeric(tstar >= 0)
      # Y1 <- 2 + mX - rho * sigma1 * lambda1(gX) + xi
      
      
      
      mX <- 2 + X %*% vecy - rho * sigma1 * lam
      #value of the parameter when eta is logistic
      pa[k] <- mean(mX + rho*sigma1* rlogis(n, 0, 1))
      if (i == 1 & j == 1) {
        el[k] <- mean(lam)
      }
    }
    param[i, j] <- mean(pa)
  }
}
elambda <- mean(el)

# main simulation
for (i in 1:(length(alphavec))) {
  for (j in 1:length(rhovec)) {
    simlow <- simulation_ui(
      FUN = generatingdata_ui_logistic, rep = rep, n = nvec[i], rho = rhovec[j],
      vecy = c(-0.5, alphavec[i]), vect = vect,
      trueValue = param[i, j], param = "Y1",
      subset = "refit", sigma_correction = 'non'
    )
    res <- results_ui(simlow$ui, param[i, j], true_b = rhovec[j] * sigma1 * elambda)
    cov_scen1_refit_sigNotCor_ours[i, j] <- res$cov[1]
    cov_scen1_refit_trueb_farrellMinusConstant[i, j] <- res$cov[2]
    
    
    #CHANGES AFTER REVISION
    bias_scen1_refit_trueb_farrellMinusConstant[i, j] <-  mean(res$tauhat_trueb) - param[i, j]
    bias_scen1_refit_sigNotCor_ours[i, j] <- mean(res$tauhat) - param[i, j]
    
    
    rm(simlow)
    gc()


    # simlow=simulation_ui(FUN=generatingdata_ui_logistic,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=c(-0.5,alphavec[i]),vect=vect,
    #                      trueValue=param[i,j],param="Y1",
    #                      subset="refit",sigma_correction='old')
    # res=results_ui(simlow$ui,param[i,j],true_b=rhovec[j]*sigma1*elambda)
    # cov_scen1_refit_sigCor[i,j]=res$cov[1]
    # rm(simlow)
    # gc()
    #
    #
    # simlow=simulation_ui(FUN=generatingdata_ui_logistic,rep=rep,n=nvec[i],rho=rhovec[j],
    #                                   vecy=c(-0.5,alphavec[i]),vect=vect,
    #                                   trueValue=param[i,j],param="Y1",
    #                      subset="noselection", sigma_correction='old')
    #   res=results_ui(simlow$ui,param[i,j],true_b=rhovec[j]*sigma1*elambda)
    #   cov_scen1_noselection_sigCor_minna[i,j]=res$cov[1]
    #   rm(simlow)
    #   gc()
    #
    #
    #
    #   simlow=simulation_ui(FUN=generatingdata_ui_logistic,rep=rep,n=nvec[i],rho=rhovec[j],
    #                        vecy=c(-0.5,alphavec[i]),vect=vect,
    #                        trueValue=param[i,j],param="Y1",
    #                        subset="noselection", sigma_correction=FALSE)
    #   res=results_ui(simlow$ui,param[i,j],true_b=rhovec[j]*sigma1*elambda)
    #   cov_scen1_noselection_sigNotCor[i,j]=res$cov[1]
    #   rm(simlow)
    #   gc()


    # simlow <- simulation_ui(
    #   FUN = generatingdata_ui_logistic, rep = rep, n = nvec[i], rho = rhovec[j],
    #   vecy = c(-0.5, alphavec[i]), vect = vect,
    #   trueValue = param[i, j], param = "Y1",
    #   subset = "double", sigma_correction='old'
    # )
    # res <- results_ui(simlow$ui, param[i, j], true_b = rhovec[j] * sigma1 * elambda)
    # cov_scen1_doubleselection_sigCor_notheory[i, j] <- res$cov[1]
    # rm(simlow)
    # gc()



    simlow <- simulation_ui(
      FUN = generatingdata_ui_logistic, rep = rep, n = nvec[i], rho = rhovec[j],
      vecy = c(-0.5, alphavec[i]), vect = vect,
      trueValue = param[i, j], param = "Y1",
      subset = "refit", sigma_correction='new'
    )
    res <- results_ui(simlow$ui, param[i, j], true_b = rhovec[j] * sigma1 * elambda)
    cov_scen1_new_sigCor_notheory[i, j] <- res$cov[1]
    
    #CHANGES AFTER REVISION
    bias_scen1_new_sigCor_notheory[i, j] <- mean(res$tauhat) -  param[i, j]
    
    
    rm(simlow)
    gc()




    #
    # simlow=simulation_ui(FUN=generatingdata_ui_logistic,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=c(-0.5,alphavec[i]),vect=vect,
    #                      trueValue=param[i,j],param="Y1",
    #                      subset="double",sigma_correction=FALSE)
    # res=results_ui(simlow$ui,param[i,j],true_b=rhovec[j]*sigma1*elambda)
    # cov_scen1_doubleselection_sigNotCor[i,j]=res$cov[1]
    # rm(simlow)
    # gc()
  }
}


#CHANGES AFTER REVISION
save.image("/Users/niloofarmoosavi/Library/CloudStorage/OneDrive-Personal/PhD/validui/code/final/Paper Two-final/scenarioOneRevisionLogistic.Rdata")

























# save.image("scenarioOne.Rdata")
#
# # load("scenarioOne.Rdata")

#CHANGED AFTER REVISION
#save.image("scenarioOne02.Rdata")
# load("scenarioOne02.Rdata")

library(xtable)
# temp=cov_scen2_refit_trueb_farrellMinusConstant
# temp <- cov_scen2_refit_sigNotCor_ours

#CHANGED AFTER REVISION
temp <- cov_scen1_refit_trueb_farrellMinusConstant
temp <- bias_scen1_refit_sigNotCor_ours
temp <- cov_scen1_new_sigCor_notheory
rownames(temp) <- nvec
colnames(temp) <- rhovec
xtable(temp)





cov_scen1_refit_trueb_farrellMinusConstant <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
cov_scen1_refit_sigNotCor_ours <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
cov_scen1_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))


