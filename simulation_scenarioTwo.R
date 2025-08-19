
# load("scenarioTwo.Rdata")

library(devtools)
# install("C:\\Users\\...\\validui\\code\\packages\\ui.causal")
# document("C:\\Users\\...\\validui\\code\\packages\\ui.causal")
# install("C:\\Users\\...\\validui\\code\\packages\\ui")
# document("C:\\Users\\...\\validui\\code\\packages\\ui")


install_github('https://github.com/stat4reg/hdim.ui.git')
library(hdim.ui)

#CHANGES AFTER REVISION
#setwd('')
#source("validui_func.R")
source("validui_func.R")



library(glmnet)
library(MASS)
library(tcltk)
library(hdm)
library("mgcv")
# Second Scenario #####
simlow <- list()
nvec <- c(1000, 5000, 20000)

#CHANGED AFTER REVISION
rhovec <- seq(0.8, 0.3, by = -0.1)
#rhovec=0.2
vecy <- c(0, -0.5)
alphavec <- c(0.12, 0.075, 0.06)

rep <- 500
# ##lowdim#########################

#####
# library(ramify)
# cov_scen2_doubleselection_sigCor_notheory=resize(cov_scen2_doubleselection_sigCor_notheory,3,6,across="columns")
# cov_scen2_refit_trueb_farrellMinusConstant=resize(cov_scen2_refit_trueb_farrellMinusConstant,3,6,across="columns")
# cov_scen2_refit_sigNotCor_ours=resize(cov_scen2_refit_sigNotCor_ours,3,6,across="columns")

cov_scen2_refit_trueb_farrellMinusConstant <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
#
# cov_scen2_noselection_sigNotCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
# cov_scen2_noselection_sigCor_minna=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
#
cov_scen2_refit_sigNotCor_ours <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
# cov_scen2_refit_sigCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
#
# cov_scen2_doubleselection_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
# cov_scen2_doubleselection_sigNotCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
cov_scen2_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))


#CHANGES AFTER REVISION
bias_scen2_refit_trueb_farrellMinusConstant <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
bias_scen2_refit_sigNotCor_ours <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
bias_scen2_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))



sigma1 <- 1

library('tcltk')
set.seed(1127)
# Finding true values of the parameter####
pa <- c()
el <- c()
n <- 10000
param <- matrix(NA, nrow = length(alphavec), ncol = length(rhovec))
elambda <- c()
for (i in 1:length(alphavec)) {
  for (j in 1:length(rhovec)) {
    for (k in 1:rep) {
      rho <- rhovec[j]
      vecy <- vecy
      vect <- c(1, alphavec[i])
      sigma1 <- 1
      xvar <- 1
      l <- max(length(vect), length(vecy))
      X <- matrix(rnorm(n * l, mean = 0, sd = xvar), n, l)
      gX <- X %*% vect
      lam <- lambda1(gX)
      mX <- 2 + X %*% vecy - rho * sigma1 * lam
      pa[k] <- mean(mX + rnorm(n, 0, 1))
      if (j == 1) {
        el[k] <- mean(lam)
      }
    }
    param[i, j] <- mean(pa)
  }
  elambda[i] <- mean(el)
}

# main simulation
for (i in 1:(length(alphavec))) {
  for (j in 1:length(rhovec)) {
    simlow <- simulation_ui(
      FUN = generatingdata_ui, rep = rep, n = nvec[i], rho = rhovec[j],
      vecy = vecy, vect = c(1, alphavec[i]),
      trueValue = param[i, j], param = "Y1",
      subset = "refit", sigma_correction = 'non'
    )
    res <- results_ui(simlow$ui, param[i, j], true_b = rhovec[j] * sigma1 * elambda[i])
    cov_scen2_refit_sigNotCor_ours[i, j] <- res$cov[1]
    cov_scen2_refit_trueb_farrellMinusConstant[i, j] <- res$cov[2]
    
    #CHANGES AFTER REVISION
    bias_scen2_refit_trueb_farrellMinusConstant[i, j] <-  mean(res$tauhat_trueb) - param[i, j]
    bias_scen2_refit_sigNotCor_ours[i, j] <- mean(res$tauhat) - param[i, j]
    
    
    
    rm(simlow)
    gc()

    #
    #   simlow=simulation_ui(FUN=generatingdata_ui,rep=rep,n=nvec[i],rho=rhovec[j],
    #                        vecy=vecy,vect=c(1,alphavec[i]),
    #                        trueValue=param[i,j],param="Y1",
    #                        subset="refit",sigma_correction = 'old')
    #   res=results_ui(simlow$ui,param[i,j],true_b=rhovec[j]*sigma1*elambda[i])
    #   cov_scen2_refit_sigCor[i,j]=res$cov[1]
    #   rm(simlow)
    #   gc()
    #
    #
    #   simlow=simulation_ui(FUN=generatingdata_ui,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=vecy,vect=c(1,alphavec[i]),
    #                      trueValue=param[i,j],param="Y1",
    #                      subset="noselection",sigma_correction = 'old')
    # res=results_ui(simlow$ui,param[i,j],true_b=rhovec[j]*sigma1*elambda[i])
    # cov_scen2_noselection_sigCor_minna[i,j]=res$cov[1]
    # rm(simlow)
    # gc()
    #
    # simlow=simulation_ui(FUN=generatingdata_ui,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=vecy,vect=c(1,alphavec[i]),
    #                      trueValue=param[i,j],param="Y1",
    #                      subset="noselection",sigma_correction = 'non')
    # res=results_ui(simlow$ui,param[i,j],true_b=rhovec[j]*sigma1*elambda[i])
    # cov_scen2_noselection_sigNotCor[i,j]=res$cov[1]
    # rm(simlow)
    # gc()

#
#     simlow <- simulation_ui(
#       FUN = generatingdata_ui, rep = rep, n = nvec[i], rho = rhovec[j],
#       vecy = vecy, vect = c(1, alphavec[i]),
#       trueValue = param[i, j], param = "Y1",
#       subset = "double", sigma_correction = 'old'
#     )
#     res <- results_ui(simlow$ui, param[i, j], true_b = rhovec[j] * sigma1 * elambda[i])
#     cov_scen2_doubleselection_sigCor_notheory[i, j] <- res$cov[1]
#     rm(simlow)
#     gc()





    simlow <- simulation_ui(
      FUN = generatingdata_ui, rep = rep, n = nvec[i], rho = rhovec[j],
      vecy = vecy, vect = c(1, alphavec[i]),
      trueValue = param[i, j], param = "Y1",
      subset = "refit", sigma_correction = 'new'
    )
    res <- results_ui(simlow$ui, param[i, j], true_b = rhovec[j] * sigma1 * elambda[i])
    cov_scen2_new_sigCor_notheory[i, j] <- res$cov[1]
    
    #CHANGES AFTER REVISION
    bias_scen2_new_sigCor_notheory[i, j] <- mean(res$tauhat) -  param[i, j]
    
    
    rm(simlow)
    gc()

    # simlow=simulation_ui(FUN=generatingdata_ui,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=vecy,vect=c(1,alphavec[i]),
    #                      trueValue=param[i,j],param="Y1",
    #                      subset="double",sigma_correction = 'non')
    # res=results_ui(simlow$ui,param[i,j],true_b=rhovec[j]*sigma1*elambda[i])
    # cov_scen2_doubleselection_sigNotCor[i,j]=res$cov[1]
    # rm(simlow)
    # gc()
  }
}



#CHANGES AFTER REVISION
save.image(".../scenarioTwoRevision.Rdata")























# save.image("scenarioTwo.Rdata")
#
# # load("scenarioTwo.Rdata")
#CHANGED AFTER REVISION
#save.image("scenarioTwo02.Rdata")
# load("scenarioTwo02.Rdata")

library(xtable)

temp <- bias_scen2_refit_trueb_farrellMinusConstant
rownames(temp) <- nvec
colnames(temp) <- rhovec
xtable(temp)

temp <- cov_scen2_noselection_sigNotCor
temp <- cov_scen2_noselection_sigCor_minna

temp <- cov_scen2_refit_sigNotCor_ours
temp <- cov_scen2_refit_sigCor

temp <- cov_scen2_doubleselection_sigNotCor
temp <- cov_scen2_doubleselection_sigCor_notheory
