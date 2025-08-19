# load("scenarioFour_final.Rdata")
library(devtools)
# install("C:\\Users\\...\\validui\\code\\packages\\ui.causal")
# document("C:\\Users\\...\\validui\\code\\packages\\ui.causal")
# install("C:\\Users\\...\\validui\\code\\packages\\ui")
# document("C:\\Users\\...\\validui\\code\\packages\\ui")



install_github('https://github.com/stat4reg/hdim.ui.git')
library(hdim.ui)

#CHANGES AFTER REVISION
#setwd()
#source("validui_func.R")
source("/Users/niloofarmoosavi/Library/CloudStorage/OneDrive-Personal/PhD/validui/code/final/Paper Two-final/validui_func.R")



library(glmnet)
library(MASS)
library(tcltk)
library(hdm)
library("mgcv")



# First Scenario #####
simlow <- list()
nvec <- c(500, 1000, 1500)

#CHANGED AFTER REVISION
rhovec <- seq(0.8, 0.3, by = -0.1)
#rhovec=0.2
rep <- 500




## lowdim#########################


#####
cov_scen4_refit_trueb_farrellMinusConstant <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
#
# cov_scen4_noselection_sigNotCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
# cov_scen4_noselection_sigCor_minna=matrix(NA,nrow = length(nvec),ncol = length(rhovec))

cov_scen4_refit_sigNotCor_ours <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
# cov_scen4_refit_sigCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))

# cov_scen4_doubleselection_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
# cov_scen4_doubleselection_sigNotCor=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
cov_scen4_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))


#CHANGES AFTER REVISION
bias_scen4_refit_trueb_farrellMinusConstant <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
bias_scen4_refit_sigNotCor_ours <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
bias_scen4_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))




sigma1 <- 1

library('tcltk')
set.seed(1127)
# Finding true values of the parameter####
pa <- c()
el <- c()
n <- 10000
# param <- matrix(NA, nrow = length(alphavec), ncol = length(rhovec))
param <- c()
elambda <- 0
# for (i in 1:length(alphavec)) {
for (j in 1:length(rhovec)) {
  for (k in 1:rep) {
    rho <- rhovec[j]
    vecy <- 0.6
    vect <- 0.3


    p <- 10
    sigma1 <- 1
    xvar <- 1
    X <- matrix(rnorm(n * p, mean = 0, sd = xvar), n, p)

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
    lam <- lambda1(gX)
    mX <-2+ X %*% b * vecy- rho*sigma1*lam



    pa[k] <- mean( mX )
    if (j == 1) {
      el[k] <- mean(lam)
    }
  }
  param[j] <- mean(pa)
  # }
}
elambda <- mean(el)

# main simulation
for (i in 1:(length(nvec))) {
  for (j in 1:length(rhovec)) {
    simlow <- simulation_ui(
      FUN = generatingdata_ui_high, rep = rep, n = nvec[i], rho = rhovec[j],
      vecy = 0.6, vect = 0.3,
      trueValue = param[j], param = "Y1",
      subset = "refit", sigma_correction = 'non'
    )
    res <- results_ui(simlow$ui, param[j],
      true_b = rhovec[j] * sigma1 * elambda, var_detect_percent = c(1:10)
    )
    cov_scen4_refit_sigNotCor_ours[i, j] <- res$cov[1]
    cov_scen4_refit_trueb_farrellMinusConstant[i, j] <- res$cov[2]
    
    #CHANGES AFTER REVISION
    bias_scen4_refit_trueb_farrellMinusConstant[i, j] <-  mean(res$tauhat_trueb) - param[j]
    bias_scen4_refit_sigNotCor_ours[i, j] <- mean(res$tauhat) - param[j]
    
    
    rm(simlow)
    gc()


    # simlow=simulation_ui(FUN=generatingdata_ui,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=c(-0.5,alphavec[i]),vect=vect,
    #                      trueValue=param[j], param="Y1",
    #                      subset="refit",sigma_correction='old')
    # res=results_ui(simlow$ui,param[j],true_b=rhovec[j]*sigma1*elambda)
    # cov_scen4_refit_sigCor[i,j]=res$cov[1]
    # rm(simlow)
    # gc()
    #
    #
    # simlow=simulation_ui(FUN=generatingdata_ui_high,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=0.5,vect=0.3,
    #                      trueValue=param[j], param="Y1",
    #                      subset="noselection", sigma_correction='old')
    # res=results_ui(simlow$ui,param[j],true_b=rhovec[j]*sigma1*elambda)
    # cov_scen4_noselection_sigCor_minna[i,j]=res$cov[1]
    # rm(simlow)
    # gc()
    #
    #
    #
    # simlow=simulation_ui(FUN=generatingdata_ui,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=c(-0.5,alphavec[i]),vect=vect,
    #                      trueValue=param[j], param="Y1",
    #                      subset="noselection", sigma_correction='non')
    # res=results_ui(simlow$ui,param[j],true_b=rhovec[j]*sigma1*elambda)
    # cov_scen4_noselection_sigNotCor[i,j]=res$cov[1]
    # rm(simlow)
    # gc()


    # simlow <- simulation_ui(
    #   FUN = generatingdata_ui_high, rep = rep, n = nvec[i], rho = rhovec[j],
    #   vecy = 0.6, vect = 0.3,
    #   trueValue = param[j], param = "Y1",
    #   subset = "double", sigma_correction='old'
    # )
    # res <- results_ui(simlow$ui, param[j], true_b = rhovec[j] * sigma1 * elambda)
    # cov_scen4_doubleselection_sigCor_notheory[i, j] <- res$cov[1]
    # rm(simlow)
    # gc()


    simlow <- simulation_ui(
      FUN = generatingdata_ui_high, rep = rep, n = nvec[i], rho = rhovec[j],
      vecy = 0.6, vect = 0.3,
      trueValue = param[j], param = "Y1",
      subset = "refit", sigma_correction = 'new'
    )
    res <- results_ui(simlow$ui, param[j],
                      true_b = rhovec[j] * sigma1 * elambda, var_detect_percent = c(1:10)
    )
    cov_scen4_new_sigCor_notheory[i, j] <- res$cov[1]
    
    #CHANGES AFTER REVISION
    bias_scen4_new_sigCor_notheory[i, j] <- mean(res$tauhat) -  param[j]
    
    
    rm(simlow)
    gc()

    # simlow=simulation_ui(FUN=generatingdata_ui,rep=rep,n=nvec[i],rho=rhovec[j],
    #                      vecy=c(-0.5,alphavec[i]),vect=vect,
    #                      trueValue=param[j], param="Y1",
    #                      subset="double",sigma_correction='non')
    # res=results_ui(simlow$ui,param[j],true_b=rhovec[j]*sigma1*elambda)
    # cov_scen4_doubleselection_sigNotCor[i,j]=res$cov[1]
    # rm(simlow)
    # gc()
  }
}



#CHANGES AFTER REVISION
save.image("/Users/.../scenarioFourRevision.Rdata")

















# # save.image("scenarioFour_final.Rdata")
#
# load("scenarioFour_final.Rdata")

#CHANGED AFTER REVISION
#save.image("scenarioFour_final02.Rdata")
# load("scenarioFour_final02.Rdata")



library(xtable)

temp <- cov_scen4_refit_trueb_farrellMinusConstant[1:3,c(1,3,5)]
rownames(temp) <- nvec
colnames(temp) <- rhovec[c(1,3,5)]
xtable(temp)

# temp=cov_scen4_noselection_sigNotCor
# temp=cov_scen4_noselection_sigCor_minna

temp <- cov_scen4_refit_sigNotCor_ours
# temp=cov_scen4_refit_sigCor

# temp=cov_scen4_doubleselection_sigNotCor
temp <- cov_scen4_new_sigCor_notheory
