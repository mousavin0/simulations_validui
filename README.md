# simulations_validui

Scenario1-3 are low dimensional cases in appendix, scenario4 is the high dimensional case.
Rdata files have the results. There are three matrices for coverages of each estimators and three for biases:



cov_scen3_refit_trueb_farrellMinusConstant=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
cov_scen3_refit_sigNotCor_ours=matrix(NA,nrow = length(nvec),ncol = length(rhovec))
cov_scen3_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))


#CHANGES AFTER REVISION
bias_scen3_refit_trueb_farrellMinusConstant <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
bias_scen3_refit_sigNotCor_ours <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))
bias_scen3_new_sigCor_notheory <- matrix(NA, nrow = length(nvec), ncol = length(rhovec))



those with _refit_trueb_farrellMinusConstant are related to the estimator tauhat - b* in the paper,
_scen3_refit_sigNotCor_ours is for tauhat - bhatrefit,
_scen3_new_sigCor_notheory is for tahhat - bhatrefit_c.



I added two generating data functions for the case where error for treatment is logistic. I will add four rcode to get the result of simulations for those dgps.
