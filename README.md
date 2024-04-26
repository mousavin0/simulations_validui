# simulations_validui

Scenario1-3 are low dimensional cases in appendix, scenario4 is the high dimensional case.
Rdata files have the results. There are three matrices for coverages of each estimators and three for biases:



cov_scen3_refit_trueb_farrellMinusConstant

cov_scen3_refit_sigNotCor_ours

cov_scen3_new_sigCor_notheory



bias_scen3_refit_trueb_farrellMinusConstant

bias_scen3_refit_sigNotCor_ours 

bias_scen3_new_sigCor_notheory 



those with _refit_trueb_farrellMinusConstant are related to the estimator tauhat - b* in the paper,
_scen3_refit_sigNotCor_ours is for tauhat - bhatrefit,
_scen3_new_sigCor_notheory is for tahhat - bhatrefit_c.



I added two generating data functions for the case where error for treatment is logistic. I will add four rcode to get the result of simulations for those dgps.
