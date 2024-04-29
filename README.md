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


# results


## Normal distribution error term for T


### High dimensional scenario

***Bstar***

% latex table generated in R 4.3.3 by xtable 1.8-4 package
% Mon Apr 29 07:50:44 2024
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrr}
  \hline
 & 0.8 & 0.7 & 0.6 & 0.5 & 0.4 & 0.3 \\ 
  \hline
500 & -0.01 & 0.04 & 0.03 & 0.02 & 0.03 & 0.02 \\ 
  1000 & 0.01 & 0.01 & 0.02 & -0.00 & 0.02 & 0.01 \\ 
  1500 & 0.00 & 0.01 & 0.01 & 0.01 & 0.01 & 0.01 \\ 
   \hline
\end{tabular}
\end{table}


***BHat***

% latex table generated in R 4.3.3 by xtable 1.8-4 package
% Mon Apr 29 07:50:44 2024
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrr}
  \hline
 & 0.8 & 0.7 & 0.6 & 0.5 & 0.4 & 0.3 \\ 
  \hline
500 & 0.03 & 0.04 & 0.02 & -0.01 & -0.00 & -0.01 \\ 
  1000 & 0.09 & 0.05 & 0.03 & -0.01 & 0.00 & -0.00 \\ 
  1500 & 0.10 & 0.06 & 0.03 & 0.01 & 0.01 & -0.00 \\ 
   \hline
\end{tabular}
\end{table}


***BHatCorrected***
% latex table generated in R 4.3.3 by xtable 1.8-4 package
% Mon Apr 29 07:50:44 2024
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrr}
  \hline
 & 0.8 & 0.7 & 0.6 & 0.5 & 0.4 & 0.3 \\ 
  \hline
500 & -0.07 & -0.05 & -0.04 & -0.05 & -0.03 & 0.00 \\ 
  1000 & -0.05 & -0.03 & -0.03 & -0.02 & -0.01 & -0.01 \\ 
  1500 & -0.02 & -0.03 & -0.02 & -0.02 & -0.01 & -0.01 \\ 
   \hline
\end{tabular}
\end{table}
