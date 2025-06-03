
#simulation_dir_pleio.R: simRes1 (supplementary)

 ind2 = ind2all[1:floor(length(ind2) * 0.5)]
        phi[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2u))
        
        alpha[ind2] = runif(length(ind2),0.01,0.03)

        ind2 = ind2all[!ind2all %in% ind2]
        alpha[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2y_td))
        



#simulation_dir_pleio2.R: simRes2 (the main stting in manuscript)

       ind2 = ind2all[1:floor(length(ind2) * 0.5)]
        phi[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2u))
        
        alpha[ind2] = rnorm(length(ind2), mean = 0.015, sd = sqrt(sigma2u))
 


#simulation_balance_pleio.R: simRes3 (supplementary)
balance pleiotropy with InSIDE assumption satisfies

         ind2all = ind2
        
        #ind2 = ind2all[1:floor(length(ind2) * 0.5)]
        #phi[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2u))
        #alpha[ind2] = runif(length(ind2),0.01,0.03)
        #ind2 = ind2all[!ind2all %in% ind2]
        
        alpha[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2y_td))
        

#simulation_balance_pleio2.R: simRes4 (supplementary)
balance pleiotropy with InSIDE assumption violates

        ind2all = ind2
                
        #ind2 = ind2all[1:floor(length(ind2) * 0.5)]
        phi[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2u))
        alpha[ind2] = runif(length(ind2),-0.03,0.03)     


#simmulation_nonlinear.R: simRes5 (supplementary) N = 10000
#simRes6: nonlinear with interaction term. N = 10000
#simRes7: small mean of normal distribution of alphaj for correlated pleiotropy. 
#simRes8: GIC simulation_dir_pleio3.R
#simRes9: different sample size. simulation_dir_pleio4.R
#simRes10: third sample under main setting. p_thr = 5e-5 for all methods. simulation_dir_pleio5.R
#simRes11: third sample under main setting. And p_thr for CARE is 5e-5. simulation_dir_pleio6.R
#simRes12: different eta. Main setting. simulation_dir_pleio7.R
#simRes13: different sample size of SNPs. simulation_dir_pleio8.R
#simRes14: nonlinear with both nonlinear y and nonlinear x. simmulation_nonlinear3.R
#simRes16: 5e-5 for other methods. simulation_dir_pleio10.R
#simRes17: additional simulation setup 1. simulation_dir_pleio11.R
#simRes18: additional simulation setup 2. simulation_dir_pleio12.R
#simRes19: additional simulation setup 3. simulation_dir_pleio13.R