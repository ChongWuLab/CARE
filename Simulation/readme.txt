
#simulation_dir_pleio.R: simRes1 (the main setting in manuscript)

 ind2 = ind2all[1:floor(length(ind2) * 0.5)]
        phi[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2u))
        
        alpha[ind2] = runif(length(ind2),0.01,0.03)

        ind2 = ind2all[!ind2all %in% ind2]
        alpha[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2y_td))
        



#simulation_dir_pleio2.R: simRes2 (supplementary)

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




