rm(list= ls())
#library(ggplot2)

summary_fun <- function(theta,N,prop_invalid,simSetting) {
    #theta = 0.03
    #propInvalid = 0.5
    #simSetting = "simRes2"
    setwd(paste0("/rsrch5/home/biostatistics/chongwulab/wzhang24/CARE/simulation_addition/",simSetting))
    #setwd(paste0("C:\\Users\\Evelyn\\OneDrive - The University of Texas Health Science Center at Houston\\Wu group\\MR\\codes\\simulation-final\\",simSetting))
    if (N == 5e+05) {
      ind1 = "N5e\\+05"
    } else {
      ind1 = paste0("N",N)
    }

    ind2 = paste0("theta",theta,"thetaU0.3prop_invalid",prop_invalid)

    care_result = list()
    care2_result = list()
    care3_result = list()


    files = list.files(getwd())
    files = files[grepl(ind1,files)]
    files = files[grepl(ind2,files)]
    
    # only using sample size = 1e5
    #files = files[grepl("N1",files)]


    for(set.ind in 1:length(files))
    {
      #load(paste("sim_results/simres_MRMethods",set.ind, "N50000pthr5e-08pi10.01theta",theta,"thetaU0.3prop_invalid0.5NxNy_ratio2_July22.Rdata",sep=""))
      load(files[set.ind])
      
      #run_time = c(run_time,runningtime)
      #settings = c(settings, setting)
      care_result = c(care_result,care_sim_result) # check
      care2_result = c(care2_result,care2_sim_result)
      care3_result = c(care3_result,care3_sim_result)
    }
    
    nsim = length(care_result)
    nsim

    care.res = matrix(NA,nsim,40)
    for(i in 1:nsim)
    {
       tmp = care_result[[i]]
       
       tmp = unlist(tmp)
       care.res[i,] = tmp
    }

    colnames(care.res) = names(tmp)
    
    
    nsim = length(care2_result)
    nsim

    care2.res = matrix(NA,nsim,40)
    for(i in 1:nsim)
    {
       tmp = care2_result[[i]]
       
       tmp = unlist(tmp)
       care2.res[i,] = tmp
    }

    colnames(care2.res) = names(tmp)

    nsim = length(care3_result)
    nsim
    care3.res = matrix(NA,nsim,40)
    for(i in 1:nsim)
    {
       tmp = care3_result[[i]]
       
       tmp = unlist(tmp)
       care3.res[i,] = tmp
    }

    colnames(care3.res) = names(tmp)
 
    cutoff = 0.05
 
    final.out = as.data.frame(matrix(NA,24,7))
    rownames(final.out) = c("CARE_MA_AIC","CARE_AIC","CARE_MA_BIC","CARE_BIC","CARE_Efron_AIC","CARE_Efron_BIC","CARE_Efron_MA_AIC","CARE_Efron_MA_BIC","CARE2_MA_AIC","CARE2_AIC","CARE2_MA_BIC","CARE2_BIC","CARE2_Efron_AIC","CARE2_Efron_BIC","CARE2_Efron_MA_AIC","CARE2_Efron_MA_BIC","CARE3_MA_AIC","CARE3_AIC","CARE3_MA_BIC","CARE3_BIC","CARE3_Efron_AIC","CARE3_Efron_BIC","CARE3_Efron_MA_AIC","CARE3_Efron_MA_BIC")
    colnames(final.out) = c("Power","Estimate","Coverage","Bias","MSE","Monte_SD","SD")
     # power
    final.out[1:8,1] = colSums(care.res[,c("MA_AIC.p","AIC.p","MA_BIC.p","BIC.p","AIC.Efron_p","BIC.Efron_p","MA_AIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care.res)[1]
    final.out[9:16,1] = colSums(care2.res[,c("MA_AIC.p","AIC.p","MA_BIC.p","BIC.p","AIC.Efron_p","BIC.Efron_p","MA_AIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care2.res)[1]
    final.out[17:24,1] = colSums(care3.res[,c("MA_AIC.p","AIC.p","MA_BIC.p","BIC.p","AIC.Efron_p","BIC.Efron_p","MA_AIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care3.res)[1]


    # Estimate
    final.out[1:8,2] = colMeans(care.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])
    final.out[9:16,2] = colMeans(care2.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])
    final.out[17:24,2] = colMeans(care3.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])

    
    # Coverage
    final.out[1,3] = sum(theta > care.res[,"MA_AIC.tilde_theta"] - 1.96 * care.res[,"MA_AIC.se"] & theta < care.res[,"MA_AIC.tilde_theta"] + 1.96 * care.res[,"MA_AIC.se"]) / dim(care.res)[1]
    final.out[2,3] = sum(theta > care.res[,"AIC.tilde_theta"] - 1.96 * care.res[,"AIC.se"] & theta < care.res[,"AIC.tilde_theta"] + 1.96 * care.res[,"AIC.se"]) / dim(care.res)[1]
    final.out[3,3] = sum(theta > care.res[,"MA_BIC.tilde_theta"] - 1.96 * care.res[,"MA_BIC.se"] & theta < care.res[,"MA_BIC.tilde_theta"] + 1.96 * care.res[,"MA_BIC.se"]) / dim(care.res)[1]
    final.out[4,3] = sum(theta > care.res[,"BIC.tilde_theta"] - 1.96 * care.res[,"BIC.se"] & theta < care.res[,"BIC.tilde_theta"] + 1.96 * care.res[,"BIC.se"]) / dim(care.res)[1]
    final.out[5,3] = sum(theta > care.res[,"AIC.tilde_theta"] - 1.96 * care.res[,"AIC.Efron_se"] & theta < care.res[,"AIC.tilde_theta"] + 1.96 * care.res[,"AIC.Efron_se"]) / dim(care.res)[1]
    final.out[6,3] = sum(theta > care.res[,"BIC.tilde_theta"] - 1.96 * care.res[,"BIC.Efron_se"] & theta < care.res[,"BIC.tilde_theta"] + 1.96 * care.res[,"BIC.Efron_se"]) / dim(care.res)[1]
    final.out[7,3] = sum(theta > care.res[,"MA_AIC.tilde_theta"] - 1.96 * care.res[,"MA_AIC.Efron_se"] & theta < care.res[,"MA_AIC.tilde_theta"] + 1.96 * care.res[,"MA_AIC.Efron_se"]) / dim(care.res)[1]
    final.out[8,3] = sum(theta > care.res[,"MA_BIC.tilde_theta"] - 1.96 * care.res[,"MA_BIC.Efron_se"] & theta < care.res[,"MA_BIC.tilde_theta"] + 1.96 * care.res[,"MA_BIC.Efron_se"]) / dim(care.res)[1]

    final.out[9,3] = sum(theta > care2.res[,"MA_AIC.tilde_theta"] - 1.96 * care2.res[,"MA_AIC.se"] & theta < care2.res[,"MA_AIC.tilde_theta"] + 1.96 * care2.res[,"MA_AIC.se"]) / dim(care2.res)[1]
    final.out[10,3] = sum(theta > care2.res[,"AIC.tilde_theta"] - 1.96 * care2.res[,"AIC.se"] & theta < care2.res[,"AIC.tilde_theta"] + 1.96 * care2.res[,"AIC.se"]) / dim(care2.res)[1]
    final.out[11,3] = sum(theta > care2.res[,"MA_BIC.tilde_theta"] - 1.96 * care2.res[,"MA_BIC.se"] & theta < care2.res[,"MA_BIC.tilde_theta"] + 1.96 * care2.res[,"MA_BIC.se"]) / dim(care2.res)[1]
    final.out[12,3] = sum(theta > care2.res[,"BIC.tilde_theta"] - 1.96 * care2.res[,"BIC.se"] & theta < care2.res[,"BIC.tilde_theta"] + 1.96 * care2.res[,"BIC.se"]) / dim(care2.res)[1]
    final.out[13,3] = sum(theta > care2.res[,"AIC.tilde_theta"] - 1.96 * care2.res[,"AIC.Efron_se"] & theta < care2.res[,"AIC.tilde_theta"] + 1.96 * care2.res[,"AIC.Efron_se"]) / dim(care2.res)[1]
    final.out[14,3] = sum(theta > care2.res[,"BIC.tilde_theta"] - 1.96 * care2.res[,"BIC.Efron_se"] & theta < care2.res[,"BIC.tilde_theta"] + 1.96 * care2.res[,"BIC.Efron_se"]) / dim(care2.res)[1]
    final.out[15,3] = sum(theta > care2.res[,"MA_AIC.tilde_theta"] - 1.96 * care2.res[,"MA_AIC.Efron_se"] & theta < care2.res[,"MA_AIC.tilde_theta"] + 1.96 * care2.res[,"MA_AIC.Efron_se"]) / dim(care2.res)[1]
    final.out[16,3] = sum(theta > care2.res[,"MA_BIC.tilde_theta"] - 1.96 * care2.res[,"MA_BIC.Efron_se"] & theta < care2.res[,"MA_BIC.tilde_theta"] + 1.96 * care2.res[,"MA_BIC.Efron_se"]) / dim(care2.res)[1]

    final.out[17,3] = sum(theta > care3.res[,"MA_AIC.tilde_theta"] - 1.96 * care3.res[,"MA_AIC.se"] & theta < care3.res[,"MA_AIC.tilde_theta"] + 1.96 * care3.res[,"MA_AIC.se"]) / dim(care3.res)[1]
    final.out[18,3] = sum(theta > care3.res[,"AIC.tilde_theta"] - 1.96 * care3.res[,"AIC.se"] & theta < care3.res[,"AIC.tilde_theta"] + 1.96 * care3.res[,"AIC.se"]) / dim(care3.res)[1]
    final.out[19,3] = sum(theta > care3.res[,"MA_BIC.tilde_theta"] - 1.96 * care3.res[,"MA_BIC.se"] & theta < care3.res[,"MA_BIC.tilde_theta"] + 1.96 * care3.res[,"MA_BIC.se"]) / dim(care3.res)[1]
    final.out[20,3] = sum(theta > care3.res[,"BIC.tilde_theta"] - 1.96 * care3.res[,"BIC.se"] & theta < care3.res[,"BIC.tilde_theta"] + 1.96 * care3.res[,"BIC.se"]) / dim(care3.res)[1]
    final.out[21,3] = sum(theta > care3.res[,"AIC.tilde_theta"] - 1.96 * care3.res[,"AIC.Efron_se"] & theta < care3.res[,"AIC.tilde_theta"] + 1.96 * care3.res[,"AIC.Efron_se"]) / dim(care3.res)[1]
    final.out[22,3] = sum(theta > care3.res[,"BIC.tilde_theta"] - 1.96 * care3.res[,"BIC.Efron_se"] & theta < care3.res[,"BIC.tilde_theta"] + 1.96 * care3.res[,"BIC.Efron_se"]) / dim(care3.res)[1]
    final.out[23,3] = sum(theta > care3.res[,"MA_AIC.tilde_theta"] - 1.96 * care3.res[,"MA_AIC.Efron_se"] & theta < care3.res[,"MA_AIC.tilde_theta"] + 1.96 * care3.res[,"MA_AIC.Efron_se"]) / dim(care3.res)[1]
    final.out[24,3] = sum(theta > care3.res[,"MA_BIC.tilde_theta"] - 1.96 * care3.res[,"MA_BIC.Efron_se"] & theta < care3.res[,"MA_BIC.tilde_theta"] + 1.96 * care3.res[,"MA_BIC.Efron_se"]) / dim(care3.res)[1]


    # bias
    final.out[1:8,4] = (theta -  colMeans(care.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])) #/ theta   
    final.out[9:16,4] = (theta -  colMeans(care2.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])) #/ theta
    final.out[17:24,4] = (theta -  colMeans(care3.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])) #/ theta
   

    # MSE
    final.out[1:8,5] = colMeans((care.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")] -theta )^2)
    final.out[9:16,5] = colMeans((care2.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")] -theta )^2)
    final.out[17:24,5] = colMeans((care3.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")] -theta )^2)
  
    # Monte SD
    final.out[1:8,6] = apply(care.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")],2,sd)
    final.out[9:16,6] = apply(care2.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")],2,sd)
    final.out[17:24,6] = apply(care3.res[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")],2,sd)
 
    # SD
    final.out[1:8,7] = colMeans(care.res[,c("MA_AIC.se","AIC.se","MA_BIC.se","BIC.se","AIC.Efron_se","BIC.Efron_se","MA_AIC.Efron_se","MA_BIC.Efron_se")])
    final.out[9:16,7] = colMeans(care2.res[,c("MA_AIC.se","AIC.se","MA_BIC.se","BIC.se","AIC.Efron_se","BIC.Efron_se","MA_AIC.Efron_se","MA_BIC.Efron_se")])
    final.out[17:24,7] = colMeans(care3.res[,c("MA_AIC.se","AIC.se","MA_BIC.se","BIC.se","AIC.Efron_se","BIC.Efron_se","MA_AIC.Efron_se","MA_BIC.Efron_se")])
  
    savefile = paste0("/rsrch5/home/biostatistics/chongwulab/wzhang24/CARE/simulation_addition/summaryRes/",simSetting,"theta",theta,"N", N, "prop_invalid",prop_invalid,".RData")
    save(final.out, file = savefile)
}





for (simSetting in c("simRes15")){
  for(theta in c(0.1, 0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.3,0.7,0.9)) {
        for (N in c(5e+05)) {
          try(summary_fun(theta,N, prop_invalid,simSetting))    
          cat("Finish theta", theta, "  Prop",prop_invalid,"\n")
        }
    }   
  }
}

for (simSetting in c("simRes20")){
  for(theta in c(0.1, 0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.5)) {
        for (N in c(5e+05)) {
          try(summary_fun(theta,N, prop_invalid,simSetting))    
          cat("Finish theta", theta, "  Prop",prop_invalid,"\n")
        }
    }   
  }
}