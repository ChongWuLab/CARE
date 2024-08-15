rm(list= ls())
#library(ggplot2)

summary_fun <- function(theta,N,prop_invalid,simSetting) {
    #theta = 0.03
    #propInvalid = 0.5
    #simSetting = "simRes2"
    setwd(paste0("/rsrch5/home/biostatistics/chongwulab/wzhang24/CARE/simulation_final_2/",simSetting))
    #setwd(paste0("C:\\Users\\Evelyn\\OneDrive - The University of Texas Health Science Center at Houston\\Wu group\\MR\\codes\\simulation-final\\",simSetting))
    if (N == 5e+05) {
      ind1 = "N5e\\+05"
    } else {
      ind1 = paste0("N",N)
    }

    ind2 = paste0("theta",theta,"thetaU0.3prop_invalid",prop_invalid)

    care_result_0.1 = list()
    care_result_0.3 = list()
    care_result_0.5 = list()
    care_result_0.7 = list()
    care_result_0.9 = list()
    #care_no_result = list()

    #ConMix_result = list()
    #MR_Lasso_result = list()
    #MRMix_result = list()
    #MR_PRESSO_result = list()
    #MRcML_result = list()
    #TSMR_result = NULL
    #MRAPSS_result = list()
    #run_time = list()
    #settings = list()

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
      care_result_0.1 = c(care_result_0.1,care_sim_result_0.1) # check
      care_result_0.3 = c(care_result_0.3,care_sim_result_0.3)
      care_result_0.5 = c(care_result_0.5,care_sim_result_0.5)
      care_result_0.7 = c(care_result_0.7,care_sim_result_0.7)
      care_result_0.9 = c(care_result_0.9,care_sim_result_0.9)
      
    }
    
    ## CARE with eta =0.5, eta_var = 0
    nsim = length(care_result_0.1)
    nsim

    care.res0.1 = matrix(NA,nsim,40)
    for(i in 1:nsim)
    {
       tmp = care_result_0.1[[i]]
       
       tmp = unlist(tmp)
       care.res0.1[i,] = tmp
    }

    colnames(care.res0.1) = names(tmp)
    
    
    nsim = length(care_result_0.3)
    nsim

    care.res0.3 = matrix(NA,nsim,40)
    for(i in 1:nsim)
    {
       tmp = care_result_0.3[[i]]
       
       tmp = unlist(tmp)
       care.res0.3[i,] = tmp
    }

    colnames(care.res0.3) = names(tmp)

    nsim = length(care_result_0.5)
    nsim
    care.res0.5 = matrix(NA,nsim,40)
    for(i in 1:nsim)
    {
       tmp = care_result_0.5[[i]]
       
       tmp = unlist(tmp)
       care.res0.5[i,] = tmp
    }
    colnames(care.res0.5) = names(tmp)

    nsim = length(care_result_0.7)
    nsim
    care.res0.7 = matrix(NA,nsim,40)
    for(i in 1:nsim)
    {
       tmp = care_result_0.7[[i]]
       
       tmp = unlist(tmp)
       care.res0.7[i,] = tmp
    }
    colnames(care.res0.7) = names(tmp)

    nsim = length(care_result_0.9)
    nsim
    care.res0.9 = matrix(NA,nsim,40)
    for(i in 1:nsim)
    {
       tmp = care_result_0.9[[i]]
       
       tmp = unlist(tmp)
       care.res0.9[i,] = tmp
    }
    colnames(care.res0.9) = names(tmp)
 
    cutoff = 0.05
 
    final.out = as.data.frame(matrix(NA,40,7))
    rownames(final.out) = c("CARE_eta0.1_MA_AIC","CARE_eta0.1_AIC","CARE_eta0.1_MA_BIC","CARE_eta0.1_BIC","CARE_eta0.1_Efron_AIC","CARE_eta0.1_Efron_BIC","CARE_eta0.1_Efron_MA_AIC","CARE_eta0.1_Efron_MA_BIC",
                            "CARE_eta0.3_MA_AIC","CARE_eta0.3_AIC","CARE_eta0.3_MA_BIC","CARE_eta0.3_BIC","CARE_eta0.3_Efron_AIC","CARE_eta0.3_Efron_BIC","CARE_eta0.3_Efron_MA_AIC","CARE_eta0.3_Efron_MA_BIC",
                            "CARE_eta0.5_MA_AIC","CARE_eta0.5_AIC","CARE_eta0.5_MA_BIC","CARE_eta0.5_BIC","CARE_eta0.5_Efron_AIC","CARE_eta0.5_Efron_BIC","CARE_eta0.5_Efron_MA_AIC","CARE_eta0.5_Efron_MA_BIC",
                            "CARE_eta0.7_MA_AIC","CARE_eta0.7_AIC","CARE_eta0.7_MA_BIC","CARE_eta0.7_BIC","CARE_eta0.7_Efron_AIC","CARE_eta0.7_Efron_BIC","CARE_eta0.7_Efron_MA_AIC","CARE_eta0.7_Efron_MA_BIC",
                            "CARE_eta0.9_MA_AIC","CARE_eta0.9_AIC","CARE_eta0.9_MA_BIC","CARE_eta0.9_BIC","CARE_eta0.9_Efron_AIC","CARE_eta0.9_Efron_BIC","CARE_eta0.9_Efron_MA_AIC","CARE_eta0.9_Efron_MA_BIC")
    colnames(final.out) = c("Power","Estimate","Coverage","Bias","MSE","Monte_SD","SD")
     # power
    final.out[1:8,1] = colSums(care.res0.1[,c("MA_AIC.p","AIC.p","MA_BIC.p","BIC.p","AIC.Efron_p","BIC.Efron_p","MA_AIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care.res0.1)[1]
    final.out[9:16,1] = colSums(care.res0.3[,c("MA_AIC.p","AIC.p","MA_BIC.p","BIC.p","AIC.Efron_p","BIC.Efron_p","MA_AIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care.res0.3)[1]
    final.out[17:24,1] = colSums(care.res0.5[,c("MA_AIC.p","AIC.p","MA_BIC.p","BIC.p","AIC.Efron_p","BIC.Efron_p","MA_AIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care.res0.5)[1]
    final.out[25:32,1] = colSums(care.res0.7[,c("MA_AIC.p","AIC.p","MA_BIC.p","BIC.p","AIC.Efron_p","BIC.Efron_p","MA_AIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care.res0.7)[1]
    final.out[33:40,1] = colSums(care.res0.9[,c("MA_AIC.p","AIC.p","MA_BIC.p","BIC.p","AIC.Efron_p","BIC.Efron_p","MA_AIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care.res0.9)[1]

    # Estimate
    final.out[1:8,2] = colMeans(care.res0.1[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])
    final.out[9:16,2] = colMeans(care.res0.3[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])
    final.out[17:24,2] = colMeans(care.res0.5[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])
    final.out[25:32,2] = colMeans(care.res0.7[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])
    final.out[33:40,2] = colMeans(care.res0.9[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])
    
    
    # Coverage
    final.out[1,3] = sum(theta > care.res0.1[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.1[,"MA_AIC.se"] & theta < care.res0.1[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.1[,"MA_AIC.se"]) / dim(care.res0.1)[1]
    final.out[2,3] = sum(theta > care.res0.1[,"AIC.tilde_theta"] - 1.96 * care.res0.1[,"AIC.se"] & theta < care.res0.1[,"AIC.tilde_theta"] + 1.96 * care.res0.1[,"AIC.se"]) / dim(care.res0.1)[1]
    final.out[3,3] = sum(theta > care.res0.1[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.1[,"MA_BIC.se"] & theta < care.res0.1[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.1[,"MA_BIC.se"]) / dim(care.res0.1)[1]
    final.out[4,3] = sum(theta > care.res0.1[,"BIC.tilde_theta"] - 1.96 * care.res0.1[,"BIC.se"] & theta < care.res0.1[,"BIC.tilde_theta"] + 1.96 * care.res0.1[,"BIC.se"]) / dim(care.res0.1)[1]
    final.out[5,3] = sum(theta > care.res0.1[,"AIC.tilde_theta"] - 1.96 * care.res0.1[,"AIC.Efron_se"] & theta < care.res0.1[,"AIC.tilde_theta"] + 1.96 * care.res0.1[,"AIC.Efron_se"]) / dim(care.res0.1)[1]
    final.out[6,3] = sum(theta > care.res0.1[,"BIC.tilde_theta"] - 1.96 * care.res0.1[,"BIC.Efron_se"] & theta < care.res0.1[,"BIC.tilde_theta"] + 1.96 * care.res0.1[,"BIC.Efron_se"]) / dim(care.res0.1)[1]
    final.out[7,3] = sum(theta > care.res0.1[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.1[,"MA_AIC.Efron_se"] & theta < care.res0.1[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.1[,"MA_AIC.Efron_se"]) / dim(care.res0.1)[1]
    final.out[8,3] = sum(theta > care.res0.1[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.1[,"MA_BIC.Efron_se"] & theta < care.res0.1[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.1[,"MA_BIC.Efron_se"]) / dim(care.res0.1)[1]

    final.out[9,3] = sum(theta > care.res0.3[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.3[,"MA_AIC.se"] & theta < care.res0.3[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.3[,"MA_AIC.se"]) / dim(care.res0.3)[1]
    final.out[10,3] = sum(theta > care.res0.3[,"AIC.tilde_theta"] - 1.96 * care.res0.3[,"AIC.se"] & theta < care.res0.3[,"AIC.tilde_theta"] + 1.96 * care.res0.3[,"AIC.se"]) / dim(care.res0.3)[1]
    final.out[11,3] = sum(theta > care.res0.3[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.3[,"MA_BIC.se"] & theta < care.res0.3[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.3[,"MA_BIC.se"]) / dim(care.res0.3)[1]
    final.out[12,3] = sum(theta > care.res0.3[,"BIC.tilde_theta"] - 1.96 * care.res0.3[,"BIC.se"] & theta < care.res0.3[,"BIC.tilde_theta"] + 1.96 * care.res0.3[,"BIC.se"]) / dim(care.res0.3)[1]
    final.out[13,3] = sum(theta > care.res0.3[,"AIC.tilde_theta"] - 1.96 * care.res0.3[,"AIC.Efron_se"] & theta < care.res0.3[,"AIC.tilde_theta"] + 1.96 * care.res0.3[,"AIC.Efron_se"]) / dim(care.res0.3)[1]
    final.out[14,3] = sum(theta > care.res0.3[,"BIC.tilde_theta"] - 1.96 * care.res0.3[,"BIC.Efron_se"] & theta < care.res0.3[,"BIC.tilde_theta"] + 1.96 * care.res0.3[,"BIC.Efron_se"]) / dim(care.res0.3)[1]
    final.out[15,3] = sum(theta > care.res0.3[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.3[,"MA_AIC.Efron_se"] & theta < care.res0.3[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.3[,"MA_AIC.Efron_se"]) / dim(care.res0.3)[1]
    final.out[16,3] = sum(theta > care.res0.3[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.3[,"MA_BIC.Efron_se"] & theta < care.res0.3[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.3[,"MA_BIC.Efron_se"]) / dim(care.res0.3)[1]

    final.out[17,3] = sum(theta > care.res0.5[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.5[,"MA_AIC.se"] & theta < care.res0.5[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.5[,"MA_AIC.se"]) / dim(care.res0.5)[1]
    final.out[18,3] = sum(theta > care.res0.5[,"AIC.tilde_theta"] - 1.96 * care.res0.5[,"AIC.se"] & theta < care.res0.5[,"AIC.tilde_theta"] + 1.96 * care.res0.5[,"AIC.se"]) / dim(care.res0.5)[1]
    final.out[19,3] = sum(theta > care.res0.5[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.5[,"MA_BIC.se"] & theta < care.res0.5[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.5[,"MA_BIC.se"]) / dim(care.res0.5)[1]
    final.out[20,3] = sum(theta > care.res0.5[,"BIC.tilde_theta"] - 1.96 * care.res0.5[,"BIC.se"] & theta < care.res0.5[,"BIC.tilde_theta"] + 1.96 * care.res0.5[,"BIC.se"]) / dim(care.res0.5)[1]
    final.out[21,3] = sum(theta > care.res0.5[,"AIC.tilde_theta"] - 1.96 * care.res0.5[,"AIC.Efron_se"] & theta < care.res0.5[,"AIC.tilde_theta"] + 1.96 * care.res0.5[,"AIC.Efron_se"]) / dim(care.res0.5)[1]
    final.out[22,3] = sum(theta > care.res0.5[,"BIC.tilde_theta"] - 1.96 * care.res0.5[,"BIC.Efron_se"] & theta < care.res0.5[,"BIC.tilde_theta"] + 1.96 * care.res0.5[,"BIC.Efron_se"]) / dim(care.res0.5)[1]
    final.out[23,3] = sum(theta > care.res0.5[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.5[,"MA_AIC.Efron_se"] & theta < care.res0.5[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.5[,"MA_AIC.Efron_se"]) / dim(care.res0.5)[1]
    final.out[24,3] = sum(theta > care.res0.5[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.5[,"MA_BIC.Efron_se"] & theta < care.res0.5[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.5[,"MA_BIC.Efron_se"]) / dim(care.res0.5)[1]

    final.out[25,3] = sum(theta > care.res0.7[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.7[,"MA_AIC.se"] & theta < care.res0.7[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.7[,"MA_AIC.se"]) / dim(care.res0.7)[1]
    final.out[26,3] = sum(theta > care.res0.7[,"AIC.tilde_theta"] - 1.96 * care.res0.7[,"AIC.se"] & theta < care.res0.7[,"AIC.tilde_theta"] + 1.96 * care.res0.7[,"AIC.se"]) / dim(care.res0.7)[1]
    final.out[27,3] = sum(theta > care.res0.7[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.7[,"MA_BIC.se"] & theta < care.res0.7[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.7[,"MA_BIC.se"]) / dim(care.res0.7)[1]
    final.out[28,3] = sum(theta > care.res0.7[,"BIC.tilde_theta"] - 1.96 * care.res0.7[,"BIC.se"] & theta < care.res0.7[,"BIC.tilde_theta"] + 1.96 * care.res0.7[,"BIC.se"]) / dim(care.res0.7)[1]
    final.out[29,3] = sum(theta > care.res0.7[,"AIC.tilde_theta"] - 1.96 * care.res0.7[,"AIC.Efron_se"] & theta < care.res0.7[,"AIC.tilde_theta"] + 1.96 * care.res0.7[,"AIC.Efron_se"]) / dim(care.res0.7)[1]
    final.out[30,3] = sum(theta > care.res0.7[,"BIC.tilde_theta"] - 1.96 * care.res0.7[,"BIC.Efron_se"] & theta < care.res0.7[,"BIC.tilde_theta"] + 1.96 * care.res0.7[,"BIC.Efron_se"]) / dim(care.res0.7)[1]
    final.out[31,3] = sum(theta > care.res0.7[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.7[,"MA_AIC.Efron_se"] & theta < care.res0.7[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.7[,"MA_AIC.Efron_se"]) / dim(care.res0.7)[1]
    final.out[32,3] = sum(theta > care.res0.7[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.7[,"MA_BIC.Efron_se"] & theta < care.res0.7[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.7[,"MA_BIC.Efron_se"]) / dim(care.res0.7)[1]

    final.out[33,3] = sum(theta > care.res0.9[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.9[,"MA_AIC.se"] & theta < care.res0.9[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.9[,"MA_AIC.se"]) / dim(care.res0.9)[1]
    final.out[34,3] = sum(theta > care.res0.9[,"AIC.tilde_theta"] - 1.96 * care.res0.9[,"AIC.se"] & theta < care.res0.9[,"AIC.tilde_theta"] + 1.96 * care.res0.9[,"AIC.se"]) / dim(care.res0.9)[1]
    final.out[35,3] = sum(theta > care.res0.9[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.9[,"MA_BIC.se"] & theta < care.res0.9[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.9[,"MA_BIC.se"]) / dim(care.res0.9)[1]
    final.out[36,3] = sum(theta > care.res0.9[,"BIC.tilde_theta"] - 1.96 * care.res0.9[,"BIC.se"] & theta < care.res0.9[,"BIC.tilde_theta"] + 1.96 * care.res0.9[,"BIC.se"]) / dim(care.res0.9)[1]
    final.out[37,3] = sum(theta > care.res0.9[,"AIC.tilde_theta"] - 1.96 * care.res0.9[,"AIC.Efron_se"] & theta < care.res0.9[,"AIC.tilde_theta"] + 1.96 * care.res0.9[,"AIC.Efron_se"]) / dim(care.res0.9)[1]
    final.out[38,3] = sum(theta > care.res0.9[,"BIC.tilde_theta"] - 1.96 * care.res0.9[,"BIC.Efron_se"] & theta < care.res0.9[,"BIC.tilde_theta"] + 1.96 * care.res0.9[,"BIC.Efron_se"]) / dim(care.res0.9)[1]
    final.out[39,3] = sum(theta > care.res0.9[,"MA_AIC.tilde_theta"] - 1.96 * care.res0.9[,"MA_AIC.Efron_se"] & theta < care.res0.9[,"MA_AIC.tilde_theta"] + 1.96 * care.res0.9[,"MA_AIC.Efron_se"]) / dim(care.res0.9)[1]
    final.out[40,3] = sum(theta > care.res0.9[,"MA_BIC.tilde_theta"] - 1.96 * care.res0.9[,"MA_BIC.Efron_se"] & theta < care.res0.9[,"MA_BIC.tilde_theta"] + 1.96 * care.res0.9[,"MA_BIC.Efron_se"]) / dim(care.res0.9)[1]

  
    # bias
    final.out[1:8,4] = (theta -  colMeans(care.res0.1[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])) #/ theta   
    final.out[9:16,4] = (theta -  colMeans(care.res0.3[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])) #/ theta
    final.out[17:24,4] = (theta -  colMeans(care.res0.5[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])) #/ theta
    final.out[25:32,4] = (theta -  colMeans(care.res0.7[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])) #/ theta
    final.out[33:40,4] = (theta -  colMeans(care.res0.9[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")])) #/ theta


    # MSE
    final.out[1:8,5] = colMeans((care.res0.1[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")] -theta )^2)
    final.out[9:16,5] = colMeans((care.res0.3[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")] -theta )^2)
    final.out[17:24,5] = colMeans((care.res0.5[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")] -theta )^2)
    final.out[25:32,5] = colMeans((care.res0.7[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")] -theta )^2)
    final.out[33:40,5] = colMeans((care.res0.9[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")] -theta )^2)

    # Monte SD
    final.out[1:8,6] = apply(care.res0.1[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")],2,sd)
    final.out[9:16,6] = apply(care.res0.3[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")],2,sd)
    final.out[17:24,6] = apply(care.res0.5[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")],2,sd)
    final.out[25:32,6] = apply(care.res0.7[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")],2,sd)
    final.out[33:40,6] = apply(care.res0.9[,c("MA_AIC.tilde_theta","AIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","AIC.tilde_theta","BIC.tilde_theta","MA_AIC.tilde_theta","MA_BIC.tilde_theta")],2,sd)

    # SD
    final.out[1:8,7] = colMeans(care.res0.1[,c("MA_AIC.se","AIC.se","MA_BIC.se","BIC.se","AIC.Efron_se","BIC.Efron_se","MA_AIC.Efron_se","MA_BIC.Efron_se")])
    final.out[9:16,7] = colMeans(care.res0.3[,c("MA_AIC.se","AIC.se","MA_BIC.se","BIC.se","AIC.Efron_se","BIC.Efron_se","MA_AIC.Efron_se","MA_BIC.Efron_se")])
    final.out[17:24,7] = colMeans(care.res0.5[,c("MA_AIC.se","AIC.se","MA_BIC.se","BIC.se","AIC.Efron_se","BIC.Efron_se","MA_AIC.Efron_se","MA_BIC.Efron_se")])
    final.out[25:32,7] = colMeans(care.res0.7[,c("MA_AIC.se","AIC.se","MA_BIC.se","BIC.se","AIC.Efron_se","BIC.Efron_se","MA_AIC.Efron_se","MA_BIC.Efron_se")])
    final.out[33:40,7] = colMeans(care.res0.9[,c("MA_AIC.se","AIC.se","MA_BIC.se","BIC.se","AIC.Efron_se","BIC.Efron_se","MA_AIC.Efron_se","MA_BIC.Efron_se")])

    savefile = paste0("/rsrch5/home/biostatistics/chongwulab/wzhang24/CARE/simulation_final_2/summaryRes/",simSetting,"theta",theta,"N", N, "prop_invalid",prop_invalid,".RData")
    save(final.out, file = savefile)
}


for (simSetting in c("simRes12")){
  for(theta in c(0.1, 0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.5)) {
        N = 5e+05
        try(summary_fun(theta,N, prop_invalid,simSetting))    
        cat("Finish theta", theta, "  Prop",prop_invalid,"\n")
    }
  }
}