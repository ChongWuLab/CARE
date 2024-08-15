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

    care_result = list()
    care_no_result = list()

    ConMix_result = list()
    MR_Lasso_result = list()
    MRMix_result = list()
    #MR_PRESSO_result = list()
    MRcML_result = list()
    TSMR_result = NULL
    MRAPSS_result = list()
    run_time = list()
    settings = list()

    files = list.files(getwd())
    files = files[grepl(ind1,files)]
    files = files[grepl(ind2,files)]
    
    # only using sample size = 1e5
    #files = files[grepl("N1",files)]


    for(set.ind in 1:length(files))
    {
      #load(paste("sim_results/simres_MRMethods",set.ind, "N50000pthr5e-08pi10.01theta",theta,"thetaU0.3prop_invalid0.5NxNy_ratio2_July22.Rdata",sep=""))
      load(files[set.ind])
      
      run_time = c(run_time,runningtime)
      settings = c(settings, setting)
      care_result = c(care_result,care_sim_result) # check
      
      care_no_result = c(care_no_result, CAREno_sim_result)
      ConMix_result =
        c(ConMix_result,ContMix_sim_result)
      MR_Lasso_result =
        c(MR_Lasso_result,MRLasso_sim_Result)
      MRMix_result =
        c(MRMix_result,MRMix_sim_result)
      # MR_PRESSO_result =
      #   c(MR_PRESSO_result,MRPresso_sim_result)
      MRcML_result =
        c(MRcML_result,MRcML_sim_result)
      TSMR_result =
        rbind(TSMR_result,TwoSampleMR_sim_result)
      MRAPSS_result = c(MRAPSS_result, MRAPSS_sim_result)
    }
    
    ## CARE with eta =0.5, eta_var = 0
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
    
    
    ## CARE with eta =0.5, eta_var = 0, no winner's curse correction
    nsim = length(care_no_result)
    nsim

    care2.res = matrix(NA,nsim,32)
    for(i in 1:nsim)
    {
       tmp = care_no_result[[i]]
       tmp = tmp[1]
       
       tmp = unlist(tmp)
       care2.res[i,] = tmp
    }

    tmpname = names(tmp)
    colnames(care2.res) = gsub("res.","",tmpname)
    #colSums(care.res<0.05) / nsim
    #colMeans(care.res)
    #quantile(care2.res[,"MA_BIC.tilde_theta"])
    #which(care2.res[,"MA_BIC.tilde_theta"]>1)
    #quantile(care2.res[,"GIC.hat_theta"])
    #quantile(care2.res[,"BIC.tilde_theta"])
    #indxtmp = which(care2.res[,"GIC.hat_theta"]>1)
    #care2.res[indxtmp,]
    
    
    # summarize the setting results
    #setting2 = matrix(NA, nsim, length(settings[[i]]))
    
    #for(i in 1:nsim) {
   #     setting2[i,] = settings[[i]]
    #}
   # colnames(setting2) = names(settings[[1]])
    
    #setting2 = colMeans(setting2)
    
    #setting = c(setting2,tmp[29:38])
    
    # running time
    #run_time2 = matrix(NA, nsim, 13)
    run_time2 = matrix(NA, nsim, 14)
    for(i in 1:nsim) {
        run_time2[i,] = run_time[[i]]
    }
    
    colnames(run_time2) = names(run_time[[i]])
    run_time = colMeans(run_time2)

    ### MRcML, MRmix, ConMix
    nsim = length(MRcML_result)
    MRcML = matrix(NA,nsim,9)
    colnames(MRcML) = c("MRcML_MA_BIC_beta","MRcML_MA_BIC_se","MRcML_MA_BIC_p","DP_MA_BIC_beta","DP_MA_BIC_se","DP_MA_BIC_p","MRmix_beta","MRmix_se","MRmix_p")
    MRcML = as.data.frame(MRcML)

    for(i in 1:nsim)
    {
        MRcML[i,] =c(MRcML_result[[i]]$MA_BIC_theta,MRcML_result[[i]]$MA_BIC_se,MRcML_result[[i]]$MA_BIC_p,MRcML_result[[i]]$MA_BIC_DP_theta,MRcML_result[[i]]$MA_BIC_DP_se,MRcML_result[[i]]$MA_BIC_DP_p,MRMix_result[[i]]$theta,MRMix_result[[i]]$SE_theta,MRMix_result[[i]]$pvalue_theta)
    }
    
    
    ### Two Sample MR
    TSMR_result = as.data.frame(TSMR_result)
    colnames(TSMR_result) = c("IVW_beta","IVW_se","Egger_beta","Egger_se","Weighted_median_beta","Weighted_median_se","Weighted_mode_beta","Weighted_mode_se","RAPS_beta","RAPS_se")
  
    TSMR_result$pval_IVW = pnorm(-abs(TSMR_result[,1] / TSMR_result[,2]))*2
    TSMR_result$pval_EGGER = pnorm(-abs(TSMR_result[,3] / TSMR_result[,4]))*2
    TSMR_result$pval_WEIGHTED_MEDIAN = pnorm(-abs(TSMR_result[,5] / TSMR_result[,6]))*2
    TSMR_result$pval_WEIGHTED_MODE = pnorm(-abs(TSMR_result[,7] / TSMR_result[,8]))*2
    TSMR_result$pval_RAPS = pnorm(-abs(TSMR_result[,9] / TSMR_result[,10]))*2


    
    ###  ConMix
    nsim = length(ConMix_result)
    ConMix = matrix(NA,nsim,5)
    colnames(ConMix) = c("ConMix_beta","ConMix_se","ConMix_p","ConMix_CILower","ConMix_CIUpper")

    for(i in 1:nsim)
    {
        ConMix[i,] = c(ConMix_result[i]$ConMix_result@Estimate,ConMix_result[i]$ConMix_result@Psi,ConMix_result[i]$ConMix_result@Pvalue,min(ConMix_result[i]$ConMix_result@CILower),max(ConMix_result[i]$ConMix_result@CIUpper))
    }
    ConMix = ConMix[!is.na(ConMix[,1]),]


    ### MR_Lasso
    nsim = length(MR_Lasso_result)

    MRLasso = as.data.frame(matrix(NA,nsim,3))
    colnames(MRLasso) = c("MRLasso_beta","MRLasso_se","MRLasso_p")
    for(i in 1:nsim)
    {
      if(is.null(MR_Lasso_result[i]$MR_Lasso_result))
      {
        #pval_MR_Lasso = c(pval_MR_Lasso,NA)
        next
      }
      MRLasso[i,] = c(MR_Lasso_result[i]$MR_Lasso_result@Estimate,MR_Lasso_result[i]$MR_Lasso_result@StdError, MR_Lasso_result[i]$MR_Lasso_result@Pvalue)
    }
    MRLasso = MRLasso[!is.na(MRLasso[,1]),]

    ### MRAPSS
    clean_list <- Filter(Negate(is.null), MRAPSS_result)
    nsim = length(clean_list)
    MRAPSS.res = matrix(NA,nsim,4)
    storage.mode(MRAPSS.res) <- "numeric"
    colnames(MRAPSS.res) = c('beta','se','p','sigma.sq')
    for (i in 1:(nsim)) {
      MRAPSS.res[i,] = as.numeric(c(clean_list[[i]]$beta, 
                                    clean_list[[i]]$beta.se, 
                                    clean_list[[i]]$pvalue, 
                                    clean_list[[i]]$sigma.sq))
    }
    
    cutoff = 0.05
    # final.out = as.data.frame(matrix(NA,21,7))
    # colnames(final.out) = c("Power","Estimate","Coverage","Bias","MSE","Monte_SD","SD")
    # 
    # rownames(final.out) = c("MRcML","MRcML_DP","MRmix","IVW","Egger","Weighted_Median","Weighted_Mode","RAPS","ConMix","PRESSO","MRLasso","CARE_MA_GIC","CARE_GIC","CARE_MA_BIC","CARE_BIC","CARE_Efron_GIC","CARE_Efron_BIC","CARE_Efron_MA_GIC","CARE_Efron_MA_BIC","CARE_no_BIC","CARE_no_BIC_Efron")
    final.out = as.data.frame(matrix(NA,22,7))
    rownames(final.out) = c("MRcML","MRcML_DP","MRmix","IVW","Egger","Weighted_Median","Weighted_Mode","RAPS","ConMix","PRESSO","MRLasso","CARE_MA_GIC","CARE_GIC","CARE_MA_BIC","CARE_BIC","CARE_Efron_GIC","CARE_Efron_BIC","CARE_Efron_MA_GIC","CARE_Efron_MA_BIC","CARE_no_BIC","CARE_no_BIC_Efron","MRAPSS")
    colnames(final.out) = c("Power","Estimate","Coverage","Bias","MSE","Monte_SD","SD")
     # power
    final.out[1,1] = sum(MRcML$MRcML_MA_BIC_p < cutoff)/dim(MRcML)[1]
    final.out[2,1] = sum(MRcML$DP_MA_BIC_p < cutoff)/dim(MRcML)[1]
    final.out[3,1] = sum(MRcML$MRmix_p < cutoff)/dim(MRcML)[1]

    final.out[4:8,1] = colSums(TSMR_result[,11:15] < cutoff) / dim(TSMR_result)[1]
    final.out[9,1] = sum(ConMix[,3] < cutoff) / dim(ConMix)[1]
    #final.out[10,1] = sum(PRESSO[,3] < cutoff) / dim(PRESSO)[1]
    final.out[11,1] = sum(MRLasso[,3] < cutoff) / dim(MRLasso)[1]

    final.out[12:19,1] = colSums(care.res[,c("MA_GIC.p","GIC.p","MA_BIC.p","BIC.p","GIC.Efron_p","BIC.Efron_p","MA_GIC.Efron_p","MA_BIC.Efron_p")] < cutoff) / dim(care.res)[1] #c("MA_GIC.p","GIC.p","MA_BIC.p","BIC.p")

    final.out[20:21,1] = colSums(care2.res[,c("BIC.p","BIC.Efron_p")] < cutoff) / dim(care2.res)[1]
    final.out[22,1] = sum(MRAPSS.res[,"p"] < cutoff) / dim(MRAPSS.res)[1]

    # Estimate
    final.out[1,2] = mean(MRcML$MRcML_MA_BIC_beta)
    final.out[2,2] = mean(MRcML$DP_MA_BIC_beta)
    final.out[3,2] = mean(MRcML$MRmix_beta)

    final.out[4:8,2] = colMeans(TSMR_result[,c(1,3,5,7,9)])
    final.out[9,2] = mean(ConMix[,1])
    #final.out[10,2] = mean(PRESSO[,1])
    final.out[11,2] = mean(MRLasso[,1])

    final.out[12:19,2] = colMeans(care.res[,c("MA_GIC.tilde_theta","GIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","GIC.tilde_theta","BIC.tilde_theta","MA_GIC.tilde_theta","GIC.tilde_theta")])

    final.out[20:21,2] = colMeans(care2.res[,c("BIC.tilde_theta","BIC.tilde_theta")])
    
    final.out[22,2] = mean(MRAPSS.res[,"beta"])
    
    
    # Coverage
    final.out[1,3] = sum(theta > MRcML$MRcML_MA_BIC_beta - 1.96 * MRcML$MRcML_MA_BIC_se & theta < MRcML$MRcML_MA_BIC_beta + 1.96 * MRcML$MRcML_MA_BIC_se) / dim(MRcML)[1]
    final.out[2,3] =  sum(theta > MRcML$DP_MA_BIC_beta - 1.96 * MRcML$DP_MA_BIC_se & theta < MRcML$DP_MA_BIC_beta + 1.96 * MRcML$DP_MA_BIC_se) / dim(MRcML)[1]
    final.out[3,3] =  sum(theta > MRcML$MRmix_beta - 1.96 * MRcML$MRmix_se & theta < MRcML$MRmix_beta + 1.96 * MRcML$MRmix_se) / dim(MRcML)[1]

    final.out[4,3] = sum(theta > TSMR_result$IVW_beta - 1.96 * TSMR_result$IVW_se & theta < TSMR_result$IVW_beta + 1.96 * TSMR_result$IVW_se) / dim(TSMR_result)[1]
    final.out[5,3] = sum(theta > TSMR_result$Egger_beta - 1.96 * TSMR_result$Egger_se & theta < TSMR_result$Egger_beta + 1.96 * TSMR_result$Egger_se) / dim(TSMR_result)[1]
    final.out[6,3] = sum(theta > TSMR_result$Weighted_median_beta - 1.96 * TSMR_result$Weighted_median_se & theta < TSMR_result$Weighted_median_beta + 1.96 * TSMR_result$Weighted_median_se) / dim(TSMR_result)[1]
    final.out[7,3] = sum(theta > TSMR_result$Weighted_mode_beta - 1.96 * TSMR_result$Weighted_mode_se & theta < TSMR_result$Weighted_mode_beta + 1.96 * TSMR_result$Weighted_mode_se) / dim(TSMR_result)[1]
    final.out[8,3] = sum(theta > TSMR_result$RAPS_beta - 1.96 * TSMR_result$RAPS_se & theta < TSMR_result$RAPS_beta + 1.96 * TSMR_result$RAPS_se) / dim(TSMR_result)[1]

    final.out[9,3] = sum(theta > ConMix[,4] & theta < ConMix[,5]) / dim(ConMix)[1]
    #final.out[10,3] = sum(theta > PRESSO[,1] - 1.96 * PRESSO[,2] & theta < PRESSO[,1] + 1.96 * PRESSO[,2]) / dim(PRESSO)[1]
    final.out[11,3] = sum(theta > MRLasso[,1] - 1.96 * MRLasso[,2] & theta < MRLasso[,1] + 1.96 * MRLasso[,2]) / dim(MRLasso)[1]

    final.out[12,3] = sum(theta > care.res[,"MA_GIC.tilde_theta"] - 1.96 * care.res[,"MA_GIC.se"] & theta < care.res[,"MA_GIC.tilde_theta"] + 1.96 * care.res[,"MA_GIC.se"]) / dim(care.res)[1]

    final.out[13,3] = sum(theta > care.res[,"GIC.tilde_theta"] - 1.96 * care.res[,"GIC.se"] & theta < care.res[,"GIC.tilde_theta"] + 1.96 * care.res[,"GIC.se"]) / dim(care.res)[1]

    final.out[14,3] = sum(theta > care.res[,"MA_BIC.tilde_theta"] - 1.96 * care.res[,"MA_BIC.se"] & theta < care.res[,"MA_BIC.tilde_theta"] + 1.96 * care.res[,"MA_BIC.se"]) / dim(care.res)[1]

    final.out[15,3] = sum(theta > care.res[,"BIC.tilde_theta"] - 1.96 * care.res[,"BIC.se"] & theta < care.res[,"BIC.tilde_theta"] + 1.96 * care.res[,"BIC.se"]) / dim(care.res)[1]

    
    final.out[16,3] = sum(theta > care.res[,"GIC.tilde_theta"] - 1.96 * care.res[,"GIC.Efron_se"] & theta < care.res[,"GIC.tilde_theta"] + 1.96 * care.res[,"GIC.Efron_se"]) / dim(care.res)[1]

    final.out[17,3] = sum(theta > care.res[,"BIC.tilde_theta"] - 1.96 * care.res[,"BIC.Efron_se"] & theta < care.res[,"BIC.tilde_theta"] + 1.96 * care.res[,"BIC.Efron_se"]) / dim(care.res)[1]

    final.out[18,3] = sum(theta > care.res[,"MA_GIC.tilde_theta"] - 1.96 * care.res[,"MA_GIC.Efron_se"] & theta < care.res[,"MA_GIC.tilde_theta"] + 1.96 * care.res[,"MA_GIC.Efron_se"]) / dim(care.res)[1]

    final.out[19,3] = sum(theta > care.res[,"MA_BIC.tilde_theta"] - 1.96 * care.res[,"MA_BIC.Efron_se"] & theta < care.res[,"MA_BIC.tilde_theta"] + 1.96 * care.res[,"MA_BIC.Efron_se"]) / dim(care.res)[1]

    final.out[20,3] = sum(theta > care2.res[,"BIC.tilde_theta"] - 1.96 * care2.res[,"BIC.se"] & theta < care2.res[,"BIC.tilde_theta"] + 1.96 * care2.res[,"BIC.se"]) / dim(care2.res)[1]

    final.out[21,3] = sum(theta > care2.res[,"BIC.tilde_theta"] - 1.96 * care2.res[,"BIC.Efron_se"] & theta < care2.res[,"BIC.tilde_theta"] + 1.96 * care2.res[,"BIC.Efron_se"]) / dim(care2.res)[1]
    
    final.out[22, 3] <- sum(theta > MRAPSS.res[, 1] - 1.96 * MRAPSS.res[, 2] &
                              theta < MRAPSS.res[, 1] + 1.96 * MRAPSS.res[, 2]) / dim(MRAPSS.res)[1]
    
    
    # bias
    final.out[1,4] = (theta - mean(MRcML$MRcML_MA_BIC_beta)) #/ theta
    final.out[2,4] = (theta - mean(MRcML$DP_MA_BIC_beta)) #/ theta
    final.out[3,4] = (theta - mean(MRcML$MRmix_beta)) #/ theta

    final.out[4:8,4] = c(theta - colMeans(TSMR_result[,c(1,3,5,7,9)])) #/ theta
    final.out[9,4] = (theta - mean(ConMix[,1])) #/ theta
    #final.out[10,4] = (theta - mean(PRESSO[,1])) #/ theta
    final.out[11,4] = (theta - mean(MRLasso[,1])) #/ theta


    final.out[12:19,4] = (theta -  colMeans(care.res[,c("MA_GIC.tilde_theta","GIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","GIC.tilde_theta","BIC.tilde_theta","MA_GIC.tilde_theta","GIC.tilde_theta")])) #/ theta
    
    final.out[20:21,4] = (theta -  colMeans(care2.res[,c("BIC.tilde_theta","BIC.tilde_theta")])) #/ theta
    final.out[22,4] = (theta - mean(MRAPSS.res[,1]))

    # MSE
    final.out[1,5] = mean((MRcML$MRcML_MA_BIC_beta -theta)^2)
    final.out[2,5] = mean((MRcML$DP_MA_BIC_beta -theta)^2)
    final.out[3,5] = mean((MRcML$MRmix_beta -theta)^2)

    final.out[4:8,5] = colMeans((TSMR_result[,c(1,3,5,7,9)]-theta)^2)
    final.out[9,5] = mean((ConMix[,1]-theta)^2)
    #final.out[10,5] = mean((PRESSO[,1]-theta)^2)
    final.out[11,5] = mean((MRLasso[,1]-theta)^2)

    final.out[12:19,5] = colMeans((care.res[,c("MA_GIC.tilde_theta","GIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","GIC.tilde_theta","BIC.tilde_theta","MA_GIC.tilde_theta","GIC.tilde_theta")] -theta )^2)

    final.out[20:21,5] = colMeans((care2.res[,c("BIC.tilde_theta","BIC.tilde_theta")] -theta )^2)
    final.out[22, 5] <- mean((MRAPSS.res[, 1] - theta)^2)
    # Monte SD
    final.out[1,6] = sd(MRcML$MRcML_MA_BIC_beta)
    final.out[2,6] = sd(MRcML$DP_MA_BIC_beta)
    final.out[3,6] = sd(MRcML$MRmix_beta)

    final.out[4:8,6] = apply(TSMR_result[,c(1,3,5,7,9)],2,sd)
    final.out[9,6] = sd(ConMix[,1])
    #final.out[10,6] = sd(PRESSO[,1])
    final.out[11,6] = sd(MRLasso[,1])

    final.out[12:19,6] = apply(care.res[,c("MA_GIC.tilde_theta","GIC.tilde_theta","MA_BIC.tilde_theta","BIC.tilde_theta","GIC.tilde_theta","BIC.tilde_theta","MA_GIC.tilde_theta","GIC.tilde_theta")],2,sd)

    final.out[20:21,6] = apply(care2.res[,c("BIC.tilde_theta","BIC.tilde_theta")],2,sd)
    final.out[22,6] <- sd(MRAPSS.res[,1])
    # SD
    final.out[1,7] = mean(MRcML$MRcML_MA_BIC_se)
    final.out[2,7] = mean(MRcML$DP_MA_BIC_se)
    final.out[3,7] = mean(MRcML$MRmix_se)

    final.out[4:8,7] = colMeans(TSMR_result[,c(2,4,6,8,10)])
    final.out[9,7] = mean(ConMix[,"ConMix_se"])
    #final.out[10,7] = mean(PRESSO[,"Presso_se"])
    final.out[11,7] = mean(MRLasso[,"MRLasso_se"])

    final.out[12:19,7] = colMeans(care.res[,c("MA_GIC.se","GIC.se","MA_BIC.se","BIC.se","GIC.Efron_se","BIC.Efron_se","MA_GIC.Efron_se","MA_BIC.Efron_se")])
    
    final.out[20:21,7] = colMeans(care2.res[,c("BIC.se","BIC.Efron_se")])
    final.out[22,7] = mean(MRAPSS.res[,2])
    savefile = paste0("/rsrch5/home/biostatistics/chongwulab/wzhang24/CARE/simulation_final_2/summaryRes/",simSetting,"theta",theta,"N", N, "prop_invalid",prop_invalid,".RData")
    save(final.out,setting,run_time, file = savefile)
}



for (simSetting in c("simRes1","simRes2","simRes3","simRes4")){
  for(theta in c(0.1,0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.3, 0.5, 0.7)) {
        try(summary_fun(theta,prop_invalid,simSetting))
        
        cat("Finish theta", theta, " Prop",prop_invalid,"\n")
    }
  }
}

for (simSetting in c("simRes5")){
  for(theta in c(0.1,0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.5)) {
        try(summary_fun(theta,prop_invalid,simSetting))
        
        cat("Finish theta", theta, " Prop",prop_invalid,"\n")
    }
  }
}

for (simSetting in c("simRes2")){
  for(theta in c(0.1,0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (N in c("5000", "10000", "50000")) {
      for (prop_invalid in c(0.5)) {
          try(summary_fun(theta,N, prop_invalid,simSetting))
          
          cat("Finish theta", theta, " ", N, " Prop",prop_invalid,"\n")
      }
    }
  }
}

for (simSetting in c("simRes2")){
  for(theta in c(0.1,0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (N in c(5e+05)) {
      for (prop_invalid in c(0.9)) {
          try(summary_fun(theta,N, prop_invalid,simSetting))      
          cat("Finish theta", theta, " ", N, " Prop",prop_invalid,"\n")
      }
    }
  }
}

for (simSetting in c("simRes6")){
  for(theta in c(0.1, 0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.5)) {
        N = 10000
        try(summary_fun(theta,N, prop_invalid,simSetting))    
        cat("Finish theta", theta, "  Prop",prop_invalid,"\n")
    }
    
  }
}

for (simSetting in c("simRes7")){
  for(theta in c(0.1, 0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.5)) {
        N = 5e+05
        try(summary_fun(theta,N, prop_invalid,simSetting))    
        cat("Finish theta", theta, "  Prop",prop_invalid,"\n")
    }
    
  }
}

for (simSetting in c("simRes8")){
  for(theta in c(0.1, 0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.5)) {
        N = 5e+05
        try(summary_fun(theta,N, prop_invalid,simSetting))    
        cat("Finish theta", theta, "  Prop",prop_invalid,"\n")
    }
    
  }
}

for (simSetting in c("simRes8")){
  for(theta in c(0.1, 0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
    for (prop_invalid in c(0.3)) {
        N = 10000
        try(summary_fun(theta,N, prop_invalid,simSetting))    
        cat("Finish theta", theta, "  Prop",prop_invalid,"\n")
    }
  }
}