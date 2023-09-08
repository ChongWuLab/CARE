rm(list= ls())
#library(ggplot2)

summary_fun <- function(theta,propInvalid,simSetting) {
    
    #theta = 0.03
    #propInvalid = 0.3
    #simSetting = "simRes1"
    #setwd(paste0("/gpfs/research/chongwu/Chong/CARE/Submit/simulation-final/",simSetting))
    setwd(paste0("C:\\Users\\Evelyn\\Dropbox\\CARE\\simulation_final_2\\", simSetting))
    ind = paste0("theta",theta,"thetaU0.3prop_invalid",propInvalid)

    care_result = list()
    care_no_result = list()

    ConMix_result = list()
    MR_Lasso_result = list()
    MRMix_result = list()
    #MR_PRESSO_result = list()
    MRAPSS_result = list()
    MRcML_result = list()
    TSMR_result = NULL
    run_time = list()
    settings = list()
    files = list.files(getwd())

    files = files[grepl(ind,files)]
    
    # only using sample size = 1e5
    #files = files[grepl("N1",files)]


    for(set.ind in 1:length(files))
    {
      load(files[set.ind])
      
      run_time = c(run_time,runningtime)
      settings = c(settings, setting)
    }
    
    nsim = length(run_time)
    # running time
    run_time2 = matrix(NA, nsim, 14)
    
    for(i in 1:nsim) {
        run_time2[i,] = run_time[[i]]
    }
    
    colnames(run_time2) = names(run_time[[i]])
    run_time = colMeans(run_time2)

    #savefile = paste0("/gpfs/research/chongwu/Chong/CARE/Submit/simulation-final/summaryRes/res_runtime_",simSetting,"theta",theta,"prop_invalid",propInvalid,".RData")
    savefile = paste0("C:\\Users\\Evelyn\\Dropbox\\CARE\\simulation_final_2\\summaryRes\\res_runtime_",simSetting,"theta",theta,"prop_invalid",propInvalid,".RData")
    save(setting,run_time,run_time2, file = savefile)
}

for (simSetting in c("simRes1")){
  for(theta in c(0.1,0.07, 0.03, 0, -0.03, -0.07, -0.1)) {
      for (propInvalid in c(0.3, 0.5, 0.7)) {
          try(summary_fun(theta,propInvalid,simSetting))
          cat("Finish theta", theta, " Prop",propInvalid,"\n")
      }
  }
}

