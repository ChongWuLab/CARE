require(TwoSampleMR)
require(data.table)
require(dplyr)


dir = "CARE/covid19/smmry_res"
setwd(dir)
files = list.files(dir)


#combine the results
CARE.res = as.data.frame(matrix(NA,length(files),33))


CARE.theta = 8 #"BIC.tilde_theta"
CARE.se = 9 #"BIC.se"
CARE.p = 10 #"BIC.p"

#CARE.se = 13 #"BIC.se" Efron
#CARE.p = 14 #"BIC.p" Efron

for(j in 1:length(files)) {
    
    tryCatch({

    file = files[j]
    load(file)

    name = gsub(".RData","",file)

    CARE.res[j,1] = name
    
    tmp = CARE.pruning[[1]]

    tmp = tmp[!is.na(tmp[,CARE.p]),,drop=FALSE]
    if(dim(tmp)[1] > 0) {
        CARE.res[j,2:33] = tmp[1,]
        
    } else {
        cat(j,"wrong check\n")
    }
    }, error = function(e) {
        cat("Probematic: ",j,"\n")
        cat("ERROR :", conditionMessage(e), "\n")
    })
    
}

colnames(CARE.res) = c("Name",colnames(tmp))


CARE.res = CARE.res[!is.na(CARE.res[,1]),]
CARE.res$exposure.id = gsub("_.*","",CARE.res[,1])
CARE.res$outcome.id = CARE.res[,1]
CARE.res$exposure.name = NA
CARE.res$outcome.name = NA



library(ieugwasr)
info = gwasinfo()
info = as.data.frame(info)

for(j in 1:dim(CARE.res)[1]) {
    tryCatch({
        
    CARE.res[j,"exposure.name"] = info[info[,1]==CARE.res[j,"exposure.id"],"trait"]
    CARE.res[j,"outcome.id"] = gsub(paste0(CARE.res[j,"exposure.id"],"_"),"",CARE.res[j,"outcome.id"])
    
    #CARE.res[j,"outcome.name"] = info[info[,1]==CARE.res[j,"outcome.id"],"trait"]

    }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
    })
}



# other competing methods:
MR.res = as.data.frame(matrix(NA,length(files),43))

colnames(MR.res) = c("MRcML_MA_BIC_beta","MRcML_MA_BIC_se","MRcML_MA_BIC_p","DP_MA_BIC_beta","DP_MA_BIC_se","DP_MA_BIC_p","ComMix_beta","ConMix_se","ConMix_p","MRLasso_beta","MRLasso_se","MRLasso_p","IVW_beta","IVW_se","IVW_p","Egger_beta","Egger_se","Egger_p","Weighted_median_beta","Weighted_median_se","Weighted_median_p","Weighted_mode_beta","Weighted_mode_se","Weighted_mode_p","RAPS_beta","RAPS_se","RAPS_p","MRmix_beta","MRmix_se","MRmix_p","MRAPSS_beta","MRAPSS_se","MRAPSS_p","CARE_MA_BIC_beta","CARE_MA_BIC_se","CARE_MA_BIC_p","CARE_BIC_beta","CARE_BIC_se","CARE_BIC_p","IVW_fixed_beta","IVW_fixed_se","IVW_fixed_p","Trait")

for(j in 1:length(files)) {
    
    file = files[j]
    load(file)
    name = gsub(".RData","",file)
    MR.res[j,43] = name
    MR = MR.clumping
    MR = MR$res
    MRAPSS = MRAPSS.pruning$MRAPSS_res
    if(is.null(MR)) {
        next
    }
    
    tmp = unlist(MR[[1]])
    
    MR.res[j,1:6] = unlist(MR[[1]])[c("MA_BIC_theta","MA_BIC_se","MA_BIC_p","MA_BIC_DP_theta","MA_BIC_DP_se","MA_BIC_DP_p")]
    
    #ConMix
    MR.res[j,7:9] = c(MR[[2]]@Estimate, MR[[2]]@Psi, MR[[2]]@Pvalue)
    
    
    #MRLasso
    if(is.null(MR[[3]])) {
        
    } else {
        MR.res[j,10:12] = c(MR[[3]]@Estimate,MR[[3]]@StdError, MR[[3]]@Pvalue)
    }
    
    
    #
    TSMR_result = MR[[9]]
    
    tmp = MR[[9]]
    names(TSMR_result) = c("IVW_beta","IVW_se","Egger_beta","Egger_se","Weighted_median_beta","Weighted_median_se","Weighted_mode_beta","Weighted_mode_se","RAPS_beta","RAPS_se")
    TSMR_result = as.data.frame(TSMR_result)
    TSMR_result = t(TSMR_result)
    
    #pIVW
    MR.res[j,13:15] = c(tmp[1:2],pnorm(-abs(TSMR_result[,1] / TSMR_result[,2])) *2)
    
    
    #pval_EGGER
    MR.res[j,16:18] = c(tmp[3:4],pnorm(-abs(TSMR_result[,3] / TSMR_result[,4]))*2)

    
    #pval_WEIGHTED_MEDIAN
    MR.res[j,19:21] = c(tmp[5:6], pnorm(-abs(TSMR_result[,5] / TSMR_result[,6]))*2)
    
    #pval_WEIGHTED_MODE
    MR.res[j,22:24] = c(tmp[7:8], pnorm(-abs(TSMR_result[,7] / TSMR_result[,8]))*2)
    
    
    #pval_RAPS
    MR.res[j,25:27] = c(tmp[9:10], pnorm(-abs(TSMR_result[,9] / TSMR_result[,10]))*2)

    #MRmix
    tmp = unlist(MR[[10]])
    MR.res[j,28:30] = tmp[c(1,4,6)]
    
    #MR-APSS
    MR.res[j, 31:33] = c(MRAPSS$beta, MRAPSS$beta.se, as.numeric(MRAPSS$pvalue))
    
    #CARE:
    tmp = MR[[11]]
    tmp = unlist(tmp[[1]][[1]])
    
    # need to change this
    MR.res[j,34:36] = tmp[c("MA_BIC.tilde_theta","MA_BIC.se","MA_BIC.p")]
    MR.res[j,37:39] = tmp[c("BIC.tilde_theta","BIC.se","BIC.p")]
    
    # IVW fixed
    tmp = MR[[12]]
    tmp = unlist(tmp)
    
    MR.res[j,40:42] = tmp[1:3]
}


colSums(MR.res<0.05/100,na.rm=T)
colSums(CARE.res<0.05/100,na.rm=T)

CARE.res = CARE.res[!is.na(CARE.res[,1]),]
MR.res = MR.res[indx,]
MR.res = MR.res[!is.na(MR.res[,1]),]
# add the trait information
save(CARE.res,MR.res,file = "CARE/covid19/res_Ana/COVID19_summarized.RData")
