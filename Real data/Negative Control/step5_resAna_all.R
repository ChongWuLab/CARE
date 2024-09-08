require(TwoSampleMR)
require(data.table)
require(dplyr)

dir = '/rsrch5/scratch/biostatistics/wzhang24/CARE/negativecontrol/smmry_res'
source("Sfig1_qq_plot_inflated_type1.R")
setwd(dir)
files = list.files(dir,pattern='.RData')

#files2 = files[!(grepl("ieu-a-1001",files) | grepl("ieu-a-996",files) | grepl("ieu-a-836",files) | grepl("ieu-a-302",files) | grepl("ukb-d-1747_6",files))]
#less than 3 million
#combine the results
final.res = as.data.frame(matrix(NA,length(files),33))

final.res2 = as.data.frame(matrix(NA,length(files),16))

# ACAT function:
ACAT<-function(Pvals,weights=NULL,is.check=TRUE){
    Pvals<-as.matrix(Pvals)
    if (is.check){
        #### check if there is NA
        if (sum(is.na(Pvals))>0){
            stop("Cannot have NAs in the p-values!")
        }
        #### check if Pvals are between 0 and 1
        if ((sum(Pvals<0)+sum(Pvals>1))>0){
            stop("P-values must be between 0 and 1!")
        }
        #### check if there are pvals that are either exactly 0 or 1.
        is.zero<-(colSums(Pvals==0)>=1)
        is.one<-(colSums(Pvals==1)>=1)
        if (sum((is.zero+is.one)==2)>0){
            stop("Cannot have both 0 and 1 p-values in the same column!")
        }

        if (sum(is.zero)>0){
            warning("There are p-values that are exactly 0!")
        }
        if (sum(is.one)>0){
            warning("There are p-values that are exactly 1!")
        }

    }
    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(weights)){
        is.weights.null<-TRUE
    }else{
        is.weights.null<-FALSE
        weights<-as.matrix(weights)
        if (sum(dim(weights)!=dim(Pvals))>0){
            stop("The dimensions of weights and Pvals must be the same!")
        }else if (is.check & (sum(weights<0)>0)){
            stop("All the weights must be nonnegative!")
        }else{
            w.sum<-colSums(weights)
            if (sum(w.sum<=0)>0){
                stop("At least one weight should be positive in each column!")
            }else{
                for (j in 1:ncol(weights)){
                    weights[,j]<-weights[,j]/w.sum[j]
                }
            }
        }

    }

    #### check if there are very small non-zero p values and calcuate the cauchy statistics
    is.small<-(Pvals<1e-15)
    if (is.weights.null){
         Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
         Pvals[is.small]<-1/Pvals[is.small]/pi
         cct.stat<-colMeans(Pvals)
    }else{
         Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
         Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
         cct.stat<-colSums(Pvals)
    }
    #### return the ACAT p value(s).
    pval<-pcauchy(cct.stat,lower.tail = F)
    return(pval)
}



for(j in 1:length(files)) {
    
    tryCatch({

    file = files[j]
    load(file)
    
    name = gsub(".RData","",file)
    final.res[j,1:32] = CARE.pruning[[1]][1,]
    final.res[j,33] = name
    
    colnames(final.res) = c(colnames(CARE.pruning[[1]]),"Trait")
    
    final.res2[j,1] = name
    
    tmp = CARE.pruning[[1]]
    
    beta = mean(tmp[,"AIC.tilde_theta"],na.rm=TRUE)
    sd = sqrt(mean(tmp[,"AIC.Efron_se"]^2, na.rm=TRUE))
    
    #tmp = tmp[!is.na(tmp[,"AIC.Efron_p"]),]
    tmp = na.omit(tmp)
    tmp = as.matrix(tmp)
    final.res2[j,2] = beta
    final.res2[j,3] = sd
    final.res2[j,4] = ACAT(tmp[,"MA_BIC.Efron_p"])
    final.res2[j,5] = mean(tmp[,"IV.AIC"])
    final.res2[j,6] = ACAT(tmp[,"BIC.Efron_p"])
    final.res2[j,7] = ACAT(tmp[,"BIC.p"])
    final.res2[j,8] = ACAT(tmp[,"BIC.hat_p"])
    final.res2[j,9] = ACAT(tmp[,"MA_AIC.Efron_p"])

    final.res2[j,10] = ACAT(tmp[,"AIC.Efron_p"])
    final.res2[j,11] = ACAT(tmp[,"AIC.p"])
    final.res2[j,12:16] = tmp[,"AIC.Efron_p"]
    }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
    })
}

final.res2 = final.res2[!is.na(final.res2[,6]),]
#dim(final.res2)
#colSums(final.res2<0.0001,na.rm=T)
#colSums(final.res2<0.01,na.rm=T)
#colSums(final.res2<0.05,na.rm=T)


#final.res2[final.res2[,6]<0.001,]
#colSums(final.res<0.001,na.rm=T)


# ebi-a-GCST006250; Intelligence
# ebi-a-GCST006572; Cognitive performance

# remove the results from the following:
# ieu-a-814, ieu-a-1058; ieu-a-85; ieu-a-86; ieu-a-90; ieu-a-91

indx = final.res[,"Trait"]
#indx = grepl("ieu-a-814",indx)  | grepl("ieu-a-90",indx) | grepl("ieu-a-85",indx) | grepl("ieu-a-86",indx) |grepl("ieu-a-91",indx)|grepl("ieu-a-92",indx) | grepl("ieu-a-1239",indx) | grepl("ieu-a-836",indx) | grepl("ebi-a-GCST006098",indx) | grepl("ebi-a-GCST006097",indx)

indx =  grepl("ieu-a-1001",indx) | grepl("ieu-a-996",indx) | grepl("ieu-a-836",indx)  | grepl("ieu-a-302",indx) | grepl("ukb-d-1747_6",indx) | grepl("ieu-a-60",indx)
#small number of snps
#ieu-a-298, ieu-a-44, ieu-a-45, ieu-a-801, 

indx = !indx

#ieu-a-814 Ischaemic stroke due to too small sample size: |grepl("ieu-a-1058",indx)

#ieu-a-1058 Celiac disease # the number of SNPs is too small, delete this

#ebi-a-GCST006098, ebi-a-GCST006097: all based on UKB; delete them


#final.res = final.res[indx,]
CARE.res = final.res2[indx,]
#CARE.res = final.res2[indx,]
CARE.res = CARE.res[!is.na(CARE.res[,1]),]

#check smallest p-values
CARE.res[CARE.res<1e-6] = 1e-6
pvals =  CARE.res[,6]
orders = order(pvals)
files2 = files[!(grepl("ieu-a-1001",files) | grepl("ieu-a-996",files) | grepl("ieu-a-836",files) | grepl("ieu-a-302",files) | grepl("ukb-d-1747_6",files))]
files2[orders]
-log10(pvals[orders][1])
tmp = CARE.res[,1]
tmp = gsub("_.*","",tmp)


#quantile(final.res2[,4])
# other competing methods:
MR.res = as.data.frame(matrix(NA,length(files),34))

colnames(MR.res) = c("MRcML_MA_BIC_beta","MRcML_MA_BIC_se","MRcML_MA_BIC_p","DP_MA_BIC_beta","DP_MA_BIC_se","DP_MA_BIC_p","ComMix_beta","ConMix_se","ConMix_p","MRLasso_beta","MRLasso_se","MRLasso_p","IVW_beta","IVW_se","IVW_p","Egger_beta","Egger_se","Egger_p","Weighted_median_beta","Weighted_median_se","Weighted_median_p","Weighted_mode_beta","Weighted_mode_se","Weighted_mode_p","RAPS_beta","RAPS_se","RAPS_p","MRmix_beta","MRmix_se","MRmix_p","MRAPSS_beta","MRAPSS_se","MRAPSS_p","Trait")
for(j in 1:length(files)) {
    
    file = files[j]
    load(file)
    name = gsub(".RData","",file)
    MR.res[j,34] = name
    MR = MR.pruning$res
    MRAPSS = MRAPSS.pruning$MRAPSS_res
    if(is.null(MR.pruning)) {
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
    MR.res[j,13:15] = c(tmp[1:2],pnorm(-abs(TSMR_result[,1] / TSMR_result[,2]))*2)
    
    
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
}


MR.res =  MR.res[!is.na( MR.res[,1]),]

#MR.res[MR.res[,"Trait"] == "ebi-a-GCST006250_ukb-d-1747_5",]

#indx = MR.res[,"Trait"]
#indx = grepl("ieu-a-814",indx) |grepl("ieu-a-1058",indx) | grepl("ieu-a-85",indx) | grepl("ieu-a-86",indx) | grepl("ieu-a-90",indx) |grepl("ieu-a-91",indx)|grepl("ieu-a-92",indx) | grepl("ieu-a-1239",indx) | grepl("ieu-a-836",indx) | grepl("ebi-a-GCST006098",indx) | grepl("ebi-a-GCST006097",indx)
#indx =  grepl("ieu-a-1239",indx) | grepl("ieu-a-1058",indx) | grepl("ebi-a-GCST006250",indx)  | grepl("ebi-a-GCST006572",indx)

#indx = !indx
indx = MR.res[,"Trait"]
#indx = grepl("ieu-a-814",indx)  | grepl("ieu-a-90",indx) | grepl("ieu-a-85",indx) | grepl("ieu-a-86",indx) |grepl("ieu-a-91",indx)|grepl("ieu-a-92",indx) | grepl("ieu-a-1239",indx) | grepl("ieu-a-836",indx) | grepl("ebi-a-GCST006098",indx) | grepl("ebi-a-GCST006097",indx)

indx =  grepl("ieu-a-1001",indx) | grepl("ieu-a-996",indx) | grepl("ieu-a-836",indx)  | grepl("ieu-a-302",indx) | grepl("ukb-d-1747_6",indx) | grepl("ieu-a-60",indx)
indx = !indx
MR.res = MR.res[indx,]
quantile(MR.res[,"MRmix_p"])
quantile(MR.res[,"Egger_p"])
quantile(MR.res[,"IVW_p"])
quantile(MR.res[,"MRcML_MA_BIC_p"])
quantile(MR.res[,"DP_MA_BIC_p"])

MR.res[MR.res[,"DP_MA_BIC_p"]<0.001,]
MR.res[MR.res[,"MRmix_p"]<0.001,]


CARE.res[CARE.res<1e-6] = 1e-6
MR.res[MR.res<1e-6] = 1e-6

#Efron_p
plot.res = list("CARE" = CARE.res[,6])

jpeg("sfig1a_pruning_efron_p_6.jpg",width=6,height=5,units="in",res=300)
qqunif.plot(plot.res, auto.key = list(corner = c(.95, .05)))
dev.off()


plot.res = list("cML-DP" = MR.res[,"DP_MA_BIC_p"],"IVW" = MR.res[,"IVW_p"],"MR-APSS" = as.numeric(MR.res[,'MRAPSS_p']))
jpeg("sfig1b.jpg", width=6,height=5,units="in",res=300)
qqunif.plot(plot.res, auto.key = list(corner = c(.95, .05)), thin.obs.places = 10,thin.exp.places = 10)
dev.off()

plot.res = list("MR-mix" = MR.res[,"MRmix_p"],"MR-Egger" = MR.res[,"Egger_p"],
                "RAPS" = MR.res[,"RAPS_p"], "ContMix" = MR.res[,"ConMix_p"]
                )
jpeg("sfig1c.jpg",width=6,height=5,units="in",res=300)
qqunif.plot(plot.res, auto.key = list(corner = c(.95, .05)))
dev.off()

plot.res = list("cML" = MR.res[,"MRcML_MA_BIC_p"],"Weighted-Median"= MR.res[,"Weighted_median_p"],
                "Weighted-Mode"= MR.res[,"Weighted_mode_p"],
                "MR-Lasso" = MR.res[!is.na(MR.res[,"MRLasso_p"]),"MRLasso_p"])
jpeg("sfig1d.jpg",width=6,height=5,units="in",res=300)
qqunif.plot(plot.res,auto.key=list(corner=c(.95,.05)))
dev.off()



###CARE with clumping
final.res2 = as.data.frame(matrix(NA,length(files),16))
for(j in 1:length(files)) {
  
  tryCatch({
    
    file = files[j]
    load(file)
    
    name = gsub(".RData","",file)
    final.res[j,1:32] = CARE.clumping[[1]][1,]
    final.res[j,33] = name
    
    colnames(final.res) = c(colnames(CARE.clumping[[1]]),"Trait")
    
    final.res2[j,1] = name
    
    tmp = CARE.clumping[[1]]
    
    beta = mean(tmp[,"AIC.tilde_theta"],na.rm=TRUE)
    sd = sqrt(mean(tmp[,"AIC.Efron_se"]^2, na.rm=TRUE))
    
    #tmp = tmp[!is.na(tmp[,"AIC.Efron_p"]),]
    tmp = na.omit(tmp)
    tmp = as.matrix(tmp)
    final.res2[j,2] = beta
    final.res2[j,3] = sd
    final.res2[j,4] = ACAT(tmp[,"MA_BIC.Efron_p"])
    final.res2[j,5] = mean(tmp[,"IV.AIC"])
    final.res2[j,6] = ACAT(tmp[,"BIC.Efron_p"])
    final.res2[j,7] = ACAT(tmp[,"BIC.p"])
    final.res2[j,8] = ACAT(tmp[,"BIC.hat_p"])
    final.res2[j,9] = ACAT(tmp[,"MA_AIC.Efron_p"])
    
    final.res2[j,10] = ACAT(tmp[,"AIC.Efron_p"])
    final.res2[j,11] = ACAT(tmp[,"AIC.p"])
    final.res2[j,12:16] = tmp[,"AIC.Efron_p"]
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

final.res2 = final.res2[!is.na(final.res2[,6]),]
indx = final.res[,"Trait"]
#indx = grepl("ieu-a-814",indx)  | grepl("ieu-a-90",indx) | grepl("ieu-a-85",indx) | grepl("ieu-a-86",indx) |grepl("ieu-a-91",indx)|grepl("ieu-a-92",indx) | grepl("ieu-a-1239",indx) | grepl("ieu-a-836",indx) | grepl("ebi-a-GCST006098",indx) | grepl("ebi-a-GCST006097",indx)
indx =  grepl("ieu-a-1001",indx) | grepl("ieu-a-996",indx) | grepl("ieu-a-836",indx)  | grepl("ieu-a-302",indx) | grepl("ukb-d-1747_6",indx) | grepl("ieu-a-60",indx)

indx = !indx

#ieu-a-814 Ischaemic stroke due to too small sample size: |grepl("ieu-a-1058",indx)

#ieu-a-1058 Celiac disease # the number of SNPs is too small, delete this

#ebi-a-GCST006098, ebi-a-GCST006097: all based on UKB; delete them


#final.res = final.res[indx,]
CARE.res = final.res2[indx,]
#CARE.res = final.res2[indx,]
CARE.res = CARE.res[!is.na(CARE.res[,1]),]



plot.res = list("CARE" = CARE.res[,6])

jpeg("supfig_NegativeControl_b.jpg",width=6,height=5,units="in",res=300)
qqunif.plot(plot.res, auto.key = list(corner = c(.95, .05)))
dev.off()


### Fixed-effect IVW
files2 = files[!(grepl("ieu-a-1001",files) | grepl("ieu-a-996",files) | grepl("ieu-a-836",files) | grepl("ieu-a-302",files) | grepl("ieu-a-60",files) | grepl("ukb-d-1747_6",files))]
#les
fe.IVW = list()
for (j in 1:length(files2)){
  file = files2[j]
  load(file)
  tmp = MR.pruning$res
  p = as.numeric(tmp$MR_IVW_fe$pval)
  fe.IVW = c(fe.IVW,p)
}
fe.IVW = unlist(fe.IVW)
fe.IVW[fe.IVW<1e-6] = 1e-6

plot.res = list("Fixed-effect IVW" = fe.IVW)
jpeg("supfig_NegativeControl_a.jpg",width=6,height=5,units="in",res=300)
qqunif.plot(plot.res, auto.key = list(corner = c(.95, .05)))
dev.off()