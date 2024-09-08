require(TwoSampleMR)
require(data.table)
require(dplyr)

dir = ""
setwd(dir)
files = list.files(dir)

files = files[grepl("B2_all_eur_V6",files) | grepl("C2_all_eur_V6",files)]
#combine the results
final.res = as.data.frame(matrix(NA,length(files),33))

final.res2 = as.data.frame(matrix(NA,length(files),15))

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
    sd = sqrt(mean(tmp[,"AIC.Efron_se"]^2,na.rm=TRUE))
    
    #tmp = tmp[!is.na(tmp[,"AIC.Efron_p"]),]
    tmp = na.omit(tmp)
    tmp = as.matrix(tmp)
    final.res2[j,2] = beta
    final.res2[j,3] = sd
    final.res2[j,4] = ACAT(tmp[,"MA_BIC.Efron_p"])
    final.res2[j,5] = mean(tmp[,"IV.AIC"])
    final.res2[j,6] = ACAT(tmp[,"BIC.Efron_p"])
    final.res2[j,7] = ACAT(tmp[,"BIC.p"])
    final.res2[j,8] = ACAT(tmp[,"MA_AIC.Efron_p"])

    final.res2[j,9] = ACAT(tmp[,"AIC.Efron_p"])
    final.res2[j,10] = ACAT(tmp[,"AIC.p"])
    final.res2[j,11] = ACAT(tmp[1:5,"AIC.Efron_p"])

    final.res2[j,12:15] = tmp[1:4,"AIC.Efron_p"]
    }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
    })
}

final.res2 = final.res2[!is.na(final.res2[,9]),]
dim(final.res2)
final.res2[final.res2[,9]<0.001,]
final.res2[final.res2[,6]<0.001,]
final.res2[final.res2[,8]<0.001,]

colSums(final.res2<0.001)

final.res = final.res[!is.na(final.res[,"AIC.Efron_p"]),]
colSums(final.res < 0.001)

final.res = final.res[final.res[,"AIC.Efron_p"]<0.001,]
final.res[final.res[,"AIC.Efron_p"]<0.001,]


final.res2 = final.res2[!is.na(final.res2[,1]),]
final.res2$exposure.id = gsub("_.*","",final.res2[,1])
final.res2$outcome.id = final.res2[,1]
final.res2$exposure.name = NA
final.res2$outcome.name = NA
final.res2$adjusted.p = p.adjust(final.res2[,6])

#tmp = final.res2[final.res2[,"adjusted.p"]<0.5,]


library(ieugwasr)
info = gwasinfo()
info = as.data.frame(info)

for(j in 1:dim(final.res2)[1]) {
    tryCatch({

    final.res2[j,"exposure.name"] = info[info[,1]==final.res2[j,"exposure.id"],"trait"]
    }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
    })
}

final.res2[final.res2[,10]<0.001,]


tmp = final.res2[,c("exposure.id","exposure.name")]
tmp = tmp[!duplicated(tmp[,1]),]
rownames(tmp) = tmp[,1]


remove.duplicated = c("ebi-a-GCST002222","ieu-a-300","ebi-a-GCST002223", "ieu-a-299", "ebi-a-GCST003116","ebi-a-GCST003769","ebi-a-GCST003770","ebi-a-GCST005327","ebi-a-GCST006686","ebi-a-GCST004599",paste0("ebi-a-GCST00",4600:4634),"ebi-a-GCST006907", "ebi-a-GCST006908","ebi-a-GCST006909","ebi-a-GCST006910",paste0("ebi-a-GCST0069",41:52),"ebi-a-GCST007090","ebi-a-GCST007091","ebi-a-GCST009722","ebi-a-GCST009760","ebi-a-GCST009761","ebi-a-GCST009762","ebi-a-GCST009763","ebi-a-GCST010776","ebi-a-GCST010777","ebi-a-GCST010778","ebi-a-GCST010779","ebi-a-GCST010780","ebi-a-GCST010783","ebi-a-GCST90000026","ebi-a-GCST90000027","ebi-a-GCST90000616","ebi-a-GCST90002412","ieu-a-1000","ebi-a-GCST005038","ukb-b-533","ukb-b-1747","ebi-a-GCST007557","ebi-a-GCST006802","ieu-a-2","ukb-b-19953","ebi-a-GCST005920","ebi-a-GCST005921","ebi-a-GCST005923","ebi-a-GCST006696","ebi-a-GCST006697","ebi-a-GCST006698","ebi-a-GCST006699","ebi-a-GCST006670","ebi-a-GCST006671","ebi-a-GCST006672","ebi-a-GCST006804","ieu-a-1090","ebi-a-GCST005840","ebi-a-GCST005841","ebi-a-GCST005842","ebi-a-GCST005843","ieu-a-1109","ebi-a-GCST006478","ieu-a-1007","ieu-a-814","ieu-a-965","ebi-a-GCST006368","ebi-a-GCST006099","ebi-a-GCST006100","ebi-a-GCST006061","ieu-a-1058","ieu-a-1239")

tmp =tmp[!tmp[,1] %in% remove.duplicated,]
dim(tmp)

remove.disease = c("ebi-a-GCST004988","ebi-a-GCST005195","ebi-a-GCST005232","ebi-a-GCST005838","ebi-a-GCST005902","ebi-a-GCST005903","ebi-a-GCST005904","ebi-a-GCST006061","ebi-a-GCST006085","ebi-a-GCST006414","ebi-a-GCST006464","ebi-a-GCST006867","ebi-a-GCST006906","ebi-a-GCST006940","ebi-a-GCST009541","ebi-a-GCST010723","ieu-a-1054","ieu-a-1102","ieu-a-1185","ieu-a-1186","ieu-a-12","ieu-a-22","ieu-a-31","ieu-a-44","ieu-a-798","ieu-a-833","ieu-a-966","ieu-a-967","ieu-a-970","ieu-a-975","ieu-a-996","ieu-b-18","ieu-b-2","ieu-b-7","ieu-b-7","ieu-a-90","ieu-a-91","ieu-a-92")
tmp =tmp[!tmp[,1] %in% remove.disease,]

B2 = final.res2[grepl("B2",final.res2[,1]),]
B2 = B2[B2$exposure.id %in% tmp[,1],]
B2$adjusted.p = p.adjust(B2[,6])
B2[B2$adjusted.p<0.1,]

B2[B2[,6]<0.05,]


final.res2[final.res2[,6]<0.05,]




C2 = final.res2[grepl("C2",final.res2[,1]),]
C2 = C2[C2$exposure.id %in% tmp[,1],]
C2$adjusted.p = p.adjust(B2[,6])
C2[C2$adjusted.p<0.05,]
C2[C2[,6]<0.05,]

library(sgof)
BH(B2[,6], alpha = 0.05)



library(ggplot2)




plot.dat = B2[B2[,6]<0.05,]

plot.dat = plot.dat[order(plot.dat[,6]),]

oddratio = exp(plot.dat[,2])

Pval = plot.dat[,6]
SE = abs(plot.dat[,2]) / qnorm(Pval/2,lower.tail=FALSE)

lower = exp(plot.dat[,2] - 1.96 * SE)
higher = exp(plot.dat[,2] + 1.96 * SE)

df <- data.frame(yAxis = dim(plot.dat)[1]:1,
                 boxOdds = oddratio,
                  boxCILow = lower,
                 boxCIHigh = higher,
Pval = Pval
                 )
boxLabels = plot.dat$exposure.name

p <- ggplot(df, aes(x = boxOdds, y = reorder(boxLabels,-Pval))) +
    geom_vline(aes(xintercept = 1), size = .5, linetype = "longdash") +
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height =
                    .2, color = "gray50") +
    geom_point(size = 3.5, color = "orange") +
    theme_bw() + scale_x_continuous(breaks = seq(0.6, 2, 0.2), labels = seq(0.6, 2, 0.2),
limits = c(0.65,1.9)) + theme(panel.grid.minor = element_blank(),text = element_text(size = 20)) +
    ylab("") +
    xlab("Odds ratio")




Pvalorder = Pval
plot.dat = C2[C2[,"exposure.name"] %in%boxLabels,]

plot.dat = plot.dat[order(plot.dat[,6]),]

oddratio = exp(plot.dat[,2])

Pval = plot.dat[,6]
SE = abs(plot.dat[,2]) / qnorm(Pval/2,lower.tail=FALSE)

lower = exp(plot.dat[,2] - 1.96 * SE)
higher = exp(plot.dat[,2] + 1.96 * SE)

df <- data.frame(yAxis = dim(plot.dat)[1]:1,
                 boxOdds = oddratio,
                  boxCILow = lower,
                 boxCIHigh = higher,
Pval = Pvalorder
                 )
boxLabels = plot.dat$exposure.name

p <- ggplot(df, aes(x = boxOdds, y = reorder(boxLabels,-Pval))) +
    geom_vline(aes(xintercept = 1), size = .5, linetype = "longdash") +
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height =
                    .2, color = "gray50") +
    geom_point(size = 3.5, color = "orange") +
    coord_trans(x = scales:::exp_trans(10))  +
    theme_bw()+
    theme(panel.grid.minor = element_blank(),text = element_text(size = 20)) +
    ylab("") +
    xlab("Odds ratio")






final.res2[final.res2[,6]<0.01,c(1:3,6)]

tmp = final.res2[grepl("B2_all_eur_V6",final.res2[,1]),]
tmp = as.data.frame(tmp)
tmp$adjusted = p.adjust(tmp[,6])
tmp[tmp$adjusted<0.2,]

final.res2[grepl("ebi-a-GCST010776",final.res2[,1]),]


indx = final.res[,"Trait"]

indx =  grepl("ieu-a-1239",indx) | grepl("ieu-a-1058",indx) #| grepl("ieu-a-1096",indx)  | grepl("ieu-a-814",indx) | grepl("ieu-a-836",indx)

indx = !indx


final.res = final.res[indx,]
CARE.res = final.res2[indx,]
CARE.res = CARE.res[!is.na(CARE.res[,1]),]

tmp = CARE.res[,1]
tmp = gsub("_.*","",tmp)


quantile(final.res2[,4])
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

MR.res[MR.res[,"IVW_p"]<0.001,]

indx = MR.res[,"Trait"]
#indx = grepl("ieu-a-814",indx) |grepl("ieu-a-1058",indx) | grepl("ieu-a-85",indx) | grepl("ieu-a-86",indx) | grepl("ieu-a-90",indx) |grepl("ieu-a-91",indx)|grepl("ieu-a-92",indx) | grepl("ieu-a-1239",indx) | grepl("ieu-a-836",indx) | grepl("ebi-a-GCST006098",indx) | grepl("ebi-a-GCST006097",indx)
indx =  grepl("ieu-a-1239",indx) | grepl("ieu-a-1058",indx) #| grepl("ieu-a-1096",indx)  | grepl("ieu-a-814",indx) | grepl("ieu-a-836",indx)

indx = !indx

MR.res = MR.res[indx,]
quantile(MR.res[,"MRmix_p"])
quantile(MR.res[,"Egger_p"])
quantile(MR.res[,"IVW_p"])
quantile(MR.res[,"MRcML_MA_BIC_p"])
quantile(MR.res[,"DP_MA_BIC_p"])

MR.res[MR.res[,"DP_MA_BIC_p"]<0.001,]
MR.res[MR.res[,"MRmix_p"]<0.001,]


CARE.res[CARE.res<1e-5] = 1e-5
MR.res[MR.res<1e-5] = 1e-5
plot.res = list("CARE" = CARE.res[,6])

pdf("sfig1a.pdf")
qqunif.plot(plot.res, auto.key = list(corner = c(.95, .05)))
dev.off()


plot.res = list("cML" = MR.res[,"MRcML_MA_BIC_p"],"cML-DP" = MR.res[,"DP_MA_BIC_p"],"IVW" = MR.res[,"IVW_p"],"Weighted-Median"= MR.res[,"Weighted_median_p"],"Weighted-Mode"= MR.res[,"Weighted_mode_p"])
pdf("sfig1b.pdf")
qqunif.plot(plot.res, auto.key = list(corner = c(.95, .05)), thin.obs.places = 10,thin.exp.places = 10)
dev.off()

plot.res = list("MR-mix" = MR.res[,"MRmix_p"],"MR-Egger" = MR.res[,"Egger_p"],"RAPS" = MR.res[,"RAPS_p"], "ContMix" = MR.res[,"ConMix_p"],"MR-Lasso" = MR.res[!is.na(MR.res[,"MRLasso_p"]),"MRLasso_p"])
pdf("sfig1c.pdf")
qqunif.plot(plot.res, auto.key = list(corner = c(.95, .05)))
dev.off()
