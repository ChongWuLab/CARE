main.dir = getwd()

require(TwoSampleMR)
require(data.table)
require(dplyr)
require(nleqslv)
require(readr)

require(mr.divw)
require(nleqslv)
library(mvtnorm)


library(MendelianRandomization)
library(glmnet)
library(mr.raps)
library(MRMix)
library(MRAPSS)

library(TwoSampleMR)
library(MRcML)
library(truncnorm)

source("step1_process_GWAS.R")
source("step3_CARE.R")
source("step4_otherMR.R")

source("mr_raps_own.R")
source("mr_lasso_own.R")

library(Rcpp)
library(RcppArmadillo)

sourceCpp("CARE_support_measurement_overlap.cpp")
sourceCpp("CARE_support_measurement_with_CD2.cpp")
sourceCpp("CARE_support.cpp")
source("cML_support.R")
source("CARE_support.R")

#job.id  = 11
savedir = paste(main.dir,"/smmry_res/",sep="")

# finished jobs:
finished.jobs = list.files(savedir)



args <- commandArgs(TRUE)
outcome.name <- as.character(args[[1]])
exposure.name <- as.character(args[[2]])

#outcome.name = 'ukb-d-1747_6'
#exposure.name = 'ieu-a-996'

# prepare the data

    
save.file = paste0(savedir,exposure.name,"_",outcome.name,".RData")

save.filejudge = paste0(exposure.name,"_",outcome.name,".RData")
if(save.filejudge %in% finished.jobs) {
    cat("This job has already been finished previously. \n")
} else {
    # create results directory:
    #exposure.name = "ebi-a-GCST005920"
    #outcome.name = "ukb-d-1747_1"
    workdir = paste0(main.dir,"/prepared_sumstats/",exposure.name,"_",outcome.name)
    
    setwd(workdir)
    
    # CARE results:
    CARE.clumping = CARE_MRbase(exposure.name,outcome.name, refpanel=NULL, etamean = 0.5, sel.pthr = 5e-5, proxies = "FALSE", nRun = 5,indep.method = "clumping", rerand = TRUE,biascorrect = "rerand")
    print("finish CARE clumping")
    CARE.pruning.direct = CARE_MRbase(exposure.name,outcome.name, refpanel=NULL, etamean = 0.5, sel.pthr = 5e-5, proxies = "FALSE", nRun = 5,indep.method = "prunning", rerand = FALSE,biascorrect = "direct")
    print("finish CARE pruning direct")
    CARE.pruning = CARE_MRbase(exposure.name,outcome.name, refpanel=NULL, etamean = 0.5, sel.pthr = 5e-5, proxies = "FALSE", nRun = 5,indep.method = "prunning", rerand = TRUE,biascorrect = "rerand")
    print("finish CARE pruning")
    #CARE.pruning.org = CARE_MRbase(exposure.name,outcome.name, cov,refpanel=NULL, etamean = 0.5, sel.pthr = 5e-5, proxies = "FALSE", nRun = 3,indep.method = "prunning.org", rerand = TRUE,biascorrect = "rerand")

    # competing methods results: (add results for competing methods)
    MR.clumping = MR_MRbase(exposure.name,outcome.name, sel.pthr = 5e-8, refpanel=NULL,indep.method = "clumping", rerand = FALSE)
    print("finish MR clumping")
    MR.pruning = MR_MRbase(exposure.name,outcome.name, sel.pthr = 5e-8, refpanel=NULL, indep.method = "prunning", rerand = FALSE)
    print("finish MR pruning")
    #MRAPSS
    MRAPSS.clumping = MR_APSS(exposure.name,outcome.name,sel.pthr=5e-5, refpanel=NULL, indep.method = 'clumping', rerand = FALSE)
    print("finish MRAPSS clumping")
    MRAPSS.pruning = MR_APSS(exposure.name,outcome.name,sel.pthr=5e-5, refpanel=NULL, indep.method = 'prunning', rerand = FALSE)
    print("finish MRAPSS pruning")
    # save results
    save.file = paste0(savedir,exposure.name,"_",outcome.name,".RData")
    save(CARE.clumping, CARE.pruning.direct, CARE.pruning, MR.clumping, MR.pruning, MRAPSS.clumping, MRAPSS.pruning,file = save.file)
}
    

    



# exposure.id = exposure.name
# outcome.id = outcome.name
# overlap.mat = cov
# refpanel=NULL
# etamean = 0.5
# sel.pthr = 5e-5
# proxies = "FALSE"
# nRun = 5

