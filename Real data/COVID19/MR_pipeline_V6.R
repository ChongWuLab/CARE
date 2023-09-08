#!/usr/bin/env Rscript
#slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
#job.id2 <- as.numeric(slurm_arrayid)
args <- commandArgs(TRUE)
outcome.name <- as.character(args[[2]])
exposure.name <- as.character(args[[1]])

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
#source("step2_LDSC.R")
source("step3_CARE.R")
source("CARE_support.R")
source("step4_otherMR.R")

source("mr_raps_own.R")
source("mr_lasso_own.R")

library(Rcpp)
library(RcppArmadillo)

sourceCpp("CARE_support_measurement_overlap.cpp")
sourceCpp("CARE_support_measurement_with_CD2.cpp")
sourceCpp("CARE_support.cpp")

library(stringr)
#job.id  = 11
savedir = "/rsrch5/scratch/biostatistics/wzhang24/CARE/covid19/smmry_res/"

# finished jobs:
finsihed.jobs = list.files(savedir)

#COVID.dir = "/gpfs/research/chongwu/shared/summary_statistics/COVID19/release6/"
COVID.dir = "/rsrch5/scratch/biostatistics/wzhang24/CARE/summarystats/covid19/release7/"
#outcome.sumstatsAll = c("COVID19_HGI_A2_ALL_leave_23andme_20210607.b37.txt","COVID19_HGI_B1_ALL_leave_23andme_20210607.b37.txt","COVID19_HGI_B2_ALL_leave_23andme_20210607.b37.txt","COVID19_HGI_C2_ALL_leave_23andme_20210607.b37.txt")
#outcome.nameAll = c("A2_V6","B1_V6","B2_V6","C2_V6")

#for(job.id in (job.id2*2 -1):(job.id2 * 2)) {

#job.id = job.id2
#job.id = 12
# prepare the data
#exposure_id = readRDS("/gpfs/research/chongwu/Chong/CARE/NegativeControl/IEU_GWAS_info.rds")
ieu_gwas_info = readRDS("/rsrch5/scratch/biostatistics/wzhang24/CARE/negativecontrol/IEU_GWAS_info.rds")
smmry_dir = "/rsrch5/scratch/biostatistics/wzhang24/CARE/summarystats/"


outcome.sumstats = paste0(COVID.dir,outcome.name,".tsv.gz")
regex = '(?<=HGI_)[^_]+'
outcome.name = str_extract(outcome.name,regex)


exposure.sumstats = paste0(smmry_dir,"converted/",exposure.name,".txt")
exposure.n = subset(ieu_gwas_info,id==exposure.name)$sample_size



save.file = paste0(savedir,exposure.name,"_",outcome.name,".RData")

save.filejudge = paste0(exposure.name,"_",outcome.name,".RData")
if(save.filejudge %in% finsihed.jobs) {
    cat("This job has already been finished previously. \n")
} else {
    # create results directory:
    #exposure.name = "ieu-a-44"
    #outcome.name = "ukb-d-1747_4"
    workdir = paste0(getwd(),"/prepared_sumstats/",exposure.name,"_",outcome.name)
    
    dir.create(workdir)
    
    setwd(workdir)
    # prepare the GWAS data
    exist_files = list.files(workdir,pattern=".txt")
    if (paste0(exposure.name,"_GWAS.txt") %in% exist_files){
      print(paste0('Exposure ',exposure.name,' already exists. Continue next step.'))
    } else {
      read_GWAS(exposure.sumstats,exposure.name,exposure.n)
      }
    
    if (paste0(outcome.name,"_GWAS.txt") %in% exist_files){
      print(paste0('Outcome ',outcome.name,' already exists. Continue next step.'))
    } else {
      read_GWAS(sumstats = outcome.sumstats,outname = outcome.name,Ninput = NULL, A1 = "ALT", A2 = "REF",SNP = "rsid",CHR = "#CHR", BETA = "all_inv_var_meta_beta", SE = "all_inv_var_meta_sebeta", MAF = "all_meta_AF", Pval = "all_inv_var_meta_p", N = NULL, nColumns = 1:18, nControl = "all_inv_var_meta_controls", nCase = "all_inv_var_meta_cases")
      }
    # step 2: construct overlapping matrix
    #ldsc.res = construct_mat(exposure.name,outcome.name)
    
    #cov = ldsc.res[[1]]
    # CARE results:
    CARE.clumping = CARE_MRbase(exposure.name,outcome.name,indep.method='clumping',rerand=TRUE,biascorrect="rerand")
    CARE.pruning = CARE_MRbase(exposure.name,outcome.name,indep.method='prunning',rerand=TRUE,biascorrect="rerand")
    # competing methods results: (add results for competing methods)
    MR.clumping = MR_MRbase(exposure.name,outcome.name,indep.method='clumping',rerand=FALSE)
    MR.pruning = MR_MRbase(exposure.name,outcome.name,indep.method='prunning',rerand=FALSE)
    MRAPSS.clumping = MR_APSS(exposure.name,outcome.name,indep.method='clumping',rerand=FALSE)
    MRAPSS.pruning = MR_APSS(exposure.name,outcome.name,indep.method='prunning',rerand=FALSE)
    # save results
    save.file = paste0(savedir,exposure.name,"_",outcome.name,".RData")
    save(CARE.clumping,CARE.pruning, MR.clumping, MR.pruning, MRAPSS.clumping, MRAPSS.pruning, file = save.file)
}

#}


exposure.id = exposure.name
outcome.id = outcome.name
overlap.mat = cov
refpanel=NULL
etamean = 0.5
sel.pthr = 5e-5
proxies = "FALSE"
nRun = 5

