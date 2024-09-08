require(mr.divw)
require(nleqslv)


library(MendelianRandomization)
library(glmnet)
library(mr.raps)
library(MRMix)
#library(MRPRESSO)

library(TwoSampleMR)
library(MRcML)
library(truncnorm)
library(MRAPSS)

library(mvtnorm)

#devtools::install_github("xue-hr/MRcML")
#devtools::install_github("gqi/MRMix")
#library(devtools)
#library(withr)
#withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/4.0/", install_github("xue-hr/MRcML"))
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")
#install_github("tye27/mr.divw")
#create output:
save_datdir = getwd()
save_datdir = paste(save_datdir,"/simRes12/",sep="")
system(paste("mkdir -p ",save_datdir,sep=""))


source("mr_raps_own.R")
source("mr_lasso_own.R")

library(Rcpp)
library(RcppArmadillo)

sourceCpp("CARE_support_measurement_with_CD2.cpp")
sourceCpp("CARE_support_measurement_overlap.cpp")
sourceCpp("cML_support.cpp")

source("CARE_support.R")
source("CARE_support2.R") #no correction, to check how bias correction helps

source("cML_support.R")

args <- commandArgs(TRUE)
indx1 <- as.numeric(args[[1]])
thetaU <- (as.character(args[[2]]))
Nin <- (as.character(args[[3]]))
PropInvalidIn <- (as.character(args[[4]]))

job.id <- as.numeric(args[[5]])


if(thetaU == "thetaU1") {
    indx2 = 1
} else {
    indx2 = 2
}

if(Nin == "N1") {
    indx3 = 1
} else if (Nin == "N2") {
    indx3 = 2
} else if (Nin == "N3") {
    indx3 = 3
} else if (Nin == "N5") {
    indx3 = 5
} else if (Nin == "N6") {
    indx3 = 6
}

if(PropInvalidIn == "Prop1") {
    indx4 = 1
} else if (PropInvalidIn == "Prop2") {
    indx4 = 2
}  else if (PropInvalidIn == "Prop3") {
    indx4 = 3
}

#thetavec = c(0.2,0.1, 0.05, 0, -0.05, -0.1, -0.2)
thetavec = c(0.1, 0.07, 0.03, 0, -0.03,-0.07, -0.1)
thetaUvec = c(0.3, 0.5)
Nvec = c(5e4, 8e4, 1e5, 1.5e5, 2.5e5, 5e5, 1e6) # 1:7
prop_invalid_vec = c(0.3, 0.5, 0.7)

#temp = as.integer(commandArgs(trailingOnly = TRUE))
temp = c(indx1,indx2,indx3,indx4)###

#temp = c(3,1,6,1)#
theta = thetavec[temp[1]] # True causal effect from X to Y
thetaU = thetaUx = thetaUvec[temp[2]] # Effect of the confounder on Y/X
N = Nvec[temp[3]] # Sample size for exposure X
prop_invalid = prop_invalid_vec[temp[4]] # Proportion of invalid IVs

pthr = 5e-8 # p-value threshold for instrument selection
pthr2 = 5e-5
NxNy_ratio = 1 # Ratio of sample sizes for X and Y
M = 2e5 # Total number of independent SNPs representing the common variants in the genome

# Model parameters for effect size distribution
pi1=0.02*(1-prop_invalid); pi3=0.01
pi2=0.02*prop_invalid;
sigma2x = sigma2y = 1e-5; sigma2u = 1e-5; sigma2x_td = sigma2y_td = 1e-5 - thetaU*thetaUx*sigma2u

print(paste("N", N, "pthr", pthr, "pi1", pi1, "theta", theta, "thetaU", thetaU, "prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio))

nx = N; ny = N/NxNy_ratio

care_sim_result_0.1 = list()
care_sim_result_0.3 = list()
care_sim_result_0.5 = list()
care_sim_result_0.7 = list()
care_sim_result_0.9 = list()

runningtime = list()
setting = list()
set.ind = job.id
simulation.ind.set = (((set.ind-1)*50+1):(set.ind*50))

outj = 1
#file = paste("simRes11/simres_MRMethods",set.ind,"N", N, "pthr", pthr,"pi1", pi1, "theta", theta, "thetaU", thetaU,"prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio,".Rdata",sep="")
#load(file)

for(sim.ind in simulation.ind.set)
{
    cat('simulation:',sim.ind,'\n')
    
    #sim.ind = 41
    numIV = 0
    numIV2 = 0
    tmpj = 0
    while(numIV <3 | numIV2 < 4) {
        
        set.seed(sim.ind + 10000 *tmpj)
        
        # Generate indices of causal SNPs
        ind1 = sample(M, round(M*pi1))
        causalsnps = ind1
        ind2 = sample(setdiff(1:M,causalsnps), round(M*pi2))
        causalsnps = c(causalsnps,ind2)
        ind3 = sample(setdiff(1:M,causalsnps), round(M*pi3))
        causalsnps = c(causalsnps,ind3)
        
        
        # Simulate effect size
        gamma = phi = alpha = rep(0,M)
        
        gamma[ind1] = rnorm(length(ind1), mean = 0, sd = sqrt(sigma2x))
        gamma[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2x_td))
        
        ind2all = ind2
        
        ind2 = ind2all[1:floor(length(ind2) * 0.5)]
        phi[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2u))
        
        alpha[ind2] = runif(length(ind2),0.01,0.03)

        ind2 = ind2all[!ind2all %in% ind2]
        alpha[ind2] = rnorm(length(ind2), mean = 0, sd = sqrt(sigma2y_td))
        
        ind2 = ind2all
        
        alpha[ind3] = rnorm(length(ind3), mean = 0, sd = sqrt(sigma2y)) #0.005
        
        # Generate summary statistics directly from summary-level model implied by individual-level model
        betax = gamma + thetaUx*phi
        betay = alpha + theta*betax + thetaU*phi
            
        
        betahat_x = betax + rnorm(M, mean = 0, sd = sqrt(1/nx))
        betahat_y = betay + rnorm(M, mean = 0, sd = sqrt(1/ny))
        betahat_xgold = betax + rnorm(M, mean = 0, sd = sqrt(1/nx))
        
        se_x = rep(sqrt(1/nx),M)
        se_y = rep(sqrt(1/ny),M)
        se_xgold = rep(sqrt(1/nx),M)
        
        ind_filter = which(2*pnorm(-sqrt(nx)*abs(betahat_x))<pthr)
        ind_filter_2 = which(2*pnorm(-sqrt(nx)*abs(betahat_x))<pthr2)
        
        numIV = length(ind_filter)
        numIV2 = length(ind_filter_2)
        weakIV = which(2*pnorm(-sqrt(nx)*abs(betahat_x))<5e-8 & 2*pnorm(-sqrt(nx)*abs(betahat_x))>5e-10)
        weakIV = length(weakIV)
        # calculate the statistics for this simulation setting:
        hertx = sum(betax^2)
        herty = sum(betay^2)
        
        # F statistics
        kappa = sum(betax^2/se_x^2)/length(betax)
        kappa.sel = sum(betax[ind_filter]^2/se_x[ind_filter]^2)/length(betax[ind_filter])
        
        F = sum(betahat_x[ind_filter]^2/se_x[ind_filter]^2)/length(betahat_x[ind_filter]) - 1
        varX = sum(betax[ind_filter]^2)
        varY = sum(betay[ind_filter]^2)
        
        sim.setting = c(numIV,weakIV,hertx,herty,kappa,kappa.sel,F,varX,varY)
        
        names(sim.setting) = c("nIV","nWeakIV","hertX","hertY","Kapp","KappSel","F","varX","varY")

        tmpj = tmpj + 1
        cat('numIV:', numIV,'numIV2:',numIV2,'\n')
        
    }
    
    sim.setting
    # save settings:
    setting = c(setting,list(setting = sim.setting))
    sum(ind_filter %in% ind1)
    sum(ind_filter %in% ind2all)


    sum(ind_filter %in% ind3)
    
    validIV = paste0("IV",ind1)
    invalidIV = paste0("IV",ind2)
    
    run.time = NULL
    nrep = 2000

    care.result0.1 = CARE2_boot(gamma.exp = betahat_x, gamma.out = betahat_y, se.exp = se_x, se.out = se_y, nx = nx,ny = ny,pthr = 5e-5, nrep = nrep, random_seeds = 0, correct.method = "rerand",etamean = 0.1,random_start = 0, ind_filter = ind_filter,betax,betay,validIV,invalidIV)
    care_sim_result_0.1[[outj]] = care.result0.1

    

    care.result0.3 = CARE2_boot(gamma.exp = betahat_x, gamma.out = betahat_y, se.exp = se_x, se.out = se_y, nx = nx,ny = ny,pthr = 5e-5, nrep = nrep, random_seeds = 0, correct.method = "rerand",etamean = 0.3,random_start = 0, ind_filter = ind_filter,betax,betay,validIV,invalidIV)
    care_sim_result_0.3[[outj]] = care.result0.3

    

    care.result0.5 = CARE2_boot(gamma.exp = betahat_x, gamma.out = betahat_y, se.exp = se_x, se.out = se_y, nx = nx,ny = ny,pthr = 5e-5, nrep = nrep, random_seeds = 0, correct.method = "rerand",etamean = 0.5,random_start = 0, ind_filter = ind_filter,betax,betay,validIV,invalidIV)
    care_sim_result_0.5[[outj]] = care.result0.5


    care.result0.7 = CARE2_boot(gamma.exp = betahat_x, gamma.out = betahat_y, se.exp = se_x, se.out = se_y, nx = nx,ny = ny,pthr = 5e-5, nrep = nrep, random_seeds = 0, correct.method = "rerand",etamean = 0.7,random_start = 0, ind_filter = ind_filter,betax,betay,validIV,invalidIV)
    care_sim_result_0.7[[outj]] = care.result0.7



    care.result0.9 = CARE2_boot(gamma.exp = betahat_x, gamma.out = betahat_y, se.exp = se_x, se.out = se_y, nx = nx,ny = ny,pthr = 5e-5, nrep = nrep, random_seeds = 0, correct.method = "rerand",etamean = 0.9,random_start = 0, ind_filter = ind_filter,betax,betay,validIV,invalidIV)
    care_sim_result_0.9[[outj]] = care.result0.9


    outj = outj + 1
    cat("Finish reran1 ",outj,"\n" )

}


save(list = c("care_sim_result_0.1",
"care_sim_result_0.3",
"care_sim_result_0.5",
"care_sim_result_0.7",
"care_sim_result_0.9"
),
file = paste(save_datdir,"simres_MRMethods",set.ind,"N", N, "pthr", pthr,
"pi1", pi1, "theta", theta, "thetaU", thetaU,
"prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio,".Rdata",sep=""))
