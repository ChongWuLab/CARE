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
#library(devtools)
#library(withr)
#withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/4.0/", install_github("xue-hr/MRcML"))

#create output:
save_datdir = getwd()
save_datdir = paste(save_datdir,"/simRes2/",sep="")
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
indx1 <- (eval(parse(text = args[[1]])))
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
} else if (Nin == "N4") {
    indx3 = 4
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
} else if (PropInvalidIn == "Prop4") {
    indx4 = 4
}

#job.id = 1; indx1 = 1; indx2 = 1; indx3 = 5; indx4 = 4
thetavec = c(0.1, 0.07, 0.03, 0, -0.03,-0.07, -0.1)
thetaUvec = c(0.3, 0.5)
Nvec = c(5000, 1e4, 5e4, 1e5, 5e5, 1e6)
prop_invalid_vec = c(0.3, 0.5, 0.7, 0.9)

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
sigma2x = 1e-5
sigma2y = 1e-5; sigma2u = 1e-5
sigma2x_td = 1e-5 - thetaU*thetaUx*sigma2u
sigma2y_td = 1e-5 - thetaU*thetaUx*sigma2u

print(paste("N", N, "pthr", pthr, "pi1", pi1, "theta", theta, "thetaU", thetaU, "prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio))

nx = N; ny = N/NxNy_ratio

care_sim_result = list()
care2_sim_result = list()
care3_sim_result = list()
CAREno_sim_result = list()

MRcML_sim_result = list()
ContMix_sim_result = list()
MRLasso_sim_Result = list()
TwoSampleMR_sim_result = NULL
MRMix_sim_result = list()
#MRPresso_sim_result = list()
MRAPSS_sim_result = list()

runningtime = list()
setting = list()
set.ind = job.id
simulation.ind.set = (((set.ind-1)*50+1):(set.ind*50))

outj = 1
#file = paste("simRes11/simres_MRMethods",set.ind,"N", N, "pthr", pthr,"pi1", pi1, "theta", theta, "thetaU", thetaU,"prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio,".Rdata",sep="")
#load(file)

for(sim.ind in simulation.ind.set)
{
    cat(sim.ind,"\n")
    
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
        
        #alpha[ind2] = rnorm(length(ind2), mean = 0.015, sd = sqrt(sigma2u))
        alpha[ind2] = rnorm(length(ind2), mean = 0.015, sd = sqrt(sigma2u)) #try mean = 0.02
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
    
    start.time = proc.time()[3]
    
    care.result = CARE2_boot(gamma.exp = betahat_x, gamma.out = betahat_y, se.exp = se_x, se.out = se_y, nx = nx,ny = ny,pthr = 5e-5, nrep = nrep, random_seeds = 0, correct.method = "rerand",etamean = 0.5,random_start = 0, ind_filter = ind_filter,betax,betay,validIV,invalidIV)
    care_sim_result[[outj]] = care.result

    run.time = c(run.time,proc.time()[3] - start.time)
    
    
    start.time = proc.time()[3]

    run.time = c(run.time,proc.time()[3] - start.time)
    
    start.time = proc.time()[3]


    run.time = c(run.time,proc.time()[3] - start.time)

    outj = outj + 1
    cat("Finish reran1 ",outj,"\n" )

    names(run.time) = c("CARE1","CARE2","CARE3")
    if (numIV > 2) # do not run other methods
    {
        b_exp = betahat_x[ind_filter]
        b_out = betahat_y[ind_filter]
        se_exp = se_x[ind_filter]
        se_out = se_y[ind_filter]
        
        ### MRcML
        start.time = proc.time()[3]

        res_MRcML = mr_cML_DP(b_exp,
        b_out,
        se_exp,
        se_out,
        random_start = 0,
        random_start_pert = 0,
        num_pert = 200,
        n = min(nx,ny)
        )
        run.time = c(run.time,proc.time()[3] - start.time)

        
        MRcML_sim_result = c(MRcML_sim_result,
        list(res_MRcML = res_MRcML))
        
        start.time = proc.time()[3]

        # Contamination Mixture
        set.seed(1)
        ConMix_result = tryCatch(mr_conmix(mr_input(bx = b_exp,
        bxse = se_exp,
        by = b_out,
        byse = se_out),
        CIMin = -0.5,
        CIMax = 0.5,
        CIStep = 0.0001),
        error=function(e){cat("ERROR :",
            conditionMessage(e),
        "\n")})
        ContMix_sim_result = c(ContMix_sim_result,
        list(ConMix_result = ConMix_result))
        run.time = c(run.time,proc.time()[3] - start.time)

        start.time = proc.time()[3]

        # MR-Lasso
        set.seed(1)
        MR_Lasso_result = tryCatch(mr_lasso_own(mr_input(bx = b_exp,
        bxse = se_exp,
        by = b_out,
        byse = se_out)
        ),
        error=function(e){cat("ERROR :",
            conditionMessage(e),
        "\n")})
        MRLasso_sim_Result = c(MRLasso_sim_Result,
        list(MR_Lasso_result = MR_Lasso_result))
        run.time = c(run.time,proc.time()[3] - start.time)

        # perform TwoSampleMR
        start.time = proc.time()[3]

        set.seed(1)
        MR_IVW = TwoSampleMR::mr_ivw(b_exp = b_exp,b_out = b_out,
        se_exp = se_exp,se_out = se_out)
        run.time = c(run.time,proc.time()[3] - start.time)

        start.time = proc.time()[3]

        set.seed(1)
        MR_EGGER = mr_egger_regression(b_exp = b_exp,b_out = b_out,
        se_exp = se_exp,se_out = se_out)
        run.time = c(run.time,proc.time()[3] - start.time)

        start.time = proc.time()[3]

        set.seed(1)
        MR_WEIGHTED_MEDIAN = mr_weighted_median(b_exp = b_exp,b_out = b_out,
        se_exp = se_exp,se_out = se_out)
        run.time = c(run.time,proc.time()[3] - start.time)

        start.time = proc.time()[3]

        set.seed(1)
        MR_WEIGHTED_MODE = mr_weighted_mode(b_exp = b_exp,b_out = b_out,
        se_exp = se_exp,se_out = se_out)
        ### mr raps para2
        run.time = c(run.time,proc.time()[3] - start.time)

        start.time = proc.time()[3]

        set.seed(1)
        MR_RAPS_para2 = NULL
        while (is.null(MR_RAPS_para2))
        {
            MR_RAPS_para2 = tryCatch(mr.raps.overdispersed.robust(b_exp, b_out, se_exp, se_out, loss.function = "huber", k = 1.345, initialization = c("l2"),suppress.warning = FALSE, niter = 20, tol = .Machine$double.eps^0.5),
            error=function(e){cat("ERROR :",
                conditionMessage(e),
            "\n")})
        }
        run.time = c(run.time,proc.time()[3] - start.time)

        
        TSMR_T1toT2 = c(MR_IVW$b, MR_IVW$se,
        MR_EGGER$b, MR_EGGER$se,
        MR_WEIGHTED_MEDIAN$b, MR_WEIGHTED_MEDIAN$se,
        MR_WEIGHTED_MODE$b, MR_WEIGHTED_MODE$se,
        MR_RAPS_para2$beta.hat,MR_RAPS_para2$beta.se)
        TwoSampleMR_sim_result = rbind(TwoSampleMR_sim_result,
        TSMR_T1toT2)
        # perform MRMix
        
        start.time = proc.time()[3]

        set.seed(1)
        MRMix_data_std = MRMix::standardize(betahat_x = b_exp,
        betahat_y = b_out,
        sx = se_exp,
        sy = se_out,
        xtype = "continuous",
        ytype = "continuous",
        nx = nx,
        ny = ny,
        MAF = NULL)
        MRMix_res = MRMix(MRMix_data_std$betahat_x_std,
        MRMix_data_std$betahat_y_std,
        MRMix_data_std$sx_std,
        MRMix_data_std$sy_std,
        profile = FALSE)
        MRMix_sim_result = c(MRMix_sim_result,
        list(MRMix_res = MRMix_res))
        run.time = c(run.time,proc.time()[3] - start.time)

        start.time = proc.time()[3]

        # MR-PRESSO
        #set.seed(1)
        #Data_Presso = data.frame(b_exp = b_exp,
        #b_out = b_out,
        #se_exp = se_exp,
        #se_out = se_out)
        
        #PRESSO_result = tryCatch(mr_presso(BetaOutcome = "b_out", BetaExposure = "b_exp", SdOutcome = "se_out", SdExposure = "se_exp", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = Data_Presso, NbDistribution = 500, SignifThreshold = 0.05),
        #error=function(e){cat("ERROR :", conditionMessage(e),"\n")})
        #MRPresso_sim_result = c(MRPresso_sim_result,
        #list(PRESSO_result = PRESSO_result))
        
        #CAREno
        overlap.mat = diag(2)
        CARE_res = CARE2_boot_overlap_nocorrect(b_exp, b_out, se_exp, se_out, nx,ny, overlap.mat, nrep = 2000)
        CAREno_sim_result = c(CAREno_sim_result, list(CARE_res = CARE_res))
        
        run.time = c(run.time,proc.time()[3] - start.time)

        #MRAPSS
        b_exp = betahat_x[ind_filter_2]
        b_out = betahat_y[ind_filter_2]
        se_exp = se_x[ind_filter_2]
        se_out = se_y[ind_filter_2]
        start.time = proc.time()[3]
        MRdata <- data.frame(b.exp = b_exp,
                             b.out = b_out,
                             se.exp = se_exp,
                             se.out = se_out,
                             Threshold = 1)
        C = diag(2)
        r_g = cor(betahat_x[!seq_along(betahat_x) %in% ind_filter_2],betahat_y[!seq_along(betahat_x) %in% ind_filter_2])
        rho_g = r_g * sqrt(hertx) * sqrt(herty)
        Omega = matrix(c(hertx,rho_g,rho_g,herty)/M,2,2)
        MRAPSS_res = MRAPSS(MRdata,
                            exposure="X",
                            outcome= "Y",
                            C = C,
                            Omega =  Omega ,
                            Cor.SelectionBias = T)
        MRAPSS_sim_result = c(MRAPSS_sim_result, list(MRAPSS_res = MRAPSS_res))
        run.time = c(run.time,proc.time()[3] - start.time)
        
        names(run.time) = c("CARE1","CARE2","CARE3","MRcML","ConMix","MR-Lasso","IVW","Egger","Median","Mode","RAPAS","MRmix","CAREno","MRAPSS")
        
        runningtime = c(runningtime,list(run.time = run.time))
        
        cat("Finish iteration ",sim.ind,"\n")
    } else {
        MRcML_sim_result = c(MRcML_sim_result,
        list(res_MRcML = NULL))
        ContMix_sim_result = c(ContMix_sim_result,
        list(ConMix_result = NULL))
        MRLasso_sim_Result = c(MRLasso_sim_Result,
        list(MR_Lasso_result = NULL))
        TwoSampleMR_sim_result = rbind(TwoSampleMR_sim_result,
        rep(0,10))
        MRMix_sim_result = c(MRMix_sim_result,
        list(MRMix_res = NULL))
        #MRPresso_sim_result = c(MRPresso_sim_result,
        #list(PRESSO_result = NULL))
        CAREno_sim_result = c(CAREno_sim_result,
                              list(CARE_res = NULL))
        MRAPSS_sim_result = c(MRAPSS_sim_result,
        list(MRAPSS_result = NULL))
    }
}


save(list = c("care_sim_result",
"care2_sim_result",
"care3_sim_result",
"CAREno_sim_result",
"MRcML_sim_result",
"ContMix_sim_result",
"MRLasso_sim_Result",
"TwoSampleMR_sim_result",
"MRMix_sim_result",
"MRAPSS_sim_result",
"runningtime",
"setting"),
file = paste(save_datdir,"simres_MRMethods",set.ind,"N", N, "pthr", pthr,
"pi1", pi1, "theta", theta, "thetaU", thetaU,
"prop_invalid", prop_invalid, "NxNy_ratio", NxNy_ratio,".Rdata",sep=""))
