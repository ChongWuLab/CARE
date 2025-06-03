require(nleqslv)

# rerandomization idea
rerandomization_bias <- function(gamma.exp,gamma.out,se.exp,se.out,pthr,etamean = 0.5,random_seeds = 999) {
    
    
    set.seed(random_seeds)
    
    ###########################################
    # Step 1: Select significant SNPs:      ###
    ###########################################
    Z = gamma.exp / se.exp
    C_sel = qnorm(pthr/2,lower.tail = FALSE)
    
    W = rnorm(length(gamma.exp), 0, etamean) # generate random noise to make the selection and estimation independent
    
    ind_filter = which(abs(gamma.exp/se.exp + W) >= C_sel & abs(gamma.exp/se.exp) >= C_sel - etamean) # we do not want to use SNPs that are selected just purely by chance
    
    if(length(ind_filter) <= 2) { #stop when only two IVs
        warning("Less than three IVs have been selected. We do not recommend using CARE in such scenarios.")
        return(NULL)
    }
    
    numIV_sel = length(ind_filter)
    gamma.exp_sel = gamma.exp[ind_filter]
    gamma.out_sel = gamma.out[ind_filter]
    se.exp_sel = se.exp[ind_filter]
    se.out_sel = se.out[ind_filter]
        
    
    #######################################################################
    # Step 2. Construct the unbiased carved estimator (also the UMVUE)  ###
    #######################################################################
    alpha1 = (-C_sel - gamma.exp_sel/se.exp_sel) / etamean
    alpha2 = (C_sel - gamma.exp_sel/se.exp_sel) / etamean
    gamma.carve = gamma.exp_sel - (se.exp_sel/etamean) * ( (dnorm(alpha2) - dnorm(alpha1)) / (pnorm(alpha1) + 1 - pnorm(alpha2)) )
    
    sigma2.carve = (1 + etamean^2 ) * se.exp_sel^2 # use the upper bound to estimate variance
    
    gamma.exp_sel = gamma.carve
    se.exp_sel = sqrt(sigma2.carve)
    
    ## save the output
    out = list(gamma.exp = gamma.exp_sel,gamma.out = gamma.out_sel, se.exp = se.exp_sel, se.out = se.out_sel,ind_filter = ind_filter)
    
    return(out)
}



#CARE2: based on the measurement error idea (no overlapped samples)
CARE2_boot <- function(gamma.exp, gamma.out, se.exp, se.out, nx,ny, pthr = 1e-5,nrep = 1000, random_seeds = 0, correct.method = "rerand", etamean = 0.5,random_start = 0,ind_filter = NULL, betax,betay,validIV,invalidIV)
{
    
    # select IV and correct winner curse bias
    if(correct.method == "rerand") {
        dat = rerandomization_bias(gamma.exp,gamma.out,se.exp,se.out,pthr,etamean = etamean,random_seeds = random_seeds)
        
        gamma.exp_sel = dat$gamma.exp
        gamma.out_sel = dat$gamma.out
        se.exp_sel = dat$se.exp
        se.out_sel = dat$se.out
        ind_filter = dat$ind_filter
        
    } else if (correct.method == "gold") {
        gamma.exp_sel = gamma.exp[ind_filter]
        gamma.out_sel = gamma.out[ind_filter]
        se.exp_sel = se.exp[ind_filter]
        se.out_sel = se.out[ind_filter]
    }
    
    corrected.mean = mean(abs(gamma.exp_sel - betax[ind_filter]))
    corrected.var = var(gamma.exp_sel - betax[ind_filter])
    
    corrected.metric = c(corrected.mean,corrected.var)
    names(corrected.metric) = c("Bias","Bias_var")
    
    allIV = paste0("IV",ind_filter)
    

    nIV = length(ind_filter)

    
    F = sum(gamma.exp[ind_filter]^2 / se.exp[ind_filter]^2 )/length(gamma.exp[ind_filter]) - 1
    varX = sum(betax[ind_filter]^2)
    varY = sum(betay[ind_filter]^2)
    
    sim.setting = c(nIV,F,varX,varY) # the statistics for this simulation setting
    names(sim.setting) = c("nIV","F","varX","varY")
    
    
    se2.out_sel = se.out_sel^2
    se2.exp_sel = se.exp_sel^2
    
    ##############################################################
    ### Bootstrap to consider post-selection bias
    ##############################################################
    
    wAll = rmultinom(nrep, size = nIV, prob = rep(1/nIV,nIV))
    
    Kvec = colSums(wAll!=0)
    KvecList = list()
    for(j in 1:dim(wAll)[2]) {
        tmp = floor(seq(0,Kvec[j] - 2,length.out = 50))
        tmp = unique(tmp)
        KvecList[[j]] = tmp
    }
    
    MRcML_result = mr_cMLC_CD_boot_efron(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = KvecList,wAll = wAll, random_start = random_start, maxit = 1000,  n = min(nx,ny), nrep = nrep)
    
    wALLmean = rowSums(wAll) / nrep
    
    wAll2 = wAll - matrix(rep(wALLmean,nrep),dim(wAll)[1],dim(wAll)[2],byrow = FALSE)
    
    t = MRcML_result[,1]
    meant = mean(t)
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    MA_se_efron = sqrt(sum(wAll3^2))
    
    
    t = MRcML_result[,2]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    BIC_se_efron = sqrt(sum(wAll3^2))
    
    
    t = MRcML_result[,3]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    MA_AIC_se_efron = sqrt(sum(wAll3^2))
    
    t = MRcML_result[,4]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    AIC_se_efron = sqrt(sum(wAll3^2))
    
    
    MA_theta = mean(MRcML_result[,1])
    MA_se = sqrt(var(MRcML_result[,1]))
    MA_p = pnorm(abs(MA_theta/MA_se), lower.tail = FALSE) * 2
    
    BIC_theta = mean(MRcML_result[,2])
    BIC_se = sqrt(var(MRcML_result[,2]))
    BIC_p = pnorm(abs(BIC_theta/BIC_se), lower.tail = FALSE) * 2
    
    MA_AIC_theta = mean(MRcML_result[,3])
    MA_AIC_se = sqrt(var(MRcML_result[,3]))
    MA_AIC_p = pnorm(abs(MA_AIC_theta/MA_AIC_se), lower.tail = FALSE) * 2
    
    AIC_theta = mean(MRcML_result[,4])
    AIC_se = sqrt(var(MRcML_result[,4]))
    AIC_p = pnorm(abs(AIC_theta/AIC_se), lower.tail = FALSE) * 2
    
    
    w = rep(1,length(gamma.exp_sel))
    se2.exp_sel = se.exp_sel^2
    se2.out_sel = se.out_sel^2
    
    
    tmp = floor(seq(0,length(gamma.exp_sel) - 2,length.out = 50))
    K_vec = unique(tmp)
    
    mrCMLraps = mr_cMLC_CD3(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = K_vec, w= w, random_start = random_start, maxit = 1000,  n = min(nx,ny))
    
    theta2 = mrCMLraps[[1]]
    
    MA_p2 = pnorm(abs(theta2[1]/MA_se), lower.tail = FALSE) * 2
    BIC_p2 = pnorm(abs(theta2[2]/BIC_se), lower.tail = FALSE) * 2
    
    MA_AIC_p2 = pnorm(abs(theta2[3]/MA_AIC_se), lower.tail = FALSE) * 2
    AIC_p2 = pnorm(abs(theta2[4]/AIC_se), lower.tail = FALSE) * 2
    
    
    MA_p3 = pnorm(abs(MA_theta/MA_se_efron), lower.tail = FALSE) * 2
    BIC_p3 = pnorm(abs(BIC_theta/BIC_se_efron), lower.tail = FALSE) * 2
    
    MA_AIC_p3 = pnorm(abs(MA_AIC_theta/MA_AIC_se_efron), lower.tail = FALSE) * 2
    AIC_p3 = pnorm(abs(AIC_theta/AIC_se_efron), lower.tail = FALSE) * 2
    
    
    MA_BIC = c(MA_theta,MA_se,MA_p,theta2[1],MA_p2,MA_se_efron,MA_p3)
    names(MA_BIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    BIC = c(BIC_theta,BIC_se,BIC_p,theta2[2],BIC_p2,BIC_se_efron,BIC_p3)
    names(BIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    MA_AIC = c(MA_AIC_theta,MA_AIC_se,MA_AIC_p,theta2[3],MA_AIC_p2,MA_AIC_se_efron,MA_AIC_p3)
    names(MA_AIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    AIC = c(AIC_theta,AIC_se,AIC_p,theta2[4],AIC_p2,AIC_se_efron,AIC_p3)
    names(AIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    ## Number of invalid IVs
    BIC_IV = mrCMLraps[[2]]
    BIC_IV = allIV[BIC_IV==0]
    
    AIC_IV = mrCMLraps[[3]]
    AIC_IV = allIV[AIC_IV==0]
    
    BIC_valid = sum(BIC_IV %in% validIV)
    AIC_valid = sum(AIC_IV %in% validIV)
    
    BIC_invalid = sum(BIC_IV %in% invalidIV)
    AIC_invalid = sum(AIC_IV %in% invalidIV)
    
    IV = c(length(BIC_IV),BIC_valid,BIC_invalid,length(AIC_IV),AIC_valid,AIC_invalid)
    names(IV) = c("BIC","BIC_valid","BIC_invalid","AIC","AIC_valid","AIC_invalid")
    
    
    out = list(MA_BIC=MA_BIC, BIC=BIC, MA_AIC = MA_AIC, AIC = AIC, setting = sim.setting, IV = IV,corrected.metric= corrected.metric)
    
    return(out)
}


#CARE2: based on the measurement error idea with overlapped samples
CARE2_boot_overlap <- function(gamma.exp, gamma.out, se.exp, se.out, nx,ny, overlap.mat, pthr = 1e-5,nrep = 1000, random_seeds = 0, correct.method = "rerand", etamean = 0.5, random_start = 0,ind_filter = NULL, betax,betay,validIV,invalidIV)
{
    
    if(correct.method == "rerand") {
        dat = rerandomization_bias(gamma.exp,gamma.out,se.exp,se.out,pthr,etamean = etamean,random_seeds = random_seeds)
        
        gamma.exp_sel = dat$gamma.exp
        gamma.out_sel = dat$gamma.out
        se.exp_sel = dat$se.exp
        se.out_sel = dat$se.out
        ind_filter = dat$ind_filter
        
    } else if (correct.method == "gold") {
        gamma.exp_sel = gamma.exp[ind_filter]
        gamma.out_sel = gamma.out[ind_filter]
        se.exp_sel = se.exp[ind_filter]
        se.out_sel = se.out[ind_filter]
    }
    
    allIV = paste0("IV",ind_filter)
    
    # generate w
    nIV = length(ind_filter)

    F = sum(gamma.exp[ind_filter]^2 / se.exp[ind_filter]^2 )/length(gamma.exp[ind_filter]) - 1
    varX = sum(betax[ind_filter]^2)
    varY = sum(betay[ind_filter]^2)
    
    sim.setting = c(nIV,F,varX,varY)
    names(sim.setting) = c("nIV","F","varX","varY")
    
    
    c1 = overlap.mat[1,1]
    c2 = overlap.mat[2,2]
    rho = overlap.mat[1,2]
    
    se2.out_sel = c2 * se.out_sel^2
    se2.exp_sel = c1 * se.exp_sel^2
    
    se.out_sel = sqrt(se2.out_sel)
    se.exp_sel = sqrt(se2.exp_sel)
    
    wAll = rmultinom(nrep, size = nIV, prob = rep(1/nIV,nIV))
    
    colSums(wAll!=0)
    
    
    Kvec = colSums(wAll!=0)
    KvecList = list()
    for(j in 1:dim(wAll)[2]) {
        tmp = floor(seq(0,Kvec[j] - 2,length.out = 50))
        tmp = unique(tmp)
        KvecList[[j]] = tmp
    }
    
    MRcML_result = cMLCOverlap_CD_boot_efron(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = KvecList,wAll = wAll, random_start = random_start, maxit = 1000,  n = min(nx,ny), nrep = nrep, se_exp = se.exp_sel, se_out = se.out_sel, rho = rho) # coordinate descent type algorithm
    
    wALLmean = rowSums(wAll) / nrep
    
    wAll2 = wAll - matrix(rep(wALLmean,nrep),dim(wAll)[1],dim(wAll)[2],byrow = FALSE)
    
    t = MRcML_result[,1]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    MA_se_efron = sqrt(sum(wAll3^2))
    
    t = MRcML_result[,2]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    BIC_se_efron = sqrt(sum(wAll3^2))
    
    
    t = MRcML_result[,3]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    MA_AIC_se_efron = sqrt(sum(wAll3^2))
    
    t = MRcML_result[,4]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    AIC_se_efron = sqrt(sum(wAll3^2))
    
    
    MA_theta = mean(MRcML_result[,1])
    MA_se = sqrt(var(MRcML_result[,1]))
    MA_p = pnorm(abs(MA_theta/MA_se), lower.tail = FALSE) * 2
    
    BIC_theta = mean(MRcML_result[,2])
    BIC_se = sqrt(var(MRcML_result[,2]))
    BIC_p = pnorm(abs(BIC_theta/BIC_se), lower.tail = FALSE) * 2
    
    MA_AIC_theta = mean(MRcML_result[,3])
    MA_AIC_se = sqrt(var(MRcML_result[,3]))
    MA_AIC_p = pnorm(abs(MA_AIC_theta/MA_AIC_se), lower.tail = FALSE) * 2
    
    AIC_theta = mean(MRcML_result[,4])
    AIC_se = sqrt(var(MRcML_result[,4]))
    AIC_p = pnorm(abs(AIC_theta/AIC_se), lower.tail = FALSE) * 2
    
    
    w = rep(1,length(gamma.exp_sel))
    se2.exp_sel = se.exp_sel^2
    se2.out_sel = se.out_sel^2
    
    
    tmp = floor(seq(0,length(gamma.exp_sel) - 2,length.out = 50))
    K_vec = unique(tmp)
    
    mrCMLraps = cMLCOverlap_CD3(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = K_vec, w= w, random_start = random_start, maxit = 1000,  n = min(nx,ny), se_exp = se.exp_sel, se_out = se.out_sel, rho = rho)
    
    theta2 = mrCMLraps[[1]]
    
    MA_p2 = pnorm(abs(theta2[1]/MA_se), lower.tail = FALSE) * 2
    BIC_p2 = pnorm(abs(theta2[2]/BIC_se), lower.tail = FALSE) * 2
    
    MA_AIC_p2 = pnorm(abs(theta2[3]/MA_AIC_se), lower.tail = FALSE) * 2
    AIC_p2 = pnorm(abs(theta2[4]/AIC_se), lower.tail = FALSE) * 2
    
    
    MA_p3 = pnorm(abs(MA_theta/MA_se_efron), lower.tail = FALSE) * 2
    BIC_p3 = pnorm(abs(BIC_theta/BIC_se_efron), lower.tail = FALSE) * 2
    
    MA_AIC_p3 = pnorm(abs(MA_AIC_theta/MA_AIC_se_efron), lower.tail = FALSE) * 2
    AIC_p3 = pnorm(abs(AIC_theta/AIC_se_efron), lower.tail = FALSE) * 2
    
    
    MA_BIC = c(MA_theta,MA_se,MA_p,theta2[1],MA_p2,MA_se_efron,MA_p3)
    names(MA_BIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    BIC = c(BIC_theta,BIC_se,BIC_p,theta2[2],BIC_p2,BIC_se_efron,BIC_p3)
    names(BIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    MA_AIC = c(MA_AIC_theta,MA_AIC_se,MA_AIC_p,theta2[3],MA_AIC_p2,MA_AIC_se_efron,MA_AIC_p3)
    names(MA_AIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    AIC = c(AIC_theta,AIC_se,AIC_p,theta2[4],AIC_p2,AIC_se_efron,AIC_p3)
    names(AIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    ## Number of invalid IVs
    BIC_IV = mrCMLraps[[2]]
    BIC_IV = allIV[BIC_IV==0]
    
    AIC_IV = mrCMLraps[[3]]
    AIC_IV = allIV[AIC_IV==0]
    
    BIC_valid = sum(BIC_IV %in% validIV)
    AIC_valid = sum(AIC_IV %in% validIV)
    
    BIC_invalid = sum(BIC_IV %in% invalidIV)
    AIC_invalid = sum(AIC_IV %in% invalidIV)
    
    IV = c(length(BIC_IV),BIC_valid,BIC_invalid,length(AIC_IV),AIC_valid,AIC_invalid)
    names(IV) = c("BIC","BIC_valid","BIC_invalid","AIC","AIC_valid","AIC_invalid")
    
    out = list(MA_BIC=MA_BIC, BIC=BIC, MA_AIC = MA_AIC, AIC = AIC, setting = sim.setting, IV = IV)
    
    return(out)
}



#CARE2 with Lasso algorithms: based on the measurement error idea (no overlapped samples)
CARE2_boot_newL1 <- function(gamma.exp, gamma.out, se.exp, se.out, nx,ny, pthr = 1e-5,nrep = 1000, random_seeds = 0, correct.method = "rerand", etamean = 0.5,random_start = 0,ind_filter = NULL, betax,betay,validIV,invalidIV, algorithm = "L0")
{
    
    # select IV and correct winner curse bias
    if(correct.method == "rerand") {
        dat = rerandomization_bias(gamma.exp,gamma.out,se.exp,se.out,pthr,etamean = etamean,random_seeds = random_seeds)
        
        gamma.exp_sel = dat$gamma.exp
        gamma.out_sel = dat$gamma.out
        se.exp_sel = dat$se.exp
        se.out_sel = dat$se.out
        ind_filter = dat$ind_filter
        
    } else if (correct.method == "gold") {
        gamma.exp_sel = gamma.exp[ind_filter]
        gamma.out_sel = gamma.out[ind_filter]
        se.exp_sel = se.exp[ind_filter]
        se.out_sel = se.out[ind_filter]
    }
    
    corrected.mean = mean(abs(gamma.exp_sel - betax[ind_filter]))
    corrected.var = var(gamma.exp_sel - betax[ind_filter])
    
    corrected.metric = c(corrected.mean,corrected.var)
    names(corrected.metric) = c("Bias","Bias_var")
    
    allIV = paste0("IV",ind_filter)
    

    nIV = length(ind_filter)

    
    F = sum(gamma.exp[ind_filter]^2 / se.exp[ind_filter]^2 )/length(gamma.exp[ind_filter]) - 1
    varX = sum(betax[ind_filter]^2)
    varY = sum(betay[ind_filter]^2)
    
    sim.setting = c(nIV,F,varX,varY) # the statistics for this simulation setting
    names(sim.setting) = c("nIV","F","varX","varY")
    
    
    se2.out_sel = se.out_sel^2
    se2.exp_sel = se.exp_sel^2
    
    ##############################################################
    ### Bootstrap to consider post-selection bias
    ##############################################################

    wAll = rmultinom(nrep, size = nIV, prob = rep(1/nIV,nIV))
    
    if(algorithm=="L0") {
        Kvec = colSums(wAll!=0)
        KvecList = list()
        for(j in 1:dim(wAll)[2]) {
            tmp = floor(seq(0,Kvec[j] - 2,length.out = 50))
            tmp = unique(tmp)
            KvecList[[j]] = tmp
        }
    
        MRcML_result = mr_cMLC_CD_boot_efron(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = KvecList,wAll = wAll, random_start = random_start, maxit = 1000,  n = min(nx,ny), nrep = nrep)
    } else if (algorithm == "L1") {
        Kvec = colSums(wAll!=0)
        #TODO: revise this
        KvecList = list()
        nlambda = 50
        for(j in 1:dim(wAll)[2]) {
            max_lambda = max( c(max(abs(gamma.out_sel)/se2.out_sel) , max(abs(gamma.out_sel - 2 * gamma.exp_sel)/se2.out_sel), max(abs(gamma.out_sel - 2 * gamma.exp_sel)/se2.out_sel )))
            lambda_seq <- exp(seq(log(max_lambda), log(max_lambda * 1e-4), length.out = nlambda))
            KvecList[[j]] = lambda_seq
        }
    
        MRcML_result = mr_CARE_L1_boot_efron(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = KvecList,wAll = wAll, random_start = random_start, maxit = 2000,  n = min(nx,ny), nrep = nrep)
    } else if (algorithm == "L1_algorithm2") {
        Kvec = colSums(wAll!=0)
        #TODO: revise this
        KvecList = list()
        nlambda = 50
        for(j in 1:dim(wAll)[2]) {
            max_lambda = max( c(max(abs(gamma.out_sel)/se2.out_sel) , max(abs(gamma.out_sel - 2 * gamma.exp_sel)/se2.out_sel), max(abs(gamma.out_sel - 2 * gamma.exp_sel)/se2.out_sel )))
            lambda_seq <- exp(seq(log(max_lambda), log(max_lambda * 1e-4), length.out = nlambda))
            KvecList[[j]] = lambda_seq
        }
    
        MRcML_result = mr_CARE_L2_boot_efron(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = KvecList,wAll = wAll, random_start = random_start, maxit = 2000,  n = min(nx,ny), nrep = nrep)
    }
    


        
    wALLmean = rowSums(wAll) / nrep
    
    wAll2 = wAll - matrix(rep(wALLmean,nrep),dim(wAll)[1],dim(wAll)[2],byrow = FALSE)
    
    t = MRcML_result[,1]
    meant = mean(t)
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    MA_se_efron = sqrt(sum(wAll3^2))
    
    
    t = MRcML_result[,2]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    BIC_se_efron = sqrt(sum(wAll3^2))
    
    
    t = MRcML_result[,3]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    MA_AIC_se_efron = sqrt(sum(wAll3^2))
    
    t = MRcML_result[,4]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    AIC_se_efron = sqrt(sum(wAll3^2))
    
    
    MA_theta = mean(MRcML_result[,1])
    MA_se = sqrt(var(MRcML_result[,1]))
    MA_p = pnorm(abs(MA_theta/MA_se), lower.tail = FALSE) * 2
    
    BIC_theta = mean(MRcML_result[,2])
    BIC_se = sqrt(var(MRcML_result[,2]))
    BIC_p = pnorm(abs(BIC_theta/BIC_se), lower.tail = FALSE) * 2
    
    MA_AIC_theta = mean(MRcML_result[,3])
    MA_AIC_se = sqrt(var(MRcML_result[,3]))
    MA_AIC_p = pnorm(abs(MA_AIC_theta/MA_AIC_se), lower.tail = FALSE) * 2
    
    AIC_theta = mean(MRcML_result[,4])
    AIC_se = sqrt(var(MRcML_result[,4]))
    AIC_p = pnorm(abs(AIC_theta/AIC_se), lower.tail = FALSE) * 2
    
    
    w = rep(1,length(gamma.exp_sel))
    se2.exp_sel = se.exp_sel^2
    se2.out_sel = se.out_sel^2
    
    if(algorithm=="L0") {
        tmp = floor(seq(0,length(gamma.exp_sel) - 2,length.out = 50))
        K_vec = unique(tmp)
        mrCMLraps = mr_cMLC_CD3(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = K_vec, w= w, random_start = random_start, maxit = 1000,  n = min(nx,ny))
    } else if (algorithm == "L1") {
        max_lambda = max( c(max(abs(gamma.out_sel)/se2.out_sel) , max(abs(gamma.out_sel - 2 * gamma.exp_sel)/se2.out_sel), max(abs(gamma.out_sel - 2 * gamma.exp_sel)/se2.out_sel )))
        lambda_seq <- exp(seq(log(max_lambda), log(max_lambda * 1e-4), length.out = nlambda))
        mrCMLraps = mr_L1_CD3(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = lambda_seq, w= w, random_start = random_start, maxit = 1000,  n = min(nx,ny))
    } else if (algorithm == "L1_algorithm2") {
        max_lambda = max( c(max(abs(gamma.out_sel)/se2.out_sel) , max(abs(gamma.out_sel - 2 * gamma.exp_sel)/se2.out_sel), max(abs(gamma.out_sel - 2 * gamma.exp_sel)/se2.out_sel )))
        lambda_seq <- exp(seq(log(max_lambda), log(max_lambda * 1e-4), length.out = nlambda))
        mrCMLraps = mr_L2_CD3(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = lambda_seq, w= w, random_start = random_start, maxit = 1000,  n = min(nx,ny))
    }

    theta2 = mrCMLraps[[1]]
    
    MA_p2 = pnorm(abs(theta2[1]/MA_se), lower.tail = FALSE) * 2
    BIC_p2 = pnorm(abs(theta2[2]/BIC_se), lower.tail = FALSE) * 2
    
    MA_AIC_p2 = pnorm(abs(theta2[3]/MA_AIC_se), lower.tail = FALSE) * 2
    AIC_p2 = pnorm(abs(theta2[4]/AIC_se), lower.tail = FALSE) * 2
    
    
    MA_p3 = pnorm(abs(MA_theta/MA_se_efron), lower.tail = FALSE) * 2
    BIC_p3 = pnorm(abs(BIC_theta/BIC_se_efron), lower.tail = FALSE) * 2
    
    MA_AIC_p3 = pnorm(abs(MA_AIC_theta/MA_AIC_se_efron), lower.tail = FALSE) * 2
    AIC_p3 = pnorm(abs(AIC_theta/AIC_se_efron), lower.tail = FALSE) * 2
    
    
    MA_BIC = c(MA_theta,MA_se,MA_p,theta2[1],MA_p2,MA_se_efron,MA_p3)
    names(MA_BIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    BIC = c(BIC_theta,BIC_se,BIC_p,theta2[2],BIC_p2,BIC_se_efron,BIC_p3)
    names(BIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    MA_AIC = c(MA_AIC_theta,MA_AIC_se,MA_AIC_p,theta2[3],MA_AIC_p2,MA_AIC_se_efron,MA_AIC_p3)
    names(MA_AIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    AIC = c(AIC_theta,AIC_se,AIC_p,theta2[4],AIC_p2,AIC_se_efron,AIC_p3)
    names(AIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    ## Number of invalid IVs
    BIC_IV = mrCMLraps[[2]]
    BIC_IV = allIV[BIC_IV==0]
    
    AIC_IV = mrCMLraps[[3]]
    AIC_IV = allIV[AIC_IV==0]
    
    BIC_valid = sum(BIC_IV %in% validIV)
    AIC_valid = sum(AIC_IV %in% validIV)
    
    BIC_invalid = sum(BIC_IV %in% invalidIV)
    AIC_invalid = sum(AIC_IV %in% invalidIV)
    
    IV = c(length(BIC_IV),BIC_valid,BIC_invalid,length(AIC_IV),AIC_valid,AIC_invalid)
    names(IV) = c("BIC","BIC_valid","BIC_invalid","AIC","AIC_valid","AIC_invalid")
    
    
    out = list(MA_BIC=MA_BIC, BIC=BIC, MA_AIC = MA_AIC, AIC = AIC, setting = sim.setting, IV = IV,corrected.metric= corrected.metric)
    
    return(out)
}