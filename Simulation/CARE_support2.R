require(nleqslv)

#CARE2: based on the measurement error idea (This version is for real data pipeline)
CARE2_boot_overlap_nocorrect <- function(gamma.exp_sel, gamma.out_sel, se.exp_sel, se.out_sel, nx,ny, overlap.mat, nrep = 1000, random_start = 0)
{
    # generate w
    nIV = length(gamma.exp_sel)
    
    F = sum(gamma.exp_sel^2 / se.exp_sel^2 )/length(gamma.exp_sel) - 1
    
    
    sim.setting = c(nIV,F)
    names(sim.setting) = c("nIV","F")
    
    c1 = overlap.mat[1,1]
    c2 = overlap.mat[2,2]
    rho = overlap.mat[1,2]
    
    se2.out_sel = c2 * se.out_sel^2
    se2.exp_sel = c1 * se.exp_sel^2
    
    se.out_sel = sqrt(se2.out_sel)
    se.exp_sel = sqrt(se2.exp_sel)
    
    wAll = rmultinom(nrep, size = nIV, prob = rep(1/nIV,nIV))
    
    
    Kvec = colSums(wAll!=0)
    KvecList = list()
    for(j in 1:dim(wAll)[2]) {
        tmp = floor(seq(0,Kvec[j] - 2,length.out = 50))
        tmp = unique(tmp)
        KvecList[[j]] = tmp
    }
    
    MRcML_result = cMLCOverlap_CD_boot_efron(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = KvecList,wAll = wAll, random_start = random_start, maxit = 1000,  n = min(nx,ny), nrep = nrep, se_exp = se.exp_sel, se_out = se.out_sel, rho = rho) # coordinate descent type algorithm;
    
    
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
    BIC_IV = sum(BIC_IV==0)
    
    AIC_IV = mrCMLraps[[3]]
    AIC_IV = sum(AIC_IV==0)
    
    IV = c(BIC_IV,AIC_IV)
    names(IV) = c("BIC","AIC")
    
    res = list(MA_BIC=MA_BIC, BIC=BIC, MA_AIC = MA_AIC, AIC = AIC, setting = sim.setting, IV = IV)
    
    out = list(res = res, BIC_IV = mrCMLraps[[2]], AIC_IV = mrCMLraps[[3]])
    return(out)
}




