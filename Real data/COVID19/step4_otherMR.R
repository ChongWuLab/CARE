require(TwoSampleMR)
require(data.table)
require(dplyr)

MR_MRbase <- function(exposure.id, outcome.id, sel.pthr = 5e-8, refpanel=NULL,indep.method = "clumping", rerand = FALSE)  {
    
    tmpdir = paste0(getwd(),"/tmp2/")
    dir.create(tmpdir)
    
    #exposure.id = "pheno10"
    #outcome.id = "AD_Jansene"
    if(is.null(refpanel)) {
        refpanel = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.EUR.ALLSNP.QC.CHR"
        warning("This is just for local use")
    }
    
    etamean = 0.5
    #sel.pthr = 5e-5
    #proxies = "FALSE"
    
    # read it through gwasvcf
    dat = fread(paste(exposure.id,"_GWAS.txt",sep=""))
    dat = as.data.frame(dat)
    
    outcomedat = fread(paste(outcome.id,"_GWAS.txt",sep=""))
    outcomedat = as.data.frame(outcomedat)
    
    #only consider SNPs in outcome data
    dat = dat[dat$SNP %in% outcomedat$SNP,]
    
    
    MR_datall = list()
    
    
    set.seed(100)
    #####################################################
    # Step 1.0: get the indepdent SNPs
    #####################################################
    
    pruned = NULL
    
    cat("Finish reading exposure data: ", exposure.id,"\n")
    
    for(chr.id in 1:22) {
        tmp.snp.file = paste0(tmpdir,exposure.id,"_CHR",chr.id,".txt")
        
        tmp = dat %>% filter(CHR==chr.id)
        
        W = rnorm(dim(tmp)[1], 0, etamean)
        
        # Step 2: Select significant SNPs:
        C_sel = qnorm(sel.pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
        
        if(rerand) { # recalculated p value when rerandomized
            tmpP = tmp[,"BETA"]/tmp[,"SE"] + W
            tmpP2 = pnorm(abs(tmpP/sqrt(1 + etamean^2)),lower.tail = FALSE) * 2
        } else {
            # use the original ones
            tmpP = tmp[,"BETA"]/tmp[,"SE"]
            tmpP2 = tmp[,"P"]
        }
        
        
        ind_filter = which(abs(tmpP) >= C_sel)
        tmpSE = tmp[ind_filter,"SE"]

        if(length(ind_filter)==0 ) {
            next
        }
        
        # select the SNPs
        tmp = tmp[ind_filter,]
        tmp = as.data.frame(tmp)
        
        if(rerand) {
            if(indep.method == "clumping") {
                tmp$P = tmpP2[ind_filter]
            } else {
                tmp$P = 1/(1 + exp(-tmpSE)) - 0.5
            }
        }
        
        tmp.ssfile = paste0(tmpdir,exposure.id,"_ss_CHR",chr.id,".txt")
        
        write.table(tmp,file = tmp.ssfile,quote=FALSE,row.names=FALSE,col.names=TRUE)
        
        tmp = tmp %>% select(SNP)
        
        write.table(tmp,file = tmp.snp.file,quote=FALSE,row.names=FALSE,col.names=FALSE)
        
        system(paste0("plink --silent --bfile ",refpanel,chr.id,  " --chr ", chr.id," --extract ", tmp.snp.file," --make-bed --out ",tmpdir,exposure.id,"_CHR",chr.id))
        
        
        if(indep.method == "prunning") {
            command = paste0("plink --silent --bfile ",tmpdir,exposure.id,"_CHR",chr.id, "  --indep-pairwise 10000 100 0.001 --out ",tmpdir,exposure.id,"_CHR",chr.id)
            system(command)
            
            pruned.file = paste0(tmpdir,exposure.id,"_CHR",chr.id,".prune.in")
        } else {
            #clumping
            command <- paste("plink --silent --bfile ",tmpdir,exposure.id,"_CHR",chr.id, " --clump ", tmp.ssfile, "  --clump-p1 1 --clump-kb 10000 --clump-r2 0.001  --clump-snp-field SNP --clump-field P --out ", tmpdir,exposure.id,"_CHR",chr.id, sep = "")
            system(command)
            pruned.file = paste(tmpdir,exposure.id,"_CHR",chr.id, ".clumped", sep = "")
        }

        if(file.exists(pruned.file)) {
           
            if(indep.method == "prunning") {
                pruned.snp = fread(pruned.file,header=FALSE)
                pruned.snp = as.data.frame(pruned.snp)
                pruned.snp = pruned.snp[,1]
                pruned = c(pruned,pruned.snp)
                cat("CHR", chr.id, " pruned SNP:", length(pruned.snp),"\n")
            } else {
                pruned.snp = fread(pruned.file)
                pruned.snp = as.data.frame(pruned.snp)
                
                pruned = c(pruned,pruned.snp[,"SNP"])
                cat("CHR", chr.id, " pruned SNP:", dim(pruned.snp)[1],"\n")
            }
        }
        
    }
    
    exp_dat = dat %>% filter(SNP %in% pruned)
    
    if(dim(exp_dat)[1]<3) {
        return(NULL)
    } else {
        
        exp_dat$id = exposure.id

        if(!"EAF" %in% colnames(exp_dat)) {
            exp_dat$EAF = NA
        }
        
        exp_dat = exp_dat[,c("POS","CHR","SE","BETA","N","id","SNP","A1","A2","EAF")]
        
        exp_dat$exposure = exposure.id
        
        colnames(exp_dat) = c("pos.exposure","chr.exposure","se.exposure","beta.exposure","samplesize.exposure","id.exposure","SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","exposure")
        
        
        cat("Finish extracting exposure data\n")
        
        #####################################################
        # Step 2: combine with the outcome GWAS dataset
        #####################################################
        
        outcome_dat = outcomedat[outcomedat[,"SNP"]%in% exp_dat$SNP,]
        
        outcome_dat$id.outcome = outcome.id
        outcome_dat$outcome = outcome.id
        
        if(!"EAF" %in% colnames(outcome_dat)) {
            outcome_dat$EAF = NA
        }
        outcome_dat = outcome_dat[, c("id.outcome","outcome","BETA","SE","A1","A2","P","EAF","N","SNP","CHR","POS")]
        
        colnames(outcome_dat) = c("id.outcome","outcome","beta.outcome","se.outcome","effect_allele.outcome","other_allele.outcome","pval.outcome","eaf.outcome","samplesize.outcome","SNP","chr","pos")
        
        
        harmonsze_dat = harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat)#
        harmonsze_dat = harmonsze_dat[harmonsze_dat[,"mr_keep"],]
        
        MR_datall = harmonsze_dat
        
        ########################################
        # Step 3: conduct analysis
        ########################################
        b_exp = harmonsze_dat$beta.exposure
        b_out = harmonsze_dat$beta.outcome
        se_exp = harmonsze_dat$se.exposure
        se_out = harmonsze_dat$se.outcome
        nx = floor(median(harmonsze_dat$samplesize.exposure))
        ny = floor(median(harmonsze_dat$samplesize.outcome))

        run.time = NULL
        
        ### MRcML
        start.time = proc.time()[3]

        res_MRcML = mr_cML_DP(b_exp,
        b_out,
        se_exp,
        se_out,
        random_start = 0,
        random_start_pert = 0,
        n = min(nx,ny)
        )
        run.time = c(run.time,proc.time()[3] - start.time)

        
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

        
        
        run.time = c(run.time,proc.time()[3] - start.time)

        # perform TwoSampleMR
        start.time = proc.time()[3]

        set.seed(1)
        MR_IVW = TwoSampleMR::mr_ivw(b_exp = b_exp,b_out = b_out,
        se_exp = se_exp,se_out = se_out)
        run.time = c(run.time,proc.time()[3] - start.time)

        
        MR_IVW_fe = TwoSampleMR::mr_ivw_fe(b_exp = b_exp,b_out = b_out,
        se_exp = se_exp,se_out = se_out)
        
        MR_IVW_mre = TwoSampleMR::mr_ivw_mre(b_exp = b_exp,b_out = b_out,
        se_exp = se_exp,se_out = se_out)
        
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
        
        # perform MRMix

        
        start.time = proc.time()[3]

        set.seed(1)
        MRMix_data_std = standardize(betahat_x = b_exp,
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

        
        run.time = c(run.time,proc.time()[3] - start.time)

        #start.time = proc.time()[3]

        #MR-PRESSO
        #set.seed(1)
        #Data_Presso = data.frame(b_exp = b_exp,
        #b_out = b_out,
        #se_exp = se_exp,
        #se_out = se_out)
        
        #PRESSO_result = tryCatch(mr_presso(BetaOutcome = "b_out", BetaExposure = "b_exp", SdOutcome = "se_out", SdExposure = "se_exp", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = Data_Presso, NbDistribution = 500, SignifThreshold = 0.05),
        #error=function(e){cat("ERROR :", conditionMessage(e),"\n")})
        #MRPresso_sim_result = c(MRPresso_sim_result,
        #list(PRESSO_result = PRESSO_result))
        
        #run.time = c(run.time,proc.time()[3] - start.time)
        
        # run CARE:
        start.time = proc.time()[3]
        nocorrection = CARE2_boot(b_exp, b_out, se_exp, se_out, nx,ny,nrep = 5000, algorithm = "CD",biascorrect = "no",etamean = 0.5, pthr = sel.pthr)
        
        direct = CARE2_boot(b_exp, b_out, se_exp, se_out, nx,ny, nrep = 5000, algorithm = "CD",biascorrect = "direct",etamean = 0.5, pthr = sel.pthr)
        
        res_CARE = list(nocorrection=nocorrection,direct = direct)
        
        run.time = c(run.time,proc.time()[3] - start.time)
        
        names(run.time) = c("MRcML","ConMix","MR-Lasso","IVW","Egger","Median","Mode","RAPAS","MRmix","CARE")
                            
        res = list(res_MRcML = res_MRcML,ConMix_result = ConMix_result,MR_Lasso_result = MR_Lasso_result,MR_IVW = MR_IVW,MR_EGGER = MR_EGGER,MR_WEIGHTED_MEDIAN = MR_WEIGHTED_MEDIAN,MR_WEIGHTED_MODE = MR_WEIGHTED_MODE,MR_RAPS_para2 = MR_RAPS_para2, TSMR_T1toT2 = TSMR_T1toT2,MRMix_res = MRMix_res,res_CARE = res_CARE,MR_IVW_fe = MR_IVW_fe, MR_IVW_mre = MR_IVW_mre, run.time = run.time)
        return(list(res = res, MRdat = MR_datall))
    }
}

MR_APSS <- function(exposure.id, outcome.id,sel.pthr = 5e-5, refpanel=NULL,indep.method = "clumping", rerand = FALSE)  {
  
  tmpdir = paste0(getwd(),"/tmp/")
  dir.create(tmpdir)
  
  #exposure.id = "pheno10"
  #outcome.id = "AD_Jansene"
  if(is.null(refpanel)) {
    refpanel = "/rsrch5/scratch/biostatistics/wzhang24/1000G/1000G.EUR.ALLSNP.QC.CHR"
    warning("This is just for local use")
  }
  
  etamean = 0.5
  #sel.pthr = 5e-5
  #proxies = "FALSE"
  
  # read it through gwasvcf
  dat = fread(paste(exposure.id,"_GWAS.txt",sep=""))
  dat = as.data.frame(dat)
  
  outcomedat = fread(paste(outcome.id,"_GWAS.txt",sep=""))
  outcomedat = as.data.frame(outcomedat)
  
  #only consider SNPs in outcome data
  dat = dat[dat$SNP %in% outcomedat$SNP,]
  
  
  MR_datall = list()
  
  
  set.seed(100)
  #####################################################
  # Step 1.0: get the indepdent SNPs
  #####################################################
  
  pruned = NULL
  
  cat("Finish reading exposure data: ", exposure.id,"\n")
  
  for(chr.id in 1:22) {
    tmp.snp.file = paste0(tmpdir,exposure.id,"_CHR",chr.id,"txt")
    
    tmp = dat %>% filter(CHR==chr.id)
    
    W = rnorm(dim(tmp)[1], 0, etamean)
    
    # Step 2: Select significant SNPs:
    C_sel = qnorm(sel.pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
    
    if(rerand) { # recalculated p value when rerandomized
      tmpP = tmp[,"BETA"]/tmp[,"SE"] + W
      tmpP2 = pnorm(abs(tmpP/sqrt(1 + etamean^2)),lower.tail = FALSE) * 2
    } else {
      # use the original ones
      tmpP = tmp[,"BETA"]/tmp[,"SE"]
      tmpP2 = tmp[,"P"]
    }
    
    
    ind_filter = which(abs(tmpP) >= C_sel)
    tmpSE = tmp[ind_filter,"SE"]
    
    if(length(ind_filter)==0 ) {
      next
    }
    
    # select the SNPs
    tmp = tmp[ind_filter,]
    tmp = as.data.frame(tmp)
    
    if(rerand) {
      if(indep.method == "clumping") {
        tmp$P = tmpP2[ind_filter]
      } else {
        tmp$P = 1/(1 + exp(-tmpSE)) - 0.5
      }
    }
    
    tmp.ssfile = paste0(tmpdir,exposure.id,"_ss_CHR",chr.id,".txt")
    
    write.table(tmp,file = tmp.ssfile,quote=FALSE,row.names=FALSE,col.names=TRUE)
    
    tmp = tmp %>% select(SNP)
    
    write.table(tmp,file = tmp.snp.file,quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    system(paste0("plink --silent --bfile ",refpanel,chr.id,  " --chr ", chr.id," --extract ", tmp.snp.file," --make-bed --out ",tmpdir,exposure.id,"_CHR",chr.id))
    
    
    if(indep.method == "prunning") {
      command = paste0("plink --silent --bfile ",tmpdir,exposure.id,"_CHR",chr.id, "  --indep-pairwise 10000 100 0.001 --out ",tmpdir,exposure.id,"_CHR",chr.id)
      system(command)
      
      pruned.file = paste0(tmpdir,exposure.id,"_CHR",chr.id,".prune.in")
    } else {
      #clumping
      command <- paste("plink --silent --bfile ",tmpdir,exposure.id,"_CHR",chr.id, " --clump ", tmp.ssfile, "  --clump-p1 1 --clump-kb 10000 --clump-r2 0.001  --clump-snp-field SNP --clump-field P --out ", tmpdir,exposure.id,"_CHR",chr.id, sep = "")
      system(command)
      pruned.file = paste(tmpdir,exposure.id,"_CHR",chr.id, ".clumped", sep = "")
    }
    
    if(file.exists(pruned.file)) {
      
      if(indep.method == "prunning") {
        pruned.snp = fread(pruned.file,header=FALSE)
        pruned.snp = as.data.frame(pruned.snp)
        pruned.snp = pruned.snp[,1]
        pruned = c(pruned,pruned.snp)
        cat("CHR", chr.id, " pruned SNP:", length(pruned.snp),"\n")
      } else {
        pruned.snp = fread(pruned.file)
        pruned.snp = as.data.frame(pruned.snp)
        
        pruned = c(pruned,pruned.snp[,"SNP"])
        cat("CHR", chr.id, " pruned SNP:", dim(pruned.snp)[1],"\n")
      }
    }
    
  }
  
  exp_dat = dat %>% filter(SNP %in% pruned)
  
  if(dim(exp_dat)[1]<3) {
    return(NULL)
  } else {
    
    exp_dat$id = exposure.id
    
    if(!"EAF" %in% colnames(exp_dat)) {
      exp_dat$EAF = NA
    }
    
    exp_dat = exp_dat[,c("POS","CHR","SE","BETA","N","id","SNP","A1","A2","EAF")]
    
    exp_dat$exposure = exposure.id
    
    colnames(exp_dat) = c("pos.exposure","chr.exposure","se.exposure","beta.exposure","samplesize.exposure","id.exposure","SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","exposure")
    
    
    cat("Finish extracting exposure data\n")
    
    #####################################################
    # Step 2: combine with the outcome GWAS dataset
    #####################################################
    
    outcome_dat = outcomedat[outcomedat[,"SNP"]%in% exp_dat$SNP,]
    
    outcome_dat$id.outcome = outcome.id
    outcome_dat$outcome = outcome.id
    
    if(!"EAF" %in% colnames(outcome_dat)) {
      outcome_dat$EAF = NA
    }
    outcome_dat = outcome_dat[, c("id.outcome","outcome","BETA","SE","A1","A2","P","EAF","N","SNP","CHR","POS")]
    
    colnames(outcome_dat) = c("id.outcome","outcome","beta.outcome","se.outcome","effect_allele.outcome","other_allele.outcome","pval.outcome","eaf.outcome","samplesize.outcome","SNP","chr","pos")
    
    
    harmonsze_dat = harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat)#
    harmonsze_dat = harmonsze_dat[harmonsze_dat[,"mr_keep"],]
    
    MR_datall = harmonsze_dat
    
    ########################################
    # Step 3: conduct analysis
    ########################################
    b_exp = harmonsze_dat$beta.exposure
    b_out = harmonsze_dat$beta.outcome
    se_exp = harmonsze_dat$se.exposure
    se_out = harmonsze_dat$se.outcome
    nx = floor(median(harmonsze_dat$samplesize.exposure))
    ny = floor(median(harmonsze_dat$samplesize.outcome))
    
    run.time = NULL
    
    ### MRAPSS
    start.time = proc.time()[3]
    MRdata <- data.frame(b.exp = b_exp,
                         b.out = b_out,
                         se.exp = se_exp,
                         se.out = se_out,
                         Threshold = 1)
    C=diag(2)
    Omega=matrix(0,2,2)
    MRAPSS_res = MRAPSS(MRdata,
                        exposure="X",
                        outcome= "Y",
                        C = C,
                        Omega =  Omega ,
                        Cor.SelectionBias = T)
    run.time = c(run.time,proc.time()[3] - start.time)
    
    #names(run.time) = c("MRcML","ConMix","MR-Lasso","IVW","Egger","Median","Mode","RAPAS","MRmix","CARE")
    res = list(MRAPSS_res=MRAPSS_res, run.time=run.time)
    
    return(res)
  }
}
