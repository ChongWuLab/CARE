require(TwoSampleMR)
require(data.table)
require(dplyr)

CARE_MRbase <- function(exposure.id, outcome.id, refpanel=NULL, etamean = 0.5, sel.pthr = 5e-5, proxies = "FALSE", nRun = 5,indep.method = "clumping", rerand = TRUE,biascorrect = "rerand")  {
    
    tmpdir = paste0(getwd(),"/tmp/")
    dir.create(tmpdir)
    
    if(is.null(refpanel)) {
        refpanel = "1000G.EUR.ALLSNP.QC.CHR"
        warning("This is just for local use")
    }
    
    # read it through gwasvcf
    dat = fread(paste(exposure.id,"_GWAS.txt",sep=""))
    dat = as.data.frame(dat)
    
    outcomedat = fread(paste(outcome.id,"_GWAS.txt",sep=""))
    outcomedat = as.data.frame(outcomedat)
    
    #only consider SNPs in outcome data
    dat = dat[dat$SNP %in% outcomedat$SNP,]
    
    
    res = matrix(NA, nRun,32)
    
    MR_datall = list()
    MR_IVall = list()
    RIVW = list()

    
    for(outj in 1:nRun) {
        
        set.seed(outj * 100)
        #####################################################
        # Step 1.0: get the independent SNPs
        #####################################################
        
        pruned = NULL
        
        cat("Repeat ",outj,"\n")
        cat("Finish reading exposure data: ", exposure.id,"\n")
        
        for(chr.id in 1:22) {
            tmp.snp.file = paste0(tmpdir,exposure.id,"_CHR",chr.id,".txt")
            
            tmp = dat %>% filter(CHR==chr.id)
            
            W = rnorm(dim(tmp)[1], 0, etamean)
            
            # Step 2: Select significant SNPs:
            C_sel = qnorm(sel.pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
            
            if(rerand) { # recalculated p value when rerandomized
                tmpZ2 = tmp[,"BETA"]/tmp[,"SE"] + W
                tmpZ = tmp[,"BETA"]/tmp[,"SE"]

                tmpP2 = pnorm(abs(tmpZ2/sqrt(1 + etamean^2)),lower.tail = FALSE) * 2
                ind_filter = which(abs(tmpZ) >= C_sel - etamean & abs(tmpZ2) >= C_sel)

            } else {
                # use the original ones
                tmpZ = tmp[,"BETA"]/tmp[,"SE"]
                tmpP2 = tmp[,"P"]
                ind_filter = which(abs(tmpZ) >= C_sel)

            }

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
            
            
            if(indep.method == "prunning.org") {
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
               
                if(indep.method == "prunning.org") {
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
        
        if(length(pruned) < 5) {
            next
        }
        
        exp_dat = dat %>% filter(SNP %in% pruned)
        exp_dat$P =  pnorm(abs(exp_dat$BETA/exp_dat$SE), lower.tail = F) * 2
        
        exp_dat = as.data.frame(exp_dat)

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
        outcome_dat = outcomedat[outcomedat[,"SNP"]%in% exp_dat$SNP, ,drop=FALSE]
        
        outcome_dat$id.outcome = outcome.id
        outcome_dat$outcome = outcome.id
        
        if(!"EAF" %in% colnames(outcome_dat)) {
            outcome_dat$EAF = NA
        }
        
        outcome_dat = outcome_dat[, c("id.outcome","outcome","BETA","SE","A1","A2","P","EAF","N","SNP","CHR","POS")]
        
        colnames(outcome_dat) = c("id.outcome","outcome","beta.outcome","se.outcome","effect_allele.outcome","other_allele.outcome","pval.outcome","eaf.outcome","samplesize.outcome","SNP","chr","pos")
        
        
        # double check the matching is correct
        
        harmonsze_dat = harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat)#
        harmonsze_dat = harmonsze_dat[harmonsze_dat[,"mr_keep"], ,drop=FALSE]
        
        MR_datall[[outj]] = harmonsze_dat
        
        if(dim(harmonsze_dat)[1] > 3) {
            
            ########################################
            # Step 3: conduct analysis
            ########################################
            gamma.exp_sel = harmonsze_dat$beta.exposure
            gamma.out_sel = harmonsze_dat$beta.outcome
            se.exp_sel = harmonsze_dat$se.exposure
            se.out_sel = harmonsze_dat$se.outcome
            nx = floor(median(harmonsze_dat$samplesize.exposure))
            ny = floor(median(harmonsze_dat$samplesize.outcome))
            
            
            ##################################
            # add results for RIVW; write down those into a function;
            ##################################
            RIVW.nopop = RIVW(gamma.exp_sel, gamma.out_sel, se.exp_sel, se.out_sel,etamean,sel.pthr)
            #RIVW.pop = RIVW(gamma.exp_sel, gamma.out_sel, se.exp_sel, se.out_sel,etamean,sel.pthr,overlap.mat)
            
            #RIVW.res = list(nopop = RIVW.nopop, pop = RIVW.pop)
            RIVW.res = list(nopop = RIVW.nopop)
            RIVW[[outj]] = RIVW.res
            
            ###################################
            CAREres = CARE2_boot(gamma.exp_sel, gamma.out_sel, se.exp_sel, se.out_sel, nx,ny, nrep = 5000, algorithm = "CD",random_start = 0,biascorrect = biascorrect,etamean = etamean, pthr=sel.pthr)
            
            tmpIV = list(AIC = CAREres[["AIC_IV"]], BIC = CAREres[["BIC_IV"]]) #, res2 =  CAREres[["res2"]],, Allres = CAREres[["MRcML_boot"]]
            MR_IVall[[outj]] = tmpIV
            
            # using CARE AIC
            tmpres = CAREres[[1]]
            thetahat = tmpres$AIC[1]
            
            tmpest = gamma.out_sel/gamma.exp_sel
            tmpest = tmpest[CAREres[["AIC_IV"]]==0] # only use valid IV
            tmpest = quantile(tmpest,c(0.1,0.9))
            
            
            if(thetahat>tmpest[1] & thetahat < tmpest[2]) {
                tmpres = unlist(tmpres)
                res[outj,] = tmpres
                colnames(res) = names(tmpres)

                break
            } #only save the results that are not exploding
        }
    }
    
    return(list(res = res, MRdat = MR_datall,MRIV = MR_IVall,RIVW = RIVW))
}

