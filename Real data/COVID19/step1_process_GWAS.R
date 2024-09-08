require(data.table)


read_GWAS <- function(sumstats,outname, Ninput = NULL,A1 = NULL, A2 = NULL, SNP = NULL, CHR = NULL, POS = NULL, Zscore = NULL, BETA = NULL, SE = NULL, MAF = NULL, Pval = NULL, N = NULL, nColumns = NULL, nControl = NULL, nCase = NULL) {
  
    if (grepl(".txt",sumstats)) {
        raw = fread(sumstats)
        raw = as.data.frame(raw)
    } else {
        raw = fread(sumstats)
        raw = as.data.frame(raw)
        warning("We try to read GWAS summary data other txt format. The data may not read properly")
        
        cat("Following are the first few lines of the dataset:", sep = "\n")
        print(head(raw))
    }

    if(!is.null(nColumns)) {
        raw = raw[,nColumns]
    }
    
    header.inner = colnames(raw)
    # Initialize the header
    header.inner <- tolower(header.inner)

    # SNP
    if(is.null(SNP)) {
        try.snp <- c("snp", "markername", "snpid", "rs", "rsid", "rs_number", "snps")
    } else {
        try.snp <- tolower(SNP)
    }
    header.inner[header.inner %in% try.snp] <- "SNP"

    # A1
    if(is.null(A1)) {
        try.a1 <- c("a1", "allele1", "allele_1", "effect_allele", "reference_allele", "inc_allele", "ea", "alt", "a1lele1", "al1ele1")
    } else {
        try.a1 = tolower(A1)
    }
    header.inner[header.inner %in% try.a1] <- "A1"

    # A2
    if(is.null(A2)) {
        try.a2 <- c("a2", "allele2", "allele_2", "other_allele", "non_effect_allele", "dec_allele", "nea", "ref", "a0")
    } else {
        try.a2 = tolower(A2)
    }
    header.inner[header.inner %in% try.a2] <- "A2"
    
    # Z-score
    if(is.null(Zscore)) {
        try.z <- c("zscore", "z-score", "gc_zscore", "z")
    } else {
        try.z = tolower(Zscore)
    }
    header.inner[header.inner %in% try.z] <- "Z"
    
    if(is.null(CHR)) {
        try.chromosome <- c("chrom", "ch", "chr", "chromosome","#chr")
    } else {
        try.chromosome = tolower(CHR)
    }
    header.inner[header.inner %in% try.chromosome] <- "CHR"

    # P
    if(is.null(Pval)) {
        try.p <- c("pvalue", "p_value", "pval", "p_val", "gc_pvalue", "p","all_inv_var_meta_p")
    } else {
        try.p = tolower(Pval)
    }
    header.inner[header.inner %in% try.p] <- "P"

    # Beta
    if(is.null(BETA)) {
        try.beta <- c("b", "beta", "effects", "effect","all_inv_var_meta_beta")
    } else {
        try.beta = tolower(BETA)
    }
    header.inner[header.inner %in% try.beta] <- "BETA"

    # Odds ratio
    try.or <- c("or")
    header.inner[header.inner %in% try.or] <- "ODDS_RATIO"

    # Log odds
    try.logodds <- c("log_odds", "logor", "log_or")
    header.inner[header.inner %in% try.logodds] <- "LOG_ODDS"

    # Standard error
    if(is.null(SE)) {
        try.se <- c("se", "sebeta", "beta_se","all_inv_var_meta_sebeta")

    } else {
        try.se = tolower(SE)
    }
    header.inner[header.inner %in% try.se] <- "SE"

    
    # MAF
    if(is.null(MAF)) {
        try.maf <- c("eaf", "frq", "maf","frq_u","f_u","af1","all_meta_af")
    } else {
        try.maf = tolower(MAF)
    }
    header.inner[header.inner %in% try.maf] <- "EAF"

    # POS
    if(is.null(POS)) {
        try.pos <- c("pos","bp","position")
    } else {
        try.pos = tolower(POS)
    }
    header.inner[header.inner %in% try.pos] <- "POS"

    # sample size
    if(is.null(N)) {
        try.n <- c("n","nsample","nsum")

    } else {
        try.n = tolower(N)
    }
    header.inner[header.inner %in% try.n] <- "N"

    colnames(raw) <- header.inner
    
    if(sum("N" %in% colnames(raw))==0) { #sample size column does not exist
        if(!is.null(nControl) & !is.null(nCase)) {
            raw$N = as.numeric(raw[,nControl]) + as.numeric(raw[,nCase]) # the number of controls + the number of cases
        }
    }
    
    list.coerce <- c("Z", "BETA", "ODDS_RATIO", "LOG_ODDS", "SE", "AF1","N")

    for (i in 1:length(header.inner)) {
        if (header.inner[i] %in% list.coerce) {
            if (class(raw[, header.inner[i]]) != "numeric") {
                class(raw[, header.inner[i]]) <- "numeric"
                cat(paste0("Column ", header.inner[i], " has wrong class and has been coerced to numeric."), sep = "\n")
                cat("=============================================================================================================", sep = "\n")
            }
        }
    }

    # Missing z-score?
    calculate.z <- FALSE
    if (!("Z" %in% header.inner)) {
        warning("No Z score column, we calculate one based on the information we have")
        
        if ("BETA" %in% header.inner & "SE" %in% header.inner) {
            raw['Z'] <- raw$BETA / raw$SE
            calculate.z <- TRUE
        } else if ("ODDS_RATIO" %in% header.inner & "SE" %in% header.inner) {
            raw['Z'] <- log(raw$ODDS_RATIO) / raw$SE
            calculate.z <- TRUE
        } else if ("LOG_ODDS" %in% header.inner & "SE" %in% header.inner) {
            raw['Z'] <- raw$LOG_ODDS / raw$SE
            calculate.z <- TRUE
        } else if ("BETA" %in% header.inner & "P" %in% header.inner) {
            raw['Z'] <- sign(raw$BETA) * abs(qnorm(raw$P / 2))
            calculate.z <- TRUE
        } else if ("ODDS_RATIO" %in% header.inner & "P" %in% header.inner) {
            raw['Z'] <- sign(log(raw$ODDS_RATIO)) * abs(qnorm(raw$P / 2))
            calculate.z <- TRUE
        } else if ("LOG_ODDS" %in% header.inner & "P" %in% header.inner) {
            raw['Z'] <- sign(raw$ODDS_RATIO) * abs(qnorm(raw$P / 2))
            calculate.z <- TRUE
        } else {
            stop("I can't calculate z-score based on the information I have. SAD FACE EMOJI.", sep = "\n")
        }
        
        if (sum(is.na(raw$Z)) != 0) {
            n.start <- nrow(raw)
            
            raw <- raw[!is.na(raw$Z),]
            n.end <- nrow(raw)
            cat(paste0(n.start - n.end, " rows removed for having invalid z-score!"), sep = "\n")
        }
        cat("=============================================================================================================", sep = "\n")
    }
    
    if(!("P" %in% header.inner)) {
        raw['P'] <- pnorm(abs(raw[,"Z"]),lower.tail=F) * 2
        warning("We add P value column based on Z score")
    }

    if(sum(header.inner %in% c("SNP","A1","A2","CHR")) !=4) {
        stop("We tried our best to match the colnames of the summary data. Please revise the colnames to the following format. You can also report in GitHub to make the lists more comphrensive:\n SNP: snp, markername, snpid, rs, rsid, rs_number, snps\n A1: a1, allele1, allele_1, effect_allele, reference_allele, inc_allele, ea, ref, a1lele1, al1ele1\n A2: a2, allele2, allele_2, other_allele, non_effect_allele, dec_allele, nea, alt, a0\n Z: zscore, z-score, gc_zscore, z\n CHR: chrom, ch, chr, chromosome")
    }
    
    outcol = c("CHR","POS","SNP","A1","A2","Z","BETA","SE","P")
    if("EAF" %in% colnames(raw)) {
        outcol = c(outcol,"EAF")
    }
    
    if("N" %in% colnames(raw)) {
        outcol = c(outcol,"N")
    }
    
    sumstat.orgin = raw[,outcol]
    
    rm(raw)
    
    #sumstat.orgin = raw[,outcol]

    # Prepare the data
    sumstat.orgin$A1 <- toupper(sumstat.orgin$A1)
    sumstat.orgin$A2 <- toupper(sumstat.orgin$A2)

    a1 = sumstat.orgin[, "A1"]
    a2 = sumstat.orgin[, "A2"]

    keep = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))

    keep2 = nchar(a1)==1 &nchar(a2)==1
    keep = keep & keep2
    
    cat("Remove ", sum(!keep), " strand-ambiguous SNPs \n")

    sumstat.orgin = sumstat.orgin[keep,]
    
    outGWAS = paste(outname,"_GWAS.txt",sep="")

    if("EAF" %in% colnames(sumstat.orgin)) {
        frq = sumstat.orgin$EAF
        indx =!is.na(frq)
        cat("Remove ", sum(!indx), " SNPs with missing MAF\n")

        sumstat.orgin = sumstat.orgin[indx,]
        frq = sumstat.orgin$EAF

        indx = frq < 0 | frq > 1
        jj = sum(indx)
        if(jj > 0) {
            warning(jj,"SNPs had MAF outside of [0,1]. EAF column may be mislabeled.")
            sumstat.orgin = sumstat.orgin[!indx,]

        }
        frq = sumstat.orgin$EAF
        frq = cbind(frq,1-frq)
        frq = apply(frq,1,min)
        indx = frq >= 0.01
        sumstat.orgin = sumstat.orgin[indx,]
        cat("Remove ", sum(!indx), " SNPs with MAF<0.01\n")
    }
    
    if("N" %in% colnames(sumstat.orgin)) {
        
        nsample = sumstat.orgin$N
        nmin = quantile(nsample,0.9)/3
        indx = nsample >= nmin
        
        sumstat.orgin = sumstat.orgin[indx,]
        cat("Remove ", sum(!indx), " SNPs with sample size <=",floor(nmin),"\n")
    }
    
    indx = !duplicated(sumstat.orgin$SNP)
    sumstat.orgin = sumstat.orgin[indx,]
    cat("Remove ", sum(!indx), " SNPs with duplicated rs numbers.\n")
    
    indx = rowSums(is.na(sumstat.orgin))==0
    sumstat.orgin = sumstat.orgin[indx,]
    cat("Remove ", sum(!indx), " SNPs with missing values.\n")
    
    if(!is.null(Ninput) & !"N" %in% colnames(sumstat.orgin)) {
        sumstat.orgin$N = Ninput
        cat("Add a column for sample size with n = ", Ninput, "\n")
    }
    
    write.table(sumstat.orgin,outGWAS,row.names=FALSE,quote=FALSE)

    cat("Finish process GWAS summary data\n")

}
