require(readr)


construct_mat <- function(pheno1,pheno2,ldsc.refDir=NULL) {
    
    # add LDSC path to the working directory

    #ldsc.dir = "/gpfs/research/chongwu/shared/software/ldsc/ldsc"
    #command = paste0("PATH=$PATH:",ldsc.dir)
    #system(command)
    
    #pheno1 = "AD_Jansene"
    #pheno2 = "pheno69"
    
    #
    if(is.null(ldsc.refDir)) {
        ldsc.refDir ="/gpfs/research/chongwu/shared/software/ldsc/eur_w_ld_chr/"
        warning("This is just local use. Please change this with your own directory.")
    }
    
    
    
    #export PATH=$PATH:/gpfs/research/chongwu/shared/software/ldsc/ldsc
    
    # save the results for rg, hert, and its standard deviation
    ldsc.res = as.data.frame(matrix(NA,3,2))
    colnames(ldsc.res) = c("h2","intercept")
    rownames(ldsc.res) = c(pheno1,pheno2,"cross")
    
    cov = matrix(NA,2,2)
    colnames(cov) = rownames(cov) = c(pheno1,pheno2)
    # calculate the c11 for pheno 1
    
    idp1 = paste0(pheno1,"_GWAS.txt")
    idp1tmp = paste0(pheno1,"_hapmap_processed")


    command = paste0("munge_sumstats.py --sumstats ",idp1, " --out ",idp1tmp," --merge-alleles ",ldsc.refDir, "w_hm3.snplist --ignore Z")
    system(command)
    
    output = paste0(pheno2,"_",pheno1,"_ldsc")
    idp1tmp = paste0(idp1tmp,".sumstats.gz")
    
    command <- paste("ldsc.py --h2 ", idp1tmp, "  --ref-ld-chr ",ldsc.refDir," --w-ld-chr ",ldsc.refDir," --out ", output, sep = "")

    system(command)

    tmp <- read_file(paste(output, ".log", sep = ""))
    
    tmp2 <- gsub(".*Observed scale h2: ","",tmp)
    ldsc.res[1,1] = gsub("\nLambda.*", "", tmp2)

    tmp <- gsub("\n\nRatio.*", "", tmp)
    tmp <- gsub(".*\nIntercept: ", "", tmp)
    
    ldsc.res[1,2] = gsub("\nRatio:.*", "", tmp)
    
    tmp <- gsub(" \\(.*", "", tmp)
    cov[1,1] <- as.numeric(tmp)
    
    
    # calculate c22
    idp2 = paste0(pheno2,"_GWAS.txt")
    idp2tmp = paste0(pheno2,"_hapmap_processed")

    command = paste0("munge_sumstats.py --sumstats ",idp2, " --out ",idp2tmp," --merge-alleles ",ldsc.refDir, "w_hm3.snplist --ignore Z")
    system(command)
    
    idp2tmp = paste0(idp2tmp,".sumstats.gz")

    command <- paste("ldsc.py --h2 ", idp2tmp, "  --ref-ld-chr ",ldsc.refDir," --w-ld-chr ",ldsc.refDir," --out ", output, sep = "")

    system(command)

    tmp <- read_file(paste(output, ".log", sep = ""))

    tmp2 <- gsub(".*Observed scale h2: ","",tmp)
    ldsc.res[2,1] = gsub("\nLambda.*", "", tmp2)

    
    tmp <- gsub("\n\nRatio.*", "", tmp)
    tmp <- gsub(".*\nIntercept: ", "", tmp)
    
    ldsc.res[2,2] = gsub("\nRatio:.*", "", tmp)

    
    tmp <- gsub(" \\(.*", "", tmp)
    cov[2,2] <- as.numeric(tmp)
    
    # calculate the off diagonal c12 = c21
    command <- paste("ldsc.py --rg ", idp1tmp,",",idp2tmp, "  --ref-ld-chr ",ldsc.refDir," --w-ld-chr ",ldsc.refDir," --out ", output, sep = "")
    system(command)
    
    
    tmp <- read_file(paste(output, ".log", sep = ""))

    tmp2 <- gsub(".*Genetic Correlation: ","",tmp)
    ldsc.res[3,1] = gsub("\nZ-score.*", "", tmp2)

    
    tmp <- gsub("\n\nGenetic Correlation.*", "", tmp)
    tmp <- gsub(".*Total Observed scale gencov", "", tmp)
    tmp <- gsub(".*\nIntercept: ", "", tmp)
    
    ldsc.res[3,2] = tmp
    
    tmp <- gsub(" \\(.*", "", tmp)
    cov[1,2] <- as.numeric(tmp)
    cov[2,1] = cov[1,2]
    
    #make NA to 0;
    cov[is.na(cov)] = 0
    
    out = list(cov = cov,ldsc.res = ldsc.res)
    return(out)
}
### Test
#idp1 = "Kunkle_etal_Stage1_results.txt"
#idp1tmp = "tmpout"
#command = paste0("munge_sumstats.py --sumstats ",idp1, " --out ",idp1tmp," --merge-alleles ",ldsc.refDir, "w_hm3.snplist --N 63926") #Kunkle_etal_Stage1_results.txt
#system(command)

#idp1tmp = paste0(idp1tmp,".sumstats.gz")
#output = "tmpout2"
#command <- paste("ldsc.py --h2 ", idp1tmp, "  --ref-ld-chr ",ldsc.refDir," --w-ld-chr ",ldsc.refDir," --out ", output, sep = "")

#system(command)
## finish
