#slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
#job.id <- as.numeric(slurm_arrayid)
args <- commandArgs(TRUE)
job.id <- as.numeric(args[[1]])
# pr
library(devtools)
library(withr)
#withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/4.1", remotes::install_github("MRCIEU/TwoSampleMR"))

#withr::with_libpaths(new = "/gpfs/home/cwu3/R/x86_64-redhat-linux-gnu-library/4.1", remotes::install_github("mrcieu/gwasvcf"))

#library(TwoSampleMR)

library(data.table)
library(ieugwasr)
library(gwasvcf)

suppressWarnings(suppressPackageStartupMessages({
    library(gwasvcf)
    library(VariantAnnotation)
    library(dplyr)
    library(magrittr)
}))

library(mr.divw)

set_bcftools()

setwd("/rsrch5/scratch/biostatistics/wzhang24/CARE/summarystats/")

files = list.files()

files = files[grepl(".vcf.gz",files)]
files = files[!grepl(".vcf.gz.tbi",files)]

files2 = gsub(".vcf.gz","",files)
#existfiles = readRDS("/gpfs/research/chongwu/shared/summary_statistics/IEU_GWAS/converted/finished.rds")
existfiles = list.files("/rsrch5/scratch/biostatistics/wzhang24/CARE/summarystats/converted/")

#files = files[!paste0(files2,".txt") %in%existfiles]

#outcome.id = "ieu-b-110"
#job.id = 1
outcome.id = files[job.id]

out.id = gsub(".vcf.gz","",outcome.id)
#

if(paste0(out.id,".txt") %in% existfiles) {
    cat("The GWAS files have already been processed. Finish.\n")
} else {
    #outcome.id = "ebi-a-GCST007800.vcf.gz"
    vcf <- readVcf(outcome.id)
    outcome_dat = vcf_to_granges(vcf) %>% dplyr::as_tibble()

    outcome_dat = as.data.frame(outcome_dat)

    outputfile = paste0("/rsrch5/scratch/biostatistics/wzhang24/CARE/summarystats/converted/",out.id,".txt")

    out_dat = outcome_dat[,c("seqnames","start","ID","ALT","REF","ES","SE","AF","SS")]

    colnames(out_dat) = c("CHR","POS","SNP","A1","A2","BETA","SE","EAF","N")

    # if AF are missing ignore that
    if(sum(is.na(out_dat$EAF)) > 0.5 * dim(out_dat)[1]) {
        col.name = colnames(out_dat)
        col.name = col.name[!col.name %in% "EAF"]
        out_dat = out_dat[,col.name]
    }

    if(sum(is.na(out_dat$N)) > 0.5 * dim(out_dat)[1]) {
        col.name = colnames(out_dat)
        col.name = col.name[!col.name %in% "N"]
        out_dat = out_dat[,col.name]
    }

    write.table(out_dat,file = outputfile,quote=FALSE,row.names=FALSE,col.names=TRUE)
}
