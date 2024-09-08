source("step1_process_GWAS.R")
source("step2_LDSC.R")
source("step3_CARE.R")
source("CARE_support.R")
source("step4_otherMR.R")

main.dir = getwd()

#job.id  = 11
savedir = paste(main.dir,"/smmry_res/",sep="")

# finished jobs:
finsihed.jobs = list.files(savedir)



args <- commandArgs(TRUE)
outcome.name <- as.character(args[[1]])
exposure.name <- as.character(args[[2]])



# prepare the data
tryCatch({
    
    ieu_gwas_info = readRDS("IEU_GWAS_info.rds")

    smmry_dir = "/rsrch5/scratch/biostatistics/wzhang24/CARE/summarystats/converted/"

    outcome.sumstats = paste0(smmry_dir,outcome.name,".txt")
    outcome.n = subset(ieu_gwas_info,id==outcome.name)$sample_size
    
    exposure.sumstats = paste0(smmry_dir,exposure.name,".txt")
    exposure.n = subset(ieu_gwas_info,id==exposure.name)$sample_size
    
    save.file = paste0(savedir,exposure.name,"_",outcome.name,".RData")
    
    save.filejudge = paste0(exposure.name,"_",outcome.name,".RData")
    if(save.filejudge %in% finsihed.jobs) {
        cat("This job has already been finished previously. \n")
    } else {
        # create results directory:
        #exposure.name = "ebi-a-GCST005920"
        #outcome.name = "ukb-d-1747_1"
        workdir = paste0(main.dir,"/prepared_sumstats/",exposure.name,"_",outcome.name)
        
        dir.create(workdir)
        
        setwd(workdir)
        # prepare the GWAS data
        
        read_GWAS(exposure.sumstats,exposure.name,exposure.n)
        
        read_GWAS(outcome.sumstats,outcome.name,outcome.n)#
    }
    
}, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
})
    


