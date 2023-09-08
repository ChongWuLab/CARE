require(data.table)
setwd("/rsrch5/scratch/biostatistics/wzhang24/CARE/summarystats/covid19/release7/")
files = list.files()
files = files[grepl("tsv.gz",files)]
files = files[!grepl("tsv.gz.tbi",files)]

outcome.id = files[job.id]
out.id = gsub(".tsv.gz","",outcome.id)

#convert the data
raw = fread(outcome.id)
outcome_dat = as.data.frame(raw)
outputfile = paste0("/rsrch5/scratch/biostatistics/wzhang24/CARE/summarystats/covid19/release7/converted/",out.id,".txt")
out_dat = outcome_dat[,c("#CHR","POS","SNP","ALT","REF","all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_meta_AF","all_inv_var_meta_")]
