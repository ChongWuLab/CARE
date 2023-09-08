fourth_args <- c("Prop1","Prop2","Prop3")
last_args <- 11:20

path <- "/rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/pbs/"


for (fourth_arg in fourth_args){
  for (last_arg in last_args){
    file_name <- paste0(path, "sim1_4_",fourth_arg,"_",last_arg,".lsf")
    # Create file content
    file_content <- paste0(
      "#BSUB -W 4:00\n",
      "#BSUB -o /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/log1_4_",fourth_arg,"_",last_arg, "\n",
      "#BSUB -e /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/error1_4_",fourth_arg,"_",last_arg, "\n",
      "#BSUB -cwd /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/\n",
      "#BSUB -q medium\n",
      "#BSUB -n 12\n",
      "#BSUB -M 17\n",
      "#BSUB -R rusage[mem=17]\n",
      "#BSUB -B\n",
      "#BSUB -N\n",
      "#BSUB -J simulate\n",
      "\n",
      "module load R/4.2.1\n",
      "Rscript simulation_dir_pleio.R 4 thetaU1 N6 ", fourth_arg, " ", last_arg
    )
    writeLines(file_content, file_name)
    
  }
}

for (fourth_arg in fourth_args){
  for (last_arg in last_args){
    file_name <- paste0(path, "sim2_4_",fourth_arg,"_",last_arg,".lsf")
    # Create file content
    file_content <- paste0(
      "#BSUB -W 4:00\n",
      "#BSUB -o /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/log2_4_",fourth_arg,"_",last_arg, "\n",
      "#BSUB -e /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/error2_4_",fourth_arg,"_",last_arg, "\n",
      "#BSUB -cwd /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/\n",
      "#BSUB -q medium\n",
      "#BSUB -n 12\n",
      "#BSUB -M 17\n",
      "#BSUB -R rusage[mem=17]\n",
      "#BSUB -B\n",
      "#BSUB -N\n",
      "#BSUB -J simulate\n",
      "\n",
      "module load R/4.2.1\n",
      "Rscript simulation_dir_pleio2.R 4 thetaU1 N6 ", fourth_arg, " ", last_arg
    )
    writeLines(file_content, file_name)
    
  }
}

for (fourth_arg in fourth_args){
  for (last_arg in last_args){
    file_name <- paste0(path, "sim3_4_",fourth_arg,"_",last_arg,".lsf")
    # Create file content
    file_content <- paste0(
      "#BSUB -W 4:00\n",
      "#BSUB -o /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/log3_4_",fourth_arg,"_",last_arg, "\n",
      "#BSUB -e /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/error3_4_",fourth_arg,"_",last_arg, "\n",
      "#BSUB -cwd /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/\n",
      "#BSUB -q medium\n",
      "#BSUB -n 12\n",
      "#BSUB -M 17\n",
      "#BSUB -R rusage[mem=17]\n",
      "#BSUB -B\n",
      "#BSUB -N\n",
      "#BSUB -J simulate\n",
      "\n",
      "module load R/4.2.1\n",
      "Rscript simulation_balance_pleio.R 4 thetaU1 N6 ", fourth_arg, " ", last_arg
    )
    writeLines(file_content, file_name)
    
  }
}

for (fourth_arg in fourth_args){
  for (last_arg in last_args){
    file_name <- paste0(path, "sim4_4_",fourth_arg,"_",last_arg,".lsf")
    # Create file content
    file_content <- paste0(
      "#BSUB -W 4:00\n",
      "#BSUB -o /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/log4_4_",fourth_arg,"_",last_arg, "\n",
      "#BSUB -e /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/error4_4_",fourth_arg,"_",last_arg, "\n",
      "#BSUB -cwd /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/\n",
      "#BSUB -q medium\n",
      "#BSUB -n 12\n",
      "#BSUB -M 17\n",
      "#BSUB -R rusage[mem=17]\n",
      "#BSUB -B\n",
      "#BSUB -N\n",
      "#BSUB -J simulate\n",
      "\n",
      "module load R/4.2.1\n",
      "Rscript simulation_balance_pleio2.R 4 thetaU1 N6 ", fourth_arg, " ", last_arg
    )
    writeLines(file_content, file_name)
    
  }
}
