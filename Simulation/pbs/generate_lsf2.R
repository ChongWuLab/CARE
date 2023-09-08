first_args <- 1:7
fourth_args <- c("Prop1","Prop2","Prop3")
last_args <- 1:10

path <- "/rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/pbs"

for (first_arg in first_args){
  for (fourth_arg in fourth_args){
    for (last_arg in last_args){
      file_name <- paste0(path, "/sim2_",first_arg,"_",fourth_arg,"_",last_arg,".lsf")
      # Create file content
      file_content <- paste0(
        "#BSUB -W 4:00\n",
        "#BSUB -o /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/log2_", first_arg, "_",fourth_arg,"_",last_arg, "\n",
        "#BSUB -e /rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/log/error2_", first_arg, "_",fourth_arg,"_",last_arg, "\n",
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
        "Rscript simulation_dir_pleio2.R ", first_arg, " thetaU1 N6 ", fourth_arg, " ", last_arg
      )
      writeLines(file_content, file_name)
      
    }
  }
}