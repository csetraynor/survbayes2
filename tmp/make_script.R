
for(i in 1:100){
  my_script <- paste0("#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntask-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4571mb
#SBATCH --time=1-0:00:00     
#SBATCH --output=my.stdout
#SBATCH --job-name=", "job_test", i, "1",
"#Call R
                      
module load intel impi imkl R
R CMD BATCH '--args ", i, "' parallel_test.R")
  write(my_script, paste0("submitscript", i,".txt") )
  }
