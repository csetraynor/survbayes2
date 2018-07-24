
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


for(i in 1:11){
  my_script <- paste0("#!/bin/bash
#PBS -l nodes=14:ppn=28
#PBS -l walltime=36:00:00
#PBS -l pmem=4571mb
    
module load intel/2017.2.174-GCC-6.3.0-2.27
module load impi/2017.2.174
module load R/3.4.3-X11-20170314
module load parallel GCC
    
MY_PARALLEL_OPTS=\"-N 1 --delay .2 -j $SLURM_NTASKS --joblog parallel-${SLURM_JOBID}.log\"
MY_SRUN_OPTS=\"-N 1 -n 1 --exclusive\"
MY_EXEC=\"R CMD BATCH '--args {1} ", i,"' parallel.R\"
    
parallel $MY_PARALLEL_OPTS srun $MY_SRUN_OPTS $MY_EXEC ::: {1..100}")
  write(my_script, paste0("pbs_intclust_", i,".txt") )
}
