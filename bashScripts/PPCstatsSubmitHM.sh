#!/bin/bash
#SBATCH --account=p31578  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=short  ### PARTITION (buyin, short, normal, w10001, etc)
#SBATCH --array=1-99 ## number of jobs to run "in parallel" 
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on each computer
#SBATCH --time=4:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=3G ## how much RAM do you need per CPU (this affects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name="PH\${SLURM_ARRAY_TASK_ID}" ## use the task id in the name of the job
#SBATCH --output=PH.%a.out ## use the jobid (A) and the specific job index (a) to name your log file
#SBATCH --error=PH.%a.err 


module purge all
module load matlab/r2022b

echo "getting started: ${SLURM_ARRAY_TASK_ID}"

matlab -singleCompThread -batch "start=${SLURM_ARRAY_TASK_ID}; PCCstatsWrapperHM;"
