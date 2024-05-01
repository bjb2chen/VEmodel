#!/bin/bash
#SBATCH --array=1-25
#SBATCH --nodes=1
#SBATCH --account=def-mnooijen 
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=40:00:00 
echo "Starting task $SLURM_ARRAY_TASK_ID"
INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" beluga_apr28_bilinears.txt)    #get a single line from input_list.txt according to SLURM_ARRAY_TASK_ID and assign it to "INPUT"
/home/bjb2chen/LOCAL/runG_mrsftd $INPUT 16                                                                     #process the input file

# cd $SLURM_SUBMIT_DIR

# export SLURM_CPUS_PER_TASK
# mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID


