#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=6G
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

export SLURM_CPUS_PER_TASK
mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID

/home/bjb2chen/LOCAL/runG_diab FeCO_SPK_gmcpt_C1_mult1_TZ11o18e_1st.inp 12 
