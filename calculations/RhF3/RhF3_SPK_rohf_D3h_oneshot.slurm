#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=3G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

export SLURM_CPUS_PER_TASK
mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID

/home/bjb2chen/LOCAL/runG_diab RhF3_SPK_rhf_D3h_oneshot.inp 12 
