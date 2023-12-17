#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

export SLURM_CPUS_PER_TASK
mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID

/home/bjb2chen/LOCAL/runG_diab RhF3_cat30_cct_rhf_D3h_e_NoValE.inp 3 
