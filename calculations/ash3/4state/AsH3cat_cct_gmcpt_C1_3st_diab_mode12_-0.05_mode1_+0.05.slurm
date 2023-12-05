#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

export SLURM_CPUS_PER_TASK
mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID

/home/$USER/LOCAL/runG_diab AsH3cat_cct_gmcpt_C1_3st_diab_mode12_-0.05_mode1_+0.05.inp 2 
