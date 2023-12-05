#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

export SLURM_CPUS_PER_TASK
mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID

/home/bjb2chen/LOCAL/runG_diab ash3cat_rohf_cct_gmcpt_7o7e_C1_4st.inp 4 
