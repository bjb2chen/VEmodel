#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

export SLURM_CPUS_PER_TASK
mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID

/home/bjb2chen/LOCAL/runG_diab H2O_ccd_rohf_C2v_mp2_gh.inp 12 
