#!/bin/bash
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=4G
#SBATCH --time=222:00:00

cd $SLURM_SUBMIT_DIR

export SLURM_CPUS_PER_TASK
mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID

# Read the list of input files into an array
# readarray -t INPUT_FILES < input_list.txt

# Get the input file for this array task
# INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

#echo $INPUT_FILE
#echo $INPUT_FILES

# Run your program with the input file
/home/bjb2chen/LOCAL/runG_diab FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x2q29.inp 10
/home/bjb2chen/LOCAL/runG_diab FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q12_+x1q11.inp 10
