#!/bin/bash
#SBATCH --array=1-30
#SBATCH --nodes=1
#SBATCH --account=def-mnooijen
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00

cd $SLURM_SUBMIT_DIR

export SLURM_CPUS_PER_TASK
mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID

# # Read the list of input files into an array
# readarray -t INPUT_FILES < beluga_apr28_bilinears.txt

# # Get the input file for this array task
# INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

# echo $INPUT_FILE

# Run your program with the input file
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q29_-x1q11.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q29_-x1q11.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_-x1q11.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_-x1q11.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_+x1q14.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_+x1q14.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_+x1q16.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_+x1q16.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_-x1q17.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_-x1q17.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q29_+x1q18.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q29_+x1q18.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_+x1q20.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_+x1q20.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q29_+x1q21.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q29_+x1q21.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_+x1q22.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_+x1q22.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_-x1q23.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q29_-x1q23.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q29_-x1q28.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q29_-x1q28.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_+x1q13.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_+x1q13.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_+x1q13.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_+x1q13.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q15.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q15.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_+x1q15.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_+x1q15.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q17.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q17.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_+x1q21.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_+x1q21.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_+x1q23.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_+x1q23.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q24.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q24.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q25.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q25.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_-x1q25.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_-x1q25.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_+x1q26.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_+x1q26.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q27.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q30_-x1q27.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_-x1q27.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_-x1q27.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_+x1q28.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q30_+x1q28.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q31_-x1q10.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q31_-x1q10.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q31_+x1q10.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q31_+x1q10.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q31_-x1q13.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_+x1q31_-x1q13.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q31_+x1q15.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q31_+x1q15.inp
echo $SLURM_ARRAY_TASK_ID
/home/bjb2chen/LOCAL/runG_mrsftd FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q31_-x1q7.inp 16
echo FeCO_SPK_gmcpt_C1_17st_diab_0.05_-x1q31_-x1q7.inp
echo $SLURM_ARRAY_TASK_ID
