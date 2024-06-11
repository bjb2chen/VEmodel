#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=0
#xBATCH --job-name=RhF3_Z1_H2_1000fs

python3 test_clean.py $1 $2 $3 $4 $5
#echo $1 $2 $3 $4 $5
