#!/bin/bash

filnam="h2scat_cct_gmcpt_6o7e_3st_diab"

for step in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8
do
  cat dist_mode.x | sed -e "s/step/$step/g" > dist_mode_1d.x
  chmod +x dist_mode_1d.x
  ./dist_mode_1d.x step1_h2s_cct_mp2_C2v_gh.out
  cp temp.inp  "$filnam"_as"$step".inp
# append the distorted structure to the input
  cat dist_structure >> "$filnam"_as"$step".inp 
# append " $END" to the input
  echo ' $END ' >> "$filnam"_as"$step".inp 
# append the data in diab_info to the input
  cat diab_info >> "$filnam"_as"$step".inp
# run gamess
  echo "running gamess"
#  runG_diab "$filnam"_as"$step".inp 4
# subgam.diab "$filnam"_as"$step".inp 48 0 1
  echo "done running gamess"
done
