#!/bin/bash

if [ $# -eq 0 ]
then
  echo ""
  echo "dist gamess_hess.out"
  echo ""
  exit 1
fi

hessout=$1
filnam="nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"

natoms=`grep ' TOTAL NUMBER OF ATOMS' $hessout|cut -d'=' -f2`
ndim=`echo "$natoms * 3"|bc -l`
echo "Dimension of all xyz coordinates:" $ndim
natom=`echo "$ndim / 3"|bc -l|cut -d'.' -f1`
echo "# of atoms:" $natom

ngroup=`echo "$ndim / 5"|bc -l|cut -d'.' -f1`
nleft=`echo "$ndim - 5* $ngroup"|bc -l`
echo $ngroup $nleft

declare -A nrmmod
declare -A freqcm

cat $hessout |sed -n "/FREQUENCIES IN CM/,/REFERENCE ON SAYVETZ CONDITIONS/p"|grep -A3 '^..[0-9]'|cut -c21- > mode.dat
cat $hessout |sed -n "/FREQUENCIES IN CM/,/REFERENCE ON SAYVETZ CONDITIONS/p"|grep 'FREQUENCY:'|cut -c19- > freq.dat

for igroup in $(seq 1 1 $ngroup)
do
  iniline=`echo " ( $igroup -1 ) * ( $ndim + 2 ) + 1"|bc -l`
  endline=`echo " ($iniline - 1 + $ndim)"|bc -l`
  echo "igroup=" $igroup
  echo "iniline=" $iniline
  ixyz=0
  for line in $(seq $iniline 1 $endline)
  do
#   echo $line
    ixyz=`echo "$ixyz + 1"|bc -l`
    echo $ixyz
    for icolumn in $(seq 1 1 5)
    do
      imode=`echo "($igroup - 1)*5 + $icolumn"|bc -l`
      echo -n $ixyz $imode ""
      cutini=`echo "($icolumn - 1)*12 + 1"|bc -l`
      cutfnl=`echo "($icolumn )*12 "|bc -l`
#     echo $cutini $cutfnl
      disp=`sed -n "$line p" < mode.dat |cut -c$cutini"-"$cutfnl`
#     sed -n "$line p" < mode.dat 
      nrmmod[$ixyz,$imode]=$disp
      echo ${nrmmod[$ixyz,$imode]} $disp

# Now read in frequencies in cm-1
      if [ $ixyz = 1 ]
      then
        freq=`sed -n "$igroup p" < freq.dat |cut -c$cutini"-"$cutfnl`
        freqcm[$imode]=$freq
        echo frequency: $imode ${freqcm[$imode]}
      fi
    done
  done

done





#for the leftover nleft modes
if [ $nleft != 0 ]
then
for igroup in $(seq $nleft 1 $nleft)
do
  iniline=`echo " $ngroup * ( $ndim + 2 ) + 1"|bc -l`
  endline=`echo " ($iniline - 1 + $ndim)"|bc -l`
  echo "igroup=leftover"
  echo "iniline=" $iniline
  ixyz=0
  for line in $(seq $iniline 1 $endline)
  do
#   echo $line
    ixyz=`echo "$ixyz + 1"|bc -l`
    echo $ixyz
    for icolumn in $(seq 1 1 $nleft)
    do
      imode=`echo "($ngroup)*5 + $icolumn"|bc -l`
      echo -n $ixyz $imode ""
      cutini=`echo "($icolumn - 1)*12 + 1"|bc -l`
      cutfnl=`echo "($icolumn )*12 "|bc -l`
#     echo $cutini $cutfnl
      disp=`sed -n "$line p" < mode.dat |cut -c$cutini"-"$cutfnl`
#     sed -n "$line p" < mode.dat
      nrmmod[$ixyz,$imode]=$disp
      echo ${nrmmod[$ixyz,$imode]} $disp

# Now read in frequencies in cm-1
      if [ $ixyz = 1 ]
      then
        freq=`tail -1 freq.dat |cut -c$cutini"-"$cutfnl`
        freqcm[$imode]=$freq
        echo frequency: $imode ${freqcm[$imode]}
      fi

    done
  done
done
fi

#Pring all frequencies
for imode in $(seq 1 1 $ndim)
do
  echo frequency: $imode ${freqcm[$imode]} CM-1
done

#read in reference structure
#The ref_structure has to be prepared by human-being and adopts the following format
# N           7.0   0.0000000000  -0.0000000000  -0.1693806842
# H           1.0  -0.4653267700   0.8059696078   0.2564602281
# H           1.0  -0.4653267700  -0.8059696078   0.2564602281
# H           1.0   0.9306535400   0.0000000000   0.2564602281
#The # symbol before each line shall not be in ref_structure.

declare -a atmlst
declare -a chrglst
declare -a refcoord

for iatom in $(seq 1 1 $natom)
do
  linecontent=`sed -n "$iatom p" < ref_structure`
  atmnam=`echo $linecontent|cut -d' ' -f1`
  atmlst[$iatom]=$atmnam
  chrglst[$iatom]=`echo $linecontent|cut -d' ' -f2`
  echo ${atmlst[$iatom]} ${chrglst[$iatom]}
  for ixyz in $(seq 1 1 3)
  do
    icomp=`echo "($iatom - 1)*3 + $ixyz"|bc -l`
#   echo $icomp
    ifield=`echo "$ixyz + 2"|bc -l`
    refcoord[$icomp]=`echo $linecontent|cut -d' ' -f$ifield`
    echo ${refcoord[$icomp]}
  done
done

#echo ${atmlst[@]}

declare -a distcoord_plus
declare -a distcoord_minus
declare -a distcoord_plus_x2
declare -a distcoord_minus_x2
declare -a distcoord_pp
declare -a distcoord_pm
declare -a distcoord_mp
declare -a distcoord_mm
#first copy refcoord to distcoord_plus and dist_coord_minus
#for icomp in $(seq 1 1 $ndim)
#do
#  distcoord_plus[$icomp]=${refcoord[$icomp]}
#  distcoord_minus[$icomp]=${refcoord[$icomp]}
# echo ${distcoord[$icomp]}
#done

#List of modes not considered
modes_excluded=(1 2 3 4 5 6)
nexclud=${#modes_excluded[@]}
echo "NUMBER OF EXCLUDED MODES:" "${#modes_excluded[@]}" $nexclud
echo "They are modes:" ${modes_excluded[@]}
nmodes_include=`echo " $ndim - $nexclud  "|bc -l`

declare -A modes_included
icount=0
for imode in $(seq 1 1 $ndim)
do
# check whether the mode is considered
  include=1
# The modes_excluded array starts counting by 0
  nexclud_end=`echo "$nexclud - 1"|bc -l`
  for iexclud in $(seq 0 1 $nexclud_end)
  do
    if [ $imode = ${modes_excluded[$iexclud]} ]
    then
      include=0
    fi
  done

# if mode considered
  if [ $include = 1 ]
  then

    icount=`echo " $icount + 1"|bc -l`
    modes_included[$icount]=$imode
    echo $icount ${modes_included[$icount]}
  fi
done
nmodes_included=$icount
echo "Number of Modes Included:" $nmodes_included


#Do diabatization calculation at the reference nondistorted structure.
#This calculation shall be a repetition of a calcualtion in preparing temp.inp
grep grace "$filnam"_refG.out
if [ $? -ne 0 ]
then
  echo "Run calculation at the undistorted reference structure"
  cp temp.inp "$filnam"_refG.inp
  cat ref_structure >> "$filnam"_refG.inp
  echo ' $END    ' >> "$filnam"_refG.inp
  # You may need to edit the following submission line based on your computer's configuration
  subgam.diab "$filnam"_refG.inp 4 0 1
else
  echo "Calculation at the reference structure is done."
fi

#Loop over all considered modes and do + and - displacements
#set step size in reduced dimensionless coordinates
qsize=0.05
#Set conversion constants
ha2ev=27.2113961318
wn2ev=0.000123981
wn2eh=0.00000455633
ang2br=1.889725989
amu2me=1822.88839

for kmode in $(seq 1 1 $nmodes_included)
do
  imode=${modes_included[$kmode]}


# Convert the reduced dimensionless qsize to the actual rsize in sqrt(amu)*Angs unit
  omga=${freqcm[$imode]}
  rsize=`echo "$qsize / ( sqrt($amu2me)*$ang2br*sqrt($omga*$wn2eh) ) "|bc -l`
  echo $imode $omga $rsize

# Loop over components
  for icomp in $(seq 1 1 $ndim)
  do
    coord_disp_plus=`echo " ${refcoord[$icomp]} + $rsize * ${nrmmod[$icomp,$imode]}  "|bc -l`
    coord_disp_minus=`echo " ${refcoord[$icomp]} - $rsize * ${nrmmod[$icomp,$imode]}  "|bc -l`
    distcoord_plus[$icomp]=$coord_disp_plus
    distcoord_minus[$icomp]=$coord_disp_minus
    coord_disp_plusx2=`echo " ${refcoord[$icomp]} + 2.0 * $rsize * ${nrmmod[$icomp,$imode]}  "|bc -l`
    coord_disp_minusx2=`echo " ${refcoord[$icomp]} - 2.0 * $rsize * ${nrmmod[$icomp,$imode]}  "|bc -l`
    distcoord_plus_x2[$icomp]=$coord_disp_plusx2
    distcoord_minus_x2[$icomp]=$coord_disp_minusx2
    echo $imode $icomp ${refcoord[$icomp]} ${nrmmod[$icomp,$imode]} $coord_disp_plus $coord_disp_minus ${distcoord_plus[$icomp]} ${distcoord_minus[$icomp]}
  done

  # delete existing dist_structure_plus and dist_structure_minus
  if [ -f dist_structure_plus ]
  then
    rm -f dist_structure_plus
  fi
  if [ -f dist_structure_minus ]
  then
    rm -f dist_structure_minus
  fi
  # delete existing dist_structure_plusx2 and dist_structure_minusx2
  if [ -f dist_structure_plusx2 ]
  then
    rm -f dist_structure_plusx2
  fi
  if [ -f dist_structure_minusx2 ]
  then
    rm -f dist_structure_minusx2
  fi

  # print the distorted structure
  for iatom in $(seq 1 1 $natom)
  do
    echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure_plus
    echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure_minus
    echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure_plusx2
    echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure_minusx2
    for ixyz in $(seq 1 1 3)
    do
      icomp=`echo "($iatom - 1)*3 + $ixyz"|bc -l`
      echo -n ${distcoord_plus[$icomp]} "" >> dist_structure_plus
      echo -n ${distcoord_minus[$icomp]} "" >> dist_structure_minus
      echo -n ${distcoord_plus_x2[$icomp]} "" >> dist_structure_plusx2
      echo -n ${distcoord_minus_x2[$icomp]} "" >> dist_structure_minusx2
    done
    echo "" >> dist_structure_plus
    echo "" >> dist_structure_minus
    echo "" >> dist_structure_plusx2
    echo "" >> dist_structure_minusx2
  done
  
  # The temp.inp template input file has to be prepared by human-being
  # I am just directing the distorted structures into the temp.inp to generate the actual input for gamess
  cp temp.inp "$filnam"_mode"$imode"_+"$qsize".inp
  cat dist_structure_plus >> "$filnam"_mode"$imode"_+"$qsize".inp
  echo ' $END ' >> "$filnam"_mode"$imode"_+"$qsize".inp
  cp temp.inp "$filnam"_mode"$imode"_-"$qsize".inp
  cat dist_structure_minus >> "$filnam"_mode"$imode"_-"$qsize".inp
  echo ' $END ' >> "$filnam"_mode"$imode"_-"$qsize".inp

  cp temp.inp "$filnam"_mode"$imode"_+"$qsize"x2.inp
  cat dist_structure_plusx2 >> "$filnam"_mode"$imode"_+"$qsize"x2.inp
  echo ' $END ' >> "$filnam"_mode"$imode"_+"$qsize"x2.inp
  cp temp.inp "$filnam"_mode"$imode"_-"$qsize"x2.inp
  cat dist_structure_minusx2 >> "$filnam"_mode"$imode"_-"$qsize"x2.inp
  echo ' $END ' >> "$filnam"_mode"$imode"_-"$qsize"x2.inp

  # submit gamess diabatization calculations for the + and - inputs
  #Check whether the calculation is done already.
  grep grace "$filnam"_mode"$imode"_+"$qsize".out
  if [ $? -ne 0 ]
  then
    echo "run calculations for" "$filnam"_mode"$imode"_+"$qsize"
  # You may need to edit the following submission line based on your computer's configuration
    subgam.diab "$filnam"_mode"$imode"_+"$qsize".inp 4 0 1
  else
    echo "$filnam"_mode"$imode"_+"$qsize" "is done"
  fi
  #Check whether the calculation is done already.
  grep grace "$filnam"_mode"$imode"_-"$qsize".out
  if [ $? -ne 0 ]
  then
    echo "run calculations for" "$filnam"_mode"$imode"_-"$qsize"
  # You may need to edit the following submission line based on your computer's configuration
    subgam.diab "$filnam"_mode"$imode"_-"$qsize".inp 4 0 1
  else
    echo "$filnam"_mode"$imode"_-"$qsize" "is done"
  fi
  grep grace "$filnam"_mode"$imode"_+"$qsize"x2.out
  if [ $? -ne 0 ]
  then
    echo "run calculations for" "$filnam"_mode"$imode"_+"$qsize"x2
  # You may need to edit the following submission line based on your computer's configuration
    subgam.diab "$filnam"_mode"$imode"_+"$qsize"x2.inp 4 0 1
  else
    echo "$filnam"_mode"$imode"_+"$qsize"x2 "is done"
  fi
  grep grace "$filnam"_mode"$imode"_-"$qsize"x2.out
  if [ $? -ne 0 ]
  then
    echo "run calculations for" "$filnam"_mode"$imode"_-"$qsize"x2
  # You may need to edit the following submission line based on your computer's configuration
    subgam.diab "$filnam"_mode"$imode"_-"$qsize"x2.inp 4 0 1
  else
    echo "$filnam"_mode"$imode"_-"$qsize"x2 "is done"
  fi

  #2D distortion to get bilinear vibronic coupling
  lmode_last=`echo "$kmode - 1"|bc -l`
  for lmode in $(seq 1 1 $lmode_last)
  do
    jmode=${modes_included[$lmode]}
    #Convert the reduced dimensionless qsize to the actual rsizep in sqrt(amu)*Angs unit for jmode
    omgap=${freqcm[$jmode]}
    rsizep=`echo "$qsize / ( sqrt($amu2me)*$ang2br*sqrt($omgap*$wn2eh) ) "|bc -l`
    echo $imode $jmode $rsize $rsizep

    # Loop over components
    for icomp in $(seq 1 1 $ndim)
    do
      coord_disp_pp=`echo " ${distcoord_plus[$icomp]} + $rsizep * ${nrmmod[$icomp,$jmode]}  "|bc -l`
      coord_disp_pm=`echo " ${distcoord_plus[$icomp]} - $rsizep * ${nrmmod[$icomp,$jmode]}  "|bc -l`
      coord_disp_mp=`echo " ${distcoord_minus[$icomp]} + $rsizep * ${nrmmod[$icomp,$jmode]}  "|bc -l`
      coord_disp_mm=`echo " ${distcoord_minus[$icomp]} - $rsizep * ${nrmmod[$icomp,$jmode]}  "|bc -l`
      distcoord_pp[$icomp]=$coord_disp_pp
      distcoord_pm[$icomp]=$coord_disp_pm
      distcoord_mp[$icomp]=$coord_disp_mp
      distcoord_mm[$icomp]=$coord_disp_mm
    done

    # delete existing dist_structure_pp, _pm, _mp, _mm
    if [ -f dist_structure_pp ]
    then
      rm -f dist_structure_pp
    fi
    if [ -f dist_structure_pm ]
    then
      rm -f dist_structure_pm
    fi
    if [ -f dist_structure_mp ]
    then
      rm -f dist_structure_mp
    fi
    if [ -f dist_structure_mm ]
    then
      rm -f dist_structure_mm
    fi

  # print the distorted structure
    for iatom in $(seq 1 1 $natom)
    do
      echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure_pp
      echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure_pm
      echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure_mp
      echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure_mm
      for ixyz in $(seq 1 1 3)
      do
        icomp=`echo "($iatom - 1)*3 + $ixyz"|bc -l`
        echo -n ${distcoord_pp[$icomp]} "" >> dist_structure_pp
        echo -n ${distcoord_pm[$icomp]} "" >> dist_structure_pm
        echo -n ${distcoord_mp[$icomp]} "" >> dist_structure_mp
        echo -n ${distcoord_mm[$icomp]} "" >> dist_structure_mm
      done
      echo "" >> dist_structure_pp
      echo "" >> dist_structure_pm
      echo "" >> dist_structure_mp
      echo "" >> dist_structure_mm
  done

  cp temp.inp "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize".inp
  cat dist_structure_pp >> "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize".inp
  echo ' $END ' >> "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize".inp

  cp temp.inp "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize".inp
  cat dist_structure_pm >> "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize".inp
  echo ' $END ' >> "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize".inp

  cp temp.inp "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize".inp
  cat dist_structure_mp >> "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize".inp
  echo ' $END ' >> "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize".inp

  cp temp.inp "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize".inp
  cat dist_structure_mm >> "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize".inp
  echo ' $END ' >> "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize".inp

  # submit gamess diabatization calculations for the bilinear inputs
  #Check whether the calculation is done already.
  grep grace "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize".out
  if [ $? -ne 0 ]
  then
    echo "run calculations for" "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize"
  # You may need to edit the following submission line based on your computer's configuration
    subgam.diab "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize".inp 4 0 1
  else
    echo "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize" "is done"
  fi
  #Check whether the calculation is done already.
  grep grace "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize".out
  if [ $? -ne 0 ]
  then
    echo "run calculations for" "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize"
  # You may need to edit the following submission line based on your computer's configuration
    subgam.diab "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize".inp 4 0 1
  else
    echo "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize" "is done"
  fi
  #Check whether the calculation is done already.
  grep grace "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize".out
  if [ $? -ne 0 ]
  then
    echo "run calculations for" "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize"
  # You may need to edit the following submission line based on your computer's configuration
    subgam.diab "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize".inp 4 0 1
  else
    echo "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize" "is done"
  fi
  #Check whether the calculation is done already.
  grep grace "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize".out
  if [ $? -ne 0 ]
  then
    echo "run calculations for" "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize"
  # You may need to edit the following submission line based on your computer's configuration
    subgam.diab "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize".inp 4 0 1
  else
    echo "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize" "is done"
  fi

  done #done loop over jmode


done


#Now we move on to extract vibronic coupling constants using finite difference
#and write the data in an mctdh operator file

#Remove the existing mctdh.op file
if [ -f mctdh.op ]
then
  rm -f mctdh.op
fi

#Heading for mctdh.op
echo "OP_DEFINE-SECTION" >> mctdh.op
echo "title" >> mctdh.op
nstate=`grep '# of states in CI      = ' "$filnam"_refG.out|tail -1|cut -d'=' -f2`
echo "$filnam $nstate states + $nmodes_include modes" >> mctdh.op
echo "end-title " >> mctdh.op
echo "end-op_define-section" >> mctdh.op
echo "" >> mctdh.op
echo "PARAMETER-SECTION" >> mctdh.op
echo "" >> mctdh.op

#Extract diabatic Hamiltonian in the reference structure
if [ -f "$filnam"_refG.out ]
then
  echo "#Diagonal and Off-diagonal diabatic Hamitonian elements at reference structure" >> mctdh.op
  for ist in $(seq 1 1 $nstate)
  do
    Ediab=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_refG.out|tail -1|cut -c62-|sed -e "s/ //g"`
    echo v"$ist" = "$Ediab"", ev" >> mctdh.op

    jlast=`echo " $ist -1"|bc -l`
    for jst in $(seq 1 1 $jlast)
    do
      Coup_ev=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_refG.out|tail -1|cut -c62-|sed -e "s/ //g"`
      echo v"$jst$ist" = "$Coup_ev"", ev" >> mctdh.op
    done
    echo "" >> mctdh.op
  done
  echo "" >> mctdh.op
  else
  echo "Skip extracting Hamiltonians from the non-existing " "$filnam"_refG.out
fi



for kmode in $(seq 1 1 $nmodes_included)
do
  imode=${modes_included[$kmode]}



  vibron_ev=`echo " ${freqcm[$imode]} * $wn2ev "|bc -l`
  echo "#Parameters for mode" $imode >> mctdh.op
  echo "#Vibron:" >> mctdh.op
  echo w_m"$imode"= "$vibron_ev"", ev" >> mctdh.op
  echo "" >> mctdh.op
  echo "#Linear and quadratic diagonal and off-diagonal vibronic coupling constants:" >> mctdh.op
  grep grace "$filnam"_mode"$imode"_+"$qsize".out
  grace_code_plus=$?
  grep grace "$filnam"_mode"$imode"_-"$qsize".out
  grace_code_minus=$?
  grep grace "$filnam"_mode"$imode"_+"$qsize"x2.out
  grace_code_plusx2=$?
  grep grace "$filnam"_mode"$imode"_-"$qsize"x2.out
  grace_code_minusx2=$?
  echo $imode $grace_code_plus $grace_code_minus

  if [ $grace_code_plus -eq 0 ] && [ $grace_code_minus -eq 0 ] && [ $grace_code_plusx2 -eq 0 ] && [ $grace_code_minusx2 -eq 0 ]
  then
    echo "good to extract"
#   nstate=`grep '# of states in CI      = ' "$filnam"_mode"$imode"_+"$qsize".out|tail -1|cut -d'=' -f2`
#   echo $nstate

    #Extract the diagonal and off-diagonal vibronic coupling
    for ist in $(seq 1 1 $nstate)
    do
      Ediab_au_plus=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_mode"$imode"_+"$qsize".out|tail -1|cut -c44-61`
      Ediab_au_plusx2=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_mode"$imode"_+"$qsize"x2.out|tail -1|cut -c44-61`
      Ediab_au_minus=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_mode"$imode"_-"$qsize".out|tail -1|cut -c44-61`
      Ediab_au_minusx2=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_mode"$imode"_-"$qsize"x2.out|tail -1|cut -c44-61`
      Ediab_au_0=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_refG.out|tail -1|cut -c44-61`
      linear_diag_ev=`echo " ( $Ediab_au_plus - $Ediab_au_minus ) * $ha2ev / ( 2 * $qsize ) "|bc -l`
      quadratic_diag_ev=`echo " ($Ediab_au_plusx2 + $Ediab_au_minusx2 - 2.0 * $Ediab_au_0) * $ha2ev / ( 4.0 * $qsize * $qsize ) "|bc -l`
      #We only view the difference between the actual force constant and the vibron as the quadratic diagonal coupling for the diabatic state.
      #This is because we add the vibron to the diagonal Hamiltonian matrix elements for all diabats.
      #So, the quadratic diagonal vibronic coupling is actually the correction to the vibron of the normal mode as a force constant.
      quadratic_diag_ev=`echo " $quadratic_diag_ev - $vibron_ev  "|bc -l`
      echo $ist $linear_diag_ev $quadratic_diag_ev
      echo l"$ist"_m"$imode" = "$linear_diag_ev"", ev" >> mctdh.op
      echo q"$ist"_m"$imode" = "$quadratic_diag_ev"", ev" >> mctdh.op
      #loop over jst
      jlast=`echo " $ist -1"|bc -l`
      for jst in $(seq 1 1 $jlast)
      do
        Coup_ev_plus=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_mode"$imode"_+"$qsize".out|tail -1|cut -c62-`
        Coup_ev_minus=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_mode"$imode"_-"$qsize".out|tail -1|cut -c62-`
        Coup_ev_plusx2=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_mode"$imode"_+"$qsize"x2.out|tail -1|cut -c62-`
        Coup_ev_minusx2=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_mode"$imode"_-"$qsize"x2.out|tail -1|cut -c62-`
        Coup_ev_refG=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_refG.out|tail -1|cut -c62-`
        linear_offdiag_ev=`echo " ( $Coup_ev_plus - $Coup_ev_minus ) / ( 2 * $qsize) "|bc -l`
	quadratic_offdiag_ev=`echo "( $Coup_ev_plusx2 + $Coup_ev_minusx2 - 2.0 * $Coup_ev_refG ) / ( 4.0 * $qsize * $qsize )"|bc -l`
        echo $jst $ist $linear_offdiag_ev
        echo l"$jst$ist"_m"$imode" = "$linear_offdiag_ev"", ev" >> mctdh.op
        echo q"$jst$ist"_m"$imode" = "$quadratic_offdiag_ev"", ev" >> mctdh.op
      done
      echo "" >> mctdh.op
    done

  else
    echo "not good to extract. Skipping mode" $imode "for extracting vibronic couplings"
  fi

  #Extracting bilinear vibronic coupling
  echo "#Bilinear diagonal and off-diagonal vibronic coupling constants:" >> mctdh.op
  lmode_last=`echo "$kmode - 1"|bc -l`
  for lmode in $(seq 1 1 $lmode_last)
  do
    jmode=${modes_included[$lmode]}

    #Check whether the four displacements are calculated properly
    grep grace "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize".out
    grace_code_pp=$?
    grep grace "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize".out
    grace_code_pm=$?
    grep grace "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize".out
    grace_code_mp=$?
    grep grace "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize".out
    grace_code_mm=$?

    if [ $grace_code_pp -eq 0 ] && [ $grace_code_pm -eq 0 ] && [ $grace_code_mp -eq 0 ] && [ $grace_code_mm -eq 0 ]
    then
      echo "Good to extract bilinear for modes $imode $jmode"
      for ist in $(seq 1 1 $nstate)
      do
        Ediab_au_pp=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize".out|tail -1|cut -c44-61`
        Ediab_au_pm=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize".out|tail -1|cut -c44-61`
        Ediab_au_mp=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize".out|tail -1|cut -c44-61`
        Ediab_au_mm=`grep "STATE #..* $ist.S GMC-PT-LEVEL DIABATIC ENERGY=" "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize".out|tail -1|cut -c44-61`
	bilinear_diag_ev=`echo "( $Ediab_au_pp + $Ediab_au_mm - $Ediab_au_pm - $Ediab_au_mp ) * $ha2ev / (4.0 * $qsize * $qsize )"|bc -l`

        echo $ist $bilinear_diag_ev 
        echo b"$ist"_m"$imode"_m"$jmode" = "$bilinear_diag_ev"", ev" >> mctdh.op

        #loop over jst
        jlast=`echo " $ist -1"|bc -l`
        for jst in $(seq 1 1 $jlast)
        do
          Coup_ev_pp=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_+"$qsize".out|tail -1|cut -c62-`
          Coup_ev_pm=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_mode"$imode"_+"$qsize"_mode"$jmode"_-"$qsize".out|tail -1|cut -c62-`
          Coup_ev_mp=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_+"$qsize".out|tail -1|cut -c62-`
          Coup_ev_mm=`grep "STATE #..* $jst &..* $ist.S GMC-PT-LEVEL COUPLING" "$filnam"_mode"$imode"_-"$qsize"_mode"$jmode"_-"$qsize".out|tail -1|cut -c62-`
          bilinear_offdiag_ev=`echo " ( $Coup_ev_pp + $Coup_ev_mm - $Coup_ev_pm - $Coup_ev_mp ) / (4.0 * $qsize * $qsize ) "|bc -l`

          echo $jst $ist $bilinear_offdiag_ev
          echo b"$jst$ist"_m"$imode"_m"$jmode" = "$bilinear_offdiag_ev"", ev" >> mctdh.op
        done
        echo "" >> mctdh.op
      done
    else
      echo "not good to extract. Skipping mode $imode mode $jmode for extracting bilinear vibronic couplings"
    fi

  done


done


echo "end-parameter-section" >> mctdh.op

#Now we construct the Hamiltonian section in mctdh.op
echo "-----------------------------------------" >> mctdh.op
echo "HAMILTONIAN-SECTION" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op


echo -n " modes | el " >> mctdh.op
for imode_include in $(seq 1 1 $nmodes_include)
do
  echo -n "| m${modes_included[$imode_include]} " >> mctdh.op
done 
echo "" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op
echo "# KINETIC OPERATOR FOR NORMAL MODES" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op
for imode_include in $(seq 1 1 $nmodes_include)
do
  mode_count=`echo " $imode_include + 1 "|bc -l`
 #echo $imode_include $mode_count
  echo "w_${modes_included[$imode_include]}   |$mode_count KE" >> mctdh.op
done
echo "-----------------------------------------" >> mctdh.op
echo "# HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op
for imode_include in $(seq 1 1 $nmodes_include)
do
  mode_count=`echo " $imode_include + 1 "|bc -l`
 #echo $imode_include $mode_count
  echo "0.5*w_${modes_included[$imode_include]}   |$mode_count  q^2 " >> mctdh.op
done
echo "-----------------------------------------" >> mctdh.op
echo "# ELECTRONIC COUPLING AT REFERENCE STRUCTURE" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op
for ist in $(seq 1 1 $nstate)
do
  echo "v$ist  |1 S""$ist""&""$ist" >> mctdh.op
done
for ist in $(seq 1 1 $nstate)
do
  jlast=`echo "$ist -1"|bc -l`
  for jst in $(seq 1 1 $jlast)
  do
    echo "v$jst$ist  |1 S""$jst""&""$ist" >> mctdh.op
  done
done
echo "-----------------------------------------" >> mctdh.op
echo "# LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op
for kmode in $(seq 1 1 $nmodes_include)
do
  imode=${modes_included[$kmode]}
  kmode_count=`echo "$kmode + 1"|bc -l`
  for ist in $(seq 1 1 $nstate)
  do
    echo l"$ist"_m"$imode |1" "S$ist""&""$ist |$kmode_count q" >> mctdh.op 
    echo q"$ist"_m"$imode |1" "S$ist""&""$ist |$kmode_count q^2" >> mctdh.op 
  done
done
echo "-----------------------------------------" >> mctdh.op
echo "# LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op
for kmode in $(seq 1 1 $nmodes_include)
do
  imode=${modes_included[$kmode]}
  kmode_count=`echo "$kmode + 1"|bc -l`
  for ist in $(seq 1 1 $nstate)
  do
    jlast=`echo "$ist -1"|bc -l`
    for jst in $(seq 1 1 $jlast)
    do
      echo l"$jst$ist"_m"$imode |1" "S$jst""&""$ist |$kmode_count q" >> mctdh.op
      echo q"$jst$ist"_m"$imode |1" "S$jst""&""$ist |$kmode_count q^2" >> mctdh.op
    done
  done
done
echo "-----------------------------------------" >> mctdh.op
echo "# BILINEAR DIAGONAL VIBRONIC COUPLINGS" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op
for kmode in $(seq 1 1 $nmodes_include)
do
  imode=${modes_included[$kmode]}
  kmode_count=`echo "$kmode +1"|bc -l`
  lmode_last=`echo "$kmode -1"|bc -l`
  for lmode in $(seq 1 1 $lmode_last)
  do
    jmode=${modes_included[$lmode]}
    lmode_count=`echo "$lmode +1"|bc -l`
    for ist in $(seq 1 1 $nstate)
    do
      echo b"$ist"_m"$imode"_m"$jmode |1 S""$ist""&""$ist |$lmode_count q |$kmode_count q" >> mctdh.op 
    done
  done
done
echo "-----------------------------------------" >> mctdh.op
echo "# BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS" >> mctdh.op
echo "-----------------------------------------" >> mctdh.op
for kmode in $(seq 1 1 $nmodes_include)
do
  imode=${modes_included[$kmode]}
  kmode_count=`echo "$kmode +1"|bc -l`
  lmode_last=`echo "$kmode -1"|bc -l`
  for lmode in $(seq 1 1 $lmode_last)
  do
    jmode=${modes_included[$lmode]}
    lmode_count=`echo "$lmode +1"|bc -l`
    for ist in $(seq 1 1 $nstate)
    do
      jlast=`echo "$ist -1"|bc -l`
      for jst in $(seq 1 1 $jlast)
      do
        echo b"$jst$ist"_m"$imode"_m"$jmode |1 S""$jst""&""$ist |$lmode_count q |$kmode_count q" >> mctdh.op
      done
    done
  done
done
echo "-----------------------------------------" >> mctdh.op
echo "" >> mctdh.op
echo "end-hamiltonian-section" >> mctdh.op
echo "" >> mctdh.op
echo "end-operator" >> mctdh.op
