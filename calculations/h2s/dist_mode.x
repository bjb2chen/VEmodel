#!/bin/bash

if [ $# -eq 0 ]
then
  echo ""
  echo "dist gamess_hess.out"
  echo ""
  exit 1
fi

gmsout=$1

natoms=`grep ' TOTAL NUMBER OF ATOMS' $gmsout|cut -d'=' -f2`
ndim=`echo "$natoms * 3"|bc -l`
echo "Dimension of all xyz coordinates:" $ndim
natom=`echo "$ndim / 3"|bc -l|cut -d'.' -f1`
echo "# of atoms:" $natom

ngroup=`echo "$ndim / 5"|bc -l|cut -d'.' -f1`
nleft=`echo "$ndim - 5* $ngroup"|bc -l`
echo $ngroup $nleft

declare -A nrmmod

cat $gmsout|sed -n "/FREQUENCIES IN CM/,/REFERENCE ON SAYVETZ CONDITIONS/p"|grep -A3 '^..[0-9]'|cut -c21- > mode.dat

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
    done
  done
done
fi

#read in reference structure
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

declare -a distcoord
#first copy refcoord to distcoord
for icomp in $(seq 1 1 $ndim)
do
  distcoord[$icomp]=${refcoord[$icomp]}
# echo ${distcoord[$icomp]}
done

#The distorstion list below is case-specific
distlst=(0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 step)
#defined in this way, the index of distlst ranges from 0 to ndim-1
for imode in $(seq 1 1 $ndim)
do
  imodex=`echo "$imode - 1"|bc -l`
  echo $imode ${distlst[$imodex]}
# distort along each mode, one by one
  distrt=${distlst[$imodex]}
  for icomp in $(seq 1 1 $ndim)
  do
    coord_disp=`echo "${distcoord[$icomp]} + $distrt * ${nrmmod[$icomp,$imode]}"|bc -l`
#   echo $coord_disp
    distcoord[$icomp]=$coord_disp
  done
done

#delete existing dist_structure
if [ -f dist_structure ]
then
  rm -f dist_structure
fi
#print final distorted structure
for iatom in $(seq 1 1 $natom)
do
  echo -n ${atmlst[$iatom]} ${chrglst[$iatom]} "" >> dist_structure
  for ixyz in $(seq 1 1 3)
  do
    icomp=`echo "($iatom - 1)*3 + $ixyz"|bc -l`
    echo -n ${distcoord[$icomp]} "" >> dist_structure
  done 
  echo "" >> dist_structure
done


