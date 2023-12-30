#!/bin/bash

# Syntax: combined_script.x GMS.out(or .dat)

file=$1

# Prompt the user to choose between "OPTIMIZED RHF" and "Semi-canonical MOs"
read -p $'Do you want:\n (1) OPTIMIZED RHF \n (2) OPTIMIZED ROHF \n (3) Semi-canonical MOs \n (4) MCSCF Natural Orbitals \n (5) OPTIMIZED MCSCF \n (6) DMO group \n (7) REFDET GROUP? \nYour choice: ' num

# Perform the corresponding operation based on the user's choice

# Handle OPTIMIZED RHF
if [ "$num" == "1" ]; then
  # Remove existing vec.dat
  if [ -f vec.dat ]; then
    rm -f vec.dat
  fi
  sed -n "/OPTIMIZED RHF/,/END/ p" $file > vec.dat
  echo "You have selected 1 - OPTIMIZED RHF"
  echo "OPTIMIZED RHF orbitals prepared in vec.dat file"

# Handle OPTIMIZED ROHF
elif [ "$num" == "2" ]; then
  # Remove existing vec.dat
  if [ -f vec.dat ]; then
    rm -f vec.dat
  fi
  sed -n "/OPTIMIZED ROHF/,/END/ p" $file > vec.dat
  echo "You have selected 2 - OPTIMIZED ROHF"
  echo "OPTIMIZED ROHF orbitals prepared in vec.dat file"

# Handle Semi-canonical MOs 
elif [ "$num" == "3" ]; then
  # Remove existing active_space_orbs
  if [ -f active_space_orbs ]; then
    rm -f active_space_orbs
  fi
  sed -n "/Semi-canonical MOs/,/END/ p" $file > active_space_orbs
  echo "You have selected 3 - Semi-canonical MOs"
  echo "Semi-canonical MOs prepared in active_space_orbs file"

# Handle MCSCF Natural Orbitals
elif [ "$num" == "4" ]; then
  # Remove existing nat_orbs_mcscf
  if [ -f nat_orbs_mcscf ]; then
    rm -f nat_orbs_mcscf
  fi
  sed -n "/NATURAL ORBITALS OF MCSCF/,/END/ p" $file > nat_orbs_mcscf
  echo "You have selected 4 - MCSCF Natural Orbitals"
  echo "Semi-canonical MOs prepared in nat_orbs_mcscf file"

# Handle OPTIMIZED MCSCF
elif [ "$num" == "5" ]; then
  # Remove existing nat_orbs_mcscf
  if [ -f optimized_mscsf ]; then
    rm -f optimized_mscsf
  fi
  sed -n "/OPTIMIZED MCSCF/,/END/ p" $file > optimized_mscsf
  echo "You have selected 5 - OPTIMIZED MCSCF"
  echo "OPTIMIZED MCSCF prepared in optimized_mscsf file"

# Handle DMO group  
elif [ "$num" == "6" ]; then
  echo "You have selected 6 - DMO group"
  read -p $'Please input (x) number of DMOs and (y) the initial DMO in format: x y \nYour choice:' ndmo inidmo

  sed -n "/Semi-canonical/,/END/ p" $file > ztemp
  
  norb=`grep '^... 1' ztemp |wc -l`
  
  echo "number of orbitals: $norb"
  
  nline=`wc -l ztemp|cut -d' ' -f1`
  nline_per_orb=`echo "($nline - 3)/$norb"|bc -l|cut -d'.' -f1`
  echo "Number of lines for each orbita: $nline_per_orb"
  
  ini_line=`echo "($inidmo - 1) * $nline_per_orb + 3"|bc -l`
  end_line=`echo " $ini_line + $ndmo * $nline_per_orb - 1"|bc -l`
  
  echo $ini_line $end_line
  
  sed -n "$ini_line, $end_line p" ztemp > ztemp2
  
  # Remove existing dmo.dat
  if [ -f dmo.dat ]
  then
    rm -f dmo.dat
  fi
  
  echo ' $DMO    ' >> dmo.dat
  for idmo in $(seq 1 1 $ndmo)
  do
    ini_line=`echo "($idmo -1) * $nline_per_orb" + 1|bc -l`
    end_line=`echo " $ini_line + $nline_per_orb - 1"|bc -l`
  # echo $idmo $ini_line $end_line
    if [ $idmo -lt 10 ]
    then
      sed -n "$ini_line, $end_line p" ztemp2 |sed -e "s/^../ $idmo/g" >> dmo.dat
    else
      sed -n "$ini_line, $end_line p" ztemp2 |sed -e "s/^../$idmo/g" >> dmo.dat
    fi
  done
  
  echo ' $END' >> dmo.dat
  echo 'DMO group prepared in dmo.dat file'

# Handle REFDET GROUP  
elif [ "$num" == "7" ]; then
  echo "You have selected 7 - REFDET group"
  nstate=`grep '# of states in CI      = ' $file|tail -1|cut -d'=' -f2`

  echo $nstate "states"
  
  sed -n "/DIABATIZATION FOR GMC-QDPT STATES/,/REFDET/ p" $file > ztemp
  
  nstate_plus1=`echo "$nstate + 1"|bc -l`
  echo " State #    $nstate_plus1" >> ztemp
  
  
  # Remove the existing refdet.out
  if [ -f refdet.out ]
  then
    rm -f refdet.out
  fi
  
  # Remove the existing phase.out
  if [ -f phase.out ]
  then
    rm -f phase.out
  fi
  
  echo ' $REFDET     ' >> refdet.out
  echo " "$nstate >> refdet.out
  
  nleft=`expr $nstate % 5`
  nline_for_5=`echo "$nstate / 5"|bc -l`
  #echo $nline_for_5 $nleft
  
  for ist in $(seq 1 1 $nstate)
  do
   ist_plus1=`echo "$ist + 1"|bc -l`
   sed -n "/State..* $ist  Energy/,/State..* $ist_plus1  Energy/ p" ztemp |grep '^ .......[0-9]'|cut -c1-9,22-33 > ztemp2
   ndet=`wc -l ztemp2|cut -d' ' -f1`
   echo " $ist $ndet" >> refdet.out
   cat ztemp2 >> refdet.out
   igroup_5=0
   for jst in $(seq 1 1 $nstate)
   do
    modulo_jst_5=`expr $jst % 5`
  # echo $modulo_jst_5
    elemnt=' 0.00000000E+00'
    if [ $jst -eq $ist ]
    then
      elemnt=' 1.00000000E+00'
    fi
    if [ "$modulo_jst_5" -eq 1 ]
    then
      igroup_5=`echo " $igroup_5 + 1 "|bc -l`
      if [ $jst -ne $nstate ]
      then
        echo -n " $ist  $igroup_5$elemnt" >> phase.out
      else
        echo " $ist  $igroup_5$elemnt" >> phase.out
      fi
    elif [ "$modulo_jst_5" -eq 0 ] || [ $jst -eq $nstate ]
    then
      echo "$elemnt" >> phase.out
    else
      echo -n "$elemnt" >> phase.out
    fi
   done
  done
  
  cat phase.out >> refdet.out
  echo ' $END ' >> refdet.out
  echo "REFDET group prepared in refdet.out file"
else
  echo "Invalid selection."
fi