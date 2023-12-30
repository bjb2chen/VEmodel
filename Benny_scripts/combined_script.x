#!/bin/bash

# Syntax: combined_script.x GMS.out(or .dat)

file=$1

# Prompt the user to choose between "OPTIMIZED RHF" and "Semi-canonical MOs"
read -p "Do you want (1) OPTIMIZED RHF \n (2) OPTIMIZED ROHF \n (3) Semi-canonical MOs? \n (4) DMO group \n (5) REFDET GROUP" num

# Perform the corresponding operation based on the user's choice
if [ "$num" == "1" ]; then
  # Remove existing vec.dat
  if [ -f vec.dat ]; then
    rm -f vec.dat
  fi
  sed -n "/OPTIMIZED RHF/,/END/ p" $file > vec.dat
  echo "You have selected 1"
  echo "OPTIMIZED RHF orbitals prepared in vec.dat file"
elif [ "$num" == "2" ]; then
  # Remove existing vec.dat
  if [ -f vec.dat ]; then
    rm -f vec.dat
  fi
  sed -n "/OPTIMIZED ROHF/,/END/ p" $file > vec.dat
  echo "You have selected 2"
  echo "OPTIMIZED ROHF orbitals prepared in vec.dat file"
elif [ "$num" == "3" ]; then
  # Remove existing active_space_orbs
  if [ -f active_space_orbs ]; then
    rm -f active_space_orbs
  fi
  sed -n "/Semi-canonical MOs/,/END/ p" $file > active_space_orbs
  echo "You have selected 2"
  echo "SEmi-canonical MOs prepared in active_space_orbs  file"
elif [ "$num" == "4" ]; then
  read -p "Please input (x) number of DMOs and (y) the initial DMO in format: x y" ndmo inidmo

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
elif [ "$num" == "5" ]; then
  nstate=`grep '# of states in CI      = ' $file|tail -1|cut -d'=' -f2`

  echo $nstate
  
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
else
  echo "Invalid selection."
fi