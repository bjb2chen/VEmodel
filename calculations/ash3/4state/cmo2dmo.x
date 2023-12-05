#!/bin/bash

# Syntax: cmo2dmo.x gamess.dat num_dmo ini_dmo

file=$1

ndmo=$2

inidmo=$3

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
