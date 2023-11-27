#!/bin/bash

file=$1

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
