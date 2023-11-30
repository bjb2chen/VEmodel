#!/bin/bash

# Syntax: dat2vec.x step1.dat

file=$1

# Remove existing vec.dat
if [ -f vec.dat ]
then
  rm -f vec.dat
fi

sed -n "/OPTIMIZED RHF/,/END/ p" $file > vec.dat
