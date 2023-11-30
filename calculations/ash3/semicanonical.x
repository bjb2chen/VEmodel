#!/bin/bash

# Syntax: semicanonical.x step3.dat

file=$1

# Remove existing active_space_orbs
if [ -f active_space_orbs ]
then
  rm -f active_space_orbs
fi

sed -n "/Semi-canonical MOs/,/END/ p" $file > active_space_orbs
