#!/bin/bash

# Syntax: combined_script.sh step3.dat

file=$1

# Prompt the user to choose between "OPTIMIZED RHF" and "Semi-canonical MOs"
read -p "Do you want (1) OPTIMIZED RHF or (2) Semi-canonical MOs? " num

# Perform the corresponding operation based on the user's choice
if [ "$num" == "1" ]; then
  # Remove existing vec.dat
  if [ -f vec.dat ]; then
    rm -f vec.dat
  fi
  sed -n "/OPTIMIZED RHF/,/END/ p" $file > vec.dat
  echo "You have selected 1"
  echo "vec.dat is prepared."
elif [ "$num" == "2" ]; then
  # Remove existing active_space_orbs
  if [ -f active_space_orbs ]; then
    rm -f active_space_orbs
  fi
  sed -n "/Semi-canonical MOs/,/END/ p" $file > active_space_orbs
  echo "You have selected 2"
  echo "active_space_orbs is prepared."
else
  echo "Invalid selection."
fi
