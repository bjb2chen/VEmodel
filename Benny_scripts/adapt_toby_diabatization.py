import sys
import pprint
import subprocess
import os
import shutil
import re
import json

# Function to get the number of atoms from thFe hessout file
def get_number_of_atoms(hessout):
    with open(hessout, 'r') as hess_file:
        for line in hess_file:
            if ' TOTAL NUMBER OF ATOMS' in line:
                natoms = int(line.split('=')[1])
                return natoms

# Function to extract lines between patterns in a file
def extract_lines_between_patterns(filename, start_pattern, end_pattern):
    selected_lines = []
    collecting = False

    with open(filename, 'r') as file:
        for line in file:
            if start_pattern in line:
                collecting = True
                selected_lines.append(line)
            elif end_pattern in line:
                collecting = False
            elif collecting:
                selected_lines.append(line)

    return selected_lines

def read_reference_structure(file_path):
    ##### SAMPLE REF STRUCT #######
    #read in reference structure
    #The ref_structure has to be prepared by human-being and adopts the following format
    # N           7.0   0.0000000000  -0.0000000000  -0.1693806842
    # H           1.0  -0.4653267700   0.8059696078   0.2564602281
    # H           1.0  -0.4653267700  -0.8059696078   0.2564602281
    # H           1.0   0.9306535400   0.0000000000   0.2564602281
    #The # symbol before each line shall not be in ref_structure.
    ###############################
    atmlst = {}
    chrglst = {}
    refcoord = {}

    with open(file_path, "r") as struct_file:
        lines = struct_file.readlines()

        for iatom, line in enumerate(lines):
            parts = line.split()
            atmnam = parts[0]
            chrg = parts[1]
            coords = parts[2:]

            atmlst[iatom + 1] = atmnam
            chrglst[iatom + 1] = chrg

            print(atmlst[iatom + 1], chrglst[iatom + 1])

            for ixyz, coord in enumerate(coords):
                icomp = (iatom * 3) + ixyz + 1
                refcoord[icomp] = float(coord)
                print(refcoord[icomp])

    return atmlst, chrglst, refcoord

def main():
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <hessout_file>")
        sys.exit(1)

    hessout = sys.argv[1]

    # ...

if __name__ == "__main__":
    main()