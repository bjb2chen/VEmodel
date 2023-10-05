#!/usr/bin/env python3

# system imports

import sys
import logging
import numpy as np

# third party imports

# local imports

try:
	sys.argv[1]
except Exception as e:
	print('dist gamess_hess.out')

hessout = sys.argv[1]
filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"

# Extract the number of atoms
with open(hessout, 'r') as hess_file:
    for line in hess_file:
        if ' TOTAL NUMBER OF ATOMS' in line:
            natoms = int(line.split('=')[1])
            ndim = natoms * 3
            break

print("Dimension of all xyz coordinates:", ndim)
natom = ndim // 3
print("# of atoms:", natom)

ngroup = ndim // 5
nleft = ndim % 5
print(ngroup, nleft)

# Create dictionaries to store data
nrmmod = {}
freqcm = {}

# Initialize an empty list to store the selected lines
selected_lines = []
filtered_set = []
freq_value_set = []

# Extract lines from given output file
with open(str(hessout), 'r') as file:
    # Flag to indicate whether to collect lines between patterns
    collecting = False
    # Iterate through each line in the file
    for line in file:
        # Check if the line contains the starting pattern
        if "FREQUENCIES IN CM" in line:
            collecting = True
            # Add the current line to the selected lines
            selected_lines.append(line)
        # Check if the line contains the ending pattern
        elif "REFERENCE ON SAYVETZ CONDITIONS" in line:
            collecting = False
        # If we're in the collection mode, add the line to selected lines
        elif collecting:
            selected_lines.append(line)

# Filter and extract lines that have third character starting as digit (grep -A3 '^..[0-9]')
for idx,modeline in enumerate(selected_lines):
	if len(modeline) > 3 and modeline[2].isdigit():    # IndexError if line is less than 3 characters
		filtered_set.append(selected_lines[idx][20:])  # Slice from 21st character to end of line
		filtered_set.append(selected_lines[idx+1][20:])
		filtered_set.append(selected_lines[idx+2][20:])

for idx, freqline in enumerate(selected_lines):
	if "FREQUENCY:" in freqline:
		freq_value_set.append(selected_lines[idx][18:])

#filtered_lines = [line for line in selected_lines if line.strip() and line[0:2].isdigit()]

# Write the extracted data to 'mode.dat'
with open('oct3_mode.dat', 'w') as output_file:
    output_file.writelines(filtered_set)

with open('oct3_freq.dat', 'w') as output_file:
	output_file.writelines(freq_value_set)

for igroup in range(1, ngroup + 1, 1):
	# No double whitespace line between groups, so no plus 2 for iniline
    iniline = (igroup - 1) * ndim + 1
    endline = iniline + ndim - 1
    print("igroup =", igroup)
    print("iniline =", iniline)
    ixyz = 0
    for line in range(iniline, endline + 1, 1):
        ixyz += 1
        print(ixyz)
        for icolumn in range(1, 6, 1):
            imode = (igroup - 1) * 5 + icolumn
            print(ixyz, imode, end=" ")
            cutini = (icolumn - 1) * 12
            cutfnl = icolumn * 12
            with open("oct3_mode.dat", "r") as mode_file:
                lines = mode_file.readlines()
                disp = lines[line - 1][cutini:cutfnl]
                nrmmod[ixyz, imode] = disp
                print(nrmmod[ixyz, imode], disp)

            if ixyz == 1:
                cutini = (icolumn - 1) * 12
                cutfnl = icolumn * 12
                with open("oct3_freq.dat", "r") as freq_file:
                    lines = freq_file.readlines()
                    freq = lines[igroup - 1][cutini:cutfnl].lstrip()
                    freqcm[imode] = freq
                    print("frequency:", imode, freqcm[imode])

# For the leftover nleft modes
if nleft != 0:
    for igroup in range(nleft, nleft + 1):
    	# No double whitespace line between groups, so no plus 2 for iniline
        iniline = ngroup * ndim + 1 
        endline = iniline + ndim - 1
        print("igroup=leftover")
        print("iniline =", iniline)
        ixyz = 0
        for line in range(iniline, endline + 1):
            ixyz += 1
            print(ixyz)
            for icolumn in range(1, nleft + 1):
                imode = ngroup * 5 + icolumn
                print(ixyz, imode, end=" ")
                cutini = (icolumn - 1) * 12
                cutfnl = icolumn * 12
                with open("oct3_mode.dat", "r") as mode_file:
                    lines = mode_file.readlines()
                    disp = lines[line - 1][cutini:cutfnl]
                nrmmod[ixyz, imode] = disp
                print(nrmmod[ixyz, imode], disp)

                if ixyz == 1:
                    cutini = (icolumn - 1) * 12
                    cutfnl = icolumn * 12
                    with open("oct3_freq.dat", "r") as freq_file:
                        lines = freq_file.readlines()
                        freq = lines[-1][cutini:cutfnl]
                    freqcm[imode] = freq
                    print("frequency:", imode, freqcm[imode])

# Print all frequencies
for imode in range(1, ndim + 1):
    print("frequency:", imode, freqcm[imode].lstrip(), "CM-1")

#read in reference structure
#The ref_structure has to be prepared by human-being and adopts the following format
# N           7.0   0.0000000000  -0.0000000000  -0.1693806842
# H           1.0  -0.4653267700   0.8059696078   0.2564602281
# H           1.0  -0.4653267700  -0.8059696078   0.2564602281
# H           1.0   0.9306535400   0.0000000000   0.2564602281
#The # symbol before each line shall not be in ref_structure.

atmlst = {}
chrglst = {}
refcoord = {}

# for iatom in range(1, natom, 1):
# 	linecontent = 

print('----------------------------')
# print('The following arguments were passed to this ' + str(sys.argv[0]) + ' program: ' + str(sys.argv[1:]))

#print(f'The location of target string in the output file {hessout} is: {target_idx}.')
#print(f'{data[target_index:target_index+100]}')
#print(f'The location of next equal sign is at: {equal_sign_index}')
#print(f'The location of natoms integer is at: {natoms_idx}')
#print(f'The string for natoms is: {natoms}')


# target_idx = data.find(target_str)
# eq_idx = data[target_idx:].find('=')
# newline_idx = data[target_idx:].find('\n')

# natoms = int(data[target_idx+eq_idx+1:target_idx+newline_idx].strip()) # +1 to go past the = sign.
# ndim = natoms*3
# print(f'Dimension of all xyz coordinates: {ndim}')
# natom = ndim/3
# print(f'# of atoms: {natom}')
# ngroup = ndim//5
# nleft = ndim - 5 * ngroup
# print(f'{ngroup} {nleft}')

# header_str = 'FREQUENCIES IN CM'
# footer_str = 'REFERENCE ON SAYVETZ CONDITIONS'

		# print(idx, modeline)
		# print(selected_lines[idx][19:])
		# print(selected_lines[idx+1][19:])
		# print(selected_lines[idx+2][19:])

#print(selected_lines)
#myFilter(selected_lines)
#print(filtered_set)

