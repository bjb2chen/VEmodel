#!/usr/bin/env python3

# shebang for python

import sys
import logging
import numpy as np

try:
	sys.argv[1]
except Exception as e:
	print('dist gamess_hess.out')
	# e.message, e.args

print('The following arguments were passed to this ' + str(sys.argv[0]) + ' program: ' + str(sys.argv[1:]))

hessout = sys.argv[1]
filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"
target_str = ' TOTAL NUMBER OF ATOMS'

with open(str(hessout), "r") as fp:
	data = fp.read()

target_idx = data.find(target_str)
eq_idx = data[target_idx:].find('=')
newline_idx = data[target_idx:].find('\n')

natoms = int(data[target_idx+eq_idx+1:target_idx+newline_idx].strip()) # +1 to go past the = sign.
ndim = natoms*3
print(f'Dimension of all xyz coordinates: {ndim}')
natom = ndim/3
print(f'# of atoms: {natom}')
ngroup = ndim//5
nleft = ndim - 5 * ngroup
print(f'{ngroup} {nleft}')

nrmmod = np.zeros(5)
freqcm = np.zeros(5)

header_str = 'FREQUENCIES IN CM'
footer_str = 'REFERENCE ON SAYVETZ CONDITIONS'

# Define the input file path (replace '$hessout' with the actual file path)
input_file = str(hessout)
print(f'{input_file}')

# Initialize an empty list to store the selected lines
selected_lines = []

# Open the input file and read its contents
with open(input_file, 'r') as file:
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

# Filter and extract lines that start with two digits (grep -A3 '^..[0-9]')
filtered_lines = [line for line in selected_lines if line.strip() and line[0:2].isdigit()]

filtered_set = []
freq_value_set = []
#print(selected_lines)

for idx,modeline in enumerate(selected_lines):
	if len(modeline) > 3 and modeline[2].isdigit():
		# print(idx, modeline)
		# print(selected_lines[idx][19:])
		# print(selected_lines[idx+1][19:])
		# print(selected_lines[idx+2][19:])
		filtered_set.append(selected_lines[idx][19:])
		filtered_set.append(selected_lines[idx+1][19:])
		filtered_set.append(selected_lines[idx+2][19:])

for idx, freqline in enumerate(selected_lines):
	if "FREQUENCY:" in freqline:
		freq_value_set.append(selected_lines[idx][18:])


# def myFilter(lines_to_filter):
# 	for idx, modeline in enumerate(lines_to_filter):
# 		# print(idx, modeline)
# 		if len(modeline) > 3 and modeline[2].isdigit():
# 			# print('here')
# 			filtered_set.append(lines_to_filter.pop(idx)[20:])
# 			filtered_set.append(lines_to_filter.pop(idx)[20:])
# 			filtered_set.append(lines_to_filter.pop(idx)[20:])
# 	return filtered_set

#myFilter(selected_lines)
#print(filtered_set)

# Extract characters from the 21st character to the end of each line (cut -c21-)
# extracted_data = [line[20:] for line in filtered_lines]

# Write the extracted data to 'mode.dat'
with open('oct3_mode.dat', 'w') as output_file:
    output_file.writelines(filtered_set)

with open('oct3_freq.dat', 'w') as output_file:
	output_file.writelines(freq_value_set)


print('----------------------------')
#print(f'The location of target string in the output file {hessout} is: {target_idx}.')
#print(f'{data[target_index:target_index+100]}')
#print(f'The location of next equal sign is at: {equal_sign_index}')
#print(f'The location of natoms integer is at: {natoms_idx}')
#print(f'The string for natoms is: {natoms}')


