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

print('The following arguments were passed to this ' + str(sys.argv[0])   \
													 + ' program: '       \
													 + str(sys.argv[1:]))

hessout = sys.argv[1]
filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"
target_str = ' TOTAL NUMBER OF ATOMS'

with open(str(hessout), "r") as fp:
	data = fp.read()

target_index = data.find(target_str)
equal_sign_index = data[target_index:].find('=')
natoms_index = data[target_index:].find('\n')

natoms = int(data[target_index+equal_sign_index+1:target_index+natoms_index].strip())
ndim = natoms*3
print(f'Dimension of all xyz coordinates: {ndim}')
natom = ndim/3
print(f'# of atoms: {natom}')
ngroup = ndim//5
nleft = ndim - 5 * ngroup
print(f'{ngroup} {nleft}')

nrmmod = np.zeros(5)
freqcm = np.zeros(5)


print('----------------------------')
#print(f'The location of target string in the output file {hessout} is: {target_index}.')
#print(f'{data[target_index:target_index+100]}')
#print(f'The location of next equal sign is at: {equal_sign_index}')
#print(f'The location of natoms integer is at: {natoms_index}')
#print(f'The string for natoms is: {natoms}')


