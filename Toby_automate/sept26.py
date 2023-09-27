#!/usr/bin/env python3

# shebang for python

import sys
import logging

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

## natoms=`grep ' TOTAL NUMBER OF ATOMS' $hessout|cut -d'=' -f2`