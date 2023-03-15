#!/usr/bin/env python3

## Generated via machine learning

import sys
import os
import subprocess

if len(sys.argv) == 1:
    print()
    print("dist gamess_hess.out")
    print()
    sys.exit(1)

hessout = sys.argv[1]
filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"

with open(hessout, "r") as f:
    file_data = f.read()

natoms = int(file_data.split(" TOTAL NUMBER OF ATOMS")[1].split("=")[1])
ndim = natoms * 3
print("Dimension of all xyz coordinates:", ndim)
natom = ndim // 3
print("# of atoms:", natom)

ngroup = ndim // 5
nleft = ndim - 5 * ngroup
print(ngroup, nleft)

nrmmod = {}
freqcm = {}

mode_lines = file_data.split("FREQUENCIES IN CM")[1].split("REFERENCE ON SAYVETZ CONDITIONS")[0].split("\n")
mode_data = [line[20:].strip() for line in mode_lines if line.startswith("  ")]
with open("mode.dat", "w") as f:
    f.write("\n".join(mode_data))
freq_lines = file_data.split("FREQUENCIES IN CM")[1].split("REFERENCE ON SAYVETZ CONDITIONS")[0].split("\n")
freq_data = [line[18:].strip() for line in freq_lines if line.startswith("  FREQUENCY:")]
with open("freq.dat", "w") as f:
    f.write("\n".join(freq_data))

for igroup in range(1, ngroup+1):
    iniline = (igroup - 1) * (ndim + 2) + 1
    endline = iniline - 1 + ndim
    print("igroup=", igroup)
    print("iniline=", iniline)
    ixyz = 0
    for line in range(iniline, endline+1):
        ixyz += 1
        print(ixyz)
        for icolumn in range(1, 6):
            imode = (igroup - 1) * 5 + icolumn
            cutini = (icolumn - 1) * 12 + 1
            cutfnl = icolumn * 12
            disp = file_data.split("\n")[line-1][cutini-1:cutfnl]
            nrmmod[(ixyz, imode)] = disp
            print(nrmmod[(ixyz, imode)], disp)
            if ixyz == 1:
                freq = freq_data[igroup-1][cutini-1:cutfnl]
                freqcm[imode] = freq
                print("frequency:", imode, freqcm[imode])
                

# Above code was taken from the bash script for lines 1 to 64
import subprocess

# for the leftover nleft modes
if nleft != 0:
    for igroup in range(nleft, nleft+1):
        iniline = int(ngroup * (ndim + 2) + 1)
        endline = int(iniline - 1 + ndim)
        print("igroup=leftover")
        print("iniline=", iniline)
        ixyz = 0
        for line in range(iniline, endline+1):
            ixyz += 1
            print(ixyz)
            for icolumn in range(1, nleft+1):
                imode = int(ngroup*5 + icolumn)
                cutini = int((icolumn - 1)*12 + 1)
                cutfnl = int(icolumn*12)
                disp = subprocess.check_output(["sed", "-n", f"{line}p", "<", "mode.dat"]).decode()[cutini-1:cutfnl]
                nrmmod[ixyz, imode] = float(disp)
                print(nrmmod[ixyz, imode], disp)

                # Now read in frequencies in cm-1
                if ixyz == 1:
                    freq = subprocess.check_output(["tail", "-1", "freq.dat"]).decode()[cutini-1:cutfnl]
                    freqcm[imode] = freq
                    print("frequency:", imode, freqcm[imode])

# Print all frequencies
for imode in range(1, ndim+1):
    print("frequency:", imode, freqcm[imode], "CM-1")

# read in reference structure
# The ref_structure has to be prepared by human-being and adopts the following format
# N           7.0   0.0000000000  -0.0000000000  -0.1693806842
# H           1.0  -0.4653267700   0.8059696078   0.2564602281
# H           1.0  -0.4653267700  -0.8059696078   0.2564602281
# H           1.0   0.9306535400   0.0000000000   0.2564602281
# The # symbol before each line shall not be in ref_structure.

atmlst = []
chrglst = []
refcoord = []

with open("ref_structure") as f:
    for iatom in range(1, natom+1):
        linecontent = f.readline().strip()
        atmnam = linecontent.split()[0]
        atmlst.append(atmnam)
        chrglst.append(float(linecontent.split()[1]))
        print(atmlst[iatom-1], chrglst[iatom-1])
        for ixyz in range(1, 4):
            icomp = int((iatom - 1)*3 + ixyz)
            ifield = ixyz + 2
            refcoord.append(float(linecontent.split()[ifield]))
            print(refcoord[icomp-1])

# Above code from line 65 till now (121) was taken from the bash script for lines 70 till 144

import numpy as np

# echo ${atmlst[@]}
# Not sure what `atmlst` is, so omitting this line

distcoord_plus = np.zeros(ndim)
distcoord_minus = np.zeros(ndim)
distcoord_plus_x2 = np.zeros(ndim)
distcoord_minus_x2 = np.zeros(ndim)
distcoord_pp = np.zeros((ndim, ndim))
distcoord_pm = np.zeros((ndim, ndim))
distcoord_mp = np.zeros((ndim, ndim))
distcoord_mm = np.zeros((ndim, ndim))

# first copy refcoord to distcoord_plus and distcoord_minus
for icomp in range(ndim):
    distcoord_plus[icomp] = refcoord[icomp]
    distcoord_minus[icomp] = refcoord[icomp]
    # echo ${distcoord[$icomp]}
    print(distcoord_plus[icomp])

# List of modes not considered
modes_excluded = [1, 2, 3, 4, 5, 6]
nexclud = len(modes_excluded)
print("NUMBER OF EXCLUDED MODES:", len(modes_excluded), nexclud)
print("They are modes:", modes_excluded)
nmodes_include = ndim - nexclud

modes_included = {}
icount = 0
for imode in range(1, ndim+1):
    # check whether the mode is considered
    include = 1
    # The modes_excluded array starts counting by 0
    nexclud_end = nexclud - 1
    for iexclud in range(nexclud):
        if imode == modes_excluded[iexclud]:
            include = 0
    # if mode considered
    if include == 1:
        icount += 1
        modes_included[icount] = imode
        print(icount, modes_included[icount])

nmodes_included = icount
print("Number of Modes Included:", nmodes_included)

w
