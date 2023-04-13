#!/usr/bin/env python

import sys
import subprocess
import pandas as pd
import numpy as np
import os
import math
import re


# Define function to run bash command and return stdout as string
def run_bash(cmd):
    return subprocess.check_output(cmd, shell=True, universal_newlines=True)

if len(sys.argv) == 1:
    print("\n")
    print("dist gamess_hess.out")
    print("\n")
    sys.exit(1)

hessout = sys.argv[1]
filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"

# Get the number of atoms from gamess_hess.out
natoms = int(run_bash(f"grep ' TOTAL NUMBER OF ATOMS' {hessout} | cut -d'=' -f2"))
ndim = natoms * 3
print("Dimension of all xyz coordinates:", ndim)
natom = ndim // 3
print("# of atoms:", natom)

ngroup = ndim // 5
nleft = ndim - 5 * ngroup
print(ngroup, nleft)

nrmmod = {}
freqcm = {}

# Get mode.dat and freq.dat files
mode_dat = run_bash(f"sed -n '/FREQUENCIES IN CM/,/REFERENCE ON SAYVETZ CONDITIONS/p' {hessout} | grep -A3 '^..[0-9]' | cut -c21-")
freq_dat = run_bash(f"sed -n '/FREQUENCIES IN CM/,/REFERENCE ON SAYVETZ CONDITIONS/p' {hessout} | grep 'FREQUENCY:' | cut -c19-")

# Loop through each group of modes
for igroup in range(1, ngroup + 1):
    iniline = (igroup - 1) * (ndim + 2) + 1
    endline = iniline - 1 + ndim
    print("igroup=", igroup)
    print("iniline=", iniline)
    ixyz = 0
    # Loop through each line of modes in the group
    for line in range(iniline, endline + 1):
        ixyz += 1
        print(ixyz)
        # Loop through each column of modes in the line
        for icolumn in range(1, 6):
            imode = (igroup - 1) * 5 + icolumn
            cutini = (icolumn - 1) * 12 + 1
            cutfnl = icolumn * 12
            disp = mode_dat[line - 1][cutini - 1:cutfnl]
            nrmmod[(ixyz, imode)] = disp
            print(nrmmod[(ixyz, imode)], disp)

            # Get frequency for the first column of modes in the line
            if ixyz == 1:
                freq = freq_dat.splitlines()[igroup - 1][cutini - 1:cutfnl]
                freqcm[imode] = freq
                print("frequency:", imode, freqcm[imode])


# Above lines (1-61) are from the bash script lines 1-64

if len(sys.argv) == 1:
    print("\ndist gamess_hess.out\n")
    sys.exit(1)

hessout = sys.argv[1]
filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"

# run command and store output in a variable
command = f"grep ' TOTAL NUMBER OF ATOMS' {hessout} | cut -d'=' -f2"
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output, error = process.communicate()
natoms = int(output.decode().strip())

ndim = natoms * 3
print(f"Dimension of all xyz coordinates: {ndim}")
natom = ndim // 3
print(f"# of atoms: {natom}")

ngroup = ndim // 5
nleft = ndim - 5 * ngroup
print(f"{ngroup} {nleft}")

nrmmod = {}
freqcm = {}

# write output to a file
with open("mode.dat", "w") as mode_file, open("freq.dat", "w") as freq_file:
    command = f"sed -n '/FREQUENCIES IN CM/,/REFERENCE ON SAYVETZ CONDITIONS/p' {hessout} | grep -A3 '^..[0-9]' | cut -c21-"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    mode_file.write(output.decode())

    command = f"sed -n '/FREQUENCIES IN CM/,/REFERENCE ON SAYVETZ CONDITIONS/p' {hessout} | grep 'FREQUENCY:' | cut -c19-"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    freq_file.write(output.decode())

# read from file and store data in dictionaries
with open("mode.dat", "r") as mode_file, open("freq.dat", "r") as freq_file:
    for igroup in range(1, ngroup + 1):
        iniline = (igroup - 1) * (ndim + 2) + 1
        endline = iniline - 1 + ndim
        print(f"igroup={igroup}")
        print(f"iniline={iniline}")
        ixyz = 0
        for line in range(iniline, endline + 1):
            ixyz += 1
            print(ixyz)
            for icolumn in range(1, 6):
                imode = (igroup - 1) * 5 + icolumn
                print(ixyz, imode, end=" ")
                cutini = (icolumn - 1) * 12 + 1
                cutfnl = icolumn * 12
                disp = mode_file.readlines()[line - 1][cutini - 1:cutfnl].strip()
                nrmmod[ixyz, imode] = disp
                print(nrmmod[ixyz, imode], disp)

                # read frequencies
                if ixyz == 1:
                    freq = freq_file.readlines()[igroup - 1][cutini - 1:cutfnl].strip()
                    freqcm[imode] = freq
                    print(f"frequency: {imode} {freqcm[imode]}")

for imode in range(1, ndim+1):
    print("frequency: ", imode, freqcm[imode], "CM-1")

#lines 65-115

# read in reference structure from file
ref_structure = pd.read_csv('ref_structure.txt', delim_whitespace=True, comment='#', header=None)

# extract atomic symbols, charges, and coordinates
atmlst = ref_structure.iloc[:, 0].tolist()
chrglst = ref_structure.iloc[:, 1].tolist()
refcoord = ref_structure.iloc[:, 2:].values.flatten()

natom = len(atmlst)
ndim = int(len(refcoord)/natom)

# print out atomic symbols and charges
for iatom in range(natom):
    print(atmlst[iatom], chrglst[iatom])

    # print out coordinates for each atom
    for ixyz in range(3):
        icomp = (iatom * 3) + ixyz
        print(refcoord[icomp])

# create arrays for various distance calculations
distcoord_plus = np.copy(refcoord)
distcoord_minus = np.copy(refcoord)
distcoord_plus_x2 = np.copy(refcoord)
distcoord_minus_x2 = np.copy(refcoord)
distcoord_pp = np.copy(refcoord)
distcoord_pm = np.copy(refcoord)
distcoord_mp = np.copy(refcoord)
distcoord_mm = np.copy(refcoord)

# list of modes not considered
modes_excluded = [1, 2, 3, 4, 5, 6]
nexclud = len(modes_excluded)
print("NUMBER OF EXCLUDED MODES:", nexclud)
print("They are modes:", modes_excluded)

# list of modes included
modes_included = {}
icount = 0
for imode in range(1, ndim+1):

    # check whether the mode is considered
    include = True
    for iexclud in range(nexclud):
        if imode == modes_excluded[iexclud]:
            include = False

    # if mode considered
    if include:
        icount += 1
        modes_included[icount] = imode
        print(icount, modes_included[icount])

nmodes_included = icount
print("Number of Modes Included:", nmodes_included)

#Above lines 116-196

# Do diabatization calculation at the reference nondistorted structure.
# This calculation shall be a repetition of a calculation in preparing temp.inp
output = subprocess.getoutput(f'grep grace {filnam}_refG.out')
if 'grace' not in output:
    print("Run calculation at the undistorted reference structure")
    with open(f'{filnam}_refG.inp', 'w') as f:
        f.write(open('temp.inp').read())
        f.write(open('ref_structure').read())
        f.write(' $END    \n')
    # You may need to edit the following submission line based on your computer's configuration
    os.system(f'subgam.diab {filnam}_refG.inp 4 0 1')
else:
    print("Calculation at the reference structure is done.")

# Loop over all considered modes and do + and - displacements
# set step size in reduced dimensionless coordinates
qsize = 0.05
# Set conversion constants
ha2ev = 27.2113961318
wn2ev = 0.000123981
wn2eh = 0.00000455633
ang2br = 1.889725989
amu2me = 1822.88839

for kmode in range(1, nmodes_included+1):
    imode = modes_included[kmode-1]

    # Convert the reduced dimensionless qsize to the actual rsize in sqrt(amu)*Angs unit
    omga = freqcm[imode]
    rsize = qsize / (math.sqrt(amu2me)*ang2br*math.sqrt(omga*wn2eh))
    print(imode, omga, rsize)

    # Loop over components
    for icomp in range(1, ndim+1):
        coord_disp_plus = refcoord[icomp] + rsize*nrmmod[icomp-1, imode]
        coord_disp_minus = refcoord[icomp] - rsize*nrmmod[icomp-1, imode]
        distcoord_plus[icomp-1] = coord_disp_plus
        distcoord_minus[icomp-1] = coord_disp_minus
        coord_disp_plusx2 = refcoord[icomp] + 2.0*rsize*nrmmod[icomp-1, imode]
        coord_disp_minusx2 = refcoord[icomp] - 2.0*rsize*nrmmod[icomp-1, imode]
        distcoord_plus_x2[icomp-1] = coord_disp_plusx2
        distcoord_minus_x2[icomp-1] = coord_disp_minusx2
        print(imode, icomp, refcoord[icomp], nrmmod[icomp-1, imode], coord_disp_plus, coord_disp_minus, distcoord_plus[icomp-1], distcoord_minus[icomp-1])

    # delete existing dist_structure_plus and dist_structure_minus
    if os.path.isfile('dist_structure_plus'):
        os.remove('dist_structure_plus')
    if os.path.isfile('dist_structure_minus'):
        os.remove('dist_structure_minus')
    # delete existing dist_structure_plusx2 and dist_structure_minusx2
    if os.path.isfile('dist_structure_plusx2'):
        os.remove('dist_structure_plusx2')
    if os.path.isfile('dist_structure_minusx2'):
        os.remove('dist_structure_minusx2')


#lines 199-265

# print the distorted structure
for iatom in range(1, natom+1):
    with open('dist_structure_plus', 'a') as fplus, open('dist_structure_minus', 'a') as fminus, open('dist_structure_plusx2', 'a') as fplusx2, open('dist_structure_minusx2', 'a') as fminusx2:
        fplus.write(f'{atmlst[iatom]} {chrglst[iatom]} ')
        fminus.write(f'{atmlst[iatom]} {chrglst[iatom]} ')
        fplusx2.write(f'{atmlst[iatom]} {chrglst[iatom]} ')
        fminusx2.write(f'{atmlst[iatom]} {chrglst[iatom]} ')
        for ixyz in range(1, 4):
            icomp = (iatom - 1) * 3 + ixyz
            fplus.write(f'{distcoord_plus[icomp]} ')
            fminus.write(f'{distcoord_minus[icomp]} ')
            fplusx2.write(f'{distcoord_plus_x2[icomp]} ')
            fminusx2.write(f'{distcoord_minus_x2[icomp]} ')
        fplus.write('\n')
        fminus.write('\n')
        fplusx2.write('\n')
        fminusx2.write('\n')
        
# The temp.inp template input file has to be prepared by human-being
# I am just directing the distorted structures into the temp.inp to generate the actual input for gamess
import shutil
shutil.copyfile('temp.inp', f'{filnam}_mode{imode}+{qsize}.inp')
with open(f'{filnam}_mode{imode}+{qsize}.inp', 'a') as finp:
    with open('dist_structure_plus', 'r') as fplus:
        finp.write(fplus.read())
        finp.write(' $END ')
        
shutil.copyfile('temp.inp', f'{filnam}_mode{imode}-{qsize}.inp')
with open(f'{filnam}_mode{imode}-{qsize}.inp', 'a') as finp:
    with open('dist_structure_minus', 'r') as fminus:
        finp.write(fminus.read())
        finp.write(' $END ')
        
shutil.copyfile('temp.inp', f'{filnam}_mode{imode}+{qsize}x2.inp')
with open(f'{filnam}_mode{imode}+{qsize}x2.inp', 'a') as finp:
    with open('dist_structure_plusx2', 'r') as fplusx2:
        finp.write(fplusx2.read())
        finp.write(' $END ')
        
shutil.copyfile('temp.inp', f'{filnam}_mode{imode}-{qsize}x2.inp')
with open(f'{filnam}_mode{imode}-{qsize}x2.inp', 'a') as finp:
    with open('dist_structure_minusx2', 'r') as fminusx2:
        finp.write(fminusx2.read())
        finp.write(' $END ')

#Above lines 267-302


# submit gamess diabatization calculations for the + and - inputs
# Check whether the calculation is done already.
if b'grace' not in subprocess.check_output(['grep', 'grace', f'{filnam}_mode{imode}_{qsize}.out']):
    print(f"run calculations for {filnam}_mode{imode}_{qsize}")
    # You may need to edit the following submission line based on your computer's configuration
    subprocess.run(['subgam.diab', f'{filnam}_mode{imode}_{qsize}.inp', '4', '0', '1'])
else:
    print(f"{filnam}_mode{imode}_{qsize} is done")

# Check whether the calculation is done already.
if b'grace' not in subprocess.check_output(['grep', 'grace', f'{filnam}_mode{imode}_{-qsize}.out']):
    print(f"run calculations for {filnam}_mode{imode}_{-qsize}")
    # You may need to edit the following submission line based on your computer's configuration
    subprocess.run(['subgam.diab', f'{filnam}_mode{imode}_{-qsize}.inp', '4', '0', '1'])
else:
    print(f"{filnam}_mode{imode}_{-qsize} is done")

if b'grace' not in subprocess.check_output(['grep', 'grace', f'{filnam}_mode{imode}_{qsize}x2.out']):
    print(f"run calculations for {filnam}_mode{imode}_{qsize}x2")
    # You may need to edit the following submission line based on your computer's configuration
    subprocess.run(['subgam.diab', f'{filnam}_mode{imode}_{qsize}x2.inp', '4', '0', '1'])
else:
    print(f"{filnam}_mode{imode}_{qsize}x2 is done")

if b'grace' not in subprocess.check_output(['grep', 'grace', f'{filnam}_mode{imode}_{-qsize}x2.out']):
    print(f"run calculations for {filnam}_mode{imode}_{-qsize}x2")
    # You may need to edit the following submission line based on your computer's configuration
    subprocess.run(['subgam.diab', f'{filnam}_mode{imode}_{-qsize}x2.inp', '4', '0', '1'])
else:
    print(f"{filnam}_mode{imode}_{-qsize}x2 is done")


#Above lines 304-342

# 2D distortion to get bilinear vibronic coupling
lmode_last = kmode - 1
for lmode in range(1, lmode_last+1):
    jmode = modes_included[lmode]
    
    # Convert the reduced dimensionless qsize to the actual rsizep in sqrt(amu)*Angs unit for jmode
    omgap = freqcm[jmode]
    rsizep = qsize / (np.sqrt(amu2me)*ang2br*np.sqrt(omgap*wn2eh))
    print(imode, jmode, rsize, rsizep)
    
    # Loop over components
    for icomp in range(1, ndim+1):
        coord_disp_pp = distcoord_plus[icomp] + rsizep * nrmmod[icomp, jmode]
        coord_disp_pm = distcoord_plus[icomp] - rsizep * nrmmod[icomp, jmode]
        coord_disp_mp = distcoord_minus[icomp] + rsizep * nrmmod[icomp, jmode]
        coord_disp_mm = distcoord_minus[icomp] - rsizep * nrmmod[icomp, jmode]
        distcoord_pp[icomp] = coord_disp_pp
        distcoord_pm[icomp] = coord_disp_pm
        distcoord_mp[icomp] = coord_disp_mp
        distcoord_mm[icomp] = coord_disp_mm
    
    # Delete existing dist_structure_pp, _pm, _mp, _mm
    if os.path.exists("dist_structure_pp"):
        os.remove("dist_structure_pp")
    if os.path.exists("dist_structure_pm"):
        os.remove("dist_structure_pm")
    if os.path.exists("dist_structure_mp"):
        os.remove("dist_structure_mp")
    if os.path.exists("dist_structure_mm"):
        os.remove("dist_structure_mm")
    
    # Print the distorted structure
    with open("dist_structure_pp", "w") as f_pp, \
         open("dist_structure_pm", "w") as f_pm, \
         open("dist_structure_mp", "w") as f_mp, \
         open("dist_structure_mm", "w") as f_mm:
        for iatom in range(1, natom+1):
            f_pp.write(atmlst[iatom] + " " + str(chrglst[iatom]) + " ")
            f_pm.write(atmlst[iatom] + " " + str(chrglst[iatom]) + " ")
            f_mp.write(atmlst[iatom] + " " + str(chrglst[iatom]) + " ")
            f_mm.write(atmlst[iatom] + " " + str(chrglst[iatom]) + " ")
            for ixyz in range(1, 4):
                icomp = (iatom - 1)*3 + ixyz
                f_pp.write(str(distcoord_pp[icomp]) + " ")
                f_pm.write(str(distcoord_pm[icomp]) + " ")
                f_mp.write(str(distcoord_mp[icomp]) + " ")
                f_mm.write(str(distcoord_mm[icomp]) + " ")
            f_pp.write("\n")
            f_pm.write("\n")
            f_mp.write("\n")
            f_mm.write("\n")

#Above lines 342-404

import shutil

# Create the four output filenames
out_files = [
    f"{filnam}_mode{imode}+{qsize}_mode{jmode}+{qsize}.inp",
    f"{filnam}_mode{imode}+{qsize}_mode{jmode}-{qsize}.inp",
    f"{filnam}_mode{imode}-{qsize}_mode{jmode}+{qsize}.inp",
    f"{filnam}_mode{imode}-{qsize}_mode{jmode}-{qsize}.inp",
]

# Define the four input file and output file pairs and their respective appended text
file_pairs = [
    ("temp.inp", "dist_structure_pp", out_files[0]),
    ("temp.inp", "dist_structure_pm", out_files[1]),
    ("temp.inp", "dist_structure_mp", out_files[2]),
    ("temp.inp", "dist_structure_mm", out_files[3]),
]

# Copy the input file and append text to create each output file
for i, (input_file, append_file, output_file) in enumerate(file_pairs):
    shutil.copy(input_file, output_file)
    with open(append_file, "r") as f:
        append_text = f.read()
    with open(output_file, "a") as f:
        f.write(append_text)
        f.write(" $END ")


#Above lines 406-420

import subprocess

filnam = "filename"  # replace with actual filename
qsize = "size"  # replace with actual size

for imode in range(1, 5):
    for jmode in range(1, 5):
        # Check whether the calculation is done already.
        if subprocess.call(['grep', 'grace', f'{filnam}_mode{imode}+{qsize}_mode{jmode}+{qsize}.out']) != 0:
            print(f'run calculations for {filnam}_mode{imode}+{qsize}_mode{jmode}+{qsize}')
            # You may need to edit the following submission line based on your computer's configuration
            subprocess.run(['subgam.diab', f'{filnam}_mode{imode}+{qsize}_mode{jmode}+{qsize}.inp', '4', '0', '1'])
        else:
            print(f'{filnam}_mode{imode}+{qsize}_mode{jmode}+{qsize} is done')
        # Check whether the calculation is done already.
        if subprocess.call(['grep', 'grace', f'{filnam}_mode{imode}+{qsize}_mode{jmode}-{qsize}.out']) != 0:
            print(f'run calculations for {filnam}_mode{imode}+{qsize}_mode{jmode}-{qsize}')
            # You may need to edit the following submission line based on your computer's configuration
            subprocess.run(['subgam.diab', f'{filnam}_mode{imode}+{qsize}_mode{jmode}-{qsize}.inp', '4', '0', '1'])
        else:
            print(f'{filnam}_mode{imode}+{qsize}_mode{jmode}-{qsize} is done')
        # Check whether the calculation is done already.
        if subprocess.call(['grep', 'grace', f'{filnam}_mode{imode}-{qsize}_mode{jmode}+{qsize}.out']) != 0:
            print(f'run calculations for {filnam}_mode{imode}-{qsize}_mode{jmode}+{qsize}')
            # You may need to edit the following submission line based on your computer's configuration
            subprocess.run(['subgam.diab', f'{filnam}_mode{imode}-{qsize}_mode{jmode}+{qsize}.inp', '4', '0', '1'])
        else:
            print(f'{filnam}_mode{imode}-{qsize}_mode{jmode}+{qsize} is done')
        # Check whether the calculation is done already.
        if subprocess.call(['grep', 'grace', f'{filnam}_mode{imode}-{qsize}_mode{jmode}-{qsize}.out']) != 0:
            print(f'run calculations for {filnam}_mode{imode}-{qsize}_mode{jmode}-{qsize}')
            # You may need to edit the following submission line based on your computer's configuration
            subprocess.run(['subgam.diab', f'{filnam}_mode{imode}-{qsize}_mode{jmode}-{qsize}.inp', '4', '0', '1'])
        else:
            print(f'{filnam}_mode{imode}-{qsize}_mode{jmode}-{qsize} is done')


#Above lines 420-467


if os.path.isfile('mctdh.op'):
    os.remove('mctdh.op')

# Heading for mctdh.op
with open('mctdh.op', 'a') as f:
    f.write('OP_DEFINE-SECTION\n')
    f.write('title\n')
    nstate = int(os.popen(f'grep \'# of states in CI      = \' {filnam}_refG.out|tail -1|cut -d\'=\' -f2').read())
    f.write(f'{filnam} {nstate} states + {nmodes_include} modes\n')
    f.write('end-title\n')
    f.write('end-op_define-section\n')
    f.write('\nPARAMETER-SECTION\n\n')

#Above lines 473-488

import os
import re

# Check if mctdh.op file exists and remove it if it does
if os.path.isfile("mctdh.op"):
    os.remove("mctdh.op")

# Write heading for mctdh.op
with open("mctdh.op", "a") as f:
    f.write("OP_DEFINE-SECTION\n")
    f.write("title\n")
    nstate = int(re.findall("# of states in CI\s+=\s+(\d+)", open(f"{filnam}_refG.out").read())[-1])
    f.write(f"{filnam} {nstate} states + {nmodes_include} modes\n")
    f.write("end-title\n")
    f.write("end-op_define-section\n\n")
    f.write("PARAMETER-SECTION\n\n")

# Extract diabatic Hamiltonian in the reference structure
if os.path.isfile(f"{filnam}_refG.out"):
    with open("mctdh.op", "a") as f:
        f.write("#Diagonal and Off-diagonal diabatic Hamitonian elements at reference structure\n")
        for ist in range(1, nstate+1):
            Ediab = float(re.findall(f"STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=\s+(\S+)", open(f"{filnam}_refG.out").read())[-1])
            f.write(f"v{ist} = {Ediab}, ev\n")
            for jst in range(1, ist):
                Coup_ev = float(re.findall(f"STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING\s+(\S+)", open(f"{filnam}_refG.out").read())[-1])
                f.write(f"v{jst}{ist} = {Coup_ev}, ev\n")
            f.write("\n")
        f.write("\n")
else:
    print(f"Skip extracting Hamiltonians from the non-existing {filnam}_refG.out")


#Above lines 490-510

nmodes_included = 10
modes_included = range(1, nmodes_included + 1)
freqcm = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
wn2ev = 0.000123984

with open("mctdh.op", "w") as mctdh_file:
    for kmode in modes_included:
        imode = modes_included[kmode - 1]

        vibron_ev = freqcm[imode - 1] * wn2ev
        mctdh_file.write("#Parameters for mode" + str(imode) + "\n")
        mctdh_file.write("#Vibron:\n")
        mctdh_file.write("w_m" + str(imode) + " = " + str(vibron_ev) + ", ev\n")
        mctdh_file.write("\n")
        mctdh_file.write("#Linear and quadratic diagonal and off-diagonal vibronic coupling constants:\n")
        grace_output_plus = subprocess.run(["grep", "grace", "filename_mode" + str(imode) + "_+" + str(qsize) + ".out"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        grace_code_plus = grace_output_plus.returncode
        grace_output_minus = subprocess.run(["grep", "grace", "filename_mode" + str(imode) + "_-" + str(qsize) + ".out"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        grace_code_minus = grace_output_minus.returncode
        grace_output_plusx2 = subprocess.run(["grep", "grace", "filename_mode" + str(imode) + "_+" + str(qsize) + "x2.out"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        grace_code_plusx2 = grace_output_plusx2.returncode
        grace_output_minusx2 = subprocess.run(["grep", "grace", "filename_mode" + str(imode) + "_-" + str(qsize) + "x2.out"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        grace_code_minusx2 = grace_output_minusx2.returncode
        print(imode, grace_code_plus, grace_code_minus)

        if grace_code_plus == 0 and grace_code_minus == 0 and grace_code_plusx2 == 0 and grace_code_minusx2 == 0:
            print("good to extract")
            # nstate = subprocess.run(["grep", "# of states in CI      = ", "filename_mode" + str(imode) + "_+" + str(qsize) + ".out"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.strip().split("=")[1]
            # print(nstate)

#Above lines 513-540

#Extract the diagonal and off-diagonal vibronic coupling
for ist in range(1, nstate+1):
    Ediab_au_plus = float(subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_mode{imode}_+{qsize}.out']).decode('utf-8').split('\n')[-2][43:61])
    Ediab_au_plusx2 = float(subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_mode{imode}_+{qsize}x2.out']).decode('utf-8').split('\n')[-2][43:61])
    Ediab_au_minus = float(subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_mode{imode}_-{qsize}.out']).decode('utf-8').split('\n')[-2][43:61])
    Ediab_au_minusx2 = float(subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_mode{imode}_-{qsize}x2.out']).decode('utf-8').split('\n')[-2][43:61])
    Ediab_au_0 = float(subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_refG.out']).decode('utf-8').split('\n')[-2][43:61])
    linear_diag_ev = ((Ediab_au_plus - Ediab_au_minus) * ha2ev) / (2 * qsize)
    quadratic_diag_ev = ((Ediab_au_plusx2 + Ediab_au_minusx2 - 2.0 * Ediab_au_0) * ha2ev) / (4.0 * qsize * qsize)
    #We only view the difference between the actual force constant and the vibron as the quadratic diagonal coupling for the diabatic state.
    #This is because we add the vibron to the diagonal Hamiltonian matrix elements for all diabats.
    #So, the quadratic diagonal vibronic coupling is actually the correction to the vibron of the normal mode as a force constant.
    quadratic_diag_ev = quadratic_diag_ev - vibron_ev
    print(f'{ist} {linear_diag_ev} {quadratic_diag_ev}')
    with open('mctdh.op', 'a') as f:
        f.write(f'l{ist}_m{imode} = {linear_diag_ev}, ev\n')
        f.write(f'q{ist}_m{imode} = {quadratic_diag_ev}, ev\n')
    #loop over jst
    jlast = ist - 1

#Above lines 542-560

jlast = 10
filnam = "file"
imode = 1
qsize = 0.1

with open("mctdh.op", "a") as mctdh_op:
    for jst in range(1, jlast+1):
        # grep and tail commands
        cmd_plus = f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_mode{imode}+{qsize}.out | tail -1 | cut -c62-'
        cmd_minus = f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_mode{imode}-{qsize}.out | tail -1 | cut -c62-'
        cmd_plusx2 = f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_mode{imode}+{qsize}x2.out | tail -1 | cut -c62-'
        cmd_minusx2 = f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_mode{imode}-{qsize}x2.out | tail -1 | cut -c62-'
        cmd_refG = f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_refG.out | tail -1 | cut -c62-'

        # execute grep and cut commands using subprocess
        Coup_ev_plus = subprocess.check_output(cmd_plus, shell=True, text=True).strip()
        Coup_ev_minus = subprocess.check_output(cmd_minus, shell=True, text=True).strip()
        Coup_ev_plusx2 = subprocess.check_output(cmd_plusx2, shell=True, text=True).strip()
        Coup_ev_minusx2 = subprocess.check_output(cmd_minusx2, shell=True, text=True).strip()
        Coup_ev_refG = subprocess.check_output(cmd_refG, shell=True, text=True).strip()

        # evaluate linear and quadratic off-diagonal elements
        linear_offdiag_ev = (float(Coup_ev_plus) - float(Coup_ev_minus)) / (2 * qsize)
        quadratic_offdiag_ev = (float(Coup_ev_plusx2) + float(Coup_ev_minusx2) - 2.0 * float(Coup_ev_refG)) / (4.0 * qsize * qsize)

        # write results to mctdh.op
        mctdh_op.write(f"{jst} {ist} {linear_offdiag_ev}\n")
        mctdh_op.write(f"l{jst}{ist}_m{imode} = {linear_offdiag_ev}, ev\n")
        mctdh_op.write(f"q{jst}{ist}_m{imode} = {quadratic_offdiag_ev}, ev\n")

    mctdh_op.write("\n")
#Above lines 562-572


#Should be an else statement for echo "not good to extract. Skipping mode" $imode "for extracting vibronic couplings"


#Above lines 577-579

# Extracting bilinear vibronic coupling
with open("mctdh.op", "a") as outfile:
    outfile.write("#Bilinear diagonal and off-diagonal vibronic coupling constants:\n")
    
    # Determine last mode
    lmode_last = kmode - 1
    
    # Iterate through modes
    for lmode in range(1, lmode_last + 1):
        jmode = modes_included[lmode]
        
        # Check whether the four displacements are calculated properly
        grace_code_pp = subprocess.call(["grep", "grace", "{}_mode{}_+{}_mode{}_+{}.out".format(filnam, imode, qsize, jmode, qsize)])
        grace_code_pm = subprocess.call(["grep", "grace", "{}_mode{}_+{}_mode{}_-{}.out".format(filnam, imode, qsize, jmode, qsize)])
        grace_code_mp = subprocess.call(["grep", "grace", "{}_mode{}_-{}_mode{}_+{}.out".format(filnam, imode, qsize, jmode, qsize)])
        grace_code_mm = subprocess.call(["grep", "grace", "{}_mode{}_-{}_mode{}_-{}.out".format(filnam, imode, qsize, jmode, qsize)])


#Above lines 581-596

if grace_code_pp == 0 and grace_code_pm == 0 and grace_code_mp == 0 and grace_code_mm == 0:
    print(f"Good to extract bilinear for modes {imode} {jmode}")
    
    # Iterate through states
    for ist in range(1, nstate + 1):
        # Extract diabatic energy values from output files
        Ediab_au_pp = subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out']).split(b'\n')[-2].decode()[43:60]
        Ediab_au_pm = subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out']).split(b'\n')[-2].decode()[43:60]
        Ediab_au_mp = subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out']).split(b'\n')[-2].decode()[43:60]
        Ediab_au_mm = subprocess.check_output(['grep', f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=', f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out']).split(b'\n')[-2].decode()[43:60]
        
        # Calculate bilinear diagonal coupling constant
        bilinear_diag_ev = (float(Ediab_au_pp) + float(Ediab_au_mm) - float(Ediab_au_pm) - float(Ediab_au_mp)) * ha2ev / (4.0 * qsize * qsize)
        
        print(f"{ist} {bilinear_diag_ev}")
        
        # Write to output file
        with open("mctdh.op", "a") as outfile:
            outfile.write(f"b{ist}_m{imode}_m{jmode} = {bilinear_diag_ev}, ev\n")

#Above lines 598-610

# loop over jst
jlast = int(ist) - 1
for jst in range(1, jlast+1):
    # extract coupling values
    Coup_ev_pp = subprocess.check_output(f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_mode{imode}+{qsize}_mode{jmode}+{qsize}.out|tail -1|cut -c62-', shell=True, text=True).strip()
    Coup_ev_pm = subprocess.check_output(f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_mode{imode}+{qsize}_mode{jmode}-{qsize}.out|tail -1|cut -c62-', shell=True, text=True).strip()
    Coup_ev_mp = subprocess.check_output(f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_mode{imode}-{qsize}_mode{jmode}+{qsize}.out|tail -1|cut -c62-', shell=True, text=True).strip()
    Coup_ev_mm = subprocess.check_output(f'grep "STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING" {filnam}_mode{imode}-{qsize}_mode{jmode}-{qsize}.out|tail -1|cut -c62-', shell=True, text=True).strip()

    # calculate bilinear_offdiag_ev
    qsize = float(qsize)
    bilinear_offdiag_ev = (float(Coup_ev_pp) + float(Coup_ev_mm) - float(Coup_ev_pm) - float(Coup_ev_mp)) / (4.0 * qsize * qsize)

    print(jst, ist, bilinear_offdiag_ev)
    with open('mctdh.op', 'a') as f:
        f.write(f'b{jst}{ist}_m{imode}_m{jmode} = {bilinear_offdiag_ev}, ev\n')

#Above lines 612-624


with open('mctdh.op', 'a') as f:
    f.write('\n')
#Line 625

#echo "not good to extract. Skipping mode $imode mode $jmode for extracting bilinear vibronic couplings"

#Line 628


# Now we construct the Hamiltonian section in mctdh.op
with open('mctdh.op', 'a') as f:
    f.write('\nend-parameter-section\n')
    f.write('-----------------------------------------\n')
    f.write('HAMILTONIAN-SECTION\n')
    f.write('-----------------------------------------\n')

    f.write(' modes | el ')
    for imode_include in range(1, nmodes_include+1):
        f.write(f'| m{modes_included[imode_include]} ')
    f.write('\n')
    f.write('-----------------------------------------\n')
    f.write('# KINETIC OPERATOR FOR NORMAL MODES\n')
    f.write('-----------------------------------------\n')
    for imode_include in range(1, nmodes_include+1):
        mode_count = imode_include + 1
        f.write(f'w_{modes_included[imode_include]}   |{mode_count} KE\n')
    f.write('-----------------------------------------\n')
    f.write('# HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES\n')
    f.write('-----------------------------------------\n')
    for imode_include in range(1, nmodes_include+1):
        mode_count = imode_include + 1
        f.write(f'0.5*w_{modes_included[imode_include]}   |{mode_count}  q^2\n')

#Above lines 637-668

with open('mctdh.op', 'a') as f:
    f.write('-----------------------------------------\n')
    f.write('# ELECTRONIC COUPLING AT REFERENCE STRUCTURE\n')
    f.write('-----------------------------------------\n')
    for ist in range(1, nstate+1):
        f.write(f'v{ist}  |1 S{ist}&{ist}\n')
    for ist in range(1, nstate+1):
        jlast = ist - 1
        for jst in range(1, jlast+1):
            f.write(f'v{jst}{ist}  |1 S{jst}&{ist}\n')
    f.write('-----------------------------------------\n')
    f.write('# LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS\n')
    f.write('-----------------------------------------\n')

#Above lines 669-686

with open('mctdh.op', 'a') as f:
    for kmode in range(1, nmodes_include+1):
        imode = modes_included[kmode-1]
        kmode_count = kmode + 1
        for ist in range(1, nstate+1):
            f.write(f'l{ist}_m{imode} |1 S{ist}&{ist} |{kmode_count} q\n')
            f.write(f'q{ist}_m{imode} |1 S{ist}&{ist} |{kmode_count} q^2\n')
    f.write('-----------------------------------------\n')
    f.write('# LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS\n')
    f.write('-----------------------------------------\n')
#Above lines 687-699

for kmode in range(1, nmodes_include + 1):
    imode = modes_included[kmode - 1]
    kmode_count = kmode + 1
    for ist in range(1, nstate + 1):
        jlast = ist - 1
        for jst in range(1, jlast + 1):
            l = "l{}{}_m{} |1".format(jst, ist, imode)
            s = "S{}&{} |{} q".format(jst, ist, kmode_count)
            print(l, s, file=open("mctdh.op", "a"))
            q = "q{}{}_m{} |1".format(jst, ist, imode)
            s = "S{}&{} |{} q^2".format(jst, ist, kmode_count)
            print(q, s, file=open("mctdh.op", "a"))

print("-----------------------------------------", file=open("mctdh.op", "a"))
print("# BILINEAR DIAGONAL VIBRONIC COUPLINGS", file=open("mctdh.op", "a"))
print("-----------------------------------------", file=open("mctdh.op", "a"))


#Above lines 700-716

for kmode in range(1, nmodes_include + 1):
    imode = modes_included[kmode - 1]
    kmode_count = kmode + 1
    lmode_last = kmode - 1
    for lmode in range(1, lmode_last + 1):
        jmode = modes_included[lmode - 1]
        lmode_count = lmode + 1
        for ist in range(1, nstate + 1):
            b = "b{}_m{}_m{} |1".format(ist, imode, jmode)
            s = "S{}&{} |{} q |{} q".format(ist, ist, lmode_count, kmode_count)
            print(b, s, file=open("mctdh.op", "a"))

print("-----------------------------------------", file=open("mctdh.op", "a"))
print("# BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS", file=open("mctdh.op", "a"))
print("-----------------------------------------", file=open("mctdh.op", "a"))


#Above lines 717-734

for kmode in range(1, nmodes_include + 1):
    imode = modes_included[kmode - 1]
    kmode_count = kmode + 1
    lmode_last = kmode - 1
    for lmode in range(1, lmode_last + 1):
        jmode = modes_included[lmode - 1]
        lmode_count = lmode + 1
        for ist in range(1, nstate + 1):
            jlast = ist - 1
            for jst in range(1, jlast + 1):
                b = "b{}_m{}_m{} |1".format(jst*10+ist, imode, jmode)
                s = "S{}&{} |{} q |{} q".format(jst, ist, lmode_count, kmode_count)
                print(b, s, file=open("mctdh.op", "a"))

print("-----------------------------------------", file=open("mctdh.op", "a"))
print("", file=open("mctdh.op", "a"))
print("end-hamiltonian-section", file=open("mctdh.op", "a"))
print("", file=open("mctdh.op", "a"))
print("end-operator", file=open("mctdh.op", "a"))

#Above lines 735-758