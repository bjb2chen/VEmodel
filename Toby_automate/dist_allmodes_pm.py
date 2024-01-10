import sys
import pprint
import subprocess
import os
import shutil
import re
import json

# Function to get the number of atoms from thFe hessout file
def get_number_of_atoms(hessout):
    with open(hessout, 'r', errors='replace') as hess_file:
        for line in hess_file:
            if ' TOTAL NUMBER OF ATOMS' in line:
                natoms = int(line.split('=')[1])
                return natoms

# Function to extract lines between patterns in a file
def extract_lines_between_patterns(filename, start_pattern, end_pattern):
    selected_lines = []
    collecting = False

    with open(filename, 'r', errors='replace') as file:
        for line in file:
            if start_pattern in line:
                collecting = True
                selected_lines.append(line)
            elif end_pattern in line:
                collecting = False
            elif collecting:
                selected_lines.append(line)

    return selected_lines

# Function to read frequency values from selected lines
def read_freq_values(selected_lines):
    freq_value_set = []

    for freqline in selected_lines:
        if "FREQUENCY:" in freqline:
            freq_value_set.append(freqline[18:])

    return freq_value_set

# Function to extract filtered set of lines
def read_mode_values(selected_lines):
    mode_value_set = []

    for idx, modeline in enumerate(selected_lines):
        if len(modeline) > 3 and modeline[2].isdigit():
            # mode_value_set.append(selected_lines[idx][20:])
            # mode_value_set.append(selected_lines[idx+1][20:])
            # mode_value_set.append(selected_lines[idx+2][20:])
            for i in range(3): # 3 for x,y,z
                mode_value_set.append(selected_lines[idx + i][20:])

    return mode_value_set

def process_mode_freq(natoms, ndim, ngroup, nleft):
    nrmmod = {}  # normal modes
    freqcm = {}  # frequencies in cm

    with open("mode.dat", 'r', errors='replace') as mode_file:
        lines_mode = mode_file.readlines()

    with open("freq.dat", 'r', errors='replace') as freq_file:
        lines_freq = freq_file.readlines()

    for igroup in range(1, ngroup + 1):
        iniline = (igroup - 1) * ndim + 1
        endline = iniline + ndim - 1
        print("igroup =", igroup)
        print("iniline =", iniline)
        #print("endline =", endline)
        ixyz = 0

        for line in range(iniline, endline + 1, 1):
            ixyz += 1
            print(ixyz)
            #print('my group line here', line)

            for icolumn in range(1, 6, 1):
                imode = (igroup - 1) * 5 + icolumn
                print(ixyz, imode, "")
                cutini = (icolumn - 1) * 12
                cutfnl = icolumn * 12

                disp = lines_mode[line - 1][cutini:cutfnl].lstrip()
                nrmmod[ixyz, imode] = float(disp) # no need to suppress sci notation, GMS is OK with E^...
                print(nrmmod[ixyz, imode], disp)

                if ixyz == 1:
                    cutini = (icolumn - 1) * 12
                    cutfnl = icolumn * 12
                    # sometimes imaginary freq ('I ... freqval') so +2
                    freq = lines_freq[igroup - 1][cutini+2:cutfnl].lstrip()
                    freqcm[imode] = float(freq)
                    print("frequency:", imode, freqcm[imode])

    # For the leftover nleft modes
    if nleft != 0:
        for igroup in range(nleft, nleft + 1):
            iniline = ngroup * ndim + 1 
            endline = iniline + ndim - 1
            print("igroup=leftover")
            print("iniline =", iniline)
            ixyz = 0

            for line in range(iniline, endline + 1):
                ixyz += 1
                print(ixyz)
                #print('my leftover group line here', line)

                for icolumn in range(1, nleft + 1):
                    imode = ngroup * 5 + icolumn
                    print(ixyz, imode, "")
                    cutini = (icolumn - 1) * 12
                    cutfnl = icolumn * 12

                    disp = lines_mode[line - 1][cutini:cutfnl].lstrip()
                    nrmmod[ixyz, imode] = float(disp)
                    print(nrmmod[ixyz, imode], disp)

                    if ixyz == 1:
                        cutini = (icolumn - 1) * 12
                        cutfnl = icolumn * 12
                        freq = lines_freq[-1][cutini+2:cutfnl].lstrip()
                        freqcm[imode] = float(freq)
                        print("frequency:", imode, freqcm[imode])

    # # Check if imode is in modes_excluded and exclude if necessary
    # for imode in modes_excluded:
    #     if imode in freqcm:
    #         del freqcm[imode]

    #Print all frequencies
    for imode in range(1, ndim + 1):
        print("frequency:", imode, str(freqcm[imode]), "CM-1")

    return nrmmod, freqcm

def compose_ref_structure(hessout, natoms):
    coord_lines = extract_lines_between_patterns(hessout, 'EQUILIBRIUM GEOMETRY LOCATED', 'INTERNUCLEAR DISTANCES')
    ref_file = 'py_ref_structure'

    if len(coord_lines) > 2:
        try:
            subprocess.run(['rm', '-f', ref_file])
        except Exception as e:
            print(f"Error deleting {ref_file}: {str(e)}")

        for atom_line in coord_lines[-natoms-1:-1]:
            with open(f'{ref_file}', "a") as python_ref_file:
                python_ref_file.write(atom_line)

        print(f'Sucessfully extracted equilibrium geometry from {hessout} and prepared ref_structure.')
    else:
        print(f'Unsucessful extraction of equilibrium geometry from {hessout}. Please prepare ref_structure manually.')

    return coord_lines

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

    with open(file_path, 'r', errors='replace') as struct_file:
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

def filter_modes(excluded_set, ndim):
    modes_included = {}
    print("NUMBER OF EXCLUDED MODES:", len(excluded_set))
    print("They are modes:", *excluded_set)

    # Iterate through modes to check inclusion
    for imode in range(1, ndim + 1):
        if imode not in excluded_set:
            modes_included[len(modes_included)+1] = imode
            print(len(modes_included), imode)

    print("Number of Modes Included:", len(modes_included))        

    return modes_included

def my_subgam(filnam, **kwargs):
    ncpus = kwargs.get('ncpus', 2)
    ngb = kwargs.get('ngb', 2)
    nhour = kwargs.get('nhour', 1)
    # Remove the ".inp" extension from the filename
    input_no_ext, extension = os.path.splitext(filnam)
    print(f"running calculations for {input_no_ext}")
    wd = os.getcwd()

    with open(f"{input_no_ext}.slurm", "w") as slurm_file:
        slurm_file.write("#!/bin/bash\n")
        slurm_file.write("#SBATCH --nodes=1\n")
        slurm_file.write(f"#SBATCH --ntasks={ncpus}\n")
        slurm_file.write(f"#SBATCH --mem-per-cpu={ngb}G\n")
        slurm_file.write(f"#SBATCH --time={nhour}:00:00\n")
        slurm_file.write("\n")
        slurm_file.write("cd $SLURM_SUBMIT_DIR\n")
        slurm_file.write("\n")
        slurm_file.write("export SLURM_CPUS_PER_TASK\n")
        slurm_file.write('mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID\n')
        slurm_file.write("\n")
        slurm_file.write(f"/home/$USER/LOCAL/runG_diab {input_no_ext}.inp {ncpus} \n")

    command = (
        "sbatch"
        f" {input_no_ext}.slurm"
    )

    return f"{input_no_ext}.slurm"

#Do diabatization calculation at the reference nondistorted structure.
#This calculation shall be a repetition of a calcualtion in preparing temp.inp
def refG_calc(refgeo, filnam):
    # Check if the calculation has already been run
    grace_exists = subprocess.call(["grep", "grace", f"{filnam}_refG.out"]) == 0
    #grace_exists = subprocess.run(["grep", "grace", f"{filnam}_refG.out"], shell=True).returncode == 0
    if not grace_exists:
        print("Run calculation at the undistorted reference structure")

        shutil.copy("temp.inp", f"{filnam}_refG.inp")

        with open(f"{filnam}_refG.inp", "a") as inp_file:
            with open(refgeo, 'r', errors='replace') as ref_structure:
                inp_file.write(ref_structure.read())

        with open(f"{filnam}_refG.inp", "a") as inp_file:
            inp_file.write(" $END\n")

        # Submit and run the refG calculation (you may need to customize this command based on your setup)
        #refG_job_result = subprocess.run(["./subgam.diab", f"{filnam}_refG.inp", "4", "0", "1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os.system("sbatch" + ' -W' + " " + my_subgam(f"{filnam}_refG.inp", ncpus=2, ngb=1, nhour=1)) # the wait 

        # At this point, refG calculation has completed successfully. 
        print("Calculation at the reference structure is done.")
    else:
        print("Calculation at the reference structure has already been done.")

    return

def diabatization(filnam, modes_included, **kwargs):

    distcoord_plus = {}
    distcoord_minus = {}
    distcoord_plus_x2 = {}
    distcoord_minus_x2 = {}
    distcoord_pp = {}
    distcoord_pm = {}
    distcoord_mp = {}
    distcoord_mm = {}

    freqcm = kwargs.get('freqcm')
    ndim = kwargs.get('ndim')
    refcoord = kwargs.get('refcoord')
    nrmmod = kwargs.get('nrmmod')
    natoms = kwargs.get('natoms')
    atmlst = kwargs.get('atmlst')
    chrglst = kwargs.get('chrglst')
    qsize = kwargs.get('qsize', 0.05)
    ha2ev = kwargs.get('ha2ev', 27.2113961318)
    wn2ev = kwargs.get('wn2ev', 0.000123981)
    wn2eh = kwargs.get('wn2eh', 0.00000455633)
    ang2br = kwargs.get('ang2br', 1.889725989)
    amu2me = kwargs.get('amu2me', 1822.88839)

    #Loop over all considered modes and do + and - displacements

    for kmode in range(1, len(modes_included) + 1):
        imode = modes_included[kmode]

        # Convert the reduced dimensionless qsize to the actual rsize in sqrt(amu)*Angs unit
        omga = freqcm[imode]
        rsize = qsize / (pow(amu2me, 0.5) * ang2br * pow(omga * wn2eh, 0.5))
        print(imode, omga, rsize, type(rsize))

        # Loop over components
        for icomp in range(1, ndim + 1):
            coord_disp_plus = refcoord[icomp] + rsize * nrmmod[icomp, imode]
            coord_disp_minus = refcoord[icomp] - rsize * nrmmod[icomp, imode]
            distcoord_plus[icomp] = coord_disp_plus
            distcoord_minus[icomp] = coord_disp_minus
            coord_disp_plusx2 = refcoord[icomp] + 2.0 * rsize * nrmmod[icomp, imode]
            coord_disp_minusx2 = refcoord[icomp] - 2.0 * rsize * nrmmod[icomp, imode]
            distcoord_plus_x2[icomp] = coord_disp_plusx2
            distcoord_minus_x2[icomp] = coord_disp_minusx2
            print(imode, icomp, refcoord[icomp], nrmmod[icomp, imode], coord_disp_plus, coord_disp_minus,
            distcoord_plus[icomp], distcoord_minus[icomp])

        # Delete existing dist_structure files
        for dist_file in ['dist_structure_plus', 'dist_structure_minus', 'dist_structure_plusx2', 'dist_structure_minusx2']:
            try:
                subprocess.run(['rm', '-f', dist_file])
            except Exception as e:
                print(f"Error deleting {dist_file}: {str(e)}")

        # Print the distorted structure
        for iatom in range(1, natoms + 1):
            with open('dist_structure_plus', 'a') as f_plus, \
                    open('dist_structure_minus', 'a') as f_minus, \
                    open('dist_structure_plusx2', 'a') as f_plusx2, \
                    open('dist_structure_minusx2', 'a') as f_minusx2:
                f_plus.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
                f_minus.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
                f_plusx2.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
                f_minusx2.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
                for ixyz in range(1, 4):
                    icomp = (iatom - 1) * 3 + ixyz
                    f_plus.write(f"{distcoord_plus[icomp]} ")
                    f_minus.write(f"{distcoord_minus[icomp]} ")
                    f_plusx2.write(f"{distcoord_plus_x2[icomp]} ")
                    f_minusx2.write(f"{distcoord_minus_x2[icomp]} ")
                f_plus.write('\n')
                f_minus.write('\n')
                f_plusx2.write('\n')
                f_minusx2.write('\n')
 
        # Create input files for diabatization calculations
        for displacement in ['+', '-']:
            if displacement == '+':
                p_or_m = 'plus'
            elif displacement == '-':
                p_or_m = 'minus'

            for suffix in ['', 'x2']:
                shutil.copy('temp.inp', f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.inp')
                with open(f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.inp', 'a') as inp_file:
                    with open(f'dist_structure_{p_or_m}{suffix}', 'r', errors='replace') as dist_file:
                        inp_file.write(dist_file.read())
                        inp_file.write(' $END')

                # Check if the calculation is done already
                grace1 = subprocess.run(["grep", "grace", f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.out'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                #if not os.path.exists(f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.out'):
                if grace1.returncode != 0:
                    print(f"Running calculations for {filnam}_mode{imode}_{displacement}{qsize}{suffix}")
                    try:
                        #subprocess.run(['./subgam.diab', f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.inp', '4', '0', '1'])
                        os.system("sbatch" + " " + my_subgam(f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.inp', ncpus=2, ngb=1, nhour=1))
                    except Exception as e:
                        print(f"Error running diabatization calculation: {str(e)}")
                else:
                    print(f"{filnam}_mode{imode}_{displacement}{qsize}{suffix} is done")

        # 2D distortion to get bilinear vibronic coupling
        for lmode in range(1, kmode):
            jmode = modes_included[lmode]
 
            # Convert the reduced dimensionless qsize to the actual rsizep in sqrt(amu)*Angs unit for jmode
            omgap = freqcm[jmode]
            rsizep = qsize / (pow(amu2me, 0.5) * ang2br * pow(omgap * wn2eh, 0.5))
            print(imode, jmode, rsize, rsizep)
 
            # Loop over components
            for icomp in range(1, ndim + 1):
                coord_disp_pp = distcoord_plus[icomp] + rsizep * nrmmod[icomp, jmode]
                coord_disp_pm = distcoord_plus[icomp] - rsizep * nrmmod[icomp, jmode]
                coord_disp_mp = distcoord_minus[icomp] + rsizep * nrmmod[icomp, jmode]
                coord_disp_mm = distcoord_minus[icomp] - rsizep * nrmmod[icomp, jmode]
                distcoord_pp[icomp] = coord_disp_pp
                distcoord_pm[icomp] = coord_disp_pm
                distcoord_mp[icomp] = coord_disp_mp
                distcoord_mm[icomp] = coord_disp_mm

            # Delete existing dist_structure files
            for dist_file in ['dist_structure_pp', 'dist_structure_pm', 'dist_structure_mp', 'dist_structure_mm']:
                try:
                    subprocess.run(['rm', '-f', dist_file])
                except Exception as e:
                    print(f"Error deleting {dist_file}: {str(e)}")
 
            # Print the distorted structure
            for iatom in range(1, natoms + 1):
                with open(f'dist_structure_pp', 'a') as f_pp, \
                        open(f'dist_structure_pm', 'a') as f_pm, \
                        open(f'dist_structure_mp', 'a') as f_mp, \
                        open(f'dist_structure_mm', 'a') as f_mm:
                    f_pp.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
                    f_pm.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
                    f_mp.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
                    f_mm.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
                    for ixyz in range(1, 4):
                        icomp = (iatom - 1) * 3 + ixyz
                        f_pp.write(f"{distcoord_pp[icomp]} ")
                        f_pm.write(f"{distcoord_pm[icomp]} ")
                        f_mp.write(f"{distcoord_mp[icomp]} ")
                        f_mm.write(f"{distcoord_mm[icomp]} ")
                    f_pp.write('\n')
                    f_pm.write('\n')
                    f_mp.write('\n')
                    f_mm.write('\n')
 
            # Create input files for diabatization calculations
            for displacement1 in ['+', '-']:
                for displacement2 in ['+', '-']:
                    if displacement1 == '+' and displacement2 == '+':
                        suffix1 = 'pp'
                    elif displacement1 == '+' and displacement2 == '-':
                        suffix1 = 'pm'
                    elif displacement1 == '-' and displacement2 == '+':
                        suffix1 = 'mp'
                    elif displacement1 == '-' and displacement2 == '-':
                        suffix1 = 'mm'

                    shutil.copy('temp.inp', f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp')
                    with open(f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp', 'a') as inp_file:
                        with open(f'dist_structure_{suffix1}', 'r', errors='replace') as dist_file:
                            inp_file.write(dist_file.read())
                        inp_file.write(' $END ')
 
                    # Check if the calculation is done already
                    output_filename = f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.out'
                    grace2 = subprocess.run(["grep", "grace", output_filename])
                    if grace2.returncode != 0:
                    #if not os.path.exists(output_filename):
                        print(f"Running calculations for {output_filename}!")
                        try:
                            #subprocess.run(['./subgam.diab', f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp', '4', '0', '1'])
                            os.system("sbatch" + " " + my_subgam(f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp', ncpus=2, ngb=1, nhour=1))
                        except Exception as e:
                            print(f"Error running diabatization calculation: {str(e)}")
                    else:
                        print(f"{output_filename} is already done.")

#Now we move on to extract vibronic coupling constants using finite difference
#and write the data in an mctdh operator file

def extract_diabatic_energy(file_path, pattern):
    with open(file_path, 'r', errors='replace') as file:
        for line in reversed(file.readlines()):
            match = re.search(pattern, line)
            if match:
                return float(line[44:62].strip().replace(" ", ""))

def extract_coupling_energy(file_path, pattern):
    with open(file_path, 'r', errors='replace') as file:
        for line in reversed(file.readlines()):
            match = re.search(pattern, line)
            if match:
                return float(line[62:].strip().replace(" ", ""))

def find_nstate(file_path, pattern='# of states in CI      = ', encoding="utf-8"):
    with open(file_path, 'r', encoding=encoding, errors='replace') as file:
        for line in file:
            if pattern in line:
                return int(line.split('=')[1].strip())
    return None  # Return None if the pattern is not found

def extract_same_state_transition_dipoles(selected_lines):
    same_state_transition_dipoles = {}

    for TDIPOLEline in selected_lines:
        try:
            if TDIPOLEline[0:5].strip().isnumeric():
                state1 = int(TDIPOLEline[0:5].strip())
                #print(state1)
                state2 = int(TDIPOLEline[5:10].strip())
    
                if state1 == state2:
                    x = float(TDIPOLEline[11:21].strip())
                    y = float(TDIPOLEline[22:31].strip())
                    z = float(TDIPOLEline[32:42].strip())
                    same_state_transition_dipoles[state1] = (f"{x:.6f}", f"{y:.6f}", f"{z:6f}")
        except Exception as e:
            print(f"ERror processing line: {TDIPOLEline} - {e}")

    return same_state_transition_dipoles

def mctdh(filnam, modes_included, **kwargs):
    nmodes = len(modes_included)
    freqcm = kwargs.get('freqcm')
    qsize = kwargs.get('qsize', 0.05)
    ha2ev = kwargs.get('ha2ev', 27.2113961318)
    wn2ev = kwargs.get('wn2ev', 0.000123981)
    wn2eh = kwargs.get('wn2eh', 0.00000455633)
    ang2br = kwargs.get('ang2br', 1.889725989)
    amu2me = kwargs.get('amu2me', 1822.88839)

    try:
        subprocess.run(['rm', '-f', 'mctdh.op'])
    except Exception as e:
        print(f"Error deleting {'mctdh.op'}: {str(e)}")

    #Heading for mctdh.op
    str1 = "OP_DEFINE-SECTION"
    str2 = "title"

    nstate = find_nstate(f'{filnam}_refG.out')

    str3 = "end-title "
    str4 = "end-op_define-section"
    str5 = ""
    # lines 482,483
    str6 = "PARAMETER-SECTION"
    str7 = ""
    str8 = f'{filnam} {nstate} states + ' + str(nmodes) + ' modes'
    strlst = [str1, str2, str8, str3, str4, str5, str6, str7]

    with open('mctdh.op', 'w') as mctdh_file:
        for idx in strlst:
            mctdh_file.write(idx+'\n')

    refG_exists = subprocess.call(["ls", f"{filnam}_refG.out"]) == 0
    if refG_exists:
        with open("mctdh.op", "a") as mctdh_file:
            mctdh_file.write("#Diagonal and Off-diagonal diabatic Hamiltonian elements at reference structure\n")
    
            for ist in range(1, nstate + 1):
                with open(f"{filnam}_refG.out", 'r', errors='replace') as refG_out:
                    lines = refG_out.readlines()
    
                    # Extract diabatic energy for state ist
                    Ediab = None
                    for line in reversed(lines):
                        state_pattern = re.compile(f"STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=")
                        #if ("STATE #" in line) and ("S GMC-PT-LEVEL DIABATIC ENERGY=" in line):
                        match = state_pattern.search(line)
                        if match:
                            Ediab = line[61:].strip().replace(" ", "")
                            break
    
                    mctdh_file.write(f"v{ist} = {Ediab}, ev\n")
    
                    # Extract coupling energy between state jst and ist
                    for jst in range(1, ist):
                        Coup_ev = None
                        for line in reversed(lines):
                            state_pattern = re.compile(f"STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING")
                            #if ("STATE #" in line) and ("S GMC-PT-LEVEL COUPLING" in line):
                            match = state_pattern.search(line)
                            if match:
                                Coup_ev = line[61:].strip().replace(" ", "")
                                break
    
                        mctdh_file.write(f"v{jst}{ist} = {Coup_ev}, ev\n")
    
                mctdh_file.write("\n")
        
            mctdh_file.write("\n")
    else:
        print(f"Skip extracting Hamiltonians from the non-existing {filnam}_refG.out")

    # Loop through modes
    for kmode in range(1, nmodes + 1):
        imode = modes_included[kmode]
    
        vibron_ev = freqcm[imode] * wn2ev
        with open("mctdh.op", "a") as mctdh_file:
            mctdh_file.write(f"#Parameters for mode {imode}\n")
            mctdh_file.write("#Vibron:\n")
            mctdh_file.write(f"w_m{imode} = {vibron_ev:.16f}, ev\n\n")
            mctdh_file.write("#Linear and quadratic diagonal and off-diagonal vibronic coupling constants:\n")
    
            grace_code_plus = subprocess.call(["grep", "grace", f"{filnam}_mode{imode}_+{qsize}.out"])
            grace_code_minus = subprocess.call(["grep", "grace", f"{filnam}_mode{imode}_-{qsize}.out"])
            grace_code_plusx2 = subprocess.call(["grep", "grace", f"{filnam}_mode{imode}_+{qsize}x2.out"])
            grace_code_minusx2 = subprocess.call(["grep", "grace", f"{filnam}_mode{imode}_-{qsize}x2.out"])
    
            if all(code == 0 for code in [grace_code_plus, grace_code_minus, grace_code_plusx2, grace_code_minusx2]):
                print("\n good to extract\n")
                # Extract the diagonal and off-diagonal vibronic coupling
                for ist in range(1, nstate + 1):
                    pattern = f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY='
                    # Extract Ediab_au_plus
                    Ediab_au_plus = extract_diabatic_energy(f'{filnam}_mode{imode}_+{qsize}.out', pattern)
                    # Extract Ediab_au_plusx2
                    Ediab_au_plusx2 = extract_diabatic_energy(f'{filnam}_mode{imode}_+{qsize}x2.out', pattern)
                    # Extract Ediab_au_minus
                    Ediab_au_minus = extract_diabatic_energy(f'{filnam}_mode{imode}_-{qsize}.out', pattern)
                    # Extract Ediab_au_minusx2
                    Ediab_au_minusx2 = extract_diabatic_energy(f'{filnam}_mode{imode}_-{qsize}x2.out', pattern)
                    # Extract Ediab_au_0
                    Ediab_au_0 = extract_diabatic_energy(f'{filnam}_refG.out', pattern)

                    # print("Ediab_au_plus: ", Ediab_au_plus)
                    # print("Ediab_au_plusx2: ", Ediab_au_plusx2)
                    # print("Ediab_au_minus: ", Ediab_au_minus)
                    # print("Ediab_au_minusx2: ", Ediab_au_minusx2)
                    # print("Ediab_au_0: ", Ediab_au_0)
    
                    linear_diag_ev = (Ediab_au_plus - Ediab_au_minus) * ha2ev / (2 * qsize)
                    quadratic_diag_ev = (Ediab_au_plusx2 + Ediab_au_minusx2 - 2.0 * Ediab_au_0) * ha2ev / (4.0 * qsize * qsize)
    
                    # We only view the difference between the actual force constant and the vibron
                    # as the quadratic diagonal coupling for the diabatic state.
                    quadratic_diag_ev = quadratic_diag_ev - vibron_ev
    
                    # Print and store results
                    print(f"{ist} {linear_diag_ev} {quadratic_diag_ev}, ev\n")
                    # machine accuracy is typically 16 digits
                    mctdh_file.write(f"l{ist}_m{imode} = {linear_diag_ev:.16f}, ev\n")
                    mctdh_file.write(f"q{ist}_m{imode} = {quadratic_diag_ev:.16f}, ev\n")
    
                    # # Loop over jst
                    jlast = ist - 1
                    for jst in range(1, jlast + 1):
                        pattern = f'STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING'
                        # Extract Coup_ev_plus
                        Coup_ev_plus = extract_coupling_energy(f'{filnam}_mode{imode}_+{qsize}.out', pattern)
                        # Extract Coup_ev_plusx2
                        Coup_ev_plusx2 = extract_coupling_energy(f'{filnam}_mode{imode}_+{qsize}x2.out', pattern)
                        # Extract Coup_ev_minus
                        Coup_ev_minus = extract_coupling_energy(f'{filnam}_mode{imode}_-{qsize}.out', pattern)
                        # Extract Coup_ev_minusx2
                        Coup_ev_minusx2 = extract_coupling_energy(f'{filnam}_mode{imode}_-{qsize}x2.out', pattern)
                        # Extract Coup_ev_0
                        Coup_ev_0= extract_coupling_energy(f'{filnam}_refG.out', pattern)

                        # print("Coup_ev_plus: ", Coup_ev_plus)
                        # print("Coup_ev_plusx2: ", Coup_ev_plusx2)
                        # print("Coup_ev_minus: ", Coup_ev_minus)
                        # print("Coup_ev_minusx2: ", Coup_ev_minusx2)
                        # print("Coup_ev_0: ", Coup_ev_0)
    
                        # Compute linear off-diagonal coupling
                        linear_offdiag_ev = (Coup_ev_plus - Coup_ev_minus) / (2 * qsize)
    
                        # Compute quadratic off-diagonal coupling
                        quadratic_offdiag_ev = (Coup_ev_plusx2 + Coup_ev_minusx2 - 2.0 * Coup_ev_0) / (4.0 * qsize * qsize)
    
                        # Print and store results
                        print(f"{jst} {ist} {linear_offdiag_ev}\n")
                        mctdh_file.write(f"l{jst}{ist}_m{imode} = {linear_offdiag_ev:.16f}, ev\n")
                        mctdh_file.write(f"q{jst}{ist}_m{imode} = {quadratic_offdiag_ev:.16f}, ev\n")
                    mctdh_file.write("\n")
    
            else:
                mctdh_file.write(f"not good to extract. Skipping mode {imode} for extracting vibronic couplings\n")
    
            # Extracting bilinear vibronic coupling
            mctdh_file.write("#Bilinear diagonal and off-diagonal vibronic coupling constants:\n")
            lmode_last = kmode - 1
            for lmode in range(1, lmode_last + 1):
                jmode = modes_included[lmode]

                grace_code_pp = subprocess.call(["grep", "grace", f"{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out"])
                grace_code_pm = subprocess.call(["grep", "grace", f"{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out"])
                grace_code_mp = subprocess.call(["grep", "grace", f"{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out"])
                grace_code_mm = subprocess.call(["grep", "grace", f"{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out"])

                if all(code == 0 for code in [grace_code_pp, grace_code_pm, grace_code_mp, grace_code_mm]):
                    print(f"\n Good to extract bilinear for modes {imode} {jmode} \n")

                    for ist in range(1, nstate + 1):
    
                        # Extract Ediab_au_pp
                        with open(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out', 'r', errors='replace') as grep_pp:
                            lines = grep_pp.readlines()
    
                            for line in reversed(lines):
                                state_pattern = re.compile(f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=')
                                match = state_pattern.search(line)
                                if match:
                                    Ediab_au_pp = float(line[44:62].strip().replace(" ", ""))
                                    break
    
                        # Extract Ediab_au_pm
                        with open(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out', 'r', errors='replace') as grep_pm:
                            lines = grep_pm.readlines()
                            
                            for line in reversed(lines):
                                state_pattern = re.compile(f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=')
                                match = state_pattern.search(line)
                                if match:
                                    Ediab_au_pm = float(line[44:62].strip().replace(" ", ""))
                                    break
    
                        # Extract Ediab_au_mp
                        with open(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out', 'r', errors='replace') as grep_mp:
                            lines = grep_mp.readlines()
                            
                            for line in reversed(lines):
                                state_pattern = re.compile(f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=')
                                match = state_pattern.search(line)
                                if match:
                                    Ediab_au_mp = float(line[44:62].strip().replace(" ", ""))
                                    break
    
                        # Extract Ediab_au_mm
                        with open(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out', 'r', errors='replace') as grep_mm:
                            lines = grep_mm.readlines()
                            
                            for line in reversed(lines):
                                state_pattern = re.compile(f'STATE #..* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=')
                                match = state_pattern.search(line)
                                if match:
                                    Ediab_au_mm = float(line[44:62].strip().replace(" ", ""))
                                    break
    
                        bilinear_diag_ev = ( Ediab_au_pp + Ediab_au_mm - Ediab_au_pm - Ediab_au_mp ) * ha2ev / (4.0 * qsize * qsize )
                    
                        print(f"{ist} {bilinear_diag_ev}")
                        mctdh_file.write(f"b{ist}_m{imode}_m{jmode} = {bilinear_diag_ev:.16f}, ev\n")
    
                        # # Loop over jst
                        jlast = ist - 1
                        for jst in range(1, jlast + 1):
        
                            # Extract Coup_ev_pp
                            with open(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out', 'r', errors='replace') as grep_coup_pp:
                                lines = grep_coup_pp.readlines()
        
                                for line in reversed(lines):
                                    state_pattern = re.compile(f'STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING')
                                    match = state_pattern.search(line)
                                    if match:
                                        Coup_ev_pp = float(line[62:].strip().replace(" ", ""))
                                        break
        
                            # Extract Coup_ev_pm
                            with open(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out', 'r', errors='replace') as grep_coup_pm:
                                lines = grep_coup_pm.readlines()
                                
                                for line in reversed(lines):
                                    state_pattern = re.compile(f'STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING')
                                    match = state_pattern.search(line)
                                    if match:
                                        Coup_ev_pm = float(line[62:].strip().replace(" ", ""))
                                        break
        
                            # Extract Coup_ev_mp
                            with open(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out', 'r', errors='replace') as grep_coup_mp:
                                lines = grep_coup_mp.readlines()
                                
                                for line in reversed(lines):
                                    state_pattern = re.compile(f'STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING')
                                    match = state_pattern.search(line)
                                    if match:
                                        Coup_ev_mp = float(line[62:].strip().replace(" ", ""))
                                        break
        
                            # Extract Coup_ev_mm
                            with open(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out', 'r', errors='replace') as grep_coup_mm:
                                lines = grep_coup_mm.readlines()
                                
                                for line in reversed(lines):
                                    state_pattern = re.compile(f'STATE #..* {jst} &..* {ist}.S GMC-PT-LEVEL COUPLING')
                                    match = state_pattern.search(line)
                                    if match:
                                        Coup_ev_mm = float(line[62:].strip().replace(" ", ""))
                                        break
        
                            bilinear_offdiag_ev = ( Coup_ev_pp + Coup_ev_mm - Coup_ev_pm - Coup_ev_mp ) / (4.0 * qsize * qsize )
                            
                            print(f"{jst} {ist} {bilinear_offdiag_ev}")
                            mctdh_file.write(f"b{jst}{ist}_m{imode}_m{jmode} = {bilinear_offdiag_ev:.16f}, ev\n")
                        
                        mctdh_file.write("\n")

                else:
                    print(f"not good to extract. Skipping mode {imode} mode {jmode} for extracting bilinear vibronic couplings")

    # Open mctdh.op file for writing
    with open('mctdh.op', 'a') as mctdh_file:
        mctdh_file.write("end-parameter-section\n")
        # Write the header
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("HAMILTONIAN-SECTION\n")
        mctdh_file.write("-----------------------------------------\n")
    
        # Write modes and mode labels
        mctdh_file.write(" modes | el")
        for imode_include in range(1, nmodes + 1):
            mctdh_file.write(f" | m{modes_included[imode_include]}")
        mctdh_file.write("\n")
        mctdh_file.write("-----------------------------------------\n")
    
        # Write KINETIC OPERATOR FOR NORMAL MODES
        mctdh_file.write("# KINETIC OPERATOR FOR NORMAL MODES\n")
        mctdh_file.write("-----------------------------------------\n")
        for imode_include in range(1, nmodes + 1):
            mode_count = imode_include + 1
            mctdh_file.write(f"w_m{modes_included[imode_include]}   |{mode_count} KE\n")
    
        # Write HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES\n")
        mctdh_file.write("-----------------------------------------\n")
        for imode_include in range(1, nmodes + 1):
            mode_count = imode_include + 1
            mctdh_file.write(f"0.5*w_m{modes_included[imode_include]}   |{mode_count}  q^2\n")
    
        # Write ELECTRONIC COUPLING AT REFERENCE STRUCTURE
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# ELECTRONIC COUPLING AT REFERENCE STRUCTURE\n")
        mctdh_file.write("-----------------------------------------\n")
        for ist in range(1, nstate + 1):
            mctdh_file.write(f"v{ist}  |1 S{ist}&{ist}\n")
        for ist in range(1, nstate + 1):
            jlast = ist - 1
            for jst in range(1, jlast + 1):
                mctdh_file.write(f"v{jst}{ist}  |1 S{jst}&{ist}\n")
    
        # Write LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            for ist in range(1, nstate + 1):
                mctdh_file.write(f"l{ist}_m{imode} |1 S{ist}&{ist} |{kmode_count} q\n")
                mctdh_file.write(f"q{ist}_m{imode} |1 S{ist}&{ist} |{kmode_count} q^2\n")
    
        # Write LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            for ist in range(1, nstate + 1):
                jlast = ist - 1
                for jst in range(1, jlast + 1):
                    mctdh_file.write(f"l{jst}{ist}_m{imode} |1 S{jst}&{ist} |{kmode_count} q\n")
                    mctdh_file.write(f"q{jst}{ist}_m{imode} |1 S{jst}&{ist} |{kmode_count} q^2\n")
    
        # Write BILINEAR DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# BILINEAR DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            lmode_last = kmode - 1
            for lmode in range(1, lmode_last + 1):
                jmode = modes_included[lmode]
                lmode_count = lmode + 1
                for ist in range(1, nstate + 1):
                    mctdh_file.write(f"b{ist}_m{imode}_m{jmode} |1 S{ist}&{ist} |{lmode_count} q |{kmode_count} q\n")
    
        # Write BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            lmode_last = kmode - 1
            for lmode in range(1, lmode_last + 1):
                jmode = modes_included[lmode]
                lmode_count = lmode + 1
                for ist in range(1, nstate + 1):
                    jlast = ist - 1
                    for jst in range(1, jlast + 1):
                        mctdh_file.write(f"b{jst}{ist}_m{imode}_m{jmode} |1 S{jst}&{ist} |{lmode_count} q |{kmode_count} q\n")
    
        # Close the file
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("\nend-hamiltonian-section\n\nend-operator\n")
    
    return

def main():
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <hessout_file>")
        sys.exit(1)

    hessout = sys.argv[1]
    filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"
    refgeo = "ref_structure"
    py_refgeo = "py_ref_structure"

    natoms = get_number_of_atoms(hessout)
    ndim = natoms * 3
    #ndim = 15 # note: can only test like this if you have proper hessout structure
    ngroup = ndim // 5
    nleft = ndim % 5
    # qsize = input("Enter desired qsize, default is 0.05: ")
    qsize = 0.05
    print("Dimension of all xyz coordinates:", ndim)
    print("# of atoms:", natoms)
    print(ngroup, nleft)

    #modes_excluded = [1, 2, 3, 4, 5, 6]
    modes_excluded = [1, 4, 5, 6, 7, 8, 9, 10, 11]
    #modes_excluded = [1, 2, 4, 7, 12]

    selected_lines = extract_lines_between_patterns(hessout, 
        "FREQUENCIES IN CM",
        "REFERENCE ON SAYVETZ CONDITIONS"
        )
    freq_value_set = read_freq_values(selected_lines)
    filtered_set = read_mode_values(selected_lines)
    
    with open('mode.dat', 'w') as output_file:
        output_file.writelines(filtered_set)

    with open('freq.dat', 'w') as output_file:
        output_file.writelines(freq_value_set)

    nrmmod, freqcm = process_mode_freq(natoms, ndim, ngroup, nleft)
    ref_file = compose_ref_structure(hessout, natoms)
    atmlst, chrglst, refcoord = read_reference_structure(py_refgeo)
    modes_included = filter_modes(modes_excluded, ndim)

    #Set conversion constants
    ha2ev = 27.2113961318
    wn2ev = 0.000123981
    wn2eh = 0.00000455633
    ang2br = 1.889725989
    amu2me = 1822.88839 

    repetition = refG_calc(py_refgeo, filnam)
    diabatize = diabatization(filnam, modes_included, freqcm=freqcm, ndim=ndim, refcoord=refcoord,\
                           nrmmod=nrmmod, natoms=natoms, atmlst=atmlst, chrglst=chrglst, \
                           qsize=qsize, ha2ev=ha2ev, wn2ev=wn2ev, wn2eh=wn2eh, ang2br=ang2br, amu2me=amu2me)
    make_mctdh = mctdh(filnam, modes_included, freqcm=freqcm, qsize=qsize, ha2ev=ha2ev, wn2ev=wn2ev, wn2eh=wn2eh, ang2br=ang2br, amu2me=amu2me)
    tdipole_block = extract_lines_between_patterns(f"{filnam}_refG.out",
    "TRANSITION DIPOLES BETWEEN DIABATS",
    "TRANSITION DIPOLES BETWEEN DIRECT MAX. DIABATS"
    )
    
    same_state_dipoles = extract_same_state_transition_dipoles(tdipole_block)
    pprint.pprint(tdipole_block)
    pprint.pprint(same_state_dipoles)

    # pprint.pprint(nrmmod)
    # print('---------nrm mod done-----------')
    # pprint.pprint(freqcm)
    # print('---------freqcm done-----------')
    # pprint.pprint(selected_lines)
    # print('---------selected_lines done-----------')
    # pprint.pprint(filtered_set)
    # print('---------filtered_set done-----------')
    # pprint.pprint(freq_value_set)
    # print('---------freq_value_set done-----------')
    # pprint.pprint(atmlst)
    # print('---------atmlst done-----------')
    # pprint.pprint(chrglst)
    # print('---------chrglst done-----------')
    # pprint.pprint(refcoord)
    # print('---------refcoord done-----------')
    # pprint.pprint(modes_included)
    # print('---------modes included done-----------')

    # ...

if __name__ == "__main__":
    main()