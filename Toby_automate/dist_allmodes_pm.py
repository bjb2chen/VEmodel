import sys
import pprint
import subprocess
import shutil

# Function to get the number of atoms from the HESS output file
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

    with open("oct3_mode.dat", "r") as mode_file:
        lines_mode = mode_file.readlines()

    with open("oct3_freq.dat", "r") as freq_file:
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
                print(ixyz, imode, end=" ")
                cutini = (icolumn - 1) * 12
                cutfnl = icolumn * 12

                disp = lines_mode[line - 1][cutini:cutfnl].lstrip()
                nrmmod[ixyz, imode] = float(disp) # no need to suppress sci notation, GMS is OK with E^...
                print(nrmmod[ixyz, imode], disp)

                if ixyz == 1:
                    cutini = (icolumn - 1) * 12
                    cutfnl = icolumn * 12
                    freq = lines_freq[igroup - 1][cutini:cutfnl].lstrip()
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
                    print(ixyz, imode, end=" ")
                    cutini = (icolumn - 1) * 12
                    cutfnl = icolumn * 12

                    disp = lines_mode[line - 1][cutini:cutfnl].lstrip()
                    nrmmod[ixyz, imode] = float(disp)
                    print(nrmmod[ixyz, imode], disp)

                    if ixyz == 1:
                        cutini = (icolumn - 1) * 12
                        cutfnl = icolumn * 12
                        freq = lines_freq[-1][cutini:cutfnl].lstrip()
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


#Do diabatization calculation at the reference nondistorted structure.
#This calculation shall be a repetition of a calcualtion in preparing temp.inp
def refG_calc(refgeo, filnam):
    # Check if the calculation has already been run
    grep_process = subprocess.run(["grep", "grace", f"{filnam}_refG.out"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if grep_process.returncode != 0:
        print("Run calculation at the undistorted reference structure")

        shutil.copy("temp.inp", f"{filnam}_refG.inp")

        with open(f"{filnam}_refG.inp", "a") as inp_file:
            with open(refgeo, "r") as ref_structure:
                inp_file.write(ref_structure.read())

        with open(f"{filnam}_refG.inp", "a") as inp_file:
            inp_file.write(" $END\n")

        # Run the calculation (you may need to customize this command based on your setup)
        subprocess.run(["./subgam.diab", f"{filnam}_refG.inp", "4", "0", "1"])

        # might want to do a sleep(1 min) here? What if refG.inp calculation fails?

        print("Calculation at the reference structure is done.")
    else:
        print("Calculation at the reference structure has already been done.")

    return

def diabatization(modes_included, freqcm, ndim, refcoord, nrmmod, natom, atmlst, chrglst, filnam):

    distcoord_plus = {}
    distcoord_minus = {}
    distcoord_plus_x2 = {}
    distcoord_minus_x2 = {}
    distcoord_pp = {}
    distcoord_pm = {}
    distcoord_mp = {}
    distcoord_mm = {}

    #Loop over all considered modes and do + and - displacements
    #Set step size in reduced dimensionless coordinates
    qsize = 0.05

    # Set conversion constants
    ha2ev = 27.2113961318
    wn2ev = 0.000123981
    wn2eh = 0.00000455633
    ang2br = 1.889725989
    amu2me = 1822.88839

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
        for iatom in range(1, natom + 1):
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
                    with open(f'dist_structure_{p_or_m}{suffix}', 'r') as dist_file:
                        inp_file.write(dist_file.read())
                        inp_file.write(' $END')

                # Check if the calculation is done already
                grace1 = subprocess.run(["grep", "grace", f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.out'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                #if not os.path.exists(f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.out'):
                if grace1.returncode != 0:
                    print(f"Running calculations for {filnam}_mode{imode}_{displacement}{qsize}{suffix}")
                    try:
                        subprocess.run(['./subgam.diab', f'{filnam}_mode{imode}_{displacement}{qsize}{suffix}.inp', '4', '0', '1'])
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
            for iatom in range(1, natom + 1):
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
                        with open(f'dist_structure_{suffix1}', 'r') as dist_file:
                            inp_file.write(dist_file.read())
                        inp_file.write(' $END ')
 
                    # Check if the calculation is done already
                    output_filename = f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.out'
                    grace2 = subprocess.run(["grep", "grace", output_filename])
                    if grace2.returncode != 0:
                    #if not os.path.exists(output_filename):
                        print(f"Running calculations for {output_filename}!")
                        try:
                            subprocess.run(['./subgam.diab', f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp', '4', '0', '1'])
                        except Exception as e:
                            print(f"Error running diabatization calculation: {str(e)}")
                    else:
                        print(f"{output_filename} is already done.")

#Now we move on to extract vibronic coupling constants using finite difference
#and write the data in an mctdh operator file

def mctdh(filnam, modes_included):
    try:
        subprocess.run(['rm', '-f', 'mctdh.op'])
    except Exception as e:
        print(f"Error deleting {'mctdh.op'}: {str(e)}")

    #Heading for mctdh.op
    str1 = "OP_DEFINE-SECTION"
    str2 = "title"

    # nstate=`grep '# of states in CI      = ' "$filnam"_refG.out|tail -1|cut -d'=' -f2`
    with open(f'{filnam}_refG.out', 'r') as refGout_file:
        for line in refGout_file:
            if '# of states in CI      = ' in line:
                nstate = int(line.split('=')[1]) # this will hopefully grab the last line
    
    str3 = "end-title "
    str4 = "end-op_define-section"
    str5 = ""
    # lines 482,483
    str6 = "PARAMETER-SECTION"
    str7 = ""
    str8 = f'{filnam} {nstate} states +' + str(len(modes_included)) + ' modes'
    strlst = [str1, str2, str8, str3, str4, str5, str6, str7]

    with open('mctdh.op', 'w') as mctdh_file:
        for idx in strlst:
            mctdh_file.write(idx+'\n')

    return
 
    
def main():
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <hessout_file>")
        sys.exit(1)

    hessout = sys.argv[1]
    filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"
    refgeo = "oct3_ref_structure"

    natoms = get_number_of_atoms(hessout)
    ndim = natoms * 3
    #ndim = 15 note that you can only test like this if you have proper hessout structure
    ngroup = ndim // 5
    nleft = ndim % 5
    print("Dimension of all xyz coordinates:", ndim)
    print("# of atoms:", natoms)
    print(ngroup, nleft)

    modes_excluded = [1, 2, 3, 4, 5, 6]
    #modes_excluded = [1, 2, 4, 7, 12]

    selected_lines = extract_lines_between_patterns(hessout, 
        "FREQUENCIES IN CM",
        "REFERENCE ON SAYVETZ CONDITIONS"
        )
    freq_value_set = read_freq_values(selected_lines)
    filtered_set = read_mode_values(selected_lines)
    
    with open('oct3_mode.dat', 'w') as output_file:
        output_file.writelines(filtered_set)

    with open('oct3_freq.dat', 'w') as output_file:
        output_file.writelines(freq_value_set)

    nrmmod, freqcm = process_mode_freq(natoms, ndim, ngroup, nleft)
    atmlst, chrglst, refcoord = read_reference_structure(refgeo)
    modes_included = filter_modes(modes_excluded, ndim)

    repetition = refG_calc(refgeo, filnam)
    diabatize = diabatization(modes_included, freqcm, ndim, refcoord, nrmmod, natoms, atmlst, chrglst, filnam)
    make_mctdh = mctdh(filnam, modes_included)

    qsize = 0.05
    #Set conversion constants
    ha2ev = 27.2113961318
    wn2ev = 0.000123981
    wn2eh = 0.00000455633
    ang2br = 1.889725989
    amu2me = 1822.88839 

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