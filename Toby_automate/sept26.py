import sys
import pprint

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
                nrmmod[ixyz, imode] = disp
                print(nrmmod[ixyz, imode], disp)

                if ixyz == 1:
                    cutini = (icolumn - 1) * 12
                    cutfnl = icolumn * 12
                    freq = lines_freq[igroup - 1][cutini:cutfnl].lstrip()
                    freqcm[imode] = freq
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
                    nrmmod[ixyz, imode] = disp
                    print(nrmmod[ixyz, imode], disp)

                    if ixyz == 1:
                        cutini = (icolumn - 1) * 12
                        cutfnl = icolumn * 12
                        freq = lines_freq[-1][cutini:cutfnl].lstrip()
                        freqcm[imode] = freq
                        print("frequency:", imode, freqcm[imode])

    # # Check if imode is in modes_excluded and exclude if necessary
    # for imode in modes_excluded:
    #     if imode in freqcm:
    #         del freqcm[imode]

    #Print all frequencies
    for imode in range(1, ndim + 1):
        print("frequency:", imode, freqcm[imode].lstrip(), "CM-1")

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
                refcoord[icomp] = coord
                print(refcoord[icomp])

    return atmlst, chrglst, refcoord

def filter_modes(excluded_set, ndim):
    modes_included = []
    print("NUMBER OF EXCLUDED MODES:", len(excluded_set))
    print("They are modes:", *excluded_set)

    # Iterate through modes to check inclusion
    for imode in range(1, ndim + 1):
        if imode not in excluded_set:
            modes_included.append(imode)
            print(len(modes_included), imode)

    print("Number of Modes Included:", len(modes_included))        

    return modes_included

def diabatization_placeholder():
    distcoord_plus = {}
    distcoord_minus = {}
    distcoord_plus_x2 = {}
    distcoord_minus_x2 = {}
    distcoord_pp = {}
    distcoord_pm = {}
    distcoord_mp = {}
    distcoord_mm = {}

    return

def main():
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <hessout_file>")
        sys.exit(1)

    hessout = sys.argv[1]
    filnam = "nh3cat_ccd_gmcpt_7o7e_C3vinC1_3st_diab"

    natoms = get_number_of_atoms(hessout)
    ndim = natoms * 3
    #ndim = 15 note that you can only test like this if you have proper hessout structure
    ngroup = ndim // 5
    nleft = ndim % 5
    print("Dimension of all xyz coordinates:", ndim)
    print("# of atoms:", natoms)
    print(ngroup, nleft)

    #modes_excluded = [1, 2, 3, 4, 5, 6]
    modes_excluded = [1, 2, 4, 7, 12]

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
    atmlst, chrglst, refcoord = read_reference_structure("oct3_ref_structure")
    modes_included = filter_modes(modes_excluded, ndim) 

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