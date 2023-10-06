import sys
import pprint

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
def extract_filtered_set(selected_lines):
    filtered_set = []

    for idx, modeline in enumerate(selected_lines):
        if len(modeline) > 3 and modeline[2].isdigit():
            # filtered_set.append(selected_lines[idx][20:])
            # filtered_set.append(selected_lines[idx+1][20:])
            # filtered_set.append(selected_lines[idx+2][20:])
            for i in range(3):
                filtered_set.append(selected_lines[idx + i][20:])

    return filtered_set

# Function to get the number of atoms from the HESS output file
def get_number_of_atoms(hessout):
    with open(hessout, 'r') as hess_file:
        for line in hess_file:
            if ' TOTAL NUMBER OF ATOMS' in line:
                natoms = int(line.split('=')[1])
                return natoms

def process_data(natoms, ndim, ngroup, nleft, modes_excluded):
    nrmmod = {}
    freqcm = {}

    with open("oct3_mode.dat", "r") as mode_file:
        lines_mode = mode_file.readlines()

    with open("oct3_freq.dat", "r") as freq_file:
        lines_freq = freq_file.readlines()

    for igroup in range(1, ngroup + 1):
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

                disp = lines_mode[line - 1][cutini:cutfnl]
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

                for icolumn in range(1, nleft + 1):
                    imode = ngroup * 5 + icolumn
                    print(ixyz, imode, end=" ")
                    cutini = (icolumn - 1) * 12
                    cutfnl = icolumn * 12

                    disp = lines_mode[line - 1][cutini:cutfnl]
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


def main():
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <hessout_file>")
        sys.exit(1)

    hessout = sys.argv[1]

    natoms = get_number_of_atoms(hessout)
    ndim = natoms * 3
    ndim = 15
    ngroup = ndim // 5
    nleft = ndim % 5

    modes_excluded = [1, 2, 3, 4, 5, 6]

    modes_excluded = [2, 3, 5, 7, 8, 9, 13]

    filtered_set = extract_lines_between_patterns(hessout, "FREQUENCIES IN CM", "REFERENCE ON SAYVETZ CONDITIONS")
    freq_value_set = read_freq_values(filtered_set)
    filtered_set = extract_filtered_set(filtered_set)
    
    with open('oct3_mode.dat', 'w') as output_file:
        output_file.writelines(filtered_set)

    with open('oct3_freq.dat', 'w') as output_file:
        output_file.writelines(freq_value_set)

    nrmmod, freqcm = process_data(natoms, ndim, ngroup, nleft, modes_excluded)

    pprint.pprint(nrmmod)
    print('---------nrm mod done-----------')
    pprint.pprint(freqcm)
    print('---------freqcm done-----------')
    #pprint.pprint(selected_lines)
    print('---------selected_lines done-----------')
    pprint.pprint(filtered_set)
    print('---------filtered_set done-----------')
    pprint.pprint(freq_value_set)
    print('---------freq_value_set done-----------')
    # pprint.pprint(atmlst)
    # print('---------atmlst done-----------')
    # pprint.pprint(chrglst)
    # print('---------chrglst done-----------')
    # pprint.pprint(refcoord)
    # print('---------refcoord done-----------')


    # ... (rest of your code)

if __name__ == "__main__":
    main()


# def process_data(natoms, ndim, ngroup, modes_excluded):
#     nrmmod = {}
#     freqcm = {}

#     for igroup in range(1, ngroup + 1, 1):
#         iniline = (igroup - 1) * ndim + 1
#         endline = iniline + ndim - 1

#         for ixyz in range(1, natoms + 1):
#             for icolumn in range(1, 6, 1):
#                 imode = (igroup - 1) * 5 + icolumn

#                 if imode in modes_excluded:
#                     continue

#                 cutini = (icolumn - 1) * 12
#                 cutfnl = icolumn * 12

#                 with open("oct3_mode.dat", "r") as mode_file:
#                     lines_mode = mode_file.readlines()
#                     disp = lines_mode[line - 1][cutini:cutfnl]

#                 nrmmod[ixyz, imode] = disp

#                 if ixyz == 1:
#                     cutini = (icolumn - 1) * 12
#                     cutfnl = icolumn * 12
#                     with open("oct3_freq.dat", "r") as freq_file:
#                         lines_freq = freq_file.readlines()
#                         freq = lines_freq[imode - 1][cutini:cutfnl].lstrip()
#                         freqcm[imode] = freq

#     return nrmmod, freqcm