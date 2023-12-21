import sys
import re
import pprint

# def get_number_of_atoms(hessout):
#     with open(hessout, 'r') as hess_file:
#         for line in hess_file:
#             if ' TOTAL NUMBER OF ATOMS' in line:
#                 natoms = int(line.split('=')[1])
#                 return natoms

# Function to extract lines between patterns in a file
def extract_lines_between_patterns(filename, start_pattern, end_pattern):
    selected_lines = []
    reversed_rightmost = []
    collecting = False

    with open(filename, 'r') as file:
        for line in file.readlines():
            if start_pattern in line:
                collecting = True
                selected_lines.append(line)
            elif end_pattern in line:
                collecting = False
            elif collecting:
                selected_lines.append(line)

    selected_lines.reverse()
    #pprint.pprint('New set\n')
    #pprint.pprint(selected_lines)

    for line in selected_lines:
        # if end_pattern in line:
        # end_pattern line not included in intial collection
        collecting = True
        reversed_rightmost.append(line)
        if start_pattern in line:
            collecting = False
        elif collecting:
            reversed_rightmost.append(line)

    reversed_rightmost.reverse()

    return reversed_rightmost

# def extract_diabatic_energy(file_path, pattern):
#     with open(file_path, 'r') as file:
#         for line in reversed(file.readlines()):
#             match = re.search(pattern, line)
#             if match:
#                 return float(line[44:62].strip().replace(" ", ""))

# def extract_coupling_energy(file_path, pattern):
#     with open(file_path, 'r') as file:
#         for line in reversed(file.readlines()):
#             match = re.search(pattern, line)
#             if match:
#                return float(line[62:].strip().replace(" ", ""))

def extract_DSOME(selected_lines, pattern):
    DSOME_set = {}
    # for ist in range(1, nstate + 1):
    for DSOMEline in selected_lines:
        if "STATE #" in DSOMEline:
            ist = DSOMEline[9:16] + ',' + DSOMEline[31:33]
            DSOME_set[ist] = DSOMEline[48:61].strip().replace(" ", "") + "+" + \
                             DSOMEline[63:77].strip().replace(" ", "")
    print(DSOME_set['12 & 13, 2'])
    return DSOME_set
    # with open(file_path, 'r') as file:
    #     for line in reversed(file.readlines()):
    #         match = re.search(pattern, line)
    #         if match:
    #             return line[48:61].strip().replace(" ", "") + "+" + \
    #                    line[63:77].strip().replace(" ", "")  

# def read_freq_values(selected_lines):
#     freq_value_set = []

#     for freqline in selected_lines:
#         if "FREQUENCY:" in freqline:
#             freq_value_set.append(freqline[18:])

#     return freq_value_set

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 your_script.py <SOC_out_file>")
        sys.exit(1)

    filnam = sys.argv[1]
    selected_lines = extract_lines_between_patterns(filnam,
        'DIABATIC SPIN-ORBIT MATRIX ELEMENTS',
        'SOC EIG. VALUES and VECTORS IN DIABATS (DIRECT MAX.)'
        )
    pprint.pprint(selected_lines)
    DSOME_block = extract_DSOME(selected_lines, 'DIABATIC SPIN-ORBIT MATRIX ELEMENTS')
    pprint.pprint(DSOME_block)

if __name__ == "__main__":
    main()