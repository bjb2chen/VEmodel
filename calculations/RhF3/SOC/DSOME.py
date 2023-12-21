import sys
import re
import pprint

# Function to extract lines between patterns in a file
def extract_lines_between_patterns(filename, start_pattern, end_pattern):
    selected_lines = []
    reversed_rightmost = []
    collecting = False

    with open(filename, 'r', encoding = "utf-8", errors="replace") as file:
        for line in file.readlines():
            if start_pattern in line:
                collecting = True
                selected_lines.append(line)
            elif end_pattern in line:
                collecting = False
            elif collecting:
                selected_lines.append(line)

    selected_lines.reverse()

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

def extract_DSOME(selected_lines, pattern, nstate):
    DSOME_set = {}
    full_extracted_set = {}
    summed_set_real = {}
    summed_set_imag = {}

    for DSOMEline in selected_lines:
        if "STATE #" in DSOMEline:
            ist = DSOMEline[9:12].strip().replace(" ", "")
            jst = DSOMEline[14:16].strip().replace(" ", "") 
            kst = ist + ' & ' + jst + ',' + DSOMEline[31:33]
            real = float(DSOMEline[48:61].strip().replace(" ", ""))
            imaginary = float(DSOMEline[63:75].strip().replace(" ", ""))
            DSOME_set[kst] = str(real) + "+" + str(imaginary)  # + "I"

    for left_state_idx in range(1, int(nstate)):
        for right_state_idx in range(left_state_idx+1, int(nstate)+1):
            for level_idx in range(1, 3):
                full_extracted_set[left_state_idx, right_state_idx, level_idx] = DSOME_set[f'{left_state_idx} & {right_state_idx}, {level_idx}']
            summed_set_real[left_state_idx, right_state_idx] = float(full_extracted_set[left_state_idx, right_state_idx, 1].split('+')[0]) + \
                                                               float(full_extracted_set[left_state_idx, right_state_idx, 2].split('+')[0])
            summed_set_imag[left_state_idx, right_state_idx] = float(full_extracted_set[left_state_idx, right_state_idx, 1].split('+')[1]) + \
                                                               float(full_extracted_set[left_state_idx, right_state_idx, 2].split('+')[1])
    
    return full_extracted_set, summed_set_real, summed_set_imag

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 DSOME.py <SOC_out_file> <nstate>")
        sys.exit(1)

    filnam = sys.argv[1]
    nstate = sys.argv[2]
    selected_lines = extract_lines_between_patterns(filnam,
        'DIABATIC SPIN-ORBIT MATRIX ELEMENTS',
        'SOC EIG. VALUES and VECTORS IN DIABATS (DIRECT MAX.)'
        )
    #pprint.pprint(selected_lines)
    #DSOME_block = extract_DSOME(selected_lines, 'DIABATIC SPIN-ORBIT MATRIX ELEMENTS')
    #pprint.pprint(DSOME_block)
    extracted = extract_DSOME(selected_lines, 'DIABATIC SPIN-ORBIT MATRIX ELEMENTS', nstate)
    pprint.pprint(extracted)
    #pprint.pprint(extracted)

if __name__ == "__main__":
    main()
