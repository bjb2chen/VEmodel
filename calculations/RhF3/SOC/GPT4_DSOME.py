import sys
import re
import os
import subprocess

# Function to extract lines between patterns in a file
def extract_lines_between_patterns(filename, start_pattern, end_pattern, encoding='ISO-8859-1'):
    collecting = False
    selected_lines = []
    with open(filename, 'r', encoding=encoding, errors='replace') as file:
        for line in file:
            if start_pattern in line:
                collecting = True
            elif end_pattern in line:
                collecting = False
            if collecting:
                selected_lines.append(line)
    return selected_lines[::-1]

# Function to extract all levels of DSOME
def extract_all_levels_DSOME(file_path, pattern, nstate, encoding='ISO-8859-1'):
    all_extracted_data = []
    with open(file_path, 'r', encoding=encoding, errors='replace') as file:
        lines = file.readlines()
    occurrence_indices = [loc for loc, val in enumerate(lines) if pattern in val]
    for index in occurrence_indices:
        DSOME_set = {}
        full_extracted_set = {}
        summed_set_real = {}
        summed_set_imag = {}
        subsequent_lines = lines[index:]
        state_lines = [line.strip() for line in subsequent_lines if line.startswith(" STATE #")]
        for line in state_lines:
            match = re.search(r"STATE # (\d+) & (\d+)'S.*LEVEL (\d+)-El.*=\s*([-\d.]+)\s*\+\s*([-\d.]+)I", line)
            if match:
                ist, jst, level, real, imaginary = match.groups()
                kst = f"{ist} & {jst}, {level}"
                DSOME_set[kst] = complex(float(real), float(imaginary))
        for left_state_idx in range(1, nstate):
            for right_state_idx in range(left_state_idx + 1, nstate + 1):
                for level_idx in range(1, 3):
                    key = f"{left_state_idx} & {right_state_idx}, {level_idx}"
                    if key in DSOME_set:
                        full_extracted_set[key] = DSOME_set[key]
                        if (left_state_idx, right_state_idx) not in summed_set_real:
                            summed_set_real[(left_state_idx, right_state_idx)] = DSOME_set[key].real
                            summed_set_imag[(left_state_idx, right_state_idx)] = DSOME_set[key].imag
                        else:
                            summed_set_real[(left_state_idx, right_state_idx)] += DSOME_set[key].real
                            summed_set_imag[(left_state_idx, right_state_idx)] += DSOME_set[key].imag
        all_extracted_data.append((full_extracted_set, summed_set_real, summed_set_imag))
    return all_extracted_data

# Function to create MCTDH operator files
def mctdh(filnam, modes_included, nstate, summed_set_real, summed_set_imag):
    try:
        os.remove('mctdh.op')
    except Exception as e:
        print(f"Error deleting 'mctdh.op': {str(e)}")
    strlst = [
        "OP_DEFINE-SECTION",
        "title",
        f'{filnam} {nstate} states + {len(modes_included)} modes',
        "end-title",
        "end-op_define-section",
        "",
        "PARAMETER-SECTION",
        ""
    ]
    with open('mctdh.op', 'w') as mctdh_file:
        mctdh_file.writelines(f"{line}\n" for line in strlst)
        mctdh_file.write("# SOC REAL AND IMAG COUPLINGS (OFF DIAGONAL ELEMENT)\n")
        for ist in range(1, nstate + 1):
            for jst in range(1, ist):
                mctdh_file.write(f"SOr_{jst}_{ist} = {summed_set_real[(jst, ist)]:.16f}, CM-1\n")
                mctdh_file.write(f"SOi_{jst}_{ist} = {summed_set_imag[(jst, ist)]:.16f}, CM-1\n")
        mctdh_file.write("end-parameter-section\n")
        mctdh_file.write("# ELECTRONIC COUPLING AT REFERENCE STRUCTURE\n")
        for ist in range(1, nstate + 1):
            mctdh_file.write(f"v{ist}  |1 S{ist}&{ist}\n")
        for ist in range(1, nstate + 1):
            jlast = ist - 1
            for jst in range(1, jlast + 1):
                mctdh_file.write(f"v{jst}{ist}  |1 S{jst}&{ist}\n")
        mctdh_file.write("# SOC LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        for kmode in range(1, len(modes_included) + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            for ist in range(1, nstate + 1):
                jlast = ist - 1
                for jst in range(1, jlast + 1):
                    mctdh_file.write(f"I*lSOi_{jst}_{ist}_m{imode} |1 Z{jst}&{ist}|{kmode_count} q\n")
                    mctdh_file.write(f"-I*lSOi_{jst}_{ist}_m{imode} |1 Z{ist}&{jst}|{kmode_count} q\n")
                    mctdh_file.write(f"I*qSOi_{jst}_{ist}_m{imode} |1 Z{jst}&{ist}|{kmode_count} q^2\n")
                    mctdh_file.write(f"-I*qSOi_{jst}_{ist}_m{imode} |1 Z{ist}&{jst}|{kmode_count} q^2\n")
        mctdh_file.write("# BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        for kmode in range(1, len(modes_included) + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            lmode_last = kmode - 1
            for lmode in range(1, lmode_last + 1):
                jmode = modes_included[lmode]
                lmode_count = lmode + 1
                for ist in range(1, nstate + 1):
                    jlast = ist - 1
                    for jst in range(1, jlast + 1):
                        mctdh_file.write(f"I*SOi_b{jst}{ist}_m{imode}_m{jmode} |1 Z{jst}&{ist} |{lmode_count} q |{kmode_count} q\n")
                        mctdh_file.write(f"-I*SOi_b{jst}{ist}_m{imode}_m{jmode} |1 Z{ist}&{jst} |{lmode_count} q |{kmode_count} q\n")
        mctdh_file.write("\nend-hamiltonian-section\n\nend-operator\n")

def find_nstate(file_path, pattern='# of states in CI      = ', encoding='ISO-8859-1'):
    with open(file_path, 'r', encoding=encoding, errors='replace') as file:
        for line in file:
            if pattern in line:
                return int(line.split('=')[1].strip())
    return None  # Return None if the pattern is not found

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <SOC_refG_file>")
        sys.exit(1)

    filnam = sys.argv[1]
    nstate = find_nstate(filnam)
    if nstate is None:
        print("Unable to determine the number of states. Check the file format.")
        sys.exit(1)

    pattern = 'DIABATIC SPIN-ORBIT MATRIX ELEMENTS'
    all_levels_data = extract_all_levels_DSOME(filnam, pattern, nstate)

    if all_levels_data:
        full_extracted_set, summed_set_real, summed_set_imag = all_levels_data[0]
        modes_included = {1: 7, 2: 8, 3: 9, 4: 10, 5: 11, 6: 12}
        mctdh(filnam, modes_included, nstate, summed_set_real, summed_set_imag)
    else:
        print("No data extracted or an error occurred")

if __name__ == "__main__":
    main()