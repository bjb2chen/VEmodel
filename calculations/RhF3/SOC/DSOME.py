import sys
import pprint
import subprocess
import os
import shutil
import re
import json

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
    append_I = {}

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
            append_I[left_state_idx, right_state_idx] = str(summed_set_imag[left_state_idx, right_state_idx]) + "I"
    
    return full_extracted_set, summed_set_real, summed_set_imag, append_I

def mctdh(filnam, modes_included, **kwargs):
    nmodes = len(modes_included)

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
    str8 = f'{filnam} {nstate} states + ' + str(nmodes) + ' modes'
    strlst = [str1, str2, str8, str3, str4, str5, str6, str7]

    with open('mctdh.op', 'w') as mctdh_file:
        for idx in strlst:
            mctdh_file.write(idx+'\n')

    # Write SOC ELECTRONIC COUPLING AT REFERENCE STRUCTURE
    mctdh_file.write("-----------------------------------------\n")
    mctdh_file.write("# SOC REAL AND IMAG COUPLINGS (OFF DIAGONAL ELEMENT)\n")
    mctdh_file.write("-----------------------------------------\n")
    # for ist in range(1, nstate + 1):
    #     mctdh_file.write(f"v{ist}  |1 S{ist}&{ist}\n")
    for ist in range(1, nstate + 1):
        jlast = ist - 1
        for jst in range(1, jlast + 1):
            mctdh_file.write(f"SOr_{jst}_{ist} = {summed_set_real[jst, ist]:.16f}, CM-1\n")
            mctdh_file.write(f"SOi_{jst}_{ist} = {summed_set_imag[jst, ist]:.16f}, CM-1\n")

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
    
        # # Write KINETIC OPERATOR FOR NORMAL MODES
        # mctdh_file.write("# KINETIC OPERATOR FOR NORMAL MODES\n")
        # mctdh_file.write("-----------------------------------------\n")
        # for imode_include in range(1, nmodes + 1):
        #     mode_count = imode_include + 1
        #     mctdh_file.write(f"w_{modes_included[imode_include]}   |{mode_count} KE\n")
    
        # # Write HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES
        # mctdh_file.write("-----------------------------------------\n")
        # mctdh_file.write("# HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES\n")
        # mctdh_file.write("-----------------------------------------\n")
        # for imode_include in range(1, nmodes + 1):
        #     mode_count = imode_include + 1
        #     mctdh_file.write(f"0.5*w_{modes_included[imode_include]}   |{mode_count}  q^2\n")
    
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
    
        # # Write LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS
        # mctdh_file.write("-----------------------------------------\n")
        # mctdh_file.write("# LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS\n")
        # mctdh_file.write("-----------------------------------------\n")
        # for kmode in range(1, nmodes + 1):
        #     imode = modes_included[kmode]
        #     kmode_count = kmode + 1
        #     for ist in range(1, nstate + 1):
        #         mctdh_file.write(f"l{ist}_m{imode} |1 S{ist}&{ist} |{kmode_count} q\n")
        #         mctdh_file.write(f"q{ist}_m{imode} |1 S{ist}&{ist} |{kmode_count} q^2\n")
    
        # # Write LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS
        # mctdh_file.write("-----------------------------------------\n")
        # mctdh_file.write("# LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        # mctdh_file.write("-----------------------------------------\n")
        # for kmode in range(1, nmodes + 1):
        #     imode = modes_included[kmode]
        #     kmode_count = kmode + 1
        #     for ist in range(1, nstate + 1):
        #         jlast = ist - 1
        #         for jst in range(1, jlast + 1):
        #             mctdh_file.write(f"l{jst}{ist}_m{imode} |1 S{jst}&{ist} |{kmode_count} q\n")
        #             mctdh_file.write(f"q{jst}{ist}_m{imode} |1 S{jst}&{ist} |{kmode_count} q^2\n")
    
        # # Write BILINEAR DIAGONAL VIBRONIC COUPLINGS
        # mctdh_file.write("-----------------------------------------\n")
        # mctdh_file.write("# BILINEAR DIAGONAL VIBRONIC COUPLINGS\n")
        # mctdh_file.write("-----------------------------------------\n")
        # for kmode in range(1, nmodes + 1):
        #     imode = modes_included[kmode]
        #     kmode_count = kmode + 1
        #     lmode_last = kmode - 1
        #     for lmode in range(1, lmode_last + 1):
        #         jmode = modes_included[lmode]
        #         lmode_count = lmode + 1
        #         for ist in range(1, nstate + 1):
        #             mctdh_file.write(f"b{ist}_m{imode}_m{jmode} |1 S{ist}&{ist} |{lmode_count} q |{kmode_count} q\n")
    
        # # Write BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS
        # mctdh_file.write("-----------------------------------------\n")
        # mctdh_file.write("# BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        # mctdh_file.write("-----------------------------------------\n")
        # for kmode in range(1, nmodes + 1):
        #     imode = modes_included[kmode]
        #     kmode_count = kmode + 1
        #     lmode_last = kmode - 1
        #     for lmode in range(1, lmode_last + 1):
        #         jmode = modes_included[lmode]
        #         lmode_count = lmode + 1
        #         for ist in range(1, nstate + 1):
        #             jlast = ist - 1
        #             for jst in range(1, jlast + 1):
        #                 mctdh_file.write(f"b{jst}{ist}_m{imode}_m{jmode} |1 S{jst}&{ist} |{lmode_count} q |{kmode_count} q\n")

        # Write SOC LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# SOC LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):  
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            for ist in range(1, nstate + 1):
                jlast = ist - 1
                for jst in range(1, jlast + 1):
                    #mctdh_file.write(f"l{jst}{ist}_m{imode} |1 S{jst}&{ist} |{kmode_count} q\n")
                    #mctdh_file.write(f"q{jst}{ist}_m{imode} |1 S{jst}&{ist} |{kmode_count} q^2\n")
                    mctdh_file.write(f"I*lSOi_{jst}_{ist}_m{imode} |1 Z{jst}&{ist}|{kmode_count} q\n")
                    #note to self: the I* is performing ARITHMETIC on SOr_{jst}_{ist} prepared earlier, does that mean we neeed to remove the l and _m{imode}
                    mctdh_file.write(f"-I*lSOi_{jst}_{ist}_m{imode} |1 Z{ist}&{jst}|{kmode_count} q\n")
                    mctdh_file.write(f"I*qSOi_{jst}_{ist}_m{imode} |1 Z{jst}&{ist}|{kmode_count} q^2\n")
                    mctdh_file.write(f"-I*qSOi_{jst}_{ist}_m{imode} |1 Z{ist}&{jst}|{kmode_count} q^2\n")

        # Write SOC BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS
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
                        #mctdh_file.write(f"b{jst}{ist}_m{imode}_m{jmode} |1 S{jst}&{ist} |{lmode_count} q |{kmode_count} q\n")
                        mctdh_file.write(f"I*SOi_b{jst}{ist}_m{imode}_m{jmode} |1 Z{jst}&{ist} |{lmode_count} q |{kmode_count} q\n")
                        mctdh_file.write(f"-I*SOi_b{jst}{ist}_m{imode}_m{jmode} |1 Z{ist}&{jst} |{lmode_count} q |{kmode_count} q\n")

        # Close the file
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("\nend-hamiltonian-section\n\nend-operator\n")
    
    return

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 DSOME.py <SOC_refG_file>")
        sys.exit(1)

    filnam = sys.argv[1]
    # nstate = sys.argv[2]
    with open(f'{filnam}', 'r', encoding = "utf-8", errors="replace") as refGout_file:
        for line in refGout_file:
            if '# of states in CI      = ' in line:
                nstate = int(line.split('=')[1]) # this will hopefully grab the last line

    #Set conversion constants
    qsize = 0.05
    ha2ev = 27.2113961318
    wn2ev = 0.000123981
    wn2eh = 0.00000455633
    ang2br = 1.889725989
    amu2me = 1822.88839 

    selected_lines = extract_lines_between_patterns(filnam,
        'DIABATIC SPIN-ORBIT MATRIX ELEMENTS',
        'SOC EIG. VALUES and VECTORS IN DIABATS (DIRECT MAX.)'
        )
    #pprint.pprint(selected_lines)
    #DSOME_block = extract_DSOME(selected_lines, 'DIABATIC SPIN-ORBIT MATRIX ELEMENTS')
    #pprint.pprint(DSOME_block)
    extracted = extract_DSOME(selected_lines, 'DIABATIC SPIN-ORBIT MATRIX ELEMENTS', nstate)
    #pprint.pprint(extracted)
    full_extracted_set, summed_set_real, summed_set_imag, append_I = extracted[0], extracted[1], extracted[2], extracted[3]
    pprint.pprint(full_extracted_set)
    pprint.pprint(summed_set_real)
    pprint.pprint(append_I)

    modes_included = {1: 7, 2: 8, 3: 9, 4: 10, 5: 11, 6: 12}
    make_mctdh = mctdh(filnam, modes_included, qsize=qsize, ha2ev=ha2ev, wn2ev=wn2ev, wn2eh=wn2eh, ang2br=ang2br, amu2me=amu2me, nstate=nstate, summed_set_real=summed_set_real, summed_set_imag=summed_set_imag)

if __name__ == "__main__":
    main()
