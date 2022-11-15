###################################################################################
'''

10/28/2022 homework: Extract normal modes (n_atoms*3 dimensions) matrix.

An output file containing the matrix is created by this script.
The script necessitates the normal mode frequencies from the optimization/hessian
calculation step .log file. Since GAMESS outputs all 3N frequencies, no need for
3N-6 or 3N-5, it is only 3N.

frequencies.py requires this information:
- 3N normal mode frequency list (cm^-1)
- .log output

It produces these files:
- frequencies.out

Typical usage:
python3 frequencies.py .log_file_here

'''
###################################################################################
# system imports
import os
import json
import sys

# third party imports

# local imports
import helper
from helper import find_string_in_file, readfile
from log_conf import log

# ask the user for number of atoms
#n_atoms = input("Amount of atoms in molecule: ")
#n_freqs = int(3*int(n_atoms))
#n_freqs_nonlinear = int(3*int(n_atoms)-6)
#n_freqs_linear = int(3*int(n_atoms)-5)
optfreq_filename = sys.argv[1]


def extract_harmonic_frequencies(file_path=optfreq_filename):
    '''
    Obtains all the harmonic frequencies (cm**-1) from
    frequency calculation on optimized geometry.
    '''

    #print(optfreq_filename)

    # check file exists
    # read file in as a single string
    data = readfile(file_path, name="GAMESS .log output")
    
    header_str = "NORMAL COORDINATE ANALYSIS IN THE HARMONIC APPROXIMATION"
    footer_str = "THERMOCHEMISTRY"
    natoms_str = " TOTAL NUMBER OF ATOMS                        ="
    nucrep_str = "THE NUCLEAR REPULSION ENERGY"

    header_index = data.index(header_str)
    footer_index = data.index(footer_str)
    natoms_index = data.index(natoms_str)
    nucrep_index = data.index(nucrep_str)

    print(f'Header index: {header_index}')
    print(f'Footer_index: {footer_index}')
    print(f'Natoms index: {natoms_index}')

    natoms = int(data[natoms_index+51:nucrep_index])
    n_freqs = 3*natoms

    #print(optfreq_filename)
    #print(header_str)
    #print(footer_str)
    print(f'THE NUMBER OF ATOMS IS: {natoms}')

    try:
        if (header_str not in data) or (footer_str not in data):
            log.warning(
                    f"File {file_path}\n"
                    f"does not have ({header_str}) and/or ({footer_str}) inside it.\n"
                    "Please check input file, something went wrong."
            )
            raise Exception(helper.generic_exception_string)

        # proceed as normal
        log.info(
            f"GAMESS log file does contain ({header_str=}) and ({footer_str})"
            " so it appears to be correctly formatted."
        )
        '''
        # attention to triple hyphen
        target_freq_str = "Frequencies ---"
        target_rmass_str = "Reduced masses ---"

        # obtain the frequencies
        if (target_freq_str not in data) or (target_rmass_str not in data):
            raise Exception(f"didn't find {target_rmass_str=} or {target_freq_str=}")
        '''
        target_ir_str = "IR INTENSITY:"
        target_sayvetz_str = "REFERENCE ON SAYVETZ CONDITIONS"

'''
        freq_lst = []
        # extract appropriate amount of frequencies
        while (len(freq_lst) != n_freqs):
            freq_start_idx = find_string_in_file(target_ir_str, header_index, footer_index, data, file_path)
            freq_end_idx = find_string_in_file(target_rmass_str, freq_start_idx, footer_index, data, file_path)

            for val in data[freq_start_idx:freq_end_idx].split():
                try:
                    harm_freq = float(val)
                    freq_lst.append(str(harm_freq))
                    log.debug(f"Extracted frequency value of {harm_freq} (cm**-1).")
                except:
                    # normally ValueError because 'Frequencies' and '---'
                    # cause float() to error out. pass helps prevent this.
                    pass

            # have to advance header_index else infinite loop.
            header_index = freq_end_idx

        log.info(f"Succesfully extracted the following frequencies: {freq_lst}")

    '''
    except Exception as e:
        log.warning("Cannot extract frequencies.")
        raise e

    return

def perform_step_6(automate_flag=False):
    """ x """
    log.info("Beginning STEP 6 of the procedure.")
    extract_harmonic_frequencies()

if (__name__ == '__main__'):

    perform_step_6()