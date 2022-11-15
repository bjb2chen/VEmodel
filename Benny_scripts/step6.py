###################################################################################
'''

Step six of the vibronic model process.

A coupling input file is created. This coupling input file necessitates
info from overdia.out, vertical excitation energies from step 2, and
normal mode frequencies in step 1 output.

coupling.inp requires this information:
- Number of excited states
- Number of normal modes
- Number of normal modes selected and threshold filter (e.g. 1.d-7)
- Vertical excitation energies list (eV)
- Normal mode frequency list (cm^-1)
- overdia.out

It produces these files:
- coupling.out
- operator .op file

'''
###################################################################################
# system imports
import os
import json

# third party imports

# local imports
import helper
from helper import find_string_in_file, readfile
from log_conf import log

with open('master_values.json') as fp:
    value = json.load(fp)

molecule_name = value["molecule_name"]
n_exstates = value["num_excited_states"]
n_freqs = value["num_freqs"]
header_str = "Excitation energies and oscillator strengths:"
footer_str = "Population analysis using the SCF Density."
exstate_str = "Excited State"
eV_str = "eV"

excited_filename = f'{molecule_name}_excited.log'
optfreq_filename = f'{molecule_name}_optfreq.log'

input_filename = f'{molecule_name}_coupling.inp'
output_filename = f'{molecule_name}_coupling.out'
# bash_filename = f'{molecule_name}_coupling.sh'
overdia_opfile_path = "/home/bjb2chen/CHEM494/CODES/overdia_2_opfile.e"


# you can make the same changes to this that I did for step5
def extract_vertical_energies(file_path=excited_filename):
    '''
    Takes the list of vertical excitation energies from
    TD-DFT calculation in step 2.
    '''

    # check file exists
    # read file in as a single string
    data = readfile(file_path, name="TD-DFT output")

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
            f"Excited log file does contain ({header_str=}) and ({footer_str})"
            " so it appears to be correctly formatted."
        )

        log.debug(f"Found ({header_str=}) ({footer_str=}) in {file_path}")

        # get indices
        header_index = data.index(header_str)
        footer_index = data.index(footer_str)

        # for some reason was getting local var referenced before assignment error
        # reassigning them locally worked
        exstate_str = "Excited State"
        eV_str = "eV"

        # obtain the excitation energies
        if (exstate_str not in data) or (eV_str not in data):
            raise Exception(f"didn't find {eV_str=} or {exstate_str=}")

        exstates_lst = []
        for i in range(1, n_exstates+1):
            
            # from testing, cannot simply do 'Excited State   {i}' because
            # if requesting e.g. 15 states,  'Excited State  10' is two spaces instead of 3.
            exstate_str = 'Excited State'
            ex_start_idx = find_string_in_file(exstate_str, header_index, footer_index, data, file_path)
            ex_end_idx = find_string_in_file(eV_str, ex_start_idx, footer_index, data, file_path)

            try:
                if len(data[ex_start_idx:ex_end_idx]) == 46:
                    str_list = data[ex_start_idx:ex_end_idx].split()
                    exstates_lst.append(str_list[4])
                    log.debug(f"For excited state {i}, extracted excitation energy value of {float(exstates_lst[-1])} eV.")
                    header_index = ex_end_idx        # advance header_idx
            except Exception as e:
                log.warning("Unable to isolate 'Excited States' line.")
                raise e

        log.info(f"Succesfully extracted the following excited state energies: {exstates_lst}")

    except Exception as e:
        log.warning("Cannot extract vertical excitation energies.")
        raise e

    return exstates_lst


# you can make the same changes to this that I did for step5
def extract_harmonic_frequencies(file_path=optfreq_filename):
    '''
    Obtains all the harmonic frequencies (cm**-1) from
    frequency calculation on optimized geometry in step 1.
    '''

    # check file exists
    # read file in as a single string
    data = readfile(file_path, name="Gaussian .log output")

    header_str = "Harmonic frequencies"
    footer_str = "Thermochemistry"

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
            f"Gaussian log file does contain ({header_str=}) and ({footer_str})"
            " so it appears to be correctly formatted."
        )
        log.debug(f"Found ({header_str=}) ({footer_str=}) in {file_path}")

        header_index = data.index(header_str)
        footer_index = data.index(footer_str)

        # attention to triple hyphen
        target_freq_str = "Frequencies ---"
        target_rmass_str = "Reduced masses ---"

        # obtain the frequencies
        if (target_freq_str not in data) or (target_rmass_str not in data):
            raise Exception(f"didn't find {target_rmass_str=} or {target_freq_str=}")

        freq_lst = []
        # extract appropriate amount of frequencies
        while (len(freq_lst) != n_freqs):
            freq_start_idx = find_string_in_file(target_freq_str, header_index, footer_index, data, file_path)
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

    except Exception as e:
        log.warning("Cannot extract frequencies.")
        raise e

    return freq_lst


def confirm_step_5_completed_succesfully(file_path=f'{molecule_name}_overdia.out'):
    ''' Make sure (to the best of our knowledge) the step5 completed successfully.
        If it did, 'OUT OF SUMMARY' should be last line.
    '''
    if not os.path.isfile(file_path):
        log.debug(f"The provided path is not valid:\n{file_path}\n")

    try:
        with open(file_path, 'r') as fp:
            data = fp.read()

    except OSError as e:
        log.debug(f"Couldn't open step 4 output file?\n{file_path}\n")
        raise e

    succesful_termination_string = 'OUT OF SUMMARY'

    if (succesful_termination_string not in data):
        raise Exception(
            "It appears that step 5 did not complete successfully"
            f"Please check the output file: \n{file_path}\n"
        )
    return


# you can make the same changes to this that I did for step5
def compose_coupling_input():
    line1 = f"{str(n_exstates):<25}! number of states"
    line2 = f"{str(n_freqs):<25}! number of modes"
    line3 = f"{(str(n_freqs) + '  1.d-7'):<25}! number of modes selected and threshold < are taken as zero"
    line4 = 'VERTICAL ENERGIES'
    line5 = '\n'.join(extract_vertical_energies(f"{molecule_name}_excited.log"))
    line6 = 'FREQUENCIES'
    line7 = '\n'.join(extract_harmonic_frequencies(f"{molecule_name}_optfreq.log"))
    line8 = f" '{molecule_name}_overdia.out'                 ! output of overdia"
    line9 = f" 'op_{molecule_name}{n_freqs}Q_{n_exstates}st.op'            ! op file for MCTDH (QUANTICS)"

    output = [line1, line2, line3, line4, line5, line6, line7, line8, line9]
    output = '\n'.join(output)

    try:
        with open(f'{molecule_name}_coupling.inp', 'w') as fp:
            fp.write(output)

        log.info(f"Generated the {molecule_name}_coupling.inp file")

    except Exception as e:
        log.warning(f"Cannot generate {molecule_name}_coupling.inp file.")
        raise e

    return


def perform_step_6(automate_flag=False):
    """ x """
    log.info("Beginning STEP 6 of the procedure.")
    confirm_step_5_completed_succesfully()
    log.info("Step 5 completed successfully.")
    compose_coupling_input()
    command = f'{overdia_opfile_path} < {input_filename} > {output_filename}'

    if False:  # not done yet
        extract_vertical_energies()
        extract_harmonic_frequencies()

    if automate_flag:
        return command

    else:
        os.system(command)
        log.info(
            "Please proceed and check if files "
            f"{molecule_name}_coupling.out and op_{molecule_name}{n_freqs}Q_{n_exstates}st.op "
            "have correct output.")
        return


if (__name__ == '__main__'):

    # code to run when file called directly
    with open('master_values.json') as fp:
        value = json.load(fp)

    molecule_name = value["molecule_name"]
    n_exstates = value["num_excited_states"]
    n_freqs = value["num_freqs"]

    perform_step_6()
