###################################################################################
'''

Perform excited state calculation to get vertical excitation energy levels needed
for a later step.

It requires this information:
- Gaussian optimization and frequency output file
- submit script

It produces these files:
= Gaussian excited states input file
- Gaussian excited states output .log file

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
from step1 import set_mem_and_processor

with open('master_values.json') as fp:
    value = json.load(fp)

molecule_name = value["molecule_name"]

output_file_path = f"{molecule_name}_excited.com"

header_str = "Symbolic Z-matrix:"
footer_str = "Distance matrix"
charge_str = "Charge =  "
multiplicity_str = " Multiplicity = "
opt_completed_str = "Optimization completed."
input_orientation_str = "Input orientation:"
distance_matrix_str = "Distance matrix"
standard_orientation_str = "Standard orientation"
new_line_string = '\n'

N_atoms = value["num_atoms"]
N_freqs = value["num_freqs"]
N_exstates = value["num_excited_states"]


def create_excited_input(file_path, nof_atoms, nof_states):
    '''
    Gaussian optimization/frequency calculation output as file_path
    Find the target string 'Symbolic Z-matrix:'
    Extract out the charge, multiplicity, and input orientation to create
    the excited state calculation input file.
    '''
    log.info(f"User has specified {nof_atoms} atoms in {molecule_name} molecule.")
    log.info(f"User has specified {nof_states} states for {molecule_name} molecule.")

    data = readfile(file_path, name="Gaussian log")

    try:
        header_index = find_string_in_file(header_str, 0, len(data), data, file_path)
        footer_index = find_string_in_file(footer_str, 0, len(data), data, file_path)

        log.info(
            f"Gaussian log file does contain ({header_str=}) and ({footer_str})"
            " so it appears to be correctly formatted."
        )

        # convert the string e.g. 'Charge =  0 Multiplicity = 1'

        # charge
        c_start_idx = find_string_in_file(charge_str, header_index, footer_index, data, file_path)
        c_end_idx = find_string_in_file(multiplicity_str, c_start_idx, footer_index, data, file_path)
        # get number
        charge_number = int(data[c_start_idx:c_end_idx].replace('Charge', '').replace('=', '').strip())

        # Multiplicity
        m_start_idx = find_string_in_file(multiplicity_str, header_index, footer_index, data, file_path)
        m_end_idx = find_string_in_file(new_line_string, m_start_idx, footer_index, data, file_path)
        # get number
        multiplicity_number = int(data[m_start_idx:m_end_idx].replace('Multiplicity', '').replace('=', '').strip())

        # add this charge,multiplicity info to master
        with open('master_values.json') as fp:
            json_contents = json.load(fp)
            json_contents["charge_multiplicity"] = f' {charge_number}  {multiplicity_number}'

        with open('master_values.json', 'w') as fp:
            json.dump(json_contents, fp)

        # check if optimization was completed, then find optimized geometry data
        opt_completed_index = find_string_in_file(opt_completed_str, 0, len(data), data, file_path)
        input_orientation_index = find_string_in_file(input_orientation_str, opt_completed_index, len(data), data, file_path)
        footer_index = find_string_in_file(standard_orientation_str, opt_completed_index, len(data), data, file_path)

        # proceed as normal, everything appears to be in order
        log.info(
            f"Gaussian log file does contain ({input_orientation_str=}) and ({opt_completed_str=})"
            " so it appears to be correctly formatted."
        )
        log.debug(f"Found ({input_orientation_str=}) and ({opt_completed_str=}) in {file_path}")

        # header_index-1 to find it properly
        start_idx = find_string_in_file(input_orientation_str, input_orientation_index-1, footer_index, data, file_path)
        end_idx = find_string_in_file(distance_matrix_str, start_idx, footer_index, data, file_path)

        line_list = data[start_idx:end_idx].splitlines()
        # after input orientation line, go past 4 more lines, but 2 lines above distance matrix
        # to get into the part of the matrix we want
        opt_geom_str = []
        for line in line_list[5:-2]:
            opt_geom_str.append(line[17:])     # get past center number
        opt_geom_str = "\n".join(opt_geom_str)

        output = "".join([
            f"%chk={molecule_name}_excited.chk\n",
            f"#p td=(nstates={nof_states}) b3lyp/6-311++g(d,p) scrf=check\n",
            "\n",
            f"Molecule {nof_states} excited states 6-311++g(d,p)\n",
            "\n",
            f"{charge_number} {multiplicity_number}\n",
            f"{opt_geom_str}\n",
            "\n"
        ])

        with open(output_file_path, 'w') as new_fp:
            new_fp.write(output)
        log.debug(f"Generated the excited state input file:\n{output_file_path}\n")
        return

    except Exception as e:
        log.warning("Cannot generate excited state input file.")
        log.warning(f"If {footer_str} not found, add 'p' after # in {molecule_name}_optfreq.com and rerun step1.py")
        raise e

    return


def perform_step_1b(molecule_name, nproc, mem, automate_flag=False):
    log.info("Beginning STEP 1b of the procedure.")
    create_excited_input(f'{molecule_name}_optfreq.log', N_atoms, N_exstates)

    # change mem/proc settings if needed
    set_mem_and_processor(f'{molecule_name}_excited.com', 'excited', nproc, mem)
    set_mem_and_processor('submit', 'excited', nproc, mem)

    command = ' '.join([
        'sbatch',
        f'--job-name={molecule_name}_excited',
        "--output='slurm-%j.out'",
        'submit'
    ])

    if automate_flag:
        return command
    else:
        os.system(command)
        return


if (__name__ == '__main__'):

    # code to run when file called directly
    with open('master_values.json') as fp:
        value = json.load(fp)

    molecule_name = value["molecule_name"]
    nproc = value["step1b"]["nproc"]
    mem = value["step1b"]["mem"]

    perform_step_1b(molecule_name, nproc, mem)
