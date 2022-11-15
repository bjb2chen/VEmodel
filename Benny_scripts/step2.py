###################################################################################
'''

Step two of the vibronic model process.

In this step we generate the input_prepfc file by taking the input orientations
and harmonic frequencies from the OPT+FREQ .log output.

It requires this information:
- Number of atoms
- Selected frequency scaling factor
- Choice of sorting atoms
- Input orientation matrix from .log output
- Harmonic frequencies from .log output

It produces these files:
- input_prepfc
- masses
- state_s0 (adiabatic states)

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

# for first three lines of input_prepfc
line1 = value["num_atoms"]
line2 = value["scaling_factor"]
line3 = value["sorting_atoms"]

molecule_name = value["molecule_name"]
output_file_path = "input_prepfc"

header_str = "Gaussian, Inc."
footer_str = "Normal termination of Gaussian"
input_orientation_str = "Input orientation:"
distance_matrix_str = "Distance matrix"
standard_orientation_str = "Standard orientation"
harm_freq_str = "Harmonic frequencies"

new_line_string = '\n'

log.info(f"User has specified {line1} atoms in {molecule_name} molecule.")
log.info(f"User has specified {line2} scaling factor for {molecule_name} molecule.")
log.info(f"User has specified \"{line3}\" for choice of sorting atoms in {molecule_name} molecule.")


def extract_orientation_harmonics(file_path):
    '''
    reads in a .log output file from a GAUSSIAN OPT+FREQ calculation,
    and returns nothing, but has effect of writing the input orientation and
    harmonic frequencies sections from .log to the input_prepfc file.
    '''

    # still need to do this
    output = [str(line1), str(line2), str("\'" + line3 + "\'")]

    data = readfile(file_path, name=".log output")

    try:
        header_index = find_string_in_file(header_str, 0, len(data), data, file_path)
        footer_index = find_string_in_file(footer_str, 0, len(data), data, file_path)
        input_orientation_index = find_string_in_file(input_orientation_str, 0, len(data), data, file_path)

        # proceed as normal, everything appears to be in order
        log.info(
            f"Gaussian log file does contain ({input_orientation_str=})"
            " so it appears to be correctly formatted."
        )
        log.debug(f"Found ({input_orientation_str=}) in {file_path}")

        # grab header and footer
        header_index = data.index(input_orientation_str)
        footer_index = data.index(standard_orientation_str)
        log.debug(f"Found ({input_orientation_str=}) in {header_index}")

        # find the input orientation matrix
        input_orientation_index = data.index(input_orientation_str)

        # header_index-1 to find it properly
        start_idx = find_string_in_file(input_orientation_str, header_index-1, footer_index, data, file_path)
        end_idx = find_string_in_file(distance_matrix_str, start_idx, footer_index, data, file_path)

        line_list = data[start_idx:end_idx].splitlines()
        # after input orientation line, go past 4 more lines
        # to get into the part of the matrix we want
        for line in line_list[5:-2]:
            output.append(line)

        # new header and footer
        # NOTE THAT THERE WILL BE TWO OCCURENCES OF "Normal termination of Gaussian"
        # one for opt, one for freq. new_footer_index needs the second, later occurence.
        new_header_index = data.index(harm_freq_str)
        log.debug(f"Found ({harm_freq_str=}) in {new_header_index}")
        new_footer_index = data.find(footer_str, data.find(footer_str)+1)
        log.debug(f"Found ({footer_str=}) in {new_footer_index}")

        # find first instance of Harmonic frequencies string
        # take every line 4 below it (inclusive) to end of .log file
        h_start_idx = find_string_in_file(harm_freq_str, new_header_index-1, new_footer_index, data, file_path)
        h_line_list = data[h_start_idx:].splitlines()
        for line in h_line_list[4:]:
            output.append(line)

        # join the output together
        makestr = new_line_string.join(str(line) for line in output)

        # input_prepfc here could be a file with first three lines already have content
        with open(output_file_path, 'w') as new_fp:
            new_fp.write(makestr)
        log.info(f"Generated the input_prepfc file:\n{output_file_path}\n")
        return

    except Exception as e:
        log.warning("Cannot generate input_prepfc file.")
        raise e

    return


def perform_step_2(automate_flag=False):
    """ x """
    log.info("Beginning STEP 2 of the procedure.")

    # Create the input_prepfc file
    extract_orientation_harmonics(f'{molecule_name}_optfreq.log')
    # command = f"{root_bin}{executeable_file_name} < {input_file_name} > {output_file_name}"
    command = '/home/bjb2chen/CHEM494/CODES/prep.fc.e < input_prepfc > state_s0'

    if automate_flag:
        # step 3 doesn't take any user input
        # it doesn't depend/interact with slurm/sbatch
        # the operation is identical for automated/manual operation
        return command

    else:
        os.system(command)
        log.debug("Generated the state_s0 and masses file")
        return


if (__name__ == '__main__'):

    # code to run when file called directly
    with open('master_values.json') as fp:
        value = json.load(fp)

    molecule_name = value["molecule_name"]

    perform_step_2()
