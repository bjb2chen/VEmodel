###################################################################################
'''

Step three of the vibronic model process.

In this step we generate the input_genepointder file, which will create all the
normal-mode displacement geometries (from original geometry) calculations in the
form of a single gaussian input file.
Remember that the objective here is to ultimately get the diabatic states from
adiabatic. The diabatization scheme needs to know the overlap between adiabatic
states and these geometries.

It requires this information:
- Number of atoms
- Number of frequencies
- Masses of each atom (output from step2)
- state_s0 file       (output from step2)
- Number of displacements along each mode
- Dimensionless shift

It produces these files:
- input_genepointder
- output_genepointder
- gaussian.com
'''
###################################################################################

# system imports
import os
# import os.path
# from os.path import join
import json

# third party imports

# local imports
import helper
from helper import find_string_in_file, readfile
from log_conf import log

with open('master_values.json') as fp:
    value = json.load(fp)

extracted_masses_file_path = "input_genepointder"
gene_pointder_path = "/home/bjb2chen/CHEM494/CODES/gene-pointder.e"
new_line_string = '\n'


def extract_masses(file_path):
    ''' takes each line of masses file '''
    data = readfile(file_path, name="masses output")

    # want every single line in masses file
    masses_list = data.splitlines()
    log.info(f"Masses have been extracted:\n{masses_list}")

    # first line in input_genepointder relates to number of atoms
    line1 = value["num_atoms"]

    # second line is for number of normal modes
    line2 = value["num_freqs"]

    # next section is for the atomic masses
    masses_lines = new_line_string.join(str(line) for line in masses_list)

    # third last line just is string 'state_s0'
    terminal_line1 = "\'state_s0\'"

    # ask for number of displacements
    terminal_line2 = value["num_displacements"]

    # last line is dimensionless shift
    terminal_line3 = '0.1d0'

    # put all these lines together and write them out to input_genepointder
    output = [
        f"{str(line1):<30}! number of atoms",
        f"{str(line2):<30}! number normal modes",
        masses_lines,
        f"{str(terminal_line1):<30}! state file prepared by prep.fc.e",
        f"{str(terminal_line2):<30}! number of displacements along each mode",
        f"{str(terminal_line3):<30}! dimensionless shift",
    ]

    # join the output together
    makestr = new_line_string.join(str(line) for line in output)
    log.info("Collected all necessary input info.")

    try:
        with open(extracted_masses_file_path, 'w') as fp:
            fp.write(makestr)
        log.debug(f"Generated the input_genepointder file:\n{extracted_masses_file_path}\n")

    except Exception as e:
        log.warning("Cannot generate input_genepointder file.")
        raise e

    return


def perform_step_3(automate_flag=False):
    log.info("Beginning STEP 3 of the procedure.")
    extract_masses("masses")

    input_file = extracted_masses_file_path
    output_file = "output_genepointder"
    command = (f"{gene_pointder_path} < {input_file} > {output_file}")

    if automate_flag:
        # step 3 doesn't take any user input
        # it doesn't depend/interact with slurm/sbatch
        # the operation is identical for automated/manual operation
        return command

    else:
        os.system(command)
        log.info(
            "Generated the output_genepointder file\n"
            "Generated the gaussian.com file."
        )
        return


if (__name__ == '__main__'):

    # code to run when file called directly
    with open('master_values.json') as fp:
        value = json.load(fp)

    molecule_name = value["molecule_name"]

    perform_step_3()
