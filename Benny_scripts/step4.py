###################################################################################
'''

Step four of the vibronic model process.

From the gaussian.com file, the default 6N-11 displaced geometries all have
default settings of cam-b3lyp/6-31G(d), excited states=8, nproc=1, Mem=2GW.
Reconfigure these settings to appropriate ones.

This step is expected to take quite some time, dependent how fast the gaussian
output can be produced. Performs TD-DFT on each adiabatic state.

It requires this information:
- gaussian.com (from running input_genepointder)
- formchk4script

It produces these files:
- 6N-11 .chk files
- 6N-11 .fchk files

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

# what to look for
header_str = "Gaussian, Inc."
footer_str = "Normal termination of Gaussian"
link_str = "--Link1--"
theory_basis_str = "cam-b3lyp/6-31G(d)"
exstates_str = "nstates=8,conver=6"
theory_basis_nstates_str = "#p cam-b3lyp/6-31G(d)  nosymm iop(3/33=4) td(nstates=8,conver=6) iop(9/40=5)"
charge_mult_str = "0  1"
output_file_path = "gaussian.com"

# desired constants
molecule_name = value["molecule_name"]
theory_basis_val = "b3lyp/6-311++G(d,p)"
exstates_val = value["num_excited_states"]
charge_mult_val = value["charge_multiplicity"]


def configure_theorybasis_nstates(file_path):

    try:
        with open(file_path, 'r') as fp:
            data = fp.readlines()

        for line in data:
            if line.find(link_str) != -1:
                log.debug("File should be default gaussian.com configuration, as anticipated.")
                break
            else:
                log.debug(
                    "File is either not gaussian.com or improperly formatted\n"
                    f"Could not find {link_str}."
                )
                break

        # Similar to new_step1, replace all instances
        # of wrong level of theory, basis set, etcetera settings in order to have unity between calculations
        if any(theory_basis_str in s for s in data) or \
           any(exstates_str in s for s in data) or \
           any(charge_mult_str in s for s in data):

            # perform replacements
            for line in data:
                if line.find(theory_basis_nstates_str) != -1:
                    location = data.index(line)
                    data[location] = (f'#p {theory_basis_val}  nosymm iop(3/33=4) td(nstates={exstates_val}) iop(9/40=5)' + '\n')
                    log.debug(f"Updated file: {file_path}\n \
                               level of theory, basis set, and number of excited states was corrected.")

            for line in data:
                if line.find(charge_mult_str) != -1:
                    location = data.index(line)
                    data[location] = (f"{charge_mult_val}" + '\n')
                    log.debug(f"Updated file: {file_path}\n charge and multiplicity values was corrected.")

            output = "".join(data)

            with open(file_path, 'w') as fp:
                fp.write(output)
                log.info(f"Updated file: {file_path}\n all settings configured.")

            return

    except Exception as e:
        log.warning(f"Cannot configure {file_path} file.")
        raise e

    return


def perform_step_4(molecule_name, nproc, mem, automate_flag=False):
    log.info("Beginning STEP 4 of the procedure.")
    configure_theorybasis_nstates('gaussian.com')
    # configure memory and nproc again, then submit calculation
    set_mem_and_processor('gaussian.com', 'gaussian', nproc, mem)
    set_mem_and_processor('submit', 'gaussian', nproc, mem)

    command = (
        f"sbatch"
        f" --job-name={molecule_name}_gauss"
        " --output='slurm-%j.out'"
        " submit")

    if automate_flag:
        return command
    else:
        os.system(command)
        return


if (__name__ == '__main__'):

    with open('master_values.json') as fp:
        value = json.load(fp)

    # global molecule_name
    molecule_name = value["molecule_name"]
    nproc = value["step4"]["nproc"]
    mem = value["step4"]["mem"]

    perform_step_4(molecule_name, nproc, mem)
