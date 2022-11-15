##################################################################################
'''

Step one of the vibronic model process.

In this step we modify the Gaussian input file and submit script to utilize
appropriate computational resources. Then, we submit the calculation to get
the output for optimized ground-state geometry and frequencies of the molecule.

It requires this information:
- Gaussian input .com file
- submit script

It produces these files:
- Gaussian output .log file

'''
###################################################################################

# system imports
import os
import json

# third party imports

# local imports
from log_conf import log
import helper

with open('master_values.json') as fp:
        value = json.load(fp)

# global molecule_name
molecule_name = value["molecule_name"]

# split this into two functions
def set_mem_and_processor(file_path, calc_type, nproc, mem):
    '''
    Sets the amount of memory and number of processors
    for the job in both input file and submit script.

    file_path is string name of file.
    e.g. 'submit'

    calc_type is string of calculation to be performed.
    e.g. 'optfreq', 'excited', 'displacement'
    '''

    log.info(f"User is asking for: {nproc} processors, {mem}GB of memory.")

    if not os.path.isfile(file_path):
        # if file_path is not a file
        log.debug(f"The provided path is not valid:\n{file_path}")

    # attempt to read the file
    data = helper.readfile(file_path, name="submit script")

    # first try/except block figures out if SBATCH script or not
    try:
        target_string = '#SBATCH'

        if target_string in data:
            # if True, file is submit script
            log.info("File is a SBATCH formatted")
            log.debug(f"Found {target_string} in {file_path}")

            # update file_path with user provided nproc and mem
            first_line = '#!/bin/bash'
            mem_str = f'#SBATCH --mem={mem}G'
            proc_str = f'#SBATCH --cpus-per-task={nproc}'
            load_gauss = 'module load gaussian/16'

            input_file = f"{molecule_name}_{calc_type}.com"
            output_file = f"{molecule_name}_{calc_type}.log"

            if (calc_type == 'gaussian'):
                input_file = "gaussian.com"
                output_file = "gaussian.log"

            submit_cmd = f'g16 < {input_file} >& {output_file}'

            output = '\n'.join([first_line, mem_str, proc_str, load_gauss, submit_cmd])

            with open(file_path, 'w') as new_fp:
                new_fp.write(output)

            log.info(f"Updated file: {file_path}\n with appropriate parameters")

            return

    except Exception as e:
        log.warning(
                f"File {file_path}\ndoes not have ({target_string}) inside it.\n"
                "Please check submit file, something went wrong."
            )
        log.warning("Updating Gaussian memory and processor settings in SBATCH file unsuccessful.")
        raise e

    try:
        with open(file_path, 'r') as fp:
            data = fp.readlines()

        # target_freq_str = 'freq=hpmodes'
        nproc_str = '%nproc='
        nproc_str2 = '%nprocshared='
        mem_str = '%Mem='

        # if any of mem_str or nproc_str/2 is in
        # any of the lines of strings
        # for every single line
        if any(mem_str in s for s in data) or \
           any(nproc_str in s for s in data) or \
           any(nproc_str2 in s for s in data):

            # perform replacements
            for line in data:
                if line.find(mem_str) != -1:
                    location = data.index(line)
                    data[location] = f'%Mem={mem}GB\n'
                    log.info(f"Updated file: {file_path}\n memory settings.")

            for line in data:
                location = data.index(line)

                if line.find(nproc_str) != -1:
                    data[location] = f'%nproc={nproc}\n'

                elif line.find(nproc_str2) != -1:
                    data[location] = f'%nprocshared={nproc}\n'

            output = "".join(data)

            with open(file_path, 'w') as fp:
                fp.write(output)
            log.info(f"Updated file: {file_path}\n processor settings.")

        else:  # add mem/proc settings starting from second line
            with open(file_path, 'r') as fp:
                data = fp.read()

            # find first newline character
            insert_here = data.index('\n')

            # grab the text above and below
            # the part we are replacing
            above = data[:insert_here+1]
            below = data[insert_here+1:]

            # the updated string
            replaced_part = (
                f'%nprocshared={nproc}' + '\n'
                f'%Mem={mem}' + 'GB' + '\n'
            )

            # glue them all together
            new_file_contents = above + replaced_part + below

            with open(file_path, 'w') as fp:
                fp.write(new_file_contents)
            log.info(f"Updated file: {file_path}\n processor and memory settings.")

            return

    except Exception as e:
        log.warning("Updating Gaussian memory and processor settings in Gaussian input file unsuccessful.")
        raise e

    return


def adjust_chk(file_path):
    ''' Sometimes you go from Gaussview the %chk linking
    to your N drive might cause an issue. Adjust chk output.
    '''
    with open(file_path, 'r') as fp:
        data = fp.readlines()

    if data[0].find("%chk") != -1:
        log.info(f"{file_path} file contains .chk point output line, updated {file_path} to remove it.")
        data = data[1:]
        out = "".join(data)

        with open(file_path, 'w') as fp:
            fp.write(out)

    return


def perform_step_1(molecule_name, nproc, mem, automate_flag=False):
    # code to run when file called directly
    log.info("Beginning STEP 1 of the procedure.")

    # apply the memory and processor changes to both files
    set_mem_and_processor(f'{molecule_name}_optfreq.com', 'optfreq', nproc, mem)
    set_mem_and_processor('submit', 'optfreq', nproc, mem)

    adjust_chk(f'{molecule_name}_optfreq.com')

    command = (
        "sbatch"
        f" --job-name={molecule_name}_opt_freq"
        " --output='slurm-%j.out'"
        " submit"
    )

    if automate_flag:
        return command
    else:
        os.system(command)
        return

# this conditional is useful if imported to other script and want to run selective parts of code
if (__name__ == '__main__'):

    with open('master_values.json') as fp:
        value = json.load(fp)

    # global molecule_name
    molecule_name = value["molecule_name"]
    nproc = value["step1"]["nproc"]
    mem = value["step1"]["mem"]

    perform_step_1(molecule_name, nproc, mem)
