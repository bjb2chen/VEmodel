###################################################################################
'''

Step five of the vibronic model process.

This step is the fundamental step and where the wizardry happens. From my
understanding, the adiabatic states get rotated by the unitary operator
to acheive the diabatic states (which are linear combos of adiabatic).
All of this falls under the diabatization scheme to get coupling info.

Runtime can potentially be very long.

It requires this information:
- gaussian.log output from prev step
- Number of normal modes
- Number of displacements/steps per mode (likely to be 2)
- Dimensionless shift
- Number of electronic states

It produces these files:
- overdia.out

'''
###################################################################################

# system imports
import os
import json

# third party imports

# local imports
from log_conf import log

with open('master_values.json') as fp:
    value = json.load(fp)

molecule_name = value["molecule_name"]
natoms = value["num_atoms"]
ndisplacements = value["num_displacements"]
step_4_output_file = "gaussian.log"
succesful_termination_string = "Normal termination of Gaussian"
input_filename = f'{molecule_name}_overdia.inp'
output_filename = f'{molecule_name}_overdia.out'
bash_filename = f'{molecule_name}_overdia.sh'
overdia_par_path = "/home/bjb2chen/CHEM494/CODES/./overdia-par.e"

formcheck_path = "./formchk4script"


def confirm_step_4_completed_succesfully(file_path=step_4_output_file):
    ''' Make sure (to the best of our knowledge) the last step completed successfully.

        For 2 displacements, it is 2*(3N-6)+1 = 6N-12+1
        so for water at 2 displacements, 6*3-11 = 7 calculations.
        Therefore, look for 7 'Normal termination of Gaussian' instances
        to know if last step completed successfully.
    '''
    if not os.path.isfile(file_path):
        log.debug(f"The provided path is not valid:\n{file_path}\n")

    try:
        with open(file_path, 'r') as fp:
            data = fp.read()

    except OSError as e:
        log.debug(f"Couldn't open step 4 output file?\n{file_path}\n")
        raise e

    log.info(
        f"Anticipated {ndisplacements}*(3*({natoms=})-6)+1 = "
        f"{((int(ndisplacements)*int(3*int(natoms)-6)+1))} instances "
        f"of \'{succesful_termination_string}\' in {step_4_output_file}.")

    log.info(
        f"Found {data.count(succesful_termination_string)} instances "
        f"of \'{succesful_termination_string}\' in {step_4_output_file}.")

    if (succesful_termination_string not in data) or \
       (data.count(succesful_termination_string) != (int(ndisplacements)*int(3*int(natoms)-6)+1)):
        raise Exception(
            "It appears that step 4 did not complete successfully"
            f"Please check the output file: \n{file_path}\n"
        )
    return


def create_overdia_inp():
    '''
    Creates the overdia_inp file, which is 4 lines long.
    '''
    # way 1
    arr = [
        ("\'" + step_4_output_file + "\'"),
        str(value['num_freqs']) + ' ' + str(value['num_displacements']),
        '0.1d0',
        value['num_excited_states']
    ]

    output = '\n'.join([
        f"{arr[0]:<30}! multilink gaussian log file",
        f"{arr[1]:<30s}! number of modes and number of steps per mode",
        f"{arr[2]:<30s}! dimensionless shift",
        f"{arr[3]:<30d}! number of electronic states",
    ])

    try:
        with open(input_filename, 'w') as fp:
            fp.write(output)

        log.debug(f"Generated the {input_filename} file")

    except Exception as e:
        log.warning(f"Cannot generate {input_filename} file.")
        raise e

    return


# need special parallelized submit script, can't use step1...
def overdia_sh(nproc, mem):
    '''
    Create the script for diabatization memory/processor settings.
    '''
    #nproc = input("Amount of processors needed, max 12: ")
    #mem = input("Amount of memory (in GB) needed, max 48GB: ")
    #log.info(f"User is asking for: {nproc} processors, {mem}GB of memory.")

    first_line = '#!/bin/bash'
    mem_str = f'#SBATCH --mem={mem}G'
    proc_str = f'#SBATCH --cpus-per-task={nproc}'
    load_gauss = 'module load gaussian/16'
    chmod_str = f'chmod +x {formcheck_path}'
    formchk = f'{formcheck_path}'
    parallel = 'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'

    slurm1 = 'echo "SLURM_SUBMIT_DIR:$SLURM_SUBMIT_DIR"'
    slurm2 = 'echo "SLURM_CPUS_PER_TASK:$SLURM_CPUS_PER_TASK"'
    slurm3 = 'echo "SLURM_MEM_PER_CPU:$SLURM_MEM_PER_CPU"'
    slurm4 = 'echo "SLURM_MEM_PER_NODE:$SLURM_MEM_PER_NODE"'

    cmd_line = f'{overdia_par_path} < {input_filename} > {output_filename}'

    output = '\n'.join([
        first_line, mem_str, proc_str, load_gauss, chmod_str,
        formchk, parallel, slurm1, slurm2, slurm3, slurm4, cmd_line
    ])

    with open(bash_filename, 'w') as new_fp:
        new_fp.write(output)

    log.debug(f"Generated the {bash_filename} file")

    return


def perform_step_5(molecule_name, nproc, mem, automate_flag=False):
    """ open and write to empty file the create_overdia_inp """
    log.info("Beginning STEP 5 of the procedure.")
    confirm_step_4_completed_succesfully()
    log.info("Step 4 completed successfully.")
    create_overdia_inp()

    if False:
        # run formchk4script, have to wait till gaussian calculation terminates
        os.system('module load gaussian/16')      # have to first load in gaussian
        os.system(f'chmod +x {formcheck_path}')    # ensures the script is executable
        os.system(f'{formcheck_path}')

    overdia_sh(nproc, mem)

    command = " ".join([
        'sbatch',
        f"--job-name={molecule_name}_overdia",
        "--output='slurm-%j.out'",
        f"{bash_filename}",  # the overdia submit script
    ])

    if automate_flag:
        return command
    else:
        os.system(command)
        return



if (__name__ == '__main__'):

    with open('master_values.json') as fp:
        value = json.load(fp)

    molecule_name = value["molecule_name"]
    nproc = value["step5"]["nproc"]
    mem = value["step5"]["mem"]
    
    perform_step_5(molecule_name, nproc, mem)
