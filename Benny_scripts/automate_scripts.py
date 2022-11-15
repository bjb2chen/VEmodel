##################################################################################
'''

Automating script execution.

We are attempting in this script to execute the step scripts 1 through 6
automatically. This is beneficial because we can then proceed to run the next step
immediately after one step terminates, without user interaction.

See issue and snippet in gitlab.

'''
###################################################################################

# system imports
# import os
import time
import sys
# from subprocess import CalledProcessError
import json

# third party imports

# local imports
from log_conf import log
import server


def wait_for_server():
    '''Instructs server to wait for given duration'''
    time.sleep(1)
    return


def run_step_one_a(dictionary):
    '''
    Automates step1.py through submitting the job to the server,
    using memory/processor settings values as dictated in step0.py.

    Optimizes the molecule to g.s. geometry and obtains frequencies.
    '''

    # open in append mode
    # with open('automation_output.txt', 'a') as fp:

    import step1
    command = step1.perform_step_1(
        dictionary['molecule_name'],
        dictionary['step1']['nproc'],
        dictionary['step1']['mem'],
        automate_flag=True
    )

    id_job, out, err = server.submit_job(command)
    log.info(f'Submitted job {id_job}\nStdout:{out}\nStderr{err}\n')

    log.info("Waiting 1 second")
    wait_for_server()

    return id_job


def run_step_one_b(dictionary, previous_job_id=None):
    '''
    Automates step1b.py, utilizes dictionary's memory/processor settings values.

    TD-DFT calculation to get excited state energies from g.s. geometry.
    '''

    #  not sure if this part is true anymore
    #  if job completes too quickly then dependency wont fire
    #  so each job (the .sh files) should sleep ~15-20 seconds after its done
    #  use the automate flag to handle that

    import step1b
    command = step1b.perform_step_1b(
        dictionary['molecule_name'],
        dictionary['step1b']['nproc'],
        dictionary['step1b']['mem'],
        automate_flag=True
    )

    # modify the command to be dependent on the previous step successfully completing
    if previous_job_id is not None:
        command = f"sbatch --dependency=afterok:{previous_job_id} " + ' '.join(command.split(' ')[1:])

    id_job, out, err = server.submit_job(command)
    log.info(f'Submitted job {id_job}\nStdout:{out}\nStderr{err}\n')

    return id_job


def run_step_two(dictionary, previous_job_id=None):
    '''
    Performs step2.py.

    Operation is identical for automated/manual execution.

    Does not need submission to slurm server, as step2 is just exctracting
    text and generating things like atomic masses data file. All facilitated by scripts.
    '''

    import step2
    command = step2.perform_step_2(automate_flag=True)

    # if we are automating all the steps
    if previous_job_id is not None:

        # make a new bash file to submit the job
        temp_step_two_automated_bash_path = "step_two_bash.sh"

        # write the bash file
        with open(temp_step_two_automated_bash_path, 'w') as fp:
            fp.write(f"#!/bin/bash\n{command}\n")

        command = f"sbatch --dependency=afterok:{previous_job_id} {temp_step_two_automated_bash_path}"

        id_job, out, err = server.submit_job(command)
        log.info(f'Submitted job {id_job}\nStdout:{out}\nStderr{err}\n')
        return id_job

    # if we are running each step independently
    else:
        out, err = server.run_locally(command)
        log.info(f"Ran step two\nStdout:{out}\nStderr{err}\n")
        # maybe check some other things?
        log.info("Generated the state_s0 and masses file")
        return


def run_step_three(dictionary, previous_job_id=None):
    '''
    Performs step3.py.

    Operation is identical for automated/manual execution.

    Does not need submission to slurm server, as step3 is just setting up
    (3N-6 * N - 1) number of Gaussian calculations. N is number of displacements
    as dictated by master_values.json setup in step0.py.
    '''

    import step3
    command = step3.perform_step_3(automate_flag=True)

    # if we are automating all the steps
    if previous_job_id is not None:

        # make a new bash file to submit the job
        temp_step_three_automated_bash_path = "step_three_bash.sh"

        # write the bash file
        with open(temp_step_three_automated_bash_path, 'w') as fp:
            fp.write(f"#!/bin/bash\n{command}\n")

        command = f"sbatch --dependency=afterok:{previous_job_id} {temp_step_three_automated_bash_path}"

        id_job, out, err = server.submit_job(command)
        log.info(f'Submitted job {id_job}\nStdout:{out}\nStderr{err}\n')
        return id_job

    # if we are running each step independently
    else:
        out, err = server.run_locally(command)
        log.info(f"Ran step two\nStdout:{out}\nStderr{err}\n")

        log.info(
                "Generated the output_genepointder file\n"
                "Generated the gaussian.com file."
            )

        return

    return id_job


def run_step_four(dictionary, previous_job_id=None):
    '''
    Automates step4.py, utilizes dictionary's memory/processor settings values.

    Runs all the (3N-6 * N - 1) Gaussian calculations prepared by step3
    according to memory/processor settings in master_values.json.

    This will give the many adiabatic states.
    '''

    import step4
    command = step4.perform_step_4(
        dictionary['molecule_name'],
        dictionary['step4']['nproc'],
        dictionary['step4']['mem'],
        automate_flag=True
    )

    # modify the command to be dependent on the previous step successfully completing
    if previous_job_id is not None:
        command = f"sbatch --dependency=afterok:{previous_job_id} " + ' '.join(command.split(' ')[1:])

    id_job, out, err = server.submit_job(command)
    log.info(f'Submitted job {id_job}\nStdout:{out}\nStderr{err}\n')

    return id_job


def run_step_five(dictionary, previous_job_id=None):
    '''
    Automates step5.py, utilizes dictionary's memory/processor settings values.

    Runs the diabatization scheme to maximize overlap between adiabatic states
    to achieve diabatic states, linear coupling terms.
    '''

    import step5
    command = step5.perform_step_5(
        dictionary['molecule_name'],
        dictionary['step5']['nproc'],
        dictionary['step5']['mem'],
        automate_flag=True
    )

    # modify the command to be dependent on the previous step successfully completing
    if previous_job_id is not None:
        command = f"sbatch --dependency=afterok:{previous_job_id} " + ' '.join(command.split(' ')[1:])

    id_job, out, err = server.submit_job(command)
    log.info(f'Submitted job {id_job}\nStdout:{out}\nStderr{err}\n')

    return id_job


def run_step_six(dictionary, previous_job_id=None):
    '''
    Performs step6.py.

    Operation is identical for automated/manual execution.

    Similar to step2/3, simply runs script to generate operator .op file
    for MCTDH input.
    '''

    import step6
    command = step6.perform_step_6(automate_flag=True)

    molecule_name = dictionary["molecule_name"]
    n_exstates = dictionary["num_excited_states"]
    n_freqs = dictionary["num_freqs"]

    # if we are automating all the steps
    if previous_job_id is not None:

        # make a new bash file to submit the job
        temp_step_six_automated_bash_path = "step_six_bash.sh"

        # write the bash file
        with open(temp_step_six_automated_bash_path, 'w') as fp:
            fp.write(f"#!/bin/bash\n{command}\n")

        command = f"sbatch --dependency=afterok:{previous_job_id} {temp_step_six_automated_bash_path}"

        id_job, out, err = server.submit_job(command)
        log.info(f'Submitted job {id_job}\nStdout:{out}\nStderr{err}\n')
        return id_job

    # if we are running each step independently
    else:
        out, err = server.run_locally(command)
        log.info(f"Ran step six\nStdout:{out}\nStderr{err}\n")

        log.info(
            "Please proceed and check if files "
            f"{molecule_name}_coupling.out and op_{molecule_name}{n_freqs}Q_{n_exstates}st.op "
            "have correct output."
        )

        return


def mini_example():
    ''' for testing or documentation purposes '''
    return {
            "molecule_name": "water",
            "num_atoms": 3,
            "num_excited_states": 2,
            "num_freqs": 3,
            "num_displacements": 4,
            "scaling_factor": "1.d0",
            "sorting_atoms": "no",
            "step1": {
                "nproc": 8,
                "mem": 30
            },
            "step1b": {
                "nproc": 7,
                "mem": 19
            },
            "step4": {
                "nproc": 10,
                "mem": 19
            },
            "step5": {
                "nproc": 11,
                "mem": 19
            }
        }


def automate_scripts():

    # have file paths for working directory for master_values.json
    # so like root_dir = /work/benny/molecule_project
    # we agree the full path looks like this
    # f"/work/benny/molecule_project/{molecule_name}/master_values.json"
    # f"/work/benny/molecule_project/{molecule_name}/scripts/"
    # f"/work/benny/molecule_project/{molecule_name}/output/step7/"
    # f"/work/benny/molecule_project/{molecule_name}/output/step1/"
    # f"/work/benny/molecule_project/{molecule_name}/logs/"
    # attempt to read from master_values.json
    try:
        # work_path = f"/work/benny/molecule_project/{molecule_name}/"
        # master_file_path = work_path + 'master_values.json'
        master_file_path = 'master_values.json'
        with open(master_file_path, 'r') as fp:
            master_values = json.load(fp)

    except Exception as e:
        log.debug(f"Couldn't read the master file?\n{master_file_path}\n")
        raise e

    # then pass in the dictionary
    dictionary = master_values
    # dictionary = mini_example()

    if len(sys.argv) >= 2:
        # we are assuming we passed in a number
        step = int(sys.argv[1])
        if step == 1:
            run_step_one_a(dictionary)
        elif step == 2:
            run_step_one_b(dictionary)
        elif step == 3:
            run_step_two(dictionary)
        elif step == 4:
            run_step_three(dictionary)
        elif step == 5:
            run_step_four(dictionary)
        elif step == 6:
            run_step_five(dictionary)
        elif step == 7:
            run_step_six(dictionary)
        elif step == -1:
            # now we are testing full automation
            id_job1 = run_step_one_a(dictionary)
            id_job2 = run_step_one_b(dictionary, id_job1)
            id_job3 = run_step_two(dictionary, id_job2)
            id_job4 = run_step_three(dictionary, id_job3)
            id_job5 = run_step_four(dictionary, id_job4)
            id_job6 = run_step_five(dictionary, id_job5)
            id_job7 = run_step_six(dictionary, id_job6)
        else:
            raise Exception("wrong number")
    else:
        print("Please provide a step number argument vector.")


if (__name__ == '__main__'):
    automate_scripts()
