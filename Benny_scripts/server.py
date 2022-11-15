# system imports
import subprocess
import threading
import signal
import sys
import os
import time

# third party imports

# local imports
from log_conf import log


# lock for asynchronous communication
job_state_lock = threading.Lock()
job_almost_done_flag = False


def subprocess_submit_asynch_wrapper(cmd, **kwargs):
    """ wrapper for subprocess.Popen function to allow for different implementation for different python versions"""
    if sys.version_info[:2] >= (3, 7):
        return subprocess.Popen(cmd, text=True, **kwargs)
    if (3, 5) <= sys.version_info[:2] <= (3, 7):
        return subprocess.Popen(cmd, universal_newlines=True, **kwargs)


def subprocess_run_wrapper(cmd, **kwargs):
    """ wrapper for subprocess.run function to allow for different implementation for different python versions"""
    if sys.version_info[:2] >= (3, 7):
        return subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if (3, 5) <= sys.version_info[:2] <= (3, 7):
        return subprocess.run(cmd, universal_newlines=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              **kwargs)


def get_hostname():
    """returns the hostname of the cluster of the server (from SLURM) as a string"""
    cmd = ['scontrol', 'show', 'config']
    result = subprocess_run_wrapper(cmd)

    for line in result.stdout.splitlines():
        if "ClusterName" in line:
            return line.split('=')[1].strip()
    else:
        raise Exception("Did not find ClusterName in the config data from scontrol!?")


def check_acct_state(id_job):
    """returns the recorded state of the job (from SLURM) as a string"""
    cmd = ['sacct', '-n', '-o', 'state', '-j', str(id_job)]
    result = subprocess_run_wrapper(cmd)
    return result.stdout


def check_running_state(id_job):
    """returns the running state of the job (from SLURM) as a string"""
    cmd = ['scontrol', 'show', 'job', str(id_job)]
    result = subprocess_run_wrapper(cmd)
    return result.stdout, result.stderr


def synchronize_with_job(id_job, job_type="default"):
    """synchronizes with a submitted job """
    log.lock(f"Synchronizing with job (id={id_job:})")

    # if the job is in the queue wait for the signal
    out, error_state = check_running_state(id_job)
    if error_state == '':
        log.lock("About to Enter critical section")
        with job_state_lock:
            if not job_almost_done_flag:
                # wait for the signal or the alarm
                log.lock(f"Waiting on \'{job_type:s}\' job (id={id_job:})")
                signal.sigwait([signal.SIGUSR1, signal.SIGALRM])
                log.lock("Woke up")

            log.lock("About to leave the critical section")
    else:
        raise Exception(f"Undefined behaviour, the job state was:\n {error_state:s}")

    # Wait until the job's state file reflects a successful execution, or an error with execution
    while True:
        out, error_state = check_running_state(id_job)

        # if this string is in the job's state file then it successfully executed
        if "COMPLETED" in out:
            break

        # If the "COMPLETED" string is not in the job's state file AND scontrol reports no errors
        # the most likely cause is the state file has not been updated since the job left the queue
        # therefore we should wait until the state file is updated
        elif error_state == '':
            time.sleep(5)

        # if the error_state from scontrol is not empty then some undefined behaviour occurred
        else:
            raise Exception(f"Undefined behaviour, the job state was:\n {error_state:s}")
    return


def check_slurm_output(path_root, id_job):
    """checks ouput file from slurm for errors memory issues, incorrect arguments, etc"""

    # sanity check - sleep so that we are sure slurm has finished writing the output file
    time.sleep(10)

    slurm_path = os.path.join(path_root, f"slurm-{id_job:d}.out")
    log.debug(f"Checking slurm output:\n{slurm_path:}\n")

    with open(slurm_path, "r") as source_file:
        if not ("Finished" in source_file):
            log.debug("We got here too early?")  # should this be an exception?

        elif "SIGSEGV" in source_file:
            raise MemoryError(f"Job {id_job:d} had a Segmentation fault, see file {slurm_path:s}")

        elif "Killed" in source_file:
            # most likely cause is that it ran out of memory
            if "Exceeded step memory limit" in source_file:
                raise MemoryError(f"Job {id_job:d} ran out of memory, see file {slurm_path:s}")
            raise Exception(f"Job {id_job:d} failed for an unknown reason, see file {slurm_path:s}")
        else:
            log.warning(f"Undefined execution, check file {slurm_path:s} for issues")
    return


def extract_id_job_from_output(out):
    """ returns the job id inside the str argument 'out' or None if the job id is not present
    if no job id can be found then a warning is raised"""

    id_job = None
    if isinstance(out, str) and len(out) >= 21:  # this is hardcoded - possibly should be changed
        id_job = int(out[20:])
    else:
        log.warning(f"Not sure how to extract job id from\n{out}\n")

    return id_job


def run_locally(command):
    """ we just want to use subprocess instead of os.system
    to catch any weird errrors. We don't need to use sbatch for these calculations.
    """
    result = subprocess_run_wrapper(command, shell=True)

    # some kind of basic error checking
    try:
        result.check_returncode()
    except subprocess.CalledProcesssError as e:
        print("Return code non zero, script didnt' run properly")
        print(f"Stdout:{out}\nStderr{err}\n")
        raise e

    return result.stdout, result.stderr


def submit_job(command):
    """craft the job submission command - no error checking"""

    """ submits the job to the slurm server"""
    result = subprocess_run_wrapper(command, shell=True)
    result.check_returncode()

    id_job = extract_id_job_from_output(result.stdout)

    return id_job, result.stdout, result.stderr


def assert_partition_exists(partition_name):
    """only checks that the given string is listed as a partition by sinfo"""
    cmd = ['sinfo', '-O', 'partition']
    result = subprocess_run_wrapper(cmd)
    assert partition_name is not None, "The partition string is None?"
    assert partition_name in result.stdout, f"Partition {partition_name} is not present in {result.stdout}"
    return
