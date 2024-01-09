import itertools as it
import os
import sys
from os.path import abspath, join
from project_parameters import *

execute_script = "submit.sh"


def submit_jobs():

    for param_list in it.product(*expression_list):


        submit_dir = join(work_root, dir_string.format(project_name, *param_list))

        nof_pbfs = param_list[0]

        # manually select amount of memory for high PBF runs
        gbs = {
            #'op_water3Q_4st': {20: '1', 100: '10', 200: '10'},
            #'op_nh36Q_5st': {20: '1', 100: '10', 200: '10'},
            'op_ph33Q_3st': {20: '1', 100: '10', 200: '10'},
        }[project_name].get(nof_pbfs, '3')


        run_on_high_mem = '--partition=highmem' if (False or param_list[0] == 1000) else ''

        for operate_string in ["Ex", "Ey", "Ez"]:
            
            job_name = "PBF{:d}_tf{:.0f}_{:s}_{:s}".format(*param_list, project_name, operate_string)
            input_file_name = f"{project_name}_{operate_string}.inp"

            command = (
                f"sbatch --mem={gbs}GB "
                f"{run_on_high_mem} "
                f"--job-name={job_name} "
                f"--chdir={submit_dir} "
                f"{submit_dir}/{execute_script} {input_file_name}"
            )
            
            # print(command)
            # sys.exit(0)

            os.system(command)



if __name__ == "__main__":
    submit_jobs()
