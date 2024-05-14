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
            f'{project_name}': {20: '1', 100: '10', 200: '10'},
        }[project_name].get(nof_pbfs, '3')


        cc_mctdh_time = f'--time=47:00:00' if is_compute_canada else ''

        for operate_string in range(A+1, A+2):

            job_name = "PBF{:d}_tf{:.0f}_{:s}_{:d}".format(*param_list, project_name, operate_string)
            input_file_name = f"{project_name}_init_st{operate_string}.inp"

            command = (
                f"sbatch "
                f"{cc_mctdh_time} "
		f"--nodes=1 "
		f"--ntasks=64 "
		f"--mem=0 "
                f"--job-name={job_name} "
                f"--chdir={submit_dir}/init_st{operate_string} "
                f"{submit_dir}/{execute_script} {input_file_name}"
            )

            #print(command)
            # sys.exit(0)

            os.system(command)



if __name__ == "__main__":
    submit_jobs()
