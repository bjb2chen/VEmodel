import itertools as it
import os
import sys
from os.path import abspath, join
# from project_parameters import *

execute_script = "submit.sh"

name_of_model = "RhF3"

z_order_list = [
    1,
    2,
    3,
]

h_order_list = [
    0,
    1,
    2,
]

tf_list = [
    20,
    50,
    #100,
    #200,
    #500,
]


expression_list = [
    z_order_list,
    h_order_list,
    tf_list,
]

# submit_dir = "./"
# execute_script = "run_clean.sh"

SOC_flag = True


def submit_jobs():

    for param_list in it.product(*expression_list):

        nZ, nH, tf = param_list

        if SOC_flag:
            job_name = "{}_Z{:d}_H{:d}_{:d}tf_SOC".format(name_of_model, nZ, nH, tf)
        else:
            job_name = "{}_Z{:d}_H{:d}_{:d}tf".format(name_of_model, nZ, nH, tf)

        command = (
            "sbatch "
            # f"--mem={gbs}GB "
            f"--job-name={job_name} "
            # f"--chdir={submit_dir} "
            "run_clean.sh "
            f" {nZ} {nH} {tf} {job_name} "
            f"{os.getcwd()} "
            # f"{submit_dir}/{execute_script}"
        )

        print(command)
        #sys.exit(0)

        os.system(command)


if __name__ == "__main__":
    submit_jobs()
