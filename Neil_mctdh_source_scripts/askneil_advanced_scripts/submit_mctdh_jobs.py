import itertools as it
import os
from os.path import abspath, join
from project_parameters import *

execute_script = "submit_mctdh.sh"


def submit_jobs():
    for FC_or_not, coupling_order in it.product(["FC", "vibronic"], ["constant", "linear", "quadratic"]):
        if FC_or_not == "FC":
            continue

        if coupling_order == "constant":
            continue

        if coupling_order == "quadratic":
            continue

        # debug select specific types of jobs (set True)
        # if True and (FC_or_not != "FC"):
        if False and (FC_or_not != "vibronic" or coupling_order != "quadratic"):
            continue

        name = "_".join([project_name, FC_or_not, coupling_order])
        for param_list in it.product(*expression_list):
            job_name = "BF{:d}_tf{:.0f}_{:s}".format(*param_list, name)
            submit_dir = join(work_root, dir_string.format(name, *param_list))

            nof_pbfs = param_list[0]

            gbs = {
                'h2o2': {10: '1', 100: '10', 200: '10', 300: '10', 400: '15', 500: '20'},
                'hcooh': {100: '15', 200: '15', 300: '20', 400: '20', 500: '20'},
                'furan': {50: '20', 100: '40', 200: '45', 300: '50', 400: '60', 500: '70'},
                'formamide': {100: '15', 200: '15', 300: '15', 400: '15', 500: '15'},
                'vcm': {100: '30', 200: '40', 300: '50', 400: '60', 500: '70'},
                'hexahelicene': {15: '10', 200: '40', 300: '50', 400: '60', 500: '70'}
            }[project_name][nof_pbfs]

            run_on_high_mem = '--partition=highmem' if (False or param_list[0] == 1000) else ''

            command = (
                f"sbatch --mem={gbs}GB "
                f"{run_on_high_mem} "
                f"--job-name={job_name} "
                f"--chdir={submit_dir} "
                f"{submit_dir}/{execute_script} {inp_file_name}"
            )
            # print(command)
            # sys.exit(0)
            os.system(command)


if __name__ == "__main__":
    submit_jobs()
