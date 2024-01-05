# system imports
import shutil
import itertools as it
import sys
import os

# third party imports

# local imports
from project_parameters import *

execute_script = "submit_cc.sh"


def submit_jobs():
    for FC_or_not, coupling_order in it.product(["FC", "vibronic"], ["constant", "linear", "quadratic"]):
        if FC_or_not == "vibronic" and coupling_order == "constant":
            continue

        if coupling_order != "linear":
            continue

        name = "_".join([project_name, FC_or_not, coupling_order])
        for param_list in it.product(*expression_list):
            submit_dir = join(work_root, dir_string.format(name, *param_list))
            shutil.copy("./cc_script.py", submit_dir)
            shutil.copy("./submit_cc.sh", submit_dir)
            command = f"sbatch --mem=5GB --partition=highmem --chdir={submit_dir} {execute_script}"
            print(command)
            # os.system(command)


if __name__ == "__main__":
    submit_jobs()
