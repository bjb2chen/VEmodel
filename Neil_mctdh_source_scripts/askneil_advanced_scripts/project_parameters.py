import socket
import os
from os.path import abspath, join

# 30 pbf ~ 1GB
# 100 pbf ~ 2GB
# 500 pbf ~ 8GB
# 1000 pbf ~ 15GB

# system_names = ["h2o", "ch2o", "co2", "nh3", ]
# medium_system_names = ["ch4", "n2o", "hooh", "bf3", ]

# pbfs = [10, 100, 200, 300, 400, 500]
# pbfs = [100, 200, 300, 400, 500]
# pbfs = [100, 200, 300]
# pbfs = [50, 100]
pbfs = [15, ]


# t_final = [10.0, 50.0]
t_final = [10.0, ]

expression_list = [pbfs, t_final]

# root directory
parent_project = "t_amplitudes_may2021"
# new models
# project_name, A, N = "h2o2", 6, 6
# project_name, A, N = "hcooh", 7, 9
# project_name, A, N = "furan", 9, 21
# project_name, A, N = "formamide", 7, 11
# project_name, A, N = "vcm", 7, 12
project_name, A, N = "hexahelicene", 15, 63
user_root = abspath("/ngraymon/projects/mctdh_compare")
home_root = abspath(f"/home/{user_root}/{parent_project}/{project_name}/")
work_root = abspath(f"/work/{user_root}/{parent_project}/mctdh/{project_name}/")

hostname = socket.gethostname()
server_flag = (hostname in ["nlogn", "feynman", 'dev002', ])
assert server_flag  # make sure we are on server
os.makedirs(work_root, exist_ok=True)  # make sure the root directory exists

# definitions
# the source files are located in the home directory by convention
op_file_name = f"{project_name}.op"
inp_file_name = f"{project_name}.inp"
execution_script = "submit_mctdh.sh"

root_path_op_files = home_root
os.makedirs(root_path_op_files, exist_ok=True)
src_path_original_op_file = join(home_root, "..", "model_src", op_file_name)
src_path_execution_script = join(home_root, "..", execution_script)

dir_string = "{:s}_PBF{:d}_tf{:.2f}"
