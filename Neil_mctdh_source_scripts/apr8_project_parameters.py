import socket
import os
from os.path import abspath, join

# 30 pbf ~ 1GB
# 100 pbf ~ 2GB
# 500 pbf ~ 8GB
# 1000 pbf ~ 15GB


# pbfs = [10, 30, 40, 50]
pbfs = [20,]
#tfinal = [10.0, 20.0, 30.0, ]
tfinal = [10.0, 20.0, 250.0 ]

expression_list = [pbfs, tfinal, ]


# root directory
# parent_project = "ground_state_energies"

# new models
# project_name, A, N = "h2o2", 6, 6
# project_name, A, N = "hcooh", 7, 9
# project_name, A, N = "furan", 9, 21
# project_name, A, N = "formamide", 7, 11
# project_name, A, N = "vcm", 7, 12

# project_name, A, N = "op_water3Q_4st", 4, 3
# project_name, A, N = "op_TiCl49Q_5st", 5, 9
project_name, A, N = "op_Caffeine69Q_3st", 3, 69

# \home\bjb2chen\740_project\mctdh_source_scripts\project_parameters.py
#user_root = abspath("/bjb2chen/740_project/claire_tiara_spectra/water_Claire_2")
#user_root = abspath("/bjb2chen/project_claire_tiara/ch2nh_vibronic_tiara")
user_root = abspath("/bjb2chen/CHEM494/lynne_calc_caffeine")
home_root = abspath(f"/home/{user_root}/{project_name}/")
work_root = abspath(f"/work/{user_root}/mctdh/{project_name}/")

# server_flag = (socket.gethostname() == "nlogn") or (socket.gethostname() == "feynman")
# assert server_flag  # make sure we are on server
os.makedirs(work_root, exist_ok=True)  # make sure the root directory exists

# definitions
# the source files are located in the home directory by convention
op_file_name = f"{project_name}.op"
# inp_file_name = f"{project_name}.inp"
execution_script = "submit.sh"


root_path_op_files = home_root
os.makedirs(root_path_op_files, exist_ok=True)

src_path_original_op_file = join(home_root, op_file_name)
src_path_execution_script = join('./', execution_script)

dir_string = "{:s}_PBF{:d}_tf{:.2f}"
