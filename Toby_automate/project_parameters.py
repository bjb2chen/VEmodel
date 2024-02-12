import socket
import os
from os.path import abspath, join

# 30 pbf ~ 1GB
# 100 pbf ~ 2GB
# 500 pbf ~ 8GB
# 1000 pbf ~ 15GB


# pbfs = [10, 30, 40, 50]
pbfs = [30, ]
tfinal = [100.0, ]

expression_list = [pbfs, tfinal, ]


# root directory
# parent_project = "ground_state_energies"

# -------------------------------------------------------------------------
# pick which molecule we're calculating

# project_name, A, Z = "h2o2", 6, 4
# project_name, A, Z = "hcooh", 7, 5
# project_name, A, Z = "furan", 9, 9
# project_name, A, Z = "formamide", 7, 6
# project_name, A, Z = "vcm", 7, 6
# project_name, A, Z = "op_water3Q_4st", 4, 3
project_name, A, Z = "op_H2O3Q_3st", 3, 3       # A = number of states, Z = number of atoms

# the total number of modes (including translational, rotational, vibrational)
N_tot = 3 * Z

""" The file_name parameters:
    {name}_{basis_set}_{calculation_type}_{point_group_symmetry}_{nof_excited_states}_diab
    {}_{}_{}_{}_{}st_diab
"""
name, basis_set, calculation_type, point_group_symmetry, nof_excited_states = (
    "H2Ocat", "cct", "gmcpt", "C1", 3
    # "NH3anion", "cctzla", "rohf", "C2v", 4
    # "H2Ocat", "cct", "gmcpt", "C1", 3
    # "NH3anion", "cctzla", "rohf", "C2v", 4
    # "H2Ocat", "cct", "gmcpt", "C1", 3
)

file_name = str(
    f"{name}_"
    f"{basis_set}_"
    f"{calculation_type}_"
    f"{point_group_symmetry}_"
    f"{nof_excited_states}st_"
    "diab"
)
print(f"Our {file_name=}")

filnam = file_name  # alias (remove later)

if False:
    project_name = file_name
# -------------------------------------------------------------------------

# build mode label objects
if False:  # if needed in future
    _all_modes = set(range(1, N_tot+1))
    if (select_modes := True):
        _selected_modes = set([7, 8, 9])
        _excluded_modes = _all_modes - _selected_modes
    else:
        _excluded_modes = set([1, 2, 3, 4, 5, 6, ] + [10, 11, 12, ])
        _selected_modes = _all_modes - _excluded_modes

    # sort the labels of the modes (sets are unordered)
    selected_mode_list = sorted([*_selected_modes])

else:
    selected_mode_list = [7, 8, 9, ]

#  (ASSUMES YOUR MODES ARE IN INCREASING ORDER)
assert sorted(selected_mode_list) == selected_mode_list, f"{selected_mode_list=} is not sorted"

""" This maps array indices to the label of the selected modes.
    The 1st mode is indexed by 0 in the array
    The 2nd mode is indexed by 1 in the array
    <...>
    For example if selected_mode_list is
        [7, 8, 9]
    then
    mode_map_dict = {0: 7, 1: 8, 2: 9}
"""
mode_map_dict = {k: v for k, v in enumerate(selected_mode_list)}
N = len(selected_mode_list)  # the number of modes should be

# -------------------------------------------------------------------------


# Project Paths
# user_root = abspath("/bjb2chen/gamess/vibronics/template_examples/")        # format is /user/.../*
user_root = abspath("/bjb2chen/gamess/vibronics/template_examples/H2O")        # format is /user/.../*
home_root = abspath(f"/home/{user_root}/home/{project_name}/")
work_root = abspath(f"/home/{user_root}/work/mctdh/{project_name}/")

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
# -------------------------------------------------------------------------

# set conversion constants
qsize = 0.05
ha2ev = 27.2113961318
wn2ev = 0.000123981
wn2eh = 0.00000455633
ang2br = 1.889725989
amu2me = 1822.888

# -------------------------------------------------------------------------
# if True it doesn't run any subprocess.run() commands
# instead it simply prints out what the commands would be to the terminal
dry_run = True

# -------------------------------------------------------------------------
