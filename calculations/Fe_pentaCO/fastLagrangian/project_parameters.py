import os
import socket
import itertools as it
from os.path import abspath, join

# -------------------------------------------------------------------------
#               set global project/calculation parameters
# -------------------------------------------------------------------------

# pick which molecule we're calculating
# (A = number of states, Z = number of atoms)

# project_name, A, Z = "h2o2", 6, 4
# project_name, A, Z = "hcooh", 7, 5
# project_name, A, Z = "furan", 9, 9
# project_name, A, Z = "formamide", 7, 6
# project_name, A, Z = "op_NH35Q_3st", 3, 4
project_name, A, Z = "op_FeCO3Q_5st", 5, 11

# the total number of modes (including translational, rotational, vibrational)
N_tot = 3 * Z

""" The file_name parameters:
    {name}_{basis_set}_{calculation_type}_{point_group_symmetry}_{nof_excited_states}st_diab
    {}_{}_{}_{}_{}st_diab
"""
name, basis_set, calculation_type, point_group_symmetry, nof_excited_states = (
    # "H2Ocat", "ccd", "gmcpt", "C1", 3
    #"H2Ocat", "ccd","gmcpt", "C1", 3
    "FeCOcat", "SPK", "gmcpt", "C1", 5
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

# -------------------------------------------------------------------------
filnam = file_name  # alias (remove later)
# -------------------------------------------------------------------------
#               set global project/calculation parameters
# -------------------------------------------------------------------------
selected_mode_list = {
    "H2Ocat": [7, 8, 9],
    # "NH3cat": [8, 9, 10, 11, 12, ],
    "NH3cat": [7, 8, 9, 10, 11, 12, ],
    "FeCOcat": [29, 30, 31, 32, 33, ],
    }[name]

nof_displacements_per_mode = {
    "H2Ocat": [10, 4, 2],
    # "NH3cat": [8, 3, 2, 2, 2],
    "NH3cat": [2, 2, 2, 2, 2, 2],
    "FeCOcat": [2, 2, 2, 2, 2],
}[name]

# -------------------------------------------------------------------------
#                       define all global flags
# -------------------------------------------------------------------------
SOC_flag = False
VECC_flag = False

# if True it doesn't run any subprocess.run() commands
# instead it simply prints out what the commands would be to the terminal
dry_run = False

suppress_zeros = True

# -------------------------------------------------------------------------
#                       define spectra parameters
# -------------------------------------------------------------------------
# pbfs = [10, 30, 40, 50]
pbfs = [50, ]
tfinal = [2000.0, ]

expression_list = [pbfs, tfinal, ]

# 30 pbf ~ 1GB
# 100 pbf ~ 2GB
# 500 pbf ~ 8GB
# 1000 pbf ~ 15GB

# -------------------------------------------------------------------------
#               useful indexing/mapping dictionaries
# -------------------------------------------------------------------------
""" This maps array indices to the label of the selected modes.
    The 1st mode is indexed by 0 in the array
    The 2nd mode is indexed by 1 in the array
    <...>
    For example if selected_mode_list is
        [7, 8, 9]
    then
    mode_map_dict = {0: 7, 1: 8, 2: 9}
"""
# build mode label objects

#  (ASSUMES YOUR MODES ARE IN INCREASING ORDER)
assert sorted(selected_mode_list) == selected_mode_list, f"{selected_mode_list=} is not sorted"
N = len(selected_mode_list)  # the number of modes should be

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


mode_map_dict = {k: v for k, v in enumerate(selected_mode_list)}


# -------------------------------------------------------------------------
# modes which have an order higher than 2
# all the +3/-3, +4/-4 displacements are handled by fitting instead of `mctdh()`
# this is a list of the indicies corresponding to those modes where their order is > 2
fitting_mode_idxs = [k for k in mode_map_dict.keys() if nof_displacements_per_mode[k] > 2]
N_fit = len(fitting_mode_idxs)

# -------------------------------------------------------------------------
# map (i,j)-> (n_i, n_j) numbers from selected_mode_list  (0, 1) -> (7, 8)
ij_map = {}
for i, j in it.combinations(range(N), 2):
    ij_map[(i, j)] = (mode_map_dict[i], mode_map_dict[j])

# -------------------------------------------------------------------------
""" map (n_i, n_j) -> (i,j)-  (7, 8) -> (0, 1)
returns the 0-indexed keys (0, 1) corresponding to the given literal-mode label keys (7, 8)
used commonly when accessing upper triangle indices of 2-D arrays/tensors
"""
reverse_ij_map = {}
for key, value in ij_map.items():
    reverse_ij_map[value] = key

# -------------------------------------------------------------------------
#                           Project Paths
# -------------------------------------------------------------------------
user_root = abspath(os.getcwd().replace('/scratch', '') if '/scratch' in os.getcwd() else os.getcwd()) # format is /user/.../*
home_root = abspath(f"/scratch/{user_root}/{project_name}/")
#work_root = abspath(f"/work/{user_root}/mctdh/{project_name}/")
work_root = abspath(f"/scratch/{user_root}/mctdh/{project_name}/")

# server_flag = (socket.gethostname() == "nlogn") or (socket.gethostname() == "feynman")
# assert server_flag  # make sure we are on server
os.makedirs(work_root, exist_ok=True)  # make sure the root directory exists

# -------------------------------------------------------------------------
#                   path/filename  definitions/conventions
# -------------------------------------------------------------------------
# the source files are located in the home directory by convention
op_file_name = f"{project_name}.op"
# inp_file_name = f"{project_name}.inp"
execution_script = "submit.sh"


root_path_op_files = home_root
os.makedirs(root_path_op_files, exist_ok=True)

src_path_original_op_file = join(home_root, op_file_name)
src_path_execution_script = join('./', execution_script)

# -------------------------------------------------------------------------
dir_string = "{:s}_PBF{:d}_tf{:.2f}"


# -------------------------------------------------------------------------
# define constant values
from types import SimpleNamespace

# constants regarding the calculation
gamess_const = SimpleNamespace()  # make an empty object
gamess_const.qsize = 0.05

# Quantum Mechanical constants
QM_const = SimpleNamespace()  # make an empty object
QM_const.ha2ev = 27.2113961318
QM_const.ev2ha = 1. / QM_const.ha2ev
QM_const.wn2ev = 0.000123981
QM_const.wn2eh = 0.00000455633
QM_const.ang2br = 1.889725989
QM_const.amu2me = 1822.888

# -------------------------------------------------------------------------
