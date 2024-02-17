""" x """

# system imports
import os
from os.path import splitext
import sys
import types
import subprocess
import shutil
import re
# import json
import functools
import itertools as it

# third party imports
import numpy as np
import pprint


# local packages
import project_parameters as pp
from project_parameters import *  # eventually remove this


# ---------------------------------------------------------------------------------------
_possible_mctdh_headers = [
    'Frequencies',
    'Electronic Hamitonian',  # NOTE - THIS IS NOT A SPELLING ERROR (Hamitonian)
    'Electronic transition moments',
    'Magnetic transition moments',
    'Linear Coupling Constants',
    'Diagonal Quadratic Coupling Constants',
    'Off_diagonal Quadratic Coupling Constants',
    'Cubic Coupling Constants',
    'Quartic Coupling Constants',
    'Quintic Coupling Constants',
    'Sextic Coupling Constants',
    'Bi-Cubic Constants',
    'Bi-Quartic Constants',
    'end-parameter-section',
]

# ---------------------------------------------------------------------------------------


def subprocess_run_wrapper(*args, **kwargs):
    """ Subprocess.run() returns the value from std.in """

    assert len(args) == 1
    command = args[0]
    assert isinstance(command, str) or isinstance(command, list), (
        f"Bad arguments! First arg should only be str or list you provided {type(command)}"
    )

    if False and __debug__:  # for checking
        print(args)
        print(command)
        print(kwargs)

        breakpoint()

    if pp.dry_run:
        if isinstance(command, list):
            command = " ".join(command)+'\n'

        print(command)
        # fake a return value
        return_obj = types.SimpleNamespace()
        return_obj.returncode = 0
        return return_obj
    else:
        return subprocess.run(command, **kwargs)


def subprocess_call_wrapper(*args, **kwargs):
    """ subprocess.call only returns 1/0  """

    assert len(args) == 1
    command = args[0]
    assert isinstance(command, str) or isinstance(command, list), (
        f"Bad arguments! First arg should only be str or list you provided {type(command)}"
    )

    if False and __debug__:  # for checking
        print(args)
        print(command)
        print(kwargs)
        breakpoint()

    if pp.dry_run:
        if isinstance(command, list):
            command = " ".join(command)+'\n'
        print(command)

        # fake a return value
        return_obj = types.SimpleNamespace()
        return_obj.returncode = 0
        return return_obj
    else:
        return subprocess.call(command, **kwargs)


def os_system_wrapper(*args, **kwargs):
    """ os.system is needed for submitting jobs on nlogn!
        sbatch on nlogn only works with os
    """
    assert len(args) == 1
    command = args[0]

    if False and __debug__:  # for checking
        print(args)
        print(command)
        print(kwargs)
        breakpoint()

    if pp.dry_run:
        print(f"{command}\n")
        return
    else:
        return os.system(command, **kwargs)


def upper_triangle_loop_indices(max_num, number_of_indicies=2):
    """ Returns generator function providing indices to iterate over upper-triangle elements of a tensor/matrix.
    Only guaranteed to work for 2D matrix ... should work for 3D but haven't checked for certain yet.

    see stackoverflow link below for an explanation on how to do combination loops
    https://stackoverflow.com/questions/69321896/how-to-do-dependent-nested-loop-using-itertools-product
    """
    assert number_of_indicies >= 2, "Don't use this fxn for 1D arrays / lists"
    return it.combinations(range(max_num), number_of_indicies)

# ---------------------------------------------------------------------------------------


linear_disp_keys = ["+1", "+2", "-1", "-2"]
bi_linear_disp_keys = ["++", "+-", "-+", "--"]

linear_disp_suffix = {
    "+1": 'plus',
    "+2": 'minus',
    "-1": 'plusx2',
    "-2": 'minusx2',
}

bi_linear_disp_suffix = {
    "++": 'pp',
    "+-": 'pm',
    "-+": 'mp',
    "--": 'mm',
}

# these could be lists or dicts... for now I'm just making them dicts
linear_disp_filenames = {
    k: f'dist_structure_{linear_disp_suffix[k]}'
    for k in linear_disp_keys
}

bi_linear_disp_filenames = {
    k: f'dist_structure_{bi_linear_disp_suffix[k]}'
    for k in bi_linear_disp_keys
}


def _remove_existing_distorted_structure_files(filename_dict):
    """ Delete existing distorted structure files """
    for filename in filename_dict.values():
        try:
            subprocess_run_wrapper(['rm', '-f', filename])
        except Exception as e:
            print(f"Error deleting {filename}: {str(e)}")
    return


# ---------------------------------------------------------------------------------------


def my_subgam(path, **kwargs):
    """ Create our own GAMESS job submission script
        Recording the script inside .slurm helps for recordkeeping """

    # if its easier to just change project parameters i would recommend doing
    # ncpus = pp.ncpus
    # nhour = pp.nhour
    # ngb = pp.ngb

    ncpus = kwargs.get('ncpus', 2)
    nhour = kwargs.get('nhour', 1)
    ngb = kwargs.get('ngb', 2)

    # Remove the ".inp" extension from the filename
    input_no_ext, extension = splitext(path)
    print(f"running calculations for {input_no_ext}")

    # wd = os.getcwd()

    file_contents = "".join([
        "#!/bin/bash\n",
        "#SBATCH --nodes=1\n",
        f"#SBATCH --ntasks={ncpus}\n",
        f"#SBATCH --mem-per-cpu={ngb}G\n",
        f"#SBATCH --time={nhour}:00:00\n",
        "\n",
        "cd $SLURM_SUBMIT_DIR\n",
        "\n",
        "export SLURM_CPUS_PER_TASK\n",
        'mkdir -p /home/$USER/.gamess_ascii_files/$SLURM_JOBID\n',
        f"/home/$USER/LOCAL/runG_diab {input_no_ext}.inp {ncpus} \n",
    ])

    with open(f"{input_no_ext}.slurm", "w") as slurm_file:
        slurm_file.write(file_contents)

    return f"{input_no_ext}.slurm"


def diabatization(**kwargs):
    """ Preform diabatization?

    Returns a dictionary of dictionaries `dist_coord`
    Whose keys are `displacement_keys` and `bi_linear_keys`.
    Each of those individual dictionaries store ( ... ?)

    (other info)
    """

    # -------------------------------------------------------------------------
    # unpacking (should get rid of these two lines eventually)
    # modes_included = kwargs['modes_included']
    filnam = filename = kwargs['project_filename']

    # kwargs has a bunch more items inside it
    freqcm = kwargs.get('freqcm')
    ndim = kwargs.get('ndim')
    refcoord = kwargs.get('refcoord')
    nrmmod = kwargs.get('nrmmod')
    natoms = kwargs.get('natoms')
    atmlst = kwargs.get('atmlst')
    chrglst = kwargs.get('chrglst')
    qsize = kwargs.get('qsize', 0.05)
    ha2ev = kwargs.get('ha2ev', 27.2113961318)
    wn2ev = kwargs.get('wn2ev', 0.000123981)
    wn2eh = kwargs.get('wn2eh', 0.00000455633)
    ang2br = kwargs.get('ang2br', 1.889725989)
    amu2me = kwargs.get('amu2me', 1822.88839)

    # if its easier to just change project parameters i would recommend doing
    # amu2me = pp.amu2me
    # <...> (etc.)
    # since you have to use these values everywhere I take back what I said about using kwargs
    # probably best to just use project parameters
    # -------------------------------------------------------------------------
    # Precomputed constants

    """
    it might be worth while to precompute constants here?
    for example in `_convert_qsize_to_rsize` we compute
        q / (pow(amu2me, 0.5) * ang2br * pow(omega * wn2eh, 0.5))

    this could be replace with
        q / qsize_to_rsize_conversion_factor[i]

    where we compute the factor for all modes included before hand
        qsize_to_rsize_conversion_factor[i] = [
            pow(amu2me, 0.5) * ang2br * pow(freqcm[i] * wn2eh, 0.5)
            for i in modes_included
        ]
    """

    # -------------------------------------------------------------------------

    # prepare the displaced/distored co-ordinates (the object we return)
    disp_coord = {}
    disp_coord.update({k: {} for k in linear_disp_keys})
    disp_coord.update({k: {} for k in bi_linear_disp_keys})

    # -------------------------------------------------------------------------
    # NEIL's mappings

    """ Note the indexing here!!
    For a system with 9 modes (and therefore 3 atoms, since Z * 3 = 9 = N_tot)
    there will be 9 frequencies `len(freqcm)` = 9
    So if `selected_mode_list` is [7, 8, 9]
    then `s_idx` will be [6, 7, 8]
    since `freqcm` is a length 9 array indexed at 0, `[6, 7, 8]` will return the 7th, 8th, 9th elements
    """
    s_idx = [i-1 for i in selected_mode_list]  # need to shift by 1 for 0-indexed arrays
    freq_array = freqcm[s_idx]  # freqcm[[6,7,8]] returns 7,8,9th freq values
    mode_array = nrmmod[:, s_idx]  # return 7,8,9th columns of Z*3 by Z*3 array, assume we return array
    atom_list = [*atmlst.values()]
    charge_list = [*chrglst.values()]

    # make a list for easy calculations
    ref_coord_array = np.array([*refcoord.values()])

    nof_surfaces = A  # alias
    nof_modes = N  # alias
    nof_atoms = Z  # alias

    NEW = np.newaxis  # alias

    if ((Z*3)**2 > 1e4):
        print("Maybe think about not precomputing? Check with Neil again?")
        breakpoint()
        import sys; sys.exit()
        # ---------------------------------------------------------------------
        # may want to consider indexing the array?

        # throw away the modes we don't need
        column_index = [j-1 for j in selected_mode_list]

        # anything whose dimensionality was N_tot we need to reduce
        # (throw away everything except the columns you needed)
        mode_array = np.array(mode_array[:, column_index])
        freq_array = np.array(freq_array[:, column_index])
        breakpoint()
        # ---------------------------------------------------------------------

    """ precomputing the bilinear makes a (ndim, nof_modes, nof_modes) sized numpy array
    for each of the 4 ++, +-, -+, -- bilinear keys.
    That might take too much memory... not sure?
    If it does then just disable this flag and compute the values inside the loop
    """
    precompute_bilinear = True
    # -------------------------------------------------------------------------

    def _preconvert_qsize_to_rsize(omega, q):
        """ Convert the reduced dimensionless q to the actual rsize in sqrt(amu)*Angs unit.
        Assumes that
            amu2me, ang2br, wn2eh
        are all constants imported from project_parameters.
        """
        R = q / (pow(amu2me, 0.5) * ang2br * pow(omega[:] * wn2eh, 0.5))
        return R

    def _precompute_linear_displacements(R_array, mode_array, reference, distored_coords):
        """
        `reference` is (Z*3, ) ~ x,y,and z coordinates for each atom
        `R_array` is assumed to be of dimension (N) ~ 1 value of displacement for all modes
        `mode_array` is assumed to be of dimension (Z*3, N)
        Z*3 == 9

        Use broadcasting.

        The first dimension iterates over the atoms and their xyz coordinates.
        The second dimension iterates over the normal modes.
        """
        Z, N = pp.Z, pp.N

        # the `:` of the second dimension of mode_array is what we index with i or j
        displacements = {
            "+1": reference[:, NEW] + 1.0 * R_array[NEW, :] * mode_array[:, :],
            "-1": reference[:, NEW] - 1.0 * R_array[NEW, :] * mode_array[:, :],
            "+2": reference[:, NEW] + 2.0 * R_array[NEW, :] * mode_array[:, :],
            "-2": reference[:, NEW] - 2.0 * R_array[NEW, :] * mode_array[:, :],
            #              (Z*3, 1)                 (1,  N)            (Z*3, N)
        }
        assert set(displacements.keys()) == set(linear_disp_keys), f"{linear_disp_keys=} no longer agree!"

        shape = (Z*3, N)
        for k in linear_disp_keys:
            assert displacements[k].shape == shape, f"{k=} {displacements[k].shape=} not {shape=}?"

        # store the displacements in the `distored_coords` dictionary
        for key in linear_disp_keys:
            distored_coords[key] = displacements[key]

        return  # don't need to return value, as `distored_coords` is a dictionary

    def _precompute_bilinear_displacements(R_array, mode_array, distored_coords):
        """
        `distored_coords` is assumed to be of dimension (Z*3, N)
        `R_array` is assumed to be of dimension (Z*3)
        `mode_array` is assumed to be of dimension (Z*3, N)
        Use broadcasting.

        The first dimension iterates over the atoms and their xyz coordinates.
        The second dimension iterates over the normal modes (i).
        The third dimension iterates over the normal modes (j).
        """
        Z, N = pp.Z, pp.N

        n = R_array.shape[0]
        assert n == mode_array.shape[1]

        displacements = {
            "++": distored_coords["+1"][:, :, NEW] + R_array[NEW, NEW, :] * mode_array[:, NEW, :],
            "+-": distored_coords["+1"][:, :, NEW] - R_array[NEW, NEW, :] * mode_array[:, NEW, :],
            "-+": distored_coords["-1"][:, :, NEW] + R_array[NEW, NEW, :] * mode_array[:, NEW, :],
            "--": distored_coords["-1"][:, :, NEW] - R_array[NEW, NEW, :] * mode_array[:, NEW, :],
            #                          (Z*3, N, 1)           (1,  1,   N)             (Z*3, 1, N)
        }
        assert set(displacements.keys()) == set(bi_linear_disp_keys), f"{bi_linear_disp_keys=} no longer agree!"

        shape = (Z*3, N, N)
        for k in bi_linear_disp_keys:
            assert displacements[k].shape == shape, f"{k=} {displacements[k].shape=} not {shape=}?"

        # store the displacements in the `distored_coords` dictionary
        for key in bi_linear_disp_keys:
            distored_coords[key] = displacements[key]

        return  # don't need to return value, as `distored_coords` is a dictionary

    def _save_distorted_structure(mode_idx, displaced_q, charge_list, atom_list, filename_list, key_list):
        """ save the distorted structure to the `distored_structure_filenames`
        `displaced_q` is an array of dimension (ndim, ndim)

        mode_idx can be (i, ) or (i, j, )
        so that
            d[(offset+0, *mode_idx)]
        is like
            d[a, i, j]
            or
            d[a, i]
        where the first index is picking the atom xyz
        and the second index is the first modes xyz
        and the third index is the second modex xyz
        """

        header, data = "{:<2s} {:} ", "{: .10f} " * natoms
        template_string = header + data

        for key in key_list:
            file_contents = []
            d = displaced_q[key]
            for idx_atom in range(natoms):
                offset = 3*idx_atom  # 3 xyz-coordinates

                if False and len(mode_idx) > 1:
                    print(*mode_idx, d.shape); breakpoint()

                string = template_string.format(
                    atom_list[idx_atom], charge_list[idx_atom],
                    d[(offset+0, *mode_idx)],  # x component
                    d[(offset+1, *mode_idx)],  # y component
                    d[(offset+2, *mode_idx)],  # z component
                )
                file_contents.append(string)

                if False and __debug__:  # if you need to debug the string
                    print(string)

            if True and __debug__:  # if you need to debug the file_contents
                print(f"MODE: {str(mode_idx)}", *file_contents, sep='\n')
                # breakpoint()

            # write to file
            path = filename_list[key]
            with open(path, 'w') as fp:
                fp.write("\n".join(file_contents))

        return

    def _create_linear_diabatization_input_files(i, filename, qsize):
        """ this could be improved somewhat

        The `q1_label` goes from (1 -> nof_modes+1) rather than (0 -> nof_modes).
        """
        q1_label = pp.mode_map_dict[i]  # returns the label of the mode at the array's ith index
        refG_out = f"{filename}_refG.out"

        # Check if the reference geometry calculation is done?
        grace0 = subprocess_run_wrapper(["grep", "DONE WITH MP2 ENERGY", refG_out])
        ref_geom_flag_exists = bool(grace0.returncode == 0)

        for d1 in ['+', '-']:
            p_or_m = {'+': 'plus', '-': 'minus'}[d1]
            for suffix in ['', 'x2']:
                games_filename = f'{filename}_mode{q1_label}_{d1}{qsize}{suffix}'
                distored_struct_file = f'dist_structure_{p_or_m}{suffix}'

                shutil.copy('temp.inp', games_filename+'.inp')

                # if your going to do this why not just copy the file then append to it?
                # with open(games_inp_file, 'a') as inp_file:
                #     with open(distored_struct_file, 'r', errors='replace') as src_fp:
                #         inp_file.write(src_fp.read())
                #         inp_file.write(' $END')

                # so you only write to this file to then change another file?
                # but i guess its good to have a record of these files?
                #
                with open(distored_struct_file, 'r', errors='replace') as fp:
                    data = fp.read()

                with open(games_filename+'.inp', 'a') as fp:
                    fp.write(data)  # can you just do data + ' $END' in one write?
                    fp.write('\n $END')

                grace1 = subprocess_run_wrapper(["grep", "DONE WITH MP2 ENERGY", games_filename+'.out'])
                gamess_calculation_not_run = bool(grace1.returncode != 0)

                # This means that refG completed successfully and `diabmode*.out` not completed
                if (ref_geom_flag_exists and gamess_calculation_not_run) or pp.dry_run:
                    print(f"Running calculations for {games_filename}")
                    try:
                        os_system_wrapper("sbatch" + " " + my_subgam(games_filename+'.inp', ncpus=2, ngb=1, nhour=1))
                    except Exception as e:
                        print(f"Error running diabatization calculation: {str(e)}")
                else:
                    print(f"{games_filename} is done")

        return

    def _create_bilinear_diabatization_input_files(i, j, filename, qsize):
        """ this could be improved somewhat
        """
        q1_label = pp.mode_map_dict[i]  # returns the label of the mode at the array's i-th index
        q2_label = pp.mode_map_dict[j]  # returns the label of the mode at the array's j-th index

        refG_out = f"{filename}_refG.out"

        # Check if the reference geometry calculation is done?
        grace0 = subprocess_run_wrapper(["grep", "DONE WITH MP2 ENERGY", refG_out])
        ref_geom_flag_exists = bool(grace0.returncode == 0)

        for d1, d2 in it.product(['+', '-'], ['+', '-']):
            suffix = bi_linear_disp_suffix[d1+d2]

            games_filename = f'{filename}_mode{q1_label}_{d1}{qsize}_mode{q2_label}_{d2}{qsize}'
            distored_struct_file = f'dist_structure_{suffix}'

            shutil.copy('temp.inp', games_filename+'.inp')

            # so you only write to this file to then change another file?
            # but i guess its good to have a record of these files?
            with open(distored_struct_file, 'r', errors='replace') as fp:
                data = fp.read()

            with open(games_filename+'.inp', 'a') as fp:
                fp.write(data)  # can you just do data + ' $END' in one write?
                fp.write('\n $END')

            # Check if the calculation is done already
            grace2 = subprocess_run_wrapper(["grep", "DONE WITH MP2 ENERGY", games_filename + '.out'])
            gamess_calculation_not_run = bool(grace2.returncode != 0)

            # this will never work? grace0 is not defined
            if (ref_geom_flag_exists and gamess_calculation_not_run) or pp.dry_run:
                print(f"Running calculations for {games_filename}!")
                try:
                    os_system_wrapper("sbatch" + " " + my_subgam(games_filename+'.inp', ncpus=2, ngb=1, nhour=1))
                except Exception as e:
                    print(f"Error running diabatization calculation: {str(e)}")
            else:
                print(f"{games_filename} is already done.")

        return

    def _compute_bi_linear_displacements(i, j, R_array, mode_array, distored_coords):
        """
        The first dimension iterates over the atoms and their xyz coordinates
        The second dimension iterates over the normal modes
        """

        # the `:` elements are indexed by (c in range(1, ndim + 1))
        # the `:` of the second dimension of mode_array is what we index with j
        displacements = {  # numpy arrays of length `ndim`
            "++": distored_coords["+1"][:, i] + R_array[NEW, :] * mode_array[:, j],
            "+-": distored_coords["+1"][:, i] - R_array[NEW, :] * mode_array[:, j],
            "-+": distored_coords["-1"][:, i] + R_array[NEW, :] * mode_array[:, j],
            "--": distored_coords["-1"][:, i] - R_array[NEW, :] * mode_array[:, j],
        }

        for key in bi_linear_disp_keys:
            distored_coords[key] = displacements[key]

        return

    # -------------------------------------------------------------------------

    rsize = _preconvert_qsize_to_rsize(freq_array, qsize)

    # modify disp_coord in-place
    _precompute_linear_displacements(rsize, mode_array, ref_coord_array, disp_coord)

    # if precompute_bilinear:  # modify disp_coord in-place (may take a lot of memory?)
    _precompute_bilinear_displacements(rsize, mode_array, disp_coord)

    # Loop over modes (that you selected) and do linear displacements
    for i in range(N):

        _remove_existing_distorted_structure_files(linear_disp_filenames)
        index = (i, )
        _save_distorted_structure(
            index, disp_coord, charge_list, atom_list,
            linear_disp_filenames,
            linear_disp_keys
        )
        _create_linear_diabatization_input_files(i, filnam, qsize)

        # 2D distortion to get bilinear vibronic coupling
        for j in range(0, i):

            # if not precompute_bilinear:  # modify disp_coord in-place
            # _compute_bi_linear_displacements(i, j, rsize, mode_array, disp_coord)

            _remove_existing_distorted_structure_files(bi_linear_disp_filenames)
            # index = (i, j) if precompute_bilinear else (j, )
            index = (i, j)
            _save_distorted_structure(
                index, disp_coord, charge_list, atom_list,
                bi_linear_disp_filenames,
                bi_linear_disp_keys,
            )
            _create_bilinear_diabatization_input_files(i, j, filnam, qsize)

    return disp_coord


# ---------------------------------------------------------------------------------------
# Now we move on to extract vibronic coupling constants using finite difference
# and write the data in an mctdh operator file


def extract_ground_state_energy(hessout, pattern):
    try:
        command = f'grep "{pattern}" {hessout} | tail -1 | cut -c40-'
        result = subprocess_run_wrapper(command, shell=True, text=True, capture_output=True)
        output = float(result.stdout.strip().replace(" ", ""))
        return output
    except Exception as e:
        print("Cannot find ground state energy")
        return None


def refG_extract(file_path, pattern):
    try:
        # Use subprocess.run with the direct command
        command = f'grep "{pattern}" "{file_path}" | tail -1 | cut -c62-'
        result = subprocess_run_wrapper(command, shell=True, text=True, capture_output=True)

        # If there is output, convert it to float
        try:
            output = float(result.stdout.strip().replace(" ", ""))
            return output
        except Exception as e:
            with open(file_path, 'r', errors='replace') as file:
                for line in reversed(file.readlines()):
                    match = re.search(pattern, line)
                    if match:
                        return float(line[62:].strip().replace(" ", ""))
    except subprocess.CalledProcessError:
        # Return None if there is an error
        return None


def extract_diabatic_energy(file_path, pattern):
    try:
        # Use subprocess.run with the direct command
        command = f'grep "{pattern}" "{file_path}" | tail -1 | cut -c44-61'
        result = subprocess_run_wrapper(command, shell=True, text=True, capture_output=True)

        # If there is output, convert it to float
        try:
            output = float(result.stdout.strip().replace(" ", ""))
            return output
        except Exception as e:
            with open(file_path, 'r', errors='replace') as file:
                for line in reversed(file.readlines()):
                    match = re.search(pattern, line)
                    if match:
                        return float(line[44:62].strip().replace(" ", ""))
    except subprocess.CalledProcessError:
        # Return None if there is an error
        return None


def extract_coupling_energy(file_path, pattern):
    try:
        # Use subprocess.run with the direct command
        command = f'grep "{pattern}" "{file_path}" | tail -1 | cut -c62-' # we insert different patterns into refG_extract and this function
        result = subprocess_run_wrapper(command, shell=True, text=True, capture_output=True)

        # If there is output, convert it to float
        try:
            output = float(result.stdout.strip().replace(" ", ""))
            return output
        except Exception as e:
            with open(file_path, 'r', errors='replace') as file:
                for line in reversed(file.readlines()):
                    match = re.search(pattern, line)
                    if match:
                        return float(line[62:].strip().replace(" ", ""))
    except subprocess.CalledProcessError:
        # Return None if there is an error
        return None


def find_nstate(file_path, pattern='# of states in CI      = '):
    with open(file_path, 'r', errors='replace') as file:
        for line in file:
            if pattern in line:
                return int(line.split('=')[1].strip())
    return None  # Return None if the pattern is not found


def extract_ground_to_excited_state_transition_dipoles(selected_lines):
    ground_to_excited_state_transition_dipoles = {}

    for TDIPOLEline in selected_lines:
        try:
            if TDIPOLEline[0:5].strip().isnumeric():
                state1 = int(TDIPOLEline[0:5].strip())
                #print(state1)
                state2 = int(TDIPOLEline[5:10].strip())

                if (state2 == 1) and (state1 != 1):
                    x = float(TDIPOLEline[11:21].strip())
                    y = float(TDIPOLEline[22:31].strip())
                    z = float(TDIPOLEline[32:42].strip())
                    ground_to_excited_state_transition_dipoles[state1] = (f"{x:.6f}", f"{y:.6f}", f"{z:6f}")
        except Exception as e:
            print(f"ERror processing line: {TDIPOLEline} - {e}")

    return ground_to_excited_state_transition_dipoles


def extract_DSOME(fname, nof_states):

    start_pattern = "HSO MATRIX IN DIABATIC REPRESENTATION (DIRECT MAXIMIZATION)"
    end_pattern = 'SOC EIG. VALUES and VECTORS IN DIABATS (DIRECT MAX.)'

    sed_command = f"sed -n '/{start_pattern}/,/{end_pattern}/p' {fname}"
    result = subprocess_run_wrapper(sed_command, shell=True, text=True, capture_output=True)

    line_list = result.stdout.splitlines()
    # i believe there is another way to check if the result is empty?
    # I think the `result` object has error codes or some such?

    if isinstance(line_list, list) and line_list != []:
        selected_lines = result.stdout.splitlines()
    else:
        print('Cannot use sed to extract')
        selected_lines = extract_lines_between_patterns(  # fallback to slow pythonic extraction
            f'{fname}',
            f'{start_pattern}',
            f'{end_pattern}',
        )
        print(f'Using selected lines from {fname}, opened via python')

    DSOME_set = {}
    full_extracted_set = {}
    summed_set_real = {}
    summed_set_imag = {}
    append_J = {}

    for DSOMEline in selected_lines:
        if "STATE #" in DSOMEline:
            try:
                ist = DSOMEline[9:12].strip().replace(" ", "")
                jst = DSOMEline[14:16].strip().replace(" ", "")
                kst = ist + ' & ' + jst + ',' + DSOMEline[31:33]
                real = DSOMEline[48:61].strip().replace(" ", "")
                imaginary = DSOMEline[63:75].strip().replace(" ", "")

                if '*' in real:
                    real = 0
                if '*' in imaginary:
                    imaginary = 0

                DSOME_set[kst] = complex(float(real), float(imaginary))

            except Exception as e:
                print(f"Error processing line: {DSOMEline} - {e}")

    assert isinstance(nof_states, int), "can be removed after code is fully updated"

    # see hyperlink for a nice way to do combination loops
    # https://stackoverflow.com/questions/69321896/how-to-do-dependent-nested-loop-using-itertools-product
    # I also made the optional function `upper_triangle_loop_indices`
    for l_idx, r_idx in it.combinations(range(nof_states), 2):

        for level_idx in range(1, 3):  # this should be replaced with array assignment
            full_extracted_set[l_idx, r_idx, level_idx] = DSOME_set[f'{l_idx} & {r_idx}, {level_idx}']

    if False:  # old approach
        real_sum = full_extracted_set[l_idx, r_idx, 1].real + full_extracted_set[l_idx, r_idx, 2].real
        summed_set_real[l_idx, r_idx] = real_sum

        imag_sum = full_extracted_set[l_idx, r_idx, 1].imag + full_extracted_set[l_idx, r_idx, 2].imag
        summed_set_imag[l_idx, r_idx] = imag_sum
    else:
        # can't we just add and then extract the real/imag components?
        new_number = full_extracted_set[l_idx, r_idx, 1] + full_extracted_set[l_idx, r_idx, 2]
        summed_set_real[l_idx, r_idx] = new_number.real
        summed_set_imag[l_idx, r_idx] = new_number.imag

    # how does this work??
    append_J[l_idx, r_idx] = complex(0, summed_set_imag[l_idx, r_idx])

    # return full_extracted_set, summed_set_real, summed_set_imag, append_J
    return [summed_set_real, summed_set_imag]


def mctdh(**kwargs):
    """ description of function """

    try:  # remove the previous file, as we are about to write to it
        subprocess_run_wrapper(['rm', '-f', 'mctdh.op'])
    except Exception as e:
        print(f"Error deleting {'mctdh.op'}: {str(e)}")

    # -------------------------------------------------------------------------
    # extract necessary parameters
    nstate, nmodes = A, N

    ndim = kwargs.get('ndim')  # shouldn't need this anymore?
    freqcm = kwargs.get('freqcm')
    nrmmod = kwargs.get('nrmmod')

    s_idx = [i-1 for i in selected_mode_list]  # need to shift by 1 for 0-indexed arrays
    freq_array = freqcm[s_idx]  # assume we return the array
    # mode_array = nrmmod[:, s_idx]  # assume we return the array

    dipoles = kwargs.get('dipoles')
    diabatize = kwargs.get('diabatize')
    hessout = kwargs.get('hessout')

    # constants
    qsize = kwargs.get('qsize', 0.05)
    ha2ev = kwargs.get('ha2ev', 27.2113961318)
    wn2ev = kwargs.get('wn2ev', 0.000123981)

    # -------------------------------------------------------------------------
    # prepare various lists and dictionaries

    # prepare filenames
    zeroth_filename = f'{file_name}_refG.out'

    def _make_displacement_filenames():

        N = pp.N  # make explicit for linter
        linear_displacement_filenames = {}
        bilinear_displacement_filenames = {}

        for i in range(N):
            q1_label = pp.mode_map_dict[i]
            for key in linear_disp_keys:
                d1 = key[0]  # select plus or minus?
                suffix = {'1': '', '2': 'x2'}[key[1]]  # number of +1/+2/-1/-2 maps to suffix
                linear_displacement_filenames[(key, i)] = str(
                    f'{file_name}_mode{q1_label}_{d1}{qsize}{suffix}.out'
                )

        for i, j in it.combinations(range(N), 2):
            q1_label, q2_label = pp.mode_map_dict[i], pp.mode_map_dict[j]
            for d1, d2 in it.product(['+', '-'], ['+', '-']):

                key = d1+d2  # just in case
                assert key in bi_linear_disp_keys, f"{d1+d2=} is not in {bi_linear_disp_keys=}"

                bilinear_displacement_filenames[(key, i, j)] = str(
                   f'{file_name}_mode{q1_label}_{d1}{qsize}_mode{q2_label}_{d2}{qsize}.out'
                )
        return linear_displacement_filenames, bilinear_displacement_filenames

    linear_displacement_filenames, bilinear_displacement_filenames = _make_displacement_filenames()

    # -------------------------------------------------------------------------
    # bad practice, works for now and we can refactor once we've finished figuring out the end product
    format_string = "{label:<25s}={value:>-15.6f}{units:>8s}\n"
    make_line = functools.partial(format_string.format, units=", ev")

    def build_op_section(job_title):
        """Returns a string which defines the `OP_DEFINE-SECTION` of an .op file"""
        start_op, end_op = "OP_DEFINE-SECTION", "end-op_define-section"
        start_title, end_title = "title", "end-title"
        return '\n'.join([
            start_op,
            start_title,
            f"{job_title:s}",
            end_title,
            end_op,
        ])

    # ---------------------------------------------

    def build_frequencies(frequencies):
        """Return a string containing the frequency information of a .op file."""
        return ''.join([
            make_line(label=f"w{i+1:0>2d}", value=w)
            # make_line(label=f"w{mode_map_dict[i]:0>2d}", value=w)
            for i, w in enumerate(frequencies)
        ])

    def build_E0(E0_array):
        """Return a string containing the energy information of a .op file."""

        def one_block_column_aligned():
            # iterates over columns first
            return ''.join([
                make_line(label=f"EH_s{row+1:0>2d}_s{col+1:0>2d}", value=E0_array[row, col])
                for col in range(A)
                for row in range(0, col+1)
                # if not np.isclose(E0_array[row, col], 0.0)  # if you don't want to print the zeros;
            ])

        def one_block_row_aligned():
            # iterates over rows first
            return ''.join([
                make_line(label=f"EH_s{row+1:0>2d}_s{col+1:0>2d}", value=E0_array[row, col])
                for row in range(A)
                for col in range(row, A)
                # if not np.isclose(E0_array[row, col], 0.0)  # if you don't want to print the zeros;
            ])

        def two_blocks():
            diag_block = ''.join([
                make_line(label=f"EH_s{a+1:0>2d}_s{a+1:0>2d}", value=E0_array[a, a])
                for a in range(A)
                # if not np.isclose(E0_array[row, col], 0.0)  # if you don't want to print the zeros;
            ])
            # iterates over columns first (upper triangle)
            off_diag_block = ''.join([
                make_line(label=f"EH_s{row+1:0>2d}_s{col+1:0>2d}", value=E0_array[row, col])
                for row, col in it.combinations(range(A), 2)
                # if not np.isclose(E0_array[row, col], 0.0)  # if you don't want to print the zeros;
            ])
            return diag_block + "\n" + off_diag_block

        # many possible methods
        if False:  # to compare methods turn to True
            s = one_block_column_aligned()
            print('one block, column aligned', s, sep='\n')

            s = one_block_row_aligned()
            print('one block, row aligned', s, sep='\n')

            # I think two blocks looks better personally
            s = two_blocks()
            print('two blocks', s, sep='\n')

        return two_blocks()

        if False:
            return ''.join([
                make_line(label=f"EH_s{row+1:0>2d}_s{col+1:0>2d}", value=E0_array[row, col])
                for idx, E in enumerate(E0_array)
                for col in range(A)
                for row in range(col+1)
                # if not np.isclose(E0_array[row, col], 0.0)  # if you don't want to print the zeros;
            ])

    def build_electronic_moments(model):
        """Return a string containing the electronic transition dipole moment information of a .op file."""

        # this is how many dimensions you have in your co-ordinate system to describe the dipole
        # (i.e. the dimensionality of your dipole's vector, not accounting for the number of surfaces)
        dipole_dimension = int(model[VMK.etdm].shape[0])

        assert 4 > dipole_dimension > 0, "Current hardcoded implementation does not support 4 dimensional dipoles"

        # we ignore the complex warning from casting the etdm (complex) to a float
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", np.ComplexWarning)

            string = ''

            for d in range(dipole_dimension):
                name = {0: 'Ex', 1: 'Ey', 2: 'Ez'}[d]
                string += ''.join([
                    make_line(label=f"{name}_s00_s{a+1:0>2d}", value=float(model[VMK.etdm][d, a]))
                    for a in range(model[VMK.A]-1)
                ])

            return string

    def benny_build_electronic_moments():

        block = ""
        block += hfs.format("ELECTRONIC TRANSITION DIPOLES")
        operate_lst = ["x", "y", "z"]

        for j in range(2, A+1):

            block += "".join([
                f"Ex_s00_s{j:>02d} = {dipoles[j][0]}\n"
                f"Ey_s00_s{j:>02d} = {dipoles[j][1]}\n"
                f"Ez_s00_s{j:>02d} = {dipoles[j][2]}\n"
                "\n"
            ])
            """ if every you need more than 3 dimensions?
            for xyz_idx, op in enumerate(operate_lst):
                block += f"E{x}_s00_s{ist:>02d} = {dipoles[ist][xyz_idx]}\n"
            """
        return block

    def build_magnetic_moments():
        raise Exception("Function not implemented yet")

    def build_linear_coupling(linear_terms, A, N):
        """Return a string containing the linear coupling constant information of a .op file."""
        if diagonal:
            return ''.join([
                make_line(
                    label=f"C1_s{a+1:0>2d}_s{a+1:0>2d}_v{i+1:0>2d}",
                    value=linear_terms[i, a]
                )
                for a, i in it.product(range(A), range(N))
                if not np.isclose(linear_terms[i, a], 0.0)
            ])
        else:
            return '\n'.join([
                ''.join([
                    make_line(
                        label=f"C1_s{a+1:0>2d}_s{a+1:0>2d}_v{i+1:0>2d}",
                        value=linear_terms[i, a, a]
                    )
                    for a, i in it.product(range(A), range(N))
                    if not np.isclose(linear_terms[i, a, a], 0.0)
                ]),
                ''.join([
                    make_line(
                        label=f"C1_s{a2+1:0>2d}_s{a1+1:0>2d}_v{i+1:0>2d}",
                        value=linear_terms[i, a2, a1]
                    )
                    for a1, a2, i in it.product(range(A), range(A), range(N))
                    if a1 != a2 and not np.isclose(linear_terms[i, a2, a1], 0.0)
                ]),
            ])

    """ You would put other functions here
    def build_diagonal_linear_coupling()
    def build_off_diagonal_linear_coupling()

    def build_diagonal_quadratic_coupling()
    def build_off_diagonal_quadratic_coupling()

    def build_diagonal_bilinear_coupling()
    def build_off_diagonal_bilinear_coupling()
    """

    def _SOC_Somehwere():
        # ------------------------------------------------------------
        # some SOC stuff that could be moved elsewhere?
        SOC_check = subprocess_run_wrapper(['grep', "(DIRECT MAX.)", zeroth_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        SOC_flag = SOC_check.returncode == 0

        if SOC_flag:
            # Extract DSOME_cm_0
            DSOME_cm_0 = extract_DSOME(zeroth_filename, nstate)
            DSOME_cm_0_real, DSOME_cm_0_imag = DSOME_cm_0[0], DSOME_cm_0[1]
        # ------------------------------------------------------------
        return

    heading, params = [], []
    EH, SOC, linear, quadratic, bilinear = [], [], [], [], []

    displacement_keys = ["+1", "+2", "-1", "-2"]

    format_string = "{label:<25s}={value:>-15.9f}{units:>8s}\n"
    make_line = functools.partial(format_string.format, units=", ev")

    spacer_format_string = f"# {'-':^60s} #\n"
    hfs = header_format_string = "# {:^60s} #\n" + spacer_format_string

    linear.append(hfs.format('Linear Coupling Constants'))
    quadratic.append(hfs.format('Quadratic Coupling Constants'))
    bilinear.append(hfs.format('Bilinear Coupling Constants'))

    def extract_E_diab_au():
        """ x """

        linear_displacement_filenames  # shape of (key, i)
        bilinear_displacement_filenames  # shape of (key, i, j)

        # strings used by `grep` to locate values to extract
        a_pattern = 'STATE #.* {a}.S GMC-PT-LEVEL DIABATIC ENERGY='
        ba_pattern = 'STATE #.* {b} &.* {a}.S GMC-PT-LEVEL COUPLING'

        # store a 1D array at each key
        E_diab_au = {np.zeros(N) for key in displacement_keys}

        for key in displacement_keys:
            for i in range(N):
                path = linear_displacement_filenames[(key, i)]
                E_diab_au[key][i] = extract_diabatic_energy(path, i_pattern.format(i))

        for i, j in it.combinations(range(N), 2):
            # for key in displacement_keys:
            Coup_ev = refG_extract(zeroth_filename, ji_pattern.format(j=jst, i=ist))

            # Extract Coup_ev_0
            Coup_ev_0 = extract_coupling_energy(zeroth_filename, pattern)

            Coup_ev = {}
            for key in displacement_keys:
                Coup_ev[key] = extract_diabatic_energy(displacement_filenames[key], pattern)

        return E_diab_au

    def extract_E0(hessout, A=pp.A):
        """ The energy values associated with the reference geometry.
        These energies don't depend on the modes
        """

        # strings used by `grep` to locate values to extract
        a_pattern = 'STATE #.* {col}.S GMC-PT-LEVEL DIABATIC ENERGY='
        ba_pattern = 'STATE #.* {row} &.* {col}.S GMC-PT-LEVEL COUPLING'

        # ground state of the (optimized geometry?) (includes fictitious surface? no/yes?)
        GSE_pattern = 'TOTAL ENERGY ='
        GSE = extract_ground_state_energy(hessout, GSE_pattern)

        # 1st diabat's energy
        D1E = extract_diabatic_energy(zeroth_filename, a_pattern.format(col=1))
        linear_shift = (D1E - GSE) * ha2ev
        print(
            f'The ground state energy is: {GSE} Hartree\n'
            f'Diabat #1 energy is: {D1E} Hartree\n'
            f'Linear shift value: {linear_shift} eV\n'
        )

        # A = pp.A + 1  # +1 for ground state's as fictitious first row/col
        A = pp.A  #?
        shape = (A, A)
        E0_array = np.zeros(shape)

        if False:  # debug
            def _reminder_produce_upper_triangle_indices(s):
                """ Assume a square matrix of shape (s, s).
                Prints the row, col indices that would index the upper triangle
                Just reminder code.
                """
                print("Indices using it.combinations")
                for row, col in it.combinations(range(s), 2):
                    print(f"{row=} {col=}")

                print("Indices using nested loops")
                for col in range(s):
                    for row in range(0, col):
                        print(f"{row=} {col=}")

                print("Indices assuming we start at 1")
                # assume indices start at 1
                for col in range(1, s+1):
                    for row in range(1, col):
                        print(f"{row=} {col=}")
                return

            _reminder_produce_upper_triangle_indices(A)
            breakpoint()

        if False:  # matches the bash style from Toby
            for a in range(A):
                E0_array[a, a] = refG_extract(zeroth_filename, a_pattern.format(col=a+1))
                E0_array[a, a] += linear_shift
                for b in range(a):
                    E0_array[b, a] = refG_extract(zeroth_filename, ba_pattern.format(row=b+1, col=a+1))
            # print(E0_array)

        if True:  # this makes more sense based on what I see in the file
            for a in range(A):
                E0_array[a, a] = refG_extract(zeroth_filename, a_pattern.format(col=a+1))
                E0_array[a, a] += linear_shift

            for b, a in it.combinations(range(A), 2):
                E0_array[b, a] = refG_extract(zeroth_filename, ba_pattern.format(row=b+1, col=a+1))
            # print(E0_array)

        if True:  # debug outputs (keep them out of the loops above for simplicity)
            for a, b in it.combinations(range(A), 2):
                print(f"Diabatic energy at state {a}: {E0_array[a, a]}")
                print(f"Coupling energy at state {b} & {a}: {E0_array[a, b]}")
            print(E0_array)

        return E0_array

    def extract_energies_benny():

        for kmode in range(N):  # kmode is array index of selected_mode_list
            imode = pp.mode_map_dict[kmode] #

            displacement_filenames = {
                "+1": f'{filnam}_mode{imode}_+{qsize}.out',
                "+2": f'{filnam}_mode{imode}_+{qsize}x2.out',
                "-1": f'{filnam}_mode{imode}_-{qsize}.out',
                "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
            }

            vibron_ev = freqcm[imode-1] * wn2ev
            params.append(make_line(label=f"w{imode:>02d}", value=vibron_ev))
            # params.append("\n")
            # Coupling.append("#Linear and quadratic diagonal and off-diagonal vibronic coupling constants:\n")

            grace_code = {}
            for key in displacement_keys:
                grace_code[key] = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", displacement_filenames[key]])

            """ either of these work (logic wise)
                if any(code != 0 for code in grace_code.values()):
                if not all(code == 0 for code in grace_code.values()):
            """

            if not all(code == 0 for code in grace_code.values()):
                params.append(f"not good to extract. Skipping mode {imode} for extracting vibronic couplings\n")

            else:  # otherwise we're good to extract
                print("\n good to extract\n")
                # Extract the diagonal and off-diagonal vibronic coupling
                for ist in range(1, nstate + 1):

                    def _make_diag_lin_quad(i):
                        pattern = f'STATE #.* {i}.S GMC-PT-LEVEL DIABATIC ENERGY='

                        # Ediab_au = [extract_diabatic_energy(displacement_filenames[kmode][k], pattern) for k in displacement keys]  $ one liner list comprehension
                        # Ediab_au = {k: extract_diabatic_energy(displacement_filenames[kmode][k], pattern) for k in displacement keys}  # one liner dictionary comprehension
                        Ediab_au = {}
                        for key in displacement_keys:
                            Ediab_au[key] = extract_diabatic_energy(displacement_filenames[key], pattern)

                        # Extract Ediab_au_0
                        Ediab_au_0 = extract_diabatic_energy(zeroth_filename, pattern)
                        linear_diag_ev = (Ediab_au["+1"] - Ediab_au["-1"]) * ha2ev / (2 * qsize)
                        quadratic_diag_ev = (Ediab_au["+2"] + Ediab_au["-2"] - 2.0 * Ediab_au_0) * ha2ev / (4.0 * qsize * qsize)

                        # We only view the difference between the actual force constant and the vibron
                        # as the quadratic diagonal coupling for the diabatic state.
                        quadratic_diag_ev = quadratic_diag_ev - vibron_ev

                        # Print and store results
                        print(f"State {i} Linear Diagonal: {linear_diag_ev} Quadratic Diagonal: {quadratic_diag_ev}, ev\n")

                        # machine accuracy is typically 16 digits
                        s1 = make_line(label=f"C1_s{i:>02d}_s{i:>02d}_v{imode:>02d}", value=linear_diag_ev)
                        s2 = make_line(label=f"C2_s{i:>02d}s{i:>02d}_v{imode:>02d}v{imode:>02d}", value=quadratic_diag_ev)
                        return s1, s2

                    s1, s2 = _make_diag_lin_quad(ist)
                    linear.append(s1)
                    quadratic.append(s2)

                    # # Loop over jst
                    jlast = ist - 1
                    for jst in range(1, jlast + 1):

                        def _make_offdiag_lin_quad(i, j):
                            pattern = f'STATE #.* {j} &.* {i}.S GMC-PT-LEVEL COUPLING'

                            # Extract Coup_ev_0
                            Coup_ev_0 = extract_coupling_energy(zeroth_filename, pattern)

                            Coup_ev = {}
                            for key in displacement_keys:
                                Coup_ev[key] = extract_diabatic_energy(displacement_filenames[key], pattern)

                            # Compute linear off-diagonal coupling
                            linear_offdiag_ev = (Coup_ev["+1"] - Coup_ev["-1"]) / (2 * qsize)
                            # Compute quadratic off-diagonal coupling
                            quadratic_offdiag_ev = (Coup_ev["+2"] + Coup_ev["-2"] - 2.0 * Coup_ev_0) / (4.0 * qsize * qsize)

                            # Print and store results
                            print(f"State {j} & {i} Linear Off-Diagonal: {linear_offdiag_ev}\n")
                            print(f"State {j} & {i} Quadratic Off-Diagonal: {quadratic_offdiag_ev}\n")
                            s1 = make_line(label=f"C1_s{j:>02d}_s{i:>02d}_v{imode:>02d}", value=linear_offdiag_ev)
                            s2 = make_line(label=f"C2_s{j:>02d}s{i:>02d}_v{imode:>02d}v{imode:>02d}", value=quadratic_offdiag_ev)
                            return s1, s2

                        s1, s2 = _make_diag_lin_quad(ist)
                        linear.append(s1)
                        quadratic.append(s2)

                        """ this is just representative (you can delete - just for learning purposes)
                        if False: # don't actually try to do right now
                            order_name = {1: 'Linear', 2: 'Quadratic'}
                            for i in [1, 2]:
                                _number = [linear_offdiag_ev, quadratic_offdiag_ev][i]
                                print(f"State {jst} & {ist} {order_name[i]} Off-Diagonal: {_number}\n")
                                oprder_list[i].append(make_line(label=f"C{i}_s{jst:>02d}_s{ist:>02d}_v{imode:>02d}", value=_number))
                        """


            # Extracting bilinear vibronic coupling
            # Coupling.append("#Bilinear diagonal and off-diagonal vibronic coupling constants:\n")

            for lmode in range(0, kmode):
                jmode = pp.mode_map_dict[lmode]
                breakpoint()

                bi_linear_displacement_filenames = {
                    "++": f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out',
                    "+-": f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out',
                    "-+": f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out',
                    "--": f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out',
                }

                grace_code_pp = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out"])
                grace_code_pm = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out"])
                grace_code_mp = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out"])
                grace_code_mm = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out"])

                if all(code == 0 for code in [grace_code_pp, grace_code_pm, grace_code_mp, grace_code_mm]):
                    print(f"\n Good to extract bilinear for modes {imode} {jmode} \n")
                    for ist in range(1, nstate + 1):
                        pattern = f'STATE #.* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY='

                        # do this style again?
                        # big_displacement_keys = ['++', '+-', '-+', '--']
                        # Ediab_au = {}
                        # for key in displacement_keys:
                        #     Ediab_au[key] = extract_diabatic_energy(big_displacement_filenames[key], pattern)


                        # Extract Ediab_au_pp
                        Ediab_au_pp = extract_diabatic_energy(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out', pattern)

                        # Extract Ediab_au_pm
                        Ediab_au_pm = extract_diabatic_energy(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out', pattern)

                        # Extract Ediab_au_mp
                        Ediab_au_mp = extract_diabatic_energy(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out', pattern)

                        # Extract Ediab_au_mm
                        Ediab_au_mm = extract_diabatic_energy(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out', pattern)

                        bilinear_diag_ev = (Ediab_au_pp + Ediab_au_mm - Ediab_au_pm - Ediab_au_mp ) * ha2ev / (4.0 * qsize * qsize )

                        print(f"State {ist} Bilinear Diagonal: {bilinear_diag_ev}\n")
                        bilinear.append(make_line(label=f"C1_s{ist:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}", value=bilinear_diag_ev))

                        # # Loop over jst
                        jlast = ist - 1
                        for jst in range(1, jlast + 1):
                            pattern = f'STATE #.* {jst} &.* {ist}.S GMC-PT-LEVEL COUPLING'
                            # Extract Coup_ev_pp
                            Coup_ev_pp = extract_coupling_energy(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out', pattern)
                            # Extract Coup_ev_pm
                            Coup_ev_pm = extract_coupling_energy(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out', pattern)

                            # Extract Coup_ev_mp
                            Coup_ev_mp = extract_coupling_energy(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out', pattern)

                            # Extract Coup_ev_mm
                            Coup_ev_mm = extract_coupling_energy(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out', pattern)

                            bilinear_offdiag_ev = ( Coup_ev_pp + Coup_ev_mm - Coup_ev_pm - Coup_ev_mp ) / (4.0 * qsize * qsize )

                            print(f"State {jst} & {ist} Bilinear Off-Diagonal: {bilinear_offdiag_ev}\n")
                            bilinear.append(make_line(label=f"C1_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}", value=bilinear_offdiag_ev))

                            if SOC_flag:

                                try:

                                    # Extract DSOME_cm_pp
                                    DSOME_cm_pp = extract_DSOME(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out', nstate)
                                    DSOME_cm_pp_real, DSOME_cm_pp_imag = DSOME_cm_pp[0], DSOME_cm_pp[1]

                                    # Extract DSOME_cm_pm
                                    DSOME_cm_pm = extract_DSOME(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out', nstate)
                                    DSOME_cm_pm_real, DSOME_cm_pm_imag = DSOME_cm_pm[0], DSOME_cm_pm[1]

                                    # Extract DSOME_cm_mp
                                    DSOME_cm_mp = extract_DSOME(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out', nstate)
                                    DSOME_cm_mp_real, DSOME_cm_mp_imag = DSOME_cm_mp[0], DSOME_cm_mp[1]

                                    # Extract DSOME_cm_mm
                                    DSOME_cm_mm = extract_DSOME(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out', nstate)
                                    DSOME_cm_mm_real, DSOME_cm_mm_imag = DSOME_cm_mm[0], DSOME_cm_mm[1]


                                except Exception as e:
                                    print(f"Error in SOC: {str(e)}")
                else:
                    print(f"not good to extract. Skipping mode {imode} mode {jmode} for extracting bilinear vibronic couplings")

    def extract_etdm():
        return

    def extract_linear():
        return

    def extract_quadratic():
        return

    def extract_bilinear():
        return

    def build_parameter_section(nof_states, nof_modes):
        """Returns a string which defines the `PARAMETER-SECTION` of an .op file"""
        start, end = "PARAMETER-SECTION", "end-parameter-section"
        header_string = "#{0:^47}#\n#{spacer:^47}#\n"
        make_header = functools.partial(header_string.format, spacer='-' * 45)
        # ----------------------------------------------------------
        # read in all the necessary parameters
        vibron_ev = freq_array * wn2ev  # fix later (uses global variables)
        E_diab_au = extract_E0(hessout)
        # extract_energies()
        E_moments = extract_etdm()
        # M_moments = extract_mtdm()
        lin_data = extract_linear()
        # quad_data = extract_quadratic()
        # bilin_data = extract_bilinear()

        # ----------------------------------------------------------
        # begin to assemble the parameter section
        return_list = [
            start,
            make_header('Frequencies'), build_frequencies(vibron_ev),
            make_header('Electronic Hamitonian'), build_E0(E_diab_au),
        ]
        print('\n'.join(return_list))
        breakpoint()  # works up till here

        return_list.extend([
            # make_header('Electronic transition moments'), build_electronic_moments(E_moments),
            # # make_header('Magnetic transition moments'), build_magnetic_moments(M_moments),
            make_header('Linear Coupling Constants'),  #build_linear_coupling(lin_data, nof_states, nof_modes),
        ])

        if False:  # turn on once finished functions
            return_list.extend([
                make_header('Diagonal Quadratic Coupling Constants'),
                # build_diagonal_quadratic_coupling(quad_data, A, N, flag_diagonal),
                make_header('Off_diagonal Quadratic Coupling Constants'),
                # build_off_diagonal_quadratic_coupling(quad_data, A, N, flag_diagonal),
                #
                make_header('Diagonal Bilinear Coupling Constants'),
                # build_diagonal_bilinear_coupling(bilin_data, A, N, flag_diagonal),
                make_header('Off_diagonal Bilinear Coupling Constants'),
                # build_off_diagonal_bilinear_coupling(bilin_data, A, N, flag_diagonal),
            ])

        if False:  # SOC_flag
            return_list.extend([
                make_header('SOC'),
                build_linear_coupling(X_lin, nof_states, nof_modes),
            ])

        # return_list.extend(linear)
        # return_list.extend(quadratic)
        # return_list.extend(bilinear)
        return_list.append(end)

        return '\n'.join(return_list)
    # # ----------------------------------------------------------

    def label_momentum(N):
        """Return a string containing the momentum labelling of a .op file."""
        spacer = '|'
        return '\n'.join([
            f"1.00*w{i:0>2d}{spacer:>12}{i+1:<3d}KE"
            for i in range(1, N+1)
        ]) + '\n'

    def label_position(N):
        """Return a string containing the position labelling of a .op file."""
        spacer = '|'
        return '\n'.join([
            f"0.50*w{i:0>2d}{spacer:>12}{i+1:<3d}q^2"
            for i in range(1, N+1)
        ]) + '\n'

    def label_energies(A):
        """Return a string containing the energy labelling of a .op file."""
        spacer = '|'
        return '\n'.join([
            f"EH_s{a:0>2d}_s{a:0>2d}{spacer:>15}1 S{a:d}&{a:d}"
            for a in range(1, A+1)
        ]) + '\n'

    def label_linear_coupling(linear_terms, A, N):
        """Return a string containing the linear coupling constant labelling of a .op file."""
        spacer = '|'
        if diagonal:
            return '\n'.join([
                f"C1_s{a:0>2d}_s{a:0>2d}_v{i:0>2d}{spacer:>11}1 S{a:d}&{a:d}{spacer:>4}{i+1}  q"
                for a, i in it.product(range(1, A+1), range(1, N+1))
                if not np.isclose(linear_terms[i-1, a-1], 0.0)
            ]) + '\n'
        else:
            return '\n'.join([
                f"C1_s{a:0>2d}_s{a:0>2d}_v{i:0>2d}{spacer:>11}1 S{a:d}&{a:d}{spacer:>4}{i+1}  q"
                for a, i in it.product(range(1, A+1), range(1, N+1))
                if not np.isclose(linear_terms[i-1, a-1, a-1], 0.0)
            ] + [
                ''  # creates a blank line between the (surface) diagonal and off-diagaonl linear terms
            ] + [
                f"C1_s{a2:0>2d}_s{a1:0>2d}_v{i:0>2d}{spacer:>11}1 S{a2:d}&{a1:d}{spacer:>4}{i+1}  q"
                for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
                if (a1 != a2) and not np.isclose(linear_terms[i-1, a2-1, a1-1], 0.0)
            ]) + '\n'

    """ You would put other functions here
    def label_diagonal_linear_coupling()
    def label_off_diagonal_linear_coupling()

    def label_diagonal_quadratic_coupling()
    def label_off_diagonal_quadratic_coupling()

    def label_diagonal_bilinear_coupling()
    def label_off_diagonal_bilinear_coupling()
    """

    def benny_label_SOC_section():

        spacer_format_string = f"# {'-':^60s} #\n"
        hfs = header_format_string = "# {:^60s} #\n" + spacer_format_string
        block = hfs.format("SOC FULL HAMILTONIAN SOC OFF-DIAGONAL VIBRONIC COUPLINGS")

        # prepare `make_line`
        format_string_1 = "{label:<25s}{link:<20s}\n"
        format_string_2 = "{label:<25s}{link:<20s}\n"
        format_string_3 = "{label:<25s}{link:<20s}\n"
        format_string_4 = "{label:<25s}{link:<20s}\n"

        make_line_1 = functools.partial(format_string_1.format)
        make_line_2 = functools.partial(format_string_2.format)
        make_line_3 = functools.partial(format_string_3.format)
        make_line_4 = functools.partial(format_string_4.format)

        for i, a in it.product(range(1,N+1), range(1,A+1)):
            for j, b in it.product(range(1,i+1), range(1,a+1)):
                i_label, j_label = modes_included[[i, j]]
                print(f"{i=}, {j=}, {a=}, {b=}")

                l1 = "C1_s{:>02d}_s{:>02d}_v{:>02d}".formmat(j, i, a)
                make_line(label=f"I*{l1}r", link=f"|1 Z{b}&{a} | {j+1} q")
                make_line(label=f"-I*{l1}i", link=f"|1 Z{a}&{b} | {j+1} q")

                l2 = "C2_s{:>02d}_s{:>02d}_v{:>02d}".formmat(j, i, a)
                make_line(label=f"I*{l2}r", link=f"|1 Z{b}&{a} | {j+1} q^2")
                make_line(label=f"-I*{l2}i", link=f"|1 Z{a}&{b} | {j+1} q^2")

                l3 = "C1_s{:>02d}_s{:>02d}_v{:>02d}_v{:>02d}".formmat(j, i, i_label, j_label)
                make_line(label=f"I*{l3}r", link=f"|1 Z{b}&{a} | {j+1} q | {i+1} q")
                make_line(label=f"-I*{l3}i", link=f"|1 Z{a}&{b} | {j+1} q | {i+1} q")

                l4 = "SOC_s{:>02d}_s{:>02d}_v{:>02d}_v{:>02d}".formmat(j, i, i_label, j_label)
                make_line(label=f"I*{l4}r", link=f"|1 Z{b}&{a} | {j+1} q")
                make_line(label=f"-I*{l4}i", link=f"|1 Z{a}&{b} | {j+1} q")
        return

    def build_hamiltonian_section(nof_states, nof_modes):
        """Returns a string which defines the `HAMILTONIAN-SECTION` of an .op file"""
        start, end = "HAMILTONIAN-SECTION", "end-hamiltonian-section"
        spec = ''.join([
            ' modes   |  el  |',
            ''.join([f" v{N+1:0>2d}|" for N in range(nof_modes)]),
            '\n'
        ])

        return_list = [
            start,
            spec,
            label_momentum(nof_modes),
            label_position(nof_modes),
            label_energies(nof_states),
            label_linear_coupling(_values_go_here, nof_states, nof_modes),
        ]

        # if highest_order > 1:
        #     return_list.append(label_diagonal_quadratic_coupling(model[VMK.G2], A, N))
        #     return_list.append(label_off_diagonal_quadratic_coupling(model[VMK.G2], A, N))

        # if highest_order > 2:
        #     return_list.append(label_diagonal_cubic_coupling(model[VMK.G3], A, N))

        # if highest_order > 3:
        #     return_list.append(label_diagonal_quartic_coupling(model[VMK.G4], A, N))

        if SOC_flag:
            return_list.append(benny_label_SOC_section())

        return_list.append(end)

        return '\n'.join(return_list)
    # # ----------------------------------------------------------

    def build_dipole_moments_section(nof_states, nof_modes):
        """Returns a string which defines the `HAMILTONIAN-SECTION_Ex` of an .op file"""
        start, end = "HAMILTONIAN-SECTION_Ex", "end-hamiltonian-section"
        spec = ''.join([
            ' modes   |  el  |',
            ''.join([f" v{n+1:0>2d}|" for n in range(nof_modes)]),
            '\n'
        ])

        def label_dipole_moments_energies(A):
            """Return a string containing the dipole moments energy labelling of a .op file."""
            spacer = '|'
            return '\n'.join([
                f"Ex_s00_s{a:0>2d}{spacer:>15}1 S{A:d}&{a:d}"
                for a in range(1, A)
            ]) + '\n'

        string = '\n'.join([
            start,
            spec,
            label_dipole_moments_energies(nof_states),
            end,
        ])
        return string

    def benny_build_dipole_moments_section():
        """ x """
        block = f"HAMILTONIAN-SECTION_Ex\n\n"

        # Write modes and mode labels
        mode_number_key = [modes_included[i] for i in range(N)]
        h_labels = ["modes", "el", ] + [
            f"v{s:>02d}"
            for s in mode_number_key
        ]

        block += " | ".join(h_labels) + "\n"
        block += f"{'-':47}\n"

        for j in range(N):
            block += f"1.0         |1 S{N+1}&{j+1}\n"

        block += "\n\nend-hamiltonian-section\n\n"
        return block

    # ----------------------------------------------------------

    def _write_op():

        job_title = f'{filnam} {A} states + ' + str(N) + ' modes'

        file_contents = "\n".join([
            build_op_section(job_title),
            build_parameter_section(A, N),
            #build_hamiltonian_section(A, N),
            #build_dipole_moments_section(A, N),
            "end-operator\n"
        ])

        with open('mctdh.op', 'w') as fp:
            fp.write(file_contents)

    def confirm_necessary_files_exist():

        def _confirm_linear_are_good(bad_mode=False):
            """ confirm all displacement_input_files are good for all selected modes """
            for i in range(N):
                grace_code = {}
                for key in linear_disp_keys:
                    grace_code[key] = subprocess_call_wrapper([
                        "grep", "DONE WITH MP2 ENERGY",
                        linear_displacement_filenames[(key, i)]
                    ])
                    print(f" ..... in file {linear_displacement_filenames[(key, i)]}")

                if not all(code == 0 for code in grace_code.values()):
                    mode_label = pp.mode_map_dict[i]
                    print("Linear/Quad mode {mode_label} not good to extract.\n")
                    bad_mode = True

            return bad_mode

        def _confirm_bilinear_are_good(bad_mode=False):
            """ confirm all displacement_input_files are good for all selected modes """
            for i, j in it.combinations(range(N), 2):
                print(i, j)
                grace_code = {}
                for key in bi_linear_disp_keys:
                    grace_code[key] = subprocess_call_wrapper([
                        "grep", "DONE WITH MP2 ENERGY",
                        bilinear_displacement_filenames[(key, i, j)]
                    ])
                    print(f" ..... in file {bilinear_displacement_filenames[(key, i, j)]}")

                if not all(code == 0 for code in grace_code.values()):
                    mode_label = (pp.mode_map_dict[i], pp.mode_map_dict[j])
                    print("Bilinear mode {mode_label} not good to extract.\n")
                    bad_mode = True

            return bad_mode

        def refG_file_exists():
            refG_exists = bool(subprocess_call_wrapper(["ls", zeroth_filename]) == 0)
            if not refG_exists:
                """ If refG doesn't exist, then maybe we can find the .... from other files `blah.txt`
                but also maybe they don't exist... need to check later
                """
                print(f"Skip extracting Hamiltonians from the non-existing {zeroth_filename}")
                breakpoint()
                raise Exception("need to add failback code if can't find refG")

            return refG_exists

        flag = bool(
            _confirm_linear_are_good()
            or _confirm_bilinear_are_good()
            or refG_file_exists()
        )
        return flag
    # ------------------------------------------------------------------------

    # execution starts here!!
    if not confirm_necessary_files_exist():
        print("Bad input files detected, please fix.")
        import sys; sys.exit()

    _write_op()

    breakpoint()

    # ------------------------------------------------------------------------
    import sys; sys.exit()
    # ------------------------------------------------------------------------

    def old_code_block():
        # different header style
        hfs = header_format_string = spacer_format_string + "# {:^60s} #\n" + spacer_format_string

        file_contents += _make_ETD_block(header_format_string)
        file_contents += "end-parameter-section\n"

        # ----------------------------------------------------------

        format_string = "{label:<25s}={value:>-15.9f}{units:>8s}\n"
        make_line = functools.partial(format_string.format, units=", ev")

        distcoord_plus, distcoord_minus, distcoord_plus_x2, distcoord_minus_x2 = diabatize[0], diabatize[1], diabatize[2], diabatize[3]
        distcoord_pp, distcoord_pm, distcoord_mp, distcoord_mm = diabatize[4], diabatize[5], diabatize[6], diabatize[7]

        coord_disp_plus = {}
        coord_disp_minus = {}
        coord_disp_plusx2 = {}
        coord_disp_minusx2 = {}
        coord_disp_pp = {}
        coord_disp_pm = {}
        coord_disp_mp = {}
        coord_disp_mm = {}

        for icomp in range(1, ndim + 1):
            coord_disp_plus[icomp] = distcoord_plus[icomp]
            print(f'icomp: {icomp}, coord_disp_plus: {coord_disp_plus[icomp]}')
            coord_disp_minus[icomp] = distcoord_minus[icomp]
            print(f'icomp: {icomp}, coord_disp_minus: {coord_disp_minus[icomp]}')
            coord_disp_plusx2[icomp] = distcoord_plus_x2[icomp]
            print(f'icomp: {icomp}, coord_disp_plusx2: {coord_disp_plusx2[icomp]}')
            coord_disp_minusx2[icomp] = distcoord_minus_x2[icomp]
            print(f'icomp: {icomp}, coord_disp_minusx2: {coord_disp_minusx2[icomp]}')
            coord_disp_pp[icomp] = distcoord_pp[icomp]
            print(f'icomp: {icomp}, coord_disp_pp: {coord_disp_pp[icomp]}')
            coord_disp_pm[icomp] = distcoord_pm[icomp]
            print(f'icomp: {icomp}, coord_disp_pm: {coord_disp_pm[icomp]}')
            coord_disp_mp[icomp] = distcoord_mp[icomp]
            print(f'icomp: {icomp}, coord_disp_mp: {coord_disp_mp[icomp]}')
            coord_disp_mm[icomp] = distcoord_mm[icomp]
            print(f'icomp: {icomp}, coord_disp_mm: {coord_disp_mm[icomp]}')

        spacer_format_string = f"# {'-':^60s} #\n"
        hfs = header_format_string = "# {:^60s} #\n" + spacer_format_string

        params.append(hfs.format('Frequencies'))
        linear.append(hfs.format('Linear Coupling Constants'))
        quadratic.append(hfs.format('Quadratic Coupling Constants'))
        bilinear.append(hfs.format('Bilinear Coupling Constants'))

        # Loop through modes
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]

            displacement_filenames = {
                "+1": f'{filnam}_mode{imode}_+{qsize}.out',
                "+2": f'{filnam}_mode{imode}_+{qsize}x2.out',
                "-1": f'{filnam}_mode{imode}_-{qsize}.out',
                "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
            }

            vibron_ev = freqcm[imode] * wn2ev
            Params.append(make_line(label=f"w{imode:>02d}", value=vibron_ev))
            # Params.append("\n")
            # Coupling.append("#Linear and quadratic diagonal and off-diagonal vibronic coupling constants:\n")

            grace_code = {}
            for key in displacement_keys:
                grace_code[key] = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", displacement_filenames[key]])

            """ either of these work (logic wise)
                if any(code != 0 for code in grace_code.values()):
                if not all(code == 0 for code in grace_code.values()):
            """

            if not all(code == 0 for code in grace_code.values()):
                print(f"not good to extract. Skipping mode {imode} for extracting vibronic couplings\n")

            else:  # otherwise we're good to extract
                print("\n good to extract\n")
                # Extract the diagonal and off-diagonal vibronic coupling
                for ist in range(1, nstate + 1):

                    def _make_diag_lin_quad(i):
                        pattern = f'STATE #.* {i}.S GMC-PT-LEVEL DIABATIC ENERGY='

                        # Ediab_au = [extract_diabatic_energy(displacement_filenames[kmode][k], pattern) for k in displacement keys]  $ one liner list comprehension
                        # Ediab_au = {k: extract_diabatic_energy(displacement_filenames[kmode][k], pattern) for k in displacement keys}  # one liner dictionary comprehension
                        Ediab_au = {}
                        for key in displacement_keys:
                            Ediab_au[key] = extract_diabatic_energy(displacement_filenames[key], pattern)

                        # Extract Ediab_au_0
                        Ediab_au_0 = extract_diabatic_energy(zeroth_filename, pattern)
                        linear_diag_ev = (Ediab_au["+1"] - Ediab_au["-1"]) * ha2ev / (2 * qsize)
                        quadratic_diag_ev = (Ediab_au["+2"] + Ediab_au["-2"] - 2.0 * Ediab_au_0) * ha2ev / (4.0 * qsize * qsize)

                        # We only view the difference between the actual force constant and the vibron
                        # as the quadratic diagonal coupling for the diabatic state.
                        quadratic_diag_ev = quadratic_diag_ev - vibron_ev

                        # Print and store results
                        print(f"State {i} Linear Diagonal: {linear_diag_ev} Quadratic Diagonal: {quadratic_diag_ev}, ev\n")

                        # machine accuracy is typically 16 digits
                        s1 = make_line(label=f"C1_s{i:>02d}_s{i:>02d}_v{imode:>02d}", value=linear_diag_ev)
                        s2 = make_line(label=f"C2_s{i:>02d}s{i:>02d}_v{imode:>02d}v{imode:>02d}", value=quadratic_diag_ev)
                        return s1, s2

                    s1, s2 = _make_diag_lin_quad(ist)
                    Linear.append(s1)
                    Quadratic.append(s2)

                    # # Loop over jst
                    jlast = ist - 1
                    for jst in range(1, jlast + 1):

                        def _make_offdiag_lin_quad(i, j):
                            pattern = f'STATE #.* {j} &.* {i}.S GMC-PT-LEVEL COUPLING'

                            # Extract Coup_ev_0
                            Coup_ev_0 = extract_coupling_energy(zeroth_filename, pattern)

                            Coup_ev = {}
                            for key in displacement_keys:
                                Coup_ev[key] = extract_diabatic_energy(displacement_filenames[key], pattern)

                            # Compute linear off-diagonal coupling
                            linear_offdiag_ev = (Coup_ev["+1"] - Coup_ev["-1"]) / (2 * qsize)
                            # Compute quadratic off-diagonal coupling
                            quadratic_offdiag_ev = (Coup_ev["+2"] + Coup_ev["-2"] - 2.0 * Coup_ev_0) / (4.0 * qsize * qsize)

                            # Print and store results
                            print(f"State {j} & {i} Linear Off-Diagonal: {linear_offdiag_ev}\n")
                            print(f"State {j} & {i} Quadratic Off-Diagonal: {quadratic_offdiag_ev}\n")
                            s1 = make_line(label=f"C1_s{j:>02d}_s{i:>02d}_v{imode:>02d}", value=linear_offdiag_ev)
                            s2 = make_line(label=f"C2_s{j:>02d}s{i:>02d}_v{imode:>02d}v{imode:>02d}", value=quadratic_offdiag_ev)
                            return s1, s2

                        s1, s2 = _make_diag_lin_quad(ist, jst)
                        Linear.append(s1)
                        Quadratic.append(s2)

                        """ this is just representative (you can delete - just for learning purposes)
                        if False: # don't actually try to do right now
                            order_name = {1: 'Linear', 2: 'Quadratic'}
                            for i in [1, 2]:
                                _number = [linear_offdiag_ev, quadratic_offdiag_ev][i]
                                print(f"State {jst} & {ist} {order_name[i]} Off-Diagonal: {_number}\n")
                                oprder_list[i].append(make_line(label=f"C{i}_s{jst:>02d}_s{ist:>02d}_v{imode:>02d}", value=_number))
                        """


            # Extracting bilinear vibronic coupling
            # Coupling.append("#Bilinear diagonal and off-diagonal vibronic coupling constants:\n")
            lmode_last = kmode - 1
            for lmode in range(1, lmode_last + 1):
                jmode = modes_included[lmode]

                grace_code_pp = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out"])
                grace_code_pm = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out"])
                grace_code_mp = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out"])
                grace_code_mm = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out"])

                if all(code == 0 for code in [grace_code_pp, grace_code_pm, grace_code_mp, grace_code_mm]):
                    print(f"\n Good to extract bilinear for modes {imode} {jmode} \n")
                    for ist in range(1, nstate + 1):
                        pattern = f'STATE #.* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY='

                        # do this style again?
                        # big_displacement_keys = ['++', '+-', '-+', '--']
                        # Ediab_au = {}
                        # for key in displacement_keys:
                        #     Ediab_au[key] = extract_diabatic_energy(big_displacement_filenames[key], pattern)


                        # Extract Ediab_au_pp
                        Ediab_au_pp = extract_diabatic_energy(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out', pattern)

                        # Extract Ediab_au_pm
                        Ediab_au_pm = extract_diabatic_energy(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out', pattern)

                        # Extract Ediab_au_mp
                        Ediab_au_mp = extract_diabatic_energy(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out', pattern)

                        # Extract Ediab_au_mm
                        Ediab_au_mm = extract_diabatic_energy(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out', pattern)

                        bilinear_diag_ev = (Ediab_au_pp + Ediab_au_mm - Ediab_au_pm - Ediab_au_mp ) * ha2ev / (4.0 * qsize * qsize )

                        print(f"State {ist} Bilinear Diagonal: {bilinear_diag_ev}\n")
                        Bilinear.append(make_line(label=f"C1_s{ist:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}", value=bilinear_diag_ev))

                        # # Loop over jst
                        jlast = ist - 1
                        for jst in range(1, jlast + 1):
                            pattern = f'STATE #.* {jst} &.* {ist}.S GMC-PT-LEVEL COUPLING'
                            # Extract Coup_ev_pp
                            Coup_ev_pp = extract_coupling_energy(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out', pattern)
                            # Extract Coup_ev_pm
                            Coup_ev_pm = extract_coupling_energy(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out', pattern)

                            # Extract Coup_ev_mp
                            Coup_ev_mp = extract_coupling_energy(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out', pattern)

                            # Extract Coup_ev_mm
                            Coup_ev_mm = extract_coupling_energy(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out', pattern)

                            bilinear_offdiag_ev = ( Coup_ev_pp + Coup_ev_mm - Coup_ev_pm - Coup_ev_mp ) / (4.0 * qsize * qsize )

                            print(f"State {jst} & {ist} Bilinear Off-Diagonal: {bilinear_offdiag_ev}\n")
                            Bilinear.append(make_line(label=f"C1_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}", value=bilinear_offdiag_ev))

                            if SOC_flag:

                                try:

                                    # Extract DSOME_cm_pp
                                    DSOME_cm_pp = extract_DSOME(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out', nstate)
                                    DSOME_cm_pp_real, DSOME_cm_pp_imag = DSOME_cm_pp[0], DSOME_cm_pp[1]

                                    # Extract DSOME_cm_pm
                                    DSOME_cm_pm = extract_DSOME(f'{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out', nstate)
                                    DSOME_cm_pm_real, DSOME_cm_pm_imag = DSOME_cm_pm[0], DSOME_cm_pm[1]

                                    # Extract DSOME_cm_mp
                                    DSOME_cm_mp = extract_DSOME(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out', nstate)
                                    DSOME_cm_mp_real, DSOME_cm_mp_imag = DSOME_cm_mp[0], DSOME_cm_mp[1]

                                    # Extract DSOME_cm_mm
                                    DSOME_cm_mm = extract_DSOME(f'{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out', nstate)
                                    DSOME_cm_mm_real, DSOME_cm_mm_imag = DSOME_cm_mm[0], DSOME_cm_mm[1]


                                except Exception as e:
                                    print(f"Error in SOC: {str(e)}")
                else:
                    print(f"not good to extract. Skipping mode {imode} mode {jmode} for extracting bilinear vibronic couplings")

        if SOC_flag:
            for kmode in range(1, nmodes + 1):
                imode = modes_included[kmode]

                displacement_filenames = {
                    "+1": f'{filnam}_mode{imode}_+{qsize}.out',
                    "+2": f'{filnam}_mode{imode}_+{qsize}x2.out',
                    "-1": f'{filnam}_mode{imode}_-{qsize}.out',
                    "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                    # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                    # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                    # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                    # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
                }

                DSOME_cm = {}
                for key in displacement_keys:
                    try:
                        dsome_real, dsome_imag = extract_DSOME(displacement_filenames[key], nstate)
                        DSOME_cm[k] = {'real': dsome_real, 'imag': dsome_imag}
                        # DSOME_cm[key] = [dsome_real, dsome_imag]
                    except Exception as e:
                        print(f"Error in SOC: {str(e)}")
                        pass  # keep executing, we just needed to log that there was an error
                        # breakpoint()  # if you needed to investigate the cause of the error

                DSOME_cm_plus_real = DSOME_cm["+1"]['real']
                DSOME_cm_plus_imag = DSOME_cm["+1"]['imag']
                DSOME_cm_minus_real = DSOME_cm["-1"]['real']
                DSOME_cm_minus_imag = DSOME_cm["-1"]['imag']
                DSOME_cm_plusx2_real = DSOME_cm["+2"]['real']
                DSOME_cm_plusx2_imag = DSOME_cm["+2"]['imag']
                DSOME_cm_minusx2_real = DSOME_cm["+2"]['real']
                DSOME_cm_minusx2_imag = DSOME_cm["+2"]['imag']


                lmode_last = kmode - 1
                for lmode in range(1, lmode_last + 1):
                    jmode = modes_included[lmode]
                    for ist in range(1, nstate + 1):
                        jlast = ist - 1
                        for jst in range(1, jlast + 1):

                            # intialize dictionaries to contain
                            linear_SOC_cm_real = {}
                            linear_SOC_cm_imag = {}
                            quadratic_SOC_cm_real = {}
                            quadratic_SOC_cm_imag = {}
                            bilinear_SOC_cm_real = {}
                            bilinear_SOC_cm_imag = {}
                            full_Ham_SOC_cm_real = {}
                            full_Ham_SOC_cm_imag = {}

                            # Set jst ist tuple as index
                            idx = (jst, ist)

                            # Compute linear SOC
                            DSOME_cm_plus_real[idx] *= coord_disp_plus[imode]
                            DSOME_cm_plus_imag[idx] *= coord_disp_plus[imode]
                            DSOME_cm_minus_real[idx] *= coord_disp_minus[imode]
                            DSOME_cm_minus_imag[idx] *= coord_disp_minus[imode]
                            linear_SOC_cm_real[idx] = (DSOME_cm_plus_real[idx] - DSOME_cm_minus_real[idx]) / (2 * qsize)
                            linear_SOC_cm_imag[idx] = (DSOME_cm_plus_imag[idx] - DSOME_cm_minus_imag[idx]) / (2 * qsize)

                            # Compute quadratic SOC
                            DSOME_cm_plusx2_real[idx] *= coord_disp_plusx2[imode] * coord_disp_plusx2[imode]
                            DSOME_cm_plusx2_imag[idx] *= coord_disp_plusx2[imode] * coord_disp_plusx2[imode]
                            DSOME_cm_minusx2_real[idx] *= coord_disp_minusx2[imode] * coord_disp_minusx2[imode]
                            DSOME_cm_minusx2_imag[idx] *= coord_disp_minusx2[imode] * coord_disp_minusx2[imode]
                            quadratic_SOC_cm_real[idx] = (DSOME_cm_plusx2_real[idx] + DSOME_cm_minusx2_real[idx] - 2.0 * DSOME_cm_0_real[idx]) / (4.0 * qsize * qsize)
                            quadratic_SOC_cm_imag[idx] = (DSOME_cm_plusx2_imag[idx] + DSOME_cm_minusx2_imag[idx] - 2.0 * DSOME_cm_0_imag[idx]) / (4.0 * qsize * qsize)

                            # Compute bilinear SOC
                            DSOME_cm_pp_real[idx] *= coord_disp_pp[imode] * coord_disp_pp[jmode]
                            DSOME_cm_pp_imag[idx] *= coord_disp_pp[imode] * coord_disp_pp[jmode]
                            DSOME_cm_mm_real[idx] *= coord_disp_mm[imode] * coord_disp_mm[jmode]
                            DSOME_cm_mm_imag[idx] *= coord_disp_mm[imode] * coord_disp_mm[jmode]
                            DSOME_cm_pm_real[idx] *= coord_disp_pm[imode] * coord_disp_pm[jmode]
                            DSOME_cm_pm_imag[idx] *= coord_disp_pm[imode] * coord_disp_pm[jmode]
                            DSOME_cm_mp_real[idx] *= coord_disp_mp[imode] * coord_disp_mp[jmode]
                            DSOME_cm_mp_imag[idx] *= coord_disp_mp[imode] * coord_disp_mp[jmode]
                            bilinear_SOC_cm_real[idx] = (DSOME_cm_pp_real[idx] + DSOME_cm_mm_real[idx] - DSOME_cm_pm_real[idx] - DSOME_cm_mp_real[idx] ) / (4.0 * qsize * qsize )
                            bilinear_SOC_cm_imag[idx] = (DSOME_cm_pp_imag[idx] + DSOME_cm_mm_imag[idx] - DSOME_cm_pm_imag[idx] - DSOME_cm_mp_imag[idx] ) / (4.0 * qsize * qsize )

                            # Compute full SOC
                            full_Ham_SOC_cm_real[idx] = (DSOME_cm_0_real[idx] + linear_SOC_cm_real[idx] + quadratic_SOC_cm_real[idx] + bilinear_SOC_cm_real[idx])
                            full_Ham_SOC_cm_imag[idx] = (DSOME_cm_0_imag[idx] + linear_SOC_cm_imag[idx] + quadratic_SOC_cm_imag[idx] + bilinear_SOC_cm_imag[idx])

                            # Hij^(0) + lij^(1)*x_1 + lij^(2)*x_2 + 0.5qij^(1)*x_1 ^ 2 + 0.5qij^(2)*x_2 ^ 2 + bij^(1,2) * x_1 x_2
                            # Does this mean I have to extract the atom coordinates from every file too?
                            # print(imode, icomp, refcoord[icomp], nrmmod[icomp, imode], coord_disp_plus, coord_disp_minus, distcoord_plus[icomp], distcoord_minus[icomp])
                            # Probably is distcoord_plus and distcoord_minus for x1,x2 respectively

                            # Print and store results
                            SOC.append(make_line(label=f"C1_s{jst:>02d}_s{ist:>02d}_v{imode:>02d}r", value=linear_SOC_cm_real[idx], units=', cm-1'))
                            SOC.append(make_line(label=f"C1_s{jst:>02d}_s{ist:>02d}_v{imode:>02d}i", value=linear_SOC_cm_imag[idx], units=', cm-1'))
                            SOC.append("\n")

                            SOC.append(make_line(label=f"C2_s{jst:>02d}s{ist:>02d}_v{imode:>02d}r", value=quadratic_SOC_cm_real[idx], units=', cm-1'))
                            SOC.append(make_line(label=f"C2_s{jst:>02d}s{ist:>02d}_v{imode:>02d}i", value=quadratic_SOC_cm_imag[idx], units=', cm-1'))
                            SOC.append("\n")

                            SOC.append(make_line(label=f"C1_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}r", value=bilinear_SOC_cm_real[idx], units=', cm-1'))
                            SOC.append(make_line(label=f"C1_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}i", value=bilinear_SOC_cm_imag[idx], units=', cm-1'))
                            SOC.append("\n")

                            print(f"State {jst:>02d} & {ist:>02d} SOC (real) at modes {imode:>02d} & {jmode:>02d} {full_Ham_SOC_cm_real[idx]}, cm-1\n")
                            SOC.append(make_line(label=f"SOC_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}r", value=full_Ham_SOC_cm_real[idx], units=', cm-1'))

                            print(f"State {jst:>02d} & {ist:>02d} SOC (imag) at modes {imode:>02d} & {jmode:>02d} {full_Ham_SOC_cm_imag[idx]}, cm-1\n")
                            SOC.append(make_line(label=f"SOC_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}i", value=full_Ham_SOC_cm_imag[idx], units=', cm-1'))
                            SOC.append("\n")

        file_contents = params + ['\n',]
        file_contents += linear + ['\n',]
        file_contents += quadratic + ['\n',]
        file_contents += bilinear + ['\n',]
        file_contents += SOC + ['\n',]

        # Params.append("\n")
        # Params.extend(Linear)
        # Params.append("\n")
        # Params.extend(Quadratic)
        # Params.append("\n")
        # Params.extend(Bilinear)
        # Params.append("\n")
        # Params.extend(SOC)
        # Params.append("\n")

        # different header style
        hfs = header_format_string = spacer_format_string + "# {:^60s} #\n" + spacer_format_string

        file_contents += _make_ETD_block(header_format_string)
        file_contents += "end-parameter-section\n"

        header_name_list = [
            "HAMILTONIAN-SECTION",
            "KINETIC OPERATOR FOR NORMAL MODES",
            "HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES",
            "ELECTRONIC COUPLING AT REFERENCE STRUCTURE",
            "LINEAR DIAGONAL VIBRONIC COUPLINGS",
            "LINEAR OFF-DIAGONAL VIBRONIC COUPLINGS",
            "QUADRATIC DIAGONAL VIBRONIC COUPLINGS",
            "QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS",
            "BILINEAR DIAGONAL VIBRONIC COUPLINGS",
            "BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS",
        ]

        # Open mctdh.op file for writing
        with open('mctdh.op', 'a') as mctdh_file:

            labels = [

            ]

            # this part isn't done yet ---- IN PROGRESS

            def _make_hamiltonian_section():

                block = header_format_string.format(header_name_list[0])

                # Write modes and mode labels
                mode_labels = [f"v{n:>02d}" for n in pp.selected_mode_list]
                block += "modes | el | " + " | ".join(mode_labels) + "\n"

                for i, label in enumerate(pp.selected_mode_list):
                    block += f"1.00*w{label:>02d}   |{i+2} KE\n"
                    block += f"0.50*w{label:>02d}   |{i+2} q^2\n"

                for a in range(1, A+2):
                    block += f"EH_s{a:>02d}_s{a:>02d} |1 S{a}&{a}\n"

                block += "\n"

                _list1 = [
                    "EH_s{}_s{}",
                    "C1_s{}_s{}_v{}",
                    "C1_s{}_s{}_v{}_v{}",
                    "C2_s{}_s{}_v{}_v{}",
                ]

                _list2 = [
                    "|1 S{a:}&{a:} |{i:} q",
                    "|1 S{b:}&{a:} |{i:} q",
                    "|1 S{a:}&{a:} |{i:} q^2",
                    "|1 S{a:}&{a:} |{i:} q |{j:} q",
                    "|1 S{b:}&{a:} |{i:} q |{j:} q",
                ]

            def label_linear_coupling(linear_terms, A, N):
                """Return a string containing the linear coupling constant labelling of a .op file."""
                spacer = '|'
                if diagonal:
                    return '\n'.join([
                        f"C1_s{a:0>2d}_s{a:0>2d}_v{i:0>2d}{spacer:>11}1 S{a:d}&{a:d}{spacer:>4}{i+1}  q"
                        for a, i in it.product(range(1, A+1), range(1, N+1))
                        if not np.isclose(linear_terms[i-1, a-1], 0.0)
                    ]) + '\n'
                else:
                    return '\n'.join([
                        f"C1_s{a:0>2d}_s{a:0>2d}_v{i:0>2d}{spacer:>11}1 S{a:d}&{a:d}{spacer:>4}{i+1}  q"
                        for a, i in it.product(range(1, A+1), range(1, N+1))
                    ] + [
                        ''  # creates a blank line between the (surface) diagonal and off-diagaonl linear terms
                    ] + [
                        f"C1_s{a2:0>2d}_s{a1:0>2d}_v{i:0>2d}{spacer:>11}1 S{a2:d}&{a1:d}{spacer:>4}{i+1}  q"
                        for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
                        if (a1 != a2)
                    ]) + '\n'


                for a in range(1, A+1):
                    for b in range(1, a):
                        block += f"EH_s{b:>02d}_s{a:>02d}  |1 S{b}&{a}\n"

                # Write LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS
                for i, a in it.product(range(1, N+1), range(1, A+1)):
                    i_label = mode_map_dict[i]
                    block += f"C1_s{a:>02d}_s{a:>02d}_v{i_label:>02d} |1 S{a}&{a} |{i+1} q\n"

                # Write LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS
                for i, a in it.product(range(1, N+1), range(1, A+1)):
                    i_label = mode_map_dict[i]
                    for b in range(1, a):
                        block += (
                            f"C1_s{b:>02d}_s{a:>02d}_v{i_label:>02d}"
                            f" |1 S{b}&{a} |{i+1} q\n"
                        )

                # Write LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS
                for i, a in it.product(range(1, N+1), range(1, A+1)):
                    i_label = mode_map_dict[i]
                    block += (
                        f"0.50*C2_s{a:>02d}s{a:>02d}_v{i_label:>02d}v{i_label:>02d}"
                        f" |1 S{a}&{a} |{i+1} q^2\n"
                    )

                # Write BILINEAR DIAGONAL VIBRONIC COUPLINGS
                for i, a in it.product(range(1, N+1), range(1, A+1)):
                    for j in range(1, i):
                        i_label, j_label = mode_map_dict[[i, j]]
                        block += (
                                f"C1_s{a:>02d}s{a:>02d}_v{i_label:>02d}v{j_label:>02d}"
                                f" |1 S{a}&{a} |{i+1} q |{j+1} q\n"
                            )

                # Write BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS
                for i, a in it.product(range(1, N+1), range(1, A+1)):
                    for j, b in it.product(range(1, i), range(1, a)):
                        i_label, j_label = mode_map_dict[[i, j]]
                        block += (
                                f"C1_s{b:>02d}s{a:>02d}_v{i_label:>02d}v{j_label:>02d}"
                                f" |1 S{b}&{a} |{i+1} q |{j+1} q\n"
                            )

                # ------------------------------------------------------------
                # Write KINETIC OPERATOR FOR NORMAL MODES (mostly fine)
                for imode_include in range(1, nmodes + 1):
                    mode_count = imode_include + 1
                    mctdh_file.write(f"1.00*w{modes_included[imode_include]:>02d}   |{mode_count} KE\n")

                # Write HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES
                for imode_include in range(1, nmodes + 1):
                    mode_count = imode_include + 1
                    mctdh_file.write(f"0.50*w{modes_included[imode_include]:>02d}   |{mode_count}  q^2\n")

                # Write ELECTRONIC COUPLING AT REFERENCE STRUCTURE
                for ist in range(1, nstate + 2):
                    mctdh_file.write(f"EH_s{ist:>02d}_s{ist:>02d} |1 S{ist}&{ist}\n")

                mctdh_file.write("\n")
                for ist in range(1, nstate + 1):
                    jlast = ist - 1
                    for jst in range(1, jlast + 1):
                        mctdh_file.write(f"EH_s{jst:>02d}_s{ist:>02d}  |1 S{jst}&{ist}\n")

                # Write LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS
                for kmode in range(1, nmodes + 1):
                    imode = modes_included[kmode]
                    kmode_count = kmode + 1
                    for ist in range(1, nstate + 1):
                        mctdh_file.write(f"C1_s{ist:>02d}_s{ist:>02d}_v{imode:>02d} |1 S{ist}&{ist} |{kmode_count} q\n")

                # Write LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS
                for kmode in range(1, nmodes + 1):
                    imode = modes_included[kmode]
                    kmode_count = kmode + 1
                    for ist in range(1, nstate + 1):
                        jlast = ist - 1
                        for jst in range(1, jlast + 1):
                            mctdh_file.write(f"C1_s{jst:>02d}_s{ist:>02d}_v{imode:>02d} |1 S{jst}&{ist} |{kmode_count} q\n")

                # Write LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS
                for kmode in range(1, nmodes + 1):
                    imode = modes_included[kmode]
                    kmode_count = kmode + 1
                    for ist in range(1, nstate + 1):
                        mctdh_file.write(f"0.50*C2_s{ist:>02d}s{ist:>02d}_v{imode:>02d}v{imode:>02d} |1 S{ist}&{ist} |{kmode_count} q^2\n")

                # Write LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS
                for kmode in range(1, nmodes + 1):
                    imode = modes_included[kmode]
                    kmode_count = kmode + 1
                    for ist in range(1, nstate + 1):
                        jlast = ist - 1
                        for jst in range(1, jlast + 1):
                            mctdh_file.write(f"0.50*C2_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{imode:>02d} |1 S{jst}&{ist} |{kmode_count} q^2\n")

                # Write BILINEAR DIAGONAL VIBRONIC COUPLINGS
                for kmode in range(1, nmodes + 1):
                    imode = modes_included[kmode]
                    kmode_count = kmode + 1
                    lmode_last = kmode - 1
                    for lmode in range(1, lmode_last + 1):
                        jmode = modes_included[lmode]
                        lmode_count = lmode + 1
                        for ist in range(1, nstate + 1):
                            mctdh_file.write(f"C1_s{ist:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d} |1 S{ist}&{ist} |{lmode_count} q |{kmode_count} q\n")

                # Write BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS
                for kmode in range(1, nmodes + 1):
                    imode = modes_included[kmode]
                    kmode_count = kmode + 1
                    lmode_last = kmode - 1
                    for lmode in range(1, lmode_last + 1):
                        jmode = modes_included[lmode]
                        lmode_count = lmode + 1
                        for ist in range(1, nstate + 1):
                            jlast = ist - 1
                            for jst in range(1, jlast + 1):
                                mctdh_file.write(f"C1_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d} |1 S{jst}&{ist} |{lmode_count} q |{kmode_count} q\n")


                mctdh_file.write("-----------------------------------------\n")
                mctdh_file.write("\nend-hamiltonian-section\n\n")

            if SOC_flag:  # all the code inside this IF should eventually move to a seperate function (like factored out)

                soc_extension_list_string = []

                key_order = [  # (THE ORDER IS VERY IMPORTANT)
                    'Linear-Real',
                    'Linear-Imag',
                    'Quadratic-Real',
                    'Quadratic-Imag',
                    'Bilinear-Real',
                    'Bilinear-Imag',
                    'SOC-Real',
                    'SOC-Imag',
                ]
                hamiltonian_blocks = {k: "" for k in key_order}
                # ----------------------------------------------------------------------

                # prepare `make_line`
                format_string = "{label:<25s}{link:<20s}\n"
                make_line = functools.partial(format_string.format)

                # Write FULL HAMILTONIAN SOC OFF-DIAGONAL VIBRONIC COUPLINGS

                """ there is ways to do this, but it may not be worth the effort right now
                for k, l, i, j in it.product():
                    hamiltonian_blocks['Linear-Real'] += make_line(label=f"I*C1_s{jst:>02d}_s{ist:>02d}_v{imode:>02d}r", link=f"|1 Z{jst}&{ist} | {kmode_count} q")  # noqa: E501
                    hamiltonian_blocks['Linear-Real'] += make_line(label=f"-I*C1_s{jst:>02d}_s{ist:>02d}_v{imode:>02d}i", link=f"|1 Z{ist}&{jst} | {kmode_count} q")  # noqa: E501
                for k, l, i, j in it.product():
                    hamiltonian_blocks['Quadratic-Real'] += make_line(label=f"I*C2_s{jst:>02d}s{ist:>02d}_v{imode:>02d}r", link=f"|1 Z{jst}&{ist} | {kmode_count} q^2")  # noqa: E501
                    hamiltonian_blocks['Quadratic-Real'] += make_line(label=f"-I*C2_s{jst:>02d}s{ist:>02d}_v{imode:>02d}i", link=f"|1 Z{ist}&{jst} | {kmode_count} q^2")  # noqa: E501
                <...> (and so forth)
                """

                # do the work
                for kmode in range(1, nmodes + 1):
                    imode = modes_included[kmode]
                    lmode_last = kmode - 1
                    for lmode in range(1, lmode_last + 1):
                        jmode = modes_included[lmode]
                        for ist in range(1, nstate + 1):
                            jlast = (ist - 1)
                            for jst in range(1, jlast + 1):

                                # note to self: the I* is performing ARITHMETIC on SOr_{jst}_{ist} prepared earlier, does that mean we neeed to remove the l and _m{imode}
                                hamiltonian_blocks['Linear-Real'] += make_line(label=f"I*C1_s{jst:>02d}_s{ist:>02d}_v{imode:>02d}r", link=f"|1 Z{jst}&{ist} | {kmode_count} q")  # noqa: E501
                                hamiltonian_blocks['Linear-Imag'] += make_line(label=f"-I*C1_s{jst:>02d}_s{ist:>02d}_v{imode:>02d}i", link=f"|1 Z{ist}&{jst} | {kmode_count} q")  # noqa: E501

                                hamiltonian_blocks['Quadratic-Real'] += make_line(label=f"I*C2_s{jst:>02d}s{ist:>02d}_v{imode:>02d}r", link=f"|1 Z{jst}&{ist} | {kmode_count} q^2")  # noqa: E501
                                hamiltonian_blocks['Quadratic-Imag'] += make_line(label=f"-I*C2_s{jst:>02d}s{ist:>02d}_v{imode:>02d}i", link=f"|1 Z{ist}&{jst} | {kmode_count} q^2")  # noqa: E501

                                hamiltonian_blocks['Bilinear-Real'] += make_line(label=f"I*C1_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}r", link=f"|1 Z{jst}&{ist} | {lmode_count} q |{kmode_count} q")  # noqa: E501
                                hamiltonian_blocks['Bilinear-Imag'] += make_line(label=f"-I*C1_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}i", link=f"|1 Z{ist}&{jst} | {lmode_count} q |{kmode_count} q")  # noqa: E501

                                hamiltonian_blocks['SOC-Real'] += make_line(label=f"I*SOC_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}r", link=f"|1 Z{jst}&{ist} | {kmode_count} q")  # noqa: E501
                                hamiltonian_blocks['SOC-Imag'] += make_line(label=f"-I*SOC_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{jmode:>02d}i", link=f"|1 Z{ist}&{jst} | {kmode_count} q")  # noqa: E501

                for k in key_order:  # glue the blocks together
                    output_string += hamiltonian_blocks[k] + "\n"

                print("Hey check the output string, and the hamiltonian_blocks"); breakpoint()
                mctdh_file.write(output_string)
                del hamiltonian_blocks  # don't need anymore

        if SOC_flag:
            file_contents += _make_SOC_section()

        file_contents += _make_Hamiltonian_operate_Ex_section()
        file_contents += "end-operator\n"

        with open("mctdh.op", "w") as fp:
            fp.write("".join(file_contents))

    return

# ---------------------------------------------------------------------------------------
# helper functions for `main()`


def extract_lines_between_patterns(filename, start_pattern, end_pattern, collecting=False):
    """ Function to extract lines between patterns in a file """
    selected_lines = []
    # would be good to replace this with memory mapping find or grep command?
    with open(filename, 'r', errors='replace') as file:
        for line in file:
            if start_pattern in line:
                collecting = True
                selected_lines.append(line)
            elif end_pattern in line:
                collecting = False
            elif collecting:
                selected_lines.append(line)

    return selected_lines


def read_freq_values(hessout):
    """ Function to read frequency values from selected lines """
    selected_lines = extract_lines_between_patterns(
        hessout,
        "FREQUENCIES IN CM",
        "REFERENCE ON SAYVETZ CONDITIONS"
    )

    freq_value_set = []
    for freqline in selected_lines:
        if "FREQUENCY:" in freqline:
            freq_value_set.append(freqline[18:])

    return freq_value_set


def read_mode_values(hessout):
    """ Function to extract filtered set of lines """
    selected_lines = extract_lines_between_patterns(
        hessout,
        "FREQUENCIES IN CM",
        "REFERENCE ON SAYVETZ CONDITIONS"
    )
    mode_value_set = []

    for idx, modeline in enumerate(selected_lines):
        if len(modeline) > 3 and modeline[2].isdigit():
            # mode_value_set.append(selected_lines[idx][20:])
            # mode_value_set.append(selected_lines[idx+1][20:])
            # mode_value_set.append(selected_lines[idx+2][20:])
            for i in range(3):  # 3 for x,y,z
                mode_value_set.append(selected_lines[idx + i][20:])

    return mode_value_set


def get_number_of_atoms(hessout):
    """ Function to get the number of atoms from the hessout file """

    # would be good to replace this with memory mapping find or grep command?
    with open(hessout, 'r', errors='replace') as hess_file:
        for line in hess_file:
            if ' TOTAL NUMBER OF ATOMS' in line:
                natoms = int(line.split('=')[1])
                return natoms


def _extract_freq_and_mode_from_hessian(path):
    """ x """
    freq_value_set, filtered_set = read_freq_values(path), read_mode_values(path)

    with open('mode.dat', 'w') as fp:
        fp.writelines(filtered_set)

    with open('freq.dat', 'w') as fp:
        fp.writelines(freq_value_set)

    return


def process_mode_freq(ndim, nof_cols=5, float_length=12):
    """ """
    nof_groups = ndim // nof_cols  # integer division
    nof_leftover_modes = ndim % nof_cols
    print(
        f"Dimension of all xyz coordinates: {ndim}\n"
        f"{ndim / 3} atoms, split into {nof_groups} groups with {nof_leftover_modes} left over\n",
    )

    # -------------------------------------------------------------------------
    with open("mode.dat", 'r', errors='replace') as mode_file:
        lines_mode = mode_file.readlines()

    mode_list = [[float(n) for n in line.strip().split()] for line in lines_mode]

    # glue the groups onto each other (we only care about the first ndim)
    for i in range(ndim):
        for g in range(nof_groups):
            mode_list[i].extend(mode_list[i+(ndim*(g+1))])

    # throw away the lists we don't need anymore
    mode_list = mode_list[:ndim]

    # turn into numpy array
    modes_array = np.array(mode_list)  # should be square?

    # -------------------------------------------------------------------------
    with open("freq.dat", 'r', errors='replace') as freq_file:
        lines_freq = freq_file.readlines()

    freq_list = [[float(n) for n in line.strip().replace('i', '').split()] for line in lines_freq]
    frequences = list(it.chain(*freq_list))  # single list
    freq_array = np.array(frequences)

    # freqcm = {}
    # for i in range(len(frequences)):
    #     key = pp.modes_included[i]
    #     freqcm[key] = frequences[i]

    # -------------------------------------------------------------------------
    if False and __debug__:   # print all frequencies
        string = "\n".join([f"frequency: {i} {freq_array[i]} CM-1" for i in range(ndim)])
        print(string)

    # return mode_array, freqcm
    return modes_array, freq_array


def compose_ref_structure(ref_file, hessout, nof_atoms):
    coord_lines = extract_lines_between_patterns(hessout, 'EQUILIBRIUM GEOMETRY LOCATED', 'INTERNUCLEAR DISTANCES')

    good_ref_structure = bool(len(coord_lines) > 2)

    if not good_ref_structure:
        print(f'Unsuccessful extraction of equilibrium geometry from {hessout}. Please prepare ref_structure manually.')
        breakpoint()  # do we want to continue execution, or do we actually want to stop the program?
        import sys; sys.exit()

    if good_ref_structure:

        """ we don't need to delete the file if we simply write a new file
        try:
            subprocess_run_wrapper(['rm', '-f', ref_file])
        except Exception as e:
            print(f"Error deleting {ref_file}: {str(e)}")
        """

        # the last element of coord_lines is an empty line (`\n`)
        assert coord_lines[-1] == '\n'

        # remove empty lines
        coord_lines = [l for l in coord_lines if l != '\n']

        # we want the lines at the end of the file (the last nof_atoms/Z lines)
        file_contents = "".join(coord_lines[-nof_atoms:])

        with open(ref_file, 'w') as fp:
            fp.write(file_contents)

        print(f'Successfully extracted equilibrium geometry from {hessout} and prepared ref_structure.')

    return coord_lines


def read_reference_structure(file_path, verbose=True):
    """ Read in the reference structure.
    The ref_structure has to be prepared by human-being and adopts the following format
        ##### SAMPLE REF STRUCT #######
        N           7.0   0.0000000000  -0.0000000000  -0.1693806842
        H           1.0  -0.4653267700   0.8059696078   0.2564602281
        H           1.0  -0.4653267700  -0.8059696078   0.2564602281
        H           1.0   0.9306535400   0.0000000000   0.2564602281
        ###############################

    This function
    """

    # atom_dict, charge_dict, ref_coords = {}, {}, {}
    atom_dict, charge_dict = {}, {}
    temp_coord_list = []

    with open(file_path, 'r', errors='replace') as struct_file:
        lines = struct_file.readlines()

    for i, line in enumerate(lines):
        atom_name, charge, *coords = line.split()

        assert isinstance(atom_name, str), f"{atom_name=} is not a string?"
        assert isinstance(charge, str), f"{charge=} is not a string?"
        assert isinstance(coords, list) and len(coords) == 3, f"{charge=} is not a list of length 3?"

        if verbose and __debug__: print(atom_name, charge, coords)

        atom_dict[i+1], charge_dict[i+1] = atom_name, charge
        temp_coord_list.append(coords)

    # flatten list and apply float to each element
    values = [*map(float, it.chain(*temp_coord_list))]
    keys = range(1, len(values)+1)  # numbers 1,2,3, ...

    # make the reference co-ordinate dictionary
    ref_coords = dict(zip(keys, values))

    if verbose and __debug__: print(atom_dict, charge_dict, ref_coords, sep='\n')

    return atom_dict, charge_dict, ref_coords


def refG_calc(refgeo, input_filename, output_filename):
    """
    Do diabatization calculation at the reference non-distorted structure.
    This calculation shall be a repetition of a calculation in preparing `temp.inp`.
    """

    # Check if the calculation has already been run
    grace_exists = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", output_filename]) == 0

    if grace_exists:
        print("Calculation at the reference structure has already been done.")
        return

    else:
        print("Run calculation at the undistorted reference structure")

        shutil.copy("temp.inp", input_filename)

        with open(refgeo, 'r', errors='replace') as ref_structure:
            data = ref_structure.read()

        # in this case we append the reference structure to file contents from "temp.inp"
        with open(input_filename, "a") as fp:
            fp.write(data)
            fp.write(" $END\n")

        # Submit and run the refG calculation (you may need to customize this command based on your setup)
        # refG_job_result = subprocess_run_wrapper(["./subgam.diab", input_filename, "4", "0", "1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os_system_wrapper("sbatch" + ' -W' + " " + my_subgam(input_filename, ncpus=2, ngb=1, nhour=1))  # the wait

        # At this point, refG calculation has completed successfully.
        print("Calculation at the reference structure is done.")

    return

# ---------------------------------------------------------------------------------------


def main(ref_file="ref_structure", ncols=5, **kwargs):
    """ x """
    hessian_filename = kwargs['hessian_filename']
    _extract_freq_and_mode_from_hessian(hessian_filename)

    nof_atoms = natoms = get_number_of_atoms(hessian_filename)
    ndim = nof_atoms * 3

    assert Z == nof_atoms, f"{Z=} is not {nof_atoms=}!?"
    assert N_tot == ndim, f"{N_tot=} is not {ndim=}!?"

    nrmmod, freqcm = process_mode_freq(N_tot, ncols)

    compose_ref_structure(ref_file, hessian_filename, nof_atoms)

    atmlst, chrglst, refcoord = read_reference_structure(ref_file)

    refG_calc(ref_file, kwargs['refG_in'], kwargs['refG_out'])

    # -------------------------------------------------------------------------
    diabatization_kwargs = kwargs.copy()
    diabatization_kwargs.update({
        'ndim': ndim,
        'freqcm': freqcm,
        'refcoord': refcoord,
        'nrmmod': nrmmod,
        'natoms': natoms,
        'atmlst': atmlst,
        'chrglst': chrglst,
        'qsize': pp.qsize,
        'ha2ev': pp.ha2ev,
        'wn2ev': pp.wn2ev,
        'wn2eh': pp.wn2eh,
        'ang2br': pp.ang2br,
        'amu2me': pp.amu2me
    })
    # name, modes = pp.filnam, pp.modes_included
    diabatize = diabatization(**diabatization_kwargs)
    print("Diabatization successfully modified")

    tdipole_block = extract_lines_between_patterns(
        kwargs['refG_out'],
        "TRANSITION DIPOLES BETWEEN DIABATS",
        "TRANSITION DIPOLES BETWEEN DIRECT MAX. DIABATS"
    )

    dipoles = extract_ground_to_excited_state_transition_dipoles(tdipole_block)

    pprint.pprint(tdipole_block)
    pprint.pprint(dipoles)

    mctdh_input_kwargs = kwargs.copy()
    mctdh_input_kwargs.update({
        'qsize': pp.qsize,
        'ha2ev': pp.ha2ev,
        'wn2ev': pp.wn2ev,
        'wn2eh': pp.wn2eh,
        'ang2br': pp.ang2br,
        'amu2me': pp.amu2me,
        'nof_electronic_states': pp.A,
        'ndim': ndim,
        'freqcm': freqcm,
        'nrmmod': nrmmod,
        'dipoles': dipoles,
        'diabatize': diabatize,
        'hessout': kwargs['hessian_filename'],
    })

    # name, modes = pp.filnam, pp.modes_included
    mctdh(**mctdh_input_kwargs)
    print("mctdh successfully modified"); return

    print('The run was a success!')

    shutil.copy("mctdh.op", kwargs['project_filename'] + '.op')

    if (extra_debug := False):
        print('---------nrm mod done-----------')
        pprint.pprint(nrmmod)
        print('---------freqcm done-----------')
        pprint.pprint(freqcm)
        print('---------selected_mode_list done-----------')
        pprint.pprint(selected_mode_list)
        print('---------atmlst done-----------')
        pprint.pprint(atmlst)
        print('---------chrglst done-----------')
        pprint.pprint(chrglst)
        print('---------refcoord done-----------')
        pprint.pprint(refcoord)
    # ...
    return


# ---------------------------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python your_script.py <hessout_file>")
        sys.exit(1)

    hessian_filename = sys.argv[1]  # read in name of hessian file
    kwargs = {'hessian_filename': hessian_filename}

    # everything we add to `kwargs` is 'imported' from project parameters

    kwargs.update({
        'project_filename': pp.filnam,
        'refG_in': f"{pp.filnam}_refG.inp",  # reference geometry
        'refG_out': f"{pp.filnam}_refG.out",  # reference geometry
        # 'modes_included': pp.modes_included,
    })

    main(**kwargs)
