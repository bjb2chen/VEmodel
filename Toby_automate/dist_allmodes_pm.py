""" x """

# system imports
import os
from os.path import splitext
import sys
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
import project_parameters
from project_parameters import *  # eventually remove this


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
            subprocess.run(['rm', '-f', filename])
        except Exception as e:
            print(f"Error deleting {filename}: {str(e)}")
    return

# ---------------------------------------------------------------------------------------

# this function should be moved to project_parameters.py!
# def filter_modes(excluded_set, ndim):
#     modes_included = {}
#     print("NUMBER OF EXCLUDED MODES:", len(excluded_set))
#     print("They are modes:", *excluded_set)

#     # Iterate through modes to check inclusion
#     for imode in range(1, ndim + 1):
#         if imode not in excluded_set:
#             modes_included[len(modes_included)+1] = imode
#             print(len(modes_included), imode)

#     print("Number of Modes Included:", len(modes_included))

#     return modes_included


def my_subgam(path, **kwargs):
    """ x """

    # if its easier to just change project parameters i would recommend doing
    # ncpus = project_parameters.ncpus
    # nhour = project_parameters.nhour
    # ngb = project_parameters.ngb

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

    command = (
        "sbatch"
        f" {input_no_ext}.slurm"
    )

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
    modes_included = kwargs['modes_included']
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
    # amu2me = project_parameters.amu2me
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
    freq_array = freqcm  # assume we return the array
    mode_array = nrmmod  # assume we return the array
    atom_list = [*atmlst.values()]
    charge_list = [*chrglst.values()]

    # make a list for easy calculations
    ref_coord_array = np.array([*refcoord.values()])

    # N should always be an integer!!
    # instead just use `modes_included` when you need to loop over stuff
    N = len(project_parameters.N)
    nof_modes = N  # alias

    NEW = np.newaxis  # alias

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
        `R_array` is assumed to be of dimension (ndim)
        `mode_array` is assumed to be of dimension (ndim, nof_modes)
        Use broadcasting.

        The first dimension iterates over the atoms and their xyz coordinates.
        The second dimension iterates over the normal modes.

        ndim == nof_modes (as far as I can tell?)
        """

        # the `:` of the second dimension of mode_array is what we index with i or j
        displacements = {
            "+1": reference[:, NEW] + 1.0 * R_array[NEW, :] * mode_array[:, :],
            "-1": reference[:, NEW] - 1.0 * R_array[NEW, :] * mode_array[:, :],
            "+2": reference[:, NEW] + 2.0 * R_array[NEW, :] * mode_array[:, :],
            "-2": reference[:, NEW] - 2.0 * R_array[NEW, :] * mode_array[:, :],
        }
        assert set(displacements.keys()) == set(linear_disp_keys), f"{linear_disp_keys=} no longer agree!"

        for k in linear_disp_keys:
            assert displacements[k].shape == (ndim, ndim), f"{k=} {displacements[k].shape=} not {(ndim, nof_modes)=}?"

        # store the displacements in the `distored_coords` dictionary
        for key in linear_disp_keys:
            distored_coords[key] = displacements[key]

        return  # don't need to return value, as `distored_coords` is a dictionary

    def _precompute_bilinear_displacements(R_array, mode_array, distored_coords):
        """
        `R_array` is assumed to be of dimension (ndim)
        `mode_array` is assumed to be of dimension (ndim, nof_modes, nof_modes)
        Use broadcasting.

        The first dimension iterates over the atoms and their xyz coordinates.
        The second dimension iterates over the normal modes (i).
        The third dimension iterates over the normal modes (j).

        ndim == nof_modes (as far as I can tell?)
        """

        displacements = {
            "++": distored_coords["+1"][:, :, NEW] + R_array[NEW, NEW, :] * mode_array[:, NEW, :],
            "+-": distored_coords["+1"][:, :, NEW] - R_array[NEW, NEW, :] * mode_array[:, NEW, :],
            "-+": distored_coords["-1"][:, :, NEW] + R_array[NEW, NEW, :] * mode_array[:, NEW, :],
            "--": distored_coords["-1"][:, :, NEW] - R_array[NEW, NEW, :] * mode_array[:, NEW, :],
        } 
        assert set(displacements.keys()) == set(bi_linear_disp_keys), f"{bi_linear_disp_keys=} no longer agree!"

        for k in bi_linear_disp_keys:
            shape = (ndim, nof_modes*3, nof_modes*3)  # 3 represents 3 xyz coordinates for each mode?
            assert displacements[k].shape == shape, f"{k=} {displacements[k].shape=} not {shape=}?"

        # store the displacements in the `distored_coords` dictionary
        for key in bi_linear_disp_keys:
            distored_coords[key] = displacements[key]

        return  # don't need to return value, as `distored_coords` is a dictionary

    def _save_distorted_structure(mode_idx, displaced_q, charge_list, atom_list, filename_list, key_list):
        """ save the distorted structure to the `distored_structure_filenames`
        `displaced_q` is an array of dimension (ndim, ndim)
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
            with open(path, 'a') as fp:
                fp.write("\n".join(file_contents))

        return

    def _create_linear_diabatization_input_files(filename, q1_label, qsize):
        """ this could be improved somewhat

        The `q1_label` goes from (1 -> nof_modes+1) rather than (0 -> nof_modes).
        """

        refG_out = f"{filename}_refG.out"

        # Check if the reference geometry calculation is done?
        grace0 = subprocess.run(["grep", "DONE WITH MP2 ENERGY", refG_out])
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
                with open(distored_struct_file, 'r', errors='replace') as fp:
                    data = fp.read()

                with open(games_filename+'.inp', 'a') as fp:
                    fp.write(data)  # can you just do data + ' $END' in one write?
                    fp.write('\n')
                    fp.write(' $END')

                grace1 = subprocess.run(["grep", "DONE WITH MP2 ENERGY", games_filename+'.out'])
                gamess_calculation_not_run = bool(grace1.returncode != 0)

                # This means that refG completed successfully and `diabmode*.out` not completed
                if ref_geom_flag_exists and gamess_calculation_not_run:
                    print(f"Running calculations for {games_filename}")
                    try:
                        os.system("sbatch" + " " + my_subgam(games_filename+'.inp', ncpus=2, ngb=1, nhour=1))
                    except Exception as e:
                        print(f"Error running diabatization calculation: {str(e)}")
                else:
                    print(f"{games_filename} is done")

        return

    def _create_bilinear_diabatization_input_files(filename, q1_label, q2_label, qsize):
        """ this could be improved somewhat
        """

        refG_out = f"{filename}_refG.out"

        # Check if the reference geometry calculation is done?
        grace0 = subprocess.run(["grep", "DONE WITH MP2 ENERGY", refG_out])
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
                fp.write('\n')
                fp.write(' $END')

            # Check if the calculation is done already
            grace2 = subprocess.run(["grep", "DONE WITH MP2 ENERGY", games_filename + '.out'])
            gamess_calculation_not_run = bool(grace2.returncode != 0)

            # this will never work? grace0 is not defined
            if ref_geom_flag_exists and gamess_calculation_not_run:
                print(f"Running calculations for {games_filename}!")
                try:
                    os.system("sbatch" + " " + my_subgam(games_filename+'.inp', ncpus=2, ngb=1, nhour=1))
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

    if precompute_bilinear:  # modify disp_coord in-place (may take a lot of memory?)
        _precompute_bilinear_displacements(rsize, mode_array, disp_coord)

    # Loop over modes and do linear displacements
    for i in range(nof_modes):
        q1_label = modes_included[i+1]

        _remove_existing_distorted_structure_files(linear_disp_filenames)
        index = (i, )
        _save_distorted_structure(
            index, disp_coord, charge_list, atom_list,
            linear_disp_filenames,
            linear_disp_keys
        )

        _create_linear_diabatization_input_files(filnam, q1_label, qsize)

        # 2D distortion to get bilinear vibronic coupling
        for j in range(i):
            q2_label = modes_included[j+1]

            if not precompute_bilinear:  # modify disp_coord in-place
                _compute_bi_linear_displacements(i, j, rsize, mode_array, disp_coord)

            _remove_existing_distorted_structure_files(bi_linear_disp_filenames)
            index = (i, j) if precompute_bilinear else (j, )
            _save_distorted_structure(
                index, disp_coord, charge_list, atom_list,
                bi_linear_disp_filenames,
                bi_linear_disp_keys,
            )

            _create_bilinear_diabatization_input_files(filnam, q1_label, q2_label, qsize)

    return disp_coord


# ---------------------------------------------------------------------------------------
# Now we move on to extract vibronic coupling constants using finite difference
# and write the data in an mctdh operator file


def extract_ground_state_energy(hessout, pattern):
    try:
        command = f'grep "{pattern}" {hessout} | tail -1 | cut -c40-'
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        output = float(result.stdout.strip().replace(" ", ""))
        return output
    except Exception as e:
        print("Cannot find ground state energy")
        return None


def refG_extract(file_path, pattern):
    try:
        # Use subprocess.run with the direct command
        command = f'grep "{pattern}" "{file_path}" | tail -1 | cut -c62-'
        result = subprocess.run(command, shell=True, text=True, capture_output=True)

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
        result = subprocess.run(command, shell=True, text=True, capture_output=True)

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
        command = f'grep "{pattern}" "{file_path}" | tail -1 | cut -c62-'
        result = subprocess.run(command, shell=True, text=True, capture_output=True)

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


def extract_DSOME(filnam, nstate):

    start_pattern = "HSO MATRIX IN DIABATIC REPRESENTATION (DIRECT MAXIMIZATION)"
    end_pattern = 'SOC EIG. VALUES and VECTORS IN DIABATS (DIRECT MAX.)'

    sed_command = f"sed -n '/{start_pattern}/,/{end_pattern}/p' {filnam}"
    result = subprocess.run(sed_command, shell=True, text=True, capture_output=True)

    if result.stdout.splitlines(): # True if not empty list
        selected_lines = result.stdout.splitlines()
    else:
        print('Cannot use sed to extract')
        selected_lines = extract_lines_between_patterns(
            f'{filnam}',
            f'{start_pattern}',
            f'{end_pattern}',
        )
        print(f'Using selected lines from {filnam}, opened via python')
    # except Exception as e:
    #     print('Cannot use sed to extract')
    #     print(f'Using selected lines from {filnam}, opened via python')

    # pprint.pprint(selected_lines)

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

    for left_state_idx in range(1, int(nstate)):
        for right_state_idx in range(left_state_idx+1, int(nstate)+1):
            for level_idx in range(1, 3):
                full_extracted_set[left_state_idx, right_state_idx, level_idx] = DSOME_set[f'{left_state_idx} & {right_state_idx}, {level_idx}']

            summed_set_real[left_state_idx, right_state_idx] = full_extracted_set[left_state_idx, right_state_idx, 1].real + \
                                                               full_extracted_set[left_state_idx, right_state_idx, 2].real

            summed_set_imag[left_state_idx, right_state_idx] = full_extracted_set[left_state_idx, right_state_idx, 1].imag + \
                                                               full_extracted_set[left_state_idx, right_state_idx, 2].imag

            append_J[left_state_idx, right_state_idx] = complex(0,summed_set_imag[left_state_idx, right_state_idx])

    #return full_extracted_set, summed_set_real, summed_set_imag, append_J
    return [summed_set_real, summed_set_imag]


def mctdh(**kwargs):
    """ description of function """

    try:  # remove the previous file, as we are about to write to it
        subprocess.run(['rm', '-f', 'mctdh.op'])
    except Exception as e:
        print(f"Error deleting {'mctdh.op'}: {str(e)}")

    # -------------------------------------------------------------------------

    filnam = filename = kwargs['project_filename']
    modes_included = kwargs['modes_included']

    # extract necessary parameters
    nstate = kwargs['nof_electronic_states']
    nmodes = len(modes_included)
    ndim = kwargs.get('ndim')
    freqcm = kwargs.get('freqcm')
    qsize = kwargs.get('qsize', 0.05)
    ha2ev = kwargs.get('ha2ev', 27.2113961318)
    wn2ev = kwargs.get('wn2ev', 0.000123981)
    # wn2eh = kwargs.get('wn2eh', 0.00000455633)
    # ang2br = kwargs.get('ang2br', 1.889725989)
    # amu2me = kwargs.get('amu2me', 1822.88839)
    dipoles = kwargs.get('dipoles')
    diabatize = kwargs.get('diabatize')
    hessout = kwargs.get('hessout')

    # -------------------------------------------------------------------------
    # prepare various lists and dictionaries

    Heading, Params = [], []
    EH, SOC, Linear, Quadratic, Bilinear = [], [], [], [], []

    # prepare filenames
    zeroth_filename = f'{filnam}_refG.out'

    displacement_keys = ["+1", "+2", "-1", "-2"]

    # allocate list of dictionaries
    # displacement_filenames = [{} for k in range(1, nmodes+1)]

    # # build the file names
    # for kmode in range(1, nmodes + 1):  # store the filenames
    #     imode = modes_included[kmode]
    #     displacement_filenames[kmode-1] = {
    #         "+1": f'{filnam}_mode{imode}_+{qsize}.out',
    #         "+2": f'{filnam}_mode{imode}_+{qsize}x2.out',
    #         "-1": f'{filnam}_mode{imode}_-{qsize}.out',
    #         "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
    #         # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
    #         # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
    #         # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
    #         # "-2": f'{filnam}_mode{imode}_-{qsize}x2.out',
    #     }

    # if True:
    #     displacement_filenames = [
    #         f'{filnam}_mode{imode}_+{qsize}x2.out',
    #         f'{filnam}_mode{imode}_-{qsize}.out',
    #         f'{filnam}_mode{imode}_-{qsize}x2.out',
    #     ]
    # elif True:
    #     for k in displacement_keys:
    #         a = {"+1":'+', "+2":'+', "-1", "-2"}[k]
    #         b = {"+2":'x2', "-2": 'x2'}.get(k, '')

    #         # if '+' in k:
    #         #     a = k[0]
    #         #     b = '' if '1' in k else 'x2' if '2' in k else 'ERROR'

    #         displacement_filenames.append(
    #             f'{filnam}_mode{imode}_{a:s}{qsize}{b:s}.out'
    #         )

    # -------------------------------------------------------------------------
    # preparing the header of the `mctdh.op`

    Heading.append("OP_DEFINE-SECTION")
    Heading.append("title")

    Heading.append(f'{filnam} {nstate} states + ' + str(nmodes) + ' modes')
    Heading.append("end-title ")
    Heading.append("end-op_define-section\n")
    Heading.append("PARAMETER-SECTION\n")

    format_string = "{label:<25s}={value:>-15.9f}{units:>8s}\n"
    make_line = functools.partial(format_string.format, units=", ev")

    refG_exists = subprocess.call(["ls", zeroth_filename]) == 0

    if refG_exists:

        SOC_check = subprocess.run(['grep', "(DIRECT MAX.)", zeroth_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        SOC_flag = SOC_check.returncode == 0

        if SOC_flag:
            # Extract DSOME_cm_0
            DSOME_cm_0 = extract_DSOME(zeroth_filename, nstate)
            DSOME_cm_0_real, DSOME_cm_0_imag = DSOME_cm_0[0], DSOME_cm_0[1]

        with open("mctdh.op", "a") as mctdh_file:
            EH.append("#             Electronic Hamitonian             #\n")
            EH.append("# --------------------------------------------- #\n")

            GSE = extract_ground_state_energy(hessout, 'TOTAL ENERGY =')
            D1E = extract_diabatic_energy(zeroth_filename, f'STATE #.* 1.S GMC-PT-LEVEL DIABATIC ENERGY=')
            print(f'The ground state energy is: {GSE} Hartree\n')
            print(f'Diabat #1 energy is: {D1E} Hartree\n')
            linear_shift = (D1E - GSE) * ha2ev
            print(f'Linear shift value: {linear_shift} eV\n')

            for ist in range(1, nstate + 1):

                Ediab = refG_extract(zeroth_filename, f'STATE #.* {ist}.S GMC-PT-LEVEL DIABATIC ENERGY=')

                EH.append(make_line(label=f"EH_s{ist:>02d}_s{ist:>02d}", value=Ediab+linear_shift))

                # Extract coupling energy between state jst and ist
                for jst in range(1, ist):

                    Coup_ev = refG_extract(zeroth_filename, f'STATE #.* {jst} &.* {ist}.S GMC-PT-LEVEL COUPLING')

                    EH.append(make_line(label=f"EH_s{jst:>02d}_s{ist:>02d}", value=Coup_ev))

                EH.append("\n")

            EH.append(make_line(label=f"EH_s{nstate+1:>02d}_s{nstate+1:>02d}", value=0.0))

            EH.append("\n")

        with open('mctdh.op', 'w') as mctdh_file:

            for idx in Heading:
                mctdh_file.write(idx+'\n')

            for idx in EH:
                mctdh_file.write(idx)

    else:
        print(f"Skip extracting Hamiltonians from the non-existing {zeroth_filename}")

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

    Params.append(f"#                  Frequencies                  #\n")
    Params.append(f"# --------------------------------------------- #\n")
    Linear.append(f"#           Linear Coupling Constants           #\n")
    Linear.append(f"# --------------------------------------------- #\n")
    Quadratic.append(f"#         Quadratic Coupling Constants          #\n")
    Quadratic.append(f"# --------------------------------------------- #\n")
    Bilinear.append(f"#         Bilinear Coupling Constants           #\n")
    Bilinear.append(f"# --------------------------------------------- #\n")

    # header_string = "# {header:^40s} #\n"
    # spacer_string = "# " + '-'*40 + "#\n"
    # Bilinear = header_string.format("Bilinear Coupling Constants")

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
            grace_code[key] = subprocess.call(["grep", "DONE WITH MP2 ENERGY", displacement_filenames[key]])

        """ either of these work (logic wise)
            if any(code != 0 for code in grace_code.values()):
            if not all(code == 0 for code in grace_code.values()):
        """

        if not all(code == 0 for code in grace_code.values()):
            Params.append(f"not good to extract. Skipping mode {imode} for extracting vibronic couplings\n")

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

            grace_code_pp = subprocess.call(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_+{qsize}_mode{jmode}_+{qsize}.out"])
            grace_code_pm = subprocess.call(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_+{qsize}_mode{jmode}_-{qsize}.out"])
            grace_code_mp = subprocess.call(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_-{qsize}_mode{jmode}_+{qsize}.out"])
            grace_code_mm = subprocess.call(["grep", "DONE WITH MP2 ENERGY", f"{filnam}_mode{imode}_-{qsize}_mode{jmode}_-{qsize}.out"])

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


    Params.append("\n")
    Params.extend(Linear)
    Params.append("\n")
    Params.extend(Quadratic)
    Params.append("\n")
    Params.extend(Bilinear)
    Params.append("\n")
    Params.extend(SOC)
    Params.append("\n")

    with open("mctdh.op", "a") as mctdh_file:
        for idx in Params:
            mctdh_file.write(idx)

    # Open mctdh.op file for writing
    with open('mctdh.op', 'a') as mctdh_file:

        # Write the TDIPOLE block
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# ELECTRONIC TRANSITION DIPOLES\n")
        mctdh_file.write("-----------------------------------------\n")

        for ist in range(2, nstate + 1):
        #     for idx in range(0, 3):
        #         mctdh_file.write(f"IO_1_{ist} = {dipoles[ist][idx]}")
        #         mctdh_file.write("\n")

            for idx in range(0, 3):
                operate_lst = ["x", "y", "z"]
                mctdh_file.write(f"E{operate_lst[idx]}_s00_s{ist:>02d} = {dipoles[ist][idx]}")
                mctdh_file.write("\n")
            mctdh_file.write("\n")

        mctdh_file.write("end-parameter-section\n")
        # Write the header
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("HAMILTONIAN-SECTION\n")
        mctdh_file.write("-----------------------------------------\n")

        # Write modes and mode labels
        mctdh_file.write(" modes | el")
        for imode_include in range(1, nmodes + 1):
            mctdh_file.write(f" | v{modes_included[imode_include]:>02d}")
        mctdh_file.write("\n")
        mctdh_file.write("-----------------------------------------\n")

        # Write KINETIC OPERATOR FOR NORMAL MODES
        mctdh_file.write("# KINETIC OPERATOR FOR NORMAL MODES\n")
        mctdh_file.write("-----------------------------------------\n")
        for imode_include in range(1, nmodes + 1):
            mode_count = imode_include + 1
            mctdh_file.write(f"1.00*w{modes_included[imode_include]:>02d}   |{mode_count} KE\n")

        # Write HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# HARMONIC OSCILLATOR POTENTIALS FOR NORMAL MODES\n")
        mctdh_file.write("-----------------------------------------\n")
        for imode_include in range(1, nmodes + 1):
            mode_count = imode_include + 1
            mctdh_file.write(f"0.50*w{modes_included[imode_include]:>02d}   |{mode_count}  q^2\n")

        # Write ELECTRONIC COUPLING AT REFERENCE STRUCTURE
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# ELECTRONIC COUPLING AT REFERENCE STRUCTURE\n")
        mctdh_file.write("-----------------------------------------\n")
        for ist in range(1, nstate + 2):
            mctdh_file.write(f"EH_s{ist:>02d}_s{ist:>02d} |1 S{ist}&{ist}\n")
        mctdh_file.write("\n")
        for ist in range(1, nstate + 1):
            jlast = ist - 1
            for jst in range(1, jlast + 1):
                mctdh_file.write(f"EH_s{jst:>02d}_s{ist:>02d}  |1 S{jst}&{ist}\n")

        # Write LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# LINEAR DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            for ist in range(1, nstate + 1):
                mctdh_file.write(f"C1_s{ist:>02d}_s{ist:>02d}_v{imode:>02d} |1 S{ist}&{ist} |{kmode_count} q\n")

        # Write LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# LINEAR OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            for ist in range(1, nstate + 1):
                jlast = ist - 1
                for jst in range(1, jlast + 1):
                    mctdh_file.write(f"C1_s{jst:>02d}_s{ist:>02d}_v{imode:>02d} |1 S{jst}&{ist} |{kmode_count} q\n")

        # Write LINEAR AND QUADRATIC DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# QUADRATIC DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            for ist in range(1, nstate + 1):
                mctdh_file.write(f"0.50*C2_s{ist:>02d}s{ist:>02d}_v{imode:>02d}v{imode:>02d} |1 S{ist}&{ist} |{kmode_count} q^2\n")

        # Write LINEAR AND QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# QUADRATIC OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
        for kmode in range(1, nmodes + 1):
            imode = modes_included[kmode]
            kmode_count = kmode + 1
            for ist in range(1, nstate + 1):
                jlast = ist - 1
                for jst in range(1, jlast + 1):
                    mctdh_file.write(f"0.50*C2_s{jst:>02d}s{ist:>02d}_v{imode:>02d}v{imode:>02d} |1 S{jst}&{ist} |{kmode_count} q^2\n")

        # Write BILINEAR DIAGONAL VIBRONIC COUPLINGS
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# BILINEAR DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
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
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("# BILINEAR OFF-DIAGONAL VIBRONIC COUPLINGS\n")
        mctdh_file.write("-----------------------------------------\n")
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
            mctdh_file.write("-----------------------------------------\n")
            mctdh_file.write("# SOC FULL HAMILTONIAN SOC OFF-DIAGONAL VIBRONIC COUPLINGS\n")
            mctdh_file.write("-----------------------------------------\n")

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

        # Close the file
        mctdh_file.write("-----------------------------------------\n")
        mctdh_file.write("\nend-hamiltonian-section\n\n")

        # # Write TRANSITION DIPOLE HAMILTONIAN-SECTION_E*
        # for idx in range(0, 3):
        #     operate_lst = ["x", "y", "z"]
        #     mctdh_file.write(f"HAMILTONIAN-SECTION_E{operate_lst[idx]}\n")
        #     mctdh_file.write("\n")

        #     # Write modes and mode labels
        #     mctdh_file.write(" modes | el")
        #     for imode_include in range(1, nmodes + 1):
        #         mctdh_file.write(f" | m{modes_included[imode_include]}")
        #     mctdh_file.write("\n")
        #     mctdh_file.write("-----------------------------------------\n")

        #     for ist in range(2, nstate + 1):
        #         mctdh_file.write(f"E{operate_lst[idx]}_1_{ist}         |1 S1&{ist}")
        #         mctdh_file.write("\n")
        #     mctdh_file.write("\n")
        #     mctdh_file.write("\nend-hamiltonian-section\n\n")

        # Write HAMILTONIAN-SECTION_IO
        mctdh_file.write(f"HAMILTONIAN-SECTION_Ex\n")
        mctdh_file.write("\n")

        # Write modes and mode labels
        mctdh_file.write(" modes | el")
        for imode_include in range(1, nmodes + 1):
            mctdh_file.write(f" | v{modes_included[imode_include]:>02d}")
        mctdh_file.write("\n")
        mctdh_file.write("-----------------------------------------\n")

        for ist in range(1, nstate + 1):
            mctdh_file.write(f"1.0         |1 S{nstate + 1}&{ist}")
            mctdh_file.write("\n")
        mctdh_file.write("\n")
        mctdh_file.write("\nend-hamiltonian-section\n\n")

        # Close the file
        mctdh_file.write("end-operator\n")

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
            for i in range(3): # 3 for x,y,z
                mode_value_set.append(selected_lines[idx + i][20:])

    return mode_value_set


def get_number_of_atoms(hessout):
    """ Function to get the number of atoms from thFe hessout file """

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
    #     key = project_parameters.modes_included[i]
    #     freqcm[key] = frequences[i]

    # -------------------------------------------------------------------------
    if False and __debug__:   # print all frequencies
        string = "\n".join([f"frequency: {i} {freq_array[i]} CM-1" for i in range(ndim)])
        print(string)

    # return mode_array, freqcm
    return modes_array, freq_array


def compose_ref_structure(hessout, natoms):
    coord_lines = extract_lines_between_patterns(hessout, 'EQUILIBRIUM GEOMETRY LOCATED', 'INTERNUCLEAR DISTANCES')
    ref_file = 'ref_structure'

    if len(coord_lines) > 2:
        try:
            subprocess.run(['rm', '-f', ref_file])
        except Exception as e:
            print(f"Error deleting {ref_file}: {str(e)}")

        for atom_line in coord_lines[-natoms-1:-1]:
            with open(ref_file, "a") as python_ref_file:
                python_ref_file.write(atom_line)

        print(f'Sucessfully extracted equilibrium geometry from {hessout} and prepared ref_structure.')
    else:
        print(f'Unsucessful extraction of equilibrium geometry from {hessout}. Please prepare ref_structure manually.')

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
    grace_exists = subprocess.call(["grep", "DONE WITH MP2 ENERGY", output_filename]) == 0

    if grace_exists:
        print("Calculation at the reference structure has already been done.")
        return

    else:
        print("Run calculation at the undistorted reference structure")

        shutil.copy("temp.inp", f"{filnam}_refG.inp")

        with open(input_filename, "a") as inp_file:
            with open(refgeo, 'r', errors='replace') as ref_structure:
                inp_file.write(ref_structure.read())

        with open(input_filename, "a") as inp_file:
            inp_file.write(" $END\n")

        # Submit and run the refG calculation (you may need to customize this command based on your setup)
        # refG_job_result = subprocess.run(["./subgam.diab", input_filename, "4", "0", "1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os.system("sbatch" + ' -W' + " " + my_subgam(input_filename, ncpus=2, ngb=1, nhour=1))  # the wait

        # At this point, refG calculation has completed successfully.
        print("Calculation at the reference structure is done.")

    return

# ---------------------------------------------------------------------------------------


def main(ref_file="ref_structure", ncols=5, **kwargs):
    """ x """
    hessian_filename = kwargs['hessian_filename']
    _extract_freq_and_mode_from_hessian(hessian_filename)

    natoms = get_number_of_atoms(hessian_filename)

    ndim = natoms * 3
    nrmmod, freqcm = process_mode_freq(ndim, ncols)

    compose_ref_structure(hessian_filename, natoms)

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
        'qsize': project_parameters.qsize,
        'ha2ev': project_parameters.ha2ev,
        'wn2ev': project_parameters.wn2ev,
        'wn2eh': project_parameters.wn2eh,
        'ang2br': project_parameters.ang2br,
        'amu2me': project_parameters.amu2me
    })
    # name, modes = project_parameters.filnam, project_parameters.modes_included
    diabatize = diabatization(**diabatization_kwargs)
    print("Diabatization successfully modified"); return

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
        'qsize': project_parameters.qsize,
        'ha2ev': project_parameters.ha2ev,
        'wn2ev': project_parameters.wn2ev,
        'wn2eh': project_parameters.wn2eh,
        'ang2br': project_parameters.ang2br,
        'amu2me': project_parameters.amu2me,
        'nof_electronic_states': project_parameters.A,
        'ndim': ndim,
        'freqcm': freqcm,
        'dipoles': dipoles,
        'diabatize': diabatize,
        'hessout': kwargs['hessian_filename'],
    })

    # name, modes = project_parameters.filnam, project_parameters.modes_included
    mctdh(**mctdh_input_kwargs)
    print("mctdh successfully modified"); return

    print('The run was a success!')

    shutil.copy("mctdh.op", kwargs['project_filename'] + '.op')

    if (extra_debug := False):
        pprint.pprint(nrmmod)
        print('---------nrm mod done-----------')
        pprint.pprint(freqcm)
        print('---------freqcm done-----------')
        pprint.pprint(selected_lines)
        print('---------selected_lines done-----------')
        pprint.pprint(filtered_set)
        print('---------filtered_set done-----------')
        pprint.pprint(freq_value_set)
        print('---------freq_value_set done-----------')
        pprint.pprint(atmlst)
        print('---------atmlst done-----------')
        pprint.pprint(chrglst)
        print('---------chrglst done-----------')
        pprint.pprint(refcoord)
        print('---------refcoord done-----------')
        pprint.pprint(modes_included)
        print('---------modes included done-----------')

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
        'project_filename': project_parameters.filnam,
        'refG_in': f"{project_parameters.filnam}_refG.inp",  # reference geometry
        'refG_out': f"{project_parameters.filnam}_refG.out",  # reference geometry
        'modes_included': project_parameters.modes_included,
    })

    main(**kwargs)
