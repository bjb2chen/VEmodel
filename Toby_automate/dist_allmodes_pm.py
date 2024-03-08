""" x """

# -----------------------------------------------
# system imports
import os
from os.path import splitext, join, basename, abspath
import sys
import subprocess
from types import SimpleNamespace
import shutil
import re
import itertools as it
import functools
import cProfile
import pstats

# -----------------------------------------------
# third party imports
import numpy as np
# from numpy import float64 as F64
from numpy import complex128 as C128
import pprint


# -----------------------------------------------
# local packages
import project_parameters as pp
from project_parameters import *  # eventually remove this
SOC_flag = pp.SOC_flag  # initialize globally
# ---------------------------------------------------------------------------------------


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
# functions used multiple times throughout the code


def _reminder_produce_upper_triangle_indices(s):
    """ Assume a square matrix of shape (s, s).
    Prints the row, col indices that would index the upper triangle
    Just reminder code in case the (i,j) (a,b) or (row, col) gets confusing.
    """
    print("Indices using it.combinations")
    for row, col in it.combinations(range(s), 2):
        print(f"{row=} {col=}")

    # this is just a wrapper for it.combinations
    for row, col in upper_triangle_loop_indices(s, 2):
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
        return_obj = SimpleNamespace()
        return_obj.returncode = 0
        return return_obj
    else:
        return subprocess.call(command, **kwargs)


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
        return_obj = SimpleNamespace()
        return_obj.returncode = 0
        return return_obj
    else:
        return subprocess.run(command, **kwargs)


def _delete_file_using_rmrf(path):
    """ x """
    try:  # remove the previous file, as we are about to write to it
        subprocess_run_wrapper(['rm', '-f', path])
    except Exception as e:
        print(f"Error deleting {path}: {str(e)}")
        breakpoint()
    return


def _remove_existing_distorted_structure_files(filename_dict):
    """ Delete existing distorted structure files """
    for filename in filename_dict.values():
        _delete_file_using_rmrf(filename)
    return


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
    if True:
        assert number_of_indicies >= 2, "Don't use this fxn for 1D arrays / lists"
        return it.combinations(range(max_num), number_of_indicies)
    else:
        indices = []
        if number_of_indicies == 2:
            for n1 in range(max_num):
                for n2 in range(max_num):
                    if n1 < n2:
                        indices.append((n1, n2))

        elif number_of_indicies == 3:
            for n1 in range(max_num):
                for n2 in range(max_num):
                    for n3 in range(max_num):
                        if n1 < n2 and n2 < n3:
                            indices.append((n1, n2, n3))
        else:
            raise Exception(f"Not Implemented for {number_of_indicies=} > 3")

        return indices


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


def extract_lines_between_patterns(path, start_pattern, end_pattern, collecting=False):
    """ Function to extract lines between patterns in a file """
    selected_lines = []
    # would be good to replace this with memory mapping find or grep command?
    with open(path, 'r', errors='replace') as file:
        for line in file:
            if start_pattern in line:
                collecting = True
            if end_pattern in line:
                collecting = False
                # break  # pattern should only occur in file once
            if collecting:
                selected_lines.append(line)

    # del selected_lines[slice(0, nof_line_skip)]

    return selected_lines


def _make_displacement_filenames():
    """ x """

    # import certain constants from project_parameters
    from project_parameters import file_name, N
    qsize = pp.gamess_const.qsize

    # initialize dictionaries
    linear_displacement_filenames = {}
    bi_linear_displacement_filenames = {}

    for i in range(N):
        q1_label = pp.mode_map_dict[i]
        for key in linear_disp_keys:
            sign = key[0]  # select plus or minus

            order = int(key[1])
            max_order = pp.nof_displacements_per_mode[i]
            if not (order <= max_order):
                continue  # skip this combination

            name = f'{file_name}_{qsize}_{sign}x{order}q{q1_label}'

            linear_displacement_filenames[(key, i)] = name + '.out'

    # the order j-before-i is important to match toby's
    # `_mode8_+0.05_mode7` (larger number on the left style)
    for i, j in upper_triangle_loop_indices(N, 2):
        # the labels are reversed because of the file names being lower triangle
        # q1 mapped to j NOT MISTAKE (matches toby gamess output naming scheme)
        q1_label = pp.mode_map_dict[j]
        q2_label = pp.mode_map_dict[i]
        for d1, d2 in it.product(['+', '-'], ['+', '-']):

            key = d1+d2  # just in case
            assert key in bi_linear_disp_keys, f"{d1+d2=} is not in {bi_linear_disp_keys=}"

            order = 1
            name = f'{file_name}_{qsize}'
            name += f'_{d1}x{order}q{q1_label}'
            name += f'_{d2}x{order}q{q2_label}'

            bi_linear_displacement_filenames[(key, i, j)] = name + '.out'

    return linear_displacement_filenames, bi_linear_displacement_filenames


# ---------------------------------------------------------------------------------------
# constants used throughout the code
highest_order = max(pp.nof_displacements_per_mode)

#linear_disp_keys = ["+1", "+2", "+3", "+4", "-1", "-2", "-3", "-4"]
linear_disp_keys = [f"+{n}" for n in range(1, highest_order+1)]
linear_disp_keys += [f"-{n}" for n in range(1, highest_order+1)]

bi_linear_disp_keys = ["++", "+-", "-+", "--"]

linear_temp_suffix = {f"+{n}": 'p'*n for n in range(1, highest_order+1)}
linear_temp_suffix.update({f"-{n}": 'm'*n for n in range(1, highest_order+1)})

# these could be lists or dicts... for now I'm just making them dicts
linear_temp_struct_filenames = {
    k: f'dist_structure_{linear_temp_suffix[k]}'
    for k in linear_disp_keys
}

bi_linear_temp_suffix = {
    "++": 'pp',
    "+-": 'pm',
    "-+": 'mp',
    "--": 'mm',
}
bi_linear_temp_struct_filenames = {
    k: f'dist_structure_{bi_linear_temp_suffix[k]}'
    for k in bi_linear_disp_keys
}

linear_displacement_filenames, bi_linear_displacement_filenames = _make_displacement_filenames()

# ---------------------------------------------------------------------------------------
# functions using grep to grab specific lines in files generated by GAMESS, Gaussian, MCTDH, etc.


def _extract_energy_from_gamessoutput_grep(file_path, pattern, column_specification_string, backup_line_idx):
    """
    Sometimes the output from a GAMESS calculatino is malformed and the file is a binary file instead of text.
    To make sure `grep` can still read this we need to use the `-a` command which does: (Process a binary file as if it were text)
    https://stackoverflow.com/questions/23512852/grep-binary-file-matches-how-to-get-normal-grep-output
    """
    try:
        # Use subprocess.run with the direct command
        command = f'grep -a "{pattern}" "{file_path}" | {column_specification_string}'
        result = subprocess_run_wrapper(command, shell=True, text=True, capture_output=True)
        assert result.stderr == "", f"Something wrong? {result.stderr=}\n{result}\n"

        try:
            # If there is output, convert it to float
            output = float(result.stdout.strip().replace(" ", ""))
            return output
        except Exception as e:
            if True:  # try to find out reason for grep failing
                print("Grep failed?!")
                print("Command\n", command)
                print("Result\n", result)
                print("Please try grep command manually before proceeding")
                breakpoint()

            with open(file_path, 'r', errors='replace') as file:
                for line in reversed(file.readlines()):
                    match = re.search(pattern, line)
                    if match:
                        return float(line[backup_line_idx].strip().replace(" ", ""))

    except subprocess.CalledProcessError:
        # Return None if there is an error
        return None


def extract_in_Hartrees(file_path, pattern):
    """ This picks up the energy in Hartrees.

    Lines in the file might look like this:

        --- DIABATIC ENERGIES (DIAGONAL ELEMENT) ---       HARTREE            EV
        STATE #  1'S GMC-PT-LEVEL DIABATIC ENERGY=     -75.798270171       0.000000000
        STATE #  2'S GMC-PT-LEVEL DIABATIC ENERGY=     -75.708488402       2.443086361

    This function gets the numbers from the HARTREE column
    """
    column_specification_string = "tail -1 | cut -c44-61"
    backup_line_idx = slice(44, 62)  # equivalent to `array[62:]`
    return _extract_energy_from_gamessoutput_grep(
        file_path, pattern,
        column_specification_string,
        backup_line_idx
    )


def extract_diabatic_energy(file_path, pattern):
    """ This picks up the energy in Hartrees. """
    return extract_in_Hartrees(file_path, pattern)


def extract_in_eV(file_path, pattern):
    """ This picks up the energy in eV.

    Lines in the file might look like this:

        --- DIABATIC ENERGIES (DIAGONAL ELEMENT) ---       HARTREE            EV
        STATE #  1'S GMC-PT-LEVEL DIABATIC ENERGY=     -75.798270171       0.000000000
        STATE #  2'S GMC-PT-LEVEL DIABATIC ENERGY=     -75.708488402       2.443086361

    This function gets the numbers from the EV column
    """
    if False:
        return _extract_energy_from_gamessoutput_memap(file_path, pattern)
    else:
        column_specification_string = "tail -1 | cut -c62-"
        backup_line_idx = slice(62, None)  # equivalent to `array[62:]`
        return _extract_energy_from_gamessoutput_grep(
            file_path, pattern,
            column_specification_string,
            backup_line_idx
        )


def refG_extract(file_path, pattern):
    """ This picks up the energy in eV. """
    return extract_in_eV(file_path, pattern)


def extract_ground_state_energy(hessian_path, pattern):
    """ This picks up the total ground-state-energy (GSE) in eV.
    Relevant lines in the file look like this:

             WAVEFUNCTION NORMALIZATION =       1.0000000000

                    ONE ELECTRON ENERGY =    -122.8353239571
                    TWO ELECTRON ENERGY =      37.4785649796
               NUCLEAR REPULSION ENERGY =       9.1280924730
                                          ------------------
                           TOTAL ENERGY =     -76.2286665045

     ELECTRON-ELECTRON POTENTIAL ENERGY =      37.4785649796
      NUCLEUS-ELECTRON POTENTIAL ENERGY =    -198.9743416468
       NUCLEUS-NUCLEUS POTENTIAL ENERGY =       9.1280924730
                                          ------------------
                 TOTAL POTENTIAL ENERGY =    -152.3676841942
                   TOTAL KINETIC ENERGY =      76.1390176897
                     VIRIAL RATIO (V/T) =       2.0011774359

    This function gets the value from the `TOTAL ENERGY =` line.
    """
    column_specification_string = "tail -1 | cut -c40-"
    backup_line_idx = slice(40, None)  # equivalent to `array[40:]`
    return _extract_energy_from_gamessoutput_grep(
        hessian_path, pattern,
        column_specification_string,
        backup_line_idx
    )


def extract_DSOME(path, nof_states, nof_electron_couplings=2):
    """ This extracts the off-diagonal spin-orbit-couplings.
    Relevant lines in the file look like this:

        HSO MATRIX IN DIABATIC REPRESENTATION (DIRECT MAXIMIZATION)

         -------------------------------------------
          -  DIABATIC SPIN-ORBIT MATRIX ELEMENTS -
         -------------------------------------------

         --- SOC COUPLINGS (OFF DIAGONAL ELEMENT) ---                  CM-1
         STATE #  1 &  2'S GMC-PT-LEVEL 1-El COUPLING  =   -3.234551 +    -7.552143I
         < ...... >
         STATE #  5 &  6'S GMC-PT-LEVEL 2-El COUPLING  =    5.784123 +    -0.552264I

        SOC EIG. VALUES and VECTORS IN DIABATS (DIRECT MAX.)

    This picks up the two numbers after `=` (-3.234551, -7.552143I) and matches them with the numbers (1 & 2).

    Since the couplings occur between electronic states then we know that there will be
        1->2, 1->3, etc. till 5->6  (for A=6)
    1->(2,3,4,5,6) 5 transitions
    2->(3,4,5,6) 4 transitions
    3->(4,5,6) 3 transitions
    4->(5,6) 2 transitions
    5->(6) 1 transitions
    so therefore we will have 5+4+3+2+1 = 15 transitions  (this is know as the n-th triangular number)
    BUT! since we have 1-electron and 2-electron couplings we end up with 15 * 2 = 30 transitions
    """
    max_nof_electron_couplings = 2  # FIXED VALUE, we can only have 1 or 2 electron couplings
    nof_transitions = np.sum(range(A)) * max_nof_electron_couplings

    def get_line_list():
        start_pattern = "HSO MATRIX IN DIABATIC REPRESENTATION (DIRECT MAXIMIZATION)"
        end_pattern = 'SOC EIG. VALUES and VECTORS IN DIABATS (DIRECT MAX.)'

        sed_command = f"sed -n '/{start_pattern}/,/{end_pattern}/p' {path}"
        result = subprocess_run_wrapper(sed_command, shell=True, text=True, capture_output=True)
        # i believe there is another way to check if the result is empty?
        # I think the `result` object has error codes or some such?

        line_list = result.stdout.splitlines()

        # check if extraction worked
        if not (isinstance(line_list, list) and line_list != []):
            print('Cannot use sed to extract')
            line_list = extract_lines_between_patterns(  # fallback to slow pythonic extraction
                f'{path}',
                f'{start_pattern}',
                f'{end_pattern}',
            )
            print(f'Using selected lines from {path}, opened via python')

        if False and __debug__:
            for i, l in enumerate(line_list): print(i, l); breakpoint()

        # remove unnecessary headers (first 7 lines and last 2 lines)
        s, e = 7, -2
        reduced_line_list = line_list[s:e]

        if True and __debug__:  # old style line list from sed
            assert 'STATE #' not in line_list[s-1]  # lines we don't need
            assert 'STATE #' in line_list[s]
            assert 'STATE #' in line_list[e-1]
            assert 'STATE #' not in line_list[e]  # lines we don't need

        return reduced_line_list

    selected_lines = get_line_list()
    assert len(selected_lines) == nof_transitions, (
        f"Why are there {len(selected_lines)} lines?"
        f"There should be {nof_transitions}"
        "Check\n" + "\n".join(selected_lines)
    )

    # check for stars in lines? this means Gamess calculation failed/is bad
    for line in selected_lines:
        assert '*' not in line, "You messed up! Gamess calculated failed?"

    shape = (A, A)
    spin_orbit_array = np.zeros(shape, dtype=C128)
    # spin_orbit_dict = {}

    def old_process_line(line):
        """ """
        ist = line[9:12].strip().replace(" ", "")
        jst = line[14:16].strip().replace(" ", "")
        real = line[48:61].strip().replace(" ", "")
        if '*' in real:
            real = 0

        imaginary = line[63:75].strip().replace(" ", "")
        if '*' in imaginary:
            imaginary = 0

        key = f"{ist} & {jst},{line[31:33]}"
        spin_orbit_dict[key] = complex(float(real), float(imaginary))
        return

    def new_process_lines(spin_orbit_dict, line1, line2):
        """
        a line looks like
            STATE #  5 &  6'S GMC-PT-LEVEL 2-El COUPLING  =    2.615044 +     0.462500I
        line.strip().split() looks like
            ['STATE', '#', '5', '&', "6'S", 'GMC-PT-LEVEL', '2-El', 'COUPLING', '=', '2.615044', '+', '0.462500I']
        line.strip().split('=')
            ["STATE #  5 &  6'S GMC-PT-LEVEL 2-El COUPLING  ", '    2.615044 +     0.462500I']
        line.strip().replace(' ', '').split('=')
            ["STATE#1&2'SGMC-PT-LEVEL1-ElCOUPLING", '-2.931761+-7.805363I']
        """
        label1, values1 = line1.strip().replace(' ', '').split('=')
        label2, values2 = line2.strip().replace(' ', '').split('=')

        # extract values ('-2.931761+-7.805363I')
        # real1, imag1 = values1.replace('I', 'j').split('+')
        # real2, imag2 = values2.replace('I', 'j').split('+')

        """ python uses `j` instead of `i` or `I` to denote imaginary #
        see https://stackoverflow.com/questions/40421736/python-numpy-convert-string-to-complex-number

        If we have negative imagine number then we get +-7.805363I since there is ALWAYS a `+` in the middle column,
        so replace it with just a `-` (since python's `complex()` can't handle `+-` )
        """
        # complex_string1 = values1.replace('I', 'j').replace('+-', '-')
        # complex_string2 = values2.replace('I', 'j').replace('+-', '-')
        number1 = complex(values1.replace('I', 'j').replace('+-', '-'))
        number2 = complex(values2.replace('I', 'j').replace('+-', '-'))

        # extract label ("STATE#1&2'SGMC-PT-LEVEL1-ElCOUPLING")
        a, b = map(int, label1[6:].split("'S")[0].split('&'))
        a2, b2 = map(int, label2[6:].split("'S")[0].split('&'))
        assert (a == a2) and (b == b2), f"The two lines have different a's and b's"

        spin_orbit_array[a-1, b-1] = number1 + number2
        return

    for i in range(0, nof_transitions, 2):
        line1, line2 = selected_lines[i:i+2]  # pick up two lines (1 & 2 electron pairs)
        try:
            # old_process_line(line1)
            # old_process_line(line2)
            new_process_lines(spin_orbit_array, line1, line2)
        except Exception as e:
            print(f"Error processing lines: {str(e)}\n {line1}\n {line2}")
            breakpoint()

    return spin_orbit_array

def extract_fitting_parameters(path, pattern='a.*='):
    return

# ---------------------------------------------------------------------------------------
# from memorymap_extract import extract_string_list, extract_from_file, find_byte_begin_and_end

def _example_processing_function(path, memmap):

    # if processing block 1
    begin_string, end_string = "block_1_begin", "block_1_end"
    # here is where you call extract string
    lines = extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=1)

    def process_block_1(): return
    b1_out = process_block_1(lines)

    begin_string, end_string = "block_2_begin", "block_2_end"
    # here is where you call extract string
    lines = extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=3)

    def process_block_2(): return
    b2_out = process_block_2(lines)

    stuff = [b1_out, b2_out, ]  # whatever you might need to return?

    return stuff


def process_func_1(path, memmap, pattern):

    if 'GMC-PT-LEVEL DIABATIC ENERGY' in pattern:
        ''' "STATE #.* 1.S GMC-PT-LEVEL DIABATIC ENERGY= '''
        begin_string = " - DM DIABT PT HAMILTONIAN MATRIX ELEMENTS -"
        end_string = " --- DIABATIC COUPLINGS (OFF DIAGONAL ELEMENTS)---"
        lines = extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=4)

        print("Ingested:")
        for i, l in enumerate(lines): print(f"Line {i:02d}: ", l)

        hartree_list = [float(l[6]) for l in lines]
        hartree_array = np.array(hartree_list)
        if False: print("Hartree array\n", hartree_array)

        eV_list = [float(l[7]) for l in lines]
        eV_array = np.array(eV_list)
        if False: print("eV array\n", eV_array)

        n = int(pattern.replace('.S GMC-PT-LEVEL DIABATIC ENERGY=', '').replace('STATE #.* ', ''))
        return eV_array[n-1]

    elif 'GMC-PT-LEVEL COUPLING' in pattern:
        ''' STATE #.* 1 &.* 2.S GMC-PT-LEVEL COUPLING '''
        begin_string = " --- DIABATIC COUPLINGS (OFF DIAGONAL ELEMENTS)---  HATREE            EV"
        end_string = "     - ROTATION TO DIABATIC STATES -    "
        lines = extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=1)
        breakpoint()

        print("Ingested:\n")
        for i, l in enumerate(lines): print(f"Line {i:02d}: ", l)

        hartree_list = [float(l[7]) for l in lines]
        hartree_array = np.array(hartree_list)
        print("Hartree array\n", hartree_array)

        eV_list = [float(l[8]) for l in lines]
        eV_array = np.array(eV_list)
        print("eV array\n", eV_array)

        breakpoint()
        s1, s2 = map(int, pattern.replace('.S GMC-PT-LEVEL COUPLING', '').replace('STATE #.*', '').split('&.*'))

        # match the s1/s2 with the appropriate line in `lines`
        return eV_array[n-1]

    else:
        print(pattern)
        breakpoint()

    # here is where you call extract string
    # lines = extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=1)

    return


def _extract_energy_from_gamessoutput_memap(path, pattern):
    """ use memory map to extract data  """
    string = extract_from_file(path, process_func_1, pattern)
    return string

# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------


def diabatization(**kwargs):
    """ Preform diabatization?

    Returns a dictionary of dictionaries `dist_coord`
    Whose keys are `displacement_keys` and `bi_linear_keys`.
    Each of those individual dictionaries store ( ... ?)

    (other info)
    """

    # -------------------------------------------------------------------------
    # unpacking kwargs
    all_frequencies_cm = kwargs['all_freq']  # N_tot
    all_normal_modes = kwargs['normal_modes']  # Z*3 x N_tot

    refcoord = kwargs['reference_coordinates']  # Z*3
    ref_coord_array = np.array([*refcoord.values()])

    atmlst = kwargs['atom_dict']
    atom_list = [*atmlst.values()]

    chglst = kwargs['charge_dict']
    charge_list = [*chglst.values()]

    # -------------------------------------------------------------------------
    # extract only the frequencies of the modes we selected
    s_idx = [i-1 for i in pp.selected_mode_list]  # need to shift by 1 for 0-indexed arrays

    freq_array = all_frequencies_cm[s_idx]  # s_idx=[6,7,8] returns 7,8,9th freq values
    mode_array = all_normal_modes[:, s_idx]  # return 7,8,9th columns of Z*3 x N_tot array, assume we return array

    """ Note the indexing here!!
        For a system with 9 modes (and therefore 3 atoms, since Z * 3 = 9 = N_tot)
        there will be 9 frequencies `len(all_frequencies_cm)` = 9
        So if `selected_mode_list` is [7, 8, 9]
        then `s_idx` will be [6, 7, 8]
        since `all_frequencies_cm` is a length 9 array indexed at 0, `[6, 7, 8]` will return the 7th, 8th, 9th elements
    """
    # -------------------------------------------------------------------------
    # extract constants
    wn2eh = pp.QM_const.wn2eh
    ang2br = pp.QM_const.ang2br
    amu2me = pp.QM_const.amu2me

    qsize = pp.gamess_const.qsize

    # -------------------------------------------------------------------------
    # prepare various lists and dictionaries
    file_name = pp.file_name  # this is the root file/calculation name for all files (displacements etc)

    # -------------------------------------------------------------------------

    # prepare the displaced/distorted co-ordinates (the object we return)
    disp_coord = {}
    disp_coord.update({k: {} for k in linear_disp_keys})
    disp_coord.update({k: {} for k in bi_linear_disp_keys})

    # -------------------------------------------------------------------------
    # NEIL's mappings / definitions
    from project_parameters import A, N, Z
    nof_surfaces = A  # alias
    nof_modes = N  # alias
    nof_atoms = Z  # alias

    NEW = np.newaxis  # alias

    # -------------------------------------------------------------------------
    if ((Z*3)**2 > 1e4):  # check in case precomputing will take too much memory
        print("Maybe think about not precomputing? Check with Neil again?")
        breakpoint()
        import sys; sys.exit()
        # ---------------------------------------------------------------------
        # may want to consider indexing the array?

        # throw away the modes we don't need
        column_index = [j-1 for j in pp.selected_mode_list]

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
    # precompute_bilinear = True
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
        # the `:` of the second dimension of mode_array is what we index with i or j
        """
        Z, N = pp.Z, pp.N

        for i in range(N):
            assert pp.nof_displacements_per_mode[i] >= 2, f"Mode {i+1} has order below 2, change code"

        # assume all modes have at least 2
        displacements = {
            "+1": reference[:, NEW] + 1.0 * R_array[NEW, :] * mode_array[:, :],
            "-1": reference[:, NEW] - 1.0 * R_array[NEW, :] * mode_array[:, :],
            "+2": reference[:, NEW] + 2.0 * R_array[NEW, :] * mode_array[:, :],
            "-2": reference[:, NEW] - 2.0 * R_array[NEW, :] * mode_array[:, :],
        }
        #                  (Z*3, 1)                 (1,  N)            (Z*3, N)
        #                  (9, 1)                   (1,  3)             (9, 3)
        #                  (12, 1)                  (1,  6)            (12, 6)

        # Do a loop here so that it scales for 8 points, 10 points
        # Change the filenames such that for 8 points, do x8 so that the filenames are not too long

        for i in range(N):
            max_order_of_qi = pp.nof_displacements_per_mode[i]
            for order in range(3, max_order_of_qi+1):
                displacements.update({
                    f"+{order}": reference[:] + order * R_array[i] * mode_array[:, i],
                    f"-{order}": reference[:] - order * R_array[i] * mode_array[:, i],
                })

        assert set(displacements.keys()) == set(linear_disp_keys), (f"{displacements.keys()=}\nand\n{linear_disp_keys=}\n no longer agree!")

        shape = (Z*3, N)

        for k in linear_disp_keys:
            if int(k[1]) > 2:
                assert displacements[k].shape == (Z*3,),  f"{k=} {displacements[k].shape=} not {(Z*3,)=}?"
            else:
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
        `displaced_q` is an array of dimension (ndim, N)

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
        Z = pp.Z  # number of atoms

        header, data = "{:<2s} {:} ", "{: .10g} " * 3  # x,y,z components
        template_string = header + data

        for key in key_list:
            file_contents = []
            d = displaced_q[key]
            for idx_atom in range(Z):
                offset = 3*idx_atom  # 3 xyz-coordinates

                if False and len(mode_idx) > 1:
                    print(*mode_idx, d.shape); breakpoint()

                # temporary, if we are doing extra displacements along 1 mode for linear
                # then the dimensionality is different
                if (d.ndim == 1) and ('+3' in key_list) and int(key[1]) > 2:
                    string = template_string.format(
                        atom_list[idx_atom], charge_list[idx_atom],
                        d[offset+0],  # x component
                        d[offset+1],  # y component
                        d[offset+2],  # z component
                    )
                else:
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
        """ this could be improved somewhat """

        # Check if the reference geometry calculation is done?
        refG_out = f"{filename}_refG.out"
        grace0 = subprocess_run_wrapper(["grep", "DONE WITH MP2 ENERGY", refG_out])
        ref_geom_flag_exists = bool(grace0.returncode == 0)

        for key in linear_disp_keys:

            order = int(key[1])
            max_order = pp.nof_displacements_per_mode[i]
            if not (order <= max_order):
                continue  # skip this combination

            # get the filename WITHOUT the '.out' ending
            games_filename = linear_displacement_filenames[(key, i)].replace('.out', '')

            shutil.copy('temp.inp', games_filename+'.inp')

            # so you only write to this file to then change another file?
            # but i guess its good to have a record of these files?
            distored_struct_file = linear_temp_struct_filenames[key]
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
                    output_filename = my_subgam(games_filename+'.inp', ncpus=2, ngb=1, nhour=1)
                    os_system_wrapper(f"sbatch {output_filename}")
                except Exception as e:
                    print(f"Error running diabatization calculation: {str(e)}")
            else:
                print(f"{games_filename} is done\n")

        return

    def _create_bilinear_diabatization_input_files(j, i, filename, qsize):
        """ this could be improved somewhat """

        # Check if the reference geometry calculation is done?
        refG_out = f"{filename}_refG.out"
        grace0 = subprocess_run_wrapper(["grep", "DONE WITH MP2 ENERGY", refG_out])
        ref_geom_flag_exists = bool(grace0.returncode == 0)

        for key in bi_linear_disp_keys:

            # get the filename WITHOUT the '.out' ending
            games_filename = bi_linear_displacement_filenames[(key, j, i)].replace('.out', '')

            shutil.copy('temp.inp', games_filename+'.inp')

            # so you only write to this file to then change another file?
            # but i guess its good to have a record of these files?
            distored_struct_file = bi_linear_temp_struct_filenames[key]
            with open(distored_struct_file, 'r', errors='replace') as fp:
                data = fp.read()

            with open(games_filename+'.inp', 'a') as fp:
                fp.write(data)  # can you just do data + ' $END' in one write?
                fp.write('\n $END')

            # Check if the calculation is done already
            grace2 = subprocess_run_wrapper(["grep", "DONE WITH MP2 ENERGY", games_filename+'.out'])
            gamess_calculation_not_run = bool(grace2.returncode != 0)

            # this will never work? grace0 is not defined
            if (ref_geom_flag_exists and gamess_calculation_not_run) or pp.dry_run:
                print(f"Running calculations for {games_filename}!")
                try:
                    output_filename = my_subgam(games_filename+'.inp', ncpus=2, ngb=1, nhour=1)
                    os_system_wrapper(f"sbatch {output_filename}")
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

        _remove_existing_distorted_structure_files(linear_temp_struct_filenames)
        index = (i, )
        _save_distorted_structure(
            index, disp_coord, charge_list, atom_list,
            linear_temp_struct_filenames,
            linear_disp_keys
        )
        _create_linear_diabatization_input_files(i, file_name, qsize)

        # 2D distortion to get bilinear vibronic coupling
        for j in range(0, i):

            # if not precompute_bilinear:  # modify disp_coord in-place
            # _compute_bi_linear_displacements(i, j, rsize, mode_array, disp_coord)

            _remove_existing_distorted_structure_files(bi_linear_temp_struct_filenames)
            # index = (i, j) if precompute_bilinear else (j, )
            index = (i, j)
            _save_distorted_structure(
                index, disp_coord, charge_list, atom_list,
                bi_linear_temp_struct_filenames,
                bi_linear_disp_keys,
            )
            _create_bilinear_diabatization_input_files(j, i, file_name, qsize)

    return disp_coord


# ---------------------------------------------------------------------------------------

# this function is not used anywhere? can probably delete?
def extract_coupling_energy(file_path, pattern):
    """ this picks up the energy in eV (same as refG_extract) """
    return extract_in_eV(file_path, pattern)


# this function is not used anywhere? can probably delete?
def find_nstate(file_path, pattern='# of states in CI      = '):
    with open(file_path, 'r', errors='replace') as file:
        for line in file:
            if pattern in line:
                return int(line.split('=')[1].strip())
    return None  # Return None if the pattern is not found

# ---------------------------------------------------------------------------------------


def fitting():
    """
    For nof_displacements_per_mode > 2, want to extract data for fitting.

    For A == 1:                                                     (extract this)

         1     E(REF-CI)=      -56.2499276613     E(GMC-PT2)=      -56.3710474990

    For A > 1:                                     (extract this)

    --- DIABATIC ENERGIES (DIAGONAL ELEMENT) ---       HARTREE            EV
    STATE #  1'S GMC-CI-LEVEL DIABATIC ENERGY=     -55.916713021       0.000000000
    STATE #  2'S GMC-CI-LEVEL DIABATIC ENERGY=     -55.694114507       6.057214072
    STATE #  3'S GMC-CI-LEVEL DIABATIC ENERGY=     -55.694113719       6.057235537

    gnuplot command:
    p 'fitting_3st_mode7.dat' u 1:2 w lp, "" u 1:3 w lp, "" u 1:4 w lp
    """
    A, N = pp.A, pp.N  # to stop flake8 complaining

    # these lines weren't used?
    # format_string_au = "{header:s}{xvals:<6s}{yvals:<-10.10f}{units:<15s}{file:<60s}\n"  # truncated to 10 d.p.
    # format_string_ev = "{header:s}{xvals:<6s}{yvals:<-10.7f}{units:<15s}{file:<60s}\n"  # truncated to 7 d.p.
    # make_line_au = functools.partial(format_string_au.format, units=", Hartree")
    # make_line_ev = functools.partial(format_string_ev.format, units=", eV")
    # xtemp, ytemp = "{:<4s} ", "{: 10g} " * A

    shape = (A, A)

    header_1 = ' qsize  | {0} | filepath\n'.format("Hamiltonian value units (A.U.) Hartree")
    header_2 = ' qsize  | {0} | filepath\n'.format("other columns on this side are Hamiltonian values, each col is one state, state 1 2 -> A, units eV")

    a_pattern = 'STATE #.* {col}.S GMC-PT-LEVEL DIABATIC ENERGY='

    ref_geom_path = f'{pp.file_name}_refG.out'

    d1_row_template = " {:+.2f} " + "   {: .10f}   " + "{path}\n"
    d2_row_template = " {:+.2f} " + "   {: .10f}   " * A + "{path}\n"

    qsize = pp.gamess_const.qsize
    ha2ev = pp.QM_const.ha2ev

    # initial/reference geom rows (first row)?
    # (Pdb) d2_row_template
    # ' {:+.2f}{: .10f}{: .10f}{: .10f}{path}\n'

    def _max_order(i):
        return pp.nof_displacements_per_mode[i]

    large_N = [i for i in range(N) if _max_order(i) > 2]

    # ---------------------------------------------------------------------
    def _extract_single_surface():
        """
            for every list there is a list
            each element of the list is a row in the file contents
            you can do "\n".join(file_row_list[i]) to get a single string which
            is the file contents of the i'th mode
        """
        file_row_list = [[] for i in range(N)]
        fitting = {}

        for i in large_N: # skip routine linear and quadratic only modes
            # ---------------------------------------------------------------------------
            refG_key = ('refG', i)  # special case?
            column_specification_string = "tail -1 | cut -c62-"
            backup_line_idx = slice(62, None)
            fitting[refG_key] = _extract_energy_from_gamessoutput_grep(
                ref_geom_path, 'E(GMC-PT2)',
                column_specification_string,
                backup_line_idx
            )
            file_row_list[i].append(
                d1_row_template.format(0, fitting[refG_key], path=ref_geom_path)
            )

            # ---------------------------------------------------------------------------
            # add remaining rows
            for key in linear_disp_keys:
                order = int(key[1])

                if not (order <= max_order):
                    # skip this combination when order > max_order
                    # as it does not exist e.g. [8, 3, 2 ...] means key('+4', 1) not exist
                    continue

                assert (max_order > 2), f"How is {max_order=} < 3?"

                column_specification_string = "tail -1 | cut -c62-"
                backup_line_idx = slice(62, None)

                fitting[(key, i)] = _extract_energy_from_gamessoutput_grep(
                    linear_displacement_filenames[(key, i)], 'E(GMC-PT2)',
                    column_specification_string,
                    backup_line_idx
                )

                x = int(key) * qsize
                file_row_list[i].append(
                    d1_row_template.format(x, fitting[(key, i)], path=linear_displacement_filenames[(key, i)])
                )

        return file_row_list, fitting

    def _extract_multi_surface():
        """ """
        file_row_list = [[] for i in range(N)]
        fitting = {}

        for i in large_N:  # skip routine linear and quadratic only modes
            # ---------------------------------------------------------------------------
            refG_y_vals = [0.0, ]*A
            file_row_list[i].append(
                d2_row_template.format(0.0, *refG_y_vals, path=ref_geom_path)
            )

            # ---------------------------------------------------------------------------
            max_order = pp.nof_displacements_per_mode[i]
            # add remaining rows
            for key in linear_disp_keys:
                order = int(key[1])

                if not (order <= max_order):
                    # skip this combination when order > max_order
                    # as it does not exist e.g. [8, 3, 2 ...] means key('+4', 1) not exist
                    continue

                assert (max_order > 2), f"How is {max_order=} < 3?"

                # E(GMC-PT2)
                array = np.zeros(shape)
                for a in range(A):

                    refG_fitting = extract_in_Hartrees(ref_geom_path, a_pattern.format(col=a+1))
                    fitting[(key, i, a+1)] = extract_in_Hartrees(linear_displacement_filenames[(key, i)], a_pattern.format(col=a+1))

                    # convert to eV
                    array[a, a] = (fitting[(key, i, a+1)] - refG_fitting) * ha2ev

                x = int(key) * qsize
                # y_vals = [array[a, a] for a in range(A)]
                y_vals = np.diag(array)  # there is numpy method for this already
                file_row_list[i].append(
                    d2_row_template.format(x, *y_vals, path=linear_displacement_filenames[(key, i)])
                )
            #
        return file_row_list, fitting

    if A == 1:
        row_list, fitting_dict = _extract_single_surface()
        header = header_1
    elif A > 1:
        row_list, fitting_dict = _extract_multi_surface()
        header = header_2
    else:
        raise Exception(f"How is {A=} < 1 ?")

    # ---------------------------------------------------------------------
    # sort
    for i in large_N:  # skip routine linear and quadratic only modes
        row_array = np.array(row_list[i])

        # grab the value of the first column and sort along it
        n_list = [float(row.split()[0]) for row in row_array]
        idx = np.argsort(n_list)

        sorted_rows = row_array[idx]  # -0.4 -> +0.4
        if False: sorted_rows = np.flip(sorted_rows)  # this gives +0.4 -> -0.4

        row_list[i] = sorted_rows.tolist()

    # ---------------------------------------------------------------------
    # save
    for i in large_N:  # skip routine linear and quadratic only modes

        string = header + "".join(row_list[i])
        path = f'fitting_{A}st_mode{pp.mode_map_dict[i]}.dat'
        with open(path, 'w') as fp:
            fp.write(string)

    # ---------------------------------------------------------------------
    # solve for the polynomial coefficients

    def _perform_fitting(path, i):
        """ This uses gnuplot to print out a,b,c's in the gnuplot logfile.
        See - http://www.bersch.net/gnuplot-doc/fit.html
        f(x) = a + b*x + c*x**2
        fit f(x) 'measured.dat' using 1:2 via a,b,c
        plot 'measured.dat' u 1:2, f(x)
        """
        max_order = pp.nof_displacements_per_mode[i]
        size = [1200, 800]
        fitting_coeff = ''.join([f"+a{i}*x**{i}" for i in range(1, max_order+1)])
        fitting_coeff_lst = ''.join([f",a{i}" for i in range(1, max_order+1)])

        plotting_command = '\n'.join([
            f"set terminal png size {size[0]},{size[1]}",
            f"set output '{path}.png'",
            f"set fit logfile '{path}_FIT.log'",
            f"f(x)=a0{fitting_coeff}",
            f"fit f(x) '{path}.dat' u 1:2 via a0{fitting_coeff_lst}",
            f"plot '{path}.dat' u 1:2 w p, f(x)",
        ])

        with open(f'{path}.log', 'w') as fp:
            fp.write(plotting_command)

        subprocess.run(['gnuplot', f'{path}.log'])
        # subprocess_run_wrapper(['gnuplot', f'{path}.log'])
        return

    for i in large_N:  # skip routine linear and quadratic only modes
        path = f'fitting_{A}st_mode{pp.mode_map_dict[i]}'
        _perform_fitting(path, i)

    # ---------------------------------------------------------------------
    # extract the polynomials coefficients and ....

        def _something(path, i):
            fitting_params = {}
            max_order = pp.nof_displacements_per_mode[i]
            for a in range(max_order+1):
                column_specification_string = "tail -1 | cut -c18-35"
                backup_line_idx = slice(18, 36)
                fitting_params[f'a{a}'] = _extract_energy_from_gamessoutput_grep(
                    f'{path}_FIT.log', 'a.*=',
                    column_specification_string,
                    backup_line_idx
                )

            pprint.pprint(fitting_params)

    for i in large_N:  # skip routine linear and quadratic only modes
        path = f'fitting_{A}st_mode{pp.mode_map_dict[i]}'
        _something(path, i)

    # do something
    # ---------------------------------------------------------------------

    return

# ---------------------------------------------------------------------------------------


def mctdh(op_path, hessian_path, all_frequencies_cm, A, N, **kwargs):
    """ This function creates an `*.op` which will be used by MCTDH.

    It extracts necessary information from the various GAMESS output files.
    It does some simple math (calculating products/sums).
    It then writes all the vibronic model parameters in a formatted output to `*.op`
    """

    # remove the previous file, as we are about to write to it
    _delete_file_using_rmrf(path=op_path)
    # -------------------------------------------------------------------------
    # extract only the frequencies of the modes we selected
    s_idx = [i-1 for i in pp.selected_mode_list]  # need to shift by 1 for 0-indexed arrays
    freq_array = all_frequencies_cm[s_idx]

    # -------------------------------------------------------------------------
    # extract constants
    constants = pp.QM_const
    ha2ev, wn2ev = constants.ha2ev, constants.wn2ev
    qsize = pp.gamess_const.qsize

    # -------------------------------------------------------------------------
    # prepare various lists and dictionaries
    file_name = pp.file_name  # this is the root file/calculation name for all files (displacements etc)

    # path to the output of the GAMESS reference geometry calculation
    ref_geom_path = f'{file_name}_refG.out'

    hessout = hessian_path
    # -------------------------------------------------------------------------
    # bad practice, works for now and we can refactor once we've finished figuring out the end product
    MCTDH_parameter_section_format_string = "{label:<25s}={value:>-15.9f}{units:>8s}\n"
    make_line = functools.partial(MCTDH_parameter_section_format_string.format, units=", ev")
    make_line_cm = functools.partial(MCTDH_parameter_section_format_string.format, units=", cm-1")
    make_line_au = functools.partial(MCTDH_parameter_section_format_string.format, units=", au")
    # -------------------------------------------------------------------------

    def make_op_section(job_title):
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

    # -------------------------------------------------------------------------
    # build functions (where we write the values of the terms)

    def build_frequencies(frequencies):
        """Return a string containing the frequency information of a .op file."""
        return ''.join([
            make_line(label=f"w{i+1:0>2d}", value=w)
            for i, w in enumerate(frequencies)
        ])

    def _incasedebug_E0(E0_array):
        """Return a string containing the energy information of a .op file.
        Has extra code for debugging
        """
        assert False, "only use for debug"

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
                for row, col in upper_triangle_loop_indices(A, 2)
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
            breakpoint()

        return two_blocks()

    def build_E0(E0_array):
        """Return a string containing the energy information of a .op file."""

        diag_block = ''.join([
            make_line(label=f"EH_s{a+1:0>2d}_s{a+1:0>2d}", value=E0_array[a, a])
            for a in range(A)
            if not suppress_zeros or not np.isclose(E0_array[a, a], 0.0)  # if you don't want to print the zeros;
        ])

        # iterates over columns first (upper triangle)
        off_diag_block = ''.join([
            make_line(label=f"EH_s{row+1:0>2d}_s{col+1:0>2d}", value=E0_array[row, col])
            for row, col in upper_triangle_loop_indices(A, 2)
            if not suppress_zeros or not np.isclose(E0_array[row, col], 0.0)  # if you don't want to print the zeros;
        ])

        string = diag_block + "\n" + off_diag_block

        if True:  # add fictious surface
            string += '\n'
            string += make_line(label=f"EH_s{A+1:0>2d}_s{A+1:0>2d}", value=0.0)

        return string

    def build_electronic_moments(dipoles_dict):
        """ Returns a string containing electronic transition dipole moments
            that will be written to a .op file. Takes in dipoles from extract_etdm

            We label both states using the key in the `dipoles_dict` FOR NOW (may change later)
            By default all values in *.op files are in atomic units (au) but we want to explicitly label them
            for clarity.
        """
        def make_xyz_blocks():
            block = ""
            """ Write the xyz component of each excitation as a single block """
            for key in dipoles_dict.keys():
                src, dst = key[0], key[1]  # the #'s identifying the states between which the excitation is occuring
                block += "".join([
                    make_line_au(label=f"E{op}_s{src:>02d}_s{dst:>02d}", value=dipoles_dict[key][xyz_idx])
                    for xyz_idx, op in enumerate(['x', 'y', 'z'])
                    # always print all transition dipole moments out even if zero
                ]) + "\n"

            return block

        def make_state_blocks():
            """ Write all Ex transitions as one block, then repeat with Ey, and then Ez"""
            block = ""
            for xyz_idx, op in enumerate(['x', 'y', 'z']):
                for key in dipoles_dict.keys():
                    src, dst = key[0], key[1]  # the #'s identifying the states between which the excitation is occuring
                    block += "".join([
                        make_line_au(label=f"E{op}_s{src:>02d}_s{dst:>02d}", value=dipoles_dict[key][xyz_idx])
                        # always print all transition dipole moments out even if zero
                    ])
                block += "\n"

            return block

        # just-in-case
        for key in dipoles_dict.keys():  # all keys are length 2 tuples
            assert isinstance(key, tuple) and len(key) == 2

        if True:  # xyz_blocks
            block = make_xyz_blocks()
        else:
            block = make_state_blocks()

        return block

    def build_magnetic_moments():
        raise Exception("Function not implemented yet")

    def build_linear_coupling(lin_dict, A, N):
        """Return a string containing the linear coupling constant information of a .op file.
        `lin_dict` is a dictionary
        for a `selected_mode_list = [7, 8, 9]`
            mode_map_dict[0] = 7
        and so
            lin_dict[mode_map_dict[0]]
        will return the array associated with the key `7`
        Thus
            lin_dict[i][a, a]
        is `lin_dict[i]` returning an array `x`
        and that array is indexed like so: `x[a, a]`
        """

        # make ordered-list of arrays stored in `lin_dict`
        assert len(lin_dict.keys()) == N
        linear = [lin_dict[mode_map_dict[i]] for i in range(N)]

        return '\n'.join([
            ''.join([
                make_line(
                    label=f"C1_s{a+1:0>2d}s{a+1:0>2d}_v{i+1:0>2d}",
                    value=linear[i][a, a]
                )
                for a, i in it.product(range(A), range(N))
                if not suppress_zeros or not np.isclose(linear[i][a, a], 0.0)
            ]),
            ''.join([
                make_line(
                    label=f"C1_s{a1+1:0>2d}s{a2+1:0>2d}_v{i+1:0>2d}",
                    value=linear[i][a1, a2]
                )
                for a1, a2, i in it.product(range(A), range(A), range(N))
                if (a1 < a2)
                and (not suppress_zeros or not np.isclose(linear[i][a1, a2], 0.0))
            ]),
        ])

    def build_quadratic_coupling(quad_dict, A, N):
        """Return a string containing the quadratic coupling constant information of a .op file."""

        # make ordered-list of arrays stored in `lin_data`
        assert len(quad_dict.keys()) == N
        quad = [quad_dict[mode_map_dict[i]] for i in range(N)]

        return '\n'.join([
            ''.join([
                make_line(
                    label=f"C2_s{a+1:0>2d}s{a+1:0>2d}_v{i+1:0>2d}v{i+1:0>2d}",
                    value=quad[i][a, a]
                )
                for a, i in it.product(range(A), range(N))
                if not suppress_zeros or not np.isclose(quad[i][a, a], 0.0)
            ]),
            ''.join([
                make_line(
                    label=f"C2_s{a1+1:0>2d}s{a2+1:0>2d}_v{i+1:0>2d}v{i+1:0>2d}",
                    value=quad[i][a1, a2]
                )
                for a1, a2, i in it.product(range(A), range(A), range(N))
                if (a1 < a2)
                and (not suppress_zeros or not np.isclose(quad[i][a1, a2], 0.0))
            ]),
        ])

    def build_bilinear_coupling(bi_lin_dict, A, N):
        """Return a string containing the BI-Linear coupling constant information of a .op file.

        We first make a new dictionary whose keys are  0-indexed (0, 1) based on `selected_mode_list`.

        So if `bi_lin_dict` has keys
            (7, 8)
        then `bi_lin` has the same values for a corresponding key
            (0, 1)
        """

        bi_lin = {}
        for old_key in bi_lin_dict.keys():
            new_key = reverse_ij_map[old_key]
            bi_lin[new_key] = bi_lin_dict[old_key]

        return '\n'.join([
            ''.join([
                make_line(
                    label=f"C1_s{a+1:0>2d}s{a+1:0>2d}_v{j1+1:0>2d}v{j2+1:0>2d}",
                    # value=0.0
                    value=bi_lin[(j1, j2)][a, a]
                )
                for a, j1, j2 in it.product(range(A), range(N), range(N))
                if (j1 < j2)
                and (not suppress_zeros or not np.isclose(bi_lin[(j1, j2)][a, a], 0.0))
            ]),
            ''.join([
                make_line(
                    label=f"C1_s{a1+1:0>2d}s{a2+1:0>2d}_v{j1+1:0>2d}v{j2+1:0>2d}",
                    # value=0.0
                    value=bi_lin[(j1, j2)][a1, a2]
                )
                for a1, a2, j1, j2 in it.product(range(A), range(A), range(N), range(N))
                if (a1 < a2) and (j1 < j2)
                and (not suppress_zeros or not np.isclose(bi_lin[(j1, j2)][a1, a2], 0.0))
            ]),
        ])

    def build_spinorbit_coupling(soc_dict, A, N):
        """ x """

        def build_linear_SOC(lin_dict, A, N):
            assert len(lin_dict.keys()) == N
            linear_soc = [lin_dict[mode_map_dict[i]] for i in range(N)]

            return ''.join([
                make_line_cm(
                    label=f"C1_s{a1+1:0>2d}s{a2+1:0>2d}_v{i+1:0>2d}r",
                    value=linear_soc[i][a1, a2].real
                ) + make_line_cm(
                    label=f"C1_s{a1+1:0>2d}s{a2+1:0>2d}_v{i+1:0>2d}i",
                    value=linear_soc[i][a1, a2].imag
                )
                for a1, a2, i in it.product(range(A), range(A), range(N))
                if (a1 < a2)
                and (not suppress_zeros or not np.isclose(linear_soc[i][a1, a2], 0.0))
            ]) + '\n'

        def build_quadratic_SOC(quad_dict, A, N):
            # make ordered-list of arrays stored in `lin_data`
            assert len(quad_dict.keys()) == N
            quad = [quad_dict[mode_map_dict[i]] for i in range(N)]

            return ''.join([
                make_line_cm(
                    label=f"C2_s{a1+1:0>2d}s{a2+1:0>2d}_v{i+1:0>2d}v{i+1:0>2d}r",
                    value=quad[i][a1, a2].real
                ) + make_line_cm(
                    label=f"C2_s{a1+1:0>2d}s{a2+1:0>2d}_v{i+1:0>2d}v{i+1:0>2d}i",
                    value=quad[i][a1, a2].imag
                )
                for a1, a2, i in it.product(range(A), range(A), range(N))
                if (a1 < a2)
                and (not suppress_zeros or not np.isclose(quad[i][a1, a2], 0.0))
            ]) + '\n'

        def build_BiLinear_SOC(bi_lin_dict, A, N):

            bi_lin = {}
            for old_key in bi_lin_dict.keys():
                new_key = reverse_ij_map[old_key] # (0, 1) <- (7, 8) for NH3
                bi_lin[new_key] = bi_lin_dict[old_key]

            return ''.join([
                make_line_cm(
                    label=f"C1_s{a1+1:0>2d}s{a2+1:0>2d}_v{j1+1:0>2d}v{j2+1:0>2d}r",
                    value=bi_lin[(j1, j2)][a1, a2].real
                ) + make_line_cm(
                    label=f"C1_s{a1+1:0>2d}s{a2+1:0>2d}_v{j1+1:0>2d}v{j2+1:0>2d}i",
                    value=bi_lin[(j1, j2)][a1, a2].imag
                )
                for a1, a2, j1, j2 in it.product(range(A), range(A), range(N), range(N))
                if (a1 < a2) and (j1 < j2)
                and (not suppress_zeros or not np.isclose(bi_lin[(j1, j2)][a1, a2], 0.0))
            ]) + '\n'

        def build_Total_SOC(total_dict, A, N):

            total = {}
            for old_key in total_dict.keys():
                new_key = reverse_ij_map[old_key]  # (0, 1) <- (7, 8) for NH3
                total[new_key] = total_dict[old_key]

            return ''.join([
                make_line_cm(
                    label=f"SOC_s{a1+1:0>2d}s{a2+1:0>2d}_v{j1+1:0>2d}v{j2+1:0>2d}r",
                    value=total[(j1, j2)][a1, a2].real
                ) + make_line_cm(
                    label=f"SOC_s{a1+1:0>2d}s{a2+1:0>2d}_v{j1+1:0>2d}v{j2+1:0>2d}i",
                    value=total[(j1, j2)][a1, a2].imag
                )
                for a1, a2, j1, j2 in it.product(range(A), range(A), range(N), range(N))
                if (a1 < a2) and (j1 < j2)
                and (not suppress_zeros or not np.isclose(total[(j1, j2)][a1, a2], 0.0))
            ]) + '\n'

        return '\n'.join([
            build_linear_SOC(soc_dict['Linear'], A, N),
            build_quadratic_SOC(soc_dict['Quadratic'], A, N),
            build_BiLinear_SOC(soc_dict['BiLinear'], A, N),
            build_Total_SOC(soc_dict['Total'], A, N),
        ]) + '\n'

    def make_parameter_section(model, A, N):
        """Returns a string which defines the `PARAMETER-SECTION` of an .op file"""
        start, end = "PARAMETER-SECTION", "end-parameter-section"
        # ----------------------------------------------------------
        # this is what controls the look of the header strings
        header_format_string = "#{0:^47}#\n#{spacer:^47}#\n"
        make_header = functools.partial(header_format_string.format, spacer='-' * 45)
        # ----------------------------------------------------------
        # read in all the necessary parameters
        dipoles = model['dipoles']
        vibron_ev = model['vibron eV']
        E0_array_eV = model['E0 eV']

        # ----------------------------------------------------------
        # begin to assemble the parameter section
        return_list = [
            start,
            make_header('Frequencies'), build_frequencies(vibron_ev),
            make_header('Electronic Hamitonian'), build_E0(E0_array_eV),
            make_header('Electronic transition moments'), build_electronic_moments(dipoles),
            # make_header('Magnetic transition moments'), build_magnetic_moments(M_moments),
        ]

        # -----------------------------------------------------------
        # required elements of model (although for debugging we may want to disable for fast debug)
        for key in ['Linear', 'Quadratic', 'BiLinear']:
            if key not in model.keys():
                print(f'NO {key} coupling found in {model.keys()}');
                breakpoint(); raise Exception()  # (temporarily) comment out if you don't need them

        # soft fail, may not always want spin-orbit-coupling?
        if 'SOC' not in model.keys():
            print(f'NO spin orbit coupling found in {model.keys()}')

        # -----------------------------------------------------------

        headers = {
            'Linear': 'Linear Coupling Constants',
            'Quadratic': 'Quadratic Coupling Constants',
            'BiLinear': 'Diagonal Bilinear Coupling Constants',
            'SOC': 'SOC',
        }

        # add them if present
        key = 'Linear'
        if key in model.keys():
            return_list += [make_header(headers[key]),  build_linear_coupling(model[key], A, N)]

        key = 'Quadratic'
        if key in model.keys():
            return_list += [make_header(headers[key]),  build_quadratic_coupling(model[key], A, N)]

        key = 'BiLinear'
        if key in model.keys():
            return_list += [make_header(headers[key]),  build_bilinear_coupling(model[key], A, N)]

        key = 'SOC'
        if key in model.keys():
            return_list += [make_header(headers[key]),  build_spinorbit_coupling(model[key], A, N)]

        return '\n'.join(return_list + [end, ])
    # -------------------------------------------------------------------------
    # labeling functions (where we label the Q.M. analog of the terms)

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

    def label_energies(energy, A):
        """Return a string containing the energy labelling of a .op file."""
        spacer = '|'

        diag_block = '\n'.join([
            f"EH_s{a:0>2d}_s{a:0>2d}{spacer:>15}1 S{a:d}&{a:d}"
            for a in range(1, A+1)
        ]) + '\n'

        off_diag_block = '\n'.join([
            f"EH_s{a1:0>2d}_s{a2:0>2d}{spacer:>15}1 S{a1:d}&{a2:d}"
            for a1, a2 in it.product(range(1, A+1), range(1, A+1))
            if (a1 < a2)
            and (not suppress_zeros or not np.isclose(energy[a1-1, a2-1], 0.0))
        ]) + '\n'

        string = diag_block + "\n" + off_diag_block

        if True:  # add fictious surface
            string += '\n'
            string += f"EH_s{A+1:0>2d}_s{A+1:0>2d}{spacer:>15}1 S{A+1:d}&{A+1:d}"

        return string

    def label_linear_coupling(lin_dict, A, N):
        """Return a string containing the linear coupling constant labelling of a .op file."""
        spacer = '|'

        assert len(lin_dict.keys()) == N
        linear_terms = [lin_dict[mode_map_dict[i]] for i in range(N)]

        return '\n'.join([
            (
                f"C1_s{a:0>2d}s{a:0>2d}_v{i:0>2d}"
                f"{spacer:>11}1 S{a:d}&{a:d}"
                f"{spacer:>4}{i+1}  q"
            )
            for a, i in it.product(range(1, A+1), range(1, N+1))
            if not suppress_zeros or not np.isclose(linear_terms[i-1][a-1, a-1], 0.0)
        ] + [
            ''  # creates a blank line between the (surface) diagonal and off-diagonal linear terms
        ] + [
            (
                f"C1_s{a1:0>2d}s{a2:0>2d}_v{i:0>2d}"
                f"{spacer:>11}1 S{a1:d}&{a2:d}"
                f"{spacer:>4}{i+1}  q"
            )
            for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
            if (a1 < a2)
            and (not suppress_zeros or not np.isclose(linear_terms[i-1][a1-1, a2-1], 0.0))
        ]) + '\n'

    def label_quadratic_coupling(quad_dict, A, N):
        """Return a string containing the quadratic coupling constant labelling of a .op file."""
        spacer = '|'

        assert len(quad_dict.keys()) == N
        quadratic_terms = [quad_dict[mode_map_dict[i]] for i in range(N)]

        return '\n'.join([
            (
                f"0.50*C2_s{a:0>2d}s{a:0>2d}_v{i:0>2d}v{i:0>2d}"
                f"{spacer:>9}1 S{a:d}&{a:d}"
                f"{spacer:>4}{i+1}  q^2"
            )
            for a, i in it.product(range(1, A+1), range(1, N+1))
            if not suppress_zeros or not np.isclose(quadratic_terms[i-1][a-1, a-1], 0.0)
        ] + [
                ''  # creates a blank line between the (surface) diagonal and off-diagonal linear terms
        ] + [
            (
                f"0.50*C2_s{a1:0>2d}s{a2:0>2d}_v{i:0>2d}v{i:0>2d}"
                f"{spacer:>9}1 S{a1:d}&{a2:d}"
                f"{spacer:>4}{i+1}  q^2"
            )
            for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
            if (a1 < a2)
            if not suppress_zeros or not np.isclose(quadratic_terms[i-1][a1-1, a2-1], 0.0)
        ]) + '\n'

    def label_bilinear_coupling(bi_lin_dict, A, N):
        """Return a string containing the quadratic coupling constant labelling of a .op file."""
        spacer = '|'

        bilinear_terms = {}
        for old_key in bi_lin_dict.keys():
            new_key = reverse_ij_map[old_key]
            bilinear_terms[new_key] = bi_lin_dict[old_key]

        return '\n'.join([
            (
                f"C1_s{a:0>2d}s{a:0>2d}_v{j1:0>2d}v{j2:0>2d}"
                f"{spacer:>9}1 S{a:d}&{a:d}"
                f"{spacer:>4}{j1+1}  q{spacer:>6}{j2+1}  q"
            )
            for a, j1, j2 in it.product(range(1, A+1), range(1, N+1), range(1, N+1))
            if (j1 < j2)
            and (not suppress_zeros or not np.isclose(bilinear_terms[(j1-1, j2-1)][a-1, a-1], 0.0))
        ] + [
                ''  # creates a blank line between the (surface) diagonal and off-diagonal linear terms
        ] + [
            (
                f"C1_s{a1:0>2d}s{a2:0>2d}_v{j1:0>2d}v{j2:0>2d}"
                f"{spacer:>9}1 S{a1:d}&{a2:d}"
                f"{spacer:>4}{j1+1}  q{spacer:>6}{j2+1}  q"
            )
            for a1, a2, j1, j2 in it.product(range(1, A+1), range(1, A+1), range(1, N+1), range(1, N+1))
            if (a1 < a2) and (j1 < j2)
            and (not suppress_zeros or not np.isclose(bilinear_terms[(j1-1, j2-1)][a1-1, a2-1], 0.0))
        ]) + '\n'

    def label_spinorbit_coupling(soc_dict, A, N):

        spacer_format_string = f"# {'-':^60s} #\n"
        hfs = header_format_string = "# {:^60s} #\n" + spacer_format_string
        block = hfs.format("SOC FULL HAMILTONIAN SOC OFF-DIAGONAL VIBRONIC COUPLINGS")

        """ Full lines look like this:
         I*C1_s##_s##_v##r |1 Z#&# | # q
        -I*C1_s##_s##_v##i |1 Z#&# | # q
         I*C2_s##_s##_v##r |1 Z#&# | # q^2
        -I*C2_s##_s##_v##i |1 Z#&# | # q^2
         I*C1_s##_s##_v##_v##r |1 Z#&# | # q | # q
        -I*C1_s##_s##_v##_v##i |1 Z#&# | # q | # q
         I*SOC_s##_s##_v##_v##r |1 Z#&# | # q
        -I*SOC_s##_s##_v##_v##i |1 Z#&# | # q
        """

        def label_linear_SOC(lin_dict, A, N):
            """Return a string containing the linear spin orbit coupling (SOC) terms"""
            spacer = '|'
            assert len(lin_dict.keys()) == N
            linear_terms = [lin_dict[mode_map_dict[i]] for i in range(N)]

            return '\n'.join([
                (
                    f"C1_s{a1:0>2d}s{a2:0>2d}_v{i:0>2d}r"
                    f"{spacer:>11}1 S{a1:d}&{a2:d}"
                    f"{spacer:>4}{i+1}  I*q"
                ) + '\n' + (
                    f"C1_s{a1:0>2d}s{a2:0>2d}_v{i:0>2d}i"
                    f"{spacer:>11}1 S{a1:d}&{a2:d}"
                    f"{spacer:>4}{i+1}  -I*q"
                )
                for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
                if (a1 < a2)
                and (not suppress_zeros or not np.isclose(linear_terms[i-1][a1-1, a2-1], 0.0))
            ]) + '\n'

        def label_quadratic_SOC(quad_dict, A, N):
            """Return a string containing the quadratic spin orbit coupling (SOC) terms"""
            spacer = '|'
            assert len(quad_dict.keys()) == N
            quadratic_terms = [quad_dict[mode_map_dict[i]] for i in range(N)]

            return '\n'.join([
                (
                    f"C2_s{a1:0>2d}s{a2:0>2d}_v{i:0>2d}v{i:0>2d}r"
                    f"{spacer:>9}1 S{a1:d}&{a2:d}"
                    f"{spacer:>4}{i+1}  I*q^2"
                ) + '\n' + (
                    f"C2_s{a1:0>2d}s{a2:0>2d}_v{i:0>2d}v{i:0>2d}i"
                    f"{spacer:>9}1 S{a1:d}&{a2:d}"
                    f"{spacer:>4}{i+1}  -I*q^2"
                )
                for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
                if (a1 < a2)
                if not suppress_zeros or not np.isclose(quadratic_terms[i-1][a1-1, a2-1], 0.0)
            ]) + '\n'

        def label_BiLinear_SOC(bi_lin_dict, A, N):
            """Return a string containing the BiLinear spin orbit coupling (SOC) terms"""
            spacer = '|'
            bilinear_terms = {}
            for old_key in bi_lin_dict.keys():  # i should make a backwards mapping dictionary
                new_key = reverse_ij_map[old_key]  # (0, 1) <- (7, 8) for NH3
                bilinear_terms[new_key] = bi_lin_dict[old_key]

            return '\n'.join([
                (
                    f" I*C1_s{a1:0>2d}s{a2:0>2d}_v{j1:0>2d}v{j2:0>2d}r"
                    f"{spacer:>9}1 S{a1:d}&{a2:d}"
                    f"{spacer:>4}{j1+1}  q{spacer:>6}{j2+1}  q"
                ) + '\n' + (
                    f"-I*C1_s{a1:0>2d}s{a2:0>2d}_v{j1:0>2d}v{j2:0>2d}i"
                    f"{spacer:>9}1 S{a1:d}&{a2:d}"
                    f"{spacer:>4}{j1+1}  q{spacer:>6}{j2+1}  q"
                )
                for a1, a2, j1, j2 in it.product(range(1, A+1), range(1, A+1), range(1, N+1), range(1, N+1))
                if (a1 < a2) and (j1 < j2)
                and (not suppress_zeros or not np.isclose(bilinear_terms[(j1-1, j2-1)][a1-1, a2-1], 0.0))
            ]) + '\n'

        def label_Total_SOC(total_dict, A, N):
            """Return a string containing the BiLinear spin orbit coupling (SOC) terms"""
            spacer = '|'
            total_terms = {}
            for old_key in total_dict.keys():  # i should make a backwards mapping dictionary
                new_key = reverse_ij_map[old_key] # (0, 1) <- (7, 8) for NH3
                total_terms[new_key] = total_dict[old_key]

            return '\n'.join([
                (
                    f" I*SOC_s{a1:0>2d}s{a2:0>2d}_v{j1:0>2d}v{j2:0>2d}r"
                    f"{spacer:>9}1 S{a1:d}&{a2:d}"
                    f"{spacer:>4}{j1+1}  q{spacer:>6}{j2+1}  q"
                ) + '\n' + (
                    f"-I*SOC_s{a1:0>2d}s{a2:0>2d}_v{j1:0>2d}v{j2:0>2d}i"
                    f"{spacer:>9}1 S{a1:d}&{a2:d}"
                    f"{spacer:>4}{j1+1}  q{spacer:>6}{j2+1}  q"
                )
                for a1, a2, j1, j2 in it.product(range(1, A+1), range(1, A+1), range(1, N+1), range(1, N+1))
                if (a1 < a2) and (j1 < j2)
                and (not suppress_zeros or not np.isclose(total_terms[(j1-1, j2-1)][a1-1, a2-1], 0.0))
            ]) + '\n'

        def _toby_bash_style():
            raise Exception('not implemented')

            # prepare `make_line`
            format_string_1 = "{label:<25s}{link:<20s}\n"
            format_string_2 = "{label:<25s}{link:<20s}\n"
            format_string_3 = "{label:<25s}{link:<20s}\n"
            format_string_4 = "{label:<25s}{link:<20s}\n"

            # make_line_1_real = functools.partial(format_string_1.format, link=)
            # make_line_1_imag = functools.partial(format_string_1.format, link=)
            # make_line_2 = functools.partial(format_string_2.format, link=)
            # make_line_3 = functools.partial(format_string_3.format, link=)
            # make_line_4 = functools.partial(format_string_4.format, link=)

            for i, a in it.product(range(1, N+1), range(1, A+1)):
                for j, b in it.product(range(1, i+1), range(1, a+1)):
                    i_label, j_label = reverse_ij_map[i,j]
                    print(f"{i=}, {j=}, {a=}, {b=}")

                    l1 = "C1_s{:>02d}_s{:>02d}_v{:>02d}".format(j, i, a)
                    make_line_1(label=f" I*{l1}r", link=f"|1 S{b}&{a} | {j+1} q")
                    make_line_1(label=f"-I*{l1}i", link=f"|1 S{a}&{b} | {j+1} q")

                    l2 = "C2_s{:>02d}_s{:>02d}_v{:>02d}".format(j, i, a)
                    make_line_1(label=f" I*{l2}r", link=f"|1 S{b}&{a} | {j+1} q^2")
                    make_line_1(label=f"-I*{l2}i", link=f"|1 S{a}&{b} | {j+1} q^2")

                    l3 = "C1_s{:>02d}_s{:>02d}_v{:>02d}_v{:>02d}".format(j, i, i_label, j_label)
                    make_line(label=f" I*{l3}r", link=f"|1 S{b}&{a} | {j+1} q | {i+1} q")
                    make_line(label=f"-I*{l3}i", link=f"|1 S{a}&{b} | {j+1} q | {i+1} q")

                    l4 = "SOC_s{:>02d}_s{:>02d}_v{:>02d}_v{:>02d}".format(j, i, i_label, j_label)
                    make_line(label=f" I*{l4}r", link=f"|1 S{b}&{a} | {j+1} q")
                    make_line(label=f"-I*{l4}i", link=f"|1 S{a}&{b} | {j+1} q")

        if True:
            string = "\n".join([
                label_linear_SOC(soc_dict['Linear'], A, N),
                label_quadratic_SOC(soc_dict['Quadratic'], A, N),
                label_BiLinear_SOC(soc_dict['BiLinear'], A, N),
                label_Total_SOC(soc_dict['Total'], A, N),
            ]) + '\n'

        else:  # not implemented
            string = _toby_bash_style()

        return string

    def make_hamiltonian_section(model, A, N, SOC_flag=False):
        """Returns a string which defines the `HAMILTONIAN-SECTION` of an .op file"""
        start, end = "HAMILTONIAN-SECTION", "end-hamiltonian-section"
        spec = ''.join([
            ' modes   |  el  |',
            ''.join([f" v{N+1:0>2d}|" for N in range(N)]),
            '\n'
        ])

        # -----------------------------------------------------------
        return_list = [
            start,
            spec,
            label_momentum(N),
            label_position(N),
            label_energies(model['E0 eV'], A),
        ]

        # -----------------------------------------------------------
        # required elements of model (although for debugging we may want to disable for fast debug)
        for key in ['Linear', 'Quadratic', 'BiLinear']:
            if key not in model.keys():
                print(f'NO {key} coupling found in {model.keys()}');
                breakpoint(); raise Exception()  # (temporarily) comment out if you don't need them

        # soft fail, may not always want spin-orbit-coupling?
        if 'SOC' not in model.keys():
            print(f'NO spin orbit coupling found in {model.keys()}')

        # -----------------------------------------------------------
        # add them if present
        key = 'Linear'
        if key in model.keys():
            return_list += [label_linear_coupling(model[key], A, N),]

        key = 'Quadratic'
        if key in model.keys():
            return_list += [label_quadratic_coupling(model[key], A, N)]

        key = 'BiLinear'
        if key in model.keys():
            return_list += [label_bilinear_coupling(model[key], A, N)]

        key = 'SOC'
        if key in model.keys():
            return_list += [label_spinorbit_coupling(model[key], A, N)]

        # -----------------------------------------------------------
        return_list.append(end)

        return '\n'.join(return_list)
    # -------------------------------------------------------------------------

    def neil_build_dipole_moments_section(nof_states, nof_modes):
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

    def make_operator_onto_dipole_moments_section(model, A, N):
        """  """
        block = f"\nHAMILTONIAN-SECTION_Ex\n"

        # Write modes and mode labels
        mode_number_key = [selected_mode_list[i] for i in range(N)]
        h_labels = ''.join([
            ' modes   |  el  |',
            ''.join([f" v{N+1:0>2d}|" for N in range(N)]),
            '\n'
        ])

        block += h_labels + "\n"
        block += f"{'-'*47}\n\n"

        # h_labels = ["modes", "el", ] + [
        #     f"v{s:>02d}"
        #     for s in mode_number_key
        # ]

        # block += " | ".join(h_labels) + "\n"

        for j in range(1, A+1):
            block += f"1.0         |1 S{A+1}&{j}\n"  # A+1, set ground state as fictitious +1 state

        block += "\nend-hamiltonian-section\n"

        return block

    # -------------------------------------------------------------------------
    # extraction functions

    def extract_E0(path, A=pp.A):
        """ The energy values associated with the reference geometry.
        These energies don't depend on the modes.
        Reads in the energies from the hessian file at `path`
        """

        # strings used by `grep` to locate values to extract
        a_pattern = 'STATE #.* {col}.S GMC-PT-LEVEL DIABATIC ENERGY='
        ba_pattern = 'STATE #.* {row} &.* {col}.S GMC-PT-LEVEL COUPLING'

        # ground state of the (optimized geometry?) (includes fictitious surface? no/yes?)
        GSE_pattern = 'TOTAL ENERGY ='
        try:
            GSE = extract_ground_state_energy(path, GSE_pattern)
        except Exception as e:
            print("Cannot find ground state energy! What to do?")
            print(str(e))
            print(f"{path=} {GSE_pattern=}")
            breakpoint()

        # 1st diabat's energy
        D1E = extract_diabatic_energy(ref_geom_path, a_pattern.format(col=1))
        linear_shift = (D1E - GSE) * ha2ev
        print(
            f'The ground state energy is: {GSE} Hartree\n'
            f'Diabat #1 energy is: {D1E} Hartree\n'
            f'Linear shift value: {linear_shift} eV\n'
        )


        A = pp.A
        shape = (A, A)
        E0_array = np.zeros(shape)

        if False and __debug__:  # debug
            # check this function to remind yourself how the indexing works
            _reminder_produce_upper_triangle_indices(A)

        def _toby_bash_style(E0_array):
            """ This matches the bash script looping style used by Toby.
            The outer loop is over the columns of the matrix a ~ col.
            The inner loop is over the rows of the matrix b ~ row.
            """
            for a in range(A):
                E0_array[a, a] = refG_extract(ref_geom_path, a_pattern.format(col=a+1))
                E0_array[a, a] += linear_shift
                for b in range(a):
                    E0_array[b, a] = refG_extract(ref_geom_path, ba_pattern.format(row=b+1, col=a+1))
            return

        def _row_first_style(E0_array):
            """ x """
            for a in range(A):
                E0_array[a, a] = extract_in_eV(ref_geom_path, a_pattern.format(col=a+1))
                E0_array[a, a] += linear_shift

            for a, b in upper_triangle_loop_indices(A, 2):
                E0_array[a, b] = extract_in_eV(ref_geom_path, ba_pattern.format(row=a+1, col=b+1))

            # for row, col in upper_triangle_loop_indices(A, 2):
            #     E0_array[row, col] = extract_in_eV(ref_geom_path, ba_pattern.format(row=row+1, col=col+1))

            return

        if True:  # this makes more sense based on what I see in the file
            _row_first_style(E0_array)
            # print("Row first", E0_array)

        if False:  # matches the bash style from Toby
            _toby_bash_style(E0_array)
            # print("Column first / bash style", E0_array)

        if False and __debug__:  # (keep them out of the loops above for simplicity)
            for a in range(A):
                print(f"Diabatic energy at state {a+1}: {E0_array[a, a]}")
            for a, b in upper_triangle_loop_indices(A, 2):
                print(f"Coupling energy at state {a+1} & {b+1}: {E0_array[a, b]}")
            # for row, col in upper_triangle_loop_indices(A, 2):
            #     print(f"Coupling energy at state {row+1} & {col+1}: {E0_array[row, col]}")
            print(E0_array)

        # added at the last minute (change later)
        E0_array_ev = E0_array
        E0_array_au = np.zeros(shape)

        # now also get the array in atomic units (Hartrees)
        for a in range(A):
            E0_array_au[a, a] = extract_in_Hartrees(ref_geom_path, a_pattern.format(col=a+1))

        for a, b in upper_triangle_loop_indices(A, 2):
            E0_array_au[a, b] = extract_in_Hartrees(ref_geom_path, ba_pattern.format(row=a+1, col=b+1))

        return E0_array_ev, E0_array_au

    def extract_etdm(path, verbose=False):
        """ Extracts the electronic transition dipole moments from refG.out
            It will extract tdm from diabat 1 -> diabat 2,3,...
            Returns dipoles, a dictionary of tdm values, columns are x,y,z

            e.g. {2: ('-0.000000', '0.182262', '0.000000'),
                  3: ('-0.000000', '0.000000', '0.000000')}
        """
        def extract_ground_to_excited_state_transition_dipoles(selected_lines):
            """ x """
            dipoles = {}

            for line in selected_lines:
                try:
                    state, pair, x, y, z = line.strip().split()
                    state, pair = int(state), int(pair)
                    if pair == state:  # transition from the fictitious ground state (0) -> to `state` #
                        pair = 0
                    dipoles[(pair, state)] = [float(x), float(y), float(z)]

                except Exception as e:
                    print(f"Error processing line: {line} - {e}")
                    breakpoint()

            return dipoles

        start = "TRANSITION DIPOLES BETWEEN DIABATS"
        end = "TRANSITION DIPOLES BETWEEN DIRECT MAX. DIABATS"
        tdipole_block = extract_lines_between_patterns(path, start, end)

        selected_lines = [  # remove all '\n' and blank lines
            line.strip() for line in tdipole_block
            if "#" not in line and line.strip() != ""
        ]

        selected_lines = selected_lines[2:]  # skip the 2 headers as shown below
        """     TRANSITION DIPOLES BETWEEN DIABATS (in A.U.):
                STATE PAIR      X         Y         Z               """

        dipoles = extract_ground_to_excited_state_transition_dipoles(selected_lines)

        if verbose and __debug__:
            pprint.pprint(tdipole_block)
            pprint.pprint(dipoles)

        return dipoles

    def extract_linear(array_style=True):
        """
        For some reason the diagonal (of linear and quadratic) needs to be in Hartrees???
        and the off-diagaonal in eV??
        """

        # strings used by `grep` to locate values to extract
        a_pattern = 'STATE #.* {col}.S GMC-PT-LEVEL DIABATIC ENERGY='
        ba_pattern = 'STATE #.* {row} &.* {col}.S GMC-PT-LEVEL COUPLING'

        shape = (A, A)

        def extract_energy_at_displaced_geometry(path, key):
            """ Diagonal values are in Hartrees (atomic units - au)
            Off-diagonal values are in electron Volts (eV)
            """
            array = np.zeros(shape)

            for a in range(A):
                array[a, a] = extract_in_Hartrees(path, a_pattern.format(col=a+1))

            for a, b in upper_triangle_loop_indices(A, 2):
                array[a, b] = extract_in_eV(path, ba_pattern.format(row=a+1, col=b+1))

            if False and __debug__:  # debug printing
                for a in range(A):
                    print(f"Linear energy {key=} at state {a+1}: {array[a, a]}")
                for a, b in upper_triangle_loop_indices(A, 2):
                    print(f"Linear energy {key=} at state {a+1} & {b+1}: {array[a, b]}")
                print(path); print(array); breakpoint()
            return array

        def _compute_using_array_style(i, temp_dict):
            linear_ev = (temp_dict["+1"] - temp_dict["-1"]) / (2*qsize)
            for a in range(A):  # make sure we multiply the diagonal by ha2ev
                linear_ev[a, a] *= ha2ev
            return linear_ev

        def _compute_using_forloop_style(i, temp_dict):
            linear_ev = np.zeros(shape)
            for a in range(A):
                linear_ev[a, a] = (temp_dict["+1"][a, a] - temp_dict["-1"][a, a]) * ha2ev / (2 * qsize)
                for b in range(a):
                    linear_ev[b, a] = (temp_dict["+1"][b, a] - temp_dict["-1"][b, a]) / (2 * qsize)
            return linear_ev

        linear_dictionary = {}  # store return values in here
        for i in range(N):
            # ----------------------------------------------------------
            # extract values from GAMESS files
            temp_dict = {}  # use new dictionary to store temporary arrays every single loop (i \in N)
            for key in ["+1", "-1"]:
                assert key in linear_disp_keys, f"{key=} not in {linear_disp_keys=}"
                path = linear_displacement_filenames[(key, i)]
                temp_dict[key] = extract_energy_at_displaced_geometry(path, key)
                print(i, path)

            # ----------------------------------------------------------
            # preform the math / calculation using the values we just extracted
            if array_style:
                linear_ev = _compute_using_array_style(i, temp_dict)
                if False and __debug__: print(linear_ev); breakpoint()  # debug

            else:  # do it with loops
                linear_ev = _compute_using_forloop_style(i, temp_dict)
                if False and __debug__: print(linear_ev); breakpoint()  # debug

            # ----------------------------------------------------------
            # store the value in the dictionary
            linear_dictionary[mode_map_dict[i]] = linear_ev

        return linear_dictionary

    def extract_quadratic(E0_array_eV, E0_array_au, vibron_ev, array_style=True):
        """ x """

        # strings used by `grep` to locate values to extract
        a_pattern = 'STATE #.* {col}.S GMC-PT-LEVEL DIABATIC ENERGY='
        ba_pattern = 'STATE #.* {row} &.* {col}.S GMC-PT-LEVEL COUPLING'
        shape = (A, A)

        def extract_energy_at_displaced_geometry(path, key):
            """ x """
            array = np.zeros(shape)

            for a in range(A):
                array[a, a] = extract_in_Hartrees(path, a_pattern.format(col=a+1))

            for a, b in upper_triangle_loop_indices(A, 2):
                array[a, b] = extract_in_eV(path, ba_pattern.format(row=a+1, col=b+1))

            if False and __debug__:  # debug printing
                for a in range(A):
                    print(f"Quadratic energy {key=} at state {a+1}: {array[a, a]}")
                for a, b in upper_triangle_loop_indices(A, 2):
                    print(f"Quadratic energy {key=} at state {a+1} & {b+1}: {array[a, b]}")
                print(path); print(array); breakpoint()  # debug

            return array

        def _compute_using_array_style(i, temp_dict):
            # until we change the extraction
            temp_fake_E0_array = np.copy(E0_array_eV)
            for a in range(A):  # make sure the diagonal is in Hartrees
                temp_fake_E0_array[a, a] = E0_array_au[a, a]

            # have to because they have different E0's (their energy scales are all different)
            quad_ev = temp_dict["+2"] + temp_dict["-2"] - 2.0 * temp_fake_E0_array
            for a in range(A):  # make sure we multiply the diagonal by ha2ev
                quad_ev[a, a] *= ha2ev

            quad_ev /= (2. * qsize) ** 2
            for a in range(A):  # make sure we subtract the vibron_ev on the diagonal
                quad_ev[a, a] -= vibron_ev[i]

            return quad_ev

        def _compute_using_forloop_style(i, temp_dict):
            quad_ev = np.zeros(shape)
            for a in range(A):
                quad_ev[a, a] = temp_dict["+2"][a, a] + temp_dict["-2"][a, a]
                quad_ev[a, a] -= 2.0 * E0_array_au[a, a]
                quad_ev[a, a] *= ha2ev  # convert back to eV's  (only for the diagonal)
                quad_ev[a, a] /= (2. * qsize) ** 2  # or equivalently: (4.0 * qsize * qsize)
                quad_ev[a, a] -= vibron_ev[i]
                for b in range(a):
                    quad_ev[b, a] = temp_dict["+2"][b, a] + temp_dict["-2"][b, a]
                    quad_ev[b, a] -= 2.0 * E0_array_eV[b, a]
                    quad_ev[b, a] /= (2. * qsize) ** 2   # or equivalently: (4.0 * qsize * qsize)

            return quad_ev

        quadratic_dictionary = {}  # store return values in here
        for i in range(N):
            # ----------------------------------------------------------
            # extract values from GAMESS files
            temp_dict = {}  # stores temporary arrays (they are different every single i \in N loop)
            for key in ["+2", "-2"]:
                assert key in linear_disp_keys, f"{key=} not in {linear_disp_keys=}"
                path = linear_displacement_filenames[(key, i)]
                temp_dict[key] = extract_energy_at_displaced_geometry(path, key)
                print(i, path)
            if False and __debug__: print(temp_dict["+2"], '\n', temp_dict["-2"]); breakpoint()

            # ----------------------------------------------------------
            # preform the math / calculation using the values we just extracted
            if not array_style:  # do it with for loops
                quad_ev = _compute_using_forloop_style(i, temp_dict)
                if False and __debug__: print(i, quad_ev, vibron_ev); breakpoint()
            else:
                quad_ev = _compute_using_array_style(i, temp_dict)
                if False and __debug__: print(quad_ev); breakpoint()

            # ----------------------------------------------------------
            # store the value in the dictionary
            quadratic_dictionary[pp.mode_map_dict[i]] = quad_ev

        return quadratic_dictionary

    def extract_bilinear(array_style=True):
        """ x """

        # strings used by `grep` to locate values to extract
        a_pattern = 'STATE #.* {col}.S GMC-PT-LEVEL DIABATIC ENERGY='
        ba_pattern = 'STATE #.* {row} &.* {col}.S GMC-PT-LEVEL COUPLING'

        A, N = pp.A, pp.N  # stop flake8 from complaining
        shape = (A, A)

        def extract_energy_at_displaced_geometry(path, key):
            """ x """
            array = np.zeros(shape)

            for a in range(A):
                array[a, a] = extract_in_Hartrees(path, a_pattern.format(col=a+1))

            for a, b in upper_triangle_loop_indices(A, 2):
                array[a, b] = extract_in_eV(path, ba_pattern.format(row=a+1, col=b+1))

            if False and __debug__:  # debug printing
                for a in range(A):
                    print(f"Bilinear energy {key=} at state {a+1}: {array[a, a]}")
                for a, b in upper_triangle_loop_indices(A, 2):
                    print(f"Bilinear energy {key=} at state {a+1} & {b+1}: {array[a, b]}")
                print(path); print(array); breakpoint()  # debug

            return array

        def _compute_using_array_style(temp_dict):
            """
            Remember that the diagonal of `temp_dict` are in hartrees,
            while the off-diagonal are in eV.
            """
            bilin_ev = np.zeros(shape)

            bilin_ev += temp_dict['++'] - temp_dict['+-']
            bilin_ev += temp_dict['--'] - temp_dict['-+']

            for a in range(A):  # make sure we multiply the diagonal by ha2ev
                bilin_ev[a, a] *= ha2ev

            bilin_ev /= (2. * qsize) ** 2
            return bilin_ev

        def _compute_using_forloop_style(i, j, temp_dict):
            """ bilin_ev needs to run through all temp_dict[key], so beyond '++'.
            You only get '--'' after temp_dict fully built
            """
            bilin_ev = np.zeros(shape)
            counter = 0

            for a in range(A):
                bilin_ev[a, a] += temp_dict['++'][a, a] - temp_dict['+-'][a, a]
                bilin_ev[a, a] += temp_dict['--'][a, a] - temp_dict['-+'][a, a]
                bilin_ev[a, a] *= ha2ev   # convert back to eV's (only for the diagonal)
                bilin_ev[a, a] /= (4.0 * qsize * qsize)

                for b in range(a):
                    bilin_ev[b, a] += temp_dict['++'][b, a] - temp_dict['+-'][b, a]
                    bilin_ev[b, a] += temp_dict['--'][b, a] - temp_dict['-+'][b, a]
                    bilin_ev[b, a] /= (4.0 * qsize * qsize)

                    if False and __debug__:
                        counter += 1
                        print(f'BILIN_EV at mode {i} & {j}, processing iteration #{counter}, states {b+1} & {a+1}')
                        pprint.pprint(bilin_ev)

            return bilin_ev

        bilinear_dictionary = {}  # store return values in here
        for i, j in upper_triangle_loop_indices(N, 2):

            # ----------------------------------------------------------
            # extract values from GAMESS files
            temp_dict = {}  # stores temporary arrays (different every single i,j \in N x N loop)
            for key in bi_linear_disp_keys:
                path = bi_linear_displacement_filenames[(key, i, j)]
                temp_dict[key] = extract_energy_at_displaced_geometry(path, key)
                print(i, j, path)

            if False and __debug__:
                print(f'TEMP DICT at mode {i} & {j}'); pprint.pprint(temp_dict)

            # ----------------------------------------------------------
            # preform the math / calculation using the values we just extracted
            if not array_style:  # do it with for loops
                bilin_ev = _compute_using_forloop_style(i, j, temp_dict)
                if False and __debug__: print(i, j, bilin_ev); breakpoint()
            else:
                bilin_ev = _compute_using_array_style(temp_dict)
                if False and __debug__:
                    print(f'BI-Linear at mode {i} & {j}');
                    pprint.pprint(bilin_ev);
                    breakpoint()

            # ----------------------------------------------------------
            # store the value in the dictionary
            key = ij_map[(i, j)]
            bilinear_dictionary[key] = bilin_ev

        if True and __debug__:  # optional printing of BiLinear dictionary
            print(
                '\n',
                f"{'':-^25}{'FULL BILINEAR DICTIONARY!':^30}{'':-^25}",
                '\n'.join([f"{k}\n{v}" for k, v in bilinear_dictionary.items()]),
                f"{'':-^80}",
                sep='\n'
            )

        return bilinear_dictionary

    def extract_soc(array_style=True):
        """ x """
        shape = (A, A)

        def _extract_linear_soc():
            """ x """
            def _compute_using_array_style(temp_dict):
                return (temp_dict["+1"] - temp_dict["-1"]) / (2*qsize)

            def _compute_using_forloop_style(i, temp_dict):
                linear_ev = np.zeros(shape)
                for a, b in upper_triangle_loop_indices(range(A), 2):
                    linear_ev[b, a] = (temp_dict["+1"][b, a] - temp_dict["-1"][b, a]) / (2*qsize)

                return linear_ev

            lin_dict = {}  # store return values in here
            for i in range(N):  # do the 1D spin_orbit_coupling terms

                # ----------------------------------------------------------
                # extract Linear from GAMESS files
                temp_dict = {}  # stores temporary arrays (different every single i,j \in N x N loop)
                for key in ['+1', '-1']:
                    assert key in linear_disp_keys, f"{key=} not in {linear_disp_keys=}"
                    path = linear_displacement_filenames[(key, i)]
                    temp_dict[key] = extract_DSOME(path, A)

                # ----------------------------------------------------------
                # preform the math / calculation using the values we just extracted
                if not array_style:  # do it with for loops
                    soc_ev = _compute_using_forloop_style(i, temp_dict)
                else:
                    soc_ev = _compute_using_array_style(temp_dict)

                # ----------------------------------------------------------
                # store the value in the dictionary
                lin_dict[pp.mode_map_dict[i]] = soc_ev

            print("Finished extracting linear SOC")
            return lin_dict

        def _extract_quadratic_soc(SOC_E0):
            """
            `SOC_E0` is the SOC energy values extracted from `ref_geom_path` i.e. refG
            """
            def _compute_using_array_style(temp_dict):
                quad_ev = temp_dict["+2"] + temp_dict["-2"] - 2.0 * SOC_E0
                quad_ev /= (2. * qsize) ** 2
                return quad_ev

            def _compute_using_forloop_style(i, j, temp_dict):
                quad_ev = np.zeros(shape)
                for a in range(A):
                    quad_ev[a, a] = temp_dict["+2"][a, a] + temp_dict["-2"][a, a] - 2.0 * SOC_E0[a, a]
                    quad_ev[a, a] /= (4.0 * qsize * qsize)

                    for b in range(a):
                        quad_ev[b, a] = temp_dict["+2"][b, a] + temp_dict["-2"][b, a] - 2.0 * SOC_E0[b, a]
                        quad_ev[b, a] /= (4.0 * qsize * qsize)
                return quad_ev

            quad_dict = {}  # store return values in here
            for i in range(N):
                # ----------------------------------------------------------
                # extract Quadratic SOC from GAMESS files
                temp_dict = {}  # stores temporary arrays (different every single i \in N loop)
                for key in ['+2', '-2']:
                    assert key in linear_disp_keys, f"{key=} not in {linear_disp_keys=}"
                    path = linear_displacement_filenames[(key, i)]
                    temp_dict[key] = extract_DSOME(path, A)
                    # the extracted value is a dictionary of complex numbers indexed by the surfaces

                if False and __debug__: print(temp_dict["+2"], '\n', temp_dict["-2"]); breakpoint()

                # ----------------------------------------------------------
                # preform the math / calculation using the values we just extracted
                if not array_style:  # do it with for loops
                    soc_ev = _compute_using_forloop_style(i, temp_dict)
                else:
                    soc_ev = _compute_using_array_style(temp_dict)

                # ----------------------------------------------------------
                # store the value in the dictionary
                quad_dict[pp.mode_map_dict[i]] = soc_ev

            print("Finished extracting quadratic SOC")
            return quad_dict

        def _extract_bilinear_soc():
            """ x """
            def _compute_using_array_style(temp_dict):
                bilin_ev = np.zeros(shape, dtype=C128)
                bilin_ev += temp_dict['++'] - temp_dict['+-']
                bilin_ev += temp_dict['--'] - temp_dict['-+']
                bilin_ev /= (2. * qsize) ** 2
                return bilin_ev

            def _compute_using_forloop_style(i, j, temp_dict):
                """ only use i,j for debug print statements """
                bilin_ev = np.zeros(shape, dtype=C128)
                for a in range(A):
                    bilin_ev[a, a] += temp_dict['++'][a, a] - temp_dict['+-'][a, a]
                    bilin_ev[a, a] += temp_dict['--'][a, a] - temp_dict['-+'][a, a]
                    bilin_ev[a, a] /= (4.0 * qsize * qsize)
                    for b in range(a):
                        bilin_ev[b, a] += temp_dict['++'][b, a] - temp_dict['+-'][b, a]
                        bilin_ev[b, a] += temp_dict['--'][b, a] - temp_dict['-+'][b, a]
                        bilin_ev[b, a] /= (4.0 * qsize * qsize)
                return bilin_ev

            bilin_dict = {}
            for i, j in upper_triangle_loop_indices(N, 2):

                # ----------------------------------------------------------
                # extract BiLinear SOC from GAMESS files
                temp_dict = {}  # stores temporary arrays (different every single i,j \in N x N loop)
                for key in bi_linear_disp_keys:
                    path = bi_linear_displacement_filenames[(key, i, j)]
                    temp_dict[key] = extract_DSOME(path, A)
                    # the extracted value is a dictionary of complex numbers indexed by the surfaces

                # ----------------------------------------------------------
                # preform the math / calculation using the values we just extracted
                if not array_style:  # do it with for loops
                    soc_ev = _compute_using_forloop_style(i, j, temp_dict)
                else:
                    soc_ev = _compute_using_array_style(temp_dict)

                # ----------------------------------------------------------
                # store the value in the dictionary
                key = ij_map[(i, j)]
                bilin_dict[key] = soc_ev

            print("Finished extracting BiLinear SOC")
            return bilin_dict

        def compute_total_soc(soc_dict):
            """ compute the full spin orbit couplings """
            total_dict = {}

            # methods are just different syntax

            if True:  # method 1
                C = soc_dict['constant']  # alias
                Lin = soc_dict['Linear']  # alias
                Quad = soc_dict['Quadratic']  # alias
                BiLin = soc_dict['BiLinear']  # alias

                for i, j in upper_triangle_loop_indices(N, 2):
                    i, j = ij_map[(i, j)]  # remap to 8,9 etc..
                    # adding (A,A) shape arrays together element-wise
                    total_dict[(i, j)] = C + Lin[i] + Quad[i] + BiLin[(i, j)]
            else:  # method 1
                for i, j in upper_triangle_loop_indices(N, 2):
                    i, j = ij_map[(i, j)]  # remap to 8,9 etc..
                    total_dict[(i, j)] = soc_dict['constant']
                    total_dict[(i, j)] += soc_dict['Linear'][i]
                    total_dict[(i, j)] += soc_dict['Quadratic'][i]
                    total_dict[(i, j)] += soc_dict['BiLinear'][(i, j)]

            print("Finished calculating the Total SOC")
            return total_dict

        soc_dict = {}
        soc_dict['constant'] = extract_DSOME(ref_geom_path, A)
        soc_dict['Linear'] = _extract_linear_soc()
        soc_dict['Quadratic'] = _extract_quadratic_soc(soc_dict['constant'])
        soc_dict['BiLinear'] = _extract_bilinear_soc()

        # compute the total
        soc_dict['Total'] = compute_total_soc(soc_dict)

        return soc_dict

    # ----------------------------------------------------------

    def _write_op():

        job_title = f'{filnam} {A} states + ' + str(N) + ' modes'

        # do all the extraction first

        vibron_ev = freq_array * wn2ev
        E0_array_eV, E0_array_au = extract_E0(hessian_path)
        model = {
            "vibron eV": vibron_ev,
            "E0 eV": E0_array_eV,
            "E0 au": E0_array_au,
            "dipoles": extract_etdm(ref_geom_path),
        }

        model["Linear"] = extract_linear()
        model["Quadratic"] = extract_quadratic(E0_array_eV, E0_array_au, vibron_ev)
        model["BiLinear"] = extract_bilinear()

        if pp.SOC_flag:
            try:
                model['SOC'] = extract_soc()
            except Exception as e:
                print(str(e), "\nFailed to extract SOC! Continue?")
                breakpoint()

        file_contents = "\n".join([
            make_op_section(job_title),
            make_parameter_section(model, A, N),
            make_hamiltonian_section(model, A, N),
            make_operator_onto_dipole_moments_section(model, A, N),
            "end-operator\n"
        ])

        if False: # VECC-compatible notation if False
            for i in range(N):
                new_i = mode_map_dict[i]
                file_contents = file_contents.replace(f'v{i+1:>02d}', f'v{new_i:>02d}')
                file_contents = file_contents.replace(f'w{i+1:>02d}', f'w{new_i:>02d}')

        with open('mctdh.op', 'w') as fp:
            fp.write(file_contents)

    def confirm_necessary_files_exist():

        def _confirm_linear_are_good(bad_mode=False):
            """ confirm all displacement_input_files are good for all selected modes """
            for i in range(N):
                grace_code = {}
                for key in linear_disp_keys:

                    order = int(key[1])
                    max_order = pp.nof_displacements_per_mode[i]
                    if not (order <= max_order):
                        continue  # skip this combination

                    grace_code[key] = subprocess_call_wrapper([
                        "grep", "DONE WITH MP2 ENERGY",
                        linear_displacement_filenames[(key, i)]
                    ])
                    print(f" ..... in file {linear_displacement_filenames[(key, i)]}")

                if not all(code == 0 for code in grace_code.values()):
                    mode_label = pp.mode_map_dict[i]
                    print(f"Linear/Quad mode {mode_label} not good to extract.\n")
                    bad_mode = True

            return bad_mode

        def _confirm_bilinear_are_good(bad_mode=False):
            """ confirm all displacement_input_files are good for all selected modes """
            for i, j in upper_triangle_loop_indices(N, 2):
                grace_code = {}
                for key in bi_linear_disp_keys:
                    grace_code[key] = subprocess_call_wrapper([
                        "grep", "DONE WITH MP2 ENERGY",
                        bi_linear_displacement_filenames[(key, i, j)]
                    ])
                    print(f" ..... in file {bi_linear_displacement_filenames[(key, i, j)]}")

                if not all(code == 0 for code in grace_code.values()):
                    mode_label = (pp.mode_map_dict[i], pp.mode_map_dict[j])
                    print(f"Bilinear mode {mode_label} not good to extract.\n")
                    bad_mode = True

            return bad_mode

        def refG_file_exists():
            refG_exists = bool(subprocess_call_wrapper(["ls", ref_geom_path]) == 0)
            if not refG_exists:
                """ If refG doesn't exist, then maybe we can find the .... from other files `blah.txt`
                but also maybe they don't exist... need to check later
                """
                print(f"Skip extracting Hamiltonians from the non-existing {ref_geom_path}")
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

    return

# ---------------------------------------------------------------------------------------
# helper functions for `main()`


def read_freq_values(path):
    """ Extract frequency values from specific block of text inside file.
    Assumes `path` points to hessian output file.
    """
    selected_lines = extract_lines_between_patterns(
        path,
        "FREQUENCIES IN CM",
        "REFERENCE ON SAYVETZ CONDITIONS"
    )

    freq_value_set = []
    for freqline in selected_lines:
        if "FREQUENCY:" in freqline:
            freq_value_set.append(freqline[18:])

    return freq_value_set


def read_mode_values(path):
    """ Extract filtered set of lines from specific block of text inside file.
    Assumes `path` points to hessian output file.
    """
    selected_lines = extract_lines_between_patterns(
        path,
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


def _extract_freq_and_mode_from_hessian(path):
    """ Extracts the frequency and normal mode information from the Hessian.
    Writes those values to two individual text files.
    """
    filtered_set = read_mode_values(path)
    with open('mode.dat', 'w') as fp:
        fp.writelines(filtered_set)

    freq_value_set = read_freq_values(path)
    with open('freq.dat', 'w') as fp:
        fp.writelines(freq_value_set)
    return


def get_number_of_atoms(path):
    """ Extract number of atoms from specific block of text inside file.
    Assumes `path` points to hessian output file.
    """

    # would be good to replace this with memory mapping find or grep command?
    with open(path, 'r', errors='replace') as hess_file:
        for line in hess_file:
            if ' TOTAL NUMBER OF ATOMS' in line:
                natoms = int(line.split('=')[1])
                return natoms


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

    freq_list = [[float(n) for n in line.strip().replace('I', '').split()] for line in lines_freq]
    frequences = list(it.chain(*freq_list))  # single list
    freq_array = np.array(frequences)

    # -------------------------------------------------------------------------
    if False and __debug__:   # print all frequencies
        string = "\n".join([f"frequency: {i} {freq_array[i]} CM-1" for i in range(ndim)])
        print(string)

    return modes_array, freq_array


def compose_ref_structure(ref_geom_path, hessian_path, nof_atoms):
    """ x """

    start = 'EQUILIBRIUM GEOMETRY LOCATED'
    end = 'INTERNUCLEAR DISTANCES'
    coord_lines = extract_lines_between_patterns(hessian_path, start, end)

    good_ref_structure = bool(len(coord_lines) > 2)

    if not good_ref_structure:
        print(
            f'Unsuccessful extraction of equilibrium geometry from {hessian_path}.\n'
            'Please prepare ref_structure manually.'
        )
        breakpoint()  # do we want to continue execution, or do we actually want to stop the program?
        import sys; sys.exit()

    """ we don't need to delete the file if we simply write a new file
    _delete_file_using_rmrf(ref_geom_path)
    """

    # the last element of coord_lines is an empty line (`\n`)
    assert coord_lines[-1] == '\n'

    # remove empty lines
    coord_lines = [l for l in coord_lines if l != '\n']

    # we want the lines at the end of the file (the last nof_atoms/Z lines)
    file_contents = "".join(coord_lines[-nof_atoms:])

    with open(ref_geom_path, 'w') as fp:
        fp.write(file_contents)

    print(
        f'Successfully extracted equilibrium geometry from\n {hessian_path=}\n'
        f"And prepared reference structure at\n {ref_geom_path=}."
    )

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


def refG_calc(ref_geom_path, **kwargs):
    """
    Do diabatization calculation at the reference non-distorted structure.
    This calculation shall be a repetition of a calculation in preparing `temp.inp`.
    """

    # input and output paths for the reference geometry calculation
    input_path = kwargs['refG_in']
    output_path = kwargs['refG_out']

    # Check if the calculation has already been run
    grace_exists = subprocess_call_wrapper(["grep", "DONE WITH MP2 ENERGY", output_path]) == 0
    if grace_exists:
        print("Calculation at the reference structure has already been done.")
        return
    print("About to run calculation at the undistorted reference structure")

    # create the input file for the reference geometry calculation
    # the first part of the file comes from `temp.inp`
    shutil.copy("temp.inp", input_path)

    # then we append the undistorted reference structure data onto the `input_path`
    with open(ref_geom_path, 'r', errors='replace') as ref_structure:
        data = ref_structure.read()
    with open(input_path, "a") as fp:
        fp.write(data)
        fp.write(" $END\n")

    # Finally we submit and run the refG calculation (you may need to customize this command based on your setup)
    # refG_job_result = subprocess_run_wrapper(["./subgam.diab", input_path, "4", "0", "1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    job_path = my_subgam(input_path, ncpus=2, ngb=1, nhour=1)
    os_system_wrapper(f"sbatch -W {job_path}")

    # At this point, refG calculation has completed successfully.
    print("Calculation at the reference structure is done.")

    return

# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------


def process_profiling_data(filename):
    """ temporary formatted printing of profiling data """

    p = pstats.Stats(filename)

    # prints the 20 functions with the largest cumulative runtime
    p.strip_dirs().sort_stats("cumulative").print_stats(20)

    # prints the 20 functions with the largest total runtime
    p.strip_dirs().sort_stats("tottime").print_stats(20)

    if True:  # this section is for analyzing specific functions

        # this sorts by cumulative runtime then only prints the stats for functions with the substring 'calculate' in their name
        p.strip_dirs().sort_stats("cumulative").print_stats('extract', 15)

        # `label` substring in their name
        p.strip_dirs().sort_stats("cumulative").print_stats('label', 15)

        # `build` substring in their name
        p.strip_dirs().sort_stats("cumulative").print_stats('build', 15)

        function_name = "_write_op"  # inside mctdh()
        p.strip_dirs().print_stats(function_name)  # this prints the stats for the specific function `function_name`
        p.strip_dirs().print_callers(function_name)  # this prints the functions that call `function_name`

        # # the same as above but for a different function
        function_name = "diabatization"
        p.strip_dirs().print_stats(function_name)
        p.strip_dirs().print_callers(function_name)  # `function_name` <== list of fxn

        # these two functions are called a lot BY other functions
        function_name = "_extract_energy_from_gamessoutput"
        p.strip_dirs().print_stats(function_name)
        p.strip_dirs().print_callees(function_name)  # `function_name` ==> list of fxn
        function_name = "extract_DSOME"
        p.strip_dirs().print_stats(function_name)
        p.strip_dirs().print_callees(function_name)

        function_name = "extract_DSOME"
        p.strip_dirs().print_callees(function_name)

        # # this does the same as above, but the callees are sorted by cumulative
        # function_name = "d"
        # p.strip_dirs().sort_stats("cumulative").print_callees(function_name)
    return


# ---------------------------------------------------------------------------------------
def main(ref_geom_path="ref_structure", ncols=5, **kwargs):
    """ x """
    from project_parameters import Z, N_tot

    hessian_filename = kwargs['hessian_filename']
    _extract_freq_and_mode_from_hessian(hessian_filename)
    nof_atoms = natoms = get_number_of_atoms(hessian_filename)
    assert Z == nof_atoms, f"{Z=} is not {nof_atoms=}!? Check hessian and project_parameters!"

    assert N_tot == (nof_atoms * 3), f"{N_tot=} is not {nof_atoms * 3=}!?  Check hessian and project_parameters!"

    normal_modes_array, frequencies_cm = process_mode_freq(N_tot, ncols)

    compose_ref_structure(ref_geom_path, hessian_filename, nof_atoms)

    atom_dict, charge_dict, reference_coordinates = read_reference_structure(ref_geom_path)

    refG_calc(ref_geom_path, **kwargs)

    # -------------------------------------------------------------------------
    diabatization_kwargs = kwargs.copy()
    diabatization_kwargs.update({
        'all_freq': frequencies_cm,
        'reference_coordinates': reference_coordinates,
        'normal_modes': normal_modes_array,
        'atom_dict': atom_dict,
        'charge_dict': charge_dict,
    })

    disp_coord = diabatization(**diabatization_kwargs)
    # turns out we don't actually use the `disp_coord` for anything?

    print("Diabatization successfully modified?")

    # -------------------------------------------------------------------------
    op_file_name = "mctdh.op"
    op_path = join("./", op_file_name)

    hessian_path = kwargs['hessian_filename']

    fitting()

    if A == 1:
        return

    mctdh(op_path, hessian_path, frequencies_cm, pp.A, pp.N, **kwargs)

    print(f"{op_path=} successfully modified\n")

    # -------------------------------------------------------------------------
    # copy <mctdh.op> file to <specific_file_name.op>
    src_path = op_path
    dst_path = join("./", pp.file_name + '.op')
    shutil.copy(src_path, dst_path)

    print(
        " "*4 + src_path,
        "copied to",
        " "*4 +dst_path,
        sep='\n'
    )

    # -------------------------------------------------------------------------
    if False and __debug__:
        header = f"\n{'-'*20}{{}}{'-'*20}\n"
        print_header = lambda s: print(header.format(s))

        print_header('Normal modes')
        for i in range(nrmmod.shape[0]):
            print(" "*4, i)
            pprint.pprint(nrmmod[i, :])

        print_header('Frequencies')
        for i in range(freqcm.shape[0]):
            print(" "*4, i, freqcm[i])

        print_header('Execution Parameters')
        print(f"Selected modes: {pp.selected_mode_list}")

        print("List of atoms:")
        for k, v in atom_dict.items():
            print(" "*4, k, v)

        print("List of charges:")
        for k, v in charge_dict.items():
            print(" "*4, k, v)

        print("Reference co-ordinates")
        for k, v in ref_coords.items():
            print(" "*4, k, v)

    return


# ---------------------------------------------------------------------------------------
if (__name__ == "__main__"):

    if len(sys.argv) != 2:
        print("Usage: python your_script.py <path_to_hessian_output>")
        sys.exit(1)

    hessian_filename = sys.argv[1]  # read in name of hessian file
    kwargs = {'hessian_filename': hessian_filename}

    # everything we add to `kwargs` is 'imported' from project parameters

    kwargs.update({
        'refG_in': f"{pp.filnam}_refG.inp",  # reference geometry
        'refG_out': f"{pp.filnam}_refG.out",  # reference geometry
        # 'modes_included': pp.modes_included,
    })

    # ---------------------------------------------------------------
    profiling = True  # just change this to enable profiling

    if not profiling:
        main(**kwargs)

    else:
        root = os.getcwd()
        filename = join(root, "cProfile_dist_allmodes_pm")

        if True:  # set this to false if you simply want to print out the profile stats again (without running all the code)
            cProfile.runctx(
                'main(**kwargs)',
                globals(),
                locals(),
                filename
            )

        # print the results of profiling to stdout
        # optional (can always be called by some other script)
        process_profiling_data(filename)

        if True:  # if you want to save the results to a file
            from contextlib import redirect_stdout  # to send the prints to a file
            with open(filename+'.txt', 'w') as f:
                with redirect_stdout(f):
                    process_profiling_data(filename)
    # ---------------------------------------------------------------
