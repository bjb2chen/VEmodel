""" this module reads in the relevant data from a (specific) *.op file and produces a corresponding *.json file as used in the package

The format of the *.op file which describes a vibronic model of a molecule is specific in that it is symmetric in the modes by definition as well as only lists coupling terms associated with two different modes, q_{3}q_{2} or q_{4}q_{3}^{2}, q_{12}^{4}

Bi cubic has q_{i}q_{j}^{2} and q_{i}^{2}q_{j}
Bi quartic has q_{i}q_{j}^{3}, q_{i}^{3}q_{j} and q_{i}^{2}q_{j}^{2}

All numerical values in the *.op file are represented in electron volts (eV)
"""

# system imports
import itertools as it
import functools
import warnings
import mmap
import os  # needed to check unix/windows for mmap flags

# third party imports
import numpy as np
from numpy import float64 as F64
from numpy import complex128 as C128
import parse

# local imports
from ..log_conf import log
from . import helper
from .helper import StringNotFoundError
from .vibronic_model_keys import VibronicModelKeys as VMK

headers = [
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


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------


def extract_energies(path, memmap):
    """x"""
    # find the beginning and ending of the important region
    memmap.seek(0)  # start looking from the beginning of the file
    beginString = headers[1]
    begin = helper.find_string_in_file(memmap, path, beginString)
    endString = 'Electronic transition moments'
    end = helper.find_string_in_file(memmap, path, endString)

    # go to the beginning of that region
    memmap.seek(begin)

    # skip the header
    helper.readlines(memmap, 3)

    # read all the relevant data
    byteData = memmap.read(end - memmap.tell())
    stringData = byteData.decode(encoding="utf-8")
    lines = [
        line for line in stringData.strip().splitlines() if "#" not in line and line != ""
    ]

    # the number of electronic states
    A = len(lines)

    # save the reference Hamiltonian into the energies array
    energies = np.zeros((A, A))

    # this needs to be changed to support off-diagonal energy values
    for a in range(A):  # old style
        list_of_words = lines[a].split()
        assert list_of_words[0] == f"EH_s{a+1:02}_s{a+1:02}", f"{list_of_words=}\nEH_s{a+1:02}_s{a+1:02}"
        assert list_of_words[-1] == "ev", f"{list_of_words[-1]=}"
        energies[a, a] = list_of_words[2]

    # most likely we would use parse to process the lines
    # this assumes that the first index is the row index and the second the column index
    if False:
        parse_pattern = "EH_s{a1:d}_s{a2:d}"
        p = parse.compile(parse_pattern)
        for line in lines:
            if line is None:  # skip all empties
                continue
            if "EH" not in line:  # expect EH prefix
                print(f"Line is malformed, missing EH?\n{line=}"); breakpoint()
                raise Exception()

            list_of_words = line.split()
            assert list_of_words[-1] == "ev", f"Can only handle energy units in eV not {list_of_words[-1]=}"

            r = p.parse(list_of_words[0])
            assert r != None, f"Failed to parse\n{line=}\nInstead got {r=}? Line probably doesn't match the {parse_pattern=}"
            try:
                a1, a2 = r['a1'], r['a2']
                energies[a1, a2] = float(list_of_words[2])
            except TypeError as e:
                print(r, line); breakpoint()
                raise (str(e))


    return energies, A


def extract_normal_mode_frequencies(path, memmap):
    """store output in frequency_array"""
    # find the beginning and ending of the important region
    try:
        # find the beginning and ending of the important region
        memmap.seek(0)  # start looking from the beginning of the file
        beginString = headers[0]
        begin = helper.find_string_in_file(memmap, path, beginString)
        # TODO - should check have a default check for these three:
        # > EOF
        # > --------------------------
        for endString in headers[1:]:
            try:
                end = helper.find_string_in_file(memmap, path, endString)
            except Exception as e:
                if f"It seems \"{endString:s}\" was not present in the file\n" in e.args[0]:
                    continue
            else:
                break
    except Exception as e:
        log.warning(f"Couldn't find {beginString}")
        return

    # go to the beginning of that region
    memmap.seek(begin)

    # skip headers
    helper.readlines(memmap, 2)

    # read all the relevant data
    byteData = memmap.read(end - memmap.tell())
    stringData = byteData.decode(encoding="utf-8")
    lines = [line for line in stringData.strip().splitlines() if "#" not in line and line != ""]

    # the number of normal modes
    N = len(lines)

    # extract the numbers and save them in the frequencies array
    frequencies = np.zeros(N)
    for j in range(N):
        list_of_words = lines[j].split()
        assert list_of_words[0] == f"w{j+1:02}"
        assert list_of_words[-1] == "ev"
        frequencies[j] = list_of_words[2]

    return frequencies, N


def find_byte_begin_and_end(path, memmap, begin_string, end_string_list):
    """ attempts to find the index(in bytes) of the start and end of the important region, which are denoted by the begin_string and end_string
    """
    try:
        memmap.seek(0)  # start looking from the beginning of the file
        begin = helper.find_string_in_file(memmap, path, begin_string)
        # TODO - should check have a default check for these three:
        # > EOF
        # > --------------------------
        for end_string in end_string_list:
            try:
                end = helper.find_string_in_file(memmap, path, end_string)
            except StringNotFoundError as e:
                continue
            else:
                break
        else:
            raise Exception("File is improperly formatted, couldn't find end string of {begin:s}")

    except StringNotFoundError as err:
        log.warning(f"Couldn't find {begin_string}")
        raise err

    return begin, end


def parse_biCubic_special_case(lines, coupling_terms):
    """handle the parsing of biCubic parameters which require a bit more work"""
    p = parse.compile("B3_s{a1:d}s{a2:d}_v{j1:d}v{j2:d}")

    for line in lines:
        r = p.parse(line[0])

        # if the first mode is smaller then it is q, q^2
        if r['j1'] < r['j2']:
            index_tuple = (r['j1']-1, r['j2']-1, r['j2']-1, r['a1']-1, r['a2']-1)
        # if the first mode is bigger then it is q^2, q
        elif r['j1'] > r['j2']:
            index_tuple = (r['j1']-1, r['j1']-1, r['j2']-1, r['a1']-1, r['a2']-1)

        coupling_terms[index_tuple] = line[1]
    return


def parse_biQuartic_special_case(lines, coupling_terms):
    """handle the parsing of biQuartic parameters which require a lot more work"""
    list_B4 = filter(lambda item: item[0].startswith("B4"), lines)
    pB4 = parse.compile("B4_s{a1:d}s{a2:d}_v{j1:d}v{j2:d}")

    for line in list_B4:
        r = pB4.parse(line[0])

        # if the first mode is smaller then it is q, q^3
        if r['j1'] < r['j2']:
            index_tuple = (r['j1']-1, r['j2']-1, r['j2']-1, r['j2']-1, r['a1']-1, r['a2']-1)
        # if the first mode is bigger then it is q^3, q
        elif r['j1'] > r['j2']:
            index_tuple = (r['j1']-1, r['j1']-1, r['j1']-1, r['j2']-1, r['a1']-1, r['a2']-1)

        coupling_terms[index_tuple] = line[1]

    list_A4 = filter(lambda item: item[0].startswith("A4"), lines)
    pA4 = parse.compile("A4_s{a1:d}s{a2:d}_v{j1:d}v{j2:d}")

    # A4's are always q^2 q^2
    for line in list_A4:
        r = pA4.parse(line[0])
        index_tuple = (r['j1']-1, r['j1']-1, r['j2']-1, r['j2']-1, r['a1']-1, r['a2']-1)
        coupling_terms[index_tuple] = line[1]
    return


def parse_lines(lines, coupling_terms, order=None):
    """extracts the numerical values of linear coupling terms from a list of strings"""
    assert order is not None, "Must provide an order parameter!"

    # If an empty list then those terms don't exist
    if not lines:
        return

    if order == 'C1':
        p = parse.compile("C1_s{a1:d}_s{a2:d}_v{j:d}")
        make_index_tuple = lambda r: (r['j']-1, r['a1']-1, r['a2']-1)
    elif order == 'C2':
        p = parse.compile("C2_s{a1:d}s{a2:d}_v{j1:d}v{j2:d}")
        make_index_tuple = lambda r: (r['j1']-1, r['j2']-1, r['a1']-1, r['a2']-1)
    elif order == 'C3':
        p = parse.compile("C3_s{a1:d}_s{a2:d}_v{j:d}")
        make_index_tuple = lambda r: (r['j']-1, r['j']-1, r['j']-1, r['a1']-1, r['a2']-1)
    elif order == 'C4':
        p = parse.compile("C4_s{a1:d}_s{a2:d}_v{j:d}")
        make_index_tuple = lambda r: (r['j']-1, r['j']-1, r['j']-1, r['j']-1, r['a1']-1, r['a2']-1)
    elif order == 'C5':
        p = parse.compile("C5_s{a1:d}_s{a2:d}_v{j:d}")
        make_index_tuple = lambda r: (r['j']-1, r['j']-1, r['j']-1, r['j']-1, r['j']-1,
                                      r['a1']-1, r['a2']-1)
    elif order == 'C6':
        p = parse.compile("C6_s{a1:d}_s{a2:d}_v{j:d}")
        make_index_tuple = lambda r: (r['j']-1, r['j']-1, r['j']-1, r['j']-1, r['j']-1, r['j']-1,
                                      r['a1']-1, r['a2']-1)
    elif order == 'B3':
        print("Warning, check biCubic results")
        return parse_biCubic_special_case(lines, coupling_terms)
    elif order == 'B4':
        print("Warning, check biQuartic results")
        return parse_biQuartic_special_case(lines, coupling_terms)
    else:
        raise Exception("Order parameter can only take on the values {C1,C2,C3,C4,C5,C6,B3,B4,A4}")

    for line in lines:
        if line is None:
            continue
        r = p.parse(line[0])
        assert r != None, f"Failed to parse\n{line=}\nInstead got {r=}? Line probably doesn't match any parse patterns"
        try:
            index_tuple = make_index_tuple(r)
        except TypeError as e:
            print(r, line); breakpoint()
            raise (str(e))
        coupling_terms[index_tuple] = line[1]

    return


def extract_string_list(path, memmap, header_index=2):
    """abstract function for extracting a list of strings for each line in a specific section of the *.op file corresponding to the memmap object"""

    begin_string = headers[header_index]
    end_string_list = headers[header_index + 1:]
    try:
        begin, end = find_byte_begin_and_end(path, memmap, begin_string, end_string_list)
    except StringNotFoundError as err:
        return list()

    # go to the beginning of that region
    memmap.seek(begin)

    # skip headers
    helper.readlines(memmap, 2)

    # get the important bytes
    byteData = memmap.read(end - memmap.tell())

    # decode it
    stringData = byteData.decode(encoding="utf-8")

    # clean it up
    stringData = stringData.strip().replace('=', '').replace(', ev', '')

    # return the list of strings (not including empty strings and comments)
    return [line.split() for line in stringData.splitlines() if "#" not in line and line != ""]


def extract_electronic_transition_dipole_moments(path, memmap, dipole_moments):
    """calls extract_string_list() with appropriate parameters
    so as to fill the dipole moments array with values from the *.op file

    `dipole_moments` is an array of shape (C, A+1)
        where C is the number of co-ordinates representing each atom (x,y,z => C = 3)
        and the second dimension is the number of surfaces PLUS 1 (representing a fictitious ground state)
    """
    idx = headers.index('Electronic transition moments')
    lines = extract_string_list(path, memmap, header_index=idx)

    for line in lines:
        if 'au' in line:
            print(lines)
            #r = input("Does not support reading in units of A.U. / Hartrees\n Proceed ignoring units?")
            #if r.lower() in ['y', 'yes']:
            #    break
            #else:
            #    raise Exception()

    # define dimensions
    C_dim = coordinate_dimension = dipole_moments.shape[0]
    A_dim = excited_electronic_dimension = dipole_moments.shape[1] - 1

    # double check
    nof_lines_expected = C_dim * A_dim
    assert len(lines) == nof_lines_expected, (
        f"The file is malformed check the path {path}.\n"
        f"We expected {nof_lines_expected=} lines but instead got {len(lines)}."
        f"\nWe should have {excited_electronic_dimension=} values for each of the {coordinate_dimension=} dimensions"
        "\nDid you not specify the correct number of dimensions?"
    )

    # read in TDM values
    for c in range(coordinate_dimension):
        for a in range(excited_electronic_dimension):
            label, number, *_ = lines[a + (c*A_dim)]
            assert ('Ex' in label) or ('Ey' in label) or ('Ez' in label), f"{label=}{number=} does not contain Ex/Ey/Ez ?!\n{lines=}"
            assert label[2:] == f"_s00_s{a+1:02}", f"Label is malformed? We expect \n{label[0:2]}_s00_s{a+1:02}\nbut got\n{label=}"
            dipole_moments[c, a] = complex(number)

    return


def extract_magnetic_transition_dipole_moments(path, memmap, dipole_moments):
    """calls extract_string_list() with appropriate parameters
    so as to fill the dipole moments array with values from the *.op file

    `dipole_moments` is an array of shape (C, A+1)
        where C is the number of co-ordinates representing each atom (x,y,z => C = 3)
        and the second dimension is the number of surfaces PLUS 1 (representing a fictitious ground state)
    """
    idx = headers.index('Magnetic transition moments')
    lines = extract_string_list(path, memmap, header_index=idx)

    for line in lines:
        if 'au' in line:
            print(lines)
            #raise Exception("Does not support reading in units of A.U. / Hartrees")

    # define dimensions
    C_dim = coordinate_dimension = dipole_moments.shape[0]
    A_dim = excited_electronic_dimension = dipole_moments.shape[1] - 1

    # double check
    nof_lines_expected = C_dim * A_dim
    assert len(lines) == nof_lines_expected, (
        f"The file is malformed check the path {path}.\n"
        f"We expected {nof_lines_expected=} lines but instead got {len(lines)}."
        f"\nWe should have {excited_electronic_dimension=} values for each of the {coordinate_dimension=} dimensions"
        "\nDid you not specify the correct number of dimensions?"
    )

    for c in range(coordinate_dimension):
        for a in range(excited_electronic_dimension):
            label, number, *_ = lines[a + (c*A_dim)]
            assert ('Mx' in label) or ('My' in label) or ('Mz' in label), f"{label=}{number=} does not contain Mx/My/Mz ?!\n{lines=}"
            assert label[2:] == f"_s00_s{a+1:02}", f"Label is malformed? We expect \n{label[0:2]}_s00_s{a+1:02}\nbut got\n{label=}"
            dipole_moments[c, a] = complex(number)

    return


def extract_linear_couplings(path, memmap, linear):
    """calls extract_string_list() with appropriate parameters so as to fill the linear coupling term array with values from the *.op file"""
    idx = headers.index('Linear Coupling Constants')
    lines = extract_string_list(path, memmap, header_index=idx)
    parse_lines(lines, linear, order='C1')
    return


def extract_quadratic_couplings(path, memmap, quadratic):
    """calls extract_string_list() with appropriate parameters so as to fill the quadratic coupling term array with values from the *.op file"""
    idx = headers.index('Diagonal Quadratic Coupling Constants')
    lines = extract_string_list(path, memmap, header_index=idx)
    parse_lines(lines, quadratic, order='C2')
    return


def extract_offdiagonal_quadratic_couplings(path, memmap, quadratic):
    """calls extract_string_list() with appropriate parameters so as to fill the quadratic coupling term array with values from the *.op file"""
    idx = headers.index('Off_diagonal Quadratic Coupling Constants')
    lines = extract_string_list(path, memmap, header_index=idx)
    parse_lines(lines, quadratic, order='C2')
    return


def extract_cubic_couplings(path, memmap, cubic):
    """calls extract_string_list() with appropriate parameters so as to fill the cubic coupling term array with values from the *.op file"""
    idx = headers.index('Cubic Coupling Constants')
    lines = extract_string_list(path, memmap, header_index=idx)
    parse_lines(lines, cubic, order='C3')
    return


def extract_bicubic_couplings(path, memmap, bicubic):
    """calls extract_string_list() with appropriate parameters so as to fill the bicubic coupling term array with values from the *.op file"""
    idx = headers.index('Bi-Cubic Constants')
    lines = extract_string_list(path, memmap, header_index=idx)
    parse_lines(lines, bicubic, order='B3')
    return


def extract_quartic_couplings(path, memmap, quartic):
    """calls extract_string_list() with appropriate parameters so as to fill the quartic coupling term array with values from the *.op file"""
    idx = headers.index('Quartic Coupling Constants')
    lines = extract_string_list(path, memmap, header_index=idx)
    parse_lines(lines, quartic, order='C4')
    return


def extract_biquartic_couplings(path, memmap, biquartic):
    """calls extract_string_list() with appropriate parameters so as to fill the biquartic coupling term array with values from the *.op file"""
    idx = headers.index('Bi-Quartic Constants')
    lines = extract_string_list(path, memmap, header_index=idx)
    parse_lines(lines, biquartic, order='B4')
    return


def surface_symmetric(*coupling_terms):
    """Returns True if all the provided coupling terms are symmetric in the surfaces, otherwise False."""

    transposition_dictionary = {
        # linear
        1: ([0, 2, 1]),
        # quadratic
        2: ([0, 1, 3, 2]),
        # cubic
        3: ([0, 1, 2, 4, 3]),
        # quartic
        4: ([0, 1, 2, 3, 5, 4]),
    }

    for term in coupling_terms:

        if term is None:
            continue

        order = len(term.shape) - 2
        assert 4 >= order >= 1, "Incorrect dimensions in provided numpy array"
        surfaces = transposition_dictionary[order]

        # check for symmetry in surfaces
        if not np.allclose(term, term.transpose(surfaces)):
            return False

    return True


def surface_symmetrize_coupling_terms(Modes, States, *coupling_terms):
    """Makes sure the four coupling terms passed in are symmetric (along the electronic surfaces)
    if they are symmetric then it simply returns without doing anything
    otherwise it checks to see if it needs to copy the lower triangle to the upper triangle."""

    for term in coupling_terms:

        if term is None:
            continue

        if surface_symmetric(term):
            continue

        # note that this overcounts the diagonal by a factor of 2
        for a, b in it.product(States, States):
            term[..., a, b] += term[..., b, a]

        # so now we correct for that factor of 2
        for a in States:
            term[..., a, a] /= 2.0


def double_quadratic_terms(number_of_modes, States, quadratic_terms):
    """If the quadratic terms are symmetric in the modes, or have an upper triangle of all zeros then we multiple all quadratic factors by 2 to account for the 1/2 factor in our mathematical definition.
    This option leaves the elements as they were, and does not attempt to symmetrize the model.
    """

    upper_triangle_idx = np.triu_indices(number_of_modes, k=1)

    for a, b in it.product(States, States):
        if not np.all(quadratic_terms[(*upper_triangle_idx, a, b)] == 0.0):
            raise Exception(
                f"The upper triangle at a({a}), b({b}) is not all zeros, "
                "therefore we cannot just double the quadratic terms\n"
            )
    else:
        quadratic_terms[:] *= 2.0


def mode_symmetrize_from_lower_triangle_quadratic_terms(number_of_modes, States, quadratic_terms):
    """If the quadratic terms are zero in the upper triangle then we copy the lower triangle to the upper triangle and multiply the diagonal terms by 2.
    This assumes the values are in the lower triangle.
    Note that the order of upper and lower triangle indices are slightly different so we need to reverse.
    """

    upper_triangle_idx = np.triu_indices(number_of_modes, k=1)
    lower_triangle_idx = np.tril_indices(number_of_modes, k=-1)
    diagonal_idx = np.diag_indices(number_of_modes)

    for a, b in it.product(States, States):
        if not np.all(quadratic_terms[(*upper_triangle_idx, a, b)] == 0.0):
            raise Exception(
                f"The upper triangle at a({a}), b({b}) is not all zeros, "
                "therefore we cannot symmetrize along the modes\n"
            )
    else:
        for a, b in it.product(States, States):
            quadratic_terms[(*reversed(lower_triangle_idx), a, b)] = quadratic_terms[(*lower_triangle_idx, a, b)]


def mode_symmetrize_from_upper_triangle_quadratic_terms(number_of_modes, States, quadratic_terms):
    """If the quadratic terms are zero in the upper triangle then we copy the lower triangle to the upper triangle and multiply the diagonal terms by 2.
    This assumes the values are in the lower triangle.
    Note that the order of upper and lower triangle indices are slightly different so we need to reverse.
    """

    upper_triangle_idx = np.triu_indices(number_of_modes, k=1)
    lower_triangle_idx = np.tril_indices(number_of_modes, k=-1)

    for a, b in it.product(States, States):
        if not np.all(quadratic_terms[(*lower_triangle_idx, a, b)] == 0.0):
            raise Exception(
                f"The lower triangle at a({a}), b({b}) is not all zeros, "
                "therefore we cannot symmetrize along the modes\n"
            )
    else:
        for a, b in it.product(States, States):
            quadratic_terms[(*reversed(upper_triangle_idx), a, b)] = quadratic_terms[(*upper_triangle_idx, a, b)]


def read_model_op_file(
        path_file_op,
        surface_symmetrize=True,
        double_quadratic=False,
        symmetrize_quadratic=False,
        highest_order=1,
        get_transition_dipole_moment=True,
        dimension_of_dipole_moments=1
):
    """Reads/parses molecule_vibron.op file and returns a dictionary in the standard format defined in the package."""

    # declare the arrays used to store the model's parameters
    excitation_energies = None
    frequencies = None
    linear = None
    quadratic = None
    cubic = None
    quartic = None

    # we will overwrite these default values
    A = N = States = Modes = 0
    size = {}

    helper.verify_file_exists(path_file_op)

    # store the energy offsets, and all the coupling terms
    with open(path_file_op, "r+b") as source_file:

        """ Turns out the memory protection arguments for mmap are different
        for windows/unix,  so this is how we handle that
        """
        running_on_windows = bool(os.name == 'nt')
        read_write_prot_kwarg = {}
        if running_on_windows:
            read_write_prot_kwarg["access"] = mmap.ACCESS_READ
        else:
            read_write_prot_kwarg["prot"] = mmap.PROT_READ

        # access the file using memory map for efficiency
        with mmap.mmap(source_file.fileno(), 0, **read_write_prot_kwarg) as mm:

            frequencies, N = extract_normal_mode_frequencies(path_file_op, mm)
            excitation_energies, A = extract_energies(path_file_op, mm)
            States = range(A)
            Modes = range(N)

            size = {
                'N': (N),
                'A': (A),
                'AA': (A, A),
                'NAA': (N, A, A),
                'NNAA': (N, N, A, A),
                'NNNAA': (N, N, N, A, A),
                'NNNNAA': (N, N, N, N, A, A),
            }

            # Assert/Initialize array sizes
            assert frequencies.shape[0] == size['N'], "Incorrect array dimensions"
            assert excitation_energies.shape == size['AA'], "Incorrect array dimensions"

            # predefine these arrays now that we know the values of N and A are
            assert 4 > dimension_of_dipole_moments > 0, 'not supported for dim > 4 or dim == 0'
            electronic_dipole_moments = np.zeros((dimension_of_dipole_moments, A), dtype=C128)
            magnetic_dipole_moments = np.zeros((dimension_of_dipole_moments, A), dtype=C128)

            linear, quadratic, cubic, quartic = None, None, None, None

            if highest_order >= 1:
                linear = np.zeros(size['NAA'], dtype=F64)
            if highest_order >= 2:
                quadratic = np.zeros(size['NNAA'], dtype=F64)
            if highest_order >= 3:
                cubic = np.zeros(size['NNNAA'], dtype=F64)
            if highest_order >= 4:
                quartic = np.zeros(size['NNNNAA'], dtype=F64)

            # special case
            if get_transition_dipole_moment:
                extract_electronic_transition_dipole_moments(path_file_op, mm, electronic_dipole_moments)
                extract_magnetic_transition_dipole_moments(path_file_op, mm, magnetic_dipole_moments)

            # read in the rest of the parameters
            if highest_order >= 1:
                extract_linear_couplings(path_file_op, mm, linear)
            if highest_order >= 2:
                extract_quadratic_couplings(path_file_op, mm, quadratic)
            if highest_order >= 2:
                extract_offdiagonal_quadratic_couplings(path_file_op, mm, quadratic)
            if highest_order >= 3:
                extract_cubic_couplings(path_file_op, mm, cubic)
            if highest_order >= 4:
                extract_quartic_couplings(path_file_op, mm, quartic)
            if highest_order >= 3:
                extract_bicubic_couplings(path_file_op, mm, cubic)
            if highest_order >= 4:
                extract_biquartic_couplings(path_file_op, mm, quartic)

    if surface_symmetrize:
        surface_symmetrize_coupling_terms(Modes, States, linear, quadratic, cubic, quartic)

        assert np.allclose(excitation_energies, excitation_energies.transpose(1, 0))
        #  assert surface_symmetric(linear, quadratic, cubic, quartic), (
        #      "Not all coupling terms are symmetric along the surface dimensions (the last 2)\n"
        #  )

    if highest_order >= 2:
        if double_quadratic:
            double_quadratic_terms(N, States, quadratic)
        elif symmetrize_quadratic:
            mode_symmetrize_from_lower_triangle_quadratic_terms(N, States, quadratic)
        else:
            log.warning("We didn't change the quadratic terms in any way, make sure this was intentional!!!")

    # and we are done
    maximal_dict = {VMK.N: N,
                    VMK.A: A,
                    VMK.E: excitation_energies,
                    VMK.w: frequencies,
                    VMK.etdm: electronic_dipole_moments,
                    VMK.mtdm: magnetic_dipole_moments,
                    }

    if highest_order >= 1:
        maximal_dict[VMK.G1] = linear
    if highest_order >= 2:
        maximal_dict[VMK.G2] = quadratic
    if highest_order >= 3:
        maximal_dict[VMK.G3] = cubic
    if highest_order >= 4:
        maximal_dict[VMK.G4] = quartic

    # if the arrays only have zeros then we might not need to store them?
    return_dict = dict((k, v) for k, v in maximal_dict.items() if not np.all(v == 0))

    return return_dict


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# bad practice, works for now and we can refactor once we've finished figuring out the end product
format_string = "{label:<25s}={value:>-15.6f}{units:>8s}\n"
make_line = functools.partial(format_string.format, units=", ev")


def build_op_section():
    """Returns a string which defines the `OP_DEFINE-SECTION` of an .op file"""

    start = "OP_DEFINE-SECTION"
    end = "end-op_define-section"

    title = "template title"
    return '\n'.join([
        start,
        "title",
        f"{title:s}",
        "end-title",
        end,
    ])


def build_frequencies(model):
    """Return a string containing the frequency information of a .op file."""
    return ''.join([
        make_line(label=f"w{idx+1:0>2d}", value=w)
        for idx, w in enumerate(model[VMK.w])
    ])


def build_energies(model, diagonal):
    """Return a string containing the energy information of a .op file."""
    if diagonal:
        return ''.join([
            make_line(label=f"EH_s{idx+1:0>2d}_s{idx+1:0>2d}", value=E)
            for idx, E in enumerate(model[VMK.E])
        ])
    else:
        return ''.join([
            make_line(label=f"EH_s{idx+1:0>2d}_s{idx+1:0>2d}", value=E)
            for idx, E in enumerate(np.diag(model[VMK.E]))
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


def build_magnetic_moments(model):
    """Return a string containing the magnetic transition dipole moment information of a .op file."""

    # this is how many dimensions you have in your co-ordinate system to describe the dipole
    # (i.e. the dimensionality of your dipole's vector, not accounting for the number of surfaces)
    dipole_dimension = int(model[VMK.mtdm].shape[0])

    assert 4 > dipole_dimension > 0, "Current hardcoded implementation does not support 4 dimensional dipoles"

    # we ignore the complex warning from casting the mtdm (complex) to a float
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", np.ComplexWarning)

        string = ''

        for d in range(dipole_dimension):
            name = {0: 'Mx', 1: 'My', 2: 'Mz'}[d]
            string += ''.join([
                make_line(label=f"{name}_s00_s{a+1:0>2d}", value=float(model[VMK.mtdm][d, a]))
                for a in range(model[VMK.A]-1)
            ])

        return string


def build_linear_coupling(linear_terms, A, N, diagonal=False):
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


def build_diagonal_quadratic_coupling(quadratic_terms, A, N, diagonal=False):
    """Return a string containing the quadratic coupling constant information of a .op file."""
    if diagonal:
        return ''.join([
            make_line(
                label=f"C2_s{a+1:0>2d}s{a+1:0>2d}_v{i+1:0>2d}v{i+1:0>2d}",
                value=quadratic_terms[i, i, a]
            )
            for a, i in it.product(range(A), range(N))
            if not np.isclose(quadratic_terms[i, i, a], 0.0)
        ])
    else:
        # other option
        # return '\n'.join([
        #     ''.join([
        #         make_line(
        #             label=f"C2_s{a+1:0>2d}s{a+1:0>2d}_v{i+1:0>2d}v{i+1:0>2d}",
        #             value=quadratic_terms[i, i, a, a]
        #         )
        #         for a, i in it.product(range(A), range(N))
        #         if not np.isclose(quadratic_terms[i, i, a, a], 0.0)
        #     ]),
        #     ''.join([
        #         make_line(
        #             label=f"C2_s{a2+1:0>2d}s{a1+1:0>2d}_v{i+1:0>2d}v{i+1:0>2d}",
        #             value=quadratic_terms[i, i, a2, a1]
        #         )
        #         for a1, a2, i in it.product(range(A), range(A), range(N))
        #         if a1 != a2 and not np.isclose(quadratic_terms[i, i, a2, a1], 0.0)
        #     ]),
        # ])
        # we use this one because it matches the provided models
        return ''.join([
            make_line(
                label=f"C2_s{a2+1:0>2d}s{a1+1:0>2d}_v{i+1:0>2d}v{i+1:0>2d}",
                value=quadratic_terms[i, i, a2, a1]
            )
            for a1, a2, i in it.product(range(A), range(A), range(N))
            if not np.isclose(quadratic_terms[i, i, a2, a1], 0.0)
        ])


def build_off_diagonal_quadratic_coupling(quadratic_terms, A, N, diagonal=False):
    """Return a string containing the quadratic coupling constant information of a .op file."""
    if diagonal:
        return ''.join([
            make_line(
                label=f"C2_s{a+1:0>2d}s{a+1:0>2d}_v{j1+1:0>2d}v{j2+1:0>2d}",
                value=quadratic_terms[j1, j2, a]
            )
            for a, j1, j2 in it.product(range(A), range(N), range(N))
            if (j1 != j2) and not np.isclose(quadratic_terms[j1, j2, a], 0.0)
        ])
    else:
        return ''.join([
            make_line(
                label=f"C2_s{a2+1:0>2d}s{a1+1:0>2d}_v{j1+1:0>2d}v{j2+1:0>2d}",
                value=quadratic_terms[j1, j2, a2, a1]
            )
            for a1, a2, j1, j2 in it.product(range(A), range(A), range(N), range(N))
            if (j1 != j2) and not np.isclose(quadratic_terms[j1, j2, a2, a1], 0.0)
        ])


def build_diagonal_cubic_coupling(cubic_terms, A, N, diagonal=False):
    """Return a string containing the quadratic coupling constant information of a .op file."""
    if diagonal:
        return ''.join([
            make_line(
                label=f"C3_s{a+1:0>2d}_s{a+1:0>2d}_v{i+1:0>2d}",
                value=cubic_terms[i, i, i, a]
            )
            for a, i in it.product(range(A), range(N))
            if not np.isclose(cubic_terms[i, i, i, a], 0.0)
        ])
    else:
        # we use this one because it matches the provided models
        return ''.join([
            make_line(
                label=f"C3_s{a2+1:0>2d}_s{a1+1:0>2d}_v{i+1:0>2d}",
                value=cubic_terms[i, i, i, a2, a1]
            )
            for a1, a2, i in it.product(range(A), range(A), range(N))
            if not np.isclose(cubic_terms[i, i, i, a2, a1], 0.0)
        ])


def build_diagonal_quartic_coupling(quartic_terms, A, N, diagonal=False):
    """Return a string containing the quadratic coupling constant information of a .op file."""
    if diagonal:
        return ''.join([
            make_line(
                label=f"C4_s{a+1:0>2d}_s{a+1:0>2d}_v{i+1:0>2d}",
                value=quartic_terms[i, i, i, i, a]
            )
            for a, i in it.product(range(A), range(N))
            if not np.isclose(quartic_terms[i, i, i, i, a], 0.0)
        ])
    else:
        # we use this one because it matches the provided models
        return ''.join([
            make_line(
                label=f"C4_s{a2+1:0>2d}_s{a1+1:0>2d}_v{i+1:0>2d}",
                value=quartic_terms[i, i, i, i, a2, a1]
            )
            for a1, a2, i in it.product(range(A), range(A), range(N))
            if not np.isclose(quartic_terms[i, i, i, i, a2, a1], 0.0)
        ])


def build_parameter_section(diagonal_model, flag_diagonal=False, highest_order=1):
    """Returns a string which defines the `PARAMETER-SECTION` of an .op file"""
    start, end = "PARAMETER-SECTION", "end-parameter-section"
    header_string = "#{:^47}#\n#{:^47}#\n"

    A, N = diagonal_model[VMK.A], diagonal_model[VMK.N]

    return_list = [
        start,
        header_string.format(headers[0], '-' * 45),
        build_frequencies(diagonal_model),
        header_string.format(headers[1], '-' * 45),
        build_energies(diagonal_model, flag_diagonal),
        header_string.format(headers[2], '-' * 45),
        build_electronic_moments(diagonal_model),
        header_string.format(headers[3], '-' * 45),
        build_magnetic_moments(diagonal_model),
        header_string.format(headers[4], '-' * 45),
        build_linear_coupling(diagonal_model[VMK.G1], A, N, flag_diagonal),
    ]

    if highest_order > 1:
        return_list.append(header_string.format(headers[5], '-' * 45))
        return_list.append(build_diagonal_quadratic_coupling(diagonal_model[VMK.G2], A, N, flag_diagonal))
        return_list.append(header_string.format(headers[6], '-' * 45))
        return_list.append(build_off_diagonal_quadratic_coupling(diagonal_model[VMK.G2], A, N, flag_diagonal))

    if highest_order > 2:
        return_list.append(header_string.format(headers[7], '-' * 45))
        return_list.append(build_diagonal_cubic_coupling(diagonal_model[VMK.G3], A, N, flag_diagonal))

    if highest_order > 3:
        return_list.append(header_string.format(headers[8], '-' * 45))
        return_list.append(build_diagonal_quartic_coupling(diagonal_model[VMK.G4], A, N, flag_diagonal))

    return_list.append(end)

    return '\n'.join(return_list)

    #   return '\n'.join([
    #       start,
    #       header_string.format(headers[0], '-' * 45),
    #       build_frequencies(diagonal_model),
    #       header_string.format(headers[1], '-' * 45),
    #       build_energies(diagonal_model, flag_diagonal),
    #       header_string.format(headers[2], '-' * 45),
    #       build_electronic_moments(diagonal_model),
    #       header_string.format(headers[3], '-' * 45),
    #       build_magnetic_moments(diagonal_model),
    #       header_string.format(headers[4], '-' * 45),
    #       build_linear_coupling(diagonal_model[VMK.G1], A, N, flag_diagonal),
    #       header_string.format(headers[5], '-' * 45),
    #       build_diagonal_quadratic_coupling(diagonal_model[VMK.G2], A, N, flag_diagonal),
    #       header_string.format('Off_diagonal Quadratic Coupling Constants', '-' * 45),
    #       header_string.format(headers[6], '-' * 45),
    #       build_off_diagonal_quadratic_coupling(diagonal_model[VMK.G2], A, N, flag_diagonal),
    #       header_string.format(headers[7], '-' * 45),
    #       build_diagonal_cubic_coupling(diagonal_model[VMK.G3], A, N, flag_diagonal),
    #       header_string.format(headers[8], '-' * 45),
    #       build_diagonal_quartic_coupling(diagonal_model[VMK.G4], A, N, flag_diagonal),
    #       end,
    #   ])


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


def label_linear_coupling(linear_terms, A, N, diagonal=False):
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


def label_diagonal_quadratic_coupling(quadratic_terms, A, N, diagonal=False):
    """Return a string containing the quadratic coupling constant labelling of a .op file."""
    spacer = '|'
    if diagonal:
        return '\n'.join([
            (
                f"C2_s{a:0>2d}s{a:0>2d}_v{i:0>2d}v{i:0>2d}{spacer:>9}1 S{a:d}&{a:d}"
                f"{spacer:>4}{i+1}  q^2"
            )
            for a, i in it.product(range(1, A+1), range(1, N+1))
            if not np.isclose(quadratic_terms[i-1, i-1, a-1], 0.0)
        ]) + '\n'
    else:
        return '\n'.join([
            (
                f"C2_s{a2:0>2d}s{a1:0>2d}_v{i:0>2d}v{i:0>2d}{spacer:>9}1 S{a2:d}&{a1:d}"
                f"{spacer:>4}{i+1}  q^2"
            )
            for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
            if not np.isclose(quadratic_terms[i-1, i-1, a2-1, a1-1], 0.0)
        ]) + '\n'


def label_off_diagonal_quadratic_coupling(quadratic_terms, A, N, diagonal=False):
    """Return a string containing the quadratic coupling constant labelling of a .op file."""
    spacer = '|'
    if diagonal:
        return '\n'.join([
            (
                f"C2_s{a:0>2d}s{a:0>2d}_v{j1:0>2d}v{j2:0>2d}{spacer:>9}1 S{a:d}&{a:d}"
                f"{spacer:>4}{j1+1}  q{spacer:>6}{j2+1}  q"
            )
            for a, j1, j2 in it.product(range(1, A+1), range(1, N+1), range(1, N+1))
            if (j1 != j2) and not np.isclose(quadratic_terms[j1-1, j2-1, a-1], 0.0)
        ]) + '\n'
    else:
        return '\n'.join([
            (
                f"C2_s{a2:0>2d}s{a1:0>2d}_v{j1:0>2d}v{j2:0>2d}{spacer:>9}1 S{a2:d}&{a1:d}"
                f"{spacer:>4}{j1+1}  q{spacer:>6}{j2+1}  q"
            )
            for a1, a2, j1, j2 in it.product(range(1, A+1), range(1, A+1), range(1, N+1), range(1, N+1))
            if (j1 != j2) and not np.isclose(quadratic_terms[j1-1, j2-1, a2-1, a1-1], 0.0)
        ]) + '\n'


def label_diagonal_cubic_coupling(cubic_terms, A, N, diagonal=False):
    """Return a string containing the cubic coupling constant labelling of a .op file."""
    spacer = '|'
    if diagonal:
        return '\n'.join([
            (
                f"C3_s{a:0>2d}_s{a:0>2d}_v{i:0>2d}{spacer:>11}1 S{a:d}&{a:d}"
                f"{spacer:>4}{i+1}  q^3"
            )
            for a, i in it.product(range(1, A+1), range(1, N+1))
            if not np.isclose(cubic_terms[i-1, i-1, i-1, a-1], 0.0)
        ]) + '\n'
    else:
        return '\n'.join([
            (
                f"C3_s{a2:0>2d}_s{a1:0>2d}_v{i:0>2d}{spacer:>11}1 S{a2:d}&{a1:d}"
                f"{spacer:>4}{i+1}  q^3"
            )
            for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
            if not np.isclose(cubic_terms[i-1, i-1, i-1, a2-1, a1-1], 0.0)
        ]) + '\n'


def label_diagonal_quartic_coupling(quartic_terms, A, N, diagonal=False):
    """Return a string containing the quartic coupling constant labelling of a .op file."""
    spacer = '|'
    if diagonal:
        return '\n'.join([
            (
                f"C4_s{a:0>2d}_s{a:0>2d}_v{i:0>2d}{spacer:>11}1 S{a:d}&{a:d}"
                f"{spacer:>4}{i+1}  q^4"
            )
            for a, i in it.product(range(1, A+1), range(1, N+1))
            if not np.isclose(quartic_terms[i-1, i-1, i-1, i-1, a-1], 0.0)
        ]) + '\n'
    else:
        return '\n'.join([
            (
                f"C4_s{a2:0>2d}_s{a1:0>2d}_v{i:0>2d}{spacer:>11}1 S{a2:d}&{a1:d}"
                f"{spacer:>4}{i+1}  q^4"
            )
            for a1, a2, i in it.product(range(1, A+1), range(1, A+1), range(1, N+1))
            if not np.isclose(quartic_terms[i-1, i-1, i-1, i-1, a2-1, a1-1], 0.0)
        ]) + '\n'


def build_hamiltonian_section(model, flag_diagonal=False, highest_order=1):
    """Returns a string which defines the `HAMILTONIAN-SECTION` of an .op file"""
    start, end = "HAMILTONIAN-SECTION", "end-hamiltonian-section"
    spec = ''.join([
        ' modes   |  el  |',
        ''.join([f" v{N+1:0>2d}|" for N in range(model[VMK.N])]),
        '\n'
    ])

    A, N = model[VMK.A], model[VMK.N]

    return_list = [
        start,
        spec,
        label_momentum(N),
        label_position(N),
        label_energies(A),
        label_linear_coupling(model[VMK.G1], A, N, flag_diagonal),
        ]

    if highest_order > 1:
        return_list.append(label_diagonal_quadratic_coupling(model[VMK.G2], A, N, flag_diagonal))
        return_list.append(label_off_diagonal_quadratic_coupling(model[VMK.G2], A, N, flag_diagonal))

    if highest_order > 2:
        return_list.append(label_diagonal_cubic_coupling(model[VMK.G3], A, N, flag_diagonal))

    if highest_order > 3:
        return_list.append(label_diagonal_quartic_coupling(model[VMK.G4], A, N, flag_diagonal))

    return_list.append(end)

    return '\n'.join(return_list)

#    string = '\n'.join([
#        start,
#        spec,
#        label_momentum(N),
#        label_position(N),
#        label_energies(A),
#        label_linear_coupling(model[VMK.G1], A, N, flag_diagonal),
#        label_diagonal_quadratic_coupling(model[VMK.G2], A, N, flag_diagonal),
#        label_off_diagonal_quadratic_coupling(model[VMK.G2], A, N, flag_diagonal),
#        label_diagonal_cubic_coupling(model[VMK.G3], A, N, flag_diagonal),
#        label_diagonal_quartic_coupling(model[VMK.G4], A, N, flag_diagonal),
#        end,
#    ])
#    return string


def build_dipole_moments_section(model, flag_diagonal=False):
    """Returns a string which defines the `HAMILTONIAN-SECTION_Ex` of an .op file"""
    start, end = "HAMILTONIAN-SECTION_Ex", "end-hamiltonian-section"
    spec = ''.join([
        ' modes   |  el  |',
        ''.join([f" v{N+1:0>2d}|" for N in range(model[VMK.N])]),
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
        label_dipole_moments_energies(model[VMK.A]),
        end,
    ])
    return string


def build_magnetic_section(model, flag_diagonal=False):
    """Returns a string which defines the `HAMILTONIAN-SECTION_Mx` of an .op file"""
    start, end = "HAMILTONIAN-SECTION_Mx", "end-hamiltonian-section"
    spec = ''.join([
        ' modes   |  el  |',
        ''.join([f" v{N+1:0>2d}|" for N in range(model[VMK.N])]),
        '\n'
    ])

    def label_magnetic_moments_energies(A):
        """Return a string containing the magnetic moments energy labelling of a .op file."""
        spacer = '|'
        return '\n'.join([
            f"Mx_s00_s{a:0>2d}{spacer:>15}1 S{A:d}&{a:d}"
            for a in range(1, A)
        ]) + '\n'

    string = '\n'.join([
        start,
        spec,
        label_magnetic_moments_energies(model[VMK.A]),
        end,
    ])
    return string


def enforce_proper_keys_are_present_in_model(model, flag_diagonal=False, highest_order=1):
    """ Force all keys to be present, fill with zeros (from `template_model`) if not present """

    from .vibronic_model_io import (
        diagonal_model_zeros_template_json_dict,
        model_zeros_template_json_dict
    )

    if flag_diagonal:
        template_model = diagonal_model_zeros_template_json_dict(model[VMK.A], model[VMK.N], highest_order)
    else:
        template_model = model_zeros_template_json_dict(model[VMK.A], model[VMK.N], highest_order)

    # check to make sure it has all the keys we require
    for n, key in enumerate(VMK.coupling_list()):

        # only include keys up to highest order
        if n+1 > highest_order:
            break

        if key not in model.keys():
            model[key] = template_model[key]

    # also force the dipole moments to exist
    for key in [VMK.etdm, VMK.mtdm]:
        if key not in model.keys():
            model[key] = template_model[key]

    return  # model is dictionary, no return needed


def generate_op_file_data(model, flag_diagonal=False, highest_order=1):
    """ Returns a string formatted for a .op file using information from `model` """
    enforce_proper_keys_are_present_in_model(model, flag_diagonal, highest_order)

    return '\n\n\n'.join([
        build_op_section(),
        build_parameter_section(model, flag_diagonal, highest_order=highest_order),
        build_hamiltonian_section(model, flag_diagonal, highest_order=highest_order),
        build_dipole_moments_section(model, flag_diagonal),
        "end-operator",
    ])
