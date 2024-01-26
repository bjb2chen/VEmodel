"""vibronic_model_io.py should handle the majority of file I/O"""

# system imports
import itertools as it
import random
import shutil
import copy
import json
import os
from os.path import isfile

# third party imports
import numpy as np
from numpy import float64 as F64
from numpy import complex128 as C128

# local imports
from .log_conf import log
from .vibronic_model_keys import VibronicModelKeys as VMK
from . import model_op

# functions that begin with an underscore are meant to be used inside this module only
# they aren't meant to be called and/or imported by other code

# Force any printing of numpy arrays to be at a set precision
# to improve reading of output/debug print statements/logs
np.set_printoptions(precision=8, suppress=True)


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Functions which define the basic structure of a model/diagonal model
# ------------------------------------------------------------------------
def model_shape_dict(A, N):
    """ returns a dictionary with the same keys as the .json file whose values are tuples representing the dimensionality of the associated value in the .json file
    Takes A - number of surfaces and N - number of modes

    """
    dictionary = {
        VMK.tdm: (A, ),
        VMK.E: (A, A),
        VMK.w: (N, ),
        VMK.G1: (N, A, A),
        VMK.G2: (N, N, A, A),
        VMK.G3: (N, N, N, A, A),
        VMK.G4: (N, N, N, N, A, A),
    }

    return dictionary


def diagonal_model_shape_dict(A, N):
    """ returns a dictionary with the same keys as the .json file whose values are tuples representing the dimensionality of the associated value in the .json file
    Takes A - number of surfaces and N - number of modes

    """
    dictionary = {
        VMK.tdm: (A, ),
        VMK.E: (A, ),
        VMK.w: (N, ),
        VMK.G1: (N, A),
        VMK.G2: (N, N, A),
        VMK.G3: (N, N, N, A),
        VMK.G4: (N, N, N, N, A),
    }

    return dictionary


def model_zeros_template_json_dict(A, N, highest_order=4):
    """ returns a dictionary that is a valid model, where all values (other than states and modes) are set to 0
    """
    shape = model_shape_dict(A, N)
    dictionary = {
        VMK.N: N,
        VMK.A: A,
        VMK.E: np.zeros(shape[VMK.E], dtype=F64),
        VMK.w: np.zeros(shape[VMK.w], dtype=F64),
        VMK.tdm: np.zeros(shape[VMK.tdm], dtype=C128),
    }

    for idx, key in enumerate(VMK.coupling_list()):
        if idx + 1 <= highest_order:
            dictionary.update({key: np.zeros(shape[key], dtype=F64)})

    return dictionary


def diagonal_model_zeros_template_json_dict(A, N, highest_order=4):
    """ returns a dictionary that is a valid diagonal model, where all values (other than states and modes) are set to 0
    """
    shape = diagonal_model_shape_dict(A, N)
    dictionary = {
        VMK.N: N,
        VMK.A: A,
        VMK.E: np.zeros(shape[VMK.E], dtype=F64),
        VMK.w: np.zeros(shape[VMK.w], dtype=F64),
        VMK.tdm: np.zeros(shape[VMK.tdm], dtype=C128),
    }

    for idx, key in enumerate(VMK.coupling_list()):
        if idx + 1 <= highest_order:
            dictionary.update({key: np.zeros(shape[key], dtype=F64)})

    return dictionary


def verify_model_parameters(kwargs):
    """make sure the provided model parameters follow the file conventions"""
    assert VMK.N in kwargs, "need the number of modes"
    assert VMK.A in kwargs, "need the number of surfaces"

    A, N = _extract_dimensions_from_dictionary(kwargs)
    shape_dict = model_shape_dict(A, N)

    for key, value in kwargs.items():
        if key in shape_dict:
            assert kwargs[key].shape == shape_dict[key], f"{key} have incorrect shape"
        else:
            log.debug(f"Found key {key} which is not present in the default dictionary")

    return


def verify_diagonal_model_parameters(kwargs):
    """make sure the provided sample parameters follow the file conventions"""
    assert VMK.N in kwargs, "need the number of modes"
    assert VMK.A in kwargs, "need the number of surfaces"

    A, N = _extract_dimensions_from_dictionary(kwargs)
    shape_dict = diagonal_model_shape_dict(A, N)

    for key, value in kwargs.items():
        if key in shape_dict:
            assert kwargs[key].shape == shape_dict[key], f"{key} have incorrect shape"
        else:
            log.debug(f"Found key {key} which is not present in the default dictionary")

    return


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Functions which return specific parameters of a model
# ------------------------------------------------------------------------
def _extract_dimensions_from_dictionary(dictionary):
    """x"""
    N = int(dictionary[VMK.N])
    A = int(dictionary[VMK.A])
    return A, N


def _extract_dimensions_from_file(path):
    """x"""
    assert isfile(path), f"invalid path:\n{path:s}"
    with open(path, mode='r', encoding='UTF8') as file:
        input_dictionary = json.loads(file.read())
    VMK.change_dictionary_keys_from_strings_to_enum_members(input_dictionary)
    return _extract_dimensions_from_dictionary(input_dictionary)


def extract_dimensions_of_model(path=None):
    """returns A, N in that order
    A is the number_of_surfaces and N is the number_of_modes
    for coupling_model.json files by using an absolute path to the file
    """
    assert path is not None, "no arguments provided"
    if type(path) is dict:
        return _extract_dimensions_from_dictionary(path)
    return _extract_dimensions_from_file(path)


def extract_dimensions_of_diagonal_model(path=None):
    """returns A, N in that order
    A is the number_of_surfaces and N is the number_of_modes
    for sampling_model.json files by using an absolute path to the file
    """
    assert path is not None, "no arguments provided"
    if type(path) is dict:
        return _extract_dimensions_from_dictionary(path)
    return _extract_dimensions_from_file(path)


def extract_maximum_order_of_model(model):
    """Returns integer `max_order`, the highest order component with non zero entries"""

    for key in reversed(VMK.key_list()):
        if (key in model) and not np.allclose(model[key], 0):
            # print('index', VMK.key_list().index(key))
            return VMK.key_list().index(key)
    else:
        return 0


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Boolean functions, for checking the truthfulness of some property
# ------------------------------------------------------------------------
def _array_is_symmetric_in_A(array):
    """ Boolean function that returns true if the provided numpy array is symmetric in the surface dimension
    where the surface dimensions (A) are by convention the last two dimensions
    this function assumes that the array is properly formatted
    """

    new_dims = list(range(array.ndim))
    # swap the last two dimensions, which by convention are the surface dimensions
    new_dims[-1], new_dims[-2] = new_dims[-2], new_dims[-1]

    return np.allclose(array, array.transpose(new_dims))


def _same_model(d1, d2):
    """ returns True if all parameters of the two dictionaries have the same dimensions
    and the same floating point numbers up to standard precision comparison
    raises an assertion error if either dictionary has incorrect parameters"""

    verify_model_parameters(d1)
    verify_model_parameters(d2)

    A1, N1 = _extract_dimensions_from_dictionary(d1)
    A2, N2 = _extract_dimensions_from_dictionary(d1)

    if A1 != A2 or N1 != N2:
        return False

    new_d1 = model_zeros_template_json_dict(A1, N1)
    new_d1.update(d1)

    new_d2 = model_zeros_template_json_dict(A2, N2)
    new_d2.update(d2)

    for key in new_d1.keys():
        if not np.allclose(new_d1[key], new_d2[key]):
            log.debug(f"These models differ for key {key}\nd1: {new_d1[key]}\nd2: {new_d2[key]}")
            return False
        # elif (key not in d2 and np.count_nonzero(d1[key]) > 0) or \
        #      (key not in d1 and np.count_nonzero(d2[key]) > 0):
        #     print(f"These models differ for key {key}\nd1: {d1[key]}\nd2: {d2[key]}")
        #     return False

    else:
        return True

    raise Exception("This line of code should not be reached!?")


def _same_diagonal_model(d1, d2):
    """ returns True if all parameters of the two dictionaries have the same dimensions
    and the same floating point numbers up to standard precision comparison
    raises an assertion error if either dictionary has incorrect parameters"""

    verify_diagonal_model_parameters(d1)
    verify_diagonal_model_parameters(d2)

    A1, N1 = _extract_dimensions_from_dictionary(d1)
    A2, N2 = _extract_dimensions_from_dictionary(d1)

    if A1 != A2 or N1 != N2:
        return False

    new_d1 = diagonal_model_zeros_template_json_dict(A1, N1)
    new_d1.update(d1)

    new_d2 = diagonal_model_zeros_template_json_dict(A2, N2)
    new_d2.update(d2)

    for key in new_d1.keys():
        if not np.allclose(new_d1[key], new_d2[key]):
            print(f"These diagonal models differ for key {key}\nd1: {new_d1[key]}\nd2: {new_d2[key]}")
            return False
    else:
        return True

    raise Exception("This line of code should not be reached!?")


def model_parameters_are_symmetric_in_surfaces(kwargs):
    """ Boolean function that returns true if the provided model's arrays are all symmetric in the surface dimension
    where the surface dimensions (A) are by convention the last two dimensions
    this function assumes that the arrays is properly formatted
    """
    verify_model_parameters(kwargs)

    for key in kwargs.keys():
        if isinstance(kwargs[key], np.ndarray) and kwargs[key].ndim >= 2:
            if not _array_is_symmetric_in_A(kwargs[key]):
                log.debug(f"{key} not symmetric in surfaces")
                return False

    return True


def model_parameters_are_symmetric_in_modes(kwargs):
    """ Boolean function that returns true if the provided model's quadratic and quartic arrays are symmetric in their mode dimensions
    this is a bit trickier than the surface dimension
    this function assumes that the arrays is properly formatted
    """
    verify_model_parameters(kwargs)

    # do the quadratic case
    key = VMK.G2
    if key in kwargs.keys():
        new_dims = list(range(kwargs[key].ndim))
        new_dims[0], new_dims[1] = new_dims[1], new_dims[0]
        if not np.allclose(kwargs[key], kwargs[key].transpose(new_dims)):
            return False

    # I don't think there is a general approach for the cubic case

    # do the quartic case
    key = VMK.G4
    if key in kwargs.keys():
        new_dims = list(range(kwargs[key].ndim))
        new_dims[0], new_dims[1] = new_dims[1], new_dims[0]
        new_dims[2], new_dims[3] = new_dims[3], new_dims[2]
        if not np.allclose(kwargs[key], kwargs[key].transpose(new_dims)):
            return False

    return True


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Functions which modify a model in some way
# ------------------------------------------------------------------------
def _make_model_diagonal(model):
    """Returns a diagonal model filled with the diagonal components of the argument `model`."""

    new_model = {}

    for key, value in model.items():
        if hasattr(value, 'shape') and len(value.shape) >= 2:
            ndims = len(value.shape)
            new_model[key] = np.diagonal(model[key], axis1=ndims - 2, axis2=ndims - 1).copy()
        else:
            new_model[key] = model[key]

    return new_model


def model_add_surface(old_model, **kwargs):
    """Modifies the provided model using kwargs and returns a copy."""

    verify_model_parameters(old_model)
    old_A = old_model[VMK.A]

    A = old_A + 1
    N = old_model[VMK.N]

    new_model = model_zeros_template_json_dict(A, N, highest_order=0)

    # energy
    for a in range(old_A):
        new_model[VMK.E][a, a] = old_model[VMK.E][a, a]

    # add energy of new surface from `kwargs`
    new_model[VMK.E][-1, -1] = kwargs[VMK.E]

    # transition dipole moment
    for a in range(old_A):
        new_model[VMK.tdm][a] = old_model[VMK.tdm][a]

    # add transition moment of new surface from `kwargs`
    new_model[VMK.tdm][-1] = kwargs[VMK.tdm]

    # frequencies don't change
    new_model[VMK.w] = old_model[VMK.w]

    # linear terms
    if VMK.G1 in kwargs:
        # not well defined how to add linear terms
        raise Exception(f"Functionality not implemented for key {VMK.G1}")
        for j in range(N):
            new_model[VMK.G1][j, :-1:] = old_model[VMK.G1]
            new_model[VMK.G1][j, -1] = kwargs[VMK.G1][j]

    # 2nd and higher order terms
    for key in [VMK.G1, VMK.G2, VMK.G3, VMK.G4]:
        if key in kwargs:
            # it is not clear how these terms would be provided and how they would be modified
            raise Exception(f"Functionality not implemented for key {key}")

    return new_model


def model_remove_ground_state(old_model, **kwargs):
    """Modifies the provided model using kwargs and returns a copy.
    Only works if provided surface is last
    """
    verify_model_parameters(old_model)
    old_A = old_model[VMK.A]

    A = old_A - 1
    N = old_model[VMK.N]

    order = extract_maximum_order_of_model(old_model)
    new_model = model_zeros_template_json_dict(A, N, highest_order=order)

    # frequencies don't change
    new_model[VMK.w] = old_model[VMK.w]

    # energy
    for a in range(A):
        new_model[VMK.tdm][a] = old_model[VMK.tdm][a]
        for b in range(A):
            new_model[VMK.E][a, b] = old_model[VMK.E][a, b]
            for key in [VMK.G1, VMK.G2, VMK.G3, VMK.G4]:
                if key in old_model:
                    new_model[key][..., a, b] = old_model[key][..., a, b]

    return new_model


def diagonal_model_add_surface(old_model, **kwargs):
    """Modifies the provided model using kwargs and returns a copy."""
    verify_diagonal_model_parameters(old_model)
    old_A = old_model[VMK.A]

    A = old_A + 1
    N = old_model[VMK.N]

    new_model = diagonal_model_zeros_template_json_dict(A, N, highest_order=1)

    # energy
    for a in range(old_A):
        new_model[VMK.E][a] = old_model[VMK.E][a]

    # add energy of new surface from `kwargs`
    new_model[VMK.E][-1] = kwargs[VMK.E]

    # transition dipole moment
    for a in range(old_A):
        new_model[VMK.tdm][a] = old_model[VMK.tdm][a]

    # add transition moment of new surface from `kwargs`
    new_model[VMK.tdm][-1] = kwargs[VMK.tdm]

    # frequencies don't change
    new_model[VMK.w] = old_model[VMK.w]

    # linear terms
    if VMK.G1 in kwargs:
        for j in range(N):
            new_model[VMK.G1][j, :-1:] = old_model[VMK.G1]
            new_model[VMK.G1][j, -1] = kwargs[VMK.G1][j]

    # 2nd and higher order terms
    for key in [VMK.G1, VMK.G2, VMK.G3, VMK.G4]:
        if key in kwargs:
            # it is not clear how these terms would be provided and how they would be modified
            raise Exception(f"Functionality not implemented for key {key}")

    return new_model


def diagonal_model_remove_surface(old_model, A_idxs):
    """Removes the surface with index A_idx, from the provided model and returns a copy."""
    verify_diagonal_model_parameters(old_model)

    # if we are given a single integer convert to a list
    if type(A_idxs) is int:
        A_idxs = [A_idxs]

    A_idxs = np.array(A_idxs)

    old_A = old_model[VMK.A]
    assert (old_A > 1) and (old_A > len(A_idxs)), \
        f"Cannot remove a surface if they are no surfaces left"
    assert max(A_idxs) < old_A, \
        f"Cannot remove the {max(A_idxs)} surface, only {old_A} surfaces exist"
    assert len(A_idxs) < old_A, \
        f"Cannot remove {len(A_idxs)} surfaces, only {old_A} surfaces exist"

    A = old_A - len(A_idxs)
    N = old_model[VMK.N]
    new_model = diagonal_model_zeros_template_json_dict(A, N, highest_order=1)

    # energy
    new_model[VMK.E] = np.delete(old_model[VMK.E], A_idxs)

    # transition moment
    new_model[VMK.tdm] = np.delete(old_model[VMK.tdm], A_idxs)

    # frequencies
    new_model[VMK.w] = old_model[VMK.w]

    # linear terms
    if VMK.G1 in old_model:
        for j in range(N):
            new_model[VMK.G1][j, :] = np.delete(old_model[VMK.G1][j, :], A_idxs)

    # print("old lins\n", old_model[VMK.G1])
    # print("new lins\n", new_model[VMK.G1])

    return new_model


def fill_offdiagonalmodes_of_model_with_zeros(model):
    """ takes a dictionary who must have values of dimensionality (..., A, A)
    and set the off-diagonal (normal modes) elements to zero

    Of course we restrict the usage of this function to models with order of 2
    as the 'off diagonal' modes for order 1 or > 3 is not well defined.
    """
    verify_model_parameters(model)
    order = extract_maximum_order_of_model(model)
    assert order <= 2, f"Ill defined operation for cubic or higher order coupling terms\nOrder was ({order:d})"
    if order <= 1:
        # if only linear terms then there are no off-diagonal mode terms
        return
    else:
        # we only treat the quadratic terms at the moment
        for i,j in it.permutations(range(model[VMK.N]), 2):
            model[VMK.G2][i, j, ...] = 0.0

    return


def fill_offdiagonalsurfaces_of_model_with_zeros(model):
    """ takes a dictionary who must have values of dimensionality (..., A, A)
    and set the off-diagonal (surface) elements to zero
    """
    verify_model_parameters(model)
    for key, value in model.items():
        if hasattr(value, 'shape') and len(value.shape) >= 2:
            for a, b in it.permutations(range(model[VMK.A]), 2):
                model[key][..., a, b] = 0.0
    return


def remove_coupling_from_model(path_source, path_destination):
    """reads in a model from path_source whose values can have dimensionality (..., A, A)
    creates a new model whose values have dimensionality (..., A) from the diagonal of the A dimension of the input model
    saves the new model to the provided path_destination
    """
    model = load_model_from_JSON(path_source)
    new_model = _make_model_diagonal(model)
    save_diagonal_model_to_JSON(path_destination, new_model)
    return


def remove_higher_order_terms(model, highest_order=4):
    for idx, key in enumerate(VMK.key_list()):
        if key in model and idx > highest_order:
            del model[key]
    return


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Functions which generate models with random values
# ------------------------------------------------------------------------
def _generate_linear_terms(linear_terms, shape, displacement, Modes):
    """ generate linear terms that are 'reasonable' """
    for i in Modes:
        upTri = np.random.uniform(-displacement[i], displacement[i], shape[VMK.E])
        # force the linear terms to be symmetric
        linear_terms[i, ...] = np.tril(upTri) + np.tril(upTri, k=-1).T
    return


def _generate_quadratic_terms(quadratic_terms, shape, displacement, Modes):
    """ generate quadratic terms that are 'reasonable' """
    for i, j in it.product(Modes, repeat=2):
        upTri = np.random.uniform(-displacement[i, j], displacement[i, j], shape[VMK.E])
        # force the quadratic terms to be symmetric
        quadratic_terms[i, j, ...] = np.tril(upTri) + np.tril(upTri, k=-1).T
        quadratic_terms[j, i, ...] = np.tril(upTri) + np.tril(upTri, k=-1).T
    return


def generate_vibronic_model_data(input_parameters=None):
    """redo this one but otherwise its fine returns e,w,l,q filled with appropriate values"""

    # default values
    paramDict = {
        'frequency_range': [0.02, 0.04],
        'energy_range': [0.0, 2.0],
        'quadratic_scaling': 0.08,
        'linear_scaling': 0.04,
        'diagonal': False,
        'numStates': 2,
        'numModes': 3,
    }

    # overwrite default values with any provided arguments
    if input_parameters is not None:
        paramDict.update(input_parameters)

    # readability
    minE, maxE = paramDict['energy_range']
    minFreq, maxFreq = paramDict['frequency_range']

    # ranges for convenience
    numModes = paramDict['numModes']
    numStates = paramDict['numStates']
    Modes = range(numModes)

    # generate the array dimensions
    shape = model_shape_dict(numStates, numModes)

    # assume we are building a coupled model
    model = model_zeros_template_json_dict(numStates, numModes)

    # generate frequencies
    model[VMK.w] = np.linspace(minFreq, maxFreq, num=numModes, endpoint=True, dtype=F64)

    # generate energy
    model[VMK.E] = np.random.uniform(minE, maxE, shape[VMK.E])
    # force the energy to be symmetric
    model[VMK.E] = np.tril(model[VMK.E]) + np.tril(model[VMK.E], k=-1).T

    # I'm not sure what appropriate parameters would be for this
    # so for now we'll just use 0.1 eV's for all the values
    model[VMK.tdm] = np.full(shape=shape[VMK.tdm], fill_value=0.1, dtype=C128)

    # calculate the linear displacement
    l_shift = paramDict['linear_scaling'] / model[VMK.w]
    _generate_linear_terms(model[VMK.G1], shape, l_shift, Modes)

    # TODO - no quadratic terms for the moment - turn back on after further testing
    # calculate the quadratic displacement
    # q_shift = np.sqrt(np.outer(frequencies, frequencies)) / paramDict['quadratic_scaling']
    # _generate_quadratic_terms(q_terms, shape, q_shift, Modes)

    # if we are building a harmonic model then zero out all off-diagonal entries
    if paramDict['diagonal']:
        d_model = diagonal_model_zeros_template_json_dict(numStates, numModes)

        d_model[VMK.E] = np.diag(model[VMK.E])
        d_model[VMK.tdm] = model[VMK.tdm]
        d_model[VMK.w] = model[VMK.w]

        for i in Modes:
            d_model[VMK.G1][i, ...] = np.diag(model[VMK.G1][i, ...])
        for i, j in it.product(Modes, repeat=2):
            d_model[VMK.G2][i, j, ...] = np.diag(model[VMK.G2][i, j, ...])

        return d_model

    assert model_parameters_are_symmetric_in_surfaces(model)
    assert model_parameters_are_symmetric_in_modes(model)

    return model


def create_random_model():
    """ returns a dictionary that is a valid model
    """
    d = {'numStates': random.randint(2, 10),
         'numModes': random.randint(2, 20),
         'quadratic_scaling': random.uniform(0.04, 0.12),
         'linear_scaling': random.uniform(0.02, 0.06),
         'diagonal': False,
         }
    return generate_vibronic_model_data(d)


def create_random_diagonal_model():
    """ returns a dictionary that is a valid diagonal model
    """
    d = {'numStates': random.randint(2, 10),
         'numModes': random.randint(2, 20),
         'quadratic_scaling': random.uniform(0.04, 0.12),
         'linear_scaling': random.uniform(0.02, 0.06),
         'diagonal': True,
         }
    return generate_vibronic_model_data(d)


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Functions which handle I/O with files ending in `.op`
# ------------------------------------------------------------------------
def read_model_op_file(
    path_file_op,
    surface_symmetrize=True,
    double_quadratic=True,
    symmetrize_quadratic=False,
):
    """ wrapper function to maintain functionality - possible remove in the future"""
    return model_op.read_model_op_file(
        path_file_op,
        surface_symmetrize=surface_symmetrize,
        double_quadratic=double_quadratic,
        symmetrize_quadratic=symmetrize_quadratic
    )


def read_raw_model_op_file(path_file_op):
    """ wrapper function to maintain functionality - possible remove in the future"""
    return model_op.read_model_op_file(
        path_file_op,
        surface_symmetrize=False,
        double_quadratic=False,
        symmetrize_quadratic=False
    )


def extract_excited_state_model_op(
    path_file_op,
    surface_symmetrize=True,
    symmetrize_quadratic=False,
    double_quadratic=False,
):
    """ wrapper function to maintain functionality - possible remove in the future"""
    model = model_op.read_model_op_file(
        path_file_op,
        surface_symmetrize=surface_symmetrize,
        double_quadratic=double_quadratic,
        symmetrize_quadratic=symmetrize_quadratic
    )
    new_model = model_remove_ground_state(model)
    return new_model


def write_diagonal_model_op_file(path_file_op, diagonal_model, half_quadratic=True):
    """Takes a diagonal model and writes it to a molecule_vibron.op files"""
    if half_quadratic:
        diagonal_model[VMK.G2] /= 2.0

    verify_diagonal_model_parameters(diagonal_model)
    formatted_data = model_op.generate_op_file_data(diagonal_model, flag_diagonal=True)
    with open(path_file_op, 'w') as source_file:
        source_file.write(formatted_data)


def write_raw_model_op_file(path_file_op, model):
    """Takes a model and writes it to a molecule_vibron.op files"""
    verify_model_parameters(model)
    formatted_data = model_op.generate_op_file_data(model)
    with open(path_file_op, 'w') as source_file:
        source_file.write(formatted_data)
    return


def write_model_op_file(path_file_op, model, half_quadratic=True):
    """Takes a model and writes it to a molecule_vibron.op files"""
    if half_quadratic:
        model[VMK.G2] /= 2.0

    write_raw_model_op_file(path_file_op, model)
    return


def create_coupling_from_op_file(dest_path, path_file_op):
    """assumes that the path_file_op is in electronic_structure"""
    model_dict = read_model_op_file(path_file_op)
    save_model_to_JSON(dest_path, model_dict)
    log.debug("Created \n{:s}\nfrom\n{:s}\n".format(dest_path, path_file_op))
    return dest_path


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Functions which handle I/O with files ending in `.json`
# ------------------------------------------------------------------------
def _save_to_JSON(path, dictionary):
    dict_copy = copy.deepcopy(dictionary)
    VMK.change_dictionary_keys_from_enum_members_to_strings(dict_copy)
    """ converts each numpy array to a list so that json can serialize them properly"""
    for key, value in list(dict_copy.items()):
        if isinstance(value, (np.ndarray, np.generic)):
            if np.count_nonzero(value) > 0:
                if key == VMK.transition_dipole_moments.value:
                    dict_copy[key] = [str(n) for n in value]
                else:
                    dict_copy[key] = value.tolist()
            else:
                del dict_copy[key]
        else:
            log.debug(f"Value {value} with Key {key} does not appear to be an ndarray")

    with open(path, mode='w', encoding='UTF8') as target_file:
        target_file.write(json.dumps(dict_copy))

    return


def save_model_to_JSON(path, dictionary):
    """ wrapper for _save_to_JSON
    calls verify_model_parameters() before calling _save_to_JSON()
    """
    verify_model_parameters(dictionary)
    log.debug(f"Saving model to {path:s}")
    _save_to_JSON(path, dictionary)
    return


def save_diagonal_model_to_JSON(path, dictionary):
    """ wrapper for _save_to_JSON
    calls verify_sample_parameters() before calling _save_to_JSON()
    """
    verify_diagonal_model_parameters(dictionary)
    log.debug(f"Saving sample to {path:s}")
    _save_to_JSON(path, dictionary)
    return


def _load_inplace_from_JSON(path, dictionary):
    """overwrites all provided values in place with the values stored in the .json file located at path"""

    with open(path, mode='r', encoding='UTF8') as file:
        input_dictionary = json.loads(file.read())

    VMK.change_dictionary_keys_from_strings_to_enum_members(input_dictionary)

    for key, value in dictionary.items():
        if isinstance(value, (np.ndarray, np.generic)):
            # this is a safer way of forcing the input arrays that have no corresponding key in the input_dictionary to have zero values
            # although this might not be necessary, it is a safer alternative at the moment
            if key not in input_dictionary:
                if key is VMK.tdm:
                    # the transition dipole moment is complex
                    # the complex numbers are stored as strings in the JSON file
                    dictionary[key].fill(complex(0.0))
                else:
                    # the rest of the model is doubles
                    dictionary[key].fill(0.0)
            elif key is VMK.tdm:
                tdm_list = list(map(complex, input_dictionary[key]))
                dictionary[key][:] = np.array(tdm_list, dtype=C128)
            else:
                dictionary[key][:] = np.array(input_dictionary[key], dtype=F64)
    return


def _load_from_JSON(path):
    """returns a dictionary filled with the values stored in the .json file located at path"""

    with open(path, mode='r', encoding='UTF8') as file:
        input_dictionary = json.loads(file.read())

    VMK.change_dictionary_keys_from_strings_to_enum_members(input_dictionary)

    for key, value in input_dictionary.items():
        if isinstance(value, list):
            if key is VMK.tdm:
                # the transition dipole moment is complex
                # the complex numbers are stored as strings in the JSON file
                value = list(map(complex, value))
                input_dictionary[key] = np.array(value, dtype=C128)
            else:
                # the rest of the model is doubles
                input_dictionary[key] = np.array(value, dtype=F64)

    # special case to always create an array of energies that are 0.0 if not provided in the .json file
    if VMK.E not in input_dictionary:
        A, N = _extract_dimensions_from_dictionary(input_dictionary)
        shape = model_shape_dict(A, N)
        input_dictionary[VMK.E] = np.zeros(shape[VMK.E], dtype=F64)

    # TODO - design decision about which arrays to fill with zeros by default?

    return input_dictionary


def load_model_from_JSON(path, dictionary=None):
    """
    if kwargs is not provided then returns a dictionary filled with the values stored in the .json file located at path

    if kwargs is provided then all values are overwritten (in place) with the values stored in the .json file located at path
    """
    log.debug(f"Loading model from {path:s}")

    # no arrays were provided so return newly created arrays after filling them with the appropriate values
    if not bool(dictionary):
        new_model_dict = _load_from_JSON(path)

        # TODO - we might want to make sure that none of the values in the dictionary have all zero values or are None

        verify_model_parameters(new_model_dict)
        return new_model_dict

    # arrays were provided so fill them with the appropriate values
    else:
        verify_model_parameters(dictionary)
        _load_inplace_from_JSON(path, dictionary)
        # check twice? might as well be cautious for the moment until test cases are written
        verify_model_parameters(dictionary)

    return


def load_diagonal_model_from_JSON(path, dictionary=None):
    """
    if kwargs is not provided then returns a dictionary filled with the values stored in the .json file located at path

    if kwargs is provided then all values are overwritten (in place) with the values stored in the .json file located at path
    """
    log.debug(f"Loading rho model (sampling model) from {path:s}")

    # no arrays were provided so return newly created arrays after filling them with the appropriate values
    if not bool(dictionary):
        new_model_dict = _load_from_JSON(path)

        # TODO - we might want to make sure that none of the values in the dictionary have all zero values or are None

        verify_diagonal_model_parameters(new_model_dict)
        return new_model_dict

    # arrays were provided so fill them with the appropriate values
    else:
        verify_diagonal_model_parameters(dictionary)
        _load_inplace_from_JSON(path, dictionary)
        # check twice? might as well be cautious for the moment until test cases are written
        verify_diagonal_model_parameters(dictionary)
    return


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# Assorted functions which don't fit in the other categories
# ------------------------------------------------------------------------
def create_harmonic_copy(path_source, path_dest):
    """wrapper function to create a harmonic model from a coupled model"""
    remove_coupling_from_model(path_source, path_dest)
    s = "Created harmonic model {:s} by removing coupling from {:s}"
    log.debug(s.format(path_dest, path_source))
    return path_dest


def create_basic_diagonal_model(path_source, path_dest):
    """wrapper function to make the simplest diagonal(sampling) model
    assumes that `path_source` points to a harmonic model
    """

    if os.path.isfile(path_dest):
        s = "Model {:s} already exists!"
        log.debug(s.format(path_dest))

    shutil.copyfile(path_source, path_dest)

    model = load_diagonal_model_from_JSON(path_dest)

    # deletes all terms of order 3 or higher
    remove_higher_order_terms(model, highest_order=2)
    # for key in VMK.key_list()[3:]:
    #     del model[key]

    save_diagonal_model_to_JSON(path_dest, model)

    s = "Created sampling model {:s} by copying {:s} and removing cubic and quartic terms"
    log.debug(s.format(path_dest, path_source))

    return path_dest


def create_deepcopy(dictionary):
    """Wrapper function for `copy.deepcopy`

    Convience function so that other code importing `vibronic_model_io` doesn't have to import copy.
    `deepcopy` is needed to create a copy of the key,value pairs of a dictionary
    """
    return copy.deepcopy(dictionary)


def create_random_orthonormal_matrix(A):
    """returns a orthonormal matrix, just a wrapper for scipy.stats.ortho_group.rvs()"""
    from scipy.stats import ortho_group
    return ortho_group.rvs(A)


def print_model(model, highest_order=None):
    """Prints all arguments of the `model` up to `highest_order` """
    print(
        f"{VMK.A.value:<20}{model[VMK.A]}",
        f"{VMK.N.value:<20}{model[VMK.N]}\n",
        sep='\n'
    )

    for key in [VMK.w, VMK.tdm, VMK.E]:
        print(f"{key.value}  {model[key].shape}\n{model[key]}\n")

    if highest_order is None:
        highest_order = extract_maximum_order_of_model(model)

    for idx, key in enumerate(VMK.coupling_list()):
        if idx + 1 <= highest_order:
            print(f"{key.value}  {model[key].shape}\n{model[key]}\n")
    return


# ------------------------------------------------------------------------
def main():
    """ currently does nothing """
    return


if (__name__ == "__main__"):
    """ x """
    main()
