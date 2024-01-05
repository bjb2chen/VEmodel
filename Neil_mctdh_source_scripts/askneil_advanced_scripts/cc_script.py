# system imports
import sys
import os
from os.path import join
import cProfile
# import pstats
# import io
import itertools as it

# third party imports
import numpy as np

# local imports
sys.path.insert(0, os.path.abspath("/home/ngraymon/public/songhao/t_amplitudes_project/t-amplitudes"))  #TEMPORARY
import project
from project.vibronic import vIO, VMK


def swap_coupling_coefficient_axes(model, coeff_order):
    """
    Currently the CC integration code expects the coefficients to have the surface dimensions first.
    When they are read in from the .op file they are the last dimensions.
    Therefore we need to shift their position.
    We do this by shifting the mode dimensions around the surface dimensions.
    """

    if coeff_order == 0:
        return  # no need to change order if their are no coupling coefficients

    print(coeff_order)
    index = VMK.key_list()[coeff_order]
    source_list = [i for i in range(coeff_order)]
    destination_list = [i for i in range(-coeff_order, 0)]
    print(source_list)
    print(destination_list)

    model[index] = np.moveaxis(model[index], source_list, destination_list)
    return


def prepare_model_for_cc_integration(model, max_order=2):
    """Removes extra parameters from .op file and reshapes the coupling coefficient tensors."""

    # remove any extra arguments we don't need (high ordered terms)
    # we want to print out if we remove any arguments while testing
    for index, key in enumerate(VMK.key_list()):
        if index > max_order and key in model:
            del model[key]
            print(f"Removed {key.name:s}")

    # we are only handling the linear and quadratic terms at the moment
    for index in [1, 2]:
        key = VMK.key_list()[index]
        if max_order >= index:
            if key not in model:
                A, N = vIO.extract_dimensions_of_model(model)
                model[key] = np.zeros(vIO.model_shape_dict(A, N)[key], dtype=float)
            swap_coupling_coefficient_axes(model, coeff_order=index)
    return


def simple_cc_run(max_order, t_final, nof_steps):
    """ temporary function while developing code """
    root_directory = os.getcwd()
    # os.walk returns a list, each element a 3-tuple, the first element is always for the root directory
    dirpaths, dirnames, filenames = next(os.walk(root_directory))
    for name in filenames:
        if ".op" in name:
            file_name = name.split('.op')[0]
            path_op = join(root_directory, name)
            break
    else:
        raise Exception("No .op file found")

    model = vIO.extract_excited_state_model_op(path_op)
    vIO.save_model_to_JSON(join(root_directory, f"{file_name}.json"), model)
    prepare_model_for_cc_integration(model, max_order=max_order)

    # ---------------------------------------------------------------------
    # debug printing
    print(f"t_final:{t_final}, nof_steps:{nof_steps}")
    print(root_directory)

    np.set_printoptions(linewidth=200, precision=16)
    vIO.print_model(model, highest_order=max_order)
    # ---------------------------------------------------------------------
    project.test_vibronic_new.generate_acf_to_file(
        model, file_name, root_directory, max_order, t_final=t_final, nof_steps=nof_steps
    )
    return


def extract_parameters(root):
    """Get the integration parameters from the folder name, and the .inp file"""

    directory_name = os.path.basename(root)
    cc_t_final = float(directory_name.split("_tf")[1])

    FC_flag = True if "_FC_" in directory_name else False
    for key,value in {"constant":0, "linear":1, "quadratic":2}.items():
        if key in directory_name:
            coupling_order = value
            break
    else:
        raise Exception("Couldn't find any keys?")

    return cc_t_final, FC_flag, coupling_order


def manage_profile_data(root):
    """Store the last 3 runs, getting rid of the oldest one each time we do another calculation.
    Copy each result to and older one, always store the latest results in cProfile_1.
    """
    spec = "cProfile_{:d}"
    dir_list = os.scandir(root)
    if os.path.exists(join(root, spec.format(3))):
        os.remove(join(root, spec.format(3)))
    for num in [2,1]:
        filename = spec.format(num)
        if filename in dir_list:
            os.system(f"mv ./{filename:s} ./{spec.format(num+1):s}")

    filename = spec.format(1)
    return filename


if (__name__ == '__main__'):

    root = os.getcwd()
    t_final, FC_flag, coupling_order = extract_parameters(root)
    nof_steps = int(t_final * 10)
    n = int(np.floor(np.log10(nof_steps)))
    r = nof_steps / pow(10, n)

    profiling = False

    if not profiling:
        simple_cc_run(coupling_order, t_final, nof_steps)
    else:
        filename = manage_profile_data(root)
        cProfile.runctx(
            'simple_cc_run(coupling_order, t_final, nof_steps)',
            globals(),
            locals(),
            filename
        )
