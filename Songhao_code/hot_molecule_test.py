# system imports
import io
import time
import os
from os.path import abspath, join, dirname, basename
import sys
import cProfile
import pstats

# third party import
import numpy as np
import json


# import the path to the package
project_dir = abspath(join(dirname(__file__), '/home/bsonghao/hot-molecule'))
sys.path.insert(0, project_dir)
inputdir = '/home/bjb2chen/VECC/vibronic_models/'
# inputdir = '/Users/pauliebao/time-dependent-vibrational-electronic-coupled-cluster-theory-for-non-adiabatic-nuclear-dynamics/original_data/vibronic_models/'
# inputdir = '/Users/pauliebao/hot-molecule/data/vibronic_models/test_models/model_on_VECC_paper/'
outputdir =  '/home/bjb2chen/VECC/hot-molecule/'

# local import
import project
from project.vibronic_model_Hamiltonian import vibronic_model_hamiltonian
from project.vibronic import vIO, VMK, vIO_wrapper

order_dict = {
    0: "constant",
    1: "linear",
    2: "quadratic",
    3: "cubic",
    4: "quartic",
}


def process_data(filename):
    """ temporary formatted printing of profiling data """
    try:
        s = io.StringIO()
        p = pstats.Stats(filename, stream = s)
        # p.strip_dirs().sort_stats("tottime").print_stats(2)
        p.strip_dirs().sort_stats("tottime").print_stats(10)
        p.strip_dirs().sort_stats("cumulative").print_stats(20)
        p.strip_dirs().sort_stats("cumulative").print_stats('contract', 100)
       # p.strip_dirs().sort_stats("cumulative").print_callees('contract_expression')

        with open(outputdir+filename+".txt", "w+") as f:
            f.write(s.getvalue())
    except Exception:
        print("cProfile data is not stored properly!")

def main():
    """main function that run TNOE simulation"""
    # Read in Hamiltonian model parameters
    # define number of vibrational model
    # name = "low_freq_model_strong_coup_vibronic_linear"
    name = "CH2O_VECC"

    integrator_flag = "RK"

    hamiltonian_truncation_order = 2

    # model = read_in_model(inputdir, name, order=1)

    model = vIO_wrapper.vibronic_input_reader(name, hamiltonian_truncation_order, inputdir)


    print("number of surfaces:{:}".format(model[VMK.A]))
    print("number of modes:{:}".format(model[VMK.N]))
    print("verticle energy (in eV):\n{:}".format(model[VMK.E]))
    print("Frequencies (in eV):\n{:}".format(model[VMK.w]))
    print("Linear coupling constants (in eV):\n{:}".format(model[VMK.G1]))
    if hamiltonian_truncation_order >=2:
        print("Quadratic coupling constants (in eV):\n{:}".format(model[VMK.G2]))

    # assert np.allclose(model[VMK.G2], np.transpose(model[VMK.G2], (1, 0, 3, 2)))

    # initialize the Hamiltonian
    model = vibronic_model_hamiltonian(model, name, truncation_order=2, FC=False, T_2_flag=True)
    # sys.exit(0)
    # calculate thermal properties using the sum over states method
    # model.sum_over_states(output_path=outputdir, basis_size=40, T_initial=2000, T_final=30, N_step=100)
    # sys.exit(0)
    # Bogoliubov transform the Hamiltonian
    model.thermal_field_transform(T_ref=2e3)
    model.reduce_H_tilde()
    # sys.exit(0)
    # run TFCC simulation
    if integrator_flag == "Euler":
        model.TFCC_integration(T_initial=1e3, T_final=3e1, N_step=3000, output_path=outputdir) #(primary 1st order Euler method)
    elif integrator_flag == "RK":
        func_string = 'model.rk45_integration(T_initial=1e4, T_final=3e1, nof_points=10000, output_path=outputdir)'
        if True:
            # conduct profiling of the main simulation code
            filename = name + "_cProfile_data"
            cProfile.runctx(
                        func_string,
                         globals(),
                         locals(),
                         filename)
            # store the profiling data
            process_data(filename)
        else:
            eval(func_string)
    else:
        pass
    return


if (__name__ == '__main__'):
    main()
