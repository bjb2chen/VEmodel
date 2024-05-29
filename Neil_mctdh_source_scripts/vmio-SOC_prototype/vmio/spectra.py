"""spectra.py should handle any I/O and processing related to producing spectra"""

# system imports
from os.path import isfile, isdir, dirname, join


# third party imports
from scipy.interpolate import interp1d
import numpy as np


# local imports
from .log_conf import log
from .vibronic import vIO, VMK
from .vibronic_hamiltonian import vibronic_hamiltonian

""" Technical note about the use of python's builtin `round` function versus numpy's `np.round` function:
Based on
https://stackoverflow.com/questions/35215161/most-efficient-way-to-map-function-over-numpy-array
and
https://numpy.org/doc/stable/reference/generated/numpy.around.html#numpy.around
it appears that python's `round` function is more accurate
and therefore preferred to use when rounding to 8 digits in `_interpolate`.
"""

# -------------------------------------------------------------------------------------------------
# helper functions
# -------------------------------------------------------------------------------------------------
def calculate_harmonic_ground_state_of_model(model):
    """Calculate the H.O. contribution to the electronic ground state"""
    harmonic_ground_state = round(0.5 * np.sum(model[VMK.w]), ndigits=9)
    return harmonic_ground_state


def calculate_harmonic_ground_state_of_op_file(path):
    """Calculate the H.O. contribution to the electronic ground state of an '.op' file"""
    model = vIO.extract_excited_state_model_op(f"./{model_name}.op")
    return calculate_harmonic_ground_state_of_model(model)


def modify_plotting_code_for_current_file(output_filename, plotting_string):
    """Attempt to overwrite the first 7 lines of `*.pl` plotting file.
    The first 7 lines are the gnuplot instructions and the rest of the lines contain the data to be plotted.
    Attempts to read the data from the file and then write back to the same file but replacing the first 7 lines with `plotting_string`.
    """
    try:
        with open(f"./{output_filename}.pl", 'r') as fp:
            data = fp.readlines()

        new_file = ''.join([plotting_string, '\n', *data[7:]])
        with open(f"./{output_filename}.pl", 'w') as fp:
            fp.write(new_file)

    except Exception as e:
        print("Something failed")
        raise e
    return


# -------------------------------------------------------------------------------------------------
# Interpolation functions
# -------------------------------------------------------------------------------------------------
def _interpolate(time, acf, nof_points):
    """ Return a tuple (`acf`,`time`) of numpy arrays `acf` is complex128, and `time` is float64.

    Interpolates the function `acf(time)` over `nof_points` equidistant points generated using `np.linspace`.
    Uses the `scipy.interpolate.interp1d` method to preform interpolation.
    """

    # some debug output
    log.debug(f"Number of data points: {len(time)}")
    log.debug(
        f"Original time is:\n"
        f"{'t0:':<5}{time[0]:12.8f}\n"
        f"{'t1:':<5}{time[1]:12.8f}\n"
        f"{'tf:':<5}{time[-1]:12.8f}\n"
    )
    assert len(time) == len(acf.real), \
        f"Length of time ({len(time)}) and acf.real ({len(acf.real)}) is not the same."

    # ---------------------------------------------------------------------------------------------
    # First we need to build the interpolation objects
    try:
        # create interpolation objects
        f_real = interp1d(time, acf.real, 'cubic')
        f_imag = interp1d(time, acf.imag, 'cubic')
    except ValueError as e:
        log.warning("Most likely the array has duplicate t values, this needs to be changed when executing")
        raise e
        if False:  # a debugging test, should remove/fix in future
            log.debug(time[-2], time[-1])
            time[-1] = round(time[-1] + 1E-7, ndigits=8)
            log.debug(time[-2], time[-1])
            f_real = interp1d(time, acf.real, 'cubic')
            f_imag = interp1d(time, acf.imag, 'cubic')

    # ---------------------------------------------------------------------------------------------
    # Then we need to generate equidistant grid points to interpolate over.

    """ We start by calculating the spacing of the grid points.
    Note: if we use `endpoint=True` then it is difficult to get good deltas
    so for now we will use `endpoint=False`.
    """
    t_init, t_final = time[0], time[-1]
    log.debug(f"Number of data points we are interpolating over = {nof_points}")
    _, dt = np.linspace(t_init, t_final, num=nof_points, endpoint=False, retstep=True)
    # we round to 8 digits since this is the limit of the precision of the `autospec84` program
    spacing = round(dt, ndigits=8)
    log.debug(f"Grid point spacing is {spacing:12.8f}")

    # ---------------------------------------------------------------------------------------------
    # Next we generate the grid points.

    """ Note: we use (`t_final`+`spacing`) because `np.arange` doesn't include the last point
    so if we want to include `t_final` (we do) then we need to make the "last point"
    one step after the point that we want to be the real "last point".
    """
    temp_grid = np.arange(t_init, t_final+spacing, step=spacing, dtype=float)

    # again we round to 8 digits
    new_time = np.array([round(x, ndigits=8) for x in temp_grid])
    log.debug(f"Grid points, start:{new_time[0]} stop:{new_time[-1]}")

    # Checking for numerical issues with the rounding and/or grid points
    # loop over every pair of grid points `x_i`, `x_i+1` to confirm they are `dt` apart
    dt, arr = abs(new_time[0] - new_time[1]), abs(new_time[0:-2] - new_time[1:-1])
    for i, a in enumerate(arr):
        assert np.isclose(dt, a)

    # ---------------------------------------------------------------------------------------------
    # Preform the interpolation.

    new_acf = np.zeros_like(new_time, dtype=complex)
    try:
        new_acf.real = f_real(new_time)
        new_acf.imag = f_imag(new_time)
    except ValueError as e:
        new_acf.real = f_real(temp_grid)
        new_acf.imag = f_imag(temp_grid)
        log.debug("Had another stupid interpolation issue")

    return new_acf, new_time


def mctdh_interpolation(cc_time, cc_acf, mctdh_t_final, mctdh_dt):
    """ Return a tuple (`cc_acf`,`cc_time`) of numpy arrays `cc_acf` is complex128, and `cc_time` is float64.

    Interpolates auto-correlation function(ACF) points over MCTDH time steps so that the two data sets can be compared.
    Returns the interpolated result as a tuple ACF(t),t.
    """

    # the number of data points MCTDH generated, we want to fit our cc data to this
    nof_mctdh_points = int(mctdh_t_final / mctdh_dt)

    # generate the new data points
    new_acf, new_time = _interpolate(cc_time, cc_acf, nof_mctdh_points)

    # normalize the results (to match MCTDH)
    normalization_factor = new_acf.real[0]
    log.debug(f"\nNormalization factor: {normalization_factor}\n")
    new_acf.real /= normalization_factor
    new_acf.imag /= normalization_factor

    return new_acf, new_time


def generate_normalized_acf_results(root, path_cc, path_mctdh, mctdh_t_final=None, mctdh_dt=None):
    """ Modify the original ACF(t) points to match MCTDH results and write them to a file. """

    def extract_from_auto_file(path_auto):
        """ Return a tuple of the final time and change in time for a MCTDH 'auto' file. """
        t, _ =  vibronic_hamiltonian.load_auto_data(path_auto)
        # we divide both by 2 because MCTDH takes advantage of symmetry in the Real Hamiltonian to output results
        # using half as many time steps as specified in the input file for the calculation.
        # However the CC integration cannot take advantage of this symmetry and needs to use the original tf and dt.
        tf = round(float(t[-1] / 2), ndigits=8)
        dt = round(float(t[1]-t[0]) / 2, ndigits=8)
        return tf, dt

    # extract the MCTDH data from the auto file if not provided through keyword arguments
    if mctdh_t_final is None and mctdh_dt is None:
        mctdh_t_final, mctdh_dt = extract_from_auto_file(join(root, path_mctdh))

    # read in the original data points we will be interpolating
    try:
        cc_time, cc_acf = vibronic_hamiltonian.load_acf_data(join(root, path_cc))
    except Exception as e:
        print(f"Cannot extract data from:\n{path_cc}\n")
        print(e)
        raise(e)

    # preform the interpolation
    new_cc_acf, new_cc_time = mctdh_interpolation(cc_time, cc_acf, mctdh_t_final, mctdh_dt)

    # save the data to disk
    # we strip '.txt' off the original path and add the suffix '_normalized.txt'
    new_path = path_cc[0:-4:] + "_normalized.txt"
    output_path = join(root, new_path)
    vibronic_hamiltonian._save_data(output_path, new_cc_time, new_cc_acf)
    return output_path


# -------------------------------------------------------------------------------------------------
# generate spectra commands for `autospec84`
# -------------------------------------------------------------------------------------------------
def _generate_pl(**kwargs):
    """Generate command string to call `autospec84` to generate spectrum data points"""

    # some basic error checking to help debugging since we're using kwargs
    # assert isdir(dirname(kwargs['output_filename'])), \
    #     f"Bad output path for spectrum file:\n{kwargs['output_filename']}\n"
    # assert isfile(kwargs['input_filename']), \
    #     f"Input filename is not a valid file!?\n{kwargs['input_filename']}\n"
    assert isinstance(kwargs['nof_points'], int)
    assert isinstance(kwargs['tau'], int)
    assert isinstance(kwargs['iexp'], int)

    command = (
        "autospec84 "
        # "-g 1 "  # to print gnuplot commands or not
        f"-o {kwargs['output_filename']:s} "
        f"-f {kwargs['input_filename']:s} "
        f"-p {kwargs['nof_points']:d} "
        # f"-EP "
        # f"-e {harmonic_ground_state} eV " # x axis shift (left/right) to account for H.O. G.S.
        f"{kwargs['left_eV']} {kwargs['right_eV']} eV "  # x axis limits (in eV)
        f"{kwargs['tau']:d} "   # tau value
        f"{kwargs['iexp']:d} "  # iexp value
    )
    return command


def generate_cc_pl(**kwargs):
    """Generate command string to call `autospec84` to generate spectrum data points for CC"""
    return _generate_pl(**kwargs)


def generate_mctdh_pl(**kwargs):
    """Generate command string to call `autospec84` to generate spectrum data points for MCTDH"""
    return _generate_pl(**kwargs)


# -------------------------------------------------------------------------------------------------
def modify_acf_file(root, path_cc, path_mctdh, mctdh_t_final=None, mctdh_step=None):
    return None
