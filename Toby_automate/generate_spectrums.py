#!python3

# system imports
import shutil
from shutil import copyfile
import itertools as it
import glob
import sys
import os
from os.path import join

# third party imports
from scipy.interpolate import interp1d
import numpy as np

# local imports
import prop_input_template
from project_parameters import *
from vibronic import vIO, VMK
#
# sys.path.insert(0, os.path.abspath("/home/ngraymon/public/songhao/t_amplitudes_project/t-amplitudes"))  # TEMPORARY
# import project
# from project.vibronic_hamiltonian import vibronic_hamiltonian


nof_points = 3000

eV_dict = {
    "h2o":  (21, 11),
    "ch2o": (19, 9.5),
    "co2":  (20, 13),
    "n2o":  (21, 11.5),
    "nh3":  (19, 9),
    #
    "h2o2":      (20.5, 8),
    "hcooh":     (20, 8),
    "furan":     (20, 13),
    "formamide": (22, 13),
    "vcm":       (18, 9.5),
    #
    "op_nh36Q_5st":  (11, 2),
    f"{project_name}": (30, 11)
}

y_dict = {
    "h2o":  (-1.5, 60),
    "ch2o": (-1.5, 50),
    "co2":  (-1.5, 100),
    "n2o":  (-1.5, 100),
    "nh3":  (-1.5, 40),
    #
    "h2o2":      (-1.5, 14),
    "hcooh":     (-1.5, 16),
    "furan":     (-1.5, 40),
    "formamide": (-1.5, 15),
    "vcm":       (-1.5, 34),
    #
    "op_nh36Q_5st":       (-1.5, 40),
    f"{project_name}": (-1.5, 100)
}

left_eV, right_EV = eV_dict[project_name]
min_y, max_y = y_dict[project_name]
# min_EV, max_EV = 40, 5
tau = 40
iexp = 1


# -------------------------------------------------------------------------------------------------
def calculate_harmonic_ground_state_of_model(model):
    """Calculate the H.O. contribution to the electronic ground state"""
    harmonic_ground_state = round(0.5 * np.sum(model[VMK.w]), ndigits=9)
    return harmonic_ground_state


def calculate_harmonic_ground_state_of_op_file(path):
    """Calculate the H.O. contribution to the electronic ground state of an '.op' file"""
    model = vIO.extract_excited_state_model_op(f"./{model_name}.op")
    return calculate_harmonic_ground_state_of_model(model)


# -------------------------------------------------------------------------------------------------
# def modify_acf_file(root, path_cc, path_mctdh, mctdh_t_final=None, mctdh_step=None):

#     def extract_from_auto_file(path_auto):
#         t, _ = vibronic_hamiltonian.load_auto_data(path_auto)
#         tf = round(float(t[-1] / 2), ndigits=8)
#         dt = round(float(t[1]-t[0]) / 2, ndigits=8)
#         return tf, dt

#     if mctdh_t_final is None and mctdh_step is None:
#         mctdh_t_final, mctdh_step = extract_from_auto_file(join(root, path_mctdh))

#     try:
#         time, acf = vibronic_hamiltonian.load_acf_data(join(root, path_cc))
#     except Exception as e:
#         print(f"Numerical issues with {path_cc}")
#         print(e)
#         return False

#     print(time.shape)
#     print(time[0], time[1], time[-2], time[-1])
#     print(len(time), len(acf))

#     # if the last two points are the same
#     try:
#         # create interpolation objects
#         f_real = interp1d(time, acf.real, 'cubic')
#         f_imag = interp1d(time, acf.imag, 'cubic')
#     except ValueError as e:
#         print("Most likely the rray has duplicate t values, this needs to be changed when executing")
#         raise e
#         # print(time[-2], time[-1])
#         # time[-1] = round(time[-1] + 1E-7, ndigits=8)
#         # print(time[-2], time[-1])
#         # f_real = interp1d(time, acf.real, 'cubic')
#         # f_imag = interp1d(time, acf.imag, 'cubic')

#     cc_t_init, cc_t_final = time[0], time[-1]
#     print(cc_t_init, cc_t_final, mctdh_t_final)
#     npoints = int(mctdh_t_final / mctdh_step)
#     print(f"Nof points = {npoints}")
#     # if we use endpoint=True then its very difficult to get good deltas
#     # so for now we will use endpoint=False
#     _, dt = np.linspace(cc_t_init, cc_t_final, num=npoints, endpoint=False, retstep=True)
#     time_step = round(dt, ndigits=8)
#     print(f"Time step is {time_step:12.8f}")
#     # then we make sure we generate uniformly spaced points
#     # we use (cc_t_final+time_step) because np.arange doesn't include the last point
#     # so if we want to include cc_t_final (we do) then we need to make the 'last point' one step after the point that we want to be the real 'last point'
#     new_x = np.arange(cc_t_init, cc_t_final+time_step, step=time_step, dtype=float)

#     print(new_x[0], new_x[-1])

#     new_time = new_x.copy()
#     for i in range(len(new_x)):
#         new_time[i] = round(new_x[i], ndigits=8)

#     new_acf = np.zeros_like(new_x, dtype=complex)
#     dt = abs(new_time[0] - new_time[1])
#     arr = abs(new_time[0:-2] - new_time[1:-1])
#     for i, a in enumerate(arr):
#         assert np.isclose(dt, a)

#     # save interpolated data to new_acf array
#     new_acf.real = f_real(new_x)
#     new_acf.imag = f_imag(new_x)

#     # normalize the results (to match MCTDH)
#     normalization_factor = new_acf.real[0]
#     print(f"\nNormalization factor: {normalization_factor}\n")
#     new_acf.real /= normalization_factor
#     new_acf.imag /= normalization_factor

#     # strip '.txt' off path and add the suffix '_normalized.txt'
#     new_path = path_cc[0:-4:]+"_normalized.txt"
#     vibronic_hamiltonian._save_data(join(root, new_path), new_time, new_acf)
#     return new_path


# -------------------------------------------------------------------------------------------------
def _generate_pl(nof_points, root_dir, output_filename, input_filename):
    """Generate command string to call `autospec84` to generate spectrum data points"""
    if "vibronic" in input_filename:
        # if "tf50" in input_filename:
        if "tf100" in input_filename:
            nof_points = 5000
        if "tf500" in input_filename:
            nof_points = 62000
        if "tf1000" in input_filename:
            nof_points = 65000

    command = (
        #
        f"cd {root_dir}; "  # temporary fix!!!!!! TEMPORARY
        #
        "autospec84 "
        # "-g 1 "  # to print gnuplot commands or not
        f"-o {output_filename:s} "
        f"-f {input_filename:s} "
        f"-p {nof_points:d} "
        # f"-EP "
        # f"-e {harmonic_ground_state} eV " # x axis shift (left/right) to account for H.O. G.S.
        f"{left_eV} {right_EV} eV "  # x axis limits (in eV)
        f"{tau:d} "   # tau value
        f"{iexp:d} "  # iexp value
    )
    return command


def generate_cc_pl(nof_points, root_dir, output_filename, input_filename):
    """Generate command string to call `autospec84` to generate spectrum data points for CC"""
    return _generate_pl(nof_points, root_dir, output_filename, input_filename)


def generate_mctdh_pl(nof_points, root_dir, output_filename, input_filename):
    """Generate command string to call `autospec84` to generate spectrum data points for MCTDH"""
    return _generate_pl(nof_points, root_dir, output_filename, input_filename)


# -------------------------------------------------------------------------------------------------
# def write_cc_mctdh_spectrum_plotting_file(configuration, *args):
#     """ a """

#     # unpack arguments
#     nof_points, root_dir, mctdh_file, cc_file, model_name, pbf, t = args
#     # print(f'{root_dir}/{mctdh_file}.pl')
#     # print(f'{root_dir}/{cc_file}.pl')
#     # doctor the file name to make it look better in the plot
#     plot_title = model_name.replace('_', ' ').replace('h2o', 'h_{2}o')

#     size = [1200, 800]
#     # size = [800, 400]

#     output_file = f'{root_dir}/both_spectrum_{model_name:s}_{nof_points:d}_PBF{pbf:d}_{int(t):d}fs_{tau:d}tau.png'

#     plotting_command = '\n'.join([
#         f"set terminal png size {size[0]},{size[1]}",
#         # f"set output './spectrum_{model_name:s}_{nof_points:d}_{t_final:d}fs_{tau:d}tau_{nof_BF}SOSBF_{mctdh_BF}{configuration}BF.png'",
#         f"set output '{output_file:s}'",
#         "set style data line", "set nologscale", "set xzeroaxis", "set xlabel 'Energy[eV]'",
#         f"set xr [ {left_eV}: {right_EV}]",
#         f"set yr [ {min_y}: {max_y}]",
#         f"set title '{plot_title:s} spectrum, n-cos = 1, tau: {tau:d}.0 1, {int(t):d}fs'",
#         f"plot \
#             '{root_dir}/{mctdh_file}.pl' using 1:3 lw 2 lc 'black' title '{configuration}',\
#             '{root_dir}/{cc_file}.pl' every 6 using 1:3 with linespoints lc 'purple' title 'CC',\
#         ",
#         # '{sos_file}.pl' using 1:3 lc 'black' title 'SOS',\
#     ])

#     path_plotting_file = f"{root_dir}/spectrum_plotting.pl"

#     # write the plotting commands to a file
#     with open(path_plotting_file, 'w') as fp:
#         fp.write(plotting_command)

#     return path_plotting_file, output_file


def write_spectrum_plotting_file(configuration, *args):
    """ a """

    # unpack arguments
    nof_points, root_dir, mctdh_file, model_name, pbf, t, operate_string = args
    # print(f'{root_dir}/{mctdh_file}.pl')
    # doctor the file name to make it look better in the plot
    plot_title = model_name.replace('_', ' ').replace('h2o2', 'h_{2}o_{2}')
    plot_title = model_name.replace('_', ' ').replace('h2o', 'h_{2}o')

    size = [1200, 800]
    #size = [800, 400]

    output_file = f'{root_dir}/{configuration}_spectrum_{model_name:s}_{nof_points:d}_PBF{pbf:d}_{int(t):d}fs_{tau:d}tau_init_st{operate_string}.png'

    plotting_command = '\n'.join([
        f"set terminal png size {size[0]},{size[1]}",
        # f"set output './spectrum_{model_name:s}_{nof_points:d}_{t_final:d}fs_{tau:d}tau_{nof_BF}SOSBF_{mctdh_BF}{configuration}BF.png'",
        f"set output '{output_file:s}'",
        "set style data line", "set nologscale", "set xzeroaxis", "set xlabel 'Energy[eV]'",
        # f"set xr [ 0.{left_eV}00000E+02: 0.{right_EV}0000E+02]",
        f"set xr [ {left_eV}: {right_EV}]",
        f"set yr [ {min_y}: {max_y}]",
        f"set title '{plot_title:s} spectrum, n-cos = 1, tau: {tau:d}.0 1, {int(t):d}fs'",
        f"plot \
            '{root_dir}/{mctdh_file}_init_st{operate_string}.pl' using 1:3 lw 2 lc 'black' title '{configuration} g1',\
            '{root_dir}/{mctdh_file}_init_st{operate_string}.pl' using 1:4 lw 2 lc 'red' title '{configuration} g2',\
        ",
        # '{cc_file}.pl' every 6 using 1:3 with linespoints lc 'purple' title 'CC',\
        # '{sos_file}.pl' using 1:3 lc 'black' title 'SOS',\
    ])

    path_plotting_file = f"{root_dir}/spectrum_plotting.pl"

    # write the plotting commands to a file
    with open(path_plotting_file, 'w') as fp:
        fp.write(plotting_command)

    return path_plotting_file, output_file


# def write_mctdh_pbf_spectrum_plotting_file(configuration, *args):
#     """Generate gnuplot script for comparing multiple MCTDH results for different #'s of PBF's """

#     # unpack arguments
#     nof_points, root_dir, spectrum_dir, model_name, t = args
#     # print(f'{root_dir}/{mctdh_file}.pl')
#     # doctor the file name to make it look better in the plot
#     plot_title = model_name.replace('_', ' ').replace('h2o', 'h_{2}o')

#     size = [1200, 800]
#     # size = [800, 400]

#     output_file = f'{spectrum_dir}/combined_{configuration}_spectrum_{model_name:s}_{nof_points:d}_{int(t):d}fs_{tau:d}tau.png'
#     d_30 = dir_string.format(model_name, 30, t)
#     d_100 = dir_string.format(model_name, 100, t)
#     d_300 = dir_string.format(model_name, 300, t)
#     d_500 = dir_string.format(model_name, 500, t)
#     d_1000 = dir_string.format(model_name, 1000, t)

#     plotting_command = '\n'.join([
#         f"set terminal png size {size[0]},{size[1]}",
#         # f"set output './spectrum_{model_name:s}_{nof_points:d}_{t_final:d}fs_{tau:d}tau_{nof_BF}SOSBF_{mctdh_BF}MCTDHBF.png'",
#         f"set output '{output_file:s}'",
#         "set style data line", "set nologscale", "set xzeroaxis", "set xlabel 'Energy[eV]'",
#         # f"set xr [ 0.{left_eV}00000E+02: 0.{right_EV}0000E+02]",
#         f"set xr [ {left_eV}: {right_EV}]",
#         f"set yr [ {min_y}: 30.0]",
#         f"set title '{plot_title:s} spectrum, n-cos = 1, tau: {tau:d}.0 1, {int(t):d}fs'",
#         f"plot \
#             '{root_dir}/{d_30}/{configuration}_spectrum_{d_30}.pl' using 1:3 lw 2 lc 'red' title '30 PBF',\
#             '{root_dir}/{d_100}/{configuration}_spectrum_{d_100}.pl' using 1:3 lw 2 lc 'green' title '100 PBF',\
#             '{root_dir}/{d_300}/{configuration}_spectrum_{d_300}.pl' using 1:3 lw 2 lc 'yellow' title '300 PBF',\
#             '{root_dir}/{d_500}/{configuration}_spectrum_{d_500}.pl' using 1:3 lw 2 lc 'blue' title '500 PBF',\
#             '{root_dir}/{d_1000}/{configuration}_spectrum_{d_1000}.pl' using 1:3 lw 2 lc 'black' title '1000 PBF',\
#         ",
#         # '{cc_file}.pl' every 6 using 1:3 with linespoints lc 'purple' title 'CC',\
#         # '{sos_file}.pl' using 1:3 lc 'black' title 'SOS',\
#     ])

#     path_plotting_file = f"{spectrum_dir}/spectrum_plotting.pl"

#     # write the plotting commands to a file
#     with open(path_plotting_file, 'w') as fp:
#         fp.write(plotting_command)

#     return path_plotting_file, output_file


# def write_mctdh_tf_spectrum_plotting_file(configuration, *args):
#     """Generate gnuplot script for comparing multiple MCTDH results for different lengths of propagation """

#     # unpack arguments
#     nof_points, root_dir, mctdh_file, model_name, pbf, t = args
#     print(f'{root_dir}/{mctdh_file}.pl')
#     # doctor the file name to make it look better in the plot
#     plot_title = model_name.replace('_', ' ').replace('h2o', 'h_{2}o')

#     size = [1200, 800]
#     # size = [800, 400]

#     plotting_command = '\n'.join([
#         f"set terminal png size {size[0]},{size[1]}",
#         # f"set output './spectrum_{model_name:s}_{nof_points:d}_{t_final:d}fs_{tau:d}tau_{nof_BF}SOSBF_{mctdh_BF}MCTDHBF.png'",
#         f"set output '{root_dir}/spectrum_{model_name:s}_{nof_points:d}_PBF{pbf:d}_{int(t):d}fs_{tau:d}tau.png'",
#         "set style data line", "set nologscale", "set xzeroaxis", "set xlabel 'Energy[eV]'",
#         # f"set xr [ 0.{left_eV}00000E+02: 0.{right_EV}0000E+02]",
#         f"set xr [ {left_eV}: {right_EV}]",
#         f"set yr [ {min_y}: {max_y}]",
#         f"set title '{plot_title:s} spectrum, n-cos = 1, tau: {tau:d}.0 1, {int(t):d}fs'",
#         f"plot \
#             '{root_dir}/{mctdh_file}.pl' using 1:3 lw 2 lc 'black' title '{configuration}',\
#         ",
#         # '{cc_file}.pl' every 6 using 1:3 with linespoints lc 'purple' title 'CC',\
#         # '{sos_file}.pl' using 1:3 lc 'black' title 'SOS',\
#     ])

#     path_plotting_file = "spectrum_plotting.pl"

#     # write the plotting commands to a file
#     with open(path_plotting_file, 'w') as fp:
#         fp.write(plotting_command)

#     return path_plotting_file


# def write_mctdh_coupling_spectrum_plotting_file(configuration, *args):
#     """Generate gnuplot script for comparing multiple MCTDH results between constant/linear/quadratic coupling terms """

#     # unpack arguments
#     nof_points, root_dir, mctdh_file, model_name, pbf, t = args
#     print(f'{root_dir}/{mctdh_file}.pl')
#     # doctor the file name to make it look better in the plot
#     plot_title = model_name.replace('_', ' ').replace('h2o', 'h_{2}o')

#     size = [1200, 800]
#     # size = [800, 400]

#     plotting_command = '\n'.join([
#         f"set terminal png size {size[0]},{size[1]}",
#         # f"set output './spectrum_{model_name:s}_{nof_points:d}_{t_final:d}fs_{tau:d}tau_{nof_BF}SOSBF_{mctdh_BF}MCTDHBF.png'",
#         f"set output '{root_dir}/spectrum_{model_name:s}_{nof_points:d}_PBF{pbf:d}_{int(t):d}fs_{tau:d}tau.png'",
#         "set style data line", "set nologscale", "set xzeroaxis", "set xlabel 'Energy[eV]'",
#         # f"set xr [ 0.{left_eV}00000E+02: 0.{right_EV}0000E+02]",
#         f"set xr [ {left_eV}: {right_EV}]",
#         f"set yr [ {min_y}: {max_y}]",
#         f"set title '{plot_title:s} spectrum, n-cos = 1, tau: {tau:d}.0 1, {int(t):d}fs'",
#         f"plot \
#             '{root_dir}/{mctdh_file}.pl' using 1:3 lw 2 lc 'black' title '{configuration}',\
#         ",
#         # '{cc_file}.pl' every 6 using 1:3 with linespoints lc 'purple' title 'CC',\
#         # '{sos_file}.pl' using 1:3 lc 'black' title 'SOS',\
#     ])

#     path_plotting_file = "spectrum_plotting.pl"

#     # write the plotting commands to a file
#     with open(path_plotting_file, 'w') as fp:
#         fp.write(plotting_command)

#     return path_plotting_file


# def write_acf_plotting_file(configuration, *args):
#     """ a """

#     # unpack arguments
#     # nof_points, mctdh_file, model_name = args
#     # nof_points, cc_file, mctdh_file, cc_file2, model_name = args

#     # doctor the file name to make it look better in the plot
#     plot_title = model_name.replace('_', ' ').replace('h2o', 'h_{2}o')

#     size = [1200, 800]
#     # style = "circle radius graph 0.002"
#     # style = "line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5"
#     style = "line 1 lt 2 pt 12 ps 1 pi -1"

#     plotting_command = '\n'.join([
#         f"set terminal png size {size[0]},{size[1]}",
#         f"set title 'ACF comparison of {t_final:d}fs'",
#         f"set style {style}", "set style data line",
#         # "set pointintervalbox 2",
#         f"set output './ACF_{model_name:s}_{nof_points:d}_{t_final:d}fs.png'",
#         "set nologscale", "set xzeroaxis", "set ylabel 'C(tau/hbar)'",
#         "set yr [ -1: 1]",
#         "set xr [ 0.0: 100.0]",
#         f"plot \
#             '{mctdh_file}' us 1:2 with linespoints ls 4 ps 3 lc 'red' title '{configuration} Real',\
#             '{cc_file}' us 1:2 lc 'green' title 'CC Real (interpolated)',\
#         ",
#         # '{cc_file2}' us 1:2 with linespoints ls 1 lc 'blue' title 'CC Real (raw RK45)'\
#     ])

#     path_plotting_file = "acf_plotting.pl"

#     # write the plotting commands to a file
#     with open(path_plotting_file, 'w') as fp:
#         fp.write(plotting_command)

#     return path_plotting_file


# -------------------------------------------------------------------------------------------------
def create_spectrum():
    return


def create_spectrums():
    return


# -------------------------------------------------------------------------------------------------
def mctdh_job_is_finished(path_to_mctdh_execution_folder):
    """Return True if MCTDH job is finished, otherwise return False."""

    path_to_output_file = join(path_to_mctdh_execution_folder, "output")

    if os.path.exists(path_to_output_file):
        with open(path_to_output_file, 'r') as fp:
            data = fp.read()

        return bool("Propagation was successful." in data)

    return False
# -------------------------------------------------------------------------------------------------


if __name__ == "__main__":

    spectrum_dir = join(work_root, "0_spectrums/")
    auto_dir = join(work_root, "0_auto_dir/")
    os.makedirs(spectrum_dir, exist_ok=True)
    os.makedirs(auto_dir, exist_ok=True)

    configuration = "tdh" if ("/tdh/" in work_root) else "mctdh"

    only_checking_output = False
    suppress_autosped84_ouptut = False
    generating_cc = False

    # temporarily plot combined MCTDH
    if False:
        coupling_order = "quadratic"
        FC_or_not = "vibronic"
        model_name = "_".join([project_name, FC_or_not, coupling_order])

        for t in t_final:
            plotting_file, output_file = write_mctdh_pbf_spectrum_plotting_file(
                nof_points, work_root, spectrum_dir, model_name, t
            )
            os.system(f"gnuplot {plotting_file}")

        sys.exit(0)

    failed_jobs = {}
    running_jobs = {}
    completed_jobs = {}
    not_submitted_jobs = {}

    # os.system('module load mctdh/84.16')  # doesn't work - need special module loading thing

    model_name = project_name

    for param_list in it.product(*expression_list):
        print(param_list)

        calculation_spec = dir_string.format(model_name, *param_list)
        for operate_string in range(1, A+2):
            root_dir = join(work_root, calculation_spec, f"init_st{operate_string}")
            print(f'root_dir: {root_dir}')
    
            path_mctdh_acf = join(model_name, "auto")
            #print("path_mctdh_acf", path_mctdh_acf, "\n")
    
            path_mctdh_spectrum = f"{configuration}_spectrum_{dir_string.format(model_name, *param_list)}"
    
            #print("path_mctdh_spectrum", path_mctdh_spectrum, "\n")
            # auto_path = f"./h2o_FC_{order:s}_{25:>03d}fs_{BF:>03d}BF_{spf:>02d}spf/auto"
    
            path_to_mctdh_execution_folder = join(root_dir, model_name)
    
            # if no folder, job was not even submitted
            if not os.path.isdir(path_to_mctdh_execution_folder):
                print(f"{configuration.upper()} job has not yet been submitted {model_name} Parameters: {param_list}")
                pbf, tf = param_list
                name = f"pbf{pbf}"
                if name not in not_submitted_jobs:
                    not_submitted_jobs[name] = [tf, ]
                else:
                    not_submitted_jobs[name].append(tf)
                continue
    
            # if slurm seems done then job might be finished
            elif mctdh_job_is_finished(path_to_mctdh_execution_folder):
                command = generate_mctdh_pl(
                    nof_points,
                    root_dir=root_dir,
                    output_filename=f"{path_mctdh_spectrum}_init_st{operate_string}",
                    input_filename=path_mctdh_acf
                )

                print(command)
                # if we just want to check how many jobs failed / haven't been submitted
                if only_checking_output:
                    continue
    
                if suppress_autosped84_ouptut:
                    os.system(command + ' 1> /dev/null')
                else:
                    print('\n')
                    os.system(command)
                    print('\n')
    
                print(f'Autospecd {model_name} Parameters: {param_list}')
    
                pbf, tf = param_list
                name = f"pbf{pbf}"
                if name not in completed_jobs:
                    completed_jobs[name] = [tf, ]
                else:
                    completed_jobs[name].append(tf)
    
                # ask Neil
                # if configuration == "mctdh":
                #     # move very large files to new folder
                #     old_dir = abspath(f"/work/{user_root}/{parent_project}/mctdh/{project_name}/{calculation_spec}/{model_name}/")
                #     new_dir = abspath(f"/work/{user_root}/{parent_project}/nearline/{project_name}/{calculation_spec}/")
    
                #     for file_name in ["psi", "dvr", "restart", "oper"]:
                #         target = join(old_dir, file_name)
                #         if os.path.exists(target):
                #             os.makedirs(new_dir, exist_ok=True)
                #             # print(f'mv -i {target} {new_dir}')
                #             os.system(f'mv -i {target} {new_dir}')
                #         else:
                #             print(f"No file here: {target}")
    
            else:
                # check for issue in latest slurm output
                paths = glob.glob(join(root_dir, 'slurm-*.out'))
                latest_slurm_file = sorted(paths, reverse=True)[0]
    
                if not os.path.isfile(latest_slurm_file):
                    print(f"Malformed slurm file {latest_slurm_file} Parameters: {param_list}")
                    continue
    
                else:
                    with open(latest_slurm_file, 'r') as fp:
                        data = fp.read()
    
                    pbf, tf = param_list
                    name = f"pbf{pbf}"
    
                    lowercase_data = data.lower()
                    error_string_list = ['memory limit', 'killed', 'error', 'CANCELLED', ]
    
                    for string in error_string_list:
                        if string in lowercase_data:
                            if 'memory limit' in data:
                                print(f"FAIL: {configuration.upper()} job failed due to memory limits {model_name} Parameters: {param_list}")
                            elif 'CANCELLED' in data:
                                print(f"FAIL: {configuration.upper()} job was cancelled {model_name} Parameters: {param_list}")
                            else:
                                print(f"FAIL: {configuration.upper()} job failed for some reason {model_name} Parameters: {param_list}")
    
                            # store job
                            if name not in failed_jobs:
                                failed_jobs[name] = [tf, ]
                            else:
                                failed_jobs[name].append(tf)
    
                            if configuration == "mctdh":
                                # ldir = abspath(f"/work/{user_root}/{parent_project}/mctdh/{project_name}/{calculation_spec}/{model_name}/")
                                ldir = abspath(f"/work/{user_root}/mctdh/{project_name}/{calculation_spec}/{model_name}/")
                                # delete the super large files
                                for file_name in ["psi", "dvr", "restart", "oper"]:
                                    target = join(ldir, file_name)
                                    if os.path.exists(target):
                                        print(f'rm {target}')
                                        os.system(f'rm {target}')
                            continue
    
                    else:
                        print(f"{configuration.upper()} job is not finished! {model_name} Parameters: {param_list}")
    
                        # store job
                        if name not in running_jobs:
                            running_jobs[name] = [tf, ]
                        else:
                            running_jobs[name].append(tf)
    
                    continue
    
            # copy auto files over
            src_acf_path = join(path_to_mctdh_execution_folder, 'auto')
            print(f'src_acf_path: {src_acf_path}')
            dst_acf_path = join(auto_dir, f"auto_{calculation_spec}_init_st{operate_string}")
            print(f'dst_acf_path: {dst_acf_path}')
            shutil.copy(src_acf_path, dst_acf_path)
    
            # ask neil
            # turned off
            # if generating_cc:
            #     if os.path.isfile(join(root_dir, f"ACF_CC_{model_name:s}_tf{int(param_list[1]):>03d}.txt")):
            #         path_cc_spectrum = f"cc_spectrum_{model_name:s}_{int(param_list[1]):>03d}fs"
            #         path_cc_acf = modify_acf_file(
            #             root_dir,
            #             f"ACF_CC_{model_name:s}_tf{int(param_list[1]):>03d}.txt",
            #             path_mctdh_acf
            #         )
            #         if path_cc_acf is False:
            #             continue
    
            #         command = generate_cc_pl(
            #             nof_points,
            #             root_dir=root_dir,
            #             output_filename=path_cc_spectrum,
            #             input_filename=path_cc_acf
            #         )
            #         os.system(command)
    
            #     else:
            #         print(f"\nCC job is not finished! {model_name}\nParameters: {param_list}\n")
            #         continue
    
            if True:  # plotting singular spectrums
                plotting_file, output_file = write_spectrum_plotting_file(
                    configuration, nof_points, root_dir, path_mctdh_spectrum, model_name, *param_list, operate_string
                )
                os.system(f"gnuplot {plotting_file}")
                shutil.copy(output_file, spectrum_dir)
    
            # if False:  # plotting combined spectrums
            #     plotting_file, output_file = write_cc_mctdh_spectrum_plotting_file(
            #         configuration, nof_points, root_dir, path_mctdh_spectrum, path_cc_spectrum, model_name, *param_list
            #     )
            #     os.system(f"gnuplot {plotting_file}")
            #     shutil.copy(output_file, spectrum_dir)
    
            # if False:  # plotting ACF
            #     plotting_file = write_acf_plotting_file(configuration, nof_points, new_cc_acf, auto_path)
            #     os.system(f"gnuplot {plotting_file}")
    
            # if False:  # plotting ACF vs SOS
            #     plotting_file = write_acf_sos_plotting_file(configuration, nof_points, new_cc_acf, auto_path, new_sos_acf)
            #     os.system(f"gnuplot {plotting_file}")

        if True:
            try:
                for operate_string in range(1, A+2):
                    root_dir_z = join(work_root, calculation_spec, f'init_st{operate_string}')
    
                    print(root_dir_z)
                    # Load data from 'auto_z' file if it exists
                    auto_z_path = f"{root_dir_z}/{path_mctdh_spectrum}_init_st{operate_string}.pl"
                    if os.path.exists(auto_z_path):
                        print(auto_z_path, 'exists')
                        auto_z = np.loadtxt(auto_z_path)
                    else:
                        auto_z = None
        
            except Exception as e:
                print(f"Error loading files: {e}")
        
            try:
                # Find a non-empty array among auto_x, auto_y, and auto_z
                non_empty_array = next(arr for arr in [auto_z] if arr is not None)
        
                # Extract x values (assuming the x values are the same in all files)
                x_values = non_empty_array[:, 0]
        
                # Sum the 4th column across the three files
                sum_column_4 = np.zeros_like(x_values)
        
                # Add the 4th column of each existing array to the sum
                for arr in [auto_z]:
                    if arr is not None:
                        sum_column_4 += arr[:, 3]
        
                # Create a new array with x values and summed 4th column
                auto_total = np.column_stack((x_values, sum_column_4))
        
                # Save the result to 'auto_total' file
                np.savetxt(f'{calculation_spec}_auto_total', auto_total, header='x_values  sum_column_4', fmt='%12.8f', comments='')
        
                print("Process completed. Check 'auto_total' file for the result.")
        
            except Exception as e:
                print(f"Error in creating auto_total: {e}")
    
        size = [1200, 800]
    
        output_file = f'{calculation_spec}_xyz_spectrum.png'
    
        plotting_command = '\n'.join([
            f"set terminal png size {size[0]},{size[1]}",
            # f"set output './spectrum_{model_name:s}_{nof_points:d}_{t_final:d}fs_{tau:d}tau_{nof_BF}SOSBF_{mctdh_BF}{configuration}BF.png'",
            f"set output '{output_file:s}'",
            "set style data line", "set nologscale", "set xzeroaxis", "set xlabel 'Energy[eV]'",
            # f"set xr [ 0.{left_eV}00000E+02: 0.{right_EV}0000E+02]",
            f"set xr [ {left_eV}: {right_EV}]",
            f"set yr [ {min_y}: {max_y}]",
            f"set title '{project_name:s} spectrum, n-cos = 1, tau: {tau:d}.0 1'",
            #, {int(t):d}fs'",
            f"plot '{calculation_spec}_auto_total' using 1:2 lw 2 lc 'black' title 'mctdh g1'",
            #     '{root_dir}/{mctdh_file}_{operate_string}.pl' using 1:3 lw 2 lc 'black' title '{configuration} g1',\
            #     '{root_dir}/{mctdh_file}_{operate_string}.pl' using 1:4 lw 2 lc 'red' title '{configuration} g2',\
            # ",
            # '{cc_file}.pl' every 6 using 1:3 with linespoints lc 'purple' title 'CC',\
            # '{sos_file}.pl' using 1:3 lc 'black' title 'SOS',\
        ])
    
        path_plotting_file = f"{calculation_spec}_spectrum_plotting.pl"
    
        # write the plotting commands to a file
        with open(path_plotting_file, 'w') as fp:
            fp.write(plotting_command)
    
        os.system(f"gnuplot {path_plotting_file}")
        print('Please check png:', output_file)

    print(f"{'-'*35}  Successfully Completed jobs  {'-'*35}\n")
    for key in completed_jobs.keys():
        print(f"{key}, tf:", completed_jobs[key])

    print(f"{'-'*35}  Not Submitted jobs  {'-'*35}\n")
    for key in not_submitted_jobs.keys():
        print(f"{key}, tf:", not_submitted_jobs[key])

    print(f"{'-'*35}  Failed jobs  {'-'*35}\n")
    for key in failed_jobs.keys():
        print(f"{key}, tf:", failed_jobs[key])

    print(f"{'-'*35}  Running jobs  {'-'*35}\n")
    for key in running_jobs.keys():
        print(f"{key}, tf:", running_jobs[key])
