# system imports
import os
from os.path import join
import sys
# import itertools as it
# import warnings

# third party imports
# import scipy
import numpy as np

# local imports

# import the path to the package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
import project
from project.vibronic import vIO, VMK
from project.vibronic_hamiltonian import vibronic_hamiltonian

model_name = 'h2o_FC_linear'
t_final = 50
tau = 20
nof_points = 4000
nof_BF = 50

auto_path = f"./h2o_FC_linear/auto"
# auto_path = f"./{model_name}/auto"


def modify_single_file(output_filename, plotting_string):
    # try to modify the file
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


def modify_acf_file(path, mctdh_t_final=25.0, mctdh_step=0.1):
    time, acf = vibronic_hamiltonian.load_acf_data(path)

    # h_factor = 4.135667696
    hbar_factor = 6.582119569E-1
    # acf.real *= np.exp(1.0 / hbar_factor)
    # acf.imag *= np.exp(1.0 / hbar_factor)

    print(time[0], time[-2], time[-1])

    t_init, t_final = time[0], time[-1]
    factor = round(t_final / mctdh_t_final, ndigits=8)
    print("factor", factor)
    print("factor2", (2 / factor))
    print("Normalization Factor?", acf.real[0])
    # time /= 2
    # time *= hbar_factor
    # time *= (2 / factor)  # works for original 2x2 (50fs vs 12.5fs) and 100fs vs 25fs
    print(time[0], time[-2], time[-1])

    # acf *= np.exp(1j * (time / 2) * 0.3446135)
    y_real = acf.real
    y_imag = acf.imag

    from scipy.interpolate import interp1d
    f_real = interp1d(time, y_real, 'cubic')
    f_imag = interp1d(time, y_imag, 'cubic')

    t_init, t_final = time[0], time[-1]
    print(t_init, t_final)
    # if we use endpoint=True then its very difficult to get good deltas
    # so for now we will use endpoint=False
    print(t_final, mctdh_t_final)
    factor = round(t_final / mctdh_t_final, ndigits=8)
    npoints = int(mctdh_t_final / mctdh_step)
    print(f"Factor = {factor}, Nof points = {npoints}")

    # if we use endpoint=True then its very difficult to get good deltas
    _, dt = np.linspace(t_init, t_final, num=npoints, endpoint=False, retstep=True)
    time_step = round(dt, ndigits=8)
    print(f"Time step is {time_step:12.8f}")
    # then we make sure we generate uniformly spaced points
    # we use (t_final+time_step) because np.arange doesn't include the last point
    # so if we want to include t_final (we do) then we need to make the 'last point' one step after the point that we want to be the real 'last point'
    new_x = np.arange(t_init, t_final, step=time_step, dtype=float)
    print(new_x[0], new_x[-1])

    time = new_x.copy()
    for i in range(len(new_x)):
        time[i] = round(new_x[i], ndigits=8)

    acf = np.zeros_like(new_x, dtype=complex)
    dt = abs(time[0] - time[1])
    arr = abs(time[0:-2] - time[1:-1])
    for i, a in enumerate(arr):
        assert np.isclose(dt, a)

    acf.real = f_real(new_x)
    acf.imag = f_imag(new_x)

    # acf *= np.exp(1j * (new_x / 8) * 1.0)

    normalization_factor = acf.real[0]
    print(f"\nNormalization factor: {normalization_factor}\n")
    acf.real /= normalization_factor
    acf.imag /= normalization_factor

    new_path = path[0:-4:]+"_normalized.txt"
    vibronic_hamiltonian._save_data(new_path, time, acf)
    return new_path


def generate_cc_pl(nof_points, output_filename, input_filename):
    # model = vIO.extract_excited_state_model_op(f"./{model_name}.op")
    # harmonic_ground_state = round(0.5 * np.sum(model[VMK.w]), ndigits=9)
    # print(f"Harmonic ground state: {harmonic_ground_state}")
    command = (
        "autospec84 "
        # "-g 1 "  # to print gnuplot commands or not
        f"-o {output_filename:s} "
        f"-f {input_filename:s} "
        f"-p {nof_points:d} "
        # f"-EP "
        # f"-e {harmonic_ground_state} eV "
        "21 10 eV "
        # "40 5 eV "
        f"{tau:d} "  # tau value
        "1 "  # iexp value
    )
    print(command)
    os.system(command)

    # set default value
    if "-q " not in command and "-p " not in command:
        nof_points = 1000
        print(f"Defaulting to {nof_points} points.")

    return nof_points


def generate_mctdh_pl(nof_points, output_filename, input_filename=auto_path):
    model = vIO.extract_excited_state_model_op(f"./{model_name}.op")
    harmonic_ground_state = -round(0.5 * np.sum(model[VMK.w]), ndigits=9)
    print(f"Harmonic ground state: {harmonic_ground_state}")
    # model = vIO.read_model_op_file(f"./{model_name}.op")
    # ground_state = -model[VMK.E][-1, -1]
    command = (
        "autospec84 "
        # "-g 1 "  # to print gnuplot commands or not
        f"-o {output_filename:s} "
        f"-f {input_filename:s} "
        f"-p {nof_points:d} "
        # f"-EP "
        f"-e {harmonic_ground_state} eV "
        # "-e 0.591746 eV "
        # "-e 0.344614 eV "
        "21 10 eV "
        # "40 5 eV "
        f"{tau:d} "  # tau value
        "1 "  # iexp value
    )
    print(command)
    os.system(command)
    return

# x1 = f"cc_spectrum_12fs.pl"
# x2 = f"cc_spectrum_25fs.pl"
# x3 = f"cc_spectrum_25fs_halved.pl"
# x4 = f"cc_spectrum_50fs_not_halved.pl"
# x5 = f"cc_spectrum_50fs_halved.pl"
# x6 = f"cc_spectrum.pl"


def write_spectrum_plotting_file(nof_points, cc_file, mctdh_file, sos_file):
    # plotting command
    plotting_string = '\n'.join([
        # "set terminal png size 800,400",
        "set terminal png size 1200,800",
        f"set output './aaaa_spectrum_{model_name:s}_{nof_points:d}_{t_final:d}fs_{tau:d}tau.png'",
        "set style data line",
        "set nologscale",
        "set xzeroaxis",
        "set xr [ 0.2100000E+02: 0.100000E+02]",
        "set yr [ -5: 80]",
        "set xlabel 'Energy[eV]'",
        f"set title '{model_name:s} Spectrum, n-cos = 1, tau: {tau:d}.0 1, {t_final:d}fs'",
        f"plot '{mctdh_file}.pl' using 1:3 lw 2 lc 'red' title 'MCTDH',\
         '{cc_file}.pl' every 6 using 1:3 with linespoints lc 'green' title 'CC',\
         '{sos_file}.pl' using 1:3 lc 'black' title 'SOS',\
          ",

        # f"plot '{mctdh_file}.pl' using 1:3 lw 4 lc 'red' title 'MCTDH',\
        #  '{x1}' every 6 using 1:3 with linespoints lc 'green' title 'CC halved',\
        #  '{x2}' every 4 using 1:3 with linespoints lc 'blue' title 'CC not halved',\
        #  '{x3}' every 6 using 1:3 lc 'black' title 'SOS halved',\
        #  '{x4}' every 4 using 1:3 lc 'orange' title 'SOS not halved',\
        # ",

        # f"plot '{mctdh_file}.pl' using 1:3 lw 4 lc 'red' title 'MCTDH',\
        #     '{x1}' every 6 using 1:3 with linespoints lc 'blue' title 'CC 12.5fs',\
        #     '{x2}' every 6 using 1:3 with linespoints lc 'magenta' title 'CC 25fs',\
        #     '{x3}' every 6 using 1:3 with linespoints lw 10 lc 'magenta' title 'CC 25fs (halved t)',\
        #     '{x4}' every 6 using 1:3 with linespoints lc 'green' title 'CC 50fs',\
        #     '{x5}' every 6 using 1:3 with linespoints lc 'orange' title 'CC 50fs (halved t)',\
        #     '{x6}' every 6 using 1:3 with linespoints lc 'yellow' title 'CC other',\
        # ",

        # f"plot '{mctdh_file}.pl' us 1:3 lc 'red' title 'MCTDH', '{cc_file}.pl' us 1:($3/6.582119564E-1) title 'CC' ",
        # f"plot '{cc_file}.pl' us 1:3 title 'CC' ",
    ])

    plotting_file = "spectrum_plotting.pl"
    # write the plotting commands to a file
    with open(plotting_file, 'w') as fp:
        fp.write(plotting_string)

    return plotting_file


def write_acf_plotting_file(nof_points, cc_file, mctdh_file, cc_file2):
    # plotting command
    plotting_string = '\n'.join([
        # "set style circle radius graph 0.002",
        "set style line 1 lt 2 pt 12 ps 1 pi -1",
        # "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5",
        # "set pointintervalbox 2",
        "set terminal png size 1200,800",
        f"set output './ACF_{model_name:s}_{nof_points:d}_{t_final:d}fs.png'",
        "set style data line",
        "set nologscale",
        "set xzeroaxis",
        "set yr [ -1: 1]",
        # f"set multiplot layout 2,2",
        "set title 'ACF comparison of {t_final:d}fs'",
        # top row (real)
        # "unset xtics", "unset xlabel",
        "set ylabel 'C(tau/hbar)'",
        "set xr [ 0.0: 100.0]",
        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 1 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 with linespoints ls 1 title 'CC Real' ,'{cc_file2}' us 1:2 with linespoints ls 1 title 'CC non interp Real' ",
        f"plot '{mctdh_file}' us 1:2 with linespoints ls 4 ps 3 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 lc 'green' title 'CC Real (interpolated)', '{cc_file2}' us 1:2 with linespoints ls 1 lc 'blue' title 'CC Real (raw RK45)' ",
        # #
        # "unset ytics", "unset ylabel",
        # "set xr [ 285.0: 300.0]",
        # #
        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 1 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 with linespoints ls 1 title 'CC Real' ",

        # # middle row (real)
        # "unset xtics", "unset xlabel",
        # "set ylabel 'C(tau/hbar)'",
        # "set xr [ 0.0: 15.0]",
        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 1 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 with linespoints ls 1 title 'CC Real' ",
        # #
        # "unset ytics", "unset ylabel",
        # "set xr [ 285.0: 300.0]",
        # #
        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 1 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 with linespoints ls 1 title 'CC Real' ",

        # # bottom row (imaginary)
        # "set ytics", "set xtics",
        # "set ylabel 'C(tau/hbar)'",
        # "set xlabel 'tau/hbar'",
        # "set xr [ 0.0: 15.0]",
        # #
        # f"plot '{mctdh_file}' us 1:3 with linespoints ls 1 lc 'red' title 'MCTDH Imag', '{cc_file}' us 1:3 with linespoints ls 1 title 'CC Imag' ",
        # #
        # "unset ytics", "unset ylabel",
        # "set xr [ 285.0: 300.0]",
        # f"plot '{mctdh_file}' us 1:3 with linespoints ls 1 lc 'red' title 'MCTDH Imag', '{cc_file}' us 1:3 with linespoints ls 1 title 'CC Imag' ",
        # #
        # "unset multiplot",
    ])

    plotting_file = "acf_plotting.pl"
    # write the plotting commands to a file
    with open(plotting_file, 'w') as fp:
        fp.write(plotting_string)

    return plotting_file


# a1 = f"ACF_CC_{model_name:s}_tf12_normalized_not_halved.txt"
# a2 = f"ACF_CC_{model_name:s}_tf25_normalized_not_halved.txt"
# a3 = f"ACF_CC_{model_name:s}_tf25_normalized_halved.txt"
# a4 = f"ACF_CC_{model_name:s}_tf50_normalized_not_halved.txt"
# a5 = f"ACF_CC_{model_name:s}_tf50_normalized_halved.txt"
# a6 = f"ACF_CC_{model_name:s}_tf{t_final}_normalized.txt"


def write_acf_sos_plotting_file(nof_points, cc_file, mctdh_file, sos_file):
    # plotting command
    plotting_string = '\n'.join([
        # "set style circle radius graph 0.002",
        "set style line 1 lt 2 lw 1.8 pt 12 ps 2.5 pi -1",
        "set style line 2 lt 1 lw 2.5 pt 1 ps 1 pi -1",
        # "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5",
        # "set pointintervalbox 2",
        "set terminal png size 1200,800",
        f"set output './ACF_{model_name:s}_{nof_points:d}_{t_final:d}fs.png'",
        "set style data line",
        "set nologscale",
        "set xzeroaxis",
        "set yr [ -1: 1]",
        # f"set multiplot layout 2,1",
        "set title 'ACF comparison of {t_final:d}fs'",
        # top row (real)
        # "unset xtics", "unset xlabel",
        "set ylabel 'C(tau/hbar)'",
        # "set key at graph 0.92, 1.08",
        # "set key outside",
        "set key opaque",
        f"set xr [ 0.0: {t_final*1.05:.2f}]",
        # f"set xr [ 0.0: {(t_final/2)*1.05:.2f}]",
        f"plot '{mctdh_file}' us 1:2 with linespoints ls 4 ps 3 lc 'red' title 'MCTDH Real',\
            '{cc_file}' using 1:2 with linespoints ls 1 lc 'green' title 'CC Real (interpolated)',\
            '{sos_file}' using 1:2 ls 2 lc 'black' title 'SOS Real',\
        ",
        # f"plot '{mctdh_file}' using 1:2 with linespoints ls 4 ps 3 lc 'red' title 'MCTDH Real',\
        #     '{a1}' every 6 using 1:2 with linespoints ls 1 lc 'green' title 'CC Real halved',\
        #     '{a2}' every 4 using 1:2 with linespoints ls 1 lc 'blue' title 'CC Real not halved',\
        #     '{a3}' every 6 using 1:2 ls 2 lc 'black' title 'SOS Real halved',\
        #     '{a4}' every 4 using 1:2 ls 2 lc 'orange' title 'SOS Real not halved',\
        # ",

        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 4 ps 3 lc 'red' title 'MCTDH Real',\
        #     '{a1}' using 1:2 with linespoints ls 1 lc 'blue' title 'CC Real 12.5fs',\
        #     '{a2}' using 1:2 with linespoints ls 2 lw 4 lc 'magenta' title 'CC Real 25fs',\
        #     '{a3}' using 1:2 with linespoints ls 4 lw 6 lc 'brown' title 'CC Real 25fs (halved t)',\
        #     '{a4}' using 1:2 with linespoints ls 3 lw 4 lc 'green' title 'CC Real 50fs',\
        #     '{a5}' using 1:2 with linespoints ls 1 lc 'black' title 'CC Real 50fs (halved t)',\
        #     '{a6}' using 1:2 with linespoints ls 1 lc 'yellow' title 'CC Real other',\
        # ",

        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 4 ps 3 lc 'red' title 'MCTDH Real',\
        #     '{cc_file}' using 1:2 with linespoints ls 1 lc 'green' title 'CC Real (interpolated)',\
        #     '{cc_file2}' using 1:2 with linespoints ls 1 lc 'blue' title 'CC Real (raw RK45)',\
        #     '{sos_file}' using 1:2 ls 2 lc 'black' title 'SOS Real',\
        # ",

        # f"plot '{mctdh_file}' us 1:3 with linespoints ls 4 ps 3 lc 'red' title 'MCTDH Imag',\
        #     '{cc_file}' using 1:3 with linespoints ls 1 lc 'green' title 'CC Imag (interpolated)',\
        #     '{sos_file}' using 1:3 ls 2 lc 'black' title 'SOS Imag',\
        # ",

        # #
        # "unset ytics", "unset ylabel",
        # "set xr [ 285.0: 300.0]",
        # #
        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 1 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 with linespoints ls 1 title 'CC Real' ",

        # # middle row (real)
        # "unset xtics", "unset xlabel",
        # "set ylabel 'C(tau/hbar)'",
        # "set xr [ 0.0: 15.0]",
        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 1 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 with linespoints ls 1 title 'CC Real' ",
        # #
        # "unset ytics", "unset ylabel",
        # "set xr [ 285.0: 300.0]",
        # #
        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 1 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 with linespoints ls 1 title 'CC Real' ",

        # # bottom row (imaginary)
        # "set ytics", "set xtics",
        # "set ylabel 'C(tau/hbar)'",
        # "set xlabel 'tau/hbar'",
        # "set xr [ 0.0: 15.0]",
        # #

        # #
        # "unset ytics", "unset ylabel",
        # "set xr [ 285.0: 300.0]",
        # f"plot '{mctdh_file}' us 1:3 with linespoints ls 1 lc 'red' title 'MCTDH Imag', '{cc_file}' us 1:3 with linespoints ls 1 title 'CC Imag' ",
        # #
        # "unset multiplot",
    ])

    plotting_file = "acf_plotting.pl"
    # write the plotting commands to a file
    with open(plotting_file, 'w') as fp:
        fp.write(plotting_string)

    return plotting_file


sos_spectrum_file = "sos_spectrum"
cc_spectrum_file = "cc_spectrum"
mctdh_spectrum_file = "mctdh_spectrum"

# generate data files
cc_acf_filename = f"ACF_CC_{model_name:s}_tf{t_final:d}.txt"
new_cc_acf = modify_acf_file(cc_acf_filename)
sos_acf_filename = f"ACF_SOS_{model_name:s}_{nof_BF}BF_tf{t_final:d}.txt"
new_sos_acf = modify_acf_file(sos_acf_filename)
#

nof_points = generate_cc_pl(nof_points, output_filename=cc_spectrum_file, input_filename=new_cc_acf)
generate_cc_pl(nof_points, output_filename=sos_spectrum_file, input_filename=new_sos_acf)
generate_mctdh_pl(nof_points, output_filename=mctdh_spectrum_file)


if True:  # plotting Spectrums
    plotting_file = write_spectrum_plotting_file(nof_points, cc_spectrum_file, mctdh_spectrum_file, sos_spectrum_file)
    os.system(f"gnuplot {plotting_file}")

if False:  # plotting ACF
    plotting_file = write_acf_plotting_file(nof_points, new_cc_acf, auto_path)
    os.system(f"gnuplot {plotting_file}")

if True:  # plotting ACF vs SOS
    plotting_file = write_acf_sos_plotting_file(nof_points, new_cc_acf, auto_path, new_sos_acf)
    os.system(f"gnuplot {plotting_file}")