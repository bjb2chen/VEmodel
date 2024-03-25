# system imports
import io
import time
import os
from os.path import abspath, join, dirname, basename
import sys
import cProfile
import pstats


# third party imports
import numpy as np
import runpy as rp
import matplotlib.pyplot as plt
import pandas as pd

# import the path to the package
#project_dir = abspath(join(dirname(__file__), '/home/bsonghao/t-amplitudes/'))
project_dir = abspath(join(dirname(__file__), '/home/bjb2chen/VECC/'))
sys.path.insert(0, project_dir)

import project

# local imports
from project.vibronic_hamiltonian import vibronic_hamiltonian
from project.vibronic import vIO, VMK, model_op, vIO_wrapper
import project.spectra
import project.log_conf


order_dict = {
    0: "constant",
    1: "linear",
    2: "quadratic",
    3: "cubic",
    4: "quartic",
}

root_directory = os.getcwd()
# file_name = "3st_water_mar19"
#file_name = "CoF3"
#file_name = "H2Ocat_ccd_gmcpt_C1_3st_diab"
#file_name = "mar19"
#file_name = "NH3_mar22"
#file_name = "water_xyz_uniquetdm"
#file_name = "RhF3_SOC_15st"
file_name = "NH3_XONLY"


order = 1
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

        with open(filename+".txt", "w+") as f:
        #    s = io.StringIO()
        #    p = pstats.Stats(filename)
            # p.strip_dirs().sort_stats("tottime").print_stats(2)
        #    p.strip_dirs().sort_stats("tottime").print_stats(6)
        #    p.strip_dirs().sort_stats("cumulative").print_stats(20)
        #   # p.strip_dirs().sort_stats("cumulative").print_stats('calculate', 15)
        #   # p.strip_dirs().sort_stats("cumulative").print_callees('_calculate_order_4_w_operator_optimized')
            f.write(s.getvalue())
    except Exception:
        print("cProfile data is not stored properly!")

    def print_based_on_criteria():
        # just example code, isn't called anywhere
        spec = "{title:^90}\n" + f"{'':*>90s}" + "\n{data}"
        string_list = []

        for title, criteria in [
            ("Highest # of calls", "calls"),
            ("Highest cumulative time", "cumulative"),
            # ("Highest percall time", "pcalls"),
            ("Highest tottime", "tottime")
        ]:
            s = io.StringIO()
            p = pstats.Stats(filename, stream=s)
            p.strip_dirs().sort_stats(criteria).print_stats(10)
            print(spec.format(
                title=title,
                data="\n".join(s.getvalue().splitlines()[7:])
            ))

    return


def generate_profile_data(filename):
    cProfile.run('project.test_vibronic_new.main()', filename)
    nof_BF += 1
    return


def generate_acf_data(model, file_name, order, t_final=10.0, nof_steps=int(1e4) , FC=False, compare_FCI=False):
    """ x """
    nof_BF = 40
    hamiltonian = vibronic_hamiltonian(
        model, file_name, HO_size=nof_BF, build_H=compare_FCI,
        cc_truncation_order=order, hamiltonian_truncation_order=order, FC=FC,
        Z_truncation_order=3, T_truncation_order=1, selected_surface=[]
        )# We truncation_order: flag to determine order of truncation in W amplitude
         # solve_W: flag to determine whether to include W in the ansatz
         # solve_Z: flag to determine whether to do the second similarity transformation

    # # check with SOS code
    if compare_FCI:

        hamiltonian.FCI_solution(t_final=t_final, nof_steps=nof_steps)

        sos_name = f"{file_name}_{nof_BF}BF_tf{int(t_final)}"

        hamiltonian.save_sos_data(sos_name, output_path=root_directory)

    # run coupled cluster propagation
    if True:
        filename = f"cProfile_{file_name}_{order_dict[order]}_tf{int(t_final):}"

        filename += "_ST_VECC_Z_{:}_T_{:}_22_06_24".format(hamiltonian.Z_truncation_order,\
                                                      hamiltonian.T_truncation_order)

        func_string = 'hamiltonian.rk45_integration(t_final=t_final, nof_points=nof_steps)'

       # eval(func_string)

        if True:
            cProfile.runctx(
                    func_string,
                     globals(),
                     locals(),
                     filename)


        process_data(filename)

        # plot norm of the wavefunction
       # Plot_Norm(hamiltonian.A, hamiltonian.t_cc, hamiltonian.Norm_cc)

        # plot energy expectation value
        # Plot_Energy(hamiltonian.A, hamiltonian.t_cc, hamiltonian.Energy_cc)

        # plot state population
       # Plot_State_Pop(hamiltonian.A, hamiltonian.t_cc, hamiltonian.state_pop_cc_DB, string = "Diabatic")
       # Plot_State_Pop(hamiltonian.A, hamiltonian.t_cc, hamiltonian.state_pop_cc_AB, string = "Adiabatic")

        # sys.exit(0)

    else:
        start_time = time.time()
        # hamiltonian.leap_frog_integration(t_final=t_final, nof_points=nof_steps)
        # hamiltonian.rk45_integration(ansatz="dU", t_final=t_final, nof_points=nof_steps)
        hamiltonian.rk45_integration(ansatz="dS", t_final=t_final, nof_points=nof_steps)
        # hamiltonian.rk45_integration(ansatz="new_1", t_final=t_final, nof_points=nof_steps)
        # hamiltonian.rk45_integration(ansatz="new_2", t_final=t_final, nof_points=nof_steps)
        # hamiltonian.rk45_integration(ansatz="new_3", t_final=t_final, nof_points=nof_steps)
        end_time = time.time()
        print(f"Very simple timing:   {end_time - start_time:10.4f}")
        # print(f"U_0 base took         {0.0331:10.4f}   for 50fs")
        # print(f"S_0 base took         {1769.2978:10.4f}   for 50fs with SD")
        # print(f"ansatz 1 took         {233.23:10.4f}   for 50fs with SD")
        # print(f"ansatz 1 took         {2261.2861:10.4f}   for 50fs with SDT")
        # print(f"ansatz 2 took         {57.7581:10.6f}   for 50fs with SD")
        # print(f"ansatz 3 took         {0:10.6f}   for 50fs with SD")


        # ansatz 2 failed

    plot_name = f"{file_name}_{order_dict[order]}_tf{int(t_final):}"
    try:
        if int_type == "new_4" and hamiltonian.solve_W:
            plot_name += "_W_order_{}".format(order_dict[hamiltonian.W_truncation_order])
            print(plot_name)
            if hamiltonian.solve_Z:
                plot_name += "_Z_on"
            else:
                plot_name += "_Z_off"
        elif int_type == "new_5":
            plot_name += "_ST_VECC_Z_{:}_T_{:}_EF_opt_state_cross_off".format(hamiltonian.Z_truncation_order,\
                                                      hamiltonian.T_truncation_order)
    except Exception:
        print("There is an error in name !")

    # make plot and store data
    output_path = hamiltonian.save_acf_data(file_name=plot_name, output_path=root_directory)
    # hamiltonian.plot_acf(file_name=plot_name, output_path=root_directory)

    return output_path


def gnuplot_spectrum(*args):
    """ a """

    # unpack arguments
    exec_name, normalized_path_ABS, model_name, order, t, FC = args


    left_eV, right_eV = 21.0, 8.5
    min_y, max_y = -3, 60
    nof_points = 5000
    tau = 40
    iexp = 1

    fc_string = "FC" if FC else "vibronic"

    # fourier transform
    command = (
        "#autospec84 "
        # "-g 1 "  # to print gnuplot commands or not
        f"-o {exec_name:s}.pl "
        f"-f {normalized_path_ABS:s} "
        f"-p {nof_points:d} "
        # f"-EP "
        # f"-e {harmonic_ground_state} eV " # x axis shift (left/right) to account for H.O. G.S.
        f"{left_eV} {right_eV} eV "  # x axis limits (in eV)
        f"{tau:d} "   # tau value
        f"{iexp:d} "  # iexp value
    )
    #     # "output_filename": join(root_directory, exec_name + ".pl"),
    #     "output_filename": "./" + exec_name + ".pl",
    #     "input_filename": "./" + normalized_path,
    #     "nof_points": nof_points,
    #     "left_eV": left_eV,
    #     "right_eV": right_eV,
    #     "tau": tau,
    #     "iexp": iexp,
    #     # 'g3': True
    # })
    #os.system(command)   # execute autospec84

    # same thing for MCTDH
    # command = project.spectra.generate_mctdh_pl(**{
    #     # "output_filename": join(root_directory, exec_name + ".pl"),
    #     "output_filename": f"./MCTDH_{fc_string}_{order_dict[order]}.pl",
    #     "input_filename": f"./auto_{fc_string}_{order_dict[order]}",
    #     # "output_filename": f"./SOS.pl",
    #     # "input_filename": f"./ACF_SOS_h2o_FC_40BF_tf50.txt",
    #     "nof_points": nof_points,
    #     "left_eV": left_eV,
    #     "right_eV": right_eV,
    #     "tau": tau,
    #     "iexp": iexp,
    # })
    print(command)
    #os.system(command)   # execute autospec84

    # doctor the file name to make it look better in the plot
    plot_title = model_name.replace('_', ' ').replace('h2o', 'h_{2}o')

    size = [1200, 800]

    # output_file = f'{root_directory}/cc_spectrum_{model_name:s}_{nof_points:d}_{tau:d}tau.png'
    output_file = f'{root_directory}/cc_spectrum_{model_name:s}.png'

    plotting_command = '\n'.join([
        f"set terminal png size {size[0]},{size[1]}",
        # f"set output './spectrum_{model_name:s}_{nof_points:d}_{t_final:d}fs_{tau:d}tau_{nof_BF}SOSBF_{mctdh_BF}MCTDHBF.png'",
        f"set output '{output_file:s}'",
        # "set style data line",
         "set nologscale", "set xzeroaxis", "set xlabel 'Energy[eV]'",
        f"set xr [ {left_eV}: {right_eV}]",
        f"set yr [ {min_y}: {max_y}]",
        f"set title '{plot_title:s} spectrum, n-cos = 1, tau: {tau:d}.0 1, {int(t):d}fs'",
        # "set style line 1 lt 2 lw 4 lc 'black'",
        # "set style line 2 lt 3 lw 2 lc 'red'",
        f"#plot '{root_directory}/op_NH36Q_3st_PBF30_tf250.00_auto_total' using 1:2 with lines ls 1 lc 'black' title 'MCTDH ({fc_string})'",
        f"plot \
            '{root_directory}/{exec_name}.pl' using 1:3 with lines ls 1 lc 'red' title 'CC',\
        ",
            # '{root_directory}/{exec_name}_6.pl' every 6 using 1:3 linetype 4 title '1E-6 1E-9 CC',\
            # '{root_directory}/{exec_name}_9.pl' every 6 using 1:3 linetype 3 title '1E-9 1E-12 CC',\
            # '{root_directory}/{exec_name}_12.pl' every 6 using 1:3 linetype 1 title '1E-12 1E-15 CC',\
    ])


    path_plotting_file = f"{root_directory}/spectrum_plotting.pl"

    # write the plotting commands to a file
    with open(path_plotting_file, 'w') as fp:
        fp.write(command)
        fp.writte("\n")
        fp.write(plotting_command)

    return path_plotting_file


def write_acf_plotting_file(cc_file, cc_file2, FC, t_final, nof_points=5000):
    # plotting command

    fc_string = "FC" if FC else "vibronic"

    plotting_string = '\n'.join([
        # "set style circle radius graph 0.002",
        "set style line 1 lt 2 pt 12 ps 1 pi -1",
        # "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5",
        # "set pointintervalbox 2",
        "set terminal png size 1200,800",
        f"set output './ACF_{model_name:s}_{order_dict[order]}_{t_final:.0f}fs.png'",
        "set style data line",
        "set nologscale",
        "set xzeroaxis",
        # "set yr [ -3: 3]",
        # f"set multiplot layout 2,2",
        f"set title 'ACF comparison of {order_dict[order]} {model_name:s} {t_final:.0f}fs'",
        # top row (real)
        # "unset xtics", "unset xlabel",
        "set ylabel 'C(tau/hbar)'",
        f"set xr [ 0.0: {t_final:d}.0]",
        # f"plot '{mctdh_file}' us 1:2 with linespoints ls 1 lc 'red' title 'MCTDH Real', '{cc_file}' us 1:2 with linespoints ls 1 title 'CC Real' ,'{cc_file2}' us 1:2 with linespoints ls 1 title 'CC non interp Real' ",
        f"plot \
            '{cc_file}' us 1:2 lc 'green' title 'CC Real (raw RK45)',\
            'ACF_SOS_{file_name}_FC_40BF_tf50.txt' us 1:2 lc 'red' title 'SOS Real',\
        "
            # '{cc_file2}' us 1:2 with linespoints ls 1 lc 'blue' title 'CC Real (interpolated)',\
            # '{root_directory}/auto_{fc_string}_{order_dict[order]}' us 1:2 with linespoints ls 4 ps 3 lc 'red' title 'MCTDH Real ({fc_string})',\
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

if (__name__ == '__main__'):

    t_final = 250.0
    FC = False
    model_name = f"{file_name}_FC" if FC else f"{file_name}_vibronic"


    project.log_conf.setLevelDebug()
    # ----------------------------------------------------------------
    print("We are running calculation for {:} model!".format(file_name))

    # read in model parameters 
    new_model = vIO_wrapper.vibronic_input_reader(name=file_name, hamiltonian_truncation_order=order, path="/home/bjb2chen/VECC/vibronic_models/")

    # run CC code
    output_path_ABS, output_path_ECD = generate_acf_data(new_model, model_name, order, t_final,FC=FC, compare_FCI=False)

    # interpolate for ACF(ABS)
    print("-"*40 + "\nInterpolating ABS\n" + "-"*40 + "\n")
    normalized_path_ABS = project.spectra.generate_normalized_acf_results(
            dirname(output_path_ABS),
            basename(output_path_ABS),
            None,
            mctdh_t_final = t_final*0.5,
            mctdh_dt = 0.1
        )
    # interpolate for ACF(ECD)
    print("-"*40 + "\nInterpolating ECD\n" + "-"*40 + "\n")
    normalized_path_ECD = project.spectra.generate_normalized_acf_results(
            dirname(output_path_ECD),
            basename(output_path_ECD),
            None,
            mctdh_t_final = t_final*0.5,
            mctdh_dt = 0.1
        )
    print('Do you want to plot the ABS spectra? Press c to continue.')
    breakpoint()
        # # plot
    print("-"*40 + "\nPlotting Spectrum\n" + "-"*40 + "\n")
    gnuplot_spectrum(
         f"{model_name}_{order_dict[order]}_tf{int(t_final):}",
         basename(normalized_path_ABS),
         f"{model_name}_{order_dict[order]}",
         order,
         t_final,
         FC,
     )

    #os.system(f"gnuplot spectrum_plotting.pl")


    # write_acf_plotting_file(output_path_ABS, normalized_path, FC, int(t_final))
    # os.system(f"gnuplot acf_plotting.pl")
