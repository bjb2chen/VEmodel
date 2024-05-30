# system imports
import filecmp
import time
import os
from os.path import join

# third party imports
import numpy as np
# import matplotlib.pyplot as plt

# local imports
from .vibronic_hamiltonian import vibronic_hamiltonian
from project.vibronic import vIO, VMK


def model_2_full_example():

    root_directory = join(os.getcwd(), "vibronic_models")

    path_model = join(root_directory, "model_2.json")

    model = vIO.load_model_from_JSON(path_model)

    time_one = time.time()

    # initialize the model
    hamiltonian = vibronic_hamiltonian(model, trans=False, theta=np.pi/6)

    time_two = time.time()

    # hamiltonian.FCI_solution()

    # run coupled cluster propagation
    hamiltonian.leap_frog_integration(t_term=10, nof_points=10000)

    time_three = time.time()

    # make plot and store data
    plot_path = join(os.getcwd(), 'vibronic_models')
    t_1, C_tau_1 = hamiltonian.plot_acf(file_name="model_2", output_path=plot_path)

    assert filecmp.cmp(
        join(root_directory, 'model_2_ACF_CC.txt'),
        join(root_directory, 'model_2_comparison_CC_output.txt')
    ), "Files are not identical!"

    time_four = time.time()

    print(
        "\n",
        f"Time to build model:    {time_two - time_one:14.10f} seconds",
        f"Time to init leap frog: {time_three - time_two:14.10f} seconds",
        f"Time to plot model:     {time_four - time_three:14.10f} seconds",
        f"Total Time:             {time_four - time_one:14.10f} seconds",
        sep='\n'
    )
    return


def generate_acf_data(model, file_name, root_directory, order, t_final=20.0, nof_steps=10000):
    """ x """

    # initialize the model
    hamiltonian = vibronic_hamiltonian(model, highest_order=order)

    # run coupled cluster propagation
    # hamiltonian.leap_frog_integration(t_final=t_final, nof_points=nof_steps, highest_order=order)
    hamiltonian.rk45_integration(t_final=t_final, nof_points=nof_steps, highest_order=order)

    n = int(np.floor(np.log10(nof_steps)))
    r = nof_steps / pow(10, n)

    order_name = {
        0: "constant",
        1: "linear",
        2: "quadratic",
        3: "cubic",
        4: "quartic",
    }[order]
    # plot_name = f"{file_name}_{order_name}_n{r}e{n}_tf{int(t_final):}"
    plot_name = f"{file_name}_{order_name}_tf{int(t_final):}"

    # make plot and store data
    hamiltonian.plot_acf(file_name=plot_name, output_path=root_directory)
    hamiltonian.save_acf_data(file_name=plot_name, output_path=root_directory)

    return


def generate_acf_to_file(model, file_name, root_directory, order, nof_BF, t_final, nof_steps):
    """ generate ACF data and save to file for later use """

    # initialize the model
    hamiltonian = vibronic_hamiltonian(model, highest_order=order, build_H=False, HO_size=nof_BF)

    # compute ACF
    hamiltonian.rk45_integration(t_final=t_final, nof_points=nof_steps, highest_order=order)

    # save ACF results to file
    n = int(np.floor(np.log10(nof_steps)))
    r = nof_steps / pow(10, n)
    order_name = order_name_dict[order]
    file_name = f"{file_name}_{nof_BF}BF_{order_name}_tf{int(t_final)}"
    hamiltonian.save_acf_data(file_name, output_path=root_directory)
    return


def generate_sos_to_file(model, file_name, root_directory, order, nof_BF, t_final, nof_steps):
    """ generate SOS data and save to file for later use """

    # initialize the model
    hamiltonian = vibronic_hamiltonian(model, highest_order=order, build_H=True, HO_size=nof_BF)

    # compute SOS
    hamiltonian.FCI_solution(t_final=t_final, nof_steps=nof_steps)

    # save SOS results to file
    n = int(np.floor(np.log10(nof_steps)))
    r = nof_steps / pow(10, n)
    # order_name = order_name_dict[order]
    file_name = f"{file_name}_{nof_BF}BF_tf{int(t_final)}"
    hamiltonian.save_sos_data(file_name, output_path=root_directory)
    return


def generate_smoothed_acf_to_file(model, file_name, root_directory, order, t_final, nof_steps):
    """ generate ACF data using rk45, smooth the data and save to file for later use """

    # initialize the model
    hamiltonian = vibronic_hamiltonian(model, highest_order=order)

    # compute ACF
    hamiltonian.rk45_integration(t_final=t_final, nof_points=nof_steps, highest_order=order)

    # save ACF results to file
    n = int(np.floor(np.log10(nof_steps)))
    r = nof_steps / pow(10, n)

    # save ACF results to file
    output_filename = f"{file_name}_tf{int(t_final):>03d}"
    path = hamiltonian.save_acf_data(output_filename, output_path=root_directory)
    # hamiltonian.plot_acf(file_name=output_filename, output_path=root_directory)

    if False:
        # smooth the data
        hamiltonian.smooth_acf_using_interpolation(nof_points=nof_steps)

        # save ACF results to file
        interpolated_filename = f"{file_name}_tf{int(t_final)}_interpolated"
        path = hamiltonian.save_acf_data(interpolated_filename, output_path=root_directory)
        # hamiltonian.plot_acf(file_name=interpolated_filename, output_path=root_directory)

    return


def generate_acf_sos_data(model, file_name, root_directory, order, nof_BF, t_final, nof_steps=10000):
    """ x """

    # initialize the model
    hamiltonian = vibronic_hamiltonian(model, highest_order=order, build_H=True, HO_size=nof_BF)

    # hamiltonian.load_acf_data(join(root_directory, 'ACF_CC_model_1_SOS10_quadratic_tf10.txt'))

    # check with SOS code
    hamiltonian.FCI_solution(t_final=t_final, nof_steps=nof_steps)

    # run coupled cluster propagation
    # hamiltonian.leap_frog_integration(t_final=t_final, nof_points=nof_steps, highest_order=order)
    hamiltonian.rk45_integration(t_final=t_final, nof_points=nof_steps, highest_order=order)

    n = int(np.floor(np.log10(nof_steps)))
    r = nof_steps / pow(10, n)
    order_name = order_name_dict[order]

    # plot_name = f"{file_name}_{nof_BF}BF_{order_name}_n{r}e{n}_tf{int(t_final):}"
    plot_name = f"{file_name}_{nof_BF}BF_{order_name}_tf{int(t_final):}"

    # make plot and store data
    hamiltonian.plot_acf_and_save_data(
        file_name=plot_name, output_path=root_directory, sos_flag=True
    )

    # simple average smoothing
    if False:
        x, y = hamiltonian.t_cc, hamiltonian.C_tau_cc
        print(x[0:7])

        differences = abs(x[0:-2] - x[1:-1])
        avg = np.average(differences)
        time_step = round(avg, ndigits=8)
        print(avg, time_step)

        smoothed_x = new_x = np.arange(0.0, stop=x[-1], step=time_step)
        smoothed_y = y
        print(len(smoothed_x), len(x))
        print(smoothed_x[0:5], '\n', smoothed_x[-5:], '\n', x[-1])

        hamiltonian.plot_acf(
            file_name=plot_name+'_smoothed', output_path=root_directory, sos_flag=True,
            acf_x_array=smoothed_x,
        )

    # interpolation
    if False:
        from scipy.interpolate import interp1d
        x, y = hamiltonian.t_cc, hamiltonian.C_tau_cc
        f = interp1d(x, y)
        f2 = interp1d(x, y, kind='quadratic')
        f3 = interp1d(x, y, kind='cubic')

        smoothed_x = np.linspace(0.0, t_final, num=nof_steps, endpoint=True)
        xnew = np.linspace(0, 10, num=41, endpoint=True)
        import matplotlib.pyplot as plt
        plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--', xnew, f3(xnew), '--')
        plt.legend(['data', 'linear', 'quadratic', 'cubic'], loc='best')

        # smoothed_y =
        # print(len(smoothed_x), len(x))
        # print(smoothed_x[0:5], '\n', smoothed_x[-5:], '\n', x[-1])

        # hamiltonian.plot_acf(
        #     file_name=plot_name+'_interpolated', output_path=root_directory, sos_flag=True,
        #     acf_x_array=smoothed_x,
        # )

    import sys
    sys.exit()

    path = hamiltonian.save_acf_data(
        "ACF_smoothed", time=smoothed_x, acf=smoothed_y, output_path=root_directory
    )

    return


def main():
    model_2_full_example()


if (__name__ == '__main__'):
    main()
