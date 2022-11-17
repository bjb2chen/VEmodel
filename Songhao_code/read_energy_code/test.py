from read_energy import *
import numpy as np


def main():
    file_name = "nh3_diabatization_2d_at_mode_7_8_dist"
    name = "nh3"
    number_of_state = 3
    data = read_energy(file_name, number_of_state, name)
    omega1 = 1687.8
    omega2 = 1687.8
    grid1 = np.linspace(-0.8, 0.8, 17)
    grid2 = np.linspace(-0.8, 0.8, 17)
    # read in adiabatic energy data
    # data.read_adiabatic_energy(name="nh3_diabatization_2d_at_mode_8_9_dist0.0_0.0.log")
    # read in diabatic energy data
    # data._read_energy(name="nh3_diabatization_2d_at_mode_8_9_dist0.0_0.0.log")
    # read in a 1d grid of energy data
    # data.scan_energy_1d(grid=np.linspace(-0.8, 0.8, 17), omega=omega1)
    # read in a 2d grid of energy data
    data.scan_energy_2d(grid_1=grid1, grid_2=grid2, omega_1=omega1, omega_2=omega2)



    return


if __name__ == "__main__":
    main()
