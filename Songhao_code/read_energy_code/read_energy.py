# system import
import os
import sys
import math
import copy
import itertools as it
import re
import json


# third party import
import numpy as np
from log_conf import log
import pandas as pd


class read_energy(object):
    """
    The python class readin excited state energies (both adiabats and diabatis)
    from outputs of GAMESS diabatization calculations and store the data
    for plotting
    """
    def __init__(self, file_name, number_of_state, name):
        """intialize the python class
        file_name: common file name for the GAMESS output
        number_of_state: number of state included in this calculation
        name: name of the molecule
        """
        self.file_name = file_name
        self.number_of_state = number_of_state
        self.name = name
        #conversion constant
        self.ha2ev=27.2113961318
        self.wn2ev=0.000123981
        self.wn2eh=0.00000455633
        self.ang2br=1.889725989
        self.amu2me=1822.88839

    def _read_energy(self, name):
        """
        readin diabatic excited state energies
        """
        def read_adiabatic_energy(file):
            """read in  adiabatic energy"""
            adiabatic_energy = []

            # read data fpr adiabatc energy
            for idx,line in enumerate(file):
                if re.search("TZ Sorted GMC-PT ENERGIES:", line):
                    for i in range(self.number_of_state):
                        # print(output[idx+1+i])
                        tmp = output[idx+1+i].split()
                        adiabatic_energy.append(float(tmp[1]))
            return adiabatic_energy

        def read_diabatic_potential(file):
            """read in diabatic potential"""
            diagonal_energy = []
            off_diagonal_energy = []
            for idx, line in enumerate(file):
                if re.search("- DM DIABT PT HAMILTONIAN MATRIX ELEMENTS -", line):
                    # extract diagonal elements
                    for i in range(self.number_of_state):
                        # print(file[idx+4+i])
                        tmp = file[idx+4+i].split()
                        diagonal_energy.append(float(tmp[6]))
                    # extract off-diagonal elements
                    number_of_off_diagonal_elements = int((self.number_of_state**2 - self.number_of_state) / 2)
                    for i in range(number_of_off_diagonal_elements):
                        # print(file[idx+self.number_of_state+6+i])
                        tmp = file[idx+self.number_of_state+6+i].split()
                        off_diagonal_energy.append(float(tmp[8]))

            return diagonal_energy, off_diagonal_energy

        # open the output file
        f = open(name, "r")
        # store the output file in a list of lines
        output = f.readlines()
        # close the file
        f.close()

        # read in adiabatic state energy data from GAMESS output
        adiabatic_energy = read_adiabatic_energy(output)
        # print("adiabatic sate energy(in Hartree)\n{:}".format(adiabatic_energy))

        # read in diabatic potential data from GAMESS output
        diabatic_diagonal_energy, diabatic_off_diagonal_energy = read_diabatic_potential(output)
        # print("diagonal diabatic state energy(in Hartree):\n{:}".format(diabatic_diagonal_energy))
        # print("off-diagonal diabatic state energy(in Hartree):\n{:}".format(diabatic_off_diagonal_energy))

        return adiabatic_energy, diabatic_diagonal_energy, diabatic_off_diagonal_energy

    def scan_energy_1d(self, grid, omega):
        """
        scan a 1d grid of energies
        grid: a grid of displacement
        omega: frequency of the corresponding vibrational mode (in cm-1)
        """
        log.info("###Start 1d scanning of energy data from GAMESS output### ")
        # store energies in a python dictionary
        energy_data = {}
        # converge length to atomic unit
        grid_atm = np.sqrt(self.amu2me) * self.ang2br * np.sqrt(omega * self.wn2eh) * grid
        # store grid in atomic unit into the python dictionary
        energy_data["Grid"] = grid_atm
        print("Grid in atm unit:\n{:}".format(energy_data["Grid"]))

        # initialize keys of the dictionary
        ## store diabatic off_diagonal energy keys in are list of strings
        off_diagonal_diabatic_key_list = []
        for i in range(self.number_of_state):
            diagonal_diabatic_key = "diagonal_diabatic_state_{:d}".format(i+1)
            adiabatic_key = "adiabatic_state_{:d}".format(i+1)
            # initial energies as empty lists
            energy_data[adiabatic_key] = []
            energy_data[diagonal_diabatic_key] = []
            # keys for off diagonal elements
            for j in range(i, self.number_of_state):
                if i != j:
                    off_diagonal_diabatic_key = "off_diagonal_diabatic_state_{:d}_{:d}".format(i+1, j+1)
                    energy_data[off_diagonal_diabatic_key] = []
                    off_diagonal_diabatic_key_list.append(off_diagonal_diabatic_key)

        # loop over a grid of displacement
        for delta in grid:
            file_name = "{:}_{:.1f}.log".format(self.file_name, delta)
            print("At displacement: {:.1f}".format(delta))
            adiabatic_energy, diagonal_diabatic_energy, off_diagonal_diabatic_energy = self._read_energy(file_name)

            # store adiabatic state energy into the python dictionary
            for idx, energy in enumerate(adiabatic_energy):
                adiabatic_key = "adiabatic_state_{:d}".format(idx+1)
                energy_data[adiabatic_key].append(energy)

            # store diagonal diabatic state energy into the python dictionary
            for idx, energy in enumerate(diagonal_diabatic_energy):
                # store diagonal elements
                diagonal_diabatic_key = "diagonal_diabatic_state_{:d}".format(idx+1)
                energy_data[diagonal_diabatic_key].append(energy)

            # store off-diagonal diabatic state energy into the python dictionary
            for idx, key in enumerate(off_diagonal_diabatic_key_list):
                # store off-diagonal elements
                energy_data[key].append(off_diagonal_diabatic_energy[idx])

        # convert energy data into pandas DataFrame
        df = pd.DataFrame(energy_data)
        print(df.head())

        # store the data to json
        df.to_json("{:}_1d_scan_data.json".format(self.name))

        log.info("###1d scanning of energy data from GAMESS output succesfully stored into file!### ")

        return

    def scan_energy_2d(self, grid_1, grid_2, omega_1, omega_2):
        """
        scan a 2d grid of energies
        grid_1, grid_2: grids of displacement
        omega_1, omega_2: frequency of the corresponding vibrational mode (in cm-1)
        """
        log.info("###Start 2d scanning of energy data from GAMESS output### ")
        # store energies in a python dictionary
        energy_data = {}
        # initial the grids (2d coordinates) as empty lists
        energy_data["Grid_1"] = []
        energy_data["Grid_2"] = []

        # initialize keys of the dictionary
        ## store diabatic off_diagonal energy keys in are list of strings
        off_diagonal_diabatic_key_list = []
        for i in range(self.number_of_state):
            diagonal_diabatic_key = "diagonal_diabatic_state_{:d}".format(i+1)
            adiabatic_key = "adiabatic_state_{:d}".format(i+1)
            # initial energies as empty lists
            energy_data[adiabatic_key] = []
            energy_data[diagonal_diabatic_key] = []
            # keys for off diagonal elements
            for j in range(i, self.number_of_state):
                if i != j:
                    off_diagonal_diabatic_key = "off_diagonal_diabatic_state_{:d}_{:d}".format(i+1, j+1)
                    energy_data[off_diagonal_diabatic_key] = []
                    off_diagonal_diabatic_key_list.append(off_diagonal_diabatic_key)
        # loop over a grid of displacement
        for delta_1, delta_2 in it.product(grid_1, grid_2):
            store_flag = False
            file_name = "{:}{:.1f}_{:.1f}.log".format(self.file_name, delta_1, delta_2)
            # open the output file
            f = open(file_name, "r")
            # store the output file in a list of lines
            output = f.readlines()
            # close the file
            f.close()
            # if the output is not succefully run, skip it
            for line in output:
                if re.search("grace", line):
                    store_flag = True
            if not store_flag:
                log.info("Warning: the diabatzation procedure fail at displacement {:.1f}, {:.1f}".format(delta_1, delta_2))
            if store_flag:
                # store coordinate data (convert to the atomic unit)
                energy_data["Grid_1"].append(np.sqrt(self.amu2me) * self.ang2br * np.sqrt(omega_1 * self.wn2eh) * delta_1)
                energy_data["Grid_2"].append(np.sqrt(self.amu2me) * self.ang2br * np.sqrt(omega_2 * self.wn2eh) * delta_2)
                # print("At displacement: {:.1f}, {:.1f}".format(delta_1, delta_2))

                adiabatic_energy, diagonal_diabatic_energy, off_diagonal_diabatic_energy = self._read_energy(file_name)

                # store adiabatic state energy into the python dictionary
                for idx, energy in enumerate(adiabatic_energy):
                    adiabatic_key = "adiabatic_state_{:d}".format(idx+1)
                    energy_data[adiabatic_key].append(energy)

                # store diagonal diabatic state energy into the python dictionary
                for idx, energy in enumerate(diagonal_diabatic_energy):
                    # store diagonal elements
                    diagonal_diabatic_key = "diagonal_diabatic_state_{:d}".format(idx+1)
                    energy_data[diagonal_diabatic_key].append(energy)

                # store off-diagonal diabatic state energy into the python dictionary
                for idx, key in enumerate(off_diagonal_diabatic_key_list):
                    # store off-diagonal elements
                    energy_data[key].append(off_diagonal_diabatic_energy[idx])

        # convert energy data into pandas DataFrame
        df = pd.DataFrame(energy_data)
        print(df.head())

        # store the data to json
        df.to_json("{:}_2d_scan_data.json".format(self.name))

        log.info("###2d scanning of energy data from GAMESS output succesfully stored into file!### ")
        return
