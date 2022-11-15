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

class dist_structure(object):

    """
    The python class aims to replace the bash script written by TZ
    The script takes in four files:

    ref_structure: The MP2 optimized structure (gain from step 1)

    diab_info: eigenvectors ($VEC), diabatic molecular orbitals
     ($DMOs and adiabats at the reference structure ($REFDET)
     (notice that the $END of DMO groups must align with the $REFDET otherwise the GAMESS will not do diabatization)

    temp.inp: This file provide a template initial settings of GAMESS GMCPT calculation.
    For different molecules, one need to manually change the active space and norb parameters.

    *.log: output file for GAMESS mp2 geometry optimization and frequency calculation

    and produce the GAMESS input file that performs diabatization calculation
    for a grid of displaced structure along each normal mode
    """
    def __init__(self, file_name, molecule_name):
        """define object instance
        file_name: name of GAMESS output file
        molecule_name: mane of the molecule
        self.linear_flag: a flag to determine if the molecule is linear
        """
        self.file_name = file_name
        self.molecule_name = molecule_name
        self.linear_flag = False

    def read_vibration(self):
        """read in vibration data (number of atoms, frequency, normal mode vector) from GAMESS output"""

        def extract_frequancy():
            """extract frequency from the GAMESS output data"""
            # initialize a list of frequancies
            freq_list = []
            # extract frequancy data form the GAMESS output file
            for line in self.GAMESS_output_list:
                # extract frequancies by searching the keyword in each line
                if re.search("FREQUENCY: ", line):
                    tmp = line.split()
                    for freq in tmp:
                        if freq != "FREQUENCY:":
                            freq_list.append(float(freq))
            return freq_list

        def extract_number_of_atom():
            """extract number of atoms from the GAMESS output data"""
            for line in self.GAMESS_output_list:
                # extract number of atom by searching the keyword in each line
                if re.search("TOTAL NUMBER OF ATOMS", line):
                    tmp = line.split()
                    number_of_atom = int(tmp[-1])

            return number_of_atom

        def extract_mode_vector():
            """extract normal mode vector from GAMESS output file"""
            def get_vector(line, iter_idx, coord_num):
                """
                multipulate with the normal mode vector and
                store them in numpy array to do linear algebra
                """
                starting_mode = iter_idx * 5
                for idx, member in enumerate(line):
                    if re.search("X", member) or re.search("Y", member) or re.search("Z", member):
                        for i in range(idx+1, len(line)):
                            mode_num = i - idx -1 + starting_mode
                            mode_vector[coord_num, mode_num] = float(line[i])
                return


            # initialize the matrix of normal mode vectors
            # into a 3*N by 3*N numpy array
            mode_vector = np.zeros([3*self.number_of_atom, 3*self.number_of_atom])
            iter_idx = 0

            for idx, line in enumerate(self.GAMESS_output_list):
                # search for the line that contain the keyword "IR INTENSITY:"
                if re.search("IR INTENSITY:", line):
                    # print("LINE NUMBER:{:}".format(idx))
                    # the line number between this idx + 2 and idx +2 + 3 * N
                    # stores the information for norm mode vector
                    coord_num = 0
                    for i in range(idx+2, idx + 2 + self.number_of_atom * 3):
                        # print(self.GAMESS_output_list[i])
                        line = self.GAMESS_output_list[i].split()
                        #print(line)
                        get_vector(line, iter_idx, coord_num)
                        coord_num += 1
                    iter_idx += 1

            #print(mode_vector)

            return mode_vector

        log.info("##### starting extracting vibration data from GAMESS output #####")
        # read in GAMESS output file of geometry optimization
        GAMESS_output = open(self.file_name, "r")
        # store the output in a list of string corresponding to each line
        GAMESS_output_list = GAMESS_output.readlines()
        # close the file
        GAMESS_output.close()
        # store the output file as an object instance
        self.GAMESS_output_list = GAMESS_output_list

        # extract frequcency data from the GAMESS output
        self.freq_list = extract_frequancy()
        print("Frequencies:\n {:}".format(self.freq_list))

        # extract totol number of atoms
        self.number_of_atom = extract_number_of_atom()
        print("Total number of atom: {:}".format(self.number_of_atom))

        # test if number of frequancies equals to 3 * N
        assert len(self.freq_list) == 3 * self.number_of_atom

        # extract normal mode vector from the GAMESS output
        self.normal_mode_vector = extract_mode_vector()
        print("Normal mode vector \n {:}".format(self.normal_mode_vector))

        log.info("#### successfully extracting vibration info from GAMESS output ####")

        return

    def _dist_structure(self, ref_file, displace):
        """generate GAMESS input give certain distorted structure

        ref_file: a file that stores the reference structure before the displacement
        displace: a vector of displacement

        new_structure = ref_structrue + displacement
        """
        def store_structure_in_vector(structure):
            """store the reference structure into a vector
             with 3*N dimension"""
            # string atom names in to a string of text
            atom_list = []
            # string atom charges in to a stong of text
            charge_list = []
            # initial the vector that store the structure
            vec = np.zeros(3 * self.number_of_atom)
            # store the elements of the vector
            for idx, line in enumerate(structure):
                tmp = line.split()
                atom_list.append(tmp[0])
                charge_list.append(tmp[1])
                for i in range(3):
                    vec[3*idx+i] = float(tmp[i+2])
            return vec, atom_list, charge_list

        def convert_vector_to_text(structure_vec, atom_list, charge_list):
            """convert the structure that stored in vector to string"""
            # store the structure in a list of string
            structure = []
            for idx, atom in enumerate(atom_list):
                tmp = []
                tmp.append(atom)
                tmp.append(charge_list[idx])
                for i in range(3):
                    tmp.append(str(structure_vec[idx*3+i]))
                tmp = " ".join(tmp)
                tmp += "\n"
                structure.append(tmp)
            return structure

        # read in reference structure from the file
        tmp = open(ref_file, "r")
        reference_structure = tmp.readlines()
        tmp.close()

        # store the reference structure into a 1d Array
        ref_structure_vec, atom_list, charge_list = store_structure_in_vector(reference_structure)
        # print("Reference structure presented in a 1d array:\n{:}".format(ref_structure_vec))

        # creat the displace structure ane presented in 1 1d Array
        new_structure_vec = ref_structure_vec + displace
        # print("Distorted structure presented in a 1d array:\n{:}".format(new_structure_vec))
        # print("list of atoms in the molecule:\n{:}".format(atom_list))

        # convert the displaced structure into a list of strings
        new_structure = convert_vector_to_text(new_structure_vec, atom_list, charge_list)
        # print("displace structure represented in a list of strings")
        # for line in new_structure:
            # print(line)

        return new_structure

    def _generate_GAMESS_input(self, structure, temp_file, diab_info, file_name):
        """generate GAMESS input files for diabatization
        structure: a list of strings stores the molecular structure
        temp_file: GAMESS template input file
        diab_info: a file that stores $VEC, $DMO and $REFDET groups
        """
        # open the template file
        f = open(temp_file, "r")
        # store template in a list of lines
        template = f.readlines()
        # close the template file
        f.close()
        # print("### template file")
        # for line in template:
            # print(line)

        # open diab_info file
        f = open(diab_info, "r")
        # store the diab_info file in a list of lines
        diabatization = f.readlines()
        # close the diab_info file
        f.close()
        # print("### diab_info file")
        # for line in diabatization:
            # print(line)

        # store distorted structure into a new file
        f = open(file_name, "w")
        # write template file
        f.writelines(template)

        # write structure info
        # structure = self._dist_structure(structure, np.zeros(9)) # this line is only needed for testing purpose
        f.writelines(structure)

        # write the diab_info
        f.writelines(diabatization)

        # close the file
        f.close()

        return

    def _write_bash_script_for_calculation(self, input_name, file_name="run_script_on_head_node"):
        """write the generated GAMESS input into bash script for calculations of diabatization
        input_name: name of the GAMESS input file
        file_name: name for the bash script
        """
        # open the bash file
        f = open(file_name, "a")

        # define a string to be added
        string = "/home/bjb2chen/gamess/rungms "
        string += input_name
        string += " 00 1 &> "
        string += input_name.replace(".inp", ".log")
        string += "\n"
        string += "sleep 10\n"

        # write the string into the bash file
        f.write(string)

        # close the bash file
        f.close()
        return

    def dis_structure_1d(self, ref_file, grid, temp_file, diab_info):
        """scan a 1d grid over each vibrational mode to systematically generate
        distorted structures"""
        log.info("#### Start 1d scan of a grid of displacement structure ####")
        # skip translational / rotational mode in the scan
        if self.linear_flag:
            skip = 5
        else:
            skip = 6
        # loop over each normal mode
        for mode in range(3*self.number_of_atom):
            # loop over a grid of displacement
            if mode >= skip:
                for delta in grid:
                    # define displacement vector
                    dis_vec = delta * self.normal_mode_vector[:, mode]
                    structure = self._dist_structure(ref_file, dis_vec)
                    # print("## mode: {:} displace: {:}".format(mode, delta))
                    # for line in structure:
                        # print(line)
                    name = "{:}_diabatization_1d_at_mode{:d}_dist{:.1f}.inp".format(self.molecule_name, mode, delta)
                    # write the get generated input file into bash file for calculation
                    self._write_bash_script_for_calculation(input_name=name)
                    # generate GAMESS input for diabatization
                    self._generate_GAMESS_input(structure, temp_file, diab_info, name)
        log.info("### 1d scan successfully conducted ! ###")
        return

    def dis_structure_2d(self, ref_file, grid, temp_file, diab_info):
        """scan  a 2d grid over each vibrational mode to systmatically generate
        distorted structures"""
        log.info("#### Start 2d scan of a grid of displacement structure ####")
        # skip translational / rotational mode in the scan
        if self.linear_flag:
            skip = 5
        else:
            skip = 6
        # loop over each normal mode
        for mode_1, mode_2 in it.product(range(3*self.number_of_atom), repeat=2):
            # loop over a grid of displacement
            if mode_1 >= skip and mode_2 >= skip:
                for delta_1, delta_2 in it.product(grid, repeat=2):
                    # define displacement vector
                    dis_vec = delta_1 * self.normal_mode_vector[:, mode_1] + delta_2 * self.normal_mode_vector[:, mode_2]
                    structure = self._dist_structure(ref_file, dis_vec)
                    name = "{:}_diabatization_2d_at_mode_{:d}_{:d}_dist{:.1f}_{:.1f}.inp".format(self.molecule_name, mode_1, mode_2, delta_1, delta_2)
                    # write the get generated input file into bash file for calculation
                    self._write_bash_script_for_calculation(input_name=name)
                    # generate GAMESS input for diabatization
                    self._generate_GAMESS_input(structure, temp_file, diab_info, name)
        log.info("### 2d scan successfully conducted ! ###")

        return
