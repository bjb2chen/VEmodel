#!python3

# system imports
import itertools as it
from shutil import copyfile
import os
from os.path import join

# third party imports

# local imports
import prop_input_template
from project_parameters import *


def create_input_file(source_dir, *args):
    '''Return formatted .inp file string.'''

    name, bf, tf = args

    nof_pbf = bf
    nof_spf = 1  # fix to 1 for the moment (this is tdh)
    # nof_spf = 4  # this is mctdh


    for operate_string in ["Ex", "Ey", "Ez"]:

        """
        the operator file name needs to be inside the *.inp file
        for some `system`
        we have three op files
            - system_Ex.op, system_Ey.op, system_Ez.op
        and three input files
            - system_Ex.inp, system_Ey.inp, system_Ez.inp
        """
        xyz_name = f"{project_name}_{operate_string}"

        input_string = "\n\n".join([
            prop_input_template.run_section_propagation.format(tfinal=tf/2, tout=0.5, name=name),
            prop_input_template.operator_section.format(opfile_name=f"{xyz_name:}"),
            prop_input_template.generate_basic_multi_set_spf_basis_section(nof_spf, N, A),
            prop_input_template.generate_basic_harmonic_oscillator_primative_basis_section(nof_pbf, N, A),
            prop_input_template.propagation_integrator_section.format(),
            prop_input_template.generate_basic_harmonic_oscillator_wavefunction_section(N, A, operate_string),

            "end-input\n"
        ])

        # the input file name
        inp_file_name = f"{xyz_name}.inp"

        # the path to the input file destination
        xyz_inp_path = join(source_dir, inp_file_name)

        # write the input file
        with open(xyz_inp_path, "w") as file:
            file.write(input_string)


def initalize_directories():
    '''Creates and populates all the directories where calculations will take place.'''
    
    for param_list in it.product(*expression_list):
        directory = join(work_root, dir_string.format(project_name, *param_list))
        os.makedirs(directory, exist_ok=True)
        os.makedirs(home_root, exist_ok=True)

        create_input_file(directory, project_name, *param_list)

        for operate_string in ["Ex", "Ey", "Ez"]:
            src_path_op = join(home_root, f"{project_name}_{operate_string}.op")
            dst_path_op = join(directory, f"{project_name}_{operate_string}.op")
            copyfile("mctdh.op", src_path_op)
            copyfile(src_path_op, dst_path_op)

            # sub directories
            xyz_subdirectory = join(directory, f"{operate_string}")
            os.makedirs(xyz_subdirectory, exist_ok=True)
            xyz_sub_path = join(xyz_subdirectory, f"{project_name}_{operate_string}.op")
            copyfile(src_path_op, xyz_sub_path)

            inp_file_name = f"{project_name}_{operate_string}.inp"
            xyz_inp_path = join(directory, inp_file_name)
            xyz_dst_path = join(xyz_subdirectory, inp_file_name)
            copyfile(xyz_inp_path, xyz_dst_path)

            xyz_sub_execution_script = join(xyz_subdirectory, execution_script)
            copyfile(src_path_execution_script, xyz_sub_execution_script)


        print(f"{directory = }")

        dst_path_execution_script = join(directory, execution_script)
        copyfile(src_path_execution_script, dst_path_execution_script)


if __name__ == "__main__":
    initalize_directories()
