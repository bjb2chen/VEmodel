#!python3

# system imports
import itertools as it
from shutil import copyfile
import os
from os.path import join

# third party imports

# local imports
from vibronic import vIO, VMK
import input_template
from project_parameters import *


def create_input_file(path, *args):
    '''Return formatted .inp file string.'''

    name, bf, tf = args

    nof_pbf = bf
    nof_spf = 1  # fix to 1 for the moment

    input_string = "\n\n".join([
        input_template.run_section_propagation.format(tfinal=tf/2, tout=0.1, name=name),
        input_template.operator_section.format(opfile_name=f"{name:}"),
        input_template.generate_basic_multi_set_spf_basis_section(nof_spf, N, A),
        input_template.generate_basic_harmonic_oscillator_primative_basis_section(nof_pbf, N, A),
        input_template.propagation_integrator_section.format(),
        input_template.generate_basic_harmonic_oscillator_wavefunction_section(N, A),

        "end-input\n"
    ])

    # print(input_string)
    with open(path, "w") as file:
        file.write(input_string)


def modify_hexahelicene(path):
    """ special case for hexahelicene """

    with open(path, "r") as file:
        data = file.read()

    # grab the parts we don't want to change
    before = data.split('PRIMITIVE-BASIS-SECTION')[0]
    after = data.split('end-primitive-basis-section')[1]
    middle = data.split('PRIMITIVE-BASIS-SECTION')[1].split('end-primitive-basis-section')[0]

    print(middle)
    import pdb; pdb.set_trace()


    hand_made_spf_section = f'SPF-BASIS-SECTION\nmulti-set{spfs}end-spf-basis-section\n'

    # glue the parts back together
    new_data = before + hand_made_spf_section + after

    with open(path, "w") as file:
        file.write(new_data)


def modify_input_file(path):
    """Replace the SPF-BASIS-SECTION with manually selected parameters
    """

    if project_name == "hexahelicene":
        modify_hexahelicene(path)
        return

    with open(path, "r") as file:
        data = file.read()

    # grab the parts we don't want to change
    before = data.split('SPF-BASIS-SECTION')[0]
    after = data.split('end-spf-basis-section')[1]

    if project_name == "h2o2":
        spfs = '''
      v01           =  6, 5, 6, 4, 6, 1
      v02           =  6, 3, 3, 3, 3, 1
      v03, v05      =  6, 4, 4, 4, 4, 1
      v04, v06      =  6, 3, 3, 3, 3, 1
'''
    elif project_name == "hcooh":
        spfs = '''
      v01, v08      =  5, 4, 5, 5, 5, 6, 1
      v02, v09      =  4, 3, 3, 3, 6, 3, 1
      v03, v04      =  6, 6, 5, 3, 3, 3, 1
      v05           =  6, 6, 6, 6, 6, 6, 1
      v06           =  6, 3, 5, 2, 3, 3, 1
      v07           =  4, 2, 3, 2, 3, 3, 1
'''
    elif project_name == "furan":
        spfs = '''
      v01, v11, v14, v15    =  2, 6, 2, 5, 6, 2, 6, 5, 1
      v02, v03, v16         =  2, 2, 2, 2, 2, 6, 2, 2, 1
      v04, v17, v18         =  2, 2, 4, 4, 2, 2, 2, 2, 1
      v05, v06, v19         =  2, 6, 5, 6, 6, 2, 2, 2, 1
      v07, v08, v20, v21    =  6, 5, 2, 7, 6, 2, 2, 2, 1
      v09, v12              =  2, 4, 2, 4, 2, 2, 2, 5, 1
      v10, v13              =  2, 2, 4, 3, 4, 2, 2, 5, 1
'''
    elif project_name == "formamide":
        spfs = '''
      v01, v10      =  6, 3, 4, 4, 4, 6, 1
      v02, v11      =  3, 4, 3, 3, 3, 3, 1
      v03           =  6, 3, 3, 3, 4, 4, 1
      v04           =  3, 5, 3, 3, 3, 3, 1
      v05, v06      =  5, 4, 6, 5, 6, 6, 1
      v07           =  2, 4, 3, 2, 3, 3, 1
      v08, v09      =  6, 4, 2, 2, 2, 3, 1
'''
    elif project_name == "vcm":
        spfs = '''
      v01            =  3, 3, 4, 3, 2, 3, 1
      v02, v10       =  3, 3, 6, 3, 2, 6, 1
      v04, v05       =  4, 5, 6, 3, 2, 3, 1
      v06            =  4, 6, 6, 3, 2, 3, 1
      v07, v08, v09  =  6, 3, 6, 2, 2, 2, 1
      v03, v11, v12  =  3, 6, 6, 3, 2, 6, 1
'''
    else:
        raise Exception(
            f"No hand built mctdh implemented for {project_name}. "
            "You need to modify `modify_input_file` in `populate_directories.py`."
        )

    hand_made_spf_section = f'SPF-BASIS-SECTION\nmulti-set{spfs}end-spf-basis-section\n'

    # glue the parts back together
    new_data = before + hand_made_spf_section + after

    with open(path, "w") as file:
        file.write(new_data)


def initalize_directories():
    '''Creates and populates all the directories where calculations will take place.'''

    original_model = vIO.read_raw_model_op_file(src_path_original_op_file)
    vIO.print_model(original_model, highest_order=1)

    for FC_or_not, coupling_order in it.product(["FC", "vibronic"], ["constant", "linear", "quadratic"]):
        if FC_or_not == "FC":
            continue

        if coupling_order == "constant":
            continue

        name = "_".join([project_name, FC_or_not, coupling_order])

        # modify op file
        model = vIO.create_deepcopy(original_model)

        vIO.remove_higher_order_terms(
            model, highest_order={"constant": 0, "linear": 1, "quadratic": 2}[coupling_order]
        )

        if FC_or_not == "FC":
            vIO.fill_offdiagonalsurfaces_of_model_with_zeros(model)

        # save file
        src_path_op_file = join(root_path_op_files, name + ".op")
        vIO.write_raw_model_op_file(src_path_op_file, model)

        for param_list in it.product(*expression_list):
            directory = join(work_root, dir_string.format(name, *param_list))
            os.makedirs(directory, exist_ok=True)

            input_file_path = join(directory, inp_file_name)
            print(input_file_path)
            create_input_file(input_file_path, name, *param_list)

            if '/tdh/' not in input_file_path and '/mctdh/' not in input_file_path:
                raise Exception('our check will fail')

            # we should only modify the spf definition if we are setting up mctdh calculations
            if '/mctdh/' in input_file_path:
                modify_input_file(input_file_path)

            dst_path_op = join(directory, name + ".op")
            dst_path_execution_script = join(directory, execution_script)
            copyfile(src_path_op_file, dst_path_op)
            copyfile(src_path_execution_script, dst_path_execution_script)


if __name__ == "__main__":
    initalize_directories()

# setup for furane Multi layer
# mlbasis-section

# 0> 9 9
#   # Electronic
#   1> [el]
#   # Vibrations
#   1> 5 5
# end-mlbasis-section


#       v01, v11, v14, v15    =  2, 6, 2, 5, 6, 2, 6, 5, 1
#       v02, v03, v16         =  2, 2, 2, 2, 2, 6, 2, 2, 1
#       v04, v17, v18         =  2, 2, 4, 4, 2, 2, 2, 2, 1
#       v05, v06, v19         =  2, 6, 5, 6, 6, 2, 2, 2, 1
#       v07, v08, v20, v21    =  6, 5, 2, 7, 6, 2, 2, 2, 1
#       v09, v12              =  2, 4, 2, 4, 2, 2, 2, 5, 1
#       v10, v13              =  2, 2, 4, 3, 4, 2, 2, 5, 1
