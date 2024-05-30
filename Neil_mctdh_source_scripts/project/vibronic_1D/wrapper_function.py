from .vibronic_model_keys import VibronicModelKeys as VMK
from . import vibronic_model_io as vIO
from . import model_op

def vibronic_input_reader(name, hamiltonian_truncation_order, path):
    """a wrapper function on top of Neil's I/O script to readin vibronic model parameter from standard op file"""
    raw_model = vIO.read_raw_model_op_file(
    f"{path}/{name}.op", get_transition_dipole_moment=True, surface_symmetrize=True,
            )
    # when we symmetrize the surfaces we divide the diagonal by 2 (that already happens and works correctly)
    A, N = vIO._extract_dimensions_from_dictionary(raw_model)
    if hamiltonian_truncation_order >=2:
        model_op.mode_symmetrize_from_upper_triangle_quadratic_terms(N, range(A), raw_model[VMK.G2])
            # now the raw_model has been symmetrized in modes as well

    if False:
        # divide all off-diagonal (vibrational) components by 2
        for i, j in it.product(range(N), range(N)):
            if i == j:
                continue
            raw_model[VMK.G2][i, j, :, :] /= 2

    if False:
        for k, v in raw_model.items():
            if isinstance(v, np.ndarray):
                if v.ndim > 2:
                    print(k, v)
                else:
                    print(k, v, sep='\n')
            else:
                print(k, v)

    # breakpoint()

    # at this point you have the model that can be put into VECC
    # and honestly just save it to json with
    #
    # vIO.save_model_to_JSON(f"{name}.json", raw_model)


    # model = vIO.load_model_from_JSON(f"{name}.json")

    # model = vIO.model_remove_ground_state(model)

    # swap electron / vibrational dimensions
    vIO.prepare_model_for_cc_integration(raw_model,hamiltonian_truncation_order)

    return raw_model
