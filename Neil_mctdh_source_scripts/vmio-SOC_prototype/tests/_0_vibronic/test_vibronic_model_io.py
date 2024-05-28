""" x """

# system imports
import random

# local imports
# from .. import context
from vmio.vibronic import vIO, VMK

# third party imports
import pytest
from pytest import raises
import numpy as np


@pytest.fixture()
def root(tmpdir_factory):
    return tmpdir_factory.mktemp("root")


@pytest.fixture()
def temp_file_path(root):
    p = root.join("tempfile.json")
    return str(p)


@pytest.fixture()
def A():
    return random.randint(2, 10)


@pytest.fixture()
def N():
    return random.randint(2, 10)


@pytest.fixture()
def random_model():
    return vIO.create_random_model()


@pytest.fixture()
def random_diagonal_model():
    return vIO.create_random_diagonal_model()


def test__array_is_symmetric_in_A():
    three_d = [[[1., 2.], [3., 4.]], [[5., 6.], [7., 8.]]]
    assert not vIO._array_is_symmetric_in_A(np.array(three_d))

    three_d = [[[1., 2.], [2., 4.]], [[5., 6.], [6., 8.]]]
    assert vIO._array_is_symmetric_in_A(np.array(three_d))


class Test_model_zeros_template_json_dict():

    def _confirm_model_zero_basic_parameters(self, model, A, N):

        # check for the right type, and that it isn't empty
        assert isinstance(model, dict)
        assert bool(model)

        # check that all elements have the correct value
        assert model[VMK.A] == A
        assert model[VMK.N] == N

    def _confirm_order_0_parameters(self, model, A, N):
        assert np.count_nonzero(model[VMK.E]) is 0
        assert np.count_nonzero(model[VMK.w]) is 0
        assert np.count_nonzero(model[VMK.etdm]) is 0
        assert np.count_nonzero(model[VMK.mtdm]) is 0

    def test_model_zeros_template_json_dict_default_args(self, A, N):
        model = vIO.model_zeros_template_json_dict(A, N)

        self._confirm_model_zero_basic_parameters(model, A, N)

        self._confirm_order_0_parameters(model, A, N)

        # default arguments return a 1st order model
        assert VMK.G1 in model and np.count_nonzero(model[VMK.G1]) is 0

        # other keys should not be present
        for key in VMK.coupling_list()[1:]:
            assert key not in model

    def test_model_zeros_template_json_dict_higher_order_args(self, A, N):
        for order in [2, 3, 4]:
            model = vIO.model_zeros_template_json_dict(A, N, highest_order=order)

            self._confirm_model_zero_basic_parameters(model, A, N)

            self._confirm_order_0_parameters(model, A, N)

            for key in VMK.coupling_list()[:order]:
                assert key in model and np.count_nonzero(model[key]) is 0

    # def test_model_zeros_template_json_dict_too_high_order(A, N):
    #     model = vIO.model_zeros_template_json_dict(A, N, highest_order=5)
    #     # add code to catch exception here


def test__generate_linear_terms(A, N):
    # confirm that the returned array is symmetric in surfaces
    # and that all values are in the correct range
    MAX = 1.0
    shape = vIO.model_shape_dict(A, N)
    displacement = np.random.uniform(0.1, MAX, size=N)
    linear_terms = np.zeros(shape[VMK.G1])

    vIO._generate_linear_terms(linear_terms, shape, displacement, range(N))

    assert np.all(-MAX <= linear_terms) and np.all(linear_terms <= MAX)
    assert vIO._array_is_symmetric_in_A(linear_terms)


def test__generate_quadratic_terms(A, N):
    # confirm that the returned array is symmetric in surfaces and modes
    # and that all values are in the correct range
    MAX = 1.0
    shape = vIO.model_shape_dict(A, N)
    displacement = np.random.uniform(0.1, MAX, size=(N, N))
    quadratic_terms = np.zeros(shape[VMK.G2])

    vIO._generate_quadratic_terms(quadratic_terms, shape, displacement, range(N))

    assert np.all(-MAX <= quadratic_terms) and np.all(quadratic_terms <= MAX)
    assert vIO._array_is_symmetric_in_A(quadratic_terms)


def test_generate_vibronic_model_data(A, N):

    max_E = random.uniform(1.0, 3.0)
    w = random.uniform(0.01, 0.03)
    scaling = random.uniform(0.02, 0.08)

    p = {'frequency_range': [w, 2*w],
         'energy_range': [0.0, max_E],
         'quadratic_scaling': 2*scaling,
         'linear_scaling': scaling,
         'diagonal': False,
         'numStates': A,
         'numModes': N,
         }

    model = vIO.generate_vibronic_model_data(p)
    assert bool(model)


class Test_extract_dimensions_of_model():

    def test_no_args(self):
        with raises(AssertionError) as e_info:
            vIO.extract_dimensions_of_model()
        assert str(e_info.value) == "no arguments provided"

    def test_path_only(self, temp_file_path, random_model):
        vIO.save_model_to_JSON(temp_file_path, random_model)
        A, N = vIO.extract_dimensions_of_model(path=temp_file_path)
        assert A == random_model[VMK.A]
        assert N == random_model[VMK.N]


class Test_extract_dimensions_of_diagonal_model():

    def test_no_args(self):
        with raises(AssertionError) as e_info:
            vIO.extract_dimensions_of_diagonal_model()
        assert str(e_info.value) == "no arguments provided"

    def test_path_only(self, temp_file_path, random_diagonal_model):
        vIO.save_diagonal_model_to_JSON(temp_file_path, random_diagonal_model)
        A, N = vIO.extract_dimensions_of_diagonal_model(path=temp_file_path)
        assert A == random_diagonal_model[VMK.A]
        assert N == random_diagonal_model[VMK.N]


def test__load_from_JSON_no_energy(temp_file_path, random_model):
    random_model[VMK.E].fill(0)
    vIO.save_model_to_JSON(temp_file_path, random_model)
    loaded_model = vIO._load_from_JSON(temp_file_path)
    assert vIO._same_model(random_model, loaded_model)


def test_load_sample_from_JSON_no_arrays_provided(temp_file_path, random_diagonal_model):
    vIO.save_diagonal_model_to_JSON(temp_file_path, random_diagonal_model)
    loaded_d_model = vIO.load_diagonal_model_from_JSON(temp_file_path, dictionary=None)
    assert vIO._same_diagonal_model(random_diagonal_model, loaded_d_model)
