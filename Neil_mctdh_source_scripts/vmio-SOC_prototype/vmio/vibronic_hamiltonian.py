""" Module Description

Some explanation / notes should go here

"""


# system imports
import os
from os.path import join
import itertools as it
import types
import sys

# third party imports
import scipy
# from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import numpy as np
import matplotlib as mpl; mpl.use('pdf')
import matplotlib.pyplot as plt
import parse  # used for loading data files
#
import opt_einsum as oe

# local imports
from project.log_conf import log
from project.two_mode_model import model_two_mode
from project.vibronic import vIO, VMK
from project.residual_equations import create_truncation_order_object
from project.residual_equations.symmetrize import symmetrize_tensor

# testing generated equations
from project.residual_equations import eT_zhz_eqs_H_2_P_4_T_1_exp_4_Z_1 as z_one_eqns
from project.residual_equations import eT_zhz_eqs_H_2_P_4_T_1_exp_4_Z_2 as z_two_eqns
from project.residual_equations import eT_zhz_eqs_H_2_P_4_T_1_exp_4_Z_3 as z_three_eqns

# tested modified integrator
from project.temporary_trial_solver import new_solve_ivp

# Force any printing of numpy arrays to be at a set precision
# to improve reading of output/debug print statements/logs
np.set_printoptions(precision=4, suppress=True, linewidth=200)

# factoring out the debug statements
print_force_constants = False
print_full_hamiltonian = False
print_U_matrix_and_percent_progress = True
check_hermicity_of_H_SOS = True
calculating_hot_bands = False

tab = " "*4

optimized = False


class vibronic_hamiltonian(object):
    """ vibronic_hamiltonian is python object """

    def __init__(
            self, model, model_name, build_H=False, HO_size=40, cc_truncation_order=1, hamiltonian_truncation_order=1,
            compare_with_mctdh=True, comparing_to_test_models=False,
            trans=False, theta=0, FC=False, W_truncation_order=0, solve_W=False, solve_Z=False,
            cal_Z0=False, cal_t0=True, solve_t=False, Z_truncation_order=2, T_truncation_order=1,
            selected_surface=[],
            force_use_optimize=False,
            calculate_population_flag=True,
    ):
        """
        trans: ?
        theta: only affects similarity_transformation (debugging?)
        """
        assert not comparing_to_test_models
        # model name
        self.model_name = model_name

        # number of potential energy surfaces
        self.A = model[VMK.A]

        # number of normal modes
        self.N = model[VMK.N]

        # using optimized paths only pays off for larger systems
        # or for higher order terms
        nof_dimensions = self.A + self.N

        self.op_einsum_flag = False

        if force_use_optimize or (nof_dimensions >= 15):
            self.op_einsum_flag = True
        else:
            self.op_einsum_flag = False

        self.calculate_population_flag = bool(calculate_population_flag)

        # vibronic model
        self.model = model

        # flag to determine whether it is FC model
        self.FC = FC

        # initialize truncation order of T, Z tensor
        self.Z_truncation_order = Z_truncation_order
        self.T_truncation_order = T_truncation_order

        if self.Z_truncation_order > 3:
            raise Exception('Only supports up to Z3 truncation, not Z4 or higher.')

        if self.T_truncation_order > 1:
            raise Exception('Only supports up to T1 truncation, not T2 or higher.')

        # ---------------------------------------------------------------------
        # make the truncation objects for the generated code
        # ---------------------------------------------------------------------
        # create the ground state/excited state ansatz
        self.ansatz = types.SimpleNamespace()
        self.ansatz.ground_state = True

        # create the truncation obj
        self.gen_trunc = types.SimpleNamespace()

        self.gen_trunc.h_at_least_linear = False
        self.gen_trunc.h_at_least_quadratic = False

        if hamiltonian_truncation_order >= 1:
            self.gen_trunc.h_at_least_linear = True
        if hamiltonian_truncation_order >= 2:
            self.gen_trunc.h_at_least_quadratic = True
        if hamiltonian_truncation_order >= 3:
            assert False, 'Not supported yet'

        self.gen_trunc.t_singles = False

        if T_truncation_order >= 1:
            self.gen_trunc.t_singles = True
        if T_truncation_order >= 2:
            assert False, 'Not supported yet'

        self.gen_trunc.z_at_least_linear = False
        self.gen_trunc.z_at_least_quadratic = False
        self.gen_trunc.z_at_least_cubic = False

        if Z_truncation_order >= 1:
            self.gen_trunc.z_at_least_linear = True
        if Z_truncation_order >= 2:
            self.gen_trunc.z_at_least_quadratic = True
        if Z_truncation_order >= 3:
            self.gen_trunc.z_at_least_cubic = True
        if Z_truncation_order >= 4:
            assert False, 'Not supported yet'

        # these functions probably need to be re-worked
        self.gen_trunc.confirm_at_least_singles = lambda: True
        self.gen_trunc.confirm_at_least_doubles = lambda: True
        self.gen_trunc.confirm_at_least_triples = lambda: True

        # ---------------------------------------------------------------------

        A, N = self.A, self.N

        # initialize derivative tensors
        self.dT = {
            0: np.zeros(A, dtype=complex),
        }

        if self.T_truncation_order >= 1:
            self.dT[1] = np.zeros((A, N), dtype=complex)

        if self.T_truncation_order >= 2:
            self.dT[2] = np.zeros((A, N, N), dtype=complex)

        self.dZ = {
            0: np.zeros((A, A), dtype=complex),
        }

        if self.Z_truncation_order >= 1:
            self.dZ[1] = np.zeros((A, A, N), dtype=complex)

        if self.Z_truncation_order >= 2:
            self.dZ[2] = np.zeros((A, A, N, N), dtype=complex)

        if self.Z_truncation_order >= 3:
            self.dZ[3] = np.zeros((A, A, N, N, N), dtype=complex)

        # ---------------------------------------------------------------------
        # initialize Lagrange multiplier
        self.la = np.zeros(self.A, dtype=complex)

        # initialize norm
        self.Norm = []

        # initialize state property
        if self.calculate_population_flag:
            self.State_pop_DB = []
            self.State_pop_AB = []

        # set up step size
        self.step_size = 5e-3

        # users defined selected surfaces to be decoupled
        self.selected_surface = selected_surface

        # vIO.temp_unswap_model_from_cc_integration(model, highest_order)
        # self.FC = vIO.model_is_FC(self.model)
        # vIO.prepare_model_for_cc_integration(model, highest_order)

        # create the `truncation_order_namedtuple` used in calculation residuals and W operators
        self.trunc = create_truncation_order_object(cc_truncation_order, hamiltonian_truncation_order)

        # print information related to truncation
        self._check_truncation_info()

        # electronic
        self.E_tdm = self.model[VMK.etdm].copy()
        # magnetic
        self.M_tdm = self.model[VMK.mtdm].copy()

        log.info("Electronic TDM:\n{:}".format(self.E_tdm.shape))
        log.info("Magnetic TDM:\n{:}".format(self.M_tdm.shape))

        self.HO_size, self.trans, self.theta = HO_size, trans, theta

        if print_force_constants:
            vIO.print_model(self.model)

        # save the flags for `initialize_hamiltonian`
        self.compare_with_mctdh = compare_with_mctdh
        self.comparing_to_test_models = comparing_to_test_models

        # the majority of the work
        self._initialize_hamiltonian()

        # precompute optimal einsum paths to speed up integration
        if self.op_einsum_flag:
            # self._initialize_einsum_paths()
            self._initialize_opt_einsum_paths()
            self._initialize_new_scheme_D_einsum_paths()
            self._initialize_cal_H_bar_tilde_einsum_paths()
            self._initialize_compute_C_matrix_einsum_paths()

        # these lines are only for testing against
        # Full configuration-interaction (FCI)
        if build_H:
            self._construct_debug_hamiltonian_for_FCI()

        # just double check that quadratic terms are symmetric
        for key in self.h.keys():
            rank = key[0] + key[1]
            if rank == 2:
                assert np.allclose(self.h[key], np.transpose(self.h[key], (0, 1, 3, 2))), (
                    f"{key = } is not symmetric"
                )

        return

    def _check_truncation_info(self):
        """ Some debugging and output just to make sure we didn't pick conflicting values """
        hamiltonian_order = vIO.extract_maximum_order_of_model(self.model)
        truncation_info = (
            '\n'
            f"The Hamiltonian provided contains at most {hamiltonian_order} order terms\n"
            # f"The calculation will considering up to    {self.trunc.hamiltonian_truncation_order} order terms\n"
            # f"Our Coupled Cluster calculation is: {self.trunc.cc_abbreviation()}\n")
        )

        truncation_info += f"Z truncation level is :{self.Z_truncation_order}\n"
        truncation_info += f"T truncation level is :{self.T_truncation_order}\n"

        log.info(truncation_info)

        # make sure we have sufficient CC truncation to calculate model to the order requested
        # assert self.trunc.cc_truncation_order >= self.trunc.hamiltonian_truncation_order, truncation_info
        # make sure that the Hamiltonian has sufficient terms to calculate `highest_order`
        assert hamiltonian_order >= self.trunc.hamiltonian_truncation_order, truncation_info

        # self.highest_order = highest_order
        # self.hamiltonian_order = hamiltonian_order
        return

    def _initialize_einsum_paths(self):
        """ Precompute einsum paths for the most computationally expensive parts of the integration.

        With optimized einsum paths runtime can be half as long.
        However since the paths are computed on every integration step, the longer the propagation time
        the more savings can be gained by precomputing static paths ahead of time.

        May also want to consult
        https://github.com/dgasmith/opt_einsum/blob/c826bb7df16f470a69f7bf90598fc27586209d11/docs/source/path_finding.rst
        """

        # ------------------------------------------------------------------------------------------------
        A, N = self.A, self.N

        # ostr = 'auto-hq'
        ostr = 'optimal'

        # ------------------------------------------------------------------------------------------------
        # constant z
        # ------------------------------------------------------------------------------------------------
        self.f_z_0_op_list = [
            oe.contract_expression('k,yx,xk->y', (N, ), (A, A), (A, N), optimize=ostr),
            oe.contract_expression('k,yxk,x->y', (N, ), (A, A, N), (A, ), optimize=ostr),
            oe.contract_expression('k,yxlk,xl->y', (N, ), (A, A, N, N), (A, N), optimize=ostr),
            oe.contract_expression('k,l,yxl,xk->y', (N, ), (N, ), (A, A, N), (A, N), optimize=ostr),
            oe.contract_expression('k,l,yxkl,x->y', (N, ), (N, ), (A, A, N, N), (A, ), optimize=ostr),
            oe.contract_expression('k,l,m,yxlm,xk->y', (N, ), (N, ), (N, ), (A, A, N, N), (A, N), optimize=ostr),
        ]

        if self.Z_truncation_order >= 2:
            self.f_z_0_op_list.extend([
                oe.contract_expression('k,yxl,xkl->y', (N, ), (A, A, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,yx,xkl->y', (N, ), (N, ), (A, A), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,yxml,xkm->y', (N, ), (N, ), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,m,yxm,xkl->y', (N, ), (N, ), (N, ), (A, A, N), (A, N, N), optimize=ostr),
                oe.contract_expression('i,j,k,l,yxij,xkl->y', (N, ), (N, ), (N, ), (N, ), (A, A, N, N), (A, N, N), optimize=ostr),

            ])

        # ------------------------------------------------------------------------------------------------
        # linear z
        # ------------------------------------------------------------------------------------------------
        self.f_z_i_op_list = [
            oe.contract_expression('yxi,x->yi', (A, A, N), (A,), optimize=ostr),
            oe.contract_expression('yxki,xk->yi', (A, A, N, N), (A, N), optimize=ostr),
            oe.contract_expression('yx,xi->yi', (A, A), (A, N), optimize=ostr),
            # adding terms associated e^T_dagger
            oe.contract_expression('k,yxi,xk->yi', (N,), (A, A, N), (A, N), optimize=ostr),
            oe.contract_expression('k,yxk,xi->yi', (N,), (A, A, N), (A, N), optimize=ostr),
            oe.contract_expression('k,yxki,x->yi', (N,), (A, A, N, N), (A,), optimize=ostr),
            oe.contract_expression('k,l,yxli,xk->yi', (N,), (N,), (A, A, N, N), (A, N), optimize=ostr),
            oe.contract_expression('k,l,yxkl,xi->yi', (N,), (N,), (A, A, N, N), (A, N), optimize=ostr)
        ]

        if self.Z_truncation_order >= 2:
            self.f_z_i_op_list.extend([
                oe.contract_expression('yxk,xki->yi', (A, A, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,yx,xki->yi', (N,), (A, A), (A, N, N), optimize=ostr),
                oe.contract_expression('k,yxli,xkl->yi', (N,), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,yxlk,xli->yi', (N,), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,yxl,xki->yi', (N,), (N,), (A, A, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,yxi,xkl->yi', (N,), (N,), (A, A, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,m,yxmi,xkl->yi', (N,), (N,), (N,), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,m,yxlm,xki->yi', (N,), (N,), (N,), (A, A, N, N), (A, N, N), optimize=ostr)
            ])

        # still testing different versions
        if False and self.T_truncation_order >= 2:
            self.f_z_i_op_list.extend([
                oe.contract_expression('kl,yxki,xl->yi', (N, N), (A, A, N, N), (A, N), optimize=ostr),
                oe.contract_expression('kl,yxkl,xi->yi', (N, N), (A, A, N, N), (A, N), optimize=ostr),
            ])

        # still testing different versions
        if False and self.Z_truncation_order >= 2 and self.T_truncation_order >= 2:
            self.f_z_i_op_list.extend([
                oe.contract_expression('kl,yxk,xli->yi', (N, N), (A, A, N), (A, N, N), optimize=ostr),
                oe.contract_expression('kl,yxi,xkl->yi', (N, N), (A, A, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,lm,yxli,xkm->yi', (N, ), (N, N), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,lm,yxkl,xmi->yi', (N, ), (N, N), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,lm,yxki,xlm->yi', (N, ), (N, N), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,lm,yxlm,xki->yi', (N, ), (N, N), (A, A, N, N), (A, N, N), optimize=ostr),
            ])

        # ------------------------------------------------------------------------------------------------
        # quadratic z
        # ------------------------------------------------------------------------------------------------
        self.f_z_ij_op_list = [
            oe.contract_expression('yx,xij->yij', (A, A), (A, N, N), optimize=ostr),
            oe.contract_expression('yxi,xj->yij', (A, A, N), (A, N), optimize=ostr),
            oe.contract_expression('yxij,x->yij', (A, A, N, N), (A, ), optimize=ostr),
            oe.contract_expression('yxki,xkj->yij', (A, A, N, N), (A, N, N), optimize=ostr),
            # adding terms associated e^T_dagger
            oe.contract_expression('k,yxi,xkj->yij', (N, ), (A, A, N), (A, N, N), optimize=ostr),
            oe.contract_expression('k,yxk,xij->yij', (N, ), (A, A, N), (A, N, N), optimize=ostr),
            oe.contract_expression('k,yxij,xk->yij', (N, ), (A, A, N, N), (A, N), optimize=ostr),
            oe.contract_expression('k,yxki,xj->yij', (N, ), (A, A, N, N), (A, N), optimize=ostr),
            oe.contract_expression('k,l,yxli,xkj->yij', (N, ), (N, ), (A, A, N, N), (A, N, N), optimize=ostr),
            oe.contract_expression('k,l,yxij,xkl->yij', (N, ), (N, ), (A, A, N, N), (A, N, N), optimize=ostr),
            oe.contract_expression('k,l,yxkl,xij->yij', (N, ), (N, ), (A, A, N, N), (A, N, N), optimize=ostr),
        ]

        # still testing different versions
        if False and self.T_truncation_order >= 2:
            self.f_z_ij_op_list.extend([
                oe.contract_expression('kl,yxki,xlj->yij', (N, N), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('kl,yxkj,xli->yij', (N, N), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('kl,yxij,xkl->yij', (N, N), (A, A, N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('kl,yxkl,xij->yij', (N, N), (A, A, N, N), (A, N, N), optimize=ostr),
            ])

        # linear t
        # self.f_t[1]_op_list = [
        #     oe.contract_expression('abk,bki->abi', (A, A, N), (A, N, N)),
        #     oe.contract_expression('abkl,bk,bli->abi', (A, A, N, N), (A, N), (A, N, N)),
        # ]

        # # quadratic t
        # self.f_t[2]_op_list = [
        #     oe.contract_expression('abkj,bki->abij', (A, A, N, N), (A, N, N)),
        #     oe.contract_expression('abki,bkj->abij', (A, A, N, N), (A, N, N)),
        #     oe.contract_expression('abkl,bli,bkj->abij', (A, A, N, N), (A, N, N), (A, N, N)),
        # ]

        # ------------------------------------------------------------------------------------------------
        # constant z residual
        # ------------------------------------------------------------------------------------------------
        self.f_z_0_residual_list = []

        if self.Z_truncation_order >= 2:
            self.f_z_0_residual_list.extend([
                oe.contract_expression('l,m,ylm->y', (N, ), (N, ), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,ykl->y', (N, ), (N, ), (A, N, N), optimize=ostr),
            ])

        if self.T_truncation_order >= 2:
            self.f_z_0_residual_list.extend([
                oe.contract_expression('k,l,kl,y->y', (N, ), (N, ), (N, N), (A, ), optimize=ostr),
                oe.contract_expression('k,l,m,lm,yk->y', (N, ), (N, ), (N, ), (N, N), (A, N), optimize=ostr),
            ])

        if self.T_truncation_order >= 2 and self.Z_truncation_order >= 2:
            self.f_z_0_residual_list.extend([
                oe.contract_expression('k,l,m,n,kl,ymn->y', (N, ), (N, ), (N, ), (N, ), (N, N), (A, N, N), optimize=ostr),
            ])

        # ------------------------------------------------------------------------------------------------
        # linear z residual
        # ------------------------------------------------------------------------------------------------
        self.f_z_1_residual_list = []

        if self.Z_truncation_order >= 2:
            self.f_z_1_residual_list.extend([
                oe.contract_expression('l,yli->yi', (N, ), (A, N, N), optimize=ostr),
                oe.contract_expression('k,yki->yi', (N, ), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,i,ykl->yi', (N, ), (N, ), (N, ), (A, N, N), optimize=ostr),
            ])

        if self.T_truncation_order >= 2:
            self.f_z_1_residual_list.extend([
                oe.contract_expression('k,ki,y->yi', (N, ), (N, N), (A, ), optimize=ostr),
                oe.contract_expression('k,l,ki,yl->yi', (N, ), (N, ), (N, N), (A, N), optimize=ostr),
                oe.contract_expression('k,l,kl,yi->yi', (N, ), (N, ), (N, N), (A, N), optimize=ostr),
            ])

        if self.T_truncation_order >= 2 and self.Z_truncation_order >= 2:
            self.f_z_1_residual_list.extend([
                oe.contract_expression('k,l,m,mi,ykl->yi', (N, ), (N, ), (N, ), (N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,m,lm,yki->yi', (N, ), (N, ), (N, ), (N, N), (A, N, N), optimize=ostr),
            ])

        # ------------------------------------------------------------------------------------------------
        # quadratic z residual
        # ------------------------------------------------------------------------------------------------
        self.f_z_2_residual_list = [
            oe.contract_expression('i,yj->yij', (N, ), (A, N), optimize=ostr),
        ]

        if self.Z_truncation_order >= 2:
            self.f_z_2_residual_list.extend([
                oe.contract_expression('k,i,ykj->yij', (N, ), (N, ), (A, N, N), optimize=ostr),
            ])

        if self.T_truncation_order >= 2:
            self.f_z_2_residual_list.extend([
                oe.contract_expression('k,ki,yj->yij', (N, ), (N, N), (A, N), optimize=ostr),
                oe.contract_expression('ij,y->yij', (N, N), (A, ), optimize=ostr),
                oe.contract_expression('k,ij,yk->yij', (N, ), (N, N), (A, N), optimize=ostr),
            ])

        if self.T_truncation_order >= 2 and self.Z_truncation_order >= 2:
            self.f_z_2_residual_list.extend([
                oe.contract_expression('k,l,ki,ylj->yij', (N, ), (N, ), (N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,kl,yij->yij', (N, ), (N, ), (N, N), (A, N, N), optimize=ostr),
                oe.contract_expression('k,l,ij,ykl->yij', (N, ), (N, ), (N, N), (A, N, N), optimize=ostr),
            ])

        return

    def _initialize_opt_einsum_paths(self):
        """Store optimized paths as member of `self` for use in the `eT_zhz_eqs_***` module during integration."""

        if self.Z_truncation_order == 1:
            self.all_opt_paths = z_one_eqns.compute_all_optimized_paths(self.A, self.N, self.ansatz, self.gen_trunc)
        if self.Z_truncation_order == 2:
            self.all_opt_paths = z_two_eqns.compute_all_optimized_paths(self.A, self.N, self.ansatz, self.gen_trunc)
        if self.Z_truncation_order == 3:
            self.all_opt_paths = z_three_eqns.compute_all_optimized_paths(self.A, self.N, self.ansatz, self.gen_trunc)

        return

    def _initialize_new_scheme_D_einsum_paths(self):
        """Store optimized paths as member of `self` for use in the scheme
        This function should eventually  be replaced with someone more permanent.
        But for testing purposes for now its sufficient
        """
        A, N = self.A, self.N

        self.d0_opt_paths = [
            oe.contract_expression('k,yk->y', (N,), (A, N), optimize='auto-hq'),
            oe.contract_expression('k,l,ykl->y', (N,), (N,), (A, N, N), optimize='auto-hq'),
            oe.contract_expression('k,l,m,yklm->y', (N,), (N,), (N,), (A, N, N, N), optimize='auto-hq'),
        ]

        if self.Z_truncation_order >= 1:
            self.d1_opt_paths = [
                oe.contract_expression('k,yki->yi', (N,), (A, N, N), optimize='auto-hq'),
                oe.contract_expression('k,l,ykli->yi', (N,), (N,), (A, N, N, N), optimize='auto-hq'),
            ]

        if self.Z_truncation_order >= 2:
            self.d2_opt_paths = [
                oe.contract_expression('l,ylij->yij', (N,), (A, N, N, N), optimize='auto-hq'),
            ]

        return

    def _initialize_cal_H_bar_tilde_einsum_paths(self):
        """Store optimized paths as member of `self` for use in `_cal_H_bar_tilde`
        This function should eventually be replaced with someone more permanent.
        But for optimization for calculating Fe(CO)5 it is fine.
        """
        A, N = self.A, self.N

        self.H_bar_tilde_paths = [
            oe.contract_expression('k,abk->ab', (N,), (A, A, N), optimize='auto-hq'),
            oe.contract_expression('k,l,abkl->ab', (N,), (N,), (A, A, N, N), optimize='auto-hq'),
            oe.contract_expression('k,abik->abi', (N,), (A, A, N, N), optimize='auto-hq'),
            oe.contract_expression('k,abki->abi', (N,), (A, A, N, N), optimize='auto-hq'),
        ]

        return

    def _initialize_compute_C_matrix_einsum_paths(self):
        """Store optimized paths as member of `self` for use in `_compute_C_matrix`
        This function should eventually be replaced with someone more permanent.
        But for optimization for calculating Fe(CO)5 it is fine.
        """
        A, N = self.A, self.N

        self.Cmat_Z0_opt_paths = [
            oe.contract_expression('k,xk->x', (N,), (A, N), optimize='auto-hq'),
        ]

        if self.Z_truncation_order >= 2:

            self.Cmat_Z2_opt_paths = [
                # s0
                oe.contract_expression('k,l,xkl->x', (N,), (N,), (A, N, N), optimize='auto-hq'),
                # s1
                oe.contract_expression('k,xik->xi', (N,), (A, N, N), optimize='auto-hq'),
                # no s2
            ]

        if self.Z_truncation_order >= 3:

            self.Cmat_Z3_opt_paths = [
                # s0
                oe.contract_expression('k,l,m,xklm->x', (N,), (N,), (N,), (A, N, N, N), optimize='auto-hq'),
                # s1
                oe.contract_expression('k,l,xikl->xi', (N,), (N,), (A, N, N, N), optimize='auto-hq'),
                # s2
                oe.contract_expression('k, xkij->xij', (N,), (A, N, N, N), optimize='auto-hq'),
            ]

        return

    def _initialize_hamiltonian(self, build_H=False):
        """ Express the Hamiltonian in a finite H.O. basis for comparing to FCI
        The basis size is defined by self.HO_size

        Explanation of the mapping between the code and the theory:
        H with a single annihilation operator is represented by
            - self.h[(1, 0)]
        H with a single creation operator is represented by
            - self.h[(0, 1)]
        """

        N, A = self.N, self.A  # for brevity

        # define coefficient tensors
        self.h = {
            (0, 0): np.zeros((A, A), dtype=complex),
            (0, 1): np.zeros((A, A, N), dtype=complex),
            (1, 0): np.zeros((A, A, N), dtype=complex),
            (1, 1): np.zeros((A, A, N, N), dtype=complex),
            (0, 2): np.zeros((A, A, N, N), dtype=complex),
            (2, 0): np.zeros((A, A, N, N), dtype=complex)
        }

        # read in input parameters
        # --------------------------------------------------------------------------
        # specific prefactors for comparing with MCTDH
        # this is the default mode!
        if self.compare_with_mctdh:
            log.info("parameter of the Hamiltonian are entered properly!")

            # energy
            self.h[(0, 0)] += self.model[VMK.E].copy()

            # H.O. ground state energy
            for a in range(A):
                self.h[(0, 0)][a, a] += 0.5 * np.sum(self.model[VMK.w])

            log.info("zero point energy: {:.8f} ev".format(0.5 * np.sum(self.model[VMK.w])))

            # frequencies
            for a, j in it.product(range(A), range(N)):
                self.h[(1, 1)][a, a, j, j] += self.model[VMK.w][j].copy()

            # linear terms
            if self.trunc.at_least_linear:
                self.h[(0, 1)] += self.model[VMK.G1].copy() / np.sqrt(2)
                self.h[(1, 0)] += self.model[VMK.G1].copy() / np.sqrt(2)

            # quadratic terms
            if self.trunc.at_least_quadratic:
                # quadratic correction to H.O. ground state energy
                self.h[(0, 0)] += 0.5 * np.trace(self.model[VMK.G2], axis1=2, axis2=3)

                # quadratic terms
                self.h[(1, 1)] += self.model[VMK.G2].copy()
                self.h[(2, 0)] += self.model[VMK.G2].copy()
                self.h[(0, 2)] += self.model[VMK.G2].copy()

            # if computing Frank-Condon(FC) model then
            # zero out all electronically-diagonal terms
            if self.FC:
                for a, b in it.product(range(A), repeat=2):
                    if a != b:
                        self.h[(1, 1)][a, b, :] = np.zeros([N, N])
                        self.h[(2, 0)][a, b, :] = np.zeros([N, N])
                        self.h[(0, 2)][a, b, :] = np.zeros([N, N])
                        self.h[(1, 0)][a, b, :] = np.zeros([N, ])
                        self.h[(0, 1)][a, b, :] = np.zeros([N, ])

        # specific prefactors for comparing comparing to fixed results for pytest
        elif self.comparing_to_test_models:
            # energy
            self.h[(0, 0)] += self.model[VMK.E].copy()

            # frequencies
            for a, j in it.product(range(A), range(N)):
                self.h[(1, 1)][a, a, j, j] += 4 * self.model[VMK.w][j].copy()

            # linear terms
            if self.trunc.at_least_linear:
                self.h[(0, 1)] += self.model[VMK.G1].copy()
                self.h[(1, 0)] += self.model[VMK.G1].copy()

            # quadratic terms
            if self.trunc.at_least_quadratic:
                self.h[(1, 1)] += 4 * self.model[VMK.G2].copy()
                self.h[(2, 0)] += 2 * self.model[VMK.G2].copy()
                self.h[(0, 2)] += 2 * self.model[VMK.G2].copy()

        else:
            Exception("You shouldn't reach this point")

        return

    def _construct_debug_hamiltonian_for_FCI(self):
        """ Should check if this function needs to be removed
        or reworked in the future?
        """

        shape = [self.A, self.HO_size**2, self.A, self.HO_size**2]
        H = np.zeros(shape, dtype=complex)

        for a, b in it.product(range(self.A), repeat=2):
            log.info(f"Doing SOS with ({a}, {b})")

            # setup the Hamiltonian
            model = model_two_mode(
                self.h[(0, 0)][a, b],
                self.h[(0, 1)][a, b, :],
                self.h[(1, 1)][a, b, :],
                self.h[(2, 0)][a, b, :]
            )

            # compute the SOS solution
            H[a, :, b][:] = model.sos_solution(HO_size=self.HO_size)

        # reshape into a matrix
        self.H = H.reshape(self.A * self.HO_size ** 2, self.A * self.HO_size ** 2)
        return

    def _similarity_trans(self, H_args, t_args):
        """ Apply a similarity transformation to the Hamiltonian.
        This is applied for each electronic surface.
        """

        A, N = self.A, self.N   # for readability

        """ TO DO: """
        # rewrite residue for T to avoid permutations of electronic labels
        def f_t_0(H, T):
            """return residue R_0"""

            # initialize as zero
            R = np.zeros([A, A], dtype=complex)

            # constant
            R += H[(0, 0)]

            # linear
            if self.T_truncation_order >= 1:
                R += np.einsum('abk,k->ab', H[(0, 1)], T[1])
                if self.trunc.at_least_quadratic:
                    R += 0.5 * np.einsum('abkl,k,l->ab', H[(0, 2)], T[1], T[1])

            # quadratic
            if self.T_truncation_order >= 2:
                R += 0.5 * np.einsum('abkl,kl->ab', H[(0, 2)], T[2])

            return R

        def f_t_I(H, T):
            """return residue R_I"""

            # initialize as zero
            R = np.zeros([A, A, N], dtype=complex)

            # linear
            R += H[(0, 1)]

            # quadratic
            if self.T_truncation_order >= 1:
                R += np.einsum('abik,k->abi', H[(0, 2)], T[1])

            return R

        def f_t_i(H, T):
            """return residue R_i"""

            # initialize
            R = np.zeros([A, A, N], dtype=complex)

            # non zero initial value of R
            R += H[(1, 0)]

            # linear
            if self.T_truncation_order >= 1:
                R += np.einsum('abki,k->abi', H[(1, 1)], T[1])

                # quadratic
                if self.T_truncation_order >= 2:
                    R += np.einsum('abk,ki->abi', H[(0, 1)], T[2])
                    R += np.einsum('abkl,k,li->abi', H[(0, 2)], T[1], T[2])

                    # cubic
                    if self.T_truncation_order >= 3:
                        R += np.einsum('abkl,kli->abi', H[(0, 2)], T[3])

            return R

        def f_t_Ij(H, T):
            """return residue R_Ij"""

            # initialize
            R = np.zeros([A, A, N, N], dtype=complex)

            # first term
            R += H[(1, 1)]

            # quadratic
            if self.T_truncation_order >= 2:
                R += np.einsum('abik,kj->abij', H[(0, 2)], T[2])

            return R

        def f_t_IJ(H, T):
            """return residue R_IJ"""

            # initialize as zero
            R = np.zeros([A, A, N, N], dtype=complex)

            # if self.hamiltonian_truncation_order >= 2:
            # quadratic
            if self.T_truncation_order >= 2:
                R += H[(0, 2)]
                pass

            return R

        def f_t_ij(H, T):
            """return residue R_ij"""

            # # initialize as zero
            R = np.zeros([A, A, N, N], dtype=complex)

            # if self.hamiltonian_truncation_order >= 2:

            # quadratic
            if self.T_truncation_order >= 2:
                R += H[(2, 0)]  # h term
                R += np.einsum('abkj,ki->abij', H[(1, 1)], T[2])
                R += np.einsum('abki,kj->abij', H[(1, 1)], T[2])
                R += 0.5 * np.einsum('abkl,ki,lj->abij', H[(0, 2)], T[2], T[2])
                R += 0.5 * np.einsum('abkl,kj,li->abij', H[(0, 2)], T[2], T[2])
            return R

        # compute similarity transformed Hamiltonian over e^T
        sim_h = {}
        sim_h[(0, 0)] = f_t_0(H_args, t_args)
        sim_h[(0, 1)] = f_t_I(H_args, t_args)
        sim_h[(1, 0)] = f_t_i(H_args, t_args)
        sim_h[(1, 1)] = f_t_Ij(H_args, t_args)
        sim_h[(0, 2)] = f_t_IJ(H_args, t_args)
        sim_h[(2, 0)] = f_t_ij(H_args, t_args)

        return sim_h

    def apply_time_step_conversion(self, step_size):
        """Change the units of each step to be in eV/fs"""

        # use consistent units of energy in eV, time in fs
        hbar_per_eV = 6.582119569e-16
        s_per_fs = 1e-15

        # apply the conversion
        new_step_size = step_size / hbar_per_eV * s_per_fs

        return new_step_size

    def remove_time_step_conversion(self, time_step_array):
        """Remove the units for each time step so that it can be compared to SOS"""

        # use consistent units of energy in eV, time in fs
        hbar_per_eV = 6.582119569e-16
        s_per_fs = 1e-15

        # remove the conversion
        time_step_array *= hbar_per_eV / s_per_fs

        if type(time_step_array) is not np.ndarray:
            # if we were in fact given a number then return it
            return time_step_array

        # otherwise we modified the array by-reference
        return

    def _postprocess_rk45_integration_results(self, sol, debug=False):
        """ extract the relevant information from the integrator object `sol` """
        # number of integration steps accepted by the integrator

        log.info(f"RK45 preformed {len(self.C_tau_ECD)} integration calculations.")
        log.info(f"RK45 accepted  {len(sol.t)} of those as solutions to the ode's.")
        if debug:
            log.debug(f"Distance we reached when we stopped: {sol.t_events[0]}")

        # Copy the time value arrays
        self.t_cc = sol.t.copy()
        # log.info(len(self.t_cc))
        # log.info(len(self.t_cc_store))

        # shift the time back to units of fs
        self.remove_time_step_conversion(self.t_cc)

        # initialize the arrays to store the autocorrelation function
        self.C_tau_cc_ABS = np.zeros_like(self.t_cc, dtype=complex)
        self.C_tau_cc_ECD = np.zeros_like(self.t_cc, dtype=complex)

        # log.info(len(self.C_tau_cc))
        # log.info(len(self.C_tau_cc_store))

        # only extract the values which correspond to time steps in the solution
        # since we save C(t) for all integration steps, but only some are accepted
        C_dic_ABS = {c[0]: c[1] for c in self.C_tau_ABS}
        for idx, t in enumerate(sol.t):
            self.C_tau_cc_ABS[idx] = C_dic_ABS[t]

        C_dic_ECD = {c[0]: c[1] for c in self.C_tau_ECD}
        for idx, t in enumerate(sol.t):
            self.C_tau_cc_ECD[idx] = C_dic_ECD[t]

        # store norm of the wave function
        self.Norm_cc = np.zeros([len(self.t_cc), self.A], dtype=complex)

        N_dic = {N[0]: N[1] for N in self.Norm}
        for idx, t in enumerate(sol.t):
            self.Norm_cc[idx, :] = N_dic[t]

        if self.calculate_population_flag:
            # store the diabatic populations --------------------
            self.state_pop_cc_DB = np.zeros([len(self.t_cc), self.A, self.A], dtype=complex)

            P_dic = {P[0]: P[1] for P in self.State_pop_DB}

            for idx, t in enumerate(sol.t):
                self.state_pop_cc_DB[idx, :] += P_dic[t]

            # store diabatic state population
            dic_diab_pop = {}
            dic_diab_pop["time(fs)"] = self.t_cc

            for x, y in it.product(range(self.A), repeat=2):
                name = f"diab state pop {x} init state {y}"
                dic_diab_pop[name] = self.state_pop_cc_DB[:, x, y].real

            import pandas as pd
            df = pd.DataFrame(dic_diab_pop)
            df.to_csv(f"Diabatic_state_pop_{self.model_name}.csv", index=False)

            # store the adiabatic populations --------------------
            self.state_pop_cc_AB = np.zeros([len(self.t_cc), self.A, self.A], dtype=complex)

            P_dic = {P[0]: P[1] for P in self.State_pop_AB}
            for idx, t in enumerate(sol.t):
                self.state_pop_cc_AB[idx, :] += P_dic[t]
            # ---------------------------------------------------------------------------------

        if debug:
            # only extract the values which correspond to time steps in the solution
            # since we save them for all integration steps, but only some are accepted
            if self.highest_order >= 0:
                self.dS_debug = np.zeros((len(sol.t), self.A, self.A), dtype=complex)
                dS_dic = {x[0]: x[1] for x in self.dS_store}
                for idx, t in enumerate(sol.t):
                    self.dS_debug[idx] = dS_dic[t]
            # ---------------------------------------------------------------------------------
            if self.highest_order >= 1:
                self.dt_i_debug = np.zeros((len(sol.t), self.A, self.A, self.N), dtype=complex)
                dt_i_dic = {x[0]: x[1] for x in self.dt_i_store}
                for idx, t in enumerate(sol.t):
                    self.dt_i_debug[idx] = dt_i_dic[t]
            # ---------------------------------------------------------------------------------
            if self.highest_order >= 2:
                self.dt_ij_debug = np.zeros((len(sol.t), self.A, self.A, self.N, self.N), dtype=complex)
                dt_ij_dic = {x[0]: x[1] for x in self.dt_ij_store}
                for idx, t in enumerate(sol.t):
                    self.dt_ij_debug[idx] = dt_ij_dic[t]
            # ---------------------------------------------------------------------------------

        log.debug(f"Status message: {sol.message}")  # description of termination reason
        log.debug(f"status: {sol.status}")  # -1: step failed, 0: reached end of tspan, 1: termination event occurred
        log.debug(f"Succeeded?: {sol.success}")  # bool if reached end of interval or termination event occurred

        return

    def _unravel_y_tensor(self, y_tensor):
        """ Restore the original shape of the flattened y tensor """

        A, N = self.A, self.N  # for brevity

        # all return tensors start as None
        Z = {0: None, 1: None, 2: None, 3: None}
        T = {0: None, 1: None, 2: None}

        # ------------------------------ restore z tensor ----------------------------

        # constant terms
        start_constant_slice_index = 0
        end_constant_slice_index = start_constant_slice_index + A * A
        Z[0] = np.reshape(
            y_tensor[start_constant_slice_index:end_constant_slice_index],
            newshape=(A, A)
        )

        if self.Z_truncation_order >= 1:
            # linear terms
            start_linear_slice_index = end_constant_slice_index
            end_linear_slice_index = start_linear_slice_index + A * A * N
            Z[1] = np.reshape(
                y_tensor[start_linear_slice_index: end_linear_slice_index],
                newshape=(A, A, N)
            )

        if self.Z_truncation_order >= 2:
            # quadratic terms
            start_quadratic_slice_index = end_linear_slice_index
            end_quadratic_slice_index = start_quadratic_slice_index + A * A * N * N
            Z[2] = np.reshape(
                y_tensor[start_quadratic_slice_index: end_quadratic_slice_index],
                newshape=(A, A, N, N)
            )

        if self.Z_truncation_order >= 3:
            # cubic terms
            start_cubic_slice_index = end_quadratic_slice_index
            end_cubic_slice_index = start_cubic_slice_index + A * A * N * N * N
            Z[3] = np.reshape(
                y_tensor[start_cubic_slice_index: end_cubic_slice_index],
                newshape=(A, A, N, N, N)
            )

        # ------------------------------ restore t tensor ----------------------------

        # need a bit of logic to figure out the final Z slice index

        start_constant_slice_index = end_constant_slice_index

        if self.Z_truncation_order >= 1:
            start_constant_slice_index = end_linear_slice_index

        if self.Z_truncation_order >= 2:
            start_constant_slice_index = end_quadratic_slice_index

        if self.Z_truncation_order >= 3:
            start_constant_slice_index = end_cubic_slice_index

        # constant terms
        end_constant_slice_index = start_constant_slice_index + A
        T[0] = np.reshape(
            y_tensor[start_constant_slice_index: end_constant_slice_index],
            newshape=(A, )
        )

        if self.T_truncation_order >= 1:
            # linear terms
            start_linear_slice_index = end_constant_slice_index
            end_linear_slice_index = start_linear_slice_index + A * N
            T[1] = np.reshape(
                y_tensor[start_linear_slice_index: end_linear_slice_index],
                newshape=(A, N)
            )

        if self.T_truncation_order >= 2:
            # quadratic terms
            start_quadratic_slice_index = end_linear_slice_index
            end_quadratic_slice_index = start_quadratic_slice_index + A * N * N
            T[2] = np.reshape(
                y_tensor[start_quadratic_slice_index: end_quadratic_slice_index],
                newshape=(A, N, N)
            )

        return Z, T

    def _ravel_y_tensor(self, Z, T):
        """ Flatten the `t` and `z` tensors into a 1D array """

        z_tensor_list = [Z[0].ravel(), ]
        if self.Z_truncation_order >= 1:
            z_tensor_list.append(Z[1].ravel())
        if self.Z_truncation_order >= 2:
            z_tensor_list.append(Z[2].ravel())
        if self.Z_truncation_order >= 3:
            z_tensor_list.append(Z[3].ravel())

        t_tensor_list = [T[0].ravel(), ]
        if self.T_truncation_order >= 1:
            t_tensor_list.append(T[1].ravel())
        if self.T_truncation_order >= 2:
            t_tensor_list.append(T[2].ravel())

        # the t tensor should come before the z tensor
        y_tensor = np.concatenate((*z_tensor_list, *t_tensor_list))

        return y_tensor

    def _print_integration_progress(self, time, t_final, *args):
        """ Prints to stdout every 1e4 steps or if current fs value is a multiple of (0.1 * `t_final`). """

        # unpack any args we wish to print
        Z, T = args

        self.counter += 1
        self.last_counter += 1

        fs = self.remove_time_step_conversion(time)

        at_least_1_femtosecond_has_passed_since_last_print = bool((fs-self.last_print) > 1.0)
        time_is_a_multiple_of_ten_percent_of_t_final = np.isclose(round(fs, 1) % (t_final / 10), 0.1)

        print_flag = bool(
            self.last_counter >= int(1e4)
            or
            (
                at_least_1_femtosecond_has_passed_since_last_print
                and time_is_a_multiple_of_ten_percent_of_t_final
            )
        )

        if print_flag:
            log.info(
                f"On integration step {self.counter:<8d} at {fs:>9.4f}fs\n"
                f"C(t)(ABS) = {self.C_tau_ABS[-1][1]:>9.4f}\n"
                f"C(t)(ECD) = {self.C_tau_ECD[-1][1]:>9.4f}\n"
            )

            # -------------------------------------------------------------------------------------
            t_amplitude_values_string = f'max constant T[0]  amplitude: {abs(T[0]).max()}\n'

            if self.T_truncation_order >= 1:
                t_amplitude_values_string += f'max single   T[1]  amplitude: {abs(T[1]).max()}\n'
            if self.T_truncation_order >= 2:
                t_amplitude_values_string += f'max double   T[2]  amplitude: {abs(T[2]).max()}\n'
            if self.T_truncation_order >= 3:
                t_amplitude_values_string += f'max triple   T[3]  amplitude: {abs(T[3]).max()}\n'

            log.info(t_amplitude_values_string)
            # -------------------------------------------------------------------------------------
            z_amplitude_values_string = f'max constant z_0   amplitude: {abs(Z[0]).max()}\n'

            if self.Z_truncation_order >= 1:
                z_amplitude_values_string += f'max single   z_i   amplitude: {abs(Z[1]).max()}\n'
            if self.Z_truncation_order >= 2:
                z_amplitude_values_string += f'max double   z_ij  amplitude: {abs(Z[2]).max()}\n'
            if self.Z_truncation_order >= 3:
                z_amplitude_values_string += f'max triple   z_ijk amplitude: {abs(Z[3]).max()}\n'

            log.info(z_amplitude_values_string)

            log.info(f'Norm of the wavefunction:\n{self.Norm[-1][1]}')

            if self.calculate_population_flag:
                log.info(f'Diabatic  State population:\n {self.State_pop_DB[-1][1]}')
                log.info(f'Adiabatic State population:\n {self.State_pop_AB[-1][1]}')

            self.last_print = fs
            self.last_counter = 0

        return

    def _compute_t_residual(self, H_bar_tilde, Z, Z_conj):
        """ compute t_residue from the net t_residue (Ehrenfest dynamics) """

        weight = np.einsum('y,y->', Z_conj[0], Z[0])

        # constant t residue
        dT = {
            0: complex(0.0),
            1: np.zeros(self.N, dtype=complex),
            2: np.zeros((self.N, self.N), dtype=complex)
        }

        # single t residue
        if self.T_truncation_order >= 1:
            dT[1] += np.einsum('y,yxi,x->i', Z_conj[0], H_bar_tilde[(1, 0)], Z[0]) / weight

        # double t residual_equation
        if self.T_truncation_order >= 2:
            dT[2] += np.einsum('y,yxij,x->ij', Z_conj[0], H_bar_tilde[(2, 0)], Z[0]) / weight

        return dT

    def _compute_t_residual_new(self, H_bar_tilde, C):
        """compute t from Ehrenfest parameterization using the new scheme (weight C)"""

        C_0_conj = np.conj(C[(0, 0)])
        C_0 = C[(0, 0)]
        weight = np.einsum('y,y->', C_0_conj, C_0)

        # constant t residue
        dT = {
            0: complex(0.0),
            1: np.zeros(self.N, dtype=complex),
            2: np.zeros((self.N, self.N), dtype=complex)
        }

        # single t residue
        if self.T_truncation_order >= 1:
            dT[1] += np.einsum('y,yxi,x->i', C_0_conj, H_bar_tilde[(1, 0)], C_0) / weight

        # double t residual_equation
        if self.T_truncation_order >= 2:
            dT[2] += np.einsum('y,yxij,x->ij', C_0_conj, H_bar_tilde[(2, 0)], C_0) / weight

        return dT

    def _compute_z_residual_new_scheme(self, R, Z, T_conj, dT, C):
        """ compute z compute z residue by subtracting t residue from net residue
            adopt the new scheme the introduce similarity transform
        """
        A, N = self.A, self.N

        T_conj_1 = T_conj[(0, 1)]

        def _sim_trans_dT(T_conj_1, dT):
            """ similarity transform dT"""
            output_tensor = {
                (0, 0): np.einsum('k,k->', T_conj_1, dT[1]),
                (1, 0): dT[1]
            }
            return output_tensor

        def _cal_D_0(T_conj_1, dz_1, dz_2, dz_3):
            """calculate constant D_0"""
            D_0 = np.zeros(self.A, dtype=complex)

            if not self.op_einsum_flag:
                D_0 += np.einsum('k,yk->y', T_conj_1, dz_1)
                D_0 += 0.5 * np.einsum('k,l,ykl->y', T_conj_1, T_conj_1, dz_2)
                D_0 += (1.0 / 6.0) * np.einsum('k,l,m,yklm->y', T_conj_1, T_conj_1, T_conj_1, dz_3)
            else:
                optimized_einsum = iter(self.d0_opt_paths)
                D_0 += next(optimized_einsum)(T_conj_1, dz_1)
                D_0 += 0.5 * next(optimized_einsum)(T_conj_1, T_conj_1, dz_2)
                D_0 += (1.0 / 6.0) * next(optimized_einsum)(T_conj_1, T_conj_1, T_conj_1, dz_3)

            return D_0

        def _cal_D_1(T_conj_1, dz_2, dz_3):
            """calculate single D_1"""
            D_1 = np.zeros([self.A, self.N], dtype=complex)

            if not self.op_einsum_flag:
                D_1 += np.einsum('k,yki->yi', T_conj_1, dz_2)
                D_1 += 0.5 * np.einsum('k,l,ykli->yi', T_conj_1, T_conj_1, dz_3)
            else:
                optimized_einsum = iter(self.d1_opt_paths)
                D_1 += next(optimized_einsum)(T_conj_1, dz_2)
                D_1 += 0.5 * next(optimized_einsum)(T_conj_1, T_conj_1, dz_3)

            return D_1

        def _cal_D_2(T_conj_1, dz_3):
            """calculate double D_2"""
            D_2 = np.zeros([self.A, self.N, self.N], dtype=complex)

            if not self.op_einsum_flag:
                D_2 += 0.5 * np.einsum('l,ylij->yij', T_conj_1, dz_3)
            else:
                optimized_einsum = iter(self.d2_opt_paths)
                D_2 += 0.5 * next(optimized_einsum)(T_conj_1, dz_3)

            return D_2

        def _cal_dT_Tran_C(dT_trans, C):
            """calculation dT_bar * C"""

            C_0, C_1, C_2, C_3 = C[(0, 0)], C[(1, 0)], C[(2, 0)], C[(3, 0)]

            dT_trans_0, dT_trans_1 = dT_trans[(0, 0)], dT_trans[(1, 0)]

            output_tensor = {
                (0, 0): dT_trans[(0, 0)] * C_0
            }

            if self.Z_truncation_order >= 1:
                output_tensor[(1, 0)] = np.einsum('i,y->yi', dT_trans_1, C_0)
                output_tensor[(1, 0)] += dT_trans_0 * C_1

            if self.Z_truncation_order >= 2:
                output_tensor[(2, 0)] = np.einsum('i,yj->yij', dT_trans_1, C_1)
                output_tensor[(2, 0)] += 0.5 * dT_trans_0 * C_2

            if self.Z_truncation_order >= 3:
                output_tensor[(3, 0)] = 0.5 * np.einsum('i,yjk->yijk', dT_trans_1, C_2)
                output_tensor[(3, 0)] += (1.0 / 6.0) * dT_trans_0 * C_3

            return output_tensor

        # similarity transform dT
        dT_bar = _sim_trans_dT(T_conj_1, dT)

        # calculate idT_bar * C
        X = _cal_dT_Tran_C(dT_bar, C)

        # calculate dz_3
        dz_3 = np.zeros([A, N, N, N], dtype=complex)
        if self.Z_truncation_order >= 3:
            dz_3 += R[3] - X[(3, 0)]
            dz_3 = symmetrize_tensor(self.N, dz_3, order=3)

        # calculate dz_2
        dz_2 = np.zeros([A, N, N], dtype=complex)
        if self.Z_truncation_order >= 2:
            D_2 = _cal_D_2(T_conj_1, dz_3)
            dz_2 += R[2]
            dz_2 -= X[(2, 0)]
            dz_2 -= D_2
            dz_2 = symmetrize_tensor(self.N, dz_2, order=2)

        # calculate dz_1
        dz_1 = np.zeros([A, N], dtype=complex)
        if self.Z_truncation_order >= 2:
            D_1 = _cal_D_1(T_conj_1, dz_2, dz_3)
            dz_1 += R[1].copy()
            dz_1 -= X[(1, 0)]
            dz_1 -= D_1

        # calculate dz_0
        D_0 = _cal_D_0(T_conj_1, dz_1, dz_2, dz_3)
        dz_0 = R[0].copy()
        dz_0 -= X[(0, 0)]
        dz_0 -= D_0

        output_tensor = {
            0: dz_0,
            1: dz_1,
            2: dz_2,
            3: dz_3
        }

        return output_tensor

    def _compute_z_residual(self, R, Z, T_conj, dT, z_opt_flag=False):
        """ compute z residue by subtracting t residue from net residue """

        def compute_z_0_residual(Z, T_conj, X, dT, dz_3, dz_2, dz_1, optimized_einsum_version=False):
            """ constant Z residue """

            dz_0 = np.zeros(self.A, dtype=complex)
            dz_0 += R[0] - (X * Z[0])

            dz_0 -= X * np.einsum('l,yl->y', T_conj[1], Z[1])
            dz_0 -= np.einsum('k,yk->y', T_conj[1], dz_1)

            if self.Z_truncation_order >= 2:
                dz_0 -= X * 0.5 * np.einsum('l,m,ylm->y', T_conj[1], T_conj[1], Z[2])
                dz_0 -= 0.5 * np.einsum('k,l,ykl->y', T_conj[1], T_conj[1], dz_2)

            if self.Z_truncation_order >= 3:
                dz_0 -= X * (1.0 / 6.0) * np.einsum('k,l,m,yklm->y', T_conj[1], T_conj[1], T_conj[1], Z[3])
                dz_0 -= (1.0 / 6.0) * np.einsum('k,l,m,yklm->y', T_conj[1], T_conj[1], T_conj[1], dz_3)

            if self.T_truncation_order >= 2:
                dz_0 -= 0.5 * np.einsum('k,l,kl,y->y', T_conj[1], T_conj[1], dT[2], Z[0])

                if self.Z_truncation_order >= 1:
                    dz_0 -= 0.5 * np.einsum('k,l,m,lm,yk->y', T_conj[1], T_conj[1], T_conj[1], dT[2], Z[1])

                if self.Z_truncation_order >= 2:
                    dz_0 -= 0.25 * np.einsum('k,l,m,n,kl,ymn->y', T_conj[1], T_conj[1], T_conj[1], T_conj[1], dT[2], Z[2])

            return dz_0

        def compute_z_1_residual(Z, T_conj, X, dT, dz_3, dz_2, optimized_einsum_version=False):
            """ single Z residue """

            if self.Z_truncation_order < 1:
                return np.zeros([self.A, self.N], dtype=complex)

            if self.Z_truncation_order >= 1:
                dz_i = np.zeros([self.A, self.N], dtype=complex)
                dz_i += R[1] - (X * Z[1])
            else:
                raise Exception('')

            dz_i -= np.einsum('i,y->yi', dT[1], Z[0])
            dz_i -= np.einsum('k,i,yk->yi', T_conj[1], dT[1], Z[1])

            if self.Z_truncation_order >= 2:
                dz_i -= X * np.einsum('l,yli->yi', T_conj[1], Z[2])
                dz_i -= np.einsum('k,yki->yi', T_conj[1], dz_2)
                dz_i -= 0.5 * np.einsum('k,l,i,ykl->yi', T_conj[1], T_conj[1], dT[1], Z[2])

            if self.Z_truncation_order >= 3:
                dz_i -= 0.5 * X * np.einsum('k,l,ykli->yi', T_conj[1], T_conj[1], Z[3])
                dz_i -= 0.5 * np.einsum('k,l,ykli->yi', T_conj[1], T_conj[1], dz_3)
                dz_i -= (1.0 / 6.0) * np.einsum('k,l,m,i,yklm->yi', T_conj[1], T_conj[1], T_conj[1], dT[1], Z[3])

            if self.T_truncation_order >= 2:
                dz_i -= np.einsum('k,ki,y->yi', T_conj[1], dT[2], Z[0])

                if self.Z_truncation_order >= 1:
                    dz_i -= np.einsum('k,l,ki,yl->yi', T_conj[1], T_conj[1], dT[2], Z[1])
                    dz_i -= 0.5 * np.einsum('k,l,kl,yi->yi', T_conj[1], T_conj[1], dT[2], Z[1])

                if self.Z_truncation_order >= 2:
                    dz_i -= 0.5 * np.einsum('k,l,m,mi,ykl->yi', T_conj[1], T_conj[1], T_conj[1], dT[2], Z[2])
                    dz_i -= 0.5 * np.einsum('k,l,m,lm,yki->yi', T_conj[1], T_conj[1], T_conj[1], dT[2], Z[2])

                if self.Z_truncation_order >= 3:
                    dz_i -= 0.25 * np.einsum('k,l,m,n,mn,ykli->yi', T_conj[1], T_conj[1], T_conj[1], T_conj[1], dT[2], Z[3])
                    dz_i -= (1.0 / 6.0) * np.einsum('k,l,m,n,ni,yklm->yi', T_conj[1], T_conj[1], T_conj[1], T_conj[1], dT[2], Z[3])

            return dz_i

        def compute_z_2_residual(Z, T_conj, X, dT, dz_3, optimized_einsum_version=False):
            """ double Z residue """

            if self.Z_truncation_order < 2:
                return np.zeros([self.A, self.N, self.N], dtype=complex)

            if self.Z_truncation_order >= 2:
                dz_ij = np.zeros([self.A, self.N, self.N], dtype=complex)
                # we need to multiply by 1/2 because symmetrization doubles the value
                # and we can't divide by 1/2 after symmetrization
                dz_ij += (R[2] - 0.5 * (X * Z[2]))
            else:
                raise Exception('')

            dz_ij -= np.einsum('i,yj->yij', dT[1], Z[1])

            if self.Z_truncation_order >= 2:
                dz_ij -= np.einsum('k,i,ykj->yij', T_conj[1], dT[1], Z[2])

            if self.Z_truncation_order >= 3:
                # 1/2 because it is already symmetric on dz_3
                dz_ij -= 0.5 * np.einsum('l,ylij->yij', T_conj[1], dz_3)
                # 1/2 because of prefactor
                dz_ij -= 0.5 * np.einsum('k,l,i,yklj->yij', T_conj[1], T_conj[1], dT[1], Z[3])
                # 1/2 because it is already symmetric on Z[3]
                dz_ij -= 0.5 * np.einsum('k,l,l,ykij->yij', T_conj[1], T_conj[1], dT[1], Z[3])

            if self.T_truncation_order >= 2:
                dz_ij -= np.einsum('k,ki,yj->yij', T_conj[1], dT[2], Z[1])
                # 1/2 because it is already symmetric on dT[2]
                dz_ij -= 0.5 * np.einsum('ij,y->yij', dT[2], Z[0])
                # 1/2 because it is already symmetric on dT[2]
                dz_ij -= 0.5 * np.einsum('k,ij,yk->yij', T_conj[1], dT[2], Z[1])

                if self.Z_truncation_order >= 2:
                    dz_ij -= np.einsum('k,l,ki,ylj->yij', T_conj[1], T_conj[1], dT[2], Z[2])
                    # 1/4 because of 1/2 prefactor and 1/2 from already symmetric on Z[2]
                    dz_ij -= 0.25 * np.einsum('k,l,kl,yij->yij', T_conj[1], T_conj[1], dT[2], Z[2])
                    # 1/4 because of 1/2 prefactor and 1/2 from already symmetric on dT[2]
                    dz_ij -= 0.25 * np.einsum('k,l,ij,ykl->yij', T_conj[1], T_conj[1], dT[2], Z[2])

                if self.Z_truncation_order >= 3:
                    dz_ij -= 0.5 * np.einsum('k,l,m,mi,yklj->yij', T_conj[1], T_conj[1], T_conj[1], dT[2], Z[3])
                    # 1/4 because of 1/2 prefactor and 1/2 from already symmetric on Z[3]
                    dz_ij -= 0.25 * np.einsum('k,l,m,lm,ykij->yij', T_conj[1], T_conj[1], T_conj[1], dT[2], Z[3])
                    # 1/12 because of 1/6 prefactor and 1/2 from already symmetric on dT[2]
                    dz_ij -= (1.0 / 12.0) * np.einsum('k,l,m,ij,yklm->yij', T_conj[1], T_conj[1], T_conj[1], dT[2], Z[3])

            # symmetrize
            dz_ij = symmetrize_tensor(self.N, dz_ij, order=2)

            return dz_ij

        def compute_z_3_residual(Z, T_conj, X, dT, optimized_einsum_version=False):
            """ triple Z residue """

            if self.Z_truncation_order < 3:
                return np.zeros([self.A, self.N, self.N, self.N], dtype=complex)

            if self.Z_truncation_order >= 3:
                dz_ijk = np.zeros([self.A, self.N, self.N, self.N], dtype=complex)
                # 1/6 because 1/3! from already symmetric on R[3] and Z[3]
                dz_ijk += (R[3] - (1.0 / 6.0) * (X * Z[3]))
            else:
                raise Exception('')

            if self.Z_truncation_order >= 2:
                # 1/2 because of 1/2 from already symmetric on Z[2]
                dz_ijk -= 0.5 * np.einsum('i,yjk->yijk', dT[1], Z[2])

            if self.Z_truncation_order >= 3:
                # 1/2 because of 1/2 from already symmetric on Z[3]
                dz_ijk -= 0.5 * np.einsum('l,i,yljk->yijk', T_conj[1], dT[1], Z[3])

            if self.T_truncation_order >= 2:
                # 1/2 because of 1/2 from already symmetric on dT[2]
                dz_ijk -= 0.5 * np.einsum('ij,yk->yijk', dT[2], Z[1])

                if self.Z_truncation_order >= 2:
                    # 1/2 because of 1/2 from already symmetric on Z[2]
                    dz_ijk -= 0.5 * np.einsum('l,li,yjk->yijk', T_conj[1], dT[2], Z[2])
                    # 1/2 because of 1/2 from already symmetric on dT[2]
                    dz_ijk -= 0.5 * np.einsum('l,ij,ylk->yijk', T_conj[1], dT[2], Z[2])

                if self.Z_truncation_order >= 3:
                    # 1/2 because of 1/2 from already symmetric on Z[3]
                    dz_ijk -= 0.5 * np.einsum('l,m,mi,yljk->yijk', T_conj[1], T_conj[1], dT[2], Z[3])
                    # 1/4 because of 1/2 prefactor and 1/2 from already symmetric on dT[2]
                    dz_ijk -= 0.25 * np.einsum('l,m,ij,ylmk->yijk', T_conj[1], T_conj[1], dT[2], Z[3])
                    # 1/12 because of 1/2 prefactor and 1/3! from already symmetric on Z[3]
                    dz_ijk -= (1.0 / 12.0) * np.einsum('l,m,lm,yijk->yijk', T_conj[1], T_conj[1], dT[2], Z[3])

            # symmetrize
            dz_ijk = symmetrize_tensor(self.N, dz_ijk, order=3)

            return dz_ijk

        # compute the common factor (t^k)^*(idt^k/dtau)
        X = np.einsum('k,k->', T_conj[1], dT[1])

        # they have to be in this order as Z1 depends on Z2 and Z0 depends on Z1
        dz_3 = compute_z_3_residual(Z, T_conj, X, dT, z_opt_flag)
        dz_2 = compute_z_2_residual(Z, T_conj, X, dT, dz_3, z_opt_flag)
        dz_1 = compute_z_1_residual(Z, T_conj, X, dT, dz_3, dz_2, z_opt_flag)
        dz_0 = compute_z_0_residual(Z, T_conj, X, dT, dz_3, dz_2, dz_1, z_opt_flag)

        dZ = {
            3: dz_3,
            2: dz_2,
            1: dz_1,
            0: dz_0,
        }

        return dZ

    def _f_z_0(self, H_bar, Z, T_conj, optimized_einsum_version=False):
        """ update G_0 and return dz_0/dtau """

        # initialize G_0
        R_y = np.zeros_like(Z[0], dtype=complex)  # similarity transformed G_0

        # initialize
        R_y += np.einsum('yx,x->y', H_bar[(0, 0)], Z[0])
        R_y += np.einsum('yxk,xk->y', H_bar[(0, 1)], Z[1])

        if self.Z_truncation_order >= 2:
            R_y += 0.5 * np.einsum('yxkl,xkl->y', H_bar[(0, 2)], Z[2])

        # adding terms associated with e^T_dagger
        R_y += np.einsum('k,yx,xk->y', T_conj[1], H_bar[(0, 0)], Z[1])
        R_y += np.einsum('k,yxk,x->y', T_conj[1], H_bar[(1, 0)], Z[0])
        R_y += np.einsum('k,yxlk,xl->y', T_conj[1], H_bar[(1, 1)], Z[1])
        R_y += np.einsum('k,l,yxl,xk->y', T_conj[1], T_conj[1], H_bar[(1, 0)], Z[1])
        R_y += 0.5 * np.einsum('k,l,yxkl,x->y', T_conj[1], T_conj[1], H_bar[(2, 0)], Z[0])
        R_y += 0.5 * np.einsum('k,l,m,yxlm,xk->y', T_conj[1], T_conj[1], T_conj[1], H_bar[(2, 0)], Z[1])

        if self.Z_truncation_order >= 2:
            R_y += np.einsum('k,yxl,xkl->y', T_conj[1], H_bar[(0, 1)], Z[2])
            R_y += 0.5 * np.einsum('k,l,yx,xkl->y', T_conj[1], T_conj[1], H_bar[(0, 0)], Z[2])
            R_y += np.einsum('k,l,yxml,xkm->y', T_conj[1], T_conj[1], H_bar[(1, 1)], Z[2])
            R_y += 0.5 * np.einsum('k,l,m,yxm,xkl->y', T_conj[1], T_conj[1], T_conj[1], H_bar[(1, 0)], Z[2])
            R_y += 0.25 * np.einsum('i,j,k,l,yxij,xkl->y', T_conj[1], T_conj[1], T_conj[1], T_conj[1], H_bar[(2, 0)], Z[2])

        return R_y

    def _f_z_i(self, H_bar, Z, T_conj, optimized_einsum_version=False):
        """ update G_i and return dz_i/dtau """

        # initialize as zeros
        R_y = np.zeros_like(Z[1], dtype=complex)

        R_y += np.einsum('yxi,x->yi', H_bar[(1, 0)], Z[0])
        R_y += np.einsum('yxki,xk->yi', H_bar[(1, 1)], Z[1])
        R_y += np.einsum('yx,xi->yi', H_bar[(0, 0)], Z[1])
        # adding terms associated e^T_dagger
        R_y += np.einsum('k,yxi,xk->yi', T_conj[1], H_bar[(1, 0)], Z[1])
        R_y += np.einsum('k,yxk,xi->yi', T_conj[1], H_bar[(1, 0)], Z[1])
        R_y += np.einsum('k,yxki,x->yi', T_conj[1], H_bar[(2, 0)], Z[0])
        R_y += np.einsum('k,l,yxli,xk->yi', T_conj[1], T_conj[1], H_bar[(2, 0)], Z[1])
        R_y += 0.5 * np.einsum('k,l,yxkl,xi->yi', T_conj[1], T_conj[1], H_bar[(2, 0)], Z[1])

        if self.Z_truncation_order >= 2:
            R_y += np.einsum('yxk,xki->yi', H_bar[(0, 1)], Z[2])
            R_y += np.einsum('k,yx,xki->yi', T_conj[1], H_bar[(0, 0)], Z[2])
            R_y += np.einsum('k,yxli,xkl->yi', T_conj[1], H_bar[(1, 1)], Z[2])
            R_y += np.einsum('k,yxlk,xli->yi', T_conj[1], H_bar[(1, 1)], Z[2])
            R_y += np.einsum('k,l,yxl,xki->yi', T_conj[1], T_conj[1], H_bar[(1, 0)], Z[2])
            R_y += 0.5 * np.einsum('k,l,yxi,xkl->yi', T_conj[1], T_conj[1], H_bar[(1, 0)], Z[2])
            R_y += 0.5 * np.einsum('k,l,m,yxmi,xkl->yi', T_conj[1], T_conj[1], T_conj[1], H_bar[(2, 0)], Z[2])
            R_y += 0.5 * np.einsum('k,l,m,yxlm,xki->yi', T_conj[1], T_conj[1], T_conj[1], H_bar[(2, 0)], Z[2])

        return R_y

    def _f_z_ij(self, H_bar, Z, T_conj, optimized_einsum_version=False):
        # return G_ij and dz_ij/dtau

        # initialize as zeros
        R_y = np.zeros_like(Z[2], dtype=complex)

        R_y += 0.5 * np.einsum('yx,xij->yij', H_bar[(0, 0)], Z[2])
        R_y += np.einsum('yxi,xj->yij', H_bar[(1, 0)], Z[1])
        R_y += 0.5 * np.einsum('yxij,x->yij', H_bar[(2, 0)], Z[0])
        R_y += np.einsum('yxki,xkj->yij', H_bar[(1, 1)], Z[2])
        # adding terms associated e^T_dagger
        R_y += np.einsum('k,yxi,xkj->yij', T_conj[1], H_bar[(1, 0)], Z[2])
        R_y += 0.5 * np.einsum('k,yxk,xij->yij', T_conj[1], H_bar[(1, 0)], Z[2])
        R_y += 0.5 * np.einsum('k,yxij,xk->yij', T_conj[1], H_bar[(2, 0)], Z[1])
        R_y += np.einsum('k,yxki,xj->yij', T_conj[1], H_bar[(2, 0)], Z[1])
        R_y += np.einsum('k,l,yxli,xkj->yij', T_conj[1], T_conj[1], H_bar[(2, 0)], Z[2])
        R_y += 0.25 * np.einsum('k,l,yxij,xkl->yij', T_conj[1], T_conj[1], H_bar[(2, 0)], Z[2])
        R_y += 0.25 * np.einsum('k,l,yxkl,xij->yij', T_conj[1], T_conj[1], H_bar[(2, 0)], Z[2])

        return R_y

    def calculate_ACF_ABS(self, E_tdm, U, b, cross=False):
        """ compute ACF from the CC amplitudes """

        if not cross:  # abandon cross correlation
            new_U = np.zeros_like(U)
            new_U[b] = U[b]
            U = new_U

        return np.einsum('ab,b,a ->', E_tdm, U, E_tdm[:, b])

    def calculate_ACF_ECD(self, M_tdm, E_tdm, U, b, cross=False):
        """ compute ACF from the CC amplitudes """

        if not cross:  # abandon cross correlation
            new_U = np.zeros_like(U)
            new_U[b] = U[b]
            U = new_U

        return np.einsum('ab,b,a ->', M_tdm, U, E_tdm[:, b])

    def calculate_CC_norm(self, T, Z, b):
        """ Calculate norm of the wavefunction from CC """

        A, N = self.A, self.N

        # pre-define as we use it multiple times
        T_conj = {1: np.conj(T[1])}
        if self.T_truncation_order >= 2:
            T_conj[2] = np.conj(T[2])

        # initialize C_args
        C_args = [np.zeros(A, dtype=complex) for i in range(4)]

        # compute C_n,x
        def Cal_C_0(C_0):
            """ Compute 0th order contribution """
            C_0 += Z[0]

            if self.Z_truncation_order >= 1:
                C_0 += np.einsum('k,xk->x', T_conj[1], Z[1])

            if self.Z_truncation_order >= 2:
                C_0 += 0.5 * np.einsum('k,l,xkl->x', T_conj[1], T_conj[1], Z[2])

                if self.T_truncation_order >= 2:
                    C_0 += 0.5 * np.einsum('kl,xkl->x', T_conj[2], Z[2])
            if self.Z_truncation_order >= 3:
                C_0 += 1./6. * np.einsum('k,l,m,xlkm->x', T_conj[1], T_conj[1], T_conj[1], Z[3])
            return

        def Cal_C_1(C_1):
            """ Compute 1st order contribution """

            if self.Z_truncation_order >= 1:
                C_1 += np.einsum('xi->x', Z[1])

            if self.Z_truncation_order >= 2:
                C_1 += np.einsum('k,xki->x', T_conj[1], Z[2])

            if self.Z_truncation_order >= 3:
                C_1 += 0.5 * np.einsum('k,l,xkli->x', T_conj[1], T_conj[1], Z[3])
            return

        def Cal_C_2(C_2, Identity):
            """ Compute 2nd order contribution """
            Q = np.ones([N, N], dtype=complex) - Identity

            if self.Z_truncation_order >= 2:
                C_2 += 0.5 * np.einsum('xij,ij->x', Z[2], Q)
                C_2 += 1. / np.sqrt(2) * np.einsum('xij,ij->x',  Z[2], Identity)

            if self.Z_truncation_order >= 3:
                C_2 += 0.5 * np.einsum('k,xkij,ij->x', T_conj[1], Z[3], Q)
                C_2 += 1. / np.sqrt(2) * np.einsum('k,xkij,ij->x', T_conj[1], Z[3], Identity)

        def Cal_C_3(C_3, Identity):
            """ Compute 3rd order contribution"""
            Q = np.ones([N, N], dtype=complex) - Identity
            if self.Z_truncation_order >= 3:
                C_3 += 1. / 6. * np.einsum('xijk,ij,jk,ik->x', Z[3], Q, Q, Q)
                C_3 += 1. / np.sqrt(6) * np.einsum('xijk,ij,jk->x', Z[3], Identity, Identity)
                C_3 += 1. / 6. * np.sqrt(2) * np.einsum('xijk,ij,jk->x', Z[3], Identity, Q)
                C_3 += 1. / 6. * np.sqrt(2) * np.einsum('xijk,ij,jk->x', Z[3], Q, Identity)
                C_3 += 1. / 6. * np.sqrt(2) * np.einsum('xijk,jk,ik->x', Z[3], Q, Identity)

            return

        Cal_C_0(C_args[0])
        Cal_C_1(C_args[1])

        if self.Z_truncation_order >= 2:
            Identity = np.eye(N, dtype=complex)
            Cal_C_2(C_args[2], Identity)
        if self.Z_truncation_order >= 3:
            Cal_C_3(C_args[3], Identity)

        s_0 = np.einsum('i,i->', T_conj[1], T[1])

        norm = 0.

        norm += np.sum(np.abs(C_args[0])**2)
        norm += np.sum(np.abs(C_args[1])**2)

        if self.Z_truncation_order >= 2:
            norm += np.sum(np.abs(C_args[2])**2)

        if self.Z_truncation_order >= 3:
            norm += np.sum(np.abs(C_args[3])**2)

        norm *= np.exp(s_0 + (2 * T[0].real))

        return norm, C_args, s_0

    def calculate_state_populations(self, s_0, norm, C_arg):
        """compute state population associated with a specific surface"""

        # unpack input tensors
        C_0, C_1, C_2, C_3 = C_arg

        # calculate density matrix
        D = np.zeros((self.A, self.A), dtype=complex)
        D += (
            np.einsum('p,q->pq', np.conjugate(C_0), C_0) +
            np.einsum('p,q->pq', np.conjugate(C_1), C_1) +
            np.einsum('p,q->pq', np.conjugate(C_2), C_2) +
            np.einsum('p,q->pq', np.conjugate(C_3), C_3)
        )
        D *= np.exp(s_0)

        # diagonalize density
        E, V = np.linalg.eigh(D)

        # diabatic state population
        state_pop_DB = np.array(np.diag(D), dtype=complex)

        # adiabatic state population
        state_pop_AB = np.array(E, dtype=complex)

        # normalize the populations
        state_pop_DB /= norm
        state_pop_AB /= norm

        return state_pop_DB, state_pop_AB

    def calculate_CC_norm_new(self, C_args):
        """calculate state population from C tensors directly """
        # calculate norm
        norm = 0.0
        norm += np.einsum('x,x->', np.conj(C_args[(0, 0)]), C_args[(0, 0)])
        norm += np.einsum('xi,xi->', np.conj(C_args[(1, 0)]), C_args[(1, 0)])

        if self.Z_truncation_order >= 2:
            norm += 0.5 * np.einsum('xij,xij->', np.conj(C_args[(2, 0)]), C_args[(2, 0)])
        if self.Z_truncation_order >= 3:
            norm += 1. / 6. * np.einsum('xijk,xijk->', np.conj(C_args[(3, 0)]), C_args[(3, 0)])

        return norm

    def calculate_density_matrix_new(self, norm, C_args):
        """calculate denisty matrix"""
        d_matrix = np.zeros((self.A, self.A), dtype=complex)
        d_matrix += np.einsum('x,y->xy', np.conj(C_args[(0, 0)]), C_args[(0, 0)])
        d_matrix += np.einsum('xi,yi->xy', np.conj(C_args[(1, 0)]), C_args[(1, 0)])
        if self.Z_truncation_order >= 2:
            d_matrix += 0.5 * np.einsum('xij,yij->xy', np.conj(C_args[(2, 0)]), C_args[(2, 0)])
        if self.Z_truncation_order >= 3:
            d_matrix += 1. / 6. * np.einsum('xijk,yijk->xy', np.conj(C_args[(3, 0)]), C_args[(3, 0)])

        d_matrix /= norm  # normalize density matrix

        return d_matrix

    def calculate_state_populations_new(self, density_matrix):
        """calculate state population from C tensors directly """

        # diagonalize density
        E, V = np.linalg.eigh(density_matrix)

        # diabatic state population
        state_pop_DB = np.array(np.diag(density_matrix), dtype=complex)

        # adiabatic state population
        state_pop_AB = np.array(E, dtype=complex)

        return state_pop_DB, state_pop_AB

    def _cal_H_bar_tilde(self, input_tensor, T_conj, opt_flag=False):
        """calculate the second similarity transform the the Hamiltonian"""
        A, N = self.A, self.N

        def f_s_0(O_mat):
            """return constant residue"""
            # initialize as zero
            R = np.zeros([A, A], dtype=complex)

            R += O_mat[(0, 0)]
            if opt_flag:
                optimized_einsum = iter(self.H_bar_tilde_paths)  # first two
                R += next(optimized_einsum)(T_conj[(0, 1)], O_mat[(1, 0)])
                R += 0.5 * next(optimized_einsum)(T_conj[(0, 1)], T_conj[(0, 1)], O_mat[(2, 0)])
            else:  # old
                R += np.einsum('k,abk->ab', T_conj[(0, 1)], O_mat[(1, 0)])
                R += 0.5 * np.einsum('k,l,abkl->ab', T_conj[(0, 1)], T_conj[(0, 1)], O_mat[(2, 0)])

            return R

        def f_s_I(O_mat):
            """return residue R_I"""
            R = np.zeros([A, A, N], dtype=complex)

            R += O_mat[(1, 0)]

            if opt_flag:
                optimized_einsum = self.H_bar_tilde_paths[2]
                R += optimized_einsum(T_conj[(0, 1)], O_mat[(2, 0)])
            else:  # old
                R += np.einsum('k,abik->abi', T_conj[(0, 1)], O_mat[(2, 0)])
            return R

        def f_s_i(O_mat):
            """return residue R_i"""
            R = np.zeros([A, A, N], dtype=complex)

            R += O_mat[(0, 1)]
            if opt_flag:
                optimized_einsum = self.H_bar_tilde_paths[3]
                R += optimized_einsum(T_conj[(0, 1)], O_mat[(1, 1)])
            else:  # old
                R += np.einsum('k,abki->abi', T_conj[(0, 1)], O_mat[(1, 1)])
            return R

        def f_s_Ij(O_mat):
            """return residue R_Ij"""
            return O_mat[(1, 1)]

        output_tensor = {
                (0, 0): f_s_0(input_tensor),
                (1, 0): f_s_I(input_tensor),
                (0, 1): f_s_i(input_tensor),
                (1, 1): f_s_Ij(input_tensor),
                (2, 0): input_tensor[(2, 0)],
                (0, 2): input_tensor[(0, 2)],
        }

        return output_tensor

    def _compute_C_matrix(self, input_tensor, T_conj, opt_flag=False):
        """"calculate intermediate quantity C"""
        A, N = self.A, self.N

        def f_s_0(O_mat):
            """return constant residue"""
            # initialize as zero
            R = np.zeros([A], dtype=complex)

            R += O_mat[(0, 0)]

            if opt_flag:
                optimized_einsum = self.Cmat_Z0_opt_paths[0]
                R += optimized_einsum(T_conj[(0, 1)], O_mat[(1, 0)])
            else:
                R += np.einsum('k,xk->x', T_conj[(0, 1)], O_mat[(1, 0)])

            if self.Z_truncation_order >= 2:
                if opt_flag:
                    optimized_einsum = self.Cmat_Z2_opt_paths[0]
                    R += optimized_einsum(T_conj[(0, 1)], T_conj[(0, 1)], O_mat[(2, 0)])
                else:
                    R += 0.5 * np.einsum('k,l,xkl->x', T_conj[(0, 1)], T_conj[(0, 1)], O_mat[(2, 0)])
            if self.Z_truncation_order >= 3:
                if opt_flag:
                    optimized_einsum = self.Cmat_Z3_opt_paths[0]
                    R += optimized_einsum(T_conj[(0, 1)], T_conj[(0, 1)], T_conj[(0, 1)], O_mat[(3, 0)])
                else:
                    R += 1./6. * np.einsum('k,l,m,xklm->x', T_conj[(0, 1)], T_conj[(0, 1)], T_conj[(0, 1)], O_mat[(3, 0)])

            return R

        def f_s_1(O_mat):
            """return single residue"""
            R = np.zeros([A, N], dtype=complex)

            R += O_mat[(1, 0)]

            if self.Z_truncation_order >= 2:
                if opt_flag:
                    optimized_einsum = self.Cmat_Z2_opt_paths[1]
                    R += optimized_einsum(T_conj[(0, 1)], O_mat[(2, 0)])
                else:
                    R += np.einsum('k,xik->xi', T_conj[(0, 1)], O_mat[(2, 0)])

            if self.Z_truncation_order >= 3:
                if opt_flag:
                    optimized_einsum = self.Cmat_Z3_opt_paths[1]
                    R += 0.5 * optimized_einsum(T_conj[(0, 1)], T_conj[(0, 1)], O_mat[(3, 0)])
                else:
                    R += 0.5 * np.einsum('k,l,xikl->xi', T_conj[(0, 1)], T_conj[(0, 1)], O_mat[(3, 0)])

            return R

        def f_s_2(O_mat):
            """return double residue"""
            R = np.zeros([A, N, N], dtype=complex)

            if self.Z_truncation_order >= 2:
                # no optimization for simple addition
                R += O_mat[(2, 0)]

            if self.Z_truncation_order >= 3:
                if opt_flag:
                    optimized_einsum = self.Cmat_Z3_opt_paths[2]
                    R += optimized_einsum(T_conj[(0, 1)], O_mat[(3, 0)])
                else:
                    R += np.einsum('k, xkij->xij', T_conj[(0, 1)], O_mat[(3, 0)])

            return R

        def f_s_3(O_mat):
            """return triple residue"""
            return O_mat[(3, 0)]

        output_tensor = {
            (0, 0): f_s_0(input_tensor),
            (1, 0): f_s_1(input_tensor),
            (2, 0): None,
            (3, 0): None
        }

        if self.Z_truncation_order >= 2:
            output_tensor[(2, 0)] = f_s_2(input_tensor)

        if self.Z_truncation_order >= 3:
            output_tensor[(3, 0)] = f_s_3(input_tensor)

        return output_tensor

    def rk45_solve_ivp_integration_function(self, time, y_tensor, t_final):
        """ Integration function used by `solve_ivp` integrator inside `rk45_integration` method.

        `time` is a float, the value of time for the current integration step
        `y_tensor` is an (n, k) dimensional tensor where the n dimension counts the ode's
        that we are attempting to solve and k can represent multiple time steps to block integrate over
        at the moment we do not do any block integration so k is 1
        """

        def solve_z(H_bar, Z, T, b, selected_surface, opt_flag=self.op_einsum_flag):
            """ Solve equations for Z amplitude
            conduct similarity transformation for the Hamiltonian over e^z
            """

            # create conjugation dictionaries
            Z_conj_dict = {
                0: np.conj(Z[0]),
            }
            T_conj_dict = {
                # 0: 0  # t_0 IS zero
                1: np.conj(T[1]),
            }

            if self.T_truncation_order >= 2:
                T_conj_dict[2] = np.conj(T[2])

            # residual dictionary
            residual = {0: None, 1: None, 2: None, 3: None}

            # -------------------------------------
            # special prep for the generated equations
            generated_flag = True

            if not generated_flag:
                e_string = """
                The current code only works using the generated equations.
                It is possible to change the code to work without the generated equations but
                it is not a simple task and requires understanding of how things work.
                You need to contact Neil to generate new equations using termfactory.
                """
                raise Exception(e_string)
            else:
                # map Z into the correct keys for generated code
                gen_Z = {(0, 0): Z[0]}

                if self.Z_truncation_order >= 1:
                    gen_Z[(1, 0)] = Z[1]
                    # gen_Z[(1, 1)] =
                    # gen_Z[(0, 1)] =

                if self.Z_truncation_order >= 2:
                    gen_Z[(2, 0)] = Z[2]

                if self.Z_truncation_order >= 3:
                    gen_Z[(3, 0)] = Z[3]

                # map T into the correct keys for generated code
                _special_T_conj = {
                    (0, 1): T_conj_dict[1],
                    # (0, 2): np.conj(T[2]),
                }

            new_scheme_flag = True
            # -------------------------------------
            # compute intermediate quantities
            # -------------------------------------
            if new_scheme_flag:
                # compute C matrix
                C = self._compute_C_matrix(gen_Z, _special_T_conj, opt_flag=True)

                # compute H_bar_tilde matrix
                H_bar_tilde = self._cal_H_bar_tilde(H_bar, _special_T_conj, opt_flag=True)

            # -------------------------------------
            # compute net residue
            # -------------------------------------

            generated_equation_exception_string = (
                f"Currently optimized paths are only generated for Z2 or Z3 not {self.Z_truncation_order=}\n"
                "You need to either:\n"
                f"set {generated_flag=} to false\n"
                f"set {opt_flag=} to false (note that it is forceably"
                "set to True for models with N+A >= 15 (large models)\n"
                f"contact Neil to generate more optimized paths for different Z values\n"
            )
            if generated_flag and opt_flag:
                if self.Z_truncation_order not in [2, 3]:
                    raise Exception(generated_equation_exception_string)

            # constant
            if not generated_flag:
                residual[0] = self._f_z_0(H_bar_tilde, C, T_conj_dict, opt_flag)
            else:
                if new_scheme_flag:
                    residual[0] = np.zeros(shape=(A, ), dtype=complex)

                    if not opt_flag:
                        z_three_eqns.add_m0_n0_HZ_terms(
                            residual[0], self.ansatz, self.gen_trunc,
                            _special_T_conj, H_bar_tilde, C
                        )
                    else:
                        z_three_eqns.add_m0_n0_HZ_terms_optimized(
                            residual[0], self.ansatz, self.gen_trunc,
                            _special_T_conj, H_bar_tilde, C,
                            self.all_opt_paths[(0, 0)][0]
                        )
                else:
                    if not opt_flag:
                        residual[0] = z_three_eqns.compute_m0_n0_amplitude(
                            A, N, self.ansatz, self.gen_trunc,
                            _special_T_conj, H_bar, gen_Z
                        )
                    else:
                        residual[0] = z_three_eqns.compute_m0_n0_amplitude_optimized(
                            A, N, self.ansatz, self.gen_trunc,
                            _special_T_conj, H_bar, gen_Z,
                            self.all_opt_paths[(0, 0)]
                        )

            # linear
            if self.T_truncation_order >= 1 or self.Z_truncation_order >= 1:
                if not generated_flag:
                    residual[1] = self._f_z_i(H_bar_tilde, C, T_conj_dict, opt_flag)
                else:
                    if new_scheme_flag:
                        residual[1] = np.zeros(shape=(A, N), dtype=complex)

                        if not opt_flag:
                            z_three_eqns.add_m0_n1_HZ_terms(
                                residual[1], self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar_tilde, C
                            )
                        else:
                            z_three_eqns.add_m0_n1_HZ_terms_optimized(
                                residual[1], self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar_tilde, C,
                                self.all_opt_paths[(0, 1)][0]
                            )
                    else:
                        if not opt_flag:
                            residual[1] = z_three_eqns.compute_m0_n1_amplitude(
                                A, N, self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar, gen_Z
                            )
                        else:
                            residual[1] = z_three_eqns.compute_m0_n1_amplitude_optimized(
                                A, N, self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar, gen_Z,
                                self.all_opt_paths[(0, 1)]
                            )

            # quadratic
            if self.T_truncation_order >= 2 or self.Z_truncation_order >= 2:
                if not generated_flag:
                    residual[2] = self._f_z_ij(H_bar_tilde, C, T_conj_dict, opt_flag)
                else:
                    if new_scheme_flag:
                        residual[2] = np.zeros(shape=(A, N, N), dtype=complex)

                        if not opt_flag:
                            z_three_eqns.add_m0_n2_HZ_terms(
                                residual[2], self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar_tilde, C
                            )
                        else:
                            z_three_eqns.add_m0_n2_HZ_terms_optimized(
                                residual[2], self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar_tilde, C,
                                self.all_opt_paths[(0, 2)][0]
                            )
                    else:
                        if not opt_flag:
                            residual[2] = z_three_eqns.compute_m0_n2_amplitude(
                                A, N, self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar, gen_Z
                            )
                        else:
                            residual[2] = z_three_eqns.compute_m0_n2_amplitude_optimized(
                                A, N, self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar, gen_Z,
                                self.all_opt_paths[(0, 2)]
                            )

                # symmetrize
                # residual[2] = symmetrize_tensor(self.N, residual[2], order=2)

            # cubic
            if self.T_truncation_order >= 3 or self.Z_truncation_order >= 3:
                if not generated_flag:
                    assert False, 'this is not implemented'
                    residual[3] = self._f_z_ijk(H_bar_tilde, C, T_conj_dict, opt_flag)
                else:
                    if new_scheme_flag:
                        residual[3] = np.zeros(shape=(A, N, N, N), dtype=complex)

                        if not opt_flag:
                            z_three_eqns.add_m0_n3_HZ_terms(
                                residual[3], self.ansatz, self.gen_trunc, _special_T_conj, H_bar_tilde, C
                            )
                        else:
                            z_three_eqns.add_m0_n3_HZ_terms_optimized(
                                residual[3], self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar_tilde, C, self.all_opt_paths[(0, 3)][0]
                            )
                    else:
                        if not opt_flag:
                            residual[3] = z_three_eqns.compute_m0_n3_amplitude(
                                A, N, self.ansatz, self.gen_trunc, _special_T_conj, H_bar, gen_Z
                            )
                        else:
                            residual[3] = z_three_eqns.compute_m0_n3_amplitude_optimized(
                                A, N, self.ansatz, self.gen_trunc,
                                _special_T_conj, H_bar, gen_Z, self.all_opt_paths[(0, 3)]
                            )

                # symmetrize
                # residual[3] = symmetrize_tensor(self.N, residual[3], order=3)

            # --------------------------------------------------------------------------------
            if new_scheme_flag:
                dT = self._compute_t_residual_new(H_bar_tilde, C)

            else:
                dT = self._compute_t_residual(H_bar_tilde, Z, Z_conj_dict)

            if new_scheme_flag:
                dZ = self._compute_z_residual_new_scheme(residual, Z, _special_T_conj, dT, C)

            else:
                dZ = self._compute_z_residual(residual, Z, T_conj_dict, dT, z_opt_flag=False)

            # apply conjugation?
            for key in dT.keys():
                dT[key] *= -1j

            for key in dZ.keys():
                dZ[key] *= -1j

            return dT, dZ, C

        A, N = self.A, self.N

        # restore the origin shape of t, z amplitudes from y_tensor
        Z_unraveled, T_unraveled = self._unravel_y_tensor(y_tensor)

        if False:  # for debugging
            log.info(f"time = {self.remove_time_step_conversion(time)}")

        if False:  # for debugging
            fs = self.remove_time_step_conversion(time)
            log.info(
                f"On integration step {self.counter:<8d} at {fs:>9.6e} fs\n"
                f"{tab}C(t)(ABS) = {self.C_tau_ABS[-1][1]:>9.4f}\n"
                f"{tab}C(t)(ECD) = {self.C_tau_ECD[-1][1]:>9.4f}"
            )
            log.info('-'*60 + f" {self.counter} " + '-'*60)

        # printing progression
        self._print_integration_progress(time, t_final, Z_unraveled, T_unraveled)

        # debug used for testing very large systems such as cytosine and hexahelicene
        # it takes super long to run these so we have to stop the program at a certain number
        # of integration steps
        if False and (25 < self.counter):
            fs = self.remove_time_step_conversion(time)
            log.info(
                f"On integration step {self.counter:<8d} at {fs:>9.6e} fs\n"
                f"{tab}C(t)(ABS) = {self.C_tau_ABS[-1][1]:>9.4f}\n"
                f"{tab}C(t)(ECD) = {self.C_tau_ECD[-1][1]:>9.4f}"
            )
            log.info('-'*60 + f" {self.counter} " + '-'*60)
            return  # attempt to gracefully stop integration

        # initialize new arrays
        self.norm = np.zeros(A, dtype=complex)

        if self.calculate_population_flag:
            self.state_pop_DB = np.zeros((A, A), dtype=complex)
            self.state_pop_AB = np.zeros((A, A), dtype=complex)

        # will hold the total of the C(t) over all surfaces
        C_tau_ABS = 0.0
        C_tau_ECD = 0.0

        # fill these arrays with zeros
        for key in self.dT.keys():
            self.dT[key].fill(complex(0.0))

        for key in self.dZ.keys():
            self.dZ[key].fill(complex(0.0))

        for b in range(A):

            # pack T amplitude
            t_dict = {0: T_unraveled[0][b]}

            if self.T_truncation_order >= 1:
                t_dict[1] = T_unraveled[1][b, :]

            if self.T_truncation_order >= 2:
                t_dict[2] = T_unraveled[2][b, :]

            # pack Z amplitude
            z_dict = {0: Z_unraveled[0][b]}

            if self.Z_truncation_order >= 1:
                z_dict[1] = Z_unraveled[1][b, :]

            if self.Z_truncation_order >= 2:
                z_dict[2] = Z_unraveled[2][b, :]

            if self.Z_truncation_order >= 3:
                z_dict[3] = Z_unraveled[3][b, :]

            # ## TO DO:
            # ## similarity transform the Hamiltonian over e^t
            H_bar = self._similarity_trans(self.h, t_dict)

            # ## TO DO:
            # ## similarity transform the Hamiltonian over 1+Z and evaluate Z/T residue
            dT, dZ, C = solve_z(H_bar, z_dict, t_dict, b, self.selected_surface)

            # compute ACF for ABS and ECD
            U = z_dict[0] * np.exp(t_dict[0])
            ACF_ABS = self.calculate_ACF_ABS(self.E_tdm, U, b)
            ACF_ECD = self.calculate_ACF_ECD(self.M_tdm, self.E_tdm, U, b)

            # accumulate the ACF's
            C_tau_ABS += ACF_ABS
            C_tau_ECD += ACF_ECD

            new_style_flag = True

            # calculate the norm of this surface
            if new_style_flag:
                norm = self.calculate_CC_norm_new(C)
            else:  # old style
                """we no longer calculate `C_args` when calculating the norm, but instead
                it is calculated above in `solve_z` and the variable `C` replaced `C_args`
                """
                norm, C_args, s_0 = self.calculate_CC_norm(t_dict, z_dict, b)

            # accumulate the norm
            self.norm[b] += norm

            if self.calculate_population_flag:
                if new_style_flag:  # new style
                    density_matrix = self.calculate_density_matrix_new(norm, C)
                    pop_DB, pop_AB = self.calculate_state_populations_new(density_matrix)
                else:
                    raise Exception('s_0 and C_args are not defined, need to confirm with Songhao; broken!?')
                    """
                    For some reason the `calculate_CC_norm` is commented out above, why?
                    This old style of calculating state populations cannot work without s_0 and C_args
                    """
                    pop_DB, pop_AB = self.calculate_state_populations(s_0, norm, C_args)

                self.state_pop_DB[b, :] += pop_DB
                self.state_pop_AB[b, :] += pop_AB

            # compute dZ's
            self.dZ[0][b, :] += dZ[0]
            self.dZ[1][b, :] += dZ[1]

            if self.Z_truncation_order >= 2:
                self.dZ[2][b, :] += dZ[2]

            if self.Z_truncation_order >= 3:
                self.dZ[3][b, :] += dZ[3]

            # compute dT's
            self.dT[0][b] += dT[0]
            self.dT[1][b, :] += dT[1]

            if self.T_truncation_order >= 2:
                self.dT[2][b, :] += dT[2]

        # flatten the z, t tensors into a 1D array
        delta_y_tensor = self._ravel_y_tensor(self.dZ, self.dT)

        # store ACF data
        self.C_tau_ABS.append((time, C_tau_ABS))
        self.C_tau_ECD.append((time, C_tau_ECD))

        # store norm and state pops for plot
        self.Norm.append((time, self.norm))

        if self.calculate_population_flag:
            self.State_pop_DB.append((time, self.state_pop_DB))
            self.State_pop_AB.append((time, self.state_pop_AB))

        return delta_y_tensor

    def rk45_integration(self, t_init=0., t_final=10., density=1.0, nof_points=10000, debug_flag=False):
        """ Runge-Kutta integration scheme
        """

        # ------------------------------------------------------------------------
        # initialize integration parameters
        # ------------------------------------------------------------------------
        log.info(
            "We are going to preform a RK4(5) integration with"
            f"{t_init=} {t_final=} {density=} {nof_points=} {debug_flag=}"
        )

        A, N = self.A, self.N  # to reduce line lengths, for conciseness

        # used for debugging purposes to print out the integration steps every n% of integration
        self.counter = 0
        self.last_counter = 0
        self.last_print = 0

        # this 'raw' step size is not scaled properly
        raw_step_size = (t_final - t_init) / nof_points

        # we need to convert the `raw_step_size`
        step_size = self.apply_time_step_conversion(raw_step_size)
        log.info(
            f"{'Raw step size (unitless):':<26} ({raw_step_size})\n"
            f"{'Real step size (eV/fs):':<26} ({step_size})"
        )

        # initialize accumulation tensor
        U_0 = np.eye(A, dtype=complex)

        # initialize autocorrelation list
        # We need to store the t values for each C(t) so that we can match them up at the end
        self.C_tau_ABS = [(0.0, np.einsum('ac,cd,ad->', self.E_tdm, U_0, self.E_tdm, dtype=complex))]
        self.C_tau_ECD = [(0.0, np.einsum('ac,cd,ad->', self.M_tdm, U_0, self.E_tdm, dtype=complex))]

        # initialize z amplitude
        initial_Z = {
            0: np.eye(A, dtype=complex),
            1: np.zeros([A, A, N], dtype=complex),
            2: np.zeros([A, A, N, N], dtype=complex),
            3: np.zeros([A, A, N, N, N], dtype=complex),
        }

        # initialize t amplitude
        initial_T = {
            0: np.zeros(A, dtype=complex),
            1: np.zeros([A, N], dtype=complex),
            2: np.zeros([A, N, N], dtype=complex),
        }

        # prepare the initial y_tensor
        initial_y_tensor = self._ravel_y_tensor(initial_Z, initial_T)

        # ------------------------------------------------------------------------
        # the integration function called by the `solve_ivp` integrator at each step of integration
        # ------------------------------------------------------------------------

        # specify the precision of the integrator so that the output for the test models is numerically identical
        if self.comparing_to_test_models:
            relative_tolerance = 1e-10
            absolute_tolerance = 1e-12
        else:
            # So the new "fix" does work but we need to up the tolerance parameters of the RK integrator to get agreement with MCTDH.
            # I would do at least 1e-9 and 1e-12, the issue is that the integrator simply takes too large steps based on the ode's represented by S_0
            relative_tolerance = 1e-07
            absolute_tolerance = 1e-09
        # ------------------------------------------------------------------------
        # call the integrator
        # ------------------------------------------------------------------------

        integration_function = self.rk45_solve_ivp_integration_function

        sol = new_solve_ivp(
            fun=integration_function,  # the function we are integrating
            method="RK45",  # the integration method we are using
            first_step=step_size,  # fix the initial step size
            t_span=(
                t_init,  # initial time
                self.apply_time_step_conversion(t_final),  # boundary time, integration end point
            ),
            y0=initial_y_tensor,  # initial state - shape (n, )
            args=(t_final, ),  # extra args to pass to `rk45_solve_ivp_integration_function`
            # max_step=self.step_size,  # maximum allowed step size
            rtol=relative_tolerance,  # relative tolerance
            atol=absolute_tolerance,  # absolute tolerance
            store_y_values=False,  # do not store the y values over the integration
            t_eval=None,  # store all the time values that we integrated over
            dense_output=False,  # extra debug information
            # we do not need to vectorize
            # this means to process multiple time steps inside the function `rk45_solve_ivp_integration_function`
            # it would be useful for a method which does some kind of block stepping
            vectorized=False,
        )

        # ------------------------------------------------------------------------
        # now we extract the relevant information from the integrator object `sol`
        # ------------------------------------------------------------------------
        self._postprocess_rk45_integration_results(sol, debug=debug_flag)

        return

    @staticmethod
    def _save_data(file_path, x_array, y_array):
        """ Save complex function y(x) into plain text file formatted as per MCTDH auto files

        This preforms file saving for both `save_acf_data` and `save_sos_data`
        `x_array` and `y_array` should be numpy arrays of the same length.
        `x_array` should have dtype float
        `y_array` should have dtype complex
        """

        # some simple data checking
        assert isinstance(x_array, np.ndarray), \
            f"x_array is not an numpy array it is a {type(x_array)}"
        assert isinstance(y_array, np.ndarray), \
            f"y_array is not an numpy array it is a {type(x_array)}"

        assert x_array.dtype == float, \
            f"x_array has type {x_array.dtype} instead of float"
        assert y_array.dtype == complex, \
            f"y_array has type {y_array.dtype} instead of complex"

        assert len(x_array) == len(y_array.real) == len(y_array.imag), \
            (
                f"Arrays must be the same length\n"
                f"x_array      has length: {len(x_array)}\n"
                f"y_array.real has length: {len(y_array.real)}\n"
                f"y_array.imag has length: {len(y_array.imag)}"
        )

        header = (
            f"#{'time[fs]': >12}"
            f" {'Re(autocorrel)': >23}"
            f" {'Im(autocorrel)': >19}"
            f" {'Abs(autocorrel)': >20}"
        )

        file_data = header + '\n'

        y_absolute = np.abs(y_array)

        # we assume all arrays are the same length
        for idx in range(len(x_array)):
            file_data += (
                f"{x_array[idx]: > 15.8f}"
                f" {y_array.real[idx]: > 23.14f}"
                f" {y_array.imag[idx]: > 19.14f}"
                f" {y_absolute[idx]: > 19.14f}"
                "\n"
            )

        # write data to file
        with open(file_path, 'w') as file:
            file.write(file_data)
        return

    @classmethod
    def class_save_acf_data(cls, file_name, time, acf, output_path="./"):
        """ Store provided auto correlation function (ACF) into plain text file """
        file_path = join(output_path, f"ACF_CC_{file_name}.txt")
        cls._save_data(file_path, time, acf)
        return file_path

    @classmethod
    def class_save_sos_data(cls, file_name, time, sos, output_path="./"):
        """ Store provided Sum-Over-States data (SOS) into plain text file """
        file_path = join(output_path, f"ACF_SOS_{file_name}.txt")
        cls._save_data(file_path, time, sos)
        return file_path

    def save_acf_data(self, file_name, time=None, acf_ABS=None, acf_ECD=None, output_path="./"):
        """ Store auto correlation function (ACF) into plain text file

        The file is saved to `output_path` + "ACF_CC_{file_name}.txt"
        Caller can provide arguments `time` and/or `acf` as alternatives to `self.t_cc` and `self.C_tau_cc`
        if they want to preform some data processing before saving to a file.
        """
        if time is None:
            time = self.t_cc
        if acf_ABS is None:
            acf_ABS = self.C_tau_cc_ABS
        if acf_ECD is None:
            acf_ECD = self.C_tau_cc_ECD

        file_path_ABS = join(output_path, f"ACF_ABS_CC_{file_name}.txt")
        self._save_data(file_path_ABS, time, acf_ABS)

        file_path_ECD = join(output_path, f"ACF_ECD_CC_{file_name}.txt")
        self._save_data(file_path_ECD, time, acf_ECD)
        return file_path_ABS, file_path_ECD

    def save_sos_data(self, file_name, time=None, sos=None, output_path="./"):
        """ Store SOS (for debugging/testing/profiling) into plain text file

        Caller can provide arguments `time` and/or `sos` as alternatives to `self.t_sos` and `self.C_tau_sos`
        if they want to preform some data processing before saving to a file
        """
        if time is None:
            time = self.t_sos
        if sos is None:
            sos = self.C_tau_sos

        file_path = join(output_path, f"ACF_SOS_{file_name}.txt")
        self._save_data(file_path, time, sos)
        return file_path

    def save_raw_integration_data(self, file_name, output_path="./"):
        """ store this objects ACF and SOS values to text files without any preprocessing """

        # but now we use the same file_name for both ACf and sos?
        acf_path = self.save_acf_data(file_name, output_path=output_path)
        sos_path = self.save_sos_data(file_name, output_path=output_path)
        return acf_path, sos_path

    @staticmethod
    def _load_auto_data(file_path):
        """ Returns a tuple (x, y),
            x is array of time steps
            y is array of complex values at each x
        Only works on `auto` files generated by MCTDH.
        """

        # read data from file
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
        except Exception as e:
            raise e

        # starting from the end of the file go line by line to the front of the file
        for line in lines[-1:0:-1]:
            length = len(line)
            # all lines of data are at least 77 characters long
            if length < 10:
                del lines[-1]
            elif length >= 77:
                break
            else:
                raise Exception(str(
                    "Line length is > 10 but < 77?\n"
                    "This file has something strange going on at the end. Please check.\n"
                    f"The line is {length} characters long.\n"
                    f"This is the offending line:\n{[line,]}\n"
                ))

        # the first line is a comment
        num_data_points = len(lines) - 1
        y_val = np.zeros(num_data_points, dtype=complex)
        time = np.zeros(num_data_points, dtype=float)

        p = parse.compile("{t:>15.8f} {y_real:>22.14f} {y_imag:>18.14f} {y_abs:>18.14f}")

        # extract data from lines of file
        for idx, line in enumerate(lines[1:]):

            r = p.parse(line)
            time[idx] = r['t']
            y_val.real[idx] = r['y_real']
            y_val.imag[idx] = r['y_imag']

        return time, y_val

    @staticmethod
    def load_auto_data(file_path):
        """ Wrapper function for `_load_auto_data`
        Returns a tuple (x, y),  x is array of time steps, y is array of ACF value at each x.
        Only works on `auto` files generated by MCTDH.
        """
        time, acf = vibronic_hamiltonian._load_auto_data(file_path)
        return time, acf

    @staticmethod
    def _load_data(file_path):
        """ Returns a tuple (x, y),
            x is array of time steps
            y is array of complex values at each x
        """

        # read data from file
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
        except Exception as e:
            raise e

        # starting from the end of the file go line by line to the front of the file
        for line in lines[-1:0:-1]:
            length = len(line)
            # all lines of data are at least 80 characters long
            if length < 10:
                del lines[-1]
            elif length >= 80:
                break
            else:
                raise Exception(str(
                    "Line length is > 10 but < 80?\n"
                    "This file has something strange going on at the end. Please check.\n"
                    f"The line is {length} characters long.\n"
                    f"This is the offending line:\n{[line,]}\n"
                ))

        # the first line is a comment
        num_data_points = len(lines) - 1
        y_val = np.zeros(num_data_points, dtype=complex)
        time = np.zeros(num_data_points, dtype=float)

        p = parse.compile("{t:>15.8f} {y_real:>23.14f} {y_imag:>19.14f} {y_abs:>19.14f}")

        # check to make sure parse isn't broken, as per issue #88
        r_test = p.parse(lines[1])
        assert r_test is not None, (
            "the python package parse was unable to process the first valid line in the dataset\n"
            "if your data set is properly formatted this likely indicates your using the wrong version"
            "of parse. version 1.18 should be stable, but not 1.19"
        )

        # extract data from lines of file
        for idx, line in enumerate(lines[1:]):

            r = p.parse(line)
            time[idx] = r['t']
            y_val.real[idx] = r['y_real']
            y_val.imag[idx] = r['y_imag']

        return time, y_val

    @staticmethod
    def load_acf_data(file_path):
        """ Wrapper function for `_load_data`
        Returns a tuple (x, y),  x is array of time steps, y is array of ACF value at each x
        """
        time, acf = vibronic_hamiltonian._load_data(file_path)
        return time, acf

    @staticmethod
    def load_sos_data(file_path):
        """ Wrapper function for `_load_data`
        Returns a tuple (x, y),  x is array of time steps, y is array of SOS value at each x
        """
        time, sos = vibronic_hamiltonian._load_data(file_path)
        return time, sos

    def load_inplace_acf_data(self, file_path):
        """ Wrapper function for `_load_data`
        Overwrite `self.C_tau_cc` and `self.t_cc` with the values of the file at `file_path`.
        """
        time, acf = self._load_data(file_path)
        self.C_tau_cc = acf.copy()
        self.t_cc = time.copy()
        return

    def load_inplace_sos_data(self, file_path):
        """ Wrapper function for `_load_data`
        Overwrite `self.C_tau_sos` and `self.t_sos` with the values of the file at `file_path`.
        """
        time, sos = self._load_data(file_path)
        self.C_tau_sos = sos.copy()
        self.t_sos = time.copy()
        return

    def plot_acf(
        self, file_name, sos_flag=False, output_path="./",
        plot_log=False, long_propagation_zoom=False,
        acf_y_array=None, acf_x_array=None, sos_y_array=None, sos_x_array=None,
    ):
        r""" Plot the Auto Correlation Function (ACF)

        `plot_log` is Boolean, if True an additional plot with log scale y-axis is generated.
        `long_propagation_zoom` if Boolean, if True changes x-axis to show only last n-1 femtoseconds of ACF in plot.

        The arguments `acf_y_array`, `acf_x_array`, `sos_y_array`, `sos_x_array` are provided so that the caller can provide alternative data for plotting.
        For large systems it would be advantageous to run one calculation of SOS data for a large number of basis functions and save the data.
        Then the ACF integration code could be called multiple times, and compared against the SOS data
        simply by loading the data from the file and providing it in arguments `sos_y_array`, `sos_x_array`.

        The plot is saved to a file at the path `output_path` + "ACF_{file_name}.png"
        an example: \home\username\plots\ACF_h2o_linear.png
        IF the `plot_log` boolean flag is true an additional plot with logarithmic scales is saved to `output_path` + "ACF_{file_name}_logscale.png"
        an example: \home\username\plots\ACF_h2o_linear_logscale.png
        """

        # if no alternative data is provided use the internal data for plotting
        if acf_y_array is None:
            acf_y_array = self.C_tau_cc_ABS
        if acf_x_array is None:
            acf_x_array = self.t_cc
        if sos_y_array is None and sos_flag:
            sos_y_array = self.C_tau_sos
        if sos_x_array is None and sos_flag:
            sos_x_array = self.t_sos

        # set style for all text in plot
        mpl.rc('font', **{'family': 'sans-serif', 'weight': 'normal', 'size': 20})

        fig = plt.figure(figsize=(14, 10))
        fig.suptitle(f"ACF of {file_name}", fontweight='bold')

        # if plotting SOS we need two different subplots
        if sos_flag:
            ax_real = fig.add_subplot(211)
            ax_imag = fig.add_subplot(212)
        else:
            ax_real = fig.add_subplot(111)
            ax_imag = ax_real  # just a pointer to `ax_real`

        if not sos_flag:
            ax_real.plot(
                acf_x_array, acf_y_array.real, label='Real component of CC',
                linewidth=3.0, color='black', alpha=0.7,
            )
            ax_imag.plot(
                acf_x_array, acf_y_array.imag, label='Imaginary component of CC',
                linewidth=3.0, color='black', alpha=0.7, ls='-.',
            )
        else:

            # for debugging we might also plot the SOS results
            ax_real.plot(
                sos_x_array, sos_y_array.real, label='Real component of SOS',
                linewidth=3.0, color='red', alpha=0.6,
            )
            ax_imag.plot(
                sos_x_array, sos_y_array.imag, label='Imaginary component of SOS',
                linewidth=3.0, color='red', alpha=0.6,
            )

            # tf = acf_x_array[-1]
            # npoints = int(tf*1e2)
            # from scipy.interpolate import interp1d
            # x, y = acf_x_array, acf_y_array.real
            # f1 = interp1d(x, y)
            # f2 = interp1d(x, y, kind='quadratic')
            # f3 = interp1d(x, y, kind='cubic')

            # xnew = np.linspace(0.0, tf, num=npoints, endpoint=True)
            # ax_real.plot(x, y, 'o', color='blue', linewidth=3.0, label='data')
            # ax_real.plot(xnew, f1(xnew), '-', color='orange', linewidth=3.0, label='linear')
            # ax_real.plot(xnew, f2(xnew), '--', color='green', linewidth=3.0, label='quadratic')
            # ax_real.plot(xnew, f3(xnew), '-.', color='magenta', linewidth=3.0, label='cubic')

            # x, y = acf_x_array, acf_y_array.imag
            # f1 = interp1d(x, y)
            # f2 = interp1d(x, y, kind='quadratic')
            # f3 = interp1d(x, y, kind='cubic')

            # ax_imag.plot(x, y, 'o', color='blue', linewidth=3.0, label='data')
            # ax_imag.plot(xnew, f1(xnew), '-', color='orange', linewidth=3.0, label='linear')
            # ax_imag.plot(xnew, f2(xnew), '--', color='green', linewidth=3.0, label='quadratic')
            # ax_imag.plot(xnew, f3(xnew), '-.', color='magenta', linewidth=3.0, label='cubic')

            # plot the CC after so its plotted 'on top of' the red SOS for better visibility
            ax_real.plot(
                acf_x_array, acf_y_array.real, label='Real component of CC',
                linewidth=3.0, color='black', alpha=1.0, ls='-.',
            )
            ax_imag.plot(
                acf_x_array, acf_y_array.imag, label='Imaginary component of CC',
                linewidth=3.0, color='black', alpha=1.0, ls='-.',
            )

        # set x labels
        ax_real.set_xlabel(r'''$\tau \,\slash\, \hbar$''')
        ax_imag.set_xlabel(r'''$\tau \,\slash\, \hbar$''')
        # set y labels
        ax_real.set_ylabel(r'''C($\tau \,\slash\, \hbar$)''')
        ax_imag.set_ylabel(r'''C($\tau \,\slash\, \hbar$)''')
        # create legend
        ax_real.legend()
        ax_imag.legend()

        # if using long propagation(> 100fs) we can zoom in on the end
        if long_propagation_zoom:
            x1, x2 = ax_real.get_xlim()
            x1 += abs(x1)
            x2 -= abs(x1)
            ax_real.set_xlim([x2 - 1, x2])
            ax_imag.set_xlim([x2 - 1, x2])

        # save plot
        plot_path = join(output_path, f"ACF_{file_name}.png")
        fig.savefig(plot_path)

        # save log plot
        if plot_log:
            ax_real.set_yscale('log')
            ax_imag.set_yscale('log')
            plot_path = join(output_path, f"ACF_{file_name}_logscale.png")
            fig.savefig(plot_path)

        plt.close(fig)
        return

    def plot_acf_and_save_data(self, file_name, sos_flag=False, output_path="./", **kwargs):
        """ Wrapper function for `plot_acf`, `save_acf_data`, and `save_sos_data`."""
        self.plot_acf(file_name, sos_flag, output_path, **kwargs)
        self.save_acf_data(file_name, output_path=output_path, **kwargs)
        if sos_flag:
            self.save_sos_data(file_name, output_path=output_path, **kwargs)
        return

    @classmethod
    def class_smooth_acf_using_interpolation(cls, time, acf, t_init=0., t_final=10., nof_points=10000, kind='cubic'):
        """ use interpolation to get evenly spaced points """
        f_real = interp1d(time, acf.real, kind)
        f_imag = interp1d(time, acf.imag, kind)

        # first we generate a reasonable dt step
        # if we use endpoint=True then its very difficult to get good deltas
        _, dt = np.linspace(t_init, t_final, num=nof_points, endpoint=False, retstep=True)
        time_step = round(dt, ndigits=2)
        print(f"Time step is {time_step:12.8f}")

        # then we make sure we generate uniformly spaced points
        # we use (t_final+time_step) because np.arange doesn't include the last point
        # so if we want to include t_final (we do) then we need to make the 'last point'
        # one step after the point that we want to be the real 'last point'
        new_x = np.arange(t_init, t_final+time_step, step=time_step, dtype=float)

        new_acf = np.zeros_like(new_x, dtype=complex)
        new_acf.real = f_real(new_x)
        new_acf.imag = f_imag(new_x)
        return new_x, new_acf

    def smooth_acf_using_interpolation(self, nof_points=10000, kind='cubic'):
        """ use interpolation to get evenly spaced points """
        x = self.t_cc
        y_real = self.C_tau_cc.real
        y_imag = self.C_tau_cc.imag

        f_real = interp1d(x, y_real, kind)
        f_imag = interp1d(x, y_imag, kind)

        t_init, t_final = x[0], x[-1]

        # first we generate a reasonable dt step
        # if we use endpoint=True then its very difficult to get good deltas
        _, dt = np.linspace(t_init, t_final, num=nof_points, endpoint=False, retstep=True)
        time_step = round(dt, ndigits=2)
        print(f"Time step is {time_step:12.8f}")

        # then we make sure we generate uniformly spaced points
        # we use (t_final+time_step) because np.arange doesn't include the last point
        # so if we want to include t_final (we do) then we need to make the 'last point'
        # one step after the point that we want to be the real 'last point'
        new_x = np.arange(t_init, t_final+time_step, step=time_step, dtype=float)

        self.t_cc = new_x.copy()
        self.C_tau_cc = np.zeros_like(new_x, dtype=complex)
        self.C_tau_cc.real = f_real(new_x)
        self.C_tau_cc.imag = f_imag(new_x)
        return
