# system imports
from math import factorial

# third party imports
import numpy as np
import opt_einsum as oe

# local imports
from ..log_conf import log
from .symmetrize import symmetrize_tensor
from . import residual_equations

# --------------------------------------------------------------------------- #
# ---------------------------- DEFAULT FUNCTIONS ---------------------------- #
# --------------------------------------------------------------------------- #

# ---------------------------- DISCONNECTED TERMS ---------------------------- #

def _order_1_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1)
    """
    raise Exception(
        "the first possible purely VECI/CC term is (1, 1) or (dt_i * t_i)"
        "which requires a residual of at least 2nd order"
    )

def _order_2_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1)
    """
    # unpack the `t_args` and 'dt_args'
    t_i, *unusedargs = t_args
    dt_i, *unusedargs = dt_args
    # Creating the 2nd order return array
    linked_disconnected_terms = np.zeros((A, A, N, N), dtype=complex)
    # the (1, 1) term
    linked_disconnected_terms += 1/factorial(2) * (
        np.einsum('aci, cbj->abij', dt_i, t_i) +
        np.einsum('aci, cbj->abij', t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_3_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1)
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, *unusedargs = t_args
    dt_i, dt_ij, *unusedargs = dt_args
    # Creating the 3rd order return array
    linked_disconnected_terms = np.zeros((A, A, N, N, N), dtype=complex)
    # the (2, 1) term
    linked_disconnected_terms += 1/(factorial(2) * factorial(2)) * (
        np.einsum('aci, cbjk->abijk', dt_i, t_ij) +
        np.einsum('aci, cbjk->abijk', t_i, dt_ij) +
        np.einsum('acij, cbk->abijk', dt_ij, t_i) +
        np.einsum('acij, cbk->abijk', t_ij, dt_i)
    )
    # the (1, 1, 1) term
    linked_disconnected_terms += 1/factorial(3) * (
        np.einsum('aci, cdj, dbk->abijk', dt_i, t_i, t_i) +
        np.einsum('aci, cdj, dbk->abijk', t_i, dt_i, t_i) +
        np.einsum('aci, cdj, dbk->abijk', t_i, t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_4_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1)
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, *unusedargs = dt_args
    # Creating the 4th order return array
    linked_disconnected_terms = np.zeros((A, A, N, N, N, N), dtype=complex)
    # the (3, 1) term
    linked_disconnected_terms += 1/(factorial(2) * factorial(3)) * (
        np.einsum('aci, cbjkl->abijkl', dt_i, t_ijk) +
        np.einsum('aci, cbjkl->abijkl', t_i, dt_ijk) +
        np.einsum('acijk, cbl->abijkl', dt_ijk, t_i) +
        np.einsum('acijk, cbl->abijkl', t_ijk, dt_i)
    )
    # the (2, 1, 1) term
    linked_disconnected_terms += 1/(factorial(3) * factorial(2)) * (
        np.einsum('aci, cdj, dbkl->abijkl', dt_i, t_i, t_ij) +
        np.einsum('aci, cdj, dbkl->abijkl', t_i, dt_i, t_ij) +
        np.einsum('aci, cdj, dbkl->abijkl', t_i, t_i, dt_ij) +
        np.einsum('aci, cdjk, dbl->abijkl', dt_i, t_ij, t_i) +
        np.einsum('aci, cdjk, dbl->abijkl', t_i, dt_ij, t_i) +
        np.einsum('aci, cdjk, dbl->abijkl', t_i, t_ij, dt_i) +
        np.einsum('acij, cdk, dbl->abijkl', dt_ij, t_i, t_i) +
        np.einsum('acij, cdk, dbl->abijkl', t_ij, dt_i, t_i) +
        np.einsum('acij, cdk, dbl->abijkl', t_ij, t_i, dt_i)
    )
    # the (1, 1, 1, 1) term
    linked_disconnected_terms += 1/factorial(4) * (
        np.einsum('aci, cdj, dek, ebl->abijkl', dt_i, t_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, ebl->abijkl', t_i, dt_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, ebl->abijkl', t_i, t_i, dt_i, t_i) +
        np.einsum('aci, cdj, dek, ebl->abijkl', t_i, t_i, t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_5_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1)
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, dt_ijkl, *unusedargs = dt_args
    # Creating the 5th order return array
    linked_disconnected_terms = np.zeros((A, A, N, N, N, N, N), dtype=complex)
    # the (4, 1) term
    linked_disconnected_terms += 1/(factorial(2) * factorial(4)) * (
        np.einsum('aci, cbjklm->abijklm', dt_i, t_ijkl) +
        np.einsum('aci, cbjklm->abijklm', t_i, dt_ijkl) +
        np.einsum('acijkl, cbm->abijklm', dt_ijkl, t_i) +
        np.einsum('acijkl, cbm->abijklm', t_ijkl, dt_i)
    )
    # the (3, 1, 1) term
    linked_disconnected_terms += 1/(factorial(3) * factorial(3)) * (
        np.einsum('aci, cdj, dbklm->abijklm', dt_i, t_i, t_ijk) +
        np.einsum('aci, cdj, dbklm->abijklm', t_i, dt_i, t_ijk) +
        np.einsum('aci, cdj, dbklm->abijklm', t_i, t_i, dt_ijk) +
        np.einsum('aci, cdjkl, dbm->abijklm', dt_i, t_ijk, t_i) +
        np.einsum('aci, cdjkl, dbm->abijklm', t_i, dt_ijk, t_i) +
        np.einsum('aci, cdjkl, dbm->abijklm', t_i, t_ijk, dt_i) +
        np.einsum('acijk, cdl, dbm->abijklm', dt_ijk, t_i, t_i) +
        np.einsum('acijk, cdl, dbm->abijklm', t_ijk, dt_i, t_i) +
        np.einsum('acijk, cdl, dbm->abijklm', t_ijk, t_i, dt_i)
    )
    # the (2, 1, 1, 1) term
    linked_disconnected_terms += 1/(factorial(4) * factorial(2)) * (
        np.einsum('aci, cdj, dek, eblm->abijklm', dt_i, t_i, t_i, t_ij) +
        np.einsum('aci, cdj, dek, eblm->abijklm', t_i, dt_i, t_i, t_ij) +
        np.einsum('aci, cdj, dek, eblm->abijklm', t_i, t_i, dt_i, t_ij) +
        np.einsum('aci, cdj, dek, eblm->abijklm', t_i, t_i, t_i, dt_ij) +
        np.einsum('aci, cdj, dekl, ebm->abijklm', dt_i, t_i, t_ij, t_i) +
        np.einsum('aci, cdj, dekl, ebm->abijklm', t_i, dt_i, t_ij, t_i) +
        np.einsum('aci, cdj, dekl, ebm->abijklm', t_i, t_i, dt_ij, t_i) +
        np.einsum('aci, cdj, dekl, ebm->abijklm', t_i, t_i, t_ij, dt_i) +
        np.einsum('aci, cdjk, del, ebm->abijklm', dt_i, t_ij, t_i, t_i) +
        np.einsum('aci, cdjk, del, ebm->abijklm', t_i, dt_ij, t_i, t_i) +
        np.einsum('aci, cdjk, del, ebm->abijklm', t_i, t_ij, dt_i, t_i) +
        np.einsum('aci, cdjk, del, ebm->abijklm', t_i, t_ij, t_i, dt_i) +
        np.einsum('acij, cdk, del, ebm->abijklm', dt_ij, t_i, t_i, t_i) +
        np.einsum('acij, cdk, del, ebm->abijklm', t_ij, dt_i, t_i, t_i) +
        np.einsum('acij, cdk, del, ebm->abijklm', t_ij, t_i, dt_i, t_i) +
        np.einsum('acij, cdk, del, ebm->abijklm', t_ij, t_i, t_i, dt_i)
    )
    # the (1, 1, 1, 1, 1) term
    linked_disconnected_terms += 1/factorial(5) * (
        np.einsum('aci, cdj, dek, efl, fbm->abijklm', dt_i, t_i, t_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fbm->abijklm', t_i, dt_i, t_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fbm->abijklm', t_i, t_i, dt_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fbm->abijklm', t_i, t_i, t_i, dt_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fbm->abijklm', t_i, t_i, t_i, t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_6_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1)
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, dt_ijkl, dt_ijklm, *unusedargs = dt_args
    # Creating the 6th order return array
    linked_disconnected_terms = np.zeros((A, A, N, N, N, N, N, N), dtype=complex)
    # the (5, 1) term
    linked_disconnected_terms += 1/(factorial(2) * factorial(5)) * (
        np.einsum('aci, cbjklmn->abijklmn', dt_i, t_ijklm) +
        np.einsum('aci, cbjklmn->abijklmn', t_i, dt_ijklm) +
        np.einsum('acijklm, cbn->abijklmn', dt_ijklm, t_i) +
        np.einsum('acijklm, cbn->abijklmn', t_ijklm, dt_i)
    )
    # the (4, 1, 1) term
    linked_disconnected_terms += 1/(factorial(3) * factorial(4)) * (
        np.einsum('aci, cdj, dbklmn->abijklmn', dt_i, t_i, t_ijkl) +
        np.einsum('aci, cdj, dbklmn->abijklmn', t_i, dt_i, t_ijkl) +
        np.einsum('aci, cdj, dbklmn->abijklmn', t_i, t_i, dt_ijkl) +
        np.einsum('aci, cdjklm, dbn->abijklmn', dt_i, t_ijkl, t_i) +
        np.einsum('aci, cdjklm, dbn->abijklmn', t_i, dt_ijkl, t_i) +
        np.einsum('aci, cdjklm, dbn->abijklmn', t_i, t_ijkl, dt_i) +
        np.einsum('acijkl, cdm, dbn->abijklmn', dt_ijkl, t_i, t_i) +
        np.einsum('acijkl, cdm, dbn->abijklmn', t_ijkl, dt_i, t_i) +
        np.einsum('acijkl, cdm, dbn->abijklmn', t_ijkl, t_i, dt_i)
    )
    # the (3, 1, 1, 1) term
    linked_disconnected_terms += 1/(factorial(4) * factorial(3)) * (
        np.einsum('aci, cdj, dek, eblmn->abijklmn', dt_i, t_i, t_i, t_ijk) +
        np.einsum('aci, cdj, dek, eblmn->abijklmn', t_i, dt_i, t_i, t_ijk) +
        np.einsum('aci, cdj, dek, eblmn->abijklmn', t_i, t_i, dt_i, t_ijk) +
        np.einsum('aci, cdj, dek, eblmn->abijklmn', t_i, t_i, t_i, dt_ijk) +
        np.einsum('aci, cdj, deklm, ebn->abijklmn', dt_i, t_i, t_ijk, t_i) +
        np.einsum('aci, cdj, deklm, ebn->abijklmn', t_i, dt_i, t_ijk, t_i) +
        np.einsum('aci, cdj, deklm, ebn->abijklmn', t_i, t_i, dt_ijk, t_i) +
        np.einsum('aci, cdj, deklm, ebn->abijklmn', t_i, t_i, t_ijk, dt_i) +
        np.einsum('aci, cdjkl, dem, ebn->abijklmn', dt_i, t_ijk, t_i, t_i) +
        np.einsum('aci, cdjkl, dem, ebn->abijklmn', t_i, dt_ijk, t_i, t_i) +
        np.einsum('aci, cdjkl, dem, ebn->abijklmn', t_i, t_ijk, dt_i, t_i) +
        np.einsum('aci, cdjkl, dem, ebn->abijklmn', t_i, t_ijk, t_i, dt_i) +
        np.einsum('acijk, cdl, dem, ebn->abijklmn', dt_ijk, t_i, t_i, t_i) +
        np.einsum('acijk, cdl, dem, ebn->abijklmn', t_ijk, dt_i, t_i, t_i) +
        np.einsum('acijk, cdl, dem, ebn->abijklmn', t_ijk, t_i, dt_i, t_i) +
        np.einsum('acijk, cdl, dem, ebn->abijklmn', t_ijk, t_i, t_i, dt_i)
    )
    # the (2, 1, 1, 1, 1) term
    linked_disconnected_terms += 1/(factorial(5) * factorial(2)) * (
        np.einsum('aci, cdj, dek, efl, fbmn->abijklmn', dt_i, t_i, t_i, t_i, t_ij) +
        np.einsum('aci, cdj, dek, efl, fbmn->abijklmn', t_i, dt_i, t_i, t_i, t_ij) +
        np.einsum('aci, cdj, dek, efl, fbmn->abijklmn', t_i, t_i, dt_i, t_i, t_ij) +
        np.einsum('aci, cdj, dek, efl, fbmn->abijklmn', t_i, t_i, t_i, dt_i, t_ij) +
        np.einsum('aci, cdj, dek, efl, fbmn->abijklmn', t_i, t_i, t_i, t_i, dt_ij) +
        np.einsum('aci, cdj, dek, eflm, fbn->abijklmn', dt_i, t_i, t_i, t_ij, t_i) +
        np.einsum('aci, cdj, dek, eflm, fbn->abijklmn', t_i, dt_i, t_i, t_ij, t_i) +
        np.einsum('aci, cdj, dek, eflm, fbn->abijklmn', t_i, t_i, dt_i, t_ij, t_i) +
        np.einsum('aci, cdj, dek, eflm, fbn->abijklmn', t_i, t_i, t_i, dt_ij, t_i) +
        np.einsum('aci, cdj, dek, eflm, fbn->abijklmn', t_i, t_i, t_i, t_ij, dt_i) +
        np.einsum('aci, cdj, dekl, efm, fbn->abijklmn', dt_i, t_i, t_ij, t_i, t_i) +
        np.einsum('aci, cdj, dekl, efm, fbn->abijklmn', t_i, dt_i, t_ij, t_i, t_i) +
        np.einsum('aci, cdj, dekl, efm, fbn->abijklmn', t_i, t_i, dt_ij, t_i, t_i) +
        np.einsum('aci, cdj, dekl, efm, fbn->abijklmn', t_i, t_i, t_ij, dt_i, t_i) +
        np.einsum('aci, cdj, dekl, efm, fbn->abijklmn', t_i, t_i, t_ij, t_i, dt_i) +
        np.einsum('aci, cdjk, del, efm, fbn->abijklmn', dt_i, t_ij, t_i, t_i, t_i) +
        np.einsum('aci, cdjk, del, efm, fbn->abijklmn', t_i, dt_ij, t_i, t_i, t_i) +
        np.einsum('aci, cdjk, del, efm, fbn->abijklmn', t_i, t_ij, dt_i, t_i, t_i) +
        np.einsum('aci, cdjk, del, efm, fbn->abijklmn', t_i, t_ij, t_i, dt_i, t_i) +
        np.einsum('aci, cdjk, del, efm, fbn->abijklmn', t_i, t_ij, t_i, t_i, dt_i) +
        np.einsum('acij, cdk, del, efm, fbn->abijklmn', dt_ij, t_i, t_i, t_i, t_i) +
        np.einsum('acij, cdk, del, efm, fbn->abijklmn', t_ij, dt_i, t_i, t_i, t_i) +
        np.einsum('acij, cdk, del, efm, fbn->abijklmn', t_ij, t_i, dt_i, t_i, t_i) +
        np.einsum('acij, cdk, del, efm, fbn->abijklmn', t_ij, t_i, t_i, dt_i, t_i) +
        np.einsum('acij, cdk, del, efm, fbn->abijklmn', t_ij, t_i, t_i, t_i, dt_i)
    )
    # the (1, 1, 1, 1, 1, 1) term
    linked_disconnected_terms += 1/factorial(6) * (
        np.einsum('aci, cdj, dek, efl, fgm, gbn->abijklmn', dt_i, t_i, t_i, t_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fgm, gbn->abijklmn', t_i, dt_i, t_i, t_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fgm, gbn->abijklmn', t_i, t_i, dt_i, t_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fgm, gbn->abijklmn', t_i, t_i, t_i, dt_i, t_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fgm, gbn->abijklmn', t_i, t_i, t_i, t_i, dt_i, t_i) +
        np.einsum('aci, cdj, dek, efl, fgm, gbn->abijklmn', t_i, t_i, t_i, t_i, t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_1_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    """
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a residual of at least 4th order"
    )

def _order_2_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    """
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a residual of at least 4th order"
    )

def _order_3_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    """
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a residual of at least 4th order"
    )

def _order_4_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, *unusedargs = dt_args
    # Creating the 4th order return array
    un_linked_disconnected_terms = np.zeros((A, A, N, N, N, N), dtype=complex)
        # the (2, 2) term
    un_linked_disconnected_terms += 1/(factorial(2) * factorial(2) * factorial(2)) * (
        np.einsum('acij, cbkl->abijkl', dt_ij, t_ij) +
        np.einsum('acij, cbkl->abijkl', t_ij, dt_ij)
    )

    return un_linked_disconnected_terms

def _order_5_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, dt_ijkl, *unusedargs = dt_args
    # Creating the 5th order return array
    un_linked_disconnected_terms = np.zeros((A, A, N, N, N, N, N), dtype=complex)
        # the (3, 2) term
    un_linked_disconnected_terms += 1/(factorial(2) * factorial(3) * factorial(2)) * (
        np.einsum('acij, cbklm->abijklm', dt_ij, t_ijk) +
        np.einsum('acij, cbklm->abijklm', t_ij, dt_ijk) +
        np.einsum('acijk, cblm->abijklm', dt_ijk, t_ij) +
        np.einsum('acijk, cblm->abijklm', t_ijk, dt_ij)
    )
        # the (2, 2, 1) term
    un_linked_disconnected_terms += 1/(factorial(3) * factorial(2) * factorial(2)) * (
        np.einsum('aci, cdjk, dblm->abijklm', dt_i, t_ij, t_ij) +
        np.einsum('aci, cdjk, dblm->abijklm', t_i, dt_ij, t_ij) +
        np.einsum('aci, cdjk, dblm->abijklm', t_i, t_ij, dt_ij) +
        np.einsum('acij, cdk, dblm->abijklm', dt_ij, t_i, t_ij) +
        np.einsum('acij, cdk, dblm->abijklm', t_ij, dt_i, t_ij) +
        np.einsum('acij, cdk, dblm->abijklm', t_ij, t_i, dt_ij) +
        np.einsum('acij, cdkl, dbm->abijklm', dt_ij, t_ij, t_i) +
        np.einsum('acij, cdkl, dbm->abijklm', t_ij, dt_ij, t_i) +
        np.einsum('acij, cdkl, dbm->abijklm', t_ij, t_ij, dt_i)
    )

    return un_linked_disconnected_terms

def _order_6_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, dt_ijkl, dt_ijklm, *unusedargs = dt_args
    # Creating the 6th order return array
    un_linked_disconnected_terms = np.zeros((A, A, N, N, N, N, N, N), dtype=complex)
        # the (4, 2) term
    un_linked_disconnected_terms += 1/(factorial(2) * factorial(4) * factorial(2)) * (
        np.einsum('acij, cbklmn->abijklmn', dt_ij, t_ijkl) +
        np.einsum('acij, cbklmn->abijklmn', t_ij, dt_ijkl) +
        np.einsum('acijkl, cbmn->abijklmn', dt_ijkl, t_ij) +
        np.einsum('acijkl, cbmn->abijklmn', t_ijkl, dt_ij)
    )
        # the (3, 3) term
    un_linked_disconnected_terms += 1/(factorial(2) * factorial(3) * factorial(3)) * (
        np.einsum('acijk, cblmn->abijklmn', dt_ijk, t_ijk) +
        np.einsum('acijk, cblmn->abijklmn', t_ijk, dt_ijk)
    )
        # the (3, 2, 1) term
    un_linked_disconnected_terms += 1/(factorial(3) * factorial(3) * factorial(2)) * (
        np.einsum('aci, cdjk, dblmn->abijklmn', dt_i, t_ij, t_ijk) +
        np.einsum('aci, cdjk, dblmn->abijklmn', t_i, dt_ij, t_ijk) +
        np.einsum('aci, cdjk, dblmn->abijklmn', t_i, t_ij, dt_ijk) +
        np.einsum('aci, cdjkl, dbmn->abijklmn', dt_i, t_ijk, t_ij) +
        np.einsum('aci, cdjkl, dbmn->abijklmn', t_i, dt_ijk, t_ij) +
        np.einsum('aci, cdjkl, dbmn->abijklmn', t_i, t_ijk, dt_ij) +
        np.einsum('acij, cdk, dblmn->abijklmn', dt_ij, t_i, t_ijk) +
        np.einsum('acij, cdk, dblmn->abijklmn', t_ij, dt_i, t_ijk) +
        np.einsum('acij, cdk, dblmn->abijklmn', t_ij, t_i, dt_ijk) +
        np.einsum('acij, cdklm, dbn->abijklmn', dt_ij, t_ijk, t_i) +
        np.einsum('acij, cdklm, dbn->abijklmn', t_ij, dt_ijk, t_i) +
        np.einsum('acij, cdklm, dbn->abijklmn', t_ij, t_ijk, dt_i) +
        np.einsum('acijk, cdl, dbmn->abijklmn', dt_ijk, t_i, t_ij) +
        np.einsum('acijk, cdl, dbmn->abijklmn', t_ijk, dt_i, t_ij) +
        np.einsum('acijk, cdl, dbmn->abijklmn', t_ijk, t_i, dt_ij) +
        np.einsum('acijk, cdlm, dbn->abijklmn', dt_ijk, t_ij, t_i) +
        np.einsum('acijk, cdlm, dbn->abijklmn', t_ijk, dt_ij, t_i) +
        np.einsum('acijk, cdlm, dbn->abijklmn', t_ijk, t_ij, dt_i)
    )
        # the (2, 2, 2) term
    un_linked_disconnected_terms += 1/(factorial(3) * factorial(2) * factorial(2) * factorial(2)) * (
        np.einsum('acij, cdkl, dbmn->abijklmn', dt_ij, t_ij, t_ij) +
        np.einsum('acij, cdkl, dbmn->abijklmn', t_ij, dt_ij, t_ij) +
        np.einsum('acij, cdkl, dbmn->abijklmn', t_ij, t_ij, dt_ij)
    )
        # the (2, 2, 1, 1) term
    un_linked_disconnected_terms += 1/(factorial(4) * factorial(2) * factorial(2)) * (
        np.einsum('aci, cdj, dekl, ebmn->abijklmn', dt_i, t_i, t_ij, t_ij) +
        np.einsum('aci, cdj, dekl, ebmn->abijklmn', t_i, dt_i, t_ij, t_ij) +
        np.einsum('aci, cdj, dekl, ebmn->abijklmn', t_i, t_i, dt_ij, t_ij) +
        np.einsum('aci, cdj, dekl, ebmn->abijklmn', t_i, t_i, t_ij, dt_ij) +
        np.einsum('aci, cdjk, del, ebmn->abijklmn', dt_i, t_ij, t_i, t_ij) +
        np.einsum('aci, cdjk, del, ebmn->abijklmn', t_i, dt_ij, t_i, t_ij) +
        np.einsum('aci, cdjk, del, ebmn->abijklmn', t_i, t_ij, dt_i, t_ij) +
        np.einsum('aci, cdjk, del, ebmn->abijklmn', t_i, t_ij, t_i, dt_ij) +
        np.einsum('aci, cdjk, delm, ebn->abijklmn', dt_i, t_ij, t_ij, t_i) +
        np.einsum('aci, cdjk, delm, ebn->abijklmn', t_i, dt_ij, t_ij, t_i) +
        np.einsum('aci, cdjk, delm, ebn->abijklmn', t_i, t_ij, dt_ij, t_i) +
        np.einsum('aci, cdjk, delm, ebn->abijklmn', t_i, t_ij, t_ij, dt_i) +
        np.einsum('acij, cdk, del, ebmn->abijklmn', dt_ij, t_i, t_i, t_ij) +
        np.einsum('acij, cdk, del, ebmn->abijklmn', t_ij, dt_i, t_i, t_ij) +
        np.einsum('acij, cdk, del, ebmn->abijklmn', t_ij, t_i, dt_i, t_ij) +
        np.einsum('acij, cdk, del, ebmn->abijklmn', t_ij, t_i, t_i, dt_ij) +
        np.einsum('acij, cdk, delm, ebn->abijklmn', dt_ij, t_i, t_ij, t_i) +
        np.einsum('acij, cdk, delm, ebn->abijklmn', t_ij, dt_i, t_ij, t_i) +
        np.einsum('acij, cdk, delm, ebn->abijklmn', t_ij, t_i, dt_ij, t_i) +
        np.einsum('acij, cdk, delm, ebn->abijklmn', t_ij, t_i, t_ij, dt_i) +
        np.einsum('acij, cdkl, dem, ebn->abijklmn', dt_ij, t_ij, t_i, t_i) +
        np.einsum('acij, cdkl, dem, ebn->abijklmn', t_ij, dt_ij, t_i, t_i) +
        np.einsum('acij, cdkl, dem, ebn->abijklmn', t_ij, t_ij, dt_i, t_i) +
        np.einsum('acij, cdkl, dem, ebn->abijklmn', t_ij, t_ij, t_i, dt_i)
    )

    return un_linked_disconnected_terms

# ---------------------------- dt AMPLITUDES ---------------------------- #

def _calculate_order_1_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Calculate the derivative of the 1 t-amplitude for use in the calculation of the residuals."""
    # unpack the `w_args`
    w_i, *unusedargs = w_args
    # Calculate the 1st order residual
    residual = residual_equations.calculate_order_1_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * np.einsum('aci,cb->abi', w_i, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        pass  # no linked disconnected terms for order < 2
    elif ansatz.VECC:
        pass  # no un-linked disconnected terms for order < 4

    # Symmetrize the residual operator
    dt_i = symmetrize_tensor(N, residual, order=1)
    return dt_i

def _calculate_order_2_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Calculate the derivative of the 2 t-amplitude for use in the calculation of the residuals."""
    # unpack the `w_args`
    w_i, w_ij, *unusedargs = w_args
    # Calculate the 2nd order residual
    residual = residual_equations.calculate_order_2_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * np.einsum('acij,cb->abij', w_ij, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_2_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
    elif ansatz.VECC:
        residual -= _order_2_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
        pass  # no un-linked disconnected terms for order < 4

    # Symmetrize the residual operator
    dt_ij = symmetrize_tensor(N, residual, order=2)
    return dt_ij

def _calculate_order_3_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Calculate the derivative of the 3 t-amplitude for use in the calculation of the residuals."""
    # unpack the `w_args`
    w_i, w_ij, w_ijk, *unusedargs = w_args
    # Calculate the 3rd order residual
    residual = residual_equations.calculate_order_3_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * np.einsum('acijk,cb->abijk', w_ijk, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_3_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
    elif ansatz.VECC:
        residual -= _order_3_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
        pass  # no un-linked disconnected terms for order < 4

    # Symmetrize the residual operator
    dt_ijk = symmetrize_tensor(N, residual, order=3)
    return dt_ijk

def _calculate_order_4_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Calculate the derivative of the 4 t-amplitude for use in the calculation of the residuals."""
    # unpack the `w_args`
    w_i, w_ij, w_ijk, w_ijkl, *unusedargs = w_args
    # Calculate the 4th order residual
    residual = residual_equations.calculate_order_4_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * np.einsum('acijkl,cb->abijkl', w_ijkl, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_4_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
    elif ansatz.VECC:
        residual -= _order_4_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
        residual -= _order_4_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args)

    # Symmetrize the residual operator
    dt_ijkl = symmetrize_tensor(N, residual, order=4)
    return dt_ijkl

def _calculate_order_5_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Calculate the derivative of the 5 t-amplitude for use in the calculation of the residuals."""
    # unpack the `w_args`
    w_i, w_ij, w_ijk, w_ijkl, w_ijklm, *unusedargs = w_args
    # Calculate the 5th order residual
    residual = residual_equations.calculate_order_5_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * np.einsum('acijklm,cb->abijklm', w_ijklm, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_5_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
    elif ansatz.VECC:
        residual -= _order_5_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
        residual -= _order_5_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args)

    # Symmetrize the residual operator
    dt_ijklm = symmetrize_tensor(N, residual, order=5)
    return dt_ijklm

def _calculate_order_6_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Calculate the derivative of the 6 t-amplitude for use in the calculation of the residuals."""
    # unpack the `w_args`
    w_i, w_ij, w_ijk, w_ijkl, w_ijklm, w_ijklmn, *unusedargs = w_args
    # Calculate the 6th order residual
    residual = residual_equations.calculate_order_6_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * np.einsum('acijklmn,cb->abijklmn', w_ijklmn, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_6_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
    elif ansatz.VECC:
        residual -= _order_6_linked_disconnected_terms(A, N, trunc, t_args, dt_args)
        residual -= _order_6_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args)

    # Symmetrize the residual operator
    dt_ijklmn = symmetrize_tensor(N, residual, order=6)
    return dt_ijklmn

# ---------------------------- WRAPPER FUNCTIONS ---------------------------- #

def solve_singles_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_i term (singles)"""

    if not trunc.singles:
        raise Exception(
            "It appears that singles is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_i = _calculate_order_1_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_i

def solve_doubles_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ij term (doubles)"""

    if not trunc.doubles:
        raise Exception(
            "It appears that doubles is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ij = _calculate_order_2_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ij

def solve_triples_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ijk term (triples)"""

    if not trunc.triples:
        raise Exception(
            "It appears that triples is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ijk = _calculate_order_3_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ijk

def solve_quadruples_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ijkl term (quadruples)"""

    if not trunc.quadruples:
        raise Exception(
            "It appears that quadruples is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ijkl = _calculate_order_4_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ijkl

def solve_quintuples_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ijklm term (quintuples)"""

    if not trunc.quintuples:
        raise Exception(
            "It appears that quintuples is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ijklm = _calculate_order_5_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ijklm

def solve_sextuples_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ijklmn term (sextuples)"""

    if not trunc.sextuples:
        raise Exception(
            "It appears that sextuples is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ijklmn = _calculate_order_6_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ijklmn

# --------------------------------------------------------------------------- #
# --------------------------- OPTIMIZED FUNCTIONS --------------------------- #
# --------------------------------------------------------------------------- #

# ---------------------------- DISCONNECTED TERMS ---------------------------- #

def _order_1_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    raise Exception(
        "the first possible purely VECI/CC term is (1, 1) or (dt_i * t_i)"
        "which requires a residual of at least 2nd order"
    )

def _order_2_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args` and 'dt_args'
    t_i, *unusedargs = t_args
    dt_i, *unusedargs = dt_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Creating the 2nd order return array
    linked_disconnected_terms = np.zeros((A, A, N, N), dtype=complex)
    # the (1, 1) term
    linked_disconnected_terms += 1/factorial(2) * (
        next(optimized_einsum)(dt_i, t_i) +
        next(optimized_einsum)(t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_3_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, *unusedargs = t_args
    dt_i, dt_ij, *unusedargs = dt_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Creating the 3rd order return array
    linked_disconnected_terms = np.zeros((A, A, N, N, N), dtype=complex)
    # the (2, 1) term
    linked_disconnected_terms += 1/(factorial(2) * factorial(2)) * (
        next(optimized_einsum)(dt_i, t_ij) +
        next(optimized_einsum)(t_i, dt_ij) +
        next(optimized_einsum)(dt_ij, t_i) +
        next(optimized_einsum)(t_ij, dt_i)
    )
    # the (1, 1, 1) term
    linked_disconnected_terms += 1/factorial(3) * (
        next(optimized_einsum)(dt_i, t_i, t_i) +
        next(optimized_einsum)(t_i, dt_i, t_i) +
        next(optimized_einsum)(t_i, t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_4_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, *unusedargs = dt_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Creating the 4th order return array
    linked_disconnected_terms = np.zeros((A, A, N, N, N, N), dtype=complex)
    # the (3, 1) term
    linked_disconnected_terms += 1/(factorial(2) * factorial(3)) * (
        next(optimized_einsum)(dt_i, t_ijk) +
        next(optimized_einsum)(t_i, dt_ijk) +
        next(optimized_einsum)(dt_ijk, t_i) +
        next(optimized_einsum)(t_ijk, dt_i)
    )
    # the (2, 1, 1) term
    linked_disconnected_terms += 1/(factorial(3) * factorial(2)) * (
        next(optimized_einsum)(dt_i, t_i, t_ij) +
        next(optimized_einsum)(t_i, dt_i, t_ij) +
        next(optimized_einsum)(t_i, t_i, dt_ij) +
        next(optimized_einsum)(dt_i, t_ij, t_i) +
        next(optimized_einsum)(t_i, dt_ij, t_i) +
        next(optimized_einsum)(t_i, t_ij, dt_i) +
        next(optimized_einsum)(dt_ij, t_i, t_i) +
        next(optimized_einsum)(t_ij, dt_i, t_i) +
        next(optimized_einsum)(t_ij, t_i, dt_i)
    )
    # the (1, 1, 1, 1) term
    linked_disconnected_terms += 1/factorial(4) * (
        next(optimized_einsum)(dt_i, t_i, t_i, t_i) +
        next(optimized_einsum)(t_i, dt_i, t_i, t_i) +
        next(optimized_einsum)(t_i, t_i, dt_i, t_i) +
        next(optimized_einsum)(t_i, t_i, t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_5_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, dt_ijkl, *unusedargs = dt_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Creating the 5th order return array
    linked_disconnected_terms = np.zeros((A, A, N, N, N, N, N), dtype=complex)
    # the (4, 1) term
    linked_disconnected_terms += 1/(factorial(2) * factorial(4)) * (
        next(optimized_einsum)(dt_i, t_ijkl) +
        next(optimized_einsum)(t_i, dt_ijkl) +
        next(optimized_einsum)(dt_ijkl, t_i) +
        next(optimized_einsum)(t_ijkl, dt_i)
    )
    # the (3, 1, 1) term
    linked_disconnected_terms += 1/(factorial(3) * factorial(3)) * (
        next(optimized_einsum)(dt_i, t_i, t_ijk) +
        next(optimized_einsum)(t_i, dt_i, t_ijk) +
        next(optimized_einsum)(t_i, t_i, dt_ijk) +
        next(optimized_einsum)(dt_i, t_ijk, t_i) +
        next(optimized_einsum)(t_i, dt_ijk, t_i) +
        next(optimized_einsum)(t_i, t_ijk, dt_i) +
        next(optimized_einsum)(dt_ijk, t_i, t_i) +
        next(optimized_einsum)(t_ijk, dt_i, t_i) +
        next(optimized_einsum)(t_ijk, t_i, dt_i)
    )
    # the (2, 1, 1, 1) term
    linked_disconnected_terms += 1/(factorial(4) * factorial(2)) * (
        next(optimized_einsum)(dt_i, t_i, t_i, t_ij) +
        next(optimized_einsum)(t_i, dt_i, t_i, t_ij) +
        next(optimized_einsum)(t_i, t_i, dt_i, t_ij) +
        next(optimized_einsum)(t_i, t_i, t_i, dt_ij) +
        next(optimized_einsum)(dt_i, t_i, t_ij, t_i) +
        next(optimized_einsum)(t_i, dt_i, t_ij, t_i) +
        next(optimized_einsum)(t_i, t_i, dt_ij, t_i) +
        next(optimized_einsum)(t_i, t_i, t_ij, dt_i) +
        next(optimized_einsum)(dt_i, t_ij, t_i, t_i) +
        next(optimized_einsum)(t_i, dt_ij, t_i, t_i) +
        next(optimized_einsum)(t_i, t_ij, dt_i, t_i) +
        next(optimized_einsum)(t_i, t_ij, t_i, dt_i) +
        next(optimized_einsum)(dt_ij, t_i, t_i, t_i) +
        next(optimized_einsum)(t_ij, dt_i, t_i, t_i) +
        next(optimized_einsum)(t_ij, t_i, dt_i, t_i) +
        next(optimized_einsum)(t_ij, t_i, t_i, dt_i)
    )
    # the (1, 1, 1, 1, 1) term
    linked_disconnected_terms += 1/factorial(5) * (
        next(optimized_einsum)(dt_i, t_i, t_i, t_i, t_i) +
        next(optimized_einsum)(t_i, dt_i, t_i, t_i, t_i) +
        next(optimized_einsum)(t_i, t_i, dt_i, t_i, t_i) +
        next(optimized_einsum)(t_i, t_i, t_i, dt_i, t_i) +
        next(optimized_einsum)(t_i, t_i, t_i, t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_6_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1)
    But not terms (5), (3, 2), (2, 2, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, dt_ijkl, dt_ijklm, *unusedargs = dt_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Creating the 6th order return array
    linked_disconnected_terms = np.zeros((A, A, N, N, N, N, N, N), dtype=complex)
    # the (5, 1) term
    linked_disconnected_terms += 1/(factorial(2) * factorial(5)) * (
        next(optimized_einsum)(dt_i, t_ijklm) +
        next(optimized_einsum)(t_i, dt_ijklm) +
        next(optimized_einsum)(dt_ijklm, t_i) +
        next(optimized_einsum)(t_ijklm, dt_i)
    )
    # the (4, 1, 1) term
    linked_disconnected_terms += 1/(factorial(3) * factorial(4)) * (
        next(optimized_einsum)(dt_i, t_i, t_ijkl) +
        next(optimized_einsum)(t_i, dt_i, t_ijkl) +
        next(optimized_einsum)(t_i, t_i, dt_ijkl) +
        next(optimized_einsum)(dt_i, t_ijkl, t_i) +
        next(optimized_einsum)(t_i, dt_ijkl, t_i) +
        next(optimized_einsum)(t_i, t_ijkl, dt_i) +
        next(optimized_einsum)(dt_ijkl, t_i, t_i) +
        next(optimized_einsum)(t_ijkl, dt_i, t_i) +
        next(optimized_einsum)(t_ijkl, t_i, dt_i)
    )
    # the (3, 1, 1, 1) term
    linked_disconnected_terms += 1/(factorial(4) * factorial(3)) * (
        next(optimized_einsum)(dt_i, t_i, t_i, t_ijk) +
        next(optimized_einsum)(t_i, dt_i, t_i, t_ijk) +
        next(optimized_einsum)(t_i, t_i, dt_i, t_ijk) +
        next(optimized_einsum)(t_i, t_i, t_i, dt_ijk) +
        next(optimized_einsum)(dt_i, t_i, t_ijk, t_i) +
        next(optimized_einsum)(t_i, dt_i, t_ijk, t_i) +
        next(optimized_einsum)(t_i, t_i, dt_ijk, t_i) +
        next(optimized_einsum)(t_i, t_i, t_ijk, dt_i) +
        next(optimized_einsum)(dt_i, t_ijk, t_i, t_i) +
        next(optimized_einsum)(t_i, dt_ijk, t_i, t_i) +
        next(optimized_einsum)(t_i, t_ijk, dt_i, t_i) +
        next(optimized_einsum)(t_i, t_ijk, t_i, dt_i) +
        next(optimized_einsum)(dt_ijk, t_i, t_i, t_i) +
        next(optimized_einsum)(t_ijk, dt_i, t_i, t_i) +
        next(optimized_einsum)(t_ijk, t_i, dt_i, t_i) +
        next(optimized_einsum)(t_ijk, t_i, t_i, dt_i)
    )
    # the (2, 1, 1, 1, 1) term
    linked_disconnected_terms += 1/(factorial(5) * factorial(2)) * (
        next(optimized_einsum)(dt_i, t_i, t_i, t_i, t_ij) +
        next(optimized_einsum)(t_i, dt_i, t_i, t_i, t_ij) +
        next(optimized_einsum)(t_i, t_i, dt_i, t_i, t_ij) +
        next(optimized_einsum)(t_i, t_i, t_i, dt_i, t_ij) +
        next(optimized_einsum)(t_i, t_i, t_i, t_i, dt_ij) +
        next(optimized_einsum)(dt_i, t_i, t_i, t_ij, t_i) +
        next(optimized_einsum)(t_i, dt_i, t_i, t_ij, t_i) +
        next(optimized_einsum)(t_i, t_i, dt_i, t_ij, t_i) +
        next(optimized_einsum)(t_i, t_i, t_i, dt_ij, t_i) +
        next(optimized_einsum)(t_i, t_i, t_i, t_ij, dt_i) +
        next(optimized_einsum)(dt_i, t_i, t_ij, t_i, t_i) +
        next(optimized_einsum)(t_i, dt_i, t_ij, t_i, t_i) +
        next(optimized_einsum)(t_i, t_i, dt_ij, t_i, t_i) +
        next(optimized_einsum)(t_i, t_i, t_ij, dt_i, t_i) +
        next(optimized_einsum)(t_i, t_i, t_ij, t_i, dt_i) +
        next(optimized_einsum)(dt_i, t_ij, t_i, t_i, t_i) +
        next(optimized_einsum)(t_i, dt_ij, t_i, t_i, t_i) +
        next(optimized_einsum)(t_i, t_ij, dt_i, t_i, t_i) +
        next(optimized_einsum)(t_i, t_ij, t_i, dt_i, t_i) +
        next(optimized_einsum)(t_i, t_ij, t_i, t_i, dt_i) +
        next(optimized_einsum)(dt_ij, t_i, t_i, t_i, t_i) +
        next(optimized_einsum)(t_ij, dt_i, t_i, t_i, t_i) +
        next(optimized_einsum)(t_ij, t_i, dt_i, t_i, t_i) +
        next(optimized_einsum)(t_ij, t_i, t_i, dt_i, t_i) +
        next(optimized_einsum)(t_ij, t_i, t_i, t_i, dt_i)
    )
    # the (1, 1, 1, 1, 1, 1) term
    linked_disconnected_terms += 1/factorial(6) * (
        next(optimized_einsum)(dt_i, t_i, t_i, t_i, t_i, t_i) +
        next(optimized_einsum)(t_i, dt_i, t_i, t_i, t_i, t_i) +
        next(optimized_einsum)(t_i, t_i, dt_i, t_i, t_i, t_i) +
        next(optimized_einsum)(t_i, t_i, t_i, dt_i, t_i, t_i) +
        next(optimized_einsum)(t_i, t_i, t_i, t_i, dt_i, t_i) +
        next(optimized_einsum)(t_i, t_i, t_i, t_i, t_i, dt_i)
    )

    return linked_disconnected_terms

def _order_1_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a residual of at least 4th order"
    )

def _order_2_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a residual of at least 4th order"
    )

def _order_3_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a residual of at least 4th order"
    )

def _order_4_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, *unusedargs = dt_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Creating the 4th order return array
    un_linked_disconnected_terms = np.zeros((A, A, N, N, N, N), dtype=complex)
        # the (2, 2) term
    un_linked_disconnected_terms += 1/(factorial(2) * factorial(2) * factorial(2)) * (
        next(optimized_einsum)(dt_ij, t_ij) +
        next(optimized_einsum)(t_ij, dt_ij)
    )

    return un_linked_disconnected_terms

def _order_5_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, dt_ijkl, *unusedargs = dt_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Creating the 5th order return array
    un_linked_disconnected_terms = np.zeros((A, A, N, N, N, N, N), dtype=complex)
        # the (3, 2) term
    un_linked_disconnected_terms += 1/(factorial(2) * factorial(3) * factorial(2)) * (
        next(optimized_einsum)(dt_ij, t_ijk) +
        next(optimized_einsum)(t_ij, dt_ijk) +
        next(optimized_einsum)(dt_ijk, t_ij) +
        next(optimized_einsum)(t_ijk, dt_ij)
    )
        # the (2, 2, 1) term
    un_linked_disconnected_terms += 1/(factorial(3) * factorial(2) * factorial(2)) * (
        next(optimized_einsum)(dt_i, t_ij, t_ij) +
        next(optimized_einsum)(t_i, dt_ij, t_ij) +
        next(optimized_einsum)(t_i, t_ij, dt_ij) +
        next(optimized_einsum)(dt_ij, t_i, t_ij) +
        next(optimized_einsum)(t_ij, dt_i, t_ij) +
        next(optimized_einsum)(t_ij, t_i, dt_ij) +
        next(optimized_einsum)(dt_ij, t_ij, t_i) +
        next(optimized_einsum)(t_ij, dt_ij, t_i) +
        next(optimized_einsum)(t_ij, t_ij, dt_i)
    )

    return un_linked_disconnected_terms

def _order_6_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):
    """Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.
    This means for order 5 we include terms such as (3, 2), (2, 2, 1)
    But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args` and 'dt_args'
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, *unusedargs = t_args
    dt_i, dt_ij, dt_ijk, dt_ijkl, dt_ijklm, *unusedargs = dt_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Creating the 6th order return array
    un_linked_disconnected_terms = np.zeros((A, A, N, N, N, N, N, N), dtype=complex)
        # the (4, 2) term
    un_linked_disconnected_terms += 1/(factorial(2) * factorial(4) * factorial(2)) * (
        next(optimized_einsum)(dt_ij, t_ijkl) +
        next(optimized_einsum)(t_ij, dt_ijkl) +
        next(optimized_einsum)(dt_ijkl, t_ij) +
        next(optimized_einsum)(t_ijkl, dt_ij)
    )
        # the (3, 3) term
    un_linked_disconnected_terms += 1/(factorial(2) * factorial(3) * factorial(3)) * (
        next(optimized_einsum)(dt_ijk, t_ijk) +
        next(optimized_einsum)(t_ijk, dt_ijk)
    )
        # the (3, 2, 1) term
    un_linked_disconnected_terms += 1/(factorial(3) * factorial(3) * factorial(2)) * (
        next(optimized_einsum)(dt_i, t_ij, t_ijk) +
        next(optimized_einsum)(t_i, dt_ij, t_ijk) +
        next(optimized_einsum)(t_i, t_ij, dt_ijk) +
        next(optimized_einsum)(dt_i, t_ijk, t_ij) +
        next(optimized_einsum)(t_i, dt_ijk, t_ij) +
        next(optimized_einsum)(t_i, t_ijk, dt_ij) +
        next(optimized_einsum)(dt_ij, t_i, t_ijk) +
        next(optimized_einsum)(t_ij, dt_i, t_ijk) +
        next(optimized_einsum)(t_ij, t_i, dt_ijk) +
        next(optimized_einsum)(dt_ij, t_ijk, t_i) +
        next(optimized_einsum)(t_ij, dt_ijk, t_i) +
        next(optimized_einsum)(t_ij, t_ijk, dt_i) +
        next(optimized_einsum)(dt_ijk, t_i, t_ij) +
        next(optimized_einsum)(t_ijk, dt_i, t_ij) +
        next(optimized_einsum)(t_ijk, t_i, dt_ij) +
        next(optimized_einsum)(dt_ijk, t_ij, t_i) +
        next(optimized_einsum)(t_ijk, dt_ij, t_i) +
        next(optimized_einsum)(t_ijk, t_ij, dt_i)
    )
        # the (2, 2, 2) term
    un_linked_disconnected_terms += 1/(factorial(3) * factorial(2) * factorial(2) * factorial(2)) * (
        next(optimized_einsum)(dt_ij, t_ij, t_ij) +
        next(optimized_einsum)(t_ij, dt_ij, t_ij) +
        next(optimized_einsum)(t_ij, t_ij, dt_ij)
    )
        # the (2, 2, 1, 1) term
    un_linked_disconnected_terms += 1/(factorial(4) * factorial(2) * factorial(2)) * (
        next(optimized_einsum)(dt_i, t_i, t_ij, t_ij) +
        next(optimized_einsum)(t_i, dt_i, t_ij, t_ij) +
        next(optimized_einsum)(t_i, t_i, dt_ij, t_ij) +
        next(optimized_einsum)(t_i, t_i, t_ij, dt_ij) +
        next(optimized_einsum)(dt_i, t_ij, t_i, t_ij) +
        next(optimized_einsum)(t_i, dt_ij, t_i, t_ij) +
        next(optimized_einsum)(t_i, t_ij, dt_i, t_ij) +
        next(optimized_einsum)(t_i, t_ij, t_i, dt_ij) +
        next(optimized_einsum)(dt_i, t_ij, t_ij, t_i) +
        next(optimized_einsum)(t_i, dt_ij, t_ij, t_i) +
        next(optimized_einsum)(t_i, t_ij, dt_ij, t_i) +
        next(optimized_einsum)(t_i, t_ij, t_ij, dt_i) +
        next(optimized_einsum)(dt_ij, t_i, t_i, t_ij) +
        next(optimized_einsum)(t_ij, dt_i, t_i, t_ij) +
        next(optimized_einsum)(t_ij, t_i, dt_i, t_ij) +
        next(optimized_einsum)(t_ij, t_i, t_i, dt_ij) +
        next(optimized_einsum)(dt_ij, t_i, t_ij, t_i) +
        next(optimized_einsum)(t_ij, dt_i, t_ij, t_i) +
        next(optimized_einsum)(t_ij, t_i, dt_ij, t_i) +
        next(optimized_einsum)(t_ij, t_i, t_ij, dt_i) +
        next(optimized_einsum)(dt_ij, t_ij, t_i, t_i) +
        next(optimized_einsum)(t_ij, dt_ij, t_i, t_i) +
        next(optimized_einsum)(t_ij, t_ij, dt_i, t_i) +
        next(optimized_einsum)(t_ij, t_ij, t_i, dt_i)
    )

    return un_linked_disconnected_terms

# ---------------------------- dt AMPLITUDES ---------------------------- #

def _calculate_order_1_dt_amplitude_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_epsilon, opt_path_list):
    """Calculate the derivative of the 1 t-amplitude for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `w_args`
    w_i, *unusedargs = w_args
    # Calculate the 1st order residual
    residual = residual_equations.calculate_order_1_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * opt_epsilon(w_i, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        pass  # no linked disconnected terms for order < 2
    elif ansatz.VECC:
        pass  # no un-linked disconnected terms for order < 4

    # Symmetrize the residual operator
    dt_i = symmetrize_tensor(N, residual, order=1)
    return dt_i

def _calculate_order_2_dt_amplitude_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_epsilon, opt_path_list):
    """Calculate the derivative of the 2 t-amplitude for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `w_args`
    w_i, w_ij, *unusedargs = w_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Calculate the 2nd order residual
    residual = residual_equations.calculate_order_2_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * opt_epsilon(w_ij, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_2_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
    elif ansatz.VECC:
        residual -= _order_2_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
        pass  # no un-linked disconnected terms for order < 4

    # Symmetrize the residual operator
    dt_ij = symmetrize_tensor(N, residual, order=2)
    return dt_ij

def _calculate_order_3_dt_amplitude_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_epsilon, opt_path_list):
    """Calculate the derivative of the 3 t-amplitude for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `w_args`
    w_i, w_ij, w_ijk, *unusedargs = w_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Calculate the 3rd order residual
    residual = residual_equations.calculate_order_3_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * opt_epsilon(w_ijk, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_3_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
    elif ansatz.VECC:
        residual -= _order_3_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
        pass  # no un-linked disconnected terms for order < 4

    # Symmetrize the residual operator
    dt_ijk = symmetrize_tensor(N, residual, order=3)
    return dt_ijk

def _calculate_order_4_dt_amplitude_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_epsilon, opt_path_list):
    """Calculate the derivative of the 4 t-amplitude for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `w_args`
    w_i, w_ij, w_ijk, w_ijkl, *unusedargs = w_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Calculate the 4th order residual
    residual = residual_equations.calculate_order_4_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * opt_epsilon(w_ijkl, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_4_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
    elif ansatz.VECC:
        residual -= _order_4_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
        residual -= _order_4_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)

    # Symmetrize the residual operator
    dt_ijkl = symmetrize_tensor(N, residual, order=4)
    return dt_ijkl

def _calculate_order_5_dt_amplitude_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_epsilon, opt_path_list):
    """Calculate the derivative of the 5 t-amplitude for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `w_args`
    w_i, w_ij, w_ijk, w_ijkl, w_ijklm, *unusedargs = w_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Calculate the 5th order residual
    residual = residual_equations.calculate_order_5_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * opt_epsilon(w_ijklm, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_5_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
    elif ansatz.VECC:
        residual -= _order_5_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
        residual -= _order_5_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)

    # Symmetrize the residual operator
    dt_ijklm = symmetrize_tensor(N, residual, order=5)
    return dt_ijklm

def _calculate_order_6_dt_amplitude_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_epsilon, opt_path_list):
    """Calculate the derivative of the 6 t-amplitude for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `w_args`
    w_i, w_ij, w_ijk, w_ijkl, w_ijklm, w_ijklmn, *unusedargs = w_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # Calculate the 6th order residual
    residual = residual_equations.calculate_order_6_residual(A, N, trunc, h_args, w_args)
    # subtract the epsilon term (which is R_0)
    residual -= 1/factorial(2) * opt_epsilon(w_ijklmn, epsilon)

    # subtract the disconnected terms
    if ansatz.VECI:
        pass  # veci does not include any disconnected terms
    elif ansatz.VE_MIXED:
        residual -= _order_6_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
    elif ansatz.VECC:
        residual -= _order_6_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)
        residual -= _order_6_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)

    # Symmetrize the residual operator
    dt_ijklmn = symmetrize_tensor(N, residual, order=6)
    return dt_ijklmn

# ---------------------------- WRAPPER FUNCTIONS ---------------------------- #

def solve_singles_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_i term (singles)"""

    if not trunc.singles:
        raise Exception(
            "It appears that singles is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_i = _calculate_order_1_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_i

def solve_doubles_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ij term (doubles)"""

    if not trunc.doubles:
        raise Exception(
            "It appears that doubles is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ij = _calculate_order_2_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ij

def solve_triples_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ijk term (triples)"""

    if not trunc.triples:
        raise Exception(
            "It appears that triples is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ijk = _calculate_order_3_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ijk

def solve_quadruples_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ijkl term (quadruples)"""

    if not trunc.quadruples:
        raise Exception(
            "It appears that quadruples is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ijkl = _calculate_order_4_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ijkl

def solve_quintuples_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ijklm term (quintuples)"""

    if not trunc.quintuples:
        raise Exception(
            "It appears that quintuples is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ijklm = _calculate_order_5_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ijklm

def solve_sextuples_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
    """Compute the change in the t_ijklmn term (sextuples)"""

    if not trunc.sextuples:
        raise Exception(
            "It appears that sextuples is not true, this cannot be."
            "Something went terribly wrong!!!"
        )
    dt_ijklmn = _calculate_order_6_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
    return dt_ijklmn

# ---------------------------- OPTIMIZED PATHS FUNCTION ---------------------------- #

def compute_optimized_paths(A, N, truncation):
    """Calculate optimized paths for the einsum calls up to `highest_order`."""

    order_1_list, order_2_list, order_3_list = [], [], []
    order_4_list, order_5_list, order_6_list = [], [], []

    return [None]

