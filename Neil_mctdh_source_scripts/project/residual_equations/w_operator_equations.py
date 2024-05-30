# system imports
from math import factorial

# third party imports
import numpy as np
import opt_einsum as oe

# local imports
from .symmetrize import symmetrize_tensor
from ..log_conf import log

# --------------------------------------------------------------------------- #
# ---------------------------- DEFAULT FUNCTIONS ---------------------------- #
# --------------------------------------------------------------------------- #

# ---------------------------- VECI/CC CONTRIBUTIONS ---------------------------- #

def _add_order_1_vemx_contributions(W_1, t_args, truncation):
    """Exists for error checking."""
    raise Exception(
        "the first possible purely VECI/CC term is (1, 1) or (dt_i * t_i)"
        "which requires a W operator of at least 2nd order"
    )

def _add_order_2_vemx_contributions(W_2, t_args, truncation):
    """Calculate the order 2 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    """
    # unpack the `t_args`
    t_i, *unusedargs = t_args
    # SINGLES contribution
    W_2 += 1/factorial(2) * (np.einsum('aci, cbj->abij', t_i, t_i))
    return

def _add_order_3_vemx_contributions(W_3, t_args, truncation):
    """Calculate the order 3 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    """
    # unpack the `t_args`
    t_i, t_ij, *unusedargs = t_args
    # DOUBLES contribution
    if truncation.doubles:
        W_3 += 1/(factorial(2) * factorial(2)) * (
            np.einsum('aci, cbjk->abijk', t_i, t_ij) +
            np.einsum('acij, cbk->abijk', t_ij, t_i)
        )
    # SINGLES contribution
    W_3 += 1/factorial(3) * (np.einsum('aci, cdj, dbk->abijk', t_i, t_i, t_i))
    return

def _add_order_4_vemx_contributions(W_4, t_args, truncation):
    """Calculate the order 4 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, *unusedargs = t_args
    # TRIPLES contribution
    if truncation.triples:
        W_4 += 1/(factorial(2) * factorial(3)) * (
            np.einsum('aci, cbjkl->abijkl', t_i, t_ijk) +
            np.einsum('acijk, cbl->abijkl', t_ijk, t_i)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_4 += 1/(factorial(3) * factorial(2)) * (
            np.einsum('aci, cdj, dbkl->abijkl', t_i, t_i, t_ij) +
            np.einsum('aci, cdjk, dbl->abijkl', t_i, t_ij, t_i) +
            np.einsum('acij, cdk, dbl->abijkl', t_ij, t_i, t_i)
        )
    # SINGLES contribution
    W_4 += 1/factorial(4) * (np.einsum('aci, cdj, dek, ebl->abijkl', t_i, t_i, t_i, t_i))
    return

def _add_order_5_vemx_contributions(W_5, t_args, truncation):
    """Calculate the order 5 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    # QUADRUPLES contribution
    if truncation.quadruples:
        W_5 += 1/(factorial(2) * factorial(4)) * (
            np.einsum('aci, cbjklm->abijklm', t_i, t_ijkl) +
            np.einsum('acijkl, cbm->abijklm', t_ijkl, t_i)
        )
    # TRIPLES contribution
    if truncation.triples:
        W_5 += 1/(factorial(3) * factorial(3)) * (
            np.einsum('aci, cdj, dbklm->abijklm', t_i, t_i, t_ijk) +
            np.einsum('aci, cdjkl, dbm->abijklm', t_i, t_ijk, t_i) +
            np.einsum('acijk, cdl, dbm->abijklm', t_ijk, t_i, t_i)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_5 += 1/(factorial(4) * factorial(2)) * (
            np.einsum('aci, cdj, dek, eblm->abijklm', t_i, t_i, t_i, t_ij) +
            np.einsum('aci, cdj, dekl, ebm->abijklm', t_i, t_i, t_ij, t_i) +
            np.einsum('aci, cdjk, del, ebm->abijklm', t_i, t_ij, t_i, t_i) +
            np.einsum('acij, cdk, del, ebm->abijklm', t_ij, t_i, t_i, t_i)
        )
    # SINGLES contribution
    W_5 += 1/factorial(5) * (np.einsum('aci, cdj, dek, efl, fbm->abijklm', t_i, t_i, t_i, t_i, t_i))
    return

def _add_order_6_vemx_contributions(W_6, t_args, truncation):
    """Calculate the order 6 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, *unusedargs = t_args
    # QUINTUPLES contribution
    if truncation.quintuples:
        W_6 += 1/(factorial(2) * factorial(5)) * (
            np.einsum('aci, cbjklmn->abijklmn', t_i, t_ijklm) +
            np.einsum('acijklm, cbn->abijklmn', t_ijklm, t_i)
        )
    # QUADRUPLES contribution
    if truncation.quadruples:
        W_6 += 1/(factorial(3) * factorial(4)) * (
            np.einsum('aci, cdj, dbklmn->abijklmn', t_i, t_i, t_ijkl) +
            np.einsum('aci, cdjklm, dbn->abijklmn', t_i, t_ijkl, t_i) +
            np.einsum('acijkl, cdm, dbn->abijklmn', t_ijkl, t_i, t_i)
        )
    # TRIPLES contribution
    if truncation.triples:
        W_6 += 1/(factorial(4) * factorial(3)) * (
            np.einsum('aci, cdj, dek, eblmn->abijklmn', t_i, t_i, t_i, t_ijk) +
            np.einsum('aci, cdj, deklm, ebn->abijklmn', t_i, t_i, t_ijk, t_i) +
            np.einsum('aci, cdjkl, dem, ebn->abijklmn', t_i, t_ijk, t_i, t_i) +
            np.einsum('acijk, cdl, dem, ebn->abijklmn', t_ijk, t_i, t_i, t_i)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_6 += 1/(factorial(5) * factorial(2)) * (
            np.einsum('aci, cdj, dek, efl, fbmn->abijklmn', t_i, t_i, t_i, t_i, t_ij) +
            np.einsum('aci, cdj, dek, eflm, fbn->abijklmn', t_i, t_i, t_i, t_ij, t_i) +
            np.einsum('aci, cdj, dekl, efm, fbn->abijklmn', t_i, t_i, t_ij, t_i, t_i) +
            np.einsum('aci, cdjk, del, efm, fbn->abijklmn', t_i, t_ij, t_i, t_i, t_i) +
            np.einsum('acij, cdk, del, efm, fbn->abijklmn', t_ij, t_i, t_i, t_i, t_i)
        )
    # SINGLES contribution
    W_6 += 1/factorial(6) * (np.einsum('aci, cdj, dek, efl, fgm, gbn->abijklmn', t_i, t_i, t_i, t_i, t_i, t_i))
    return

# ---------------------------- VECC CONTRIBUTIONS ---------------------------- #

def _add_order_1_vecc_contributions(W_1, t_args, truncation):
    """Exists for error checking."""
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a W operator of at least 4th order"
    )

def _add_order_2_vecc_contributions(W_2, t_args, truncation):
    """Exists for error checking."""
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a W operator of at least 4th order"
    )

def _add_order_3_vecc_contributions(W_3, t_args, truncation):
    """Exists for error checking."""
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a W operator of at least 4th order"
    )

def _add_order_4_vecc_contributions(W_4, t_args, truncation):
    """Calculate the order 4 VECC contributions to the W operator
    for use in the calculation of the residuals.
    """
    # unpack the `t_args`
    t_i, t_ij, *unusedargs = t_args
    # DOUBLES contribution
    if truncation.doubles:
        W_4 += 1/(factorial(2) * factorial(2) * factorial(2)) * (
            np.einsum('acij, cbkl->abijkl', t_ij, t_ij)
        )
    return

def _add_order_5_vecc_contributions(W_5, t_args, truncation):
    """Calculate the order 5 VECC contributions to the W operator
    for use in the calculation of the residuals.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, *unusedargs = t_args
    # TRIPLES contribution
    if truncation.triples:
        W_5 += 1/(factorial(2) * factorial(3) * factorial(2)) * (
            np.einsum('acij, cbklm->abijklm', t_ij, t_ijk) +
            np.einsum('acijk, cblm->abijklm', t_ijk, t_ij)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_5 += 1/(factorial(3) * factorial(2) * factorial(2)) * (
            np.einsum('aci, cdjk, dblm->abijklm', t_i, t_ij, t_ij) +
            np.einsum('acij, cdk, dblm->abijklm', t_ij, t_i, t_ij) +
            np.einsum('acij, cdkl, dbm->abijklm', t_ij, t_ij, t_i)
        )
    return

def _add_order_6_vecc_contributions(W_6, t_args, truncation):
    """Calculate the order 6 VECC contributions to the W operator
    for use in the calculation of the residuals.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    # QUADRUPLES contribution
    if truncation.quadruples:
        W_6 += 1/(factorial(2) * factorial(4) * factorial(2)) * (
            np.einsum('acij, cbklmn->abijklmn', t_ij, t_ijkl) +
            np.einsum('acijkl, cbmn->abijklmn', t_ijkl, t_ij)
        )
    # TRIPLES contribution
    if truncation.triples:
        W_6 += 1/(factorial(2) * factorial(3) * factorial(3)) * (
            np.einsum('acijk, cblmn->abijklmn', t_ijk, t_ijk)
        )
    # TRIPLES contribution
    if truncation.triples:
        W_6 += 1/(factorial(3) * factorial(3) * factorial(2)) * (
            np.einsum('aci, cdjk, dblmn->abijklmn', t_i, t_ij, t_ijk) +
            np.einsum('aci, cdjkl, dbmn->abijklmn', t_i, t_ijk, t_ij) +
            np.einsum('acij, cdk, dblmn->abijklmn', t_ij, t_i, t_ijk) +
            np.einsum('acij, cdklm, dbn->abijklmn', t_ij, t_ijk, t_i) +
            np.einsum('acijk, cdl, dbmn->abijklmn', t_ijk, t_i, t_ij) +
            np.einsum('acijk, cdlm, dbn->abijklmn', t_ijk, t_ij, t_i)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_6 += 1/(factorial(3) * factorial(2) * factorial(2) * factorial(2)) * (
            np.einsum('acij, cdkl, dbmn->abijklmn', t_ij, t_ij, t_ij)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_6 += 1/(factorial(4) * factorial(2) * factorial(2)) * (
            np.einsum('aci, cdj, dekl, ebmn->abijklmn', t_i, t_i, t_ij, t_ij) +
            np.einsum('aci, cdjk, del, ebmn->abijklmn', t_i, t_ij, t_i, t_ij) +
            np.einsum('aci, cdjk, delm, ebn->abijklmn', t_i, t_ij, t_ij, t_i) +
            np.einsum('acij, cdk, del, ebmn->abijklmn', t_ij, t_i, t_i, t_ij) +
            np.einsum('acij, cdk, delm, ebn->abijklmn', t_ij, t_i, t_ij, t_i) +
            np.einsum('acij, cdkl, dem, ebn->abijklmn', t_ij, t_ij, t_i, t_i)
        )
    return

# ---------------------------- W OPERATOR FUNCTIONS ---------------------------- #

def _calculate_order_1_w_operator(A, N, t_args, ansatz, truncation):
    """Calculate the order 1 W operator for use in the calculation of the residuals."""
    # unpack the `t_args`
    t_i, *unusedargs = t_args
    # Creating the 1st order W operator
    W_1 = np.zeros((A, A, N), dtype=complex)
    # Singles contribution
    W_1 += t_i
    return W_1

def _calculate_order_2_w_operator(A, N, t_args, ansatz, truncation):
    """Calculate the order 2 W operator for use in the calculation of the residuals."""
    # unpack the `t_args`
    t_i, t_ij, *unusedargs = t_args
    # Creating the 2nd order W operator
    W_2 = np.zeros((A, A, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.doubles:
        W_2 += 1/factorial(2) * t_ij
    if ansatz.VE_MIXED:
        _add_order_2_vemx_contributions(W_2, t_args, truncation)
    elif ansatz.VECC:
        _add_order_2_vemx_contributions(W_2, t_args, truncation)
        pass  # no VECC contributions for order < 4

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_2, order=2)
    return symmetric_w

def _calculate_order_3_w_operator(A, N, t_args, ansatz, truncation):
    """Calculate the order 3 W operator for use in the calculation of the residuals."""
    # unpack the `t_args`
    t_i, t_ij, t_ijk, *unusedargs = t_args
    # Creating the 3rd order W operator
    W_3 = np.zeros((A, A, N, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.triples:
        W_3 += 1/factorial(3) * t_ijk
    if ansatz.VE_MIXED:
        _add_order_3_vemx_contributions(W_3, t_args, truncation)
    elif ansatz.VECC:
        _add_order_3_vemx_contributions(W_3, t_args, truncation)
        pass  # no VECC contributions for order < 4

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_3, order=3)
    return symmetric_w

def _calculate_order_4_w_operator(A, N, t_args, ansatz, truncation):
    """Calculate the order 4 W operator for use in the calculation of the residuals."""
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    # Creating the 4th order W operator
    W_4 = np.zeros((A, A, N, N, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.quadruples:
        W_4 += 1/factorial(4) * t_ijkl
    if ansatz.VE_MIXED:
        _add_order_4_vemx_contributions(W_4, t_args, truncation)
    elif ansatz.VECC:
        _add_order_4_vemx_contributions(W_4, t_args, truncation)
        _add_order_4_vecc_contributions(W_4, t_args, truncation)

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_4, order=4)
    return symmetric_w

def _calculate_order_5_w_operator(A, N, t_args, ansatz, truncation):
    """Calculate the order 5 W operator for use in the calculation of the residuals."""
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, *unusedargs = t_args
    # Creating the 5th order W operator
    W_5 = np.zeros((A, A, N, N, N, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.quintuples:
        W_5 += 1/factorial(5) * t_ijklm
    if ansatz.VE_MIXED:
        _add_order_5_vemx_contributions(W_5, t_args, truncation)
    elif ansatz.VECC:
        _add_order_5_vemx_contributions(W_5, t_args, truncation)
        _add_order_5_vecc_contributions(W_5, t_args, truncation)

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_5, order=5)
    return symmetric_w

def _calculate_order_6_w_operator(A, N, t_args, ansatz, truncation):
    """Calculate the order 6 W operator for use in the calculation of the residuals."""
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, t_ijklmn, *unusedargs = t_args
    # Creating the 6th order W operator
    W_6 = np.zeros((A, A, N, N, N, N, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.sextuples:
        W_6 += 1/factorial(6) * t_ijklmn
    if ansatz.VE_MIXED:
        _add_order_6_vemx_contributions(W_6, t_args, truncation)
    elif ansatz.VECC:
        _add_order_6_vemx_contributions(W_6, t_args, truncation)
        _add_order_6_vecc_contributions(W_6, t_args, truncation)

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_6, order=6)
    return symmetric_w

def compute_w_operators(A, N, t_args, ansatz, truncation):
    """Compute a number of W operators depending on the level of truncation."""

    if not truncation.singles:
        raise Exception(
            "It appears that `singles` is not true, this cannot be.\n"
            "Something went terribly wrong!!!\n\n"
            f"{truncation}\n"
        )

    w_1 = _calculate_order_1_w_operator(A, N, t_args, ansatz, truncation)
    w_2 = _calculate_order_2_w_operator(A, N, t_args, ansatz, truncation)
    w_3 = _calculate_order_3_w_operator(A, N, t_args, ansatz, truncation)

    if not truncation.doubles:
        return w_1, w_2, w_3, None, None, None
    else:
        w_4 = _calculate_order_4_w_operator(A, N, t_args, ansatz, truncation)

    if not truncation.triples:
        return w_1, w_2, w_3, w_4, None, None
    else:
        w_5 = _calculate_order_5_w_operator(A, N, t_args, ansatz, truncation)

    if not truncation.quadruples:
        return w_1, w_2, w_3, w_4, w_5, None
    else:
        w_6 = _calculate_order_6_w_operator(A, N, t_args, ansatz, truncation)

    if not truncation.quintuples:
        return w_1, w_2, w_3, w_4, w_5, w_6
    else:
        raise Exception(
            "Attempting to calculate W^7 operator (quintuples)\n"
            "This is currently not implemented!!\n"
        )

# --------------------------------------------------------------------------- #
# --------------------------- OPTIMIZED FUNCTIONS --------------------------- #
# --------------------------------------------------------------------------- #

# ---------------------------- VECI/CC CONTRIBUTIONS ---------------------------- #

def _add_order_1_vemx_contributions_optimized(W_1, t_args, truncation, opt_path_list):
    """Exists for error checking."""
    raise Exception(
        "the first possible purely VECI/CC term is (1, 1) or (dt_i * t_i)"
        "which requires a W operator of at least 2nd order"
    )

def _add_order_2_vemx_contributions_optimized(W_2, t_args, truncation, opt_path_list):
    """Calculate the order 2 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, *unusedargs = t_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # SINGLES contribution
    W_2 += 1/factorial(2) * (next(optimized_einsum)(t_i, t_i))
    return

def _add_order_3_vemx_contributions_optimized(W_3, t_args, truncation, opt_path_list):
    """Calculate the order 3 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, *unusedargs = t_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # DOUBLES contribution
    if truncation.doubles:
        W_3 += 1/(factorial(2) * factorial(2)) * (
            next(optimized_einsum)(t_i, t_ij) +
            next(optimized_einsum)(t_ij, t_i)
        )
    # SINGLES contribution
    W_3 += 1/factorial(3) * (next(optimized_einsum)(t_i, t_i, t_i))
    return

def _add_order_4_vemx_contributions_optimized(W_4, t_args, truncation, opt_path_list):
    """Calculate the order 4 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, *unusedargs = t_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # TRIPLES contribution
    if truncation.triples:
        W_4 += 1/(factorial(2) * factorial(3)) * (
            next(optimized_einsum)(t_i, t_ijk) +
            next(optimized_einsum)(t_ijk, t_i)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_4 += 1/(factorial(3) * factorial(2)) * (
            next(optimized_einsum)(t_i, t_i, t_ij) +
            next(optimized_einsum)(t_i, t_ij, t_i) +
            next(optimized_einsum)(t_ij, t_i, t_i)
        )
    # SINGLES contribution
    W_4 += 1/factorial(4) * (next(optimized_einsum)(t_i, t_i, t_i, t_i))
    return

def _add_order_5_vemx_contributions_optimized(W_5, t_args, truncation, opt_path_list):
    """Calculate the order 5 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # QUADRUPLES contribution
    if truncation.quadruples:
        W_5 += 1/(factorial(2) * factorial(4)) * (
            next(optimized_einsum)(t_i, t_ijkl) +
            next(optimized_einsum)(t_ijkl, t_i)
        )
    # TRIPLES contribution
    if truncation.triples:
        W_5 += 1/(factorial(3) * factorial(3)) * (
            next(optimized_einsum)(t_i, t_i, t_ijk) +
            next(optimized_einsum)(t_i, t_ijk, t_i) +
            next(optimized_einsum)(t_ijk, t_i, t_i)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_5 += 1/(factorial(4) * factorial(2)) * (
            next(optimized_einsum)(t_i, t_i, t_i, t_ij) +
            next(optimized_einsum)(t_i, t_i, t_ij, t_i) +
            next(optimized_einsum)(t_i, t_ij, t_i, t_i) +
            next(optimized_einsum)(t_ij, t_i, t_i, t_i)
        )
    # SINGLES contribution
    W_5 += 1/factorial(5) * (next(optimized_einsum)(t_i, t_i, t_i, t_i, t_i))
    return

def _add_order_6_vemx_contributions_optimized(W_6, t_args, truncation, opt_path_list):
    """Calculate the order 6 VECI/CC (mixed) contributions to the W operator
    for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, *unusedargs = t_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # QUINTUPLES contribution
    if truncation.quintuples:
        W_6 += 1/(factorial(2) * factorial(5)) * (
            next(optimized_einsum)(t_i, t_ijklm) +
            next(optimized_einsum)(t_ijklm, t_i)
        )
    # QUADRUPLES contribution
    if truncation.quadruples:
        W_6 += 1/(factorial(3) * factorial(4)) * (
            next(optimized_einsum)(t_i, t_i, t_ijkl) +
            next(optimized_einsum)(t_i, t_ijkl, t_i) +
            next(optimized_einsum)(t_ijkl, t_i, t_i)
        )
    # TRIPLES contribution
    if truncation.triples:
        W_6 += 1/(factorial(4) * factorial(3)) * (
            next(optimized_einsum)(t_i, t_i, t_i, t_ijk) +
            next(optimized_einsum)(t_i, t_i, t_ijk, t_i) +
            next(optimized_einsum)(t_i, t_ijk, t_i, t_i) +
            next(optimized_einsum)(t_ijk, t_i, t_i, t_i)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_6 += 1/(factorial(5) * factorial(2)) * (
            next(optimized_einsum)(t_i, t_i, t_i, t_i, t_ij) +
            next(optimized_einsum)(t_i, t_i, t_i, t_ij, t_i) +
            next(optimized_einsum)(t_i, t_i, t_ij, t_i, t_i) +
            next(optimized_einsum)(t_i, t_ij, t_i, t_i, t_i) +
            next(optimized_einsum)(t_ij, t_i, t_i, t_i, t_i)
        )
    # SINGLES contribution
    W_6 += 1/factorial(6) * (next(optimized_einsum)(t_i, t_i, t_i, t_i, t_i, t_i))
    return

# ---------------------------- VECC CONTRIBUTIONS ---------------------------- #

def _add_order_1_vecc_contributions_optimized(W_1, t_args, truncation, opt_path_list):
    """Exists for error checking."""
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a W operator of at least 4th order"
    )

def _add_order_2_vecc_contributions_optimized(W_2, t_args, truncation, opt_path_list):
    """Exists for error checking."""
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a W operator of at least 4th order"
    )

def _add_order_3_vecc_contributions_optimized(W_3, t_args, truncation, opt_path_list):
    """Exists for error checking."""
    raise Exception(
        "the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"
        "which requires a W operator of at least 4th order"
    )

def _add_order_4_vecc_contributions_optimized(W_4, t_args, truncation, opt_path_list):
    """Calculate the order 4 VECC contributions to the W operator
    "for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, *unusedargs = t_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # DOUBLES contribution
    if truncation.doubles:
        W_4 += 1/(factorial(2) * factorial(2) * factorial(2)) * (
            next(optimized_einsum)(t_ij, t_ij)
        )
    return

def _add_order_5_vecc_contributions_optimized(W_5, t_args, truncation, opt_path_list):
    """Calculate the order 5 VECC contributions to the W operator
    "for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, *unusedargs = t_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # TRIPLES contribution
    if truncation.triples:
        W_5 += 1/(factorial(2) * factorial(3) * factorial(2)) * (
            next(optimized_einsum)(t_ij, t_ijk) +
            next(optimized_einsum)(t_ijk, t_ij)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_5 += 1/(factorial(3) * factorial(2) * factorial(2)) * (
            next(optimized_einsum)(t_i, t_ij, t_ij) +
            next(optimized_einsum)(t_ij, t_i, t_ij) +
            next(optimized_einsum)(t_ij, t_ij, t_i)
        )
    return

def _add_order_6_vecc_contributions_optimized(W_6, t_args, truncation, opt_path_list):
    """Calculate the order 6 VECC contributions to the W operator
    "for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    # make an iterable out of the `opt_path_list`
    optimized_einsum = iter(opt_path_list)
    # QUADRUPLES contribution
    if truncation.quadruples:
        W_6 += 1/(factorial(2) * factorial(4) * factorial(2)) * (
            next(optimized_einsum)(t_ij, t_ijkl) +
            next(optimized_einsum)(t_ijkl, t_ij)
        )
    # TRIPLES contribution
    if truncation.triples:
        W_6 += 1/(factorial(2) * factorial(3) * factorial(3)) * (
            next(optimized_einsum)(t_ijk, t_ijk)
        )
    # TRIPLES contribution
    if truncation.triples:
        W_6 += 1/(factorial(3) * factorial(3) * factorial(2)) * (
            next(optimized_einsum)(t_i, t_ij, t_ijk) +
            next(optimized_einsum)(t_i, t_ijk, t_ij) +
            next(optimized_einsum)(t_ij, t_i, t_ijk) +
            next(optimized_einsum)(t_ij, t_ijk, t_i) +
            next(optimized_einsum)(t_ijk, t_i, t_ij) +
            next(optimized_einsum)(t_ijk, t_ij, t_i)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_6 += 1/(factorial(3) * factorial(2) * factorial(2) * factorial(2)) * (
            next(optimized_einsum)(t_ij, t_ij, t_ij)
        )
    # DOUBLES contribution
    if truncation.doubles:
        W_6 += 1/(factorial(4) * factorial(2) * factorial(2)) * (
            next(optimized_einsum)(t_i, t_i, t_ij, t_ij) +
            next(optimized_einsum)(t_i, t_ij, t_i, t_ij) +
            next(optimized_einsum)(t_i, t_ij, t_ij, t_i) +
            next(optimized_einsum)(t_ij, t_i, t_i, t_ij) +
            next(optimized_einsum)(t_ij, t_i, t_ij, t_i) +
            next(optimized_einsum)(t_ij, t_ij, t_i, t_i)
        )
    return

# ---------------------------- W OPERATOR FUNCTIONS ---------------------------- #

def _calculate_order_1_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_opt_path_list, vecc_opt_path_list):
    """Calculate the order 1 W operator for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, *unusedargs = t_args
    # Creating the 1st order W operator
    W_1 = np.zeros((A, A, N), dtype=complex)
    # Singles contribution
    W_1 += t_i
    return W_1

def _calculate_order_2_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_opt_path_list, vecc_opt_path_list):
    """Calculate the order 2 W operator for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, *unusedargs = t_args
    # Creating the 2nd order W operator
    W_2 = np.zeros((A, A, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.doubles:
        W_2 += 1/factorial(2) * t_ij
    if ansatz.VE_MIXED:
        _add_order_2_vemx_contributions_optimized(W_2, t_args, truncation, vemx_opt_path_list)
    elif ansatz.VECC:
        _add_order_2_vemx_contributions_optimized(W_2, t_args, truncation, vemx_opt_path_list)
        pass  # no VECC contributions for order < 4

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_2, order=2)
    return symmetric_w

def _calculate_order_3_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_opt_path_list, vecc_opt_path_list):
    """Calculate the order 3 W operator for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, *unusedargs = t_args
    # Creating the 3rd order W operator
    W_3 = np.zeros((A, A, N, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.triples:
        W_3 += 1/factorial(3) * t_ijk
    if ansatz.VE_MIXED:
        _add_order_3_vemx_contributions_optimized(W_3, t_args, truncation, vemx_opt_path_list)
    elif ansatz.VECC:
        _add_order_3_vemx_contributions_optimized(W_3, t_args, truncation, vemx_opt_path_list)
        pass  # no VECC contributions for order < 4

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_3, order=3)
    return symmetric_w

def _calculate_order_4_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_opt_path_list, vecc_opt_path_list):
    """Calculate the order 4 W operator for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, *unusedargs = t_args
    # Creating the 4th order W operator
    W_4 = np.zeros((A, A, N, N, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.quadruples:
        W_4 += 1/factorial(4) * t_ijkl
    if ansatz.VE_MIXED:
        _add_order_4_vemx_contributions_optimized(W_4, t_args, truncation, vemx_opt_path_list)
    elif ansatz.VECC:
        _add_order_4_vemx_contributions_optimized(W_4, t_args, truncation, vemx_opt_path_list)
        _add_order_4_vecc_contributions_optimized(W_4, t_args, truncation, vecc_opt_path_list)

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_4, order=4)
    return symmetric_w

def _calculate_order_5_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_opt_path_list, vecc_opt_path_list):
    """Calculate the order 5 W operator for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, *unusedargs = t_args
    # Creating the 5th order W operator
    W_5 = np.zeros((A, A, N, N, N, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.quintuples:
        W_5 += 1/factorial(5) * t_ijklm
    if ansatz.VE_MIXED:
        _add_order_5_vemx_contributions_optimized(W_5, t_args, truncation, vemx_opt_path_list)
    elif ansatz.VECC:
        _add_order_5_vemx_contributions_optimized(W_5, t_args, truncation, vemx_opt_path_list)
        _add_order_5_vecc_contributions_optimized(W_5, t_args, truncation, vecc_opt_path_list)

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_5, order=5)
    return symmetric_w

def _calculate_order_6_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_opt_path_list, vecc_opt_path_list):
    """Calculate the order 6 W operator for use in the calculation of the residuals.
    Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.
    """
    # unpack the `t_args`
    t_i, t_ij, t_ijk, t_ijkl, t_ijklm, t_ijklmn, *unusedargs = t_args
    # Creating the 6th order W operator
    W_6 = np.zeros((A, A, N, N, N, N, N, N), dtype=complex)

    # add the VECI contribution
    if truncation.sextuples:
        W_6 += 1/factorial(6) * t_ijklmn
    if ansatz.VE_MIXED:
        _add_order_6_vemx_contributions_optimized(W_6, t_args, truncation, vemx_opt_path_list)
    elif ansatz.VECC:
        _add_order_6_vemx_contributions_optimized(W_6, t_args, truncation, vemx_opt_path_list)
        _add_order_6_vecc_contributions_optimized(W_6, t_args, truncation, vecc_opt_path_list)

    # Symmetrize the W operator
    symmetric_w = symmetrize_tensor(N, W_6, order=6)
    return symmetric_w

def compute_w_operators_optimized(A, N, t_args, ansatz, truncation, vemx_optimized_paths, vecc_optimized_paths):
    """Compute a number of W operators depending on the level of truncation."""

    if not truncation.singles:
        raise Exception(
            "It appears that `singles` is not true, this cannot be.\n"
            "Something went terribly wrong!!!\n\n"
            f"{truncation}\n"
        )

    w_1 = _calculate_order_1_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_optimized_paths[0], vecc_optimized_paths[0])
    w_2 = _calculate_order_2_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_optimized_paths[1], vecc_optimized_paths[1])
    w_3 = _calculate_order_3_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_optimized_paths[2], vecc_optimized_paths[2])

    if not truncation.doubles:
        return w_1, w_2, w_3, None, None, None
    else:
        w_4 = _calculate_order_4_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_optimized_paths[3], vecc_optimized_paths[3])

    if not truncation.triples:
        return w_1, w_2, w_3, w_4, None, None
    else:
        w_5 = _calculate_order_5_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_optimized_paths[4], vecc_optimized_paths[4])

    if not truncation.quadruples:
        return w_1, w_2, w_3, w_4, w_5, None
    else:
        w_6 = _calculate_order_6_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_optimized_paths[5], vecc_optimized_paths[5])

    if not truncation.quintuples:
        return w_1, w_2, w_3, w_4, w_5, w_6
    else:
        raise Exception(
            "Attempting to calculate W^7 operator (quintuples)\n"
            "This is currently not implemented!!\n"
        )


# ---------------------------- OPTIMIZED PATHS FUNCTION ---------------------------- #


def compute_optimized_vecc_paths(A, N, truncation):
    """Calculate optimized paths for the VECC einsum calls up to `highest_order`."""

    order_1_list, order_2_list, order_3_list = [], [], []
    order_4_list, order_5_list, order_6_list = [], [], []

    if truncation.doubles:
        if truncation.triples:
            order_4_list.extend([
            ])

        order_4_list.extend([
        ])
    else:
        log.warning("Didn't calculate optimized VECC paths of the W^4 operator")

    if truncation.triples:
        if truncation.quadruples:
            order_5_list.extend([
            ])

        order_5_list.extend([
            oe.contract_expression('aci, cdjk, dblm->abijklm', (A, A, N), (A, A, N, N), (A, A, N, N)),
            oe.contract_expression('acij, cdk, dblm->abijklm', (A, A, N, N), (A, A, N), (A, A, N, N)),
            oe.contract_expression('acij, cdkl, dbm->abijklm', (A, A, N, N), (A, A, N, N), (A, A, N)),
        ])
    else:
        log.warning("Didn't calculate optimized VECC paths of the W^5 operator")

    if truncation.quadruples:
        if truncation.quintuples:
            order_6_list.extend([
            ])

        order_6_list.extend([
            oe.contract_expression('acijk, cblmn->abijklmn', (A, A, N, N, N), (A, A, N, N, N)),
            oe.contract_expression('aci, cdjk, dblmn->abijklmn', (A, A, N), (A, A, N, N), (A, A, N, N, N)),
            oe.contract_expression('aci, cdjkl, dbmn->abijklmn', (A, A, N), (A, A, N, N, N), (A, A, N, N)),
            oe.contract_expression('acij, cdk, dblmn->abijklmn', (A, A, N, N), (A, A, N), (A, A, N, N, N)),
            oe.contract_expression('acij, cdklm, dbn->abijklmn', (A, A, N, N), (A, A, N, N, N), (A, A, N)),
            oe.contract_expression('acijk, cdl, dbmn->abijklmn', (A, A, N, N, N), (A, A, N), (A, A, N, N)),
            oe.contract_expression('acijk, cdlm, dbn->abijklmn', (A, A, N, N, N), (A, A, N, N), (A, A, N)),
            oe.contract_expression('acij, cdkl, dbmn->abijklmn', (A, A, N, N), (A, A, N, N), (A, A, N, N)),
            oe.contract_expression('aci, cdj, dekl, ebmn->abijklmn', (A, A, N), (A, A, N), (A, A, N, N), (A, A, N, N)),
            oe.contract_expression('aci, cdjk, del, ebmn->abijklmn', (A, A, N), (A, A, N, N), (A, A, N), (A, A, N, N)),
            oe.contract_expression('aci, cdjk, delm, ebn->abijklmn', (A, A, N), (A, A, N, N), (A, A, N, N), (A, A, N)),
            oe.contract_expression('acij, cdk, del, ebmn->abijklmn', (A, A, N, N), (A, A, N), (A, A, N), (A, A, N, N)),
            oe.contract_expression('acij, cdk, delm, ebn->abijklmn', (A, A, N, N), (A, A, N), (A, A, N, N), (A, A, N)),
            oe.contract_expression('acij, cdkl, dem, ebn->abijklmn', (A, A, N, N), (A, A, N, N), (A, A, N), (A, A, N)),
        ])
    else:
        log.warning("Didn't calculate optimized VECC paths of the W^6 operator")

    return [order_1_list, order_4_list, order_5_list, order_6_list]

