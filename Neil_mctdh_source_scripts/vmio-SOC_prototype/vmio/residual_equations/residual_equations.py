# system imports

# third party imports
import numpy as np

# local imports
from .symmetrize import symmetrize_tensor


def calculate_order_0_residual(A, N, truncation, h_args, w_args):
    """Calculate the 0 order residual as a function of the W operators."""
    h_ab, h_abI, h_abi, h_abIj, h_abIJ, h_abij = h_args
    w_i, *unusedargs = w_args

    R = np.zeros((A, A), dtype=complex)

    assert truncation.singles, \
        f"Cannot calculate order 0 residual for {truncation.cc_truncation_order}"

    R += 1.0 * h_ab

    R += 1.0 * np.einsum('acm,cbm->ab', h_abI, w_i)

    if truncation.quadratic:
        if w_ij is not None:
            R += (1/2) * np.einsum('acmn,cbmn->ab', h_abIJ, w_ij)

    return R

def calculate_order_1_residual(A, N, truncation, h_args, w_args):
    """Calculate the 1 order residual as a function of the W operators."""
    h_ab, h_abI, h_abi, h_abIj, h_abIJ, h_abij = h_args
    w_i, w_ij, *unusedargs = w_args

    R = np.zeros((A, A, N), dtype=complex)

    assert truncation.singles, \
        f"Cannot calculate order 1 residual for {truncation.cc_truncation_order}"

    R += 1.0 * np.einsum('ac,cbi->abi', h_ab, w_i)
        # NO
        # R += np.einsum('ac,cbi->abi',h_ab,t_i)

    if w_ij is not None:
        R += 1.0 * np.einsum('acm,cbmi->abi', h_abI, w_ij)
        # R += np.einsum('ack,cbki->abi', h_abI, t_ij)

        # R += 0.5*   (np.einsum('ackl,cdk,dbli->abi', h_abIJ, t_i, t_ij)+\
        #              np.einsum('ackl,cdli,dbk->abi', h_abIJ, t_ij, t_i))

    # NO
    # R += 0.5 * (np.einsum('ack,cdk,dbi->abi',h_abI,t_i,t_i)+\
    #             np.einsum('ack,cdi,dbk->abi',h_abI,t_i,t_i))

    if truncation.quadratic:
        if w_ijk is not None:
            R += (3/6) * np.einsum('acmn,cbmni->abi', h_abIJ, w_ijk)

            # NO
            # R += 0.25 * (np.einsum('ackl,cdkl,dbi->abi',h_abIJ,t_ij,t_i)+\
            #              np.einsum('ackl,cdi,dbkl->abi',h_abIJ,t_i,t_ij))

    R += 1.0 * h_abi  # R += h_abi

    R += 1.0 * np.einsum('acmi,cbm->abi', h_abIj, w_i)  #R = np.einsum('acki,cbk->abi', h_abIj, t_i).copy()

    return R


    # unlinked terms




    R +=  1./12. * (np.einsum('ackl,cdk,del,ebi->abi',h_abIJ,t_i,t_i,t_i)+\
                    np.einsum('ackl,cdl,dek,ebi->abi',h_abIJ,t_i,t_i,t_i)+\
                    np.einsum('ackl,cdk,dei,ebl->abi',h_abIJ,t_i,t_i,t_i)+\
                    np.einsum('ackl,cdl,dei,ebk->abi',h_abIJ,t_i,t_i,t_i)+\
                    np.einsum('ackl,cdi,dek,ebl->abi',h_abIJ,t_i,t_i,t_i)+\
                    np.einsum('ackl,cdi,del,ebk->abi',h_abIJ,t_i,t_i,t_i))





def calculate_order_2_residual(A, N, truncation, h_args, w_args):
    """Calculate the 2 order residual as a function of the W operators."""
    h_ab, h_abI, h_abi, h_abIj, h_abIJ, h_abij = h_args
    w_i, w_ij, w_ijk, *unusedargs = w_args

    R = np.zeros((A, A, N, N), dtype=complex)

    assert truncation.doubles, \
        f"Cannot calculate order 2 residual for {truncation.cc_truncation_order}"

    if w_ij is not None:
        R += (1/2) * np.einsum('ac,cbij->abij', h_ab, w_ij)
        R += 1.0 * np.einsum('acmi,cbmj->abij', h_abIj, w_ij)

    if w_ijk is not None:
        R += (3/6) * np.einsum('acm,cbmij->abij', h_abI, w_ijk)

    if truncation.quadratic:
        if w_ijkl is not None:
            R += (6/24) * np.einsum('acmn,cbmnij->abij', h_abIJ, w_ijkl)
        else:
            R += (1/2) * h_abij

    R += 1.0 * np.einsum('aci,cbj->abij', h_abi, w_i)

    return R

def calculate_order_3_residual(A, N, truncation, h_args, w_args):
    """Calculate the 3 order residual as a function of the W operators."""
    h_ab, h_abI, h_abi, h_abIj, h_abIJ, h_abij = h_args
    w_i, w_ij, w_ijk, w_ijkl, *unusedargs = w_args

    R = np.zeros((A, A, N, N, N), dtype=complex)

    assert truncation.triples, \
        f"Cannot calculate order 3 residual for {truncation.cc_truncation_order}"

    if w_ijk is not None:
        R += (1/6) * np.einsum('ac,cbijk->abijk', h_ab, w_ijk)
        R += (3/6) * np.einsum('acmi,cbmjk->abijk', h_abIj, w_ijk)

    if w_ijkl is not None:
        R += (4/24) * np.einsum('acm,cbmijk->abijk', h_abI, w_ijkl)

    if truncation.quadratic:
        if w_ijklm is not None:
            R += (10/120) * np.einsum('acmn,cbmnijk->abijk', h_abIJ, w_ijklm)
        else:
            R += (1/2) * np.einsum('acij,cbk->abijk', h_abij, w_i)

    if w_ij is not None:
        R += (1/2) * np.einsum('aci,cbjk->abijk', h_abi, w_ij)

    return R

def calculate_order_4_residual(A, N, truncation, h_args, w_args):
    """Calculate the 4 order residual as a function of the W operators."""
    h_ab, h_abI, h_abi, h_abIj, h_abIJ, h_abij = h_args
    w_i, w_ij, w_ijk, w_ijkl, w_ijklm, *unusedargs = w_args

    R = np.zeros((A, A, N, N, N, N), dtype=complex)

    assert truncation.quadruples, \
        f"Cannot calculate order 4 residual for {truncation.cc_truncation_order}"

    if w_ijkl is not None:
        R += (1/24) * np.einsum('ac,cbijkl->abijkl', h_ab, w_ijkl)
        R += (4/24) * np.einsum('acmi,cbmjkl->abijkl', h_abIj, w_ijkl)

    if w_ijklm is not None:
        R += (5/120) * np.einsum('acm,cbmijkl->abijkl', h_abI, w_ijklm)

    if truncation.quadratic:
        if w_ijklmn is not None:
            R += (15/720) * np.einsum('acmn,cbmnijkl->abijkl', h_abIJ, w_ijklmn)
        if w_ij is not None:
            R += (1/(2*2)) * np.einsum('acij,cbkl->abijkl', h_abij, w_ij)

    if w_ijk is not None:
        R += (1/6) * np.einsum('aci,cbjkl->abijkl', h_abi, w_ijk)

    return R

def calculate_order_5_residual(A, N, truncation, h_args, w_args):
    """Calculate the 5 order residual as a function of the W operators."""
    h_ab, h_abI, h_abi, h_abIj, h_abIJ, h_abij = h_args
    w_i, w_ij, w_ijk, w_ijkl, w_ijklm, w_ijklmn, *unusedargs = w_args

    R = np.zeros((A, A, N, N, N, N, N), dtype=complex)

    assert truncation.quintuples, \
        f"Cannot calculate order 5 residual for {truncation.cc_truncation_order}"

    if w_ijklm is not None:
        R += (1/120) * np.einsum('ac,cbijkl->abijkl', h_ab, w_ijklm)
        R += (5/120) * np.einsum('acmi,cbmjkl->abijkl', h_abIj, w_ijklm)

    if w_ijklmn is not None:
        R += (6/720) * np.einsum('acm,cbmijkl->abijkl', h_abI, w_ijklmn)

    if truncation.quadratic:
        if w_ijklmno is not None:
            R += (21/5040) * np.einsum('acmn,cbmnijkl->abijkl', h_abIJ, w_ijklmno)
        if w_ijk is not None:
            R += (1/(2*2)) * np.einsum('acij,cbkl->abijkl', h_abij, w_ijk)

    if w_ijkl is not None:
        R += (1/24) * np.einsum('aci,cbjkl->abijkl', h_abi, w_ijkl)

    return R

def calculate_order_6_residual(A, N, truncation, h_args, w_args):
    """Calculate the 6 order residual as a function of the W operators."""
    h_ab, h_abI, h_abi, h_abIj, h_abIJ, h_abij = h_args
    w_i, w_ij, w_ijk, w_ijkl, w_ijklm, w_ijklmn, w_ijklmno, *unusedargs = w_args

    R = np.zeros((A, A, N, N, N, N, N, N), dtype=complex)

    assert truncation.sextuples, \
        f"Cannot calculate order 6 residual for {truncation.cc_truncation_order}"

    if w_ijklmn is not None:
        R += (1/720) * np.einsum('ac,cbijkl->abijkl', h_ab, w_ijklmn)
        R += (6/720) * np.einsum('acmi,cbmjkl->abijkl', h_abIj, w_ijklmn)

    if w_ijklmno is not None:
        R += (7/5040) * np.einsum('acm,cbmijkl->abijkl', h_abI, w_ijklmno)

    if truncation.quadratic:
        if w_ijklmnop is not None:
            R += (28/40320) * np.einsum('acmn,cbmnijkl->abijkl', h_abIJ, w_ijklmnop)
        if w_ijkl is not None:
            R += (1/(2*2)) * np.einsum('acij,cbkl->abijkl', h_abij, w_ijkl)

    if w_ijklm is not None:
        R += (1/120) * np.einsum('aci,cbjkl->abijkl', h_abi, w_ijklm)

    return R
