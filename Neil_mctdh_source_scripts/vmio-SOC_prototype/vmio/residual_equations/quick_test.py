# system imports
# import math
import itertools as it
from collections import namedtuple

# third party imports
import numpy as np

# local imports
import generate_residual_equations as gre


s_operator_namedtuple = namedtuple('s_operator', ['name', 'n'])

indices = 'ijklmnqrsv'

fbar = "\\overbar{f}"
fbar2 = "\\overbar{f}^{2}"
f___ = "f"
f___2 = "f^2"


def generate_omega_operator(maximum_rank=2):
    # m is associated with 'd' creation operators and n is associated with 'b' annihilation operators
    return_list = []

    for m in range(maximum_rank + 1):              # m is the upper label
        for n in range(maximum_rank + 1 - m):      # n is the lower label
            name = ""
            if m == 0 and n == 0:
                continue
            if m >= 1:
                name += "d"*m
            if n >= 1:
                name += "b"*n

            h_operator = gre.general_operator_namedtuple(name, m, n)
            return_list.append(h_operator)

    return gre.hamiltonian_namedtuple(maximum_rank, return_list)


def generate_S_operator_no_hotband(maximum_rank=2):
    return_list = []

    for m in range(maximum_rank + 1):      # m is the upper label

        if m == 0:
            s_operator = s_operator_namedtuple("1", 0)
        else:
            s_operator = s_operator_namedtuple(f"s^{m}", m)

        return_list.append(s_operator)

    return gre.hamiltonian_namedtuple(maximum_rank, return_list)


def generate_S_operator(maximum_rank=2):
    # m is associated with 'd' creation operators and n is associated with 'b' annihilation operators
    return_list = []

    for m in range(maximum_rank + 1):              # m is the upper label
        for n in range(maximum_rank + 1 - m):      # n is the lower label
            if m == 0 and n == 0:
                continue
            elif m == 0:
                s_operator = gre.general_operator_namedtuple(f"s_{n}", 0, n)
            elif n == 0:
                s_operator = gre.general_operator_namedtuple(f"s^{m}", m, 0)
            else:
                s_operator = gre.general_operator_namedtuple(f"s^{m}_{n}", m, n)

            return_list.append(s_operator)

    return gre.hamiltonian_namedtuple(maximum_rank, return_list)


def s_0_series_term(omega, H, s, term_list, total_list):
    """ fill up the `term_list` and `total_list` for the S^0 term"""


    print("\nS^0 term\n")


    for h in H.operator_list:

        nof_creation_ops = omega.m + h.m
        nof_annhiliation_ops = omega.n + h.n

        # only terms which can pair off all operators are non zero
        if nof_creation_ops != nof_annhiliation_ops:
            continue

        # realistically the fbars only appear based on b,d order
        # fbars = "fbar" if omega.name == 'db' else ""  # omega contribution
        fbars = "fbar" if h.m >= 1 else "" # h contribution

        term = ' * '.join([fbars, h.name])
        term_list.append(term)
        total_list.append([omega.name, h.name])

        print(f"nof operators: {nof_creation_ops}  {str(total_list[-1]):<28} {str(term_list[-1]):<28}")

    return


def s_1_series_term(omega, H, s_list, term_list, total_list):
    """ fill up the `term_list` and `total_list` for the S^1 term"""
    print("\nS^1 term\n")
    for s in s_list:
        for h in H.operator_list:

            nof_creation_ops = omega.m + h.m + s.m
            nof_annhiliation_ops = omega.n + h.n + s.n

            # only terms which can pair off all operators are non zero
            if nof_creation_ops != nof_annhiliation_ops:
                continue

            # realistically the fbars only appear based on b,d order
            # fbars = "fbar" if omega.name == 'db' else ""  # omega contribution
            fbars = "fbar" if h.m >= 1 else "" # h contribution

            term = ' * '.join([fbars, h.name, s.name.replace('s', 't')])
            term_list.append(term)
            total_list.append([omega.name, h.name, s.name])
            print(f"nof operators: {nof_creation_ops}  {str(total_list[-1]):<30} {str(term_list[-1]):<28}")

    return


def _s_n_logic(omega, H, s_series_term, term_list, total_list):
    total_list.append([omega.name, h.name, s.name])
    # print(total_list[-1])

    term = []
    h_idx = indices[0:h.m]

    print(f"S contributes {s.n} creation operators and {s.m} annihilation operators")
    print(f"H contributes {h.n} creation operators and {h.m} annihilation operators")
    print(f"W contributes {omega.n} creation operators and {omega.m} annihilation operators")

    if h.m == 0 and h.n == 0:
        term.append("h_0")
        if s.m >= 1:
            s_idx = indices[0:s.m]
            term.extend([f"fbar_{s.m}", f"t^{s_idx}"])
        if s.n >= 1:
            s_idx = indices[h.n:s.n]
            term.extend([f"fbar^{s.n}", f"t_{s_idx}"])

    elif h.m == 0:
        h_idx = indices[0:h.n]
        term.append(f"h_{h_idx}")
        if s.m >= 1:
            s_idx = indices[h.m:s.m]
            term.extend([f"fbar_{s.m}", f"t^{s_idx}"])
        if s.n >= 1:
            s_idx = indices[h.n:s.n]
            term.extend([f"fbar^{s.n}", f"t_{s_idx}"])


    elif h.n == 0:
        h_idx = indices[0:h.m]
        term.extend(["fbar", f"h^{h_idx}"])
        if s.m >= 1:
            s_idx = indices[h.m:s.m]
            term.extend([f"fbar_{s.m}", f"t^{s_idx}"])
        if s.n >= 1:
            s_idx = indices[h.n:s.n]
            term.extend([f"fbar^{s.n}", f"t_{s_idx}"])

    else:
        term.extend(["fbar", f"h_{h_idx}"])
        if s.m >= 1:
            s_idx = indices[h.m:s.m]
            term.extend([f"fbar_{s.m}", f"t^{s_idx}"])
        if s.n >= 1:
            s_idx = indices[h.n:s.n]
            term.extend([f"fbar^{s.n}", f"t_{s_idx}"])



    term = " * ".join(term)


    print(f"nof operators: {nof_creation_ops}  {str(total_list[-1]):<28} {term}")

    # print(f"nof creation: {nof_creation_ops}")
    # print(f"nof annhiliation: {nof_annhiliation_ops}")
    # print(f"f bars: {0}")
    # print(f"f unbar: {0}")

    # now we count the number of creation and annihilation operators
    return


def s_n_series_term(omega, H, s_series_term, term_list, total_list):
    """ fill up the `term_list` and `total_list` for the S^n term"""
    print("\nS^n term\n")
    for product in s_series_term:

        product_m = sum([s.m for s in product])
        product_n = sum([s.n for s in product])

        for h in H.operator_list:

            nof_creation_ops = omega.m + h.m + product_m
            nof_annhiliation_ops = omega.n + h.n + product_n

            # only terms which can pair off all operators are non zero
            if nof_creation_ops != nof_annhiliation_ops:
                continue

            # realistically the fbars only appear based on b,d order
            # fbars = "fbar" if omega.name == 'db' else ""  # omega contribution
            fbars = "fbar" if h.m - product_n >= 1 else "" # h contribution

            s_list = ' '.join([s.name.replace('s', 't') for s in product])
            term = ' * '.join([fbars, h.name, s_list])

            term_list.append(term)
            total_list.append([omega.name, h.name, s_list])

            print(f"nof operators: {nof_creation_ops}  {str(total_list[-1]):<32} {str(term_list[-1]):<28}")

    return


def main(taylor_truncation_order):

    omega_operator = generate_omega_operator(taylor_truncation_order)
    print(omega_operator)


    S = generate_S_operator()

    s_taylor_expansion = [None,]*3
    s_taylor_expansion[0] = gre.general_operator_namedtuple("1", 0, 0)
    s_taylor_expansion[1] = S.operator_list
    s_taylor_expansion[2] = list(it.product(S.operator_list, S.operator_list))

    # print('0', s_taylor_expansion[0])
    # print('1', s_taylor_expansion[1])
    # print('2', s_taylor_expansion[2][0])

    # generate the Hamiltonian data
    H = gre.generate_hamiltonian_operator()
    # for h in H.operator_list:
    #     print(h)

    term_list = []
    total_list = []

    count = 0
    for omega in omega_operator.operator_list:

        if omega.name == "d":
            continue

        for s_series_term in s_taylor_expansion:
            count += 1

            # S operator is simply 1 in this case
            if isinstance(s_series_term, gre.general_operator_namedtuple):
                s_0_series_term(omega, H, s_series_term, term_list, total_list)

            elif isinstance(s_series_term, list):

                # S^1 operator, simple
                if isinstance(s_series_term[0], gre.general_operator_namedtuple):
                    s_1_series_term(omega, H, s_series_term, term_list, total_list)

                # S^n operator, most complicated
                elif len(s_series_term[0]) >= 2:
                    s_n_series_term(omega, H, s_series_term, term_list, total_list)

                else:
                    raise Exception(f"Wrong S?\n {count=}, {s_series_term=}")

            else:
                raise Exception(f"s_series_term is not list?\n {count=}, {s_series_term=}")





    # print(np.array(total_list))
    return


if (__name__ == '__main__'):

    taylor_truncation_order = 1
    main(taylor_truncation_order)
