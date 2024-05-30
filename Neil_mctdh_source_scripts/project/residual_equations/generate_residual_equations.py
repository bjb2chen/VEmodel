# system imports
import math
import itertools as it
from collections import namedtuple

# third party imports
import numpy as np

# local imports


# define
tab = " "*4

# our Hamiltonian is
# H = (h_0 + omega + h^1 + h_1) + h^1_1 + h^2 + h_2
# but we can ignore the omega when calculating residuals as we add it back in at a later point
# so we use this H = h_0 + h^1 + h_1 + h^1_1 + h^2 + h_2

# ----------------------------------------------------------------------------------------------- #
# -------------------------------  NAMED TUPLES DEFINITIONS  ------------------------------------ #
# ----------------------------------------------------------------------------------------------- #
# the building blocks for the h & w components of each residual term are stored in named tuples

# rterm_namedtuple = namedtuple('rterm_namedtuple', ['prefactor', 'h', 'w'])
# we originally defined a class so that we can overload the `__eq__` operator
# because we needed to compare the rterm tuples, however I think I might have removed that code
# so Shanmei or I should check if we even need the class anymore
class residual_term(namedtuple('residual_term', ['prefactor', 'h', 'w'])):
    __slots__ = ()
    #
    def __eq__(self, other_term):
        return bool(
            self.prefactor == other_term.prefactor and \
            np.array_equal(self.h, other_term.h) and \
            np.array_equal(self.w, other_term.w)
        )


# tuple for a general operator as described on page 1, eq 1
general_operator_namedtuple = namedtuple('operator', ['name', 'm', 'n'])

hamiltonian_namedtuple = namedtuple('hamiltonian', ['maximum_rank', 'operator_list'])

"""rather than just lists and dictionaries using namedtuples makes the code much more concise
we can write things like `h.max_i` instead of `h[0]` and the label of the member explicitly
describes what value it contains making the code more readable and user friendly """
h_namedtuple = namedtuple('h_term', ['max_i', 'max_k'])
w_namedtuple = namedtuple('w_term', ['max_i', 'max_k', 'order'])

# used in building the W operators
t_term_namedtuple = namedtuple('t_term_namedtuple', ['string', 'order', 'shape'])

# ----------------------------------------------------------------------------------------------- #
# ------------------------------------  HELPER FUNCTIONS  --------------------------------------- #
# ----------------------------------------------------------------------------------------------- #

def unique_permutations(iterable):
    """Return a sorted list of unique permutations of the items in some iterable."""
    return sorted(list(set(it.permutations(iterable))))


def build_symmetrizing_function(max_order=5, show_perm=False):

    string = ""
    string += (
        f"\ndef symmetrize_tensor(tensor, order):\n"
        f"{tab}'''Symmetrizing a tensor (the W operator or the residual) by tracing over all permutations.'''\n"
        f"{tab}X = np.zeros_like(tensor, dtype=complex)\n"
    )

    string += (
        f"{tab}if order == 0:\n"
        f"{tab}{tab}return tensor\n"
        f"{tab}if order == 1:\n"
        f"{tab}{tab}return tensor\n"
    )

    for n in range(2, max_order+1):
        string += f"{tab}if order == {n}:\n"
        for p in it.permutations(range(2, n+2)):
            string += f"{tab}{tab}X += np.transpose(tensor, {(0,1) + p})\n"

    string += f"{tab}return X\n"
    return string

    X = np.zeros_like(W_2, dtype=complex)
    X += np.transpose(R, (0, 1, 2, 3))
    X += np.transpose(R, (0, 1, 3, 2))
    W_2 = X.copy()


def print_residual_data(R_lists, term_lists, print_equations = False, print_tuples = False):
    """Print to stdout in a easily readable format the residual terms and term tuples."""
    if print_equations:
        for i, R in enumerate(R_lists):
            print(f"{'':-<30} R_{i} {'':-<30}")
            for a in R:
                print(f"{tab} - {a}")
        print(f"{'':-<65}\n{'':-<65}\n")

    if print_tuples:
        for i, terms in enumerate(term_lists):
            print(f"{'':-<30} R_{i} {'':-<30}")
            for term in terms:
                print(f"{tab} - {term}")
        print(f"{'':-<65}\n{'':-<65}\n")

    return


def _partitions(number):
    """Return partitions if n. See `https://en.wikipedia.org/wiki/Partition_(number_theory)`"""
    answer = set()
    answer.add((number,))
    for x in range(1, number):
        for y in _partitions(number - x):
            answer.add(tuple(sorted((x, ) + y, reverse=True)))


    return sorted(list(answer), reverse=True)


def generate_partitions_of_n(n):
    """Return partitions of n. Such as (5,), (4, 1), (3, 1, 1), (2, 2, 1) ... etc."""
    return _partitions(n)


def generate_mixed_partitions_of_n(n):
    """Return partitions of n that include at most one number greater than 1.
    Such as (5,), (4, 1), (3, 1, 1), (2, 1, 1, 1) ... etc, but not (3, 2) or (2, 2, 1)
    """
    return [p for p in _partitions(n) if n - max(p) + 1 == len(p)]


def genereate_connected_partitions_of_n(n):
    """Return partitions of n which are only comprised of 1's.
    Such as (1, 1), or (1, 1, 1). The max value should only ever be 1.
    """
    return tuple([1]*n)


def generate_linked_disconnected_partitions_of_n(n):
    """Return partitions of n that include at most one number greater than 1 and not `n`.
    Such as (4, 1), (3, 1, 1), (2, 1, 1, 1) ... etc, but not (5,), (3, 2), (2, 2, 1)
    """
    return [p for p in _partitions(n) if n - max(p) + 1 == len(p) and max(p) < n]


def generate_un_linked_disconnected_partitions_of_n(n):
    """Return partitions of n that represent the unlinked disconnected wave operator parts.
    Such as (3, 2), (2, 2, 1) ... etc, but not (5,), (4, 1), (3, 1, 1), (2, 1, 1, 1)
    """
    new_set = set(_partitions(n)) - set(generate_mixed_partitions_of_n(n))
    return sorted(list(new_set), reverse=True)

# ----------------------------------------------------------------------------------------------- #
# --------------------------------  GENERATING RESIDUAL DATA  ----------------------------------- #
# ----------------------------------------------------------------------------------------------- #

def generate_hamiltonian_operator(maximum_rank=2):
    """Return a `hamiltonian_namedtuple`.
    It contains an `operator_list` of namedtuples for each term that looks like equation 6 based on `maximum_rank`
    (which is the sum of the ranks (m,n)).
    Equation 6 is a Hamiltonian with `maximum_rank` of 2
    """
    return_list = []

    for m in range(maximum_rank + 1):              # m is the upper label
        for n in range(maximum_rank + 1 - m):      # n is the lower label
            if m == 0 and n == 0:
                h_operator = general_operator_namedtuple("h_0", 0, 0)
            elif m == 0:
                h_operator = general_operator_namedtuple(f"h_{n}", 0, n)
            elif n == 0:
                h_operator = general_operator_namedtuple(f"h^{m}", m, 0)
            else:
                h_operator = general_operator_namedtuple(f"h^{m}_{n}", m, n)

            return_list.append(h_operator)

    return hamiltonian_namedtuple(maximum_rank, return_list)


# ------------- constructing the prefactor -------------- #

def extract_numerator_denominator_from_string(s):
    """Return the number part of the numerator and denominator from a string (s) of a fraction w factorial in denominator part."""

    if "*" in s: # we are only trying to get the first fraction
        s = s.split('*')[0]

    # print(s)
    numer, denom = s.replace('(', '').replace(')', '').split(sep="/")
    denom = denom.split(sep='!')[0]

    # print(f"numer: {numer} denom: {denom}")
    return [int(numer), int(denom)]


def simplified_prefactor(pre):
    """Creates the simplified form of the given prefactor f."""

    arr = extract_numerator_denominator_from_string(pre)
    numerator, denominator = arr[0], arr[1]

    if pre == "*(1/2)":    # case when 1/2 is the only prefactor, delete * sign
        pre = "1/2"
    elif denominator == numerator:
        if "*(1/2)" in pre:
            if denominator == 1 or denominator == 2:  # case when (1/1!) or (2/2!) presents, which will both be recognized as 1
                pre = "(1/2)"     # 1*(1/2) = (1/2)
            else:
                pre = f"({numerator}/(2*{denominator}))"
        else:
            if denominator == 1 or denominator == 2:
                pre = ""   # use empty string to present 1
    elif "*(1/2)" in pre:  # case when 1/2 is multiplied to the prefactor
        if numerator % 2 == 0:
            pre[0] = str(numerator // 2)
            pre = pre[:-6]  # get rid of "*(1/2)"
        else:
            pre = pre[:3] + "(2*" + pre[3:-6] + ")"  # add "2*" to the front of the denominator
    return pre


def construct_prefactor(h, p, simplify_flag=False):
    """Creates the string for the prefactor of the tensor h.
    Returns condensed fraction by default.
    If `simplify_flag` is true it reduces the fraction using the gcf(greatest common factor).
    If the prefactor is equal to 1 or there is no prefactor then returns an empty string.
    """

    # is there a 'master' equations / a general approach for any number of i's / k's

    # Should be able to simplify down to 3 cases? maybe?
    # case 1 - only 1 k, any number of i's or only 1 i, any number of k's
    # case 2 - 2 or more k's AND 2 or more i's
    # case 3 - only i or only k presents

    if h.m > p:
        return ""

    prefactor = ""
    i_number_w = p - h.m
    total_number_w = p - h.m + h.n

    # special case when there is no w operator and h^2 doesn't present
    if total_number_w == 0:
        if h.m == 2:
            return "(1/2)"
        else:
            return ""

    # case 1. when there is only 1 k or i label on w operator
    if i_number_w == 1 or h.n == 1:
        prefactor = f"({total_number_w}/{total_number_w}!)"  # needs simplification: n/n! = 1/(n-1)

    # case 2. when 2 or more k's AND 2 or more i's on w
    elif h.n > 1 and i_number_w > 1:
        prefactor +=  f"({sum(x for x in range(total_number_w))}/{total_number_w}!)"

    # case 3. when only i or only k presents on w operator
    elif h.n == 0 or i_number_w == 0:
        prefactor += f"(1/{total_number_w}!)"

    # special case: when h^2 is included, 1/2 needs to be multiplies to the term
    if h.m == 2:
        prefactor += "*(1/2)"

    # if simplification is needed
    if simplify_flag:            # call simplification step
        prefactor = simplified_prefactor(prefactor)

    return prefactor


# ---- construct string labels for each residual term --- #

def construct_upper_w_label(h, p):
    """Creates the string for the "upper" labels of the tensor w.
    Returns `^{str}` if there are upper labels
    Otherwise returns an empty string
    """
    if (h.m == p and h.n == 0) or h.m > p:   # case when there is no w operator needed
        return ""

    w_label = "^{"   # if w operator exist, initialize w_label
    for a in range(h.m+1, p+1):   # add i_p
        w_label += f"i_{a},"
    for b in range(1, h.n+1):     # add k_h
        w_label += f"k_{b},"

    assert w_label != "^{", "Whoops you missed a case, check the logic!!!"

    w_label = w_label[:-1] + "}"  # delete extra comma and add close bracket
    return w_label


def construct_upper_h_label(h, p):
    """Creates the string for the "upper" labels of the tensor h.
    Returns `^{str}` if there are upper labels
    Otherwise returns an empty string
    """
    if h.m == 0 or h.m > p:  # case when there is no upper i label or no proper h operator
        return ""

    upper_h_result = "^{"  # initialize the return string if upper label exist (m!=0)
    for c in range(1, h.m+1):
        upper_h_result += f"i_{c},"

    upper_h_result = upper_h_result[:-1] + "}"  # delete extra comma and add close bracket
    return upper_h_result


def construct_lower_h_label(h, p):
    """Creates the string for the "lower" labels of the tensor h.
    Returns `_{str}` if there are lower labels
    Otherwise returns an empty string
    """
    # case when h_o presents
    if h.m == 0 and h.n == 0:
        return "_0"

    # case when h operator doesn't have lower label
    if h.n == 0 or h.m > p:
        return ""

    # initialize the return string if lower label exist (n!=0)
    lower_h_result = "_{"
    for d in range(1, h.n+1):
        lower_h_result += f"k_{d},"

    lower_h_result = lower_h_result[:-1] + "}"
    return lower_h_result


# ------- constructing individual residual terms -------- #

def generate_p_term(str_fac):
    """Generate a floating point number calculated by the fraction in the input string fac_str"""

    # check if the prefactor is 1
    if str_fac == "":
        return "1.0"

    if str_fac == "(1/2) * ":
        return "0.5"

    arr = extract_numerator_denominator_from_string(str_fac)
    numerator, denominator = arr[0], arr[1]

    if "/(2*" in str_fac:
        return ("(" + str(numerator) + "/(2*" + str(math.factorial(denominator)) + "))")
    else:
        return ("(" + str(numerator) + "/" + str(math.factorial(denominator)) + ")")


def generate_h_term(str_h):
    """Generate an h_namedtuple; contains max_i and max_k of h operator"""
    if "0" in str_h:
        # special case for h_0
        return h_namedtuple(0, 0)

    return h_namedtuple(str_h.count("i"), str_h.count("k"))


def generate_w_term(str_w):
    """Generate an w_namedtuple; contains max_i and max_k of w operator"""
    if str_w == "":
        #if there is no w operator, return [0,0,0]
        return w_namedtuple(0, 0, 0)

    max_i = str_w.count("i")
    max_k = str_w.count("k")
    return w_namedtuple(max_i, max_k, max_i+max_k)

# ------------------------------------------------------- #

def generate_residual_string_list(hamiltonian, order):
    """Return a list of strings that will represent each term in the equations 59-63; based on order === p

    the indices i go from i_1, i_2 to i_order

    For each term in the Hamiltonian (see hamiltonian.png) we use the value of order === p
    to determine how much pairing happens when projecting (see residual_equation.png)
    such that we produce a string like the one we see in equation 59
    Refer back to operator_equation.png to remind yourself of the indices for each of the terms in the Hamiltonian
    """
    return_list = []
    term_list = []

    for h_operator in hamiltonian:

        # initialize the new hamiltonian operator, the w operator and the numeric part of each hamiltonian operator
        h_result, w_result, prefactor = "h", "w", ""

        # construct each hamiltonian operator and w operator
        h_result += construct_upper_h_label(h_operator, order)
        h_result += construct_lower_h_label(h_operator, order)
        w_result += construct_upper_w_label(h_operator, order)
        prefactor = construct_prefactor(h_operator,  order, True)

        if h_result == "h":
            continue

        # make terms here
        term_list += [residual_term(generate_p_term(prefactor), generate_h_term(h_result), generate_w_term(w_result)), ]

        # add '*' symbols for the text representation
        w_result = "" if (w_result == "w") else (" * " + w_result)
        if not (prefactor == ""):
            prefactor = prefactor + " * "

        return_list += [prefactor + h_result  + w_result, ]

    return return_list, term_list


def generate_residual_data(H, max_order):
    """Return two lists of length `max_order` containing the data on the residual equations.
    The `R_lists` contains string representations of the residual equations, to be used to generate latex documents.
    The `term_lists` contains tuple representations of the residuals to be used to generate python code for use in simulation.
    `term_lists` is a list of lists, each of which contain tuples `(prefactor, h, w)` representing terms of that particular residual.
    """
    lst = [tuple(generate_residual_string_list(H, order=order)) for order in range(max_order+1)]
    R_lists = [tup[0] for tup in lst]
    term_lists = [tup[1] for tup in lst]
    return R_lists, term_lists

# ----------------------------------------------------------------------------------------------- #
# -------------------------------  GENERATING RESIDUAL LATEX  ----------------------------------- #
# ----------------------------------------------------------------------------------------------- #

# use these if we define \newcommand to map textbf to \bt and so forth
if True:
    bold_t_latex = "\\bt"
    bold_h_latex = "\\bh"
    bold_w_latex = "\\bw"
    bold_c_latex = "\\bc"
    bold_d_latex = "\\bd"
else:
    bold_t_latex = "\\textbf{t}"
    bold_h_latex = "\\textbf{h}"
    bold_w_latex = "\\textbf{w}"
    bold_c_latex = "\\textbf{c}"
    bold_d_latex = "\\textbf{d}"


# ----------------- generating VECI residual latex -------------------- #
def _generate_veci_explicit_latex_term(term_list, h_list, w_string, order):
    """Generate the latex for the residual equations including all the indices and factors."""

    for h in h_list:
        h_string = bold_h_latex
        h_string += construct_upper_h_label(h, order)
        h_string += construct_lower_h_label(h, order)

        # no W operator
        if (h.m == order and h.n == 0) or h.m > order:
            w_string = ""

        # 1 W operator
        else:
            w_string = bold_t_latex + construct_upper_w_label(h, order)


        prefactor = construct_prefactor(h,  order, True)
        if prefactor != "":
            numerator, denominator = prefactor[1:-1].split('/')
            prefactor = f"\\left(\\dfrac{{{numerator}}}{{{denominator}}}\\right)"

        # save the string representation of the term
        # print(prefactor + h_string + w_string)
        term_list.append(prefactor + h_string + w_string)
    return


def _generate_veci_condensed_latex_term(term_list, h_list, w_string):
    """Append a string representation of a summation term for a residual to `term_list`."""

    if len(h_list) == 1:
        h_string = h_list[0].name.replace('h', bold_h_latex)
    else:
        bold_h_list = [h.name.replace('h', bold_h_latex) for h in h_list]
        h_string = f"({' + '.join(bold_h_list)})"

    # save the string representation of the term
    term_list.append(h_string + w_string)
    return


def _generate_veci_latex_form(hamiltonian, order, max_order, explicit=False):
    """Generate the VECI latex for the residual equations.
    Default is to generate the condensed form, but if `explicit` flag is `True` then
    explicit form with all i and k terms and prefactors is returned.
    """
    w_range = list(range(0, order + hamiltonian.maximum_rank + 1))

    # generate a list of H terms of equal or lower order than the residual
    operators = []
    for h in hamiltonian.operator_list:
        if h.m <= order:
            operators.append(h)

    # generate each term we need to sum over
    term_list = []
    for n in w_range:

        # create the w string
        w_string = f"{bold_t_latex}^{{{n}}}" if n != 0 else ""
        h_list = []

        # collect h terms that map to w^n
        for i, h in enumerate(operators):
            if abs(h.m - h.n - order) == n:
                h_list.append(h)
                operators.pop(i)

        # if no terms we skip this W operator
        if len(h_list) == 0:
            continue

        # append all term strings to `term_list`
        if explicit:
            _generate_veci_explicit_latex_term(term_list, h_list, w_string, order)
            # _generate_mixedcc_explicit_latex_term(term_list, h_list, w_list, order, expand_order_n_w)
        else:
            _generate_veci_condensed_latex_term(term_list, h_list, w_string)
            # _generate_mixedcc_condensed_latex_term(term_list, h_list, w_list)

        # loop


    # return the full string (with prefix and all term's wrapped in plus symbols)
    if explicit:
        i_terms = ", ".join([f"i_{n}" for n in range(1, order+1)]) if order != 0 else "0"
        return_string = f"\\begin{{split}} \\textbf{{R}}_{{{i_terms}}} &=\n" + " +\n".join(term_list) + " \\end{split}"
        return return_string
    else:
        return_string = f"\\textbf{{R}}_{{{order}}} &= " + " + ".join(term_list)
        return return_string


def _generate_veci_latex(hamiltonian, max_order):
    """Generate all of the VECI latex equations.
    We use `join` to insert two backward's slashes \\ BETWEEN each line
    rather then adding them to end and having extra trailing slashes on the last line.
    The user is expected to manually copy the relevant lines from the text file into a latex file
    and generate the pdf themselves.
    """

    return_string = "" # store it all in here
    return_string += "VECI condensed latex\n"
    return_string += ' \\\\\n%\n'.join([_generate_veci_latex_form(hamiltonian, order, max_order) for order in range(max_order+1)])
    return_string += "\n"*4
    return_string += "VECI explicit latex\n"
    return_string += ' \\\\\n%\n'.join([_generate_veci_latex_form(hamiltonian, order, max_order, explicit=True) for order in range(max_order+1)])
    return return_string


# ----------------- generating VECI/CC residual latex -------------------- #
def _construct_c_term(h, power, order):
    """Creates the string for the "upper" labels of the tensor c.
    c_n is defined as S * (1/n!) * (t^1)^n.
    Returns `c^{str}` if there are upper labels
    Otherwise returns an empty string
    """
    i_start = h.m+1
    nof_k = h.n

    upper_labels = []

    for n in range(i_start, order+1):
        upper_labels.append(f"i_{n}")

    for n in range(1, nof_k+1):
        upper_labels.append(f"k_{n}")

    return f"{bold_c_latex}^{{{', '.join(upper_labels)}}}"


def _construct_c_and_t_term(h, c_power, t_power, order):
    """Creates the string for the "upper" labels of the c and t tensors.
    c_n is defined as S * (1/n!) * (t^1)^n.
    Returns `c^{str1}t^{str2}` if there are upper labels
    Otherwise returns an empty string
    """
    nof_k = h.n

    string = ""

    return ""

    for n in range(1, c_power - nof_k + 1):
        string += f"{bold_c_latex}^{{i_{n}}}"

    string += f"{bold_t_latex}^{{k_{n}}}"

    return string


def _extract_power(string):
    """Return integer representation of the power that `string` is being raised to."""
    return int(string.split('{')[1].strip('}')[0])


def _extract_ct_powers(string):
    """Return integer representation of the powers that the two terms in `string` are being raised to."""
    n1 = int(string.split('{')[1].strip('}')[0])
    n2 = int(string.split('{')[2].strip('}')[0])
    return n1, n2


def _generate_mixedcc_explicit_latex_term(term_list, h_list, w_list, order, expand_order_n_w=None):
    """Generate the latex for the residual equations including all the indices and factors."""

    for h in h_list:
        h_string = bold_h_latex
        h_string += construct_upper_h_label(h, order)
        h_string += construct_lower_h_label(h, order)


        prefactor = construct_prefactor(h,  order, True)
        if prefactor != "":
            numerator, denominator = prefactor[1:-1].split('/')
            prefactor = f"\\left(\\frac{{{numerator}}}{{{denominator}}}\\right)"
            # prefactor = f"\\frac{{{numerator}}}{{{denominator}}}"

        # no W operator
        if (h.m == order and h.n == 0) or h.m > order:
            w_string = ""
            term_list.append(prefactor + h_string + w_string)
            continue

        # 1 W operator
        elif len(w_list) == 1 and bold_w_latex in w_list[0]:

            power = _extract_power(w_list[0])

            if expand_order_n_w is None or power < expand_order_n_w:
                w_string = bold_w_latex + construct_upper_w_label(h, order)
                term_list.append(prefactor + h_string + w_string)
            else:
                t_string = bold_t_latex + construct_upper_w_label(h, order)
                term_list.append(prefactor + h_string + t_string)

            continue

        # many W operators
        else:
            w_string = ""
            for w_op in w_list:

                # these terms only appear because we are mixing in the VECC through the e^(t^1) operator
                if bold_c_latex in w_op:
                    # the term is c^x
                    if bold_t_latex not in w_op:
                        # add 1/n! to prefactor
                        power = _extract_power(w_op)
                        # we don't need to update prefactor if we keep the c term
                        # updated_prefactor = f"\\left(\\dfrac{{{numerator}}}{{{denominator}*{power}!}}\\right)"
                        term_list.append(prefactor + h_string + _construct_c_term(h, power, order))
                        continue

                    # the term is c^x * t^y
                    else:
                        # c_power, t_power = _extract_ct_powers(w_op)
                        # term_list.append(prefactor + h_string + "\\left[" + w_op + "\\right]")

                        term_list.append(h_string + "\\left[" + w_op + "\\right]")  # no prefactor for the moment

                        # we will develop the expansion code here at a later date

                        # # we have to modify the prefactor
                        # if c_power > 1:
                        #     # add 1/n! to prefactor
                        #     updated_prefactor = f"\\left(\\dfrac{{{numerator}}}{{{denominator}*{c_power}!}}\\right)"
                        #     term_list.append(updated_prefactor + h_string + _construct_c_and_t_term(h, c_power, t_power, order))
                        # # normal prefactor
                        # else:
                        #     term_list.append(prefactor + h_string + _construct_c_and_t_term(h, c_power, t_power, order))

                        continue

                # this is the condensed linked-disconnected term contribution
                elif bold_d_latex in w_op:
                    d_string = bold_d_latex + construct_upper_w_label(h, order)
                    term_list.append(prefactor + h_string + d_string)
                    continue


                # this is the VECI contribution (simple t amplitude t^y)
                elif bold_t_latex in w_op and w_op.count(bold_t_latex) == 1:
                    t_string = bold_t_latex + construct_upper_w_label(h, order)
                    term_list.append(prefactor + h_string + t_string)
                    continue
                else:
                    raise Exception()

    return


def _generate_mixedcc_condensed_latex_term(term_list, h_list, w_list):
    """Append a string representation of a summation term for a residual to `term_list`."""

    if len(h_list) == 1:
        h_string = h_list[0].name.replace('h', bold_h_latex)
    else:
        bold_h_list = [h.name.replace('h', bold_h_latex) for h in h_list]
        h_string = f"({' + '.join(bold_h_list)})"

    if len(w_list) == 1:
        w_string = w_list[0]
    else:
        w_string = f"({' + '.join(w_list)})"

    # save the string representation of the term
    term_list.append(h_string + w_string)
    return


def _generate_mixedcc_latex_form(hamiltonian, order, max_order, expand_order_n_w=3, condense_disconnected_terms=True, explicit=False):
    """Generate the VECI/CC latex for the residual equations.
    Default is to generate the condensed form, but if `explicit` flag is `True` then
    explicit form with all i and k terms and prefactors is returned.
    """
    w_range = list(range(0, order + hamiltonian.maximum_rank + 1))

    # generate a list of H terms of equal or lower order than the residual
    operators = []
    for h in hamiltonian.operator_list:
        if h.m <= order:
            operators.append(h)


    term_list = [] # store all terms in here


    # generate each term we need to sum over
    for n in w_range:

        h_list, w_list = [], []

        # fill up the `w_list`
        if expand_order_n_w is None or n < expand_order_n_w:
            # create the plain w string
            w_list.append(f"{bold_w_latex}^{{{n}}}" if n != 0 else "")

        elif n >= expand_order_n_w:
            # connected term
            w_list.append(f"{bold_t_latex}^{{{n}}}")

            # linked disconnected terms
            if condense_disconnected_terms:
                w_list.append(f"{bold_d_latex}^{{{n}}}")
            else:
                for partition in sorted(generate_linked_disconnected_partitions_of_n(n)):
                    if max(partition) == 1:
                        w_list.append(f"{bold_c_latex}^{{{len(partition)}}}")
                    else:
                        w_list.append(f"{bold_c_latex}^{{{len(partition)-1}}}" + f"{bold_t_latex}^{{{max(partition)}}}")


        # collect h terms that map to the w^n operators in `w_list`
        for i, h in enumerate(operators):
            if abs(h.m - h.n - order) == n:
                h_list.append(h)
                operators.pop(i)

        # if no terms we skip this W operator
        if len(h_list) == 0:
            continue

        # append all term strings to `term_list`
        if explicit:
            _generate_mixedcc_explicit_latex_term(term_list, h_list, w_list, order, expand_order_n_w)
        else:
            _generate_mixedcc_condensed_latex_term(term_list, h_list, w_list)

        # loop


    # return the full string (with prefix and all term's wrapped in plus symbols)
    if explicit:
        i_terms = ", ".join([f"i_{n}" for n in range(1, order+1)]) if order != 0 else "0"
        return_string = f"\\begin{{split}} \\textbf{{R}}_{{{i_terms}}} &=\n" + " +\n".join(term_list) + " \\end{split}"
        return return_string
    else:
        return_string = f"\\textbf{{R}}_{{{order}}} &= " + " + ".join(term_list)
        return return_string


def _generate_mixedcc_wave_operator_2(hamiltonian, order, max_order):
    """Generate the latex for the VECI/CC mixed wave operator w^{`order`}."""
    return_string = ""

    if order == 0:
        return f"{bold_w_latex}^{{0}} &= \\textbf{{1}}"
    elif order == 1:
        return f"{bold_w_latex}^{{1}} &= {bold_t_latex}^{{1}}"
    else:
        c_list, d_list = [], []

        c_list.append(f"{bold_c_latex}^{{{order}}}")
        for i in range(2, order):
            c_list.append(f"{bold_c_latex}^{{{order-i}}}" + f"{bold_t_latex}^{{{i}}}")
        c_list.append(f"{bold_t_latex}^{{{order}}}")

        d_list.append(f"\\textbf{{d}}^{{{order}}}")
        d_list.append(f"{bold_t_latex}^{{{order}}}")

        return_string = f"{bold_w_latex}^{{{order}}} &= " + " + ".join(c_list) + " &&=" + " + ".join(d_list)
    return return_string


def _generate_mixedcc_wave_operator_1(hamiltonian, max_order):
    """Generate the latex for all the VECI/CC mixed wave operators up to w^{`max_order`}."""
    line1, line2 = "", ""
    for order in range(max_order + 1):
        if order == 0:
            # line1 += f"{bold_w_latex}^{{0}} &= \\textbf{{1}} &"
            # line2 += f"{bold_w_latex}^{{0}} &= \\textbf{{1}} &"
            continue
        elif order == 1:
            # line1 += f"{bold_w_latex}^{{1}} &= {bold_t_latex}^{{1}} &"
            # line2 += f"{bold_w_latex}^{{1}} &= {bold_t_latex}^{{1}} &"
            continue
        else:
            # c_list is the uncompressed format
            # d_list is the compressed format
            c_list, d_list = [], []
            c_list, d_list = [], []

            # connected term
            c_list.append(f"{bold_t_latex}^{{{order}}}")
            d_list.append(f"{bold_t_latex}^{{{order}}}")

            # linked disconnected terms
            d_list.append(f"\\textbf{{d}}^{{{order}}}")
            for partition in sorted(generate_linked_disconnected_partitions_of_n(order)):
                if max(partition) == 1:
                    c_list.append(f"{bold_c_latex}^{{{len(partition)}}}")
                else:
                    c_list.append(f"{bold_c_latex}^{{{len(partition)-1}}}" + f"{bold_t_latex}^{{{max(partition)}}}")

            # record lines
            line1 += f"{bold_w_latex}^{{{order}}} &= " + " + ".join(c_list) + " & "
            line2 += f"{bold_w_latex}^{{{order}}} &= " + " + ".join(d_list) + " & "

    return line1 + ' \\\\\n' + line2


def _generate_mixedcc_latex(hamiltonian, max_order):
    """Generate all of the VECI/CC mixed latex equations.
    We use `join` to insert two backward's slashes \\ BETWEEN each line
    rather then adding them to end and having extra trailing slashes on the last line.
    The user is expected to manually copy the relevant lines from the text file into a latex file
    and generate the pdf themselves.
    """
    expand_order_n_w = 3  # expand all terms of order n or higher
    return_string = ""  # store it all in here


    return_string += "VECI/CC wave operator latex\n"
    if True:
        return_string += _generate_mixedcc_wave_operator_1(hamiltonian, max_order)
    else:
        return_string += ' \\\\\n'.join(
            [_generate_mixedcc_wave_operator_2(hamiltonian, order, max_order) for order in range(max_order+1)]
        )

    return_string += "\n"*4
    return_string += "VECI/CC mixed condensed latex\n"
    return_string += ' \\\\\n%\n'.join(
        [_generate_mixedcc_latex_form(hamiltonian, order, max_order, expand_order_n_w) for order in range(max_order+1)]
    )

    return_string += "\n"*4
    return_string += "VECI/CC mixed explicit latex\n"
    return_string += ' \\\\\n%\n'.join(
        [_generate_mixedcc_latex_form(hamiltonian, order, max_order, expand_order_n_w, explicit=True) for order in range(max_order+1)]
    )
    return return_string


# ------------------- generating VECC residual latex --------------------- #

def _generate_vecc_explicit_latex_term(term_list, h_list, w_list, order, expand_order_n_w=None):
    """Generate the latex for the residual equations including all the indices and factors."""

    for h in h_list:
        h_string = bold_h_latex
        h_string += construct_upper_h_label(h, order)
        h_string += construct_lower_h_label(h, order)


        prefactor = construct_prefactor(h,  order, True)
        if prefactor != "":
            numerator, denominator = prefactor[1:-1].split('/')
            prefactor = f"\\left(\\frac{{{numerator}}}{{{denominator}}}\\right)"
            # prefactor = f"\\frac{{{numerator}}}{{{denominator}}}"

        # no W operator
        if (h.m == order and h.n == 0) or h.m > order:
            w_string = ""
            term_list.append(prefactor + h_string + w_string)
            continue

        # 1 W operator
        elif len(w_list) == 1 and bold_w_latex in w_list[0]:

            power = _extract_power(w_list[0])

            if expand_order_n_w is None or power < expand_order_n_w:
                w_string = bold_w_latex + construct_upper_w_label(h, order)
                term_list.append(prefactor + h_string + w_string)
            else:
                t_string = bold_t_latex + construct_upper_w_label(h, order)
                term_list.append(prefactor + h_string + t_string)

            continue

        # many W operators
        else:
            w_string = ""
            for w_op in w_list:

                # these terms only appear because we are mixing in the VECC through the e^(t^1) operator
                if bold_c_latex in w_op:
                    # the term is c^x
                    if bold_t_latex not in w_op:
                        # add 1/n! to prefactor
                        power = _extract_power(w_op)
                        # we don't need to update prefactor if we keep the c term
                        # updated_prefactor = f"\\left(\\dfrac{{{numerator}}}{{{denominator}*{power}!}}\\right)"
                        term_list.append(prefactor + h_string + _construct_c_term(h, power, order))
                        continue

                    # the term is c^x * t^y
                    else:
                        # c_power, t_power = _extract_ct_powers(w_op)
                        # term_list.append(prefactor + h_string + "\\left[" + w_op + "\\right]")

                        term_list.append(h_string + "\\left[" + w_op + "\\right]")  # no prefactor for the moment

                        # we will develop the expansion code here at a later date

                        # # we have to modify the prefactor
                        # if c_power > 1:
                        #     # add 1/n! to prefactor
                        #     updated_prefactor = f"\\left(\\dfrac{{{numerator}}}{{{denominator}*{c_power}!}}\\right)"
                        #     term_list.append(updated_prefactor + h_string + _construct_c_and_t_term(h, c_power, t_power, order))
                        # # normal prefactor
                        # else:
                        #     term_list.append(prefactor + h_string + _construct_c_and_t_term(h, c_power, t_power, order))

                        continue

                # this is the condensed linked-disconnected term contribution
                elif bold_d_latex in w_op:
                    d_string = bold_d_latex + construct_upper_w_label(h, order)
                    term_list.append(prefactor + h_string + d_string)
                    continue

                # this is the VECI contribution (simple t amplitude t^y)
                elif bold_t_latex in w_op and w_op.count(bold_t_latex) == 1:
                    t_string = bold_t_latex + construct_upper_w_label(h, order)
                    term_list.append(prefactor + h_string + t_string)
                    continue

                # this is the unlinked-disconnected term (the VECC contributions)
                elif bold_t_latex in w_op and w_op.count(bold_t_latex) > 1:
                    # string = bold_t_latex + construct_upper_w_label(h, order)
                    # term_list.append(prefactor + h_string + string)
                    term_list.append(h_string + "\\left[" + w_op + "\\right]")  # no prefactor for the moment
                    continue

                else:
                    raise Exception()

    return


def _generate_vecc_condensed_latex_term(term_list, h_list, w_list):
    """Append a string representation of a summation term for a residual to `term_list`."""

    if len(h_list) == 1:
        h_string = h_list[0].name.replace('h', bold_h_latex)
    else:
        bold_h_list = [h.name.replace('h', bold_h_latex) for h in h_list]
        h_string = f"({' + '.join(bold_h_list)})"

    if len(w_list) == 1:
        w_string = w_list[0]
    else:
        w_string = f"({' + '.join(w_list)})"

    # save the string representation of the term
    term_list.append(h_string + w_string)
    return


def _generate_vecc_latex_form(hamiltonian, order, max_order, expand_order_n_w=3, condense_disconnected_terms=True, explicit=False):
    """Generate the VECC latex for the residual equations.
    Default is to generate the condensed form, but if `explicit` flag is `True` then
    explicit form with all i and k terms and prefactors is returned.
    """
    w_range = list(range(0, order + hamiltonian.maximum_rank + 1))

    # generate a list of H terms of equal or lower order than the residual
    operators = []
    for h in hamiltonian.operator_list:
        if h.m <= order:
            operators.append(h)

    term_list = []  # store all terms in here


    # generate each term we need to sum over
    for n in w_range:

        h_list, w_list = [], []

        # fill up the `w_list`
        if expand_order_n_w is None or n < expand_order_n_w:
            # create the plain w string
            w_list.append(f"{bold_w_latex}^{{{n}}}" if n != 0 else "")

        elif n >= expand_order_n_w:

            # connected term
            w_list.append(f"{bold_t_latex}^{{{n}}}")

            # linked disconnected terms
            if condense_disconnected_terms:
                w_list.append(f"{bold_d_latex}^{{{n}}}")
            else:
                for partition in sorted(generate_linked_disconnected_partitions_of_n(n)):
                    if max(partition) == 1:
                        w_list.append(f"{bold_c_latex}^{{{len(partition)}}}")
                    else:
                        w_list.append(f"{bold_c_latex}^{{{len(partition)-1}}}" + f"{bold_t_latex}^{{{max(partition)}}}")

            # un-linked disconnected terms
            for partition in sorted(generate_un_linked_disconnected_partitions_of_n(n)):
                unlinked_disconnected_term = "".join([f"{bold_t_latex}^{{{term}}}" for term in partition])
                w_list.append(unlinked_disconnected_term)

        # collect h terms that map to the w^n operators in `w_list`
        for i, h in enumerate(operators):
            if abs(h.m - h.n - order) == n:
                h_list.append(h)
                operators.pop(i)

        # if no terms we skip this W operator
        if len(h_list) == 0:
            continue

        # append all term strings to `term_list`
        if explicit:
            _generate_vecc_explicit_latex_term(term_list, h_list, w_list, order, expand_order_n_w)
        else:
            _generate_vecc_condensed_latex_term(term_list, h_list, w_list)

    # return the full string (with prefix and all term's wrapped in plus symbols)
    if explicit:
        i_terms = ", ".join([f"i_{n}" for n in range(1, order+1)]) if order != 0 else "0"
        return_string = f"\\begin{{split}} \\textbf{{R}}_{{{i_terms}}} &=\n" + " +\n".join(term_list) + " \\end{split}"
        return return_string
    else:
        # higher order terms are very long and will probably need to be split over many lines
        if order >= 4:
            return_string = f"\\begin{{split}} \\textbf{{R}}_{{{order}}} &=\n" + " +\n".join(term_list) + " \\end{split}"
            return return_string
        else:
            return_string = f"\\textbf{{R}}_{{{order}}} &= " + " + ".join(term_list)
            return return_string

    return return_string


def _generate_vecc_wave_operator(hamiltonian, max_order):
    """Generate the latex for all the VECC wave operators up to w^{`max_order`}."""
    line1, line2 = "", ""
    for order in range(max_order + 1):
        if order == 0:
            # line1 += f"{bold_w_latex}^{{0}} &= \\textbf{{1}} &"
            # line2 += f"{bold_w_latex}^{{0}} &= \\textbf{{1}} &"
            continue
        elif order == 1:
            # line1 += f"{bold_w_latex}^{{1}} &= {bold_t_latex}^{{1}} &"
            # line2 += f"{bold_w_latex}^{{1}} &= {bold_t_latex}^{{1}} &"
            continue
        else:
            # c_list is the uncompressed format
            # d_list is the compressed format
            c_list, d_list = [], []

            # connected term
            c_list.append(f"{bold_t_latex}^{{{order}}}")
            d_list.append(f"{bold_t_latex}^{{{order}}}")

            # linked disconnected terms
            d_list.append(f"\\textbf{{d}}^{{{order}}}")
            for partition in sorted(generate_linked_disconnected_partitions_of_n(order)):
                if max(partition) == 1:
                    c_list.append(f"{bold_c_latex}^{{{len(partition)}}}")
                else:
                    c_list.append(f"{bold_c_latex}^{{{len(partition)-1}}}" + f"{bold_t_latex}^{{{max(partition)}}}")

            # un-linked disconnected terms
            for partition in sorted(generate_un_linked_disconnected_partitions_of_n(order)):
                string = ""
                for term in partition:
                    string += f"{bold_t_latex}^{{{term}}}"
                c_list.append(string)
                d_list.append(string)

            # record lines
            line1 += f"{bold_w_latex}^{{{order}}} &= " + " + ".join(c_list) + " & "
            line2 += f"{bold_w_latex}^{{{order}}} &= " + " + ".join(d_list) + " & "

    return line1 + ' \\\\\n' + line2


def _generate_vecc_latex(hamiltonian, max_order):
    """Generate all of the VECI/CC mixed latex equations.
    We use `join` to insert two backward's slashes \\ BETWEEN each line
    rather then adding them to end and having extra trailing slashes on the last line.
    The user is expected to manually copy the relevant lines from the text file into a latex file
    and generate the pdf themselves.
    """
    expand_order_n_w = 3  # expand all terms of order n or higher
    return_string = ""  # store it all in here


    return_string += "VECC wave operator latex\n"
    return_string += _generate_vecc_wave_operator(hamiltonian, max_order)

    return_string += "\n"*4
    return_string += "VECC condensed latex\n"
    return_string += ' \\\\\n%\n'.join(
        [_generate_vecc_latex_form(hamiltonian, order, max_order, expand_order_n_w) for order in range(max_order+1)]
    )

    return_string += "\n"*4
    return_string += "VECC explicit latex\n"
    return_string += ' \\\\\n%\n'.join(
        [_generate_vecc_latex_form(hamiltonian, order, max_order, expand_order_n_w, explicit=True) for order in range(max_order+1)]
    )
    return return_string


# ------------------------------------------------------------------------ #
def generate_residual_equations_latex(max_residual_order, path="./generated_latex.txt"):
    """Generates and saves to a file the code to calculate the residual equations for the CC approach."""

    # generate the Hamiltonian data
    H = generate_hamiltonian_operator()
    for h in H.operator_list:
        print(h)

    # write latex code
    string = _generate_veci_latex(H, max_residual_order)

    # write latex code
    string += "\n"*4 + _generate_mixedcc_latex(H, max_residual_order)

    # write latex code
    string += "\n"*4 + _generate_vecc_latex(H, max_residual_order)

    # save data
    with open(path, 'w') as fp:
        fp.write(string)

# ----------------------------------------------------------------------------------------------- #
# -----------------------------  GENERATING RESIDUAL EQUATIONS  --------------------------------- #
# ----------------------------------------------------------------------------------------------- #
h_dict = {
    (0,0): 'h_ab',
    (0,1): 'h_abI',
    (1,0): 'h_abi',
    (1,1): 'h_abIj',
    (0,2): 'h_abIJ',
    (2,0): 'h_abij',
}

w_dict = {
    1: 'w_i',
    2: 'w_ij',
    3: 'w_ijk',
    4: 'w_ijkl',
    5: 'w_ijklm',
    6: 'w_ijklmn',
    7: 'w_ijklmno',
    8: 'w_ijklmnop',
}


i_list = ['i', 'j', 'k', 'l']
k_list = ['m', 'n']


def _generate_einsum_h_indices(term):
    """Generate the indices for the h term in the residual's einsum equation
    upper or lower k's should map to the letters 'm' and 'n'
    upper or lower i's should map to the letters 'ijkl'
    """
    h_dims = 'ac'
    h_dims += ''.join(k_list[0:term.h.max_k])  # generate mn's
    h_dims += ''.join(i_list[0:term.h.max_i])  # generate jikl's
    return h_dims


def _generate_einsum_w_indices(term):
    """Generate the indices for the W operator in the residual's einsum equation
    upper or lower k's should map to the letters 'm' and 'n'
    upper or lower i's should map to the letters 'ijkl'
    """
    w_dims = 'cb'
    w_dims += ''.join(k_list[0:term.w.max_k])  # generate mn's
    w_dims += ''.join(i_list[term.h.max_i:(term.h.max_i + term.w.max_i)])  # generate jikl's

    return w_dims


def _generate_einsum_ouput_indices(term):
    """Generate the labels for the residual output."""
    # number of normal mode labels in output should be equal to the sum of the `max_i`'s
    return 'ab' + ''.join(i_list[0:(term.h.max_i + term.w.max_i)])


def _residual_terms_einsum(term, suppress_1_prefactor=True):
    """Returns a python code in `str` format which calculates the contribution to the residual from `term`."""

    # create the prefactor
    if suppress_1_prefactor and (term.prefactor == 1.0):  # if we don't want to print prefactors that are 1
        prefactor = ""
    else:
        prefactor = str(term.prefactor) + ' * '

    # create the string naming the h python object
    h = h_dict[(term.h.max_i, term.h.max_k)]

    # if there is no W operator then we simply return the h tensor multiplied by the prefactor
    if term.w.order == 0:
        return f"R += {prefactor}{h}\n"

    # create the string naming the w python object
    w = w_dict[term.w.order]

    # generate command
    h_dims = _generate_einsum_h_indices(term)
    w_dims = _generate_einsum_w_indices(term)
    out_dims = _generate_einsum_ouput_indices(term)

    return f"R += {prefactor}np.einsum('{h_dims},{w_dims}->{out_dims}', {h}, {w})\n"


def _same_w_order_term_list(current_term, term_list):
    """Returns a list of terms in the list whose W operator order is the same as the `current_term`'s"""
    return [term for term in term_list if term.w.order == current_term.w.order]


def write_residual_function_string(residual_terms_list, order):
    """ Output a string of python code to calculate a specific residual """

    string = "" # we store the output in this string

    # this line of python code labels the terms in the `w_args` tuple
    w_args_string = ", ".join([w_dict[n] for n in range(1, order+2)] + ["*unusedargs"])
    # print(f"{w_args_string=}")

    # the function definition and initialization of the residual array
    string += (
        f"\ndef calculate_order_{order}_residual(A, N, truncation, h_args, w_args):\n"
        f'{tab}"""Calculate the {order} order residual as a function of the W operators."""\n'
        f"{tab}h_ab, h_abI, h_abi, h_abIj, h_abIJ, h_abij = h_args\n"
        f"{tab}{w_args_string} = w_args\n"
        "\n"
        f"{tab}R = np.zeros(({', '.join(['A','A',] + ['N',]*order)}), dtype=complex)\n"
    )

    # add code to assert correct truncation
    if order <= 1:
        assertion_check = f"truncation.{taylor_series_order_tag[1]}"
    else:
        assertion_check = f"truncation.{taylor_series_order_tag[order]}"

    string += (
        "\n"
        f'{tab}assert {assertion_check}, \\\n'
        f'{tab}{tab}f"Cannot calculate order {order} residual for {{truncation.cc_truncation_order}}"\n'
    )

    # a list of all terms whose h operator is quadratic (2nd or higher order)
    quadratic_terms = [term for term in residual_terms_list if (term.h.max_i >= 2) or (term.h.max_k >= 2)]

    # each time we write a term to the `string` object append it to this list
    # so we don't write any terms twice
    already_printed_list = []

    for term in residual_terms_list:

        if term in already_printed_list:
            continue

        string += '\n' # space out the einsum commands

        if term.w.order == 0 or term.w.order == 1:
            already_printed_list.append(term)
            string += f"{tab}" + _residual_terms_einsum(term)

        # if the h operator is 2nd or higher order
        elif (term.h.max_i >= 2) or (term.h.max_k >= 2):

            string += f"{tab}if truncation.quadratic:\n"

            # find all the other terms whose h operator is quadratic (2nd or higher order)
            for quad_term in quadratic_terms:
                already_printed_list.append(quad_term)
                # if the next quadratic term has a W operator which is 2nd or higher order
                if bool(quad_term.w.order >= 2):
                    string += (
                        f"{tab}{tab}if {w_dict[quad_term.w.order]} is not None:\n"
                        f"{tab}{tab}{tab}" + _residual_terms_einsum(quad_term)
                    )
                else:
                    string += (
                        f"{tab}{tab}else:\n"
                        f"{tab}{tab}{tab}" + _residual_terms_einsum(quad_term)
                    )

        # if the W operator is 2nd or higher order
        elif bool(term.w.order >= 2):
            string +=f"{tab}if {w_dict[term.w.order]} is not None:\n"

            # find all the other terms whose W operator is the same order as `term`
            for same_w_term in _same_w_order_term_list(term, residual_terms_list):
                already_printed_list.append(same_w_term)
                string += f"{tab}{tab}" + _residual_terms_einsum(same_w_term)

        else:
            raise Exception("We shouldn't reach this else!")

    string += (f"\n{tab}return R\n")
    return string


def generate_python_code_for_residual_functions(term_lists, max_order):
    """Return a string containing the python code to generate residual functions up to (and including) `max_order`.
    `term_lists` is a list of lists, each of which contain tuples `(prefactor, h, w)` representing terms of that particular residual.
    Requires the following header: `"import numpy as np"`.
    """
    lst = [write_residual_function_string(term_lists[order], order=order) for order in range(max_order+1)]
    return "".join(lst)


def generate_residual_equations_file(max_residual_order, path="./residual_equations.py"):
    """Generates and saves to a file the code to calculate the residual equations for the CC approach."""

    # generate the Hamiltonian data
    H = generate_hamiltonian_operator()
    for h in H.operator_list:
        print(h)

    # generate the residual data
    R_lists, term_lists = generate_residual_data(H.operator_list, max_order=max_residual_order)

    # if we want to print stuff for debugging
    print_residual_data(R_lists, term_lists, print_equations=False, print_tuples=False)

    # start with the import statements
    file_data = (
        "# system imports\n"
        "\n"
        "# third party imports\n"
        "import numpy as np\n"
        "\n"
        "# local imports\n"
        "from .symmetrize import symmetrize_tensor\n"
        "\n"
    )

    # write the functions to calculate the residuals
    file_data += generate_python_code_for_residual_functions(term_lists, max_order=max_residual_order)

    # save data
    with open(path, 'w') as fp:
        fp.write(file_data)

    return

# ----------------------------------------------------------------------------------------------- #
# ---------------------------  GENERATING W OPERATOR EQUATIONS  --------------------------------- #
# ----------------------------------------------------------------------------------------------- #

def _next_list(lst, n):
    ''' if items in `lst` are all one, or there is only one 1 in `lst`, add 1 to the last item and delete the first 1 in the list like:
            [1, 1, 1] -> [1, 2] and [1, 2] -> [3];
        if there is no 1 in `lst`, return the list itself;
        if there more than one 1 in `lst`, delete one of 1 and add 1 to the last and the one before the last item separately, like:
            [1, 1, 2] -> [(1,3), (2, 2)]'''
    if lst.count(1) == 0:
        return [lst, ]
    elif lst.count(1) == n or lst.count(1) == 1:
        result = lst[1:-1] + [(lst[-1]+1), ]
        return [result, ]
    else:
        result = [(lst[1:-2] + [lst[-2]+1, ] + [lst[-1], ]), (lst[1:-1] + [(lst[-1]+1), ])]
        return result


def _generate_t_lists(n):
    ''' generates a list of lists in which each item is from 1 to n and the sum of items are n, for example:
        4 -> [[4], [1, 3], [2, 2], [1, 1, 2], [1, 1, 1, 1]]'''
    first = [1]*n
    result = []
    l = [first]
    current = _next_list(l[-1], n)

    while current != [l[-1]]: # if there is no 1 in current(a list), l gets all items, bit not in correct order
        l += current
        current = _next_list(l[-1], n)

    # sort l to get the final result
    previous_max = max(l[0]) # initialize previous_max to record the maximum value in the previous list

    for sub_l in l:
        current_max = max(sub_l) # record the maximum in the first list in l

        # if the maximum item in sub-l is larger than the max in previous sub-list, add sub_l to the end of result list
        if current_max >= previous_max:
            result.append(sub_l)
            previous_max = current_max
        # if the max in sub-l is smaller than the max in previous sub-list, add sub_l to the beginning of result list
        else:
            result = [sub_l] + result
    result.reverse()
    return result


def _generate_t_terms_dictionary(n):
    ''' generates a dictionary in which keys are tuple that show the tag for t terms as number, and values are the actual t terms,
        for example: {(1, 4): ('t_i', 't_ijkl')}'''
    result = {}
    t_tag_list = ["i", "j", "k", "l", "m", "n"]
    t_num_list = _generate_t_lists(n)  # a list of lists like [1,1], [2,1,1]....
    for num_term in t_num_list:
        t_term = [] # stores t_i, t_ij, t_ijk....
        for n in num_term:
            t_str = "t_" + "".join(t_tag_list[:n]) # add i,j,k.....to "t_"
            # print(n, t_str, t_terms[n-1].string)
            t_term.append(t_str)
        result[tuple(num_term)] = tuple(t_term)
    return result


def _permutation_function(l_t, l_fix, l_perm, n):
    ''' generates a list of character combinations which will be used in the einsum calculation,
        l_t contains combination of t terms group like [["t_i", "t_ij"], ["t_ij", "t_i"]];
        l_fix contains the fixed part of character combinations which will be used in the einsum calculation like
            ["ac", "cd", "db"] for (t_i, t_i, t_ij) group;
        l_perm contains groups of characters which will be added to the fixed part later, for example:
            [[i, j], [j, i]] for [t_i, t_i] and ["ac","cb"]
        '''

    result = []
    for char_group in l_perm:
        s_result = []
        for t_group in l_t:
            ss_result = [list(l_fix)]
            i = 0
            p_sum = 0
            for t_item in t_group:
                length = len(t_item) - 2 # delete the length of "t_" part
                ss_result[0][i] += "".join(char_group[p_sum:length+p_sum])
                p_sum += length
                i += 1
            if len(l_fix) != n:
                s_result += ss_result
        if len(l_fix) == n:
            s_result += ss_result
        result += [s_result]

    #print(result)
    return result


def _write_permutations(perm_t, perm_char_list, W_array, prefactor):
    ''' generates lines for w function if permutations inside the einsum calculation is needed'''
    result = ""
    for p_group in perm_char_list:
        # list of command strings for einsum
        # cmd = ",".join(einsum_perm) + " -> " + einsum_result
        cmd_list = [", ".join(p_c) for p_c in p_group]
        # list of arguments to einsum call
        arg_list = [", ".join(perm) for perm in perm_t]
        # compile a list of the einsum calls (as strings)
        einsum_list = [f"np.einsum('{cmd_list[i]}', {arg_list[i]})" for i in range(len(p_group))]
        # join them together
        sum_of_einsums = " + ".join(einsum_list)
        # add the whole line to the main string
        result += f"{tab}{W_array} += {prefactor} * ({sum_of_einsums})\n"
    return result


def old_write_w_function_strings_show_all(order):
    """ x """

    string = "" # we store the output in this string

    num_tag = ["zero", "1st", "2nd", "3rd", "4th", "5th"]
    word_tag = ["", "Singles", "Doubles", "Triples", "Quadruples", "Quintuples"]
    einsum_surface_tags = "acdef"
    tag_str = "ijklmn"
    W_array = f"W_{order}"

    string += (
        f"\ndef calculate_order_{order}_w_operator(A, N, t_args):\n"
        f'{tab}"""Calculate the {order} order W operator for use in the calculation of the residuals."""\n'
        f"{tab}t_i, t_ij, t_ijk, t_ijkl = t_args\n"
        f"{tab}# Creating the {num_tag[order]} order W operator\n"
        f"{tab}{W_array} = np.zeros(({', '.join(['A','A',] + ['N',]*order)}), dtype=complex)\n"
    )

    # simple cases, zero and 1st order
    if order == 0:
        return f"{tab}return {W_array}\n"
    if order == 1:
        string += (
            f"{tab}# Singles contribution\n"
            f"{tab}{W_array} += t_i\n"
            f"{tab}return {W_array}\n"
        )
        return string


    t_term_dict = _generate_t_terms_dictionary(order)


    previous_max = order
    for num_group in t_term_dict.keys():

        current_max = max(num_group)
        w_prefactor = generate_w_prefactors(list(num_group))

        # Add comment
        string += f"{tab}# {word_tag[current_max]} contribution\n"

        if len(num_group) == 1: # no permutation is needed for this term
            string += f"{tab}{W_array} += {w_prefactor} * {t_term_dict[num_group][0]}\n"
            continue

        einsum_fixed = [] # contains [ac, cd, de.....]
        einsum_perm = []

        t_perm = it.permutations(t_term_dict[num_group]) # permutations of t_i, t_ij, t_ijk....

        # generate the list of tensor indices that we don't iterate over
        i = 0
        previous_sum = 0
        for n in num_group:
            if i == len(num_group)-1:
                einsum_fixed.append(einsum_surface_tags[i]+"b")
                print(f"{einsum_fixed=}")
            else:
                einsum_fixed.append(einsum_surface_tags[i:i+2])
                print(f"{einsum_fixed=}")
            i += 1

        if current_max < previous_max: # check for the new type of contribution
            previous_max = current_max
        elif current_max == previous_max: # if is the same type of contribution
            string += "\n"

        tag_perm = it.permutations(tag_str[:order]) # permutation of (i,j,k), (i,k,j).....
        t_list = list(t_perm)
        e_combination = _permutation_function(t_list, einsum_fixed, list(tag_perm),order) # a list
        string += _write_permutations(t_list, e_combination, W_array, w_prefactor)


    string += f"{tab}return {W_array}\n"

    return string

# ------------------------------------------------------- #

t_terms = [
    None,
    t_term_namedtuple("t_i", 1, "(A, A, N)"),
    t_term_namedtuple("t_ij", 2, "(A, A, N, N)"),
    t_term_namedtuple("t_ijk", 3, "(A, A, N, N, N)"),
    t_term_namedtuple("t_ijkl", 4, "(A, A, N, N, N, N)"),
    t_term_namedtuple("t_ijklm", 5, "(A, A, N, N, N, N, N)"),
    t_term_namedtuple("t_ijklmn", 6, "(A, A, N, N, N, N, N, N)"),
    t_term_namedtuple("t_ijklmno", 7, "(A, A, N, N, N, N, N, N, N)"),
    t_term_namedtuple("t_ijklmnop", 8, "(A, A, N, N, N, N, N, N, N, N)"),
]


def _generate_w_operator_prefactor(tupl):
    """Generates the prefactor for each part of W terms.
    The theory from which the prefactors arise goes as follows:
        - The Taylor series contributes a `1/factorial(length(x))`.
        - Each integer `n` in the tuple contributes `1/factorial(n)`.
    We choose not to print out the 1/factorial(1) prefactors.
    """

    # if `tupl` is (1,)
    if max(tupl) == 1 and len(tupl) == 1:
        return ""

    # if all items in `tupl` are 1, such as (1, 1, 1, ...)
    elif max(tupl) == 1:
        return f"1/factorial({len(tupl)})"

    # if there is only one item in `tupl` and it's not 1
    elif len(tupl) == 1:
        return f"1/factorial({max(tupl)})"

    # otherwise `tupl` has 2 or more terms that are not all 1
    # we need to include a factorial of the length
    # and a factorial for each term that is greater than 1
    else:
        denominator_list = [f"factorial({len(tupl)})", ]
        for n in tupl:
            if n != 1:
                denominator_list.append(f"factorial({n})")
        denominator = " * ".join(denominator_list)
        return f"1/({denominator})"


num_tag = ["zero", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"]
taylor_series_order_tag = ["", "singles", "doubles", "triples", "quadruples", "quintuples", "sextuples"]
einsum_surface_tags = "acdefghi"
tag_str = "ijklmnop"


# ------------------------------------------------------- #
def _generate_surface_index(partition):
    """Return a list of strings of length 2 [`ac`, `cd`, `db`].
    Where the first char in the list should be `a` and the last char should be `b`
    and the rest of the letters in between are in alphabetical order.
    """
    max_i = len(partition)-1
    assert max_i <= 7, "Our current `einsum_surface_tags` can't support 7th order W operators"
    return_list = [einsum_surface_tags[i:i+2] for i in range(max_i)]
    return_list.append(einsum_surface_tags[max_i] + "b")
    return return_list


def _generate_mode_index(partition, order):
    """Return a list of strings  [`ij`, `k`, `l`] representing the mode indices
    for the einsum of the W operator.
    """
    assert order <= 6, "Our current `tag_str` can't support 7th order W operators"
    combinations = unique_permutations(partition)
    # log.debug(f"{combinations=}")

    return_list = []
    for comb in combinations:
        N_list = list(tag_str[:order])
        comb_list = [''.join([N_list.pop(0) for _ in range(n) ]) for n in comb]
        return_list.append(comb_list)
    return return_list


def _w_einsum_list(partition, order):
    """Returns a list of strings. Each string represents a call to np.einsum.
    The list is filled relative to the specific `partition` that is being calculated.
    """

    # the unique permutations of the `partition` of the integer `order`
    combinations = unique_permutations(partition)
    # the einsum indices for the surface dimensions
    surf_index = _generate_surface_index(partition)
    # the einsum indices for the normal mode dimensions
    mode_index = _generate_mode_index(partition, order)

    return_list = []

    # `combinations` is a list of tuples such as (2, 2, 1, ) or (5,)
    for i, tupl in enumerate(combinations):
        # the input dimensions are two character strings representing the surface dimensions
        # plus 1 or more characters representing the normal mode dimensions
        in_dims = ", ".join([surf_index[a]+mode_index[i][a] for a in range(len(surf_index)) ])
        # the output dimension is the same for all einsum calls for a given `partition` argument
        out_dims = f"ab{tag_str[0:order]}"
        # the names of the arguments to the einsum call are stored in the list `t_terms`
        # and the objects are accessed using each integer in the tuple (2, 1) -> ("t_ij", "t_i")
        pterms = ", ".join([t_terms[n].string for n in tupl])
        # print(f"np.einsum('{in_dims}->{out_dims}', {pterms})")
        return_list.append(f"np.einsum('{in_dims}->{out_dims}', {pterms})")

    return return_list


def _optimized_w_einsum_list(partition, order, iterator_name='optimized_einsum'):
    """Returns a list of strings. Each string represents a call to np.einsum.
    The list is filled relative to the specific `partition` that is being calculated.

    This function and `_construct_w_function_definition` need to agree on the name of the
    iterator containing the optimized paths. At the moment it is named `optimized_einsum`.
    """

    # the unique permutations of the `partition` of the integer `order`
    combinations = unique_permutations(partition)
    # the einsum indices for the surface dimensions
    surf_index = _generate_surface_index(partition)
    # the einsum indices for the normal mode dimensions
    mode_index = _generate_mode_index(partition, order)

    return_list = []

    # `combinations` is a list of tuples such as (2, 2, 1, ) or (5,)
    for i, tupl in enumerate(combinations):
        # the names of the arguments to the einsum call are stored in the list `t_terms`
        # and the objects are accessed using each integer in the tuple (2, 1) -> ("t_ij", "t_i")
        pterms = ", ".join([t_terms[n].string for n in tupl])
        return_list.append(f"next({iterator_name})({pterms})")

    return return_list


# ------------------------------------------------------- #
def _construct_vemx_contributions_definition(return_string, order, opt_einsum=False, iterator_name='optimized_einsum'):
    """Return the string containing the python code to prepare the function definition and unpack the t_args.
    This function and `_optimized_w_einsum_list` need to agree on the name of the
    iterator containing the optimized paths. At the moment it is named `optimized_einsum`.
    """

    # this line of python code labels the terms in the `t_args` tuple
    t_arg_string = ", ".join([t_terms[n].string for n in range(1, order)] + ["*unusedargs"])

    W_array = f"W_{order}"  # name of the W operator

    if not opt_einsum:
        return_string += (
            f"\ndef _add_order_{order}_vemx_contributions({W_array}, t_args, truncation):\n"
            f'{tab}"""Calculate the order {order} VECI/CC (mixed) contributions to the W operator\n'
            f'{tab}for use in the calculation of the residuals.\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args`\n"
            f"{tab}{t_arg_string} = t_args\n"
        )
    else:
        return_string += (
            f"\ndef _add_order_{order}_vemx_contributions_optimized({W_array}, t_args, truncation, opt_path_list):\n"
            f'{tab}"""Calculate the order {order} VECI/CC (mixed) contributions to the W operator\n'
            f'{tab}for use in the calculation of the residuals.\n'
            f'{tab}Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args`\n"
            f"{tab}{t_arg_string} = t_args\n"
        )
        if order >= 2:
            return_string += (
                f"{tab}# make an iterable out of the `opt_path_list`\n"
                f"{tab}{iterator_name} = iter(opt_path_list)\n"
            )

    return return_string


def _generate_vemx_contributions(order, opt_einsum=False):
    """ x """
    assert order <= 6, "Can only handle up to 6th order due to `einsum_surface_tags` limit"

    return_string = "" # we store the output in this string

    # prepare the function definition, unpacking of arguments, and initialization of the W array
    return_string += _construct_vemx_contributions_definition(return_string, order, opt_einsum)

    # special exception for order less than 2
    if order < 2:
        lst = return_string.split('"""')
        exception = ( "\n"
            f"{tab}raise Exception(\n"
            f'{tab}{tab}"the first possible purely VECI/CC term is (1, 1) or (dt_i * t_i)"\n'
            f'{tab}{tab}"which requires a W operator of at least 2nd order"\n'
            f"{tab})\n"
        )
        return '\"\"\"'.join([lst[0], "Exists for error checking.", exception])

    W_array = f"W_{order}"  # name of the W operator


    # for each partition (the mathematical term) of the integer `order`
    for partition in generate_linked_disconnected_partitions_of_n(order):

        contribution_name = taylor_series_order_tag[max(partition)].upper()
        contribution_string = ""

        # Label the order of this partition's contribution
        return_string += f"{tab}# {contribution_name} contribution\n"

        # if this partition is contributing doubles or higher then we
        # add the if statement to check for correct truncation level
        # note that we will add `{tab}`s to every line in `contribution_string`
        # before adding it to return_string (to account for indentation of if statement)
        if max(partition) >= 2:
            return_string += f"{tab}if truncation.{contribution_name.lower()}:\n"

        # make the prefactor
        prefactor = _generate_w_operator_prefactor(partition)

        # order 1 case is simple
        if len(partition) == 1: # no permutation is needed for this term
            print(partition, partition[0])
            # we have to space the line correct (how many tabs)
            if max(partition) >= 2:
                return_string += f"{tab}{tab}{W_array} += {prefactor} * {t_terms[partition[0]].string}\n"
            else:
                return_string += f"{tab}{W_array} += {prefactor} * {t_terms[partition[0]].string}\n"
            # continue onto the next partition
            continue

        # compile a list of the einsum calls
        einsum_list = _optimized_w_einsum_list(partition, order) if opt_einsum else _w_einsum_list(partition, order)

        if len(partition) == order:
            # join them together
            sum_of_einsums = " + ".join(einsum_list)
            # add the whole line to the main string
            contribution_string += f"{tab}{W_array} += {prefactor} * ({sum_of_einsums})\n"
        else:
            # join them together
            sum_of_einsums = f" +\n{tab}{tab}".join(einsum_list)
            # add the whole line to the main string
            contribution_string += (
                f"{tab}{W_array} += {prefactor} * (\n"
                f"{tab}{tab}{sum_of_einsums}\n"
                f"{tab})\n"
            )

        # we need to add a `{tab}` to every line in the `contribution_string`
        # to space the lines correctly to account for the if statement
        if max(partition) >= 2:
            lines = contribution_string.splitlines()
            contribution_string = "".join([f"{tab}{line}\n" for line in lines])
        # add the contribution to the `return_string`
        return_string += contribution_string
    # end of loop
    return_string +=  f"{tab}return\n"
    return return_string


# ------------------------------------------------------- #
def _construct_vecc_contributions_definition(return_string, order, opt_einsum=False, iterator_name='optimized_einsum'):
    """Return the string containing the python code to prepare the function definition and unpack the t_args.
    This function and `_optimized_w_einsum_list` need to agree on the name of the
    iterator containing the optimized paths. At the moment it is named `optimized_einsum`.
    """

    # this line of python code labels the terms in the `t_args` tuple
    t_arg_string = ", ".join([t_terms[n].string for n in range(1, order-1)] + ["*unusedargs"])

    W_array = f"W_{order}"  # name of the W operator

    if not opt_einsum:
        return_string += (
            f"\ndef _add_order_{order}_vecc_contributions({W_array}, t_args, truncation):\n"
            f'{tab}"""Calculate the order {order} VECC contributions to the W operator\n'
            f'{tab}for use in the calculation of the residuals.\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args`\n"
            f"{tab}{t_arg_string} = t_args\n"
        )
    else:
        return_string += (
            f"\ndef _add_order_{order}_vecc_contributions_optimized({W_array}, t_args, truncation, opt_path_list):\n"
            f'{tab}"""Calculate the order {order} VECC contributions to the W operator\n'
            f'{tab}"for use in the calculation of the residuals.\n'
            f'{tab}Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args`\n"
            f"{tab}{t_arg_string} = t_args\n"
        )
        if order >= 2:
            return_string += (
                f"{tab}# make an iterable out of the `opt_path_list`\n"
                f"{tab}{iterator_name} = iter(opt_path_list)\n"
            )

    return return_string


def _generate_vecc_contributions(order, opt_einsum=False):
    """ x """
    assert order <= 6, "Can only handle up to 6th order due to `einsum_surface_tags` limit"

    return_string = "" # we store the output in this string

    # prepare the function definition, unpacking of arguments, and initialization of the W array
    return_string += _construct_vecc_contributions_definition(return_string, order, opt_einsum)

    # special exception for order less than 2
    if order < 4:
        lst = return_string.split('"""')
        exception = ( "\n"
            f"{tab}raise Exception(\n"
            f'{tab}{tab}"the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"\n'
            f'{tab}{tab}"which requires a W operator of at least 4th order"\n'
            f"{tab})\n"
        )
        return '\"\"\"'.join([lst[0], "Exists for error checking.", exception])

    W_array = f"W_{order}"  # name of the W operator

    # for each partition (the mathematical term) of the integer `order`
    for partition in generate_un_linked_disconnected_partitions_of_n(order):

        contribution_name = taylor_series_order_tag[max(partition)].upper()
        contribution_string = ""

        # Label the order of this partition's contribution
        return_string += f"{tab}# {contribution_name} contribution\n"

        # if this partition is contributing doubles or higher then we
        # add the if statement to check for correct truncation level
        # note that we will add `{tab}`s to every line in `contribution_string`
        # before adding it to return_string (to account for indentation of if statement)
        if max(partition) >= 2:
            return_string += f"{tab}if truncation.{contribution_name.lower()}:\n"

        # make the prefactor
        prefactor = _generate_w_operator_prefactor(partition)

        # order 1 case is simple
        if len(partition) == 1: # no permutation is needed for this term
            print(partition, partition[0])
            # we have to space the line correct (how many tabs)
            if max(partition) >= 2:
                return_string += f"{tab}{tab}{W_array} += {prefactor} * {t_terms[partition[0]].string}\n"
            else:
                return_string += f"{tab}{W_array} += {prefactor} * {t_terms[partition[0]].string}\n"
            # continue onto the next partition
            continue

        # compile a list of the einsum calls
        einsum_list = _optimized_w_einsum_list(partition, order) if opt_einsum else _w_einsum_list(partition, order)

        if len(partition) == order:
            # join them together
            sum_of_einsums = " + ".join(einsum_list)
            # add the whole line to the main string
            contribution_string += f"{tab}{W_array} += {prefactor} * ({sum_of_einsums})\n"
        else:
            # join them together
            sum_of_einsums = f" +\n{tab}{tab}".join(einsum_list)
            # add the whole line to the main string
            contribution_string += (
                f"{tab}{W_array} += {prefactor} * (\n"
                f"{tab}{tab}{sum_of_einsums}\n"
                f"{tab})\n"
            )

        # we need to add a `{tab}` to every line in the `contribution_string`
        # to space the lines correctly to account for the if statement
        if max(partition) >= 2:
            lines = contribution_string.splitlines()
            contribution_string = "".join([f"{tab}{line}\n" for line in lines])
        # add the contribution to the `return_string`
        return_string += contribution_string
    # end of loop
    return_string +=  f"{tab}return\n"
    return return_string


# ------------------------------------------------------- #
def _construct_w_function_definition(return_string, order, opt_einsum=False, iterator_name='optimized_einsum'):
    """Return the string containing the python code to prepare the function definition and unpack the t_args.
    This function and `_optimized_w_einsum_list` need to agree on the name of the
    iterator containing the optimized paths. At the moment it is named `optimized_einsum`.
    """

    # this line of python code labels the terms in the `t_args` tuple
    t_arg_string = ", ".join([t_terms[n].string for n in range(1, order+1)] + ["*unusedargs"])

    if not opt_einsum:
        return_string += (
            f"\ndef _calculate_order_{order}_w_operator(A, N, t_args, ansatz, truncation):\n"
            f'{tab}"""Calculate the order {order} W operator for use in the calculation of the residuals."""\n'
            f"{tab}# unpack the `t_args`\n"
            f"{tab}{t_arg_string} = t_args\n"
        )
    else:
        return_string += (
            f"\ndef _calculate_order_{order}_w_operator_optimized(A, N, t_args, ansatz, truncation, vemx_opt_path_list, vecc_opt_path_list):\n"
            f'{tab}"""Calculate the order {order} W operator for use in the calculation of the residuals.\n'
            f'{tab}Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args`\n"
            f"{tab}{t_arg_string} = t_args\n"
        )

    return return_string


def _write_w_function_strings(order, opt_einsum=False):
    """ x """
    assert order <= 6, "Can only handle up to 6th order due to `einsum_surface_tags` limit"

    return_string = "" # we store the output in this string

    # prepare the function definition, unpacking of arguments, and initialization of the W array
    return_string += _construct_w_function_definition(return_string, order, opt_einsum)

    W_array = f"W_{order}"  # name of the W operator

    # initialize the W array
    return_string += (
        f"{tab}# Creating the {num_tag[order]} order W operator\n"
        f"{tab}{W_array} = np.zeros(({', '.join(['A','A',] + ['N',]*order)}), dtype=complex)\n"
    )

    # special base cases, zero and 1st order
    if order == 0:
        raise Exception("We should not call `_write_w_function_strings` with order 0.")
    if order == 1:
        return_string += (
            f"{tab}# Singles contribution\n"
            f"{tab}{W_array} += t_i\n"
            f"{tab}return {W_array}\n"
        )
        # return, no more code is needed for order 1 case
        return return_string

    # ------------------------------------------------------- #
    return_string += "\n"  # spacing the text to make code easier to read
    # ------------------------------------------------------- #
    # generate the appropriate strings `vemx_operations` and `vecc_operations
    # which are used in the f-strings about 25 lines below

    opt_string = "_optimized" if opt_einsum else ""
    vemx_opt_path_string = ", vemx_opt_path_list" if opt_einsum else""
    vecc_opt_path_string = ", vecc_opt_path_list" if opt_einsum else""

    if order < 2:
        vemx_operations = f"{tab}{tab}pass  # no VECI/CC (mixed) contributions for order < 2\n"
        vecc_operations = f"{tab}{tab}pass  # no VECC contributions for order < 4\n"
    else:
        vemx_operations = (
            f"{tab}{tab}_add_order_{order}_vemx_contributions{opt_string}({W_array}, t_args, truncation{vemx_opt_path_string})\n"
        )
        if order < 4:
            vecc_operations = (
                f"{tab}{tab}_add_order_{order}_vemx_contributions{opt_string}({W_array}, t_args, truncation{vemx_opt_path_string})\n"
                f"{tab}{tab}pass  # no VECC contributions for order < 4\n"
            )
        else:
            vecc_operations = (
                f"{tab}{tab}_add_order_{order}_vemx_contributions{opt_string}({W_array}, t_args, truncation{vemx_opt_path_string})\n"
                f"{tab}{tab}_add_order_{order}_vecc_contributions{opt_string}({W_array}, t_args, truncation{vecc_opt_path_string})\n"
            )

    # ------------------------------------------------------- #
    veci_prefactor = _generate_w_operator_prefactor((order,))
    # finally, we add it all together
    return_string += (
        f"{tab}# add the VECI contribution\n"
        f"{tab}if truncation.{taylor_series_order_tag[order]}:\n"
        f"{tab}{tab}{W_array} += {veci_prefactor} * {t_terms[order].string}\n"
        # disconnected terms
        f"{tab}if ansatz.VE_MIXED:\n"
        f"{vemx_operations}"
        f"{tab}elif ansatz.VECC:\n"
        f"{vecc_operations}"
    )
    # ------------------------------------------------------- #
    return_string += "\n"  # spacing the text to make code easier to read
    # ------------------------------------------------------- #
    return_string += (
        f"{tab}# Symmetrize the W operator\n"
        f"{tab}symmetric_w = symmetrize_tensor(N, {W_array}, order={order})\n"
        f"{tab}return symmetric_w\n"
    )

    return return_string


# ------------------------------------------------------- #
def _write_master_w_compute_function(max_order, opt_einsum=False):
    """Write the wrapper function which `vibronic_hamiltonian.py` calls.
    This functions decides, up to what order of W operator is calculated,
    based on the truncation level inside the namedtuple `truncation`.
    """
    if not opt_einsum:
        function_string = "w_{0} = _calculate_order_{0}_w_operator(A, N, t_args, ansatz, truncation)"

        string = f'''
            def compute_w_operators(A, N, t_args, ansatz, truncation):
                """Compute a number of W operators depending on the level of truncation."""

                if not truncation.singles:
                    raise Exception(
                        "It appears that `singles` is not true, this cannot be.\\n"
                        "Something went terribly wrong!!!\\n\\n"
                        f"{{truncation}}\\n"
                    )

                {function_string.format(1)}
                {function_string.format(2)}
                {function_string.format(3)}

                if not truncation.doubles:
                    return w_1, w_2, w_3, None, None, None
                else:
                    {function_string.format(4)}

                if not truncation.triples:
                    return w_1, w_2, w_3, w_4, None, None
                else:
                    {function_string.format(5)}

                if not truncation.quadruples:
                    return w_1, w_2, w_3, w_4, w_5, None
                else:
                    {function_string.format(6)}

                if not truncation.quintuples:
                    return w_1, w_2, w_3, w_4, w_5, w_6
                else:
                    raise Exception(
                        "Attempting to calculate W^7 operator (quintuples)\\n"
                        "This is currently not implemented!!\\n"
                    )
        '''
    else:
        function_string = (
            "w_{0} = _calculate_order_{0}_w_operator_optimized("
            "A, N, t_args, ansatz, truncation, vemx_optimized_paths[{1}], vecc_optimized_paths[{1}]"
            ")"
        )

        string = f'''
            def compute_w_operators_optimized(A, N, t_args, ansatz, truncation, vemx_optimized_paths, vecc_optimized_paths):
                """Compute a number of W operators depending on the level of truncation."""

                if not truncation.singles:
                    raise Exception(
                        "It appears that `singles` is not true, this cannot be.\\n"
                        "Something went terribly wrong!!!\\n\\n"
                        f"{{truncation}}\\n"
                    )

                {function_string.format(1, 0)}
                {function_string.format(2, 1)}
                {function_string.format(3, 2)}

                if not truncation.doubles:
                    return w_1, w_2, w_3, None, None, None
                else:
                    {function_string.format(4, 3)}

                if not truncation.triples:
                    return w_1, w_2, w_3, w_4, None, None
                else:
                    {function_string.format(5, 4)}

                if not truncation.quadruples:
                    return w_1, w_2, w_3, w_4, w_5, None
                else:
                    {function_string.format(6, 5)}

                if not truncation.quintuples:
                    return w_1, w_2, w_3, w_4, w_5, w_6
                else:
                    raise Exception(
                        "Attempting to calculate W^7 operator (quintuples)\\n"
                        "This is currently not implemented!!\\n"
                    )
        '''


    # remove three indents from string block
    lines = string.splitlines()
    trimmed_string = "\n".join([line[4*3:] for line in lines])

    return trimmed_string


# ------------------------------------------------------- #
def _t_term_shape_string(order):
    """Return the string `(A, A, N, ...)` with `order` number of `N`'s."""
    return f"({', '.join(['A','A',] + ['N',]*order)})"


def _contracted_expressions(partition_list, order):
    """Return a list of each of the `oe.contract_expression` calls
    for a W operator of order `order`.
    """
    exp_list = []

    # for each partition (the mathematical term) of the integer `order`
    for partition in partition_list:

        # there is nothing to optimize for the N^th case
        # we simply add the largest t_term to the W operator
        # no multiplication required
        if len(partition) == 1:
            continue

        # the unique permutations of the `partition` of the integer `order`
        combinations = unique_permutations(partition)
        # the einsum indices for the surface dimensions
        surf_index = _generate_surface_index(partition)
        # the einsum indices for the normal mode dimensions
        mode_index = _generate_mode_index(partition, order)

        # `combinations` is a list of tuples such as (2, 2, 1, ) or (5,)
        for i, tupl in enumerate(combinations):
            # the input dimensions are two character strings representing the surface dimensions
            # plus 1 or more characters representing the normal mode dimensions
            in_dims = ", ".join([surf_index[a]+mode_index[i][a] for a in range(len(surf_index)) ])
            # the output dimension is the same for all einsum calls for a given `partition` argument
            out_dims = f"ab{tag_str[0:order]}"
            # the shape of the t_terms are (A, A, N, ...) where the number of N dimensions is
            # determined by the integer elements of each tuple `tupl`
            # so (2, 1) -> ("(A, A, N, N)", "(A, A, N")
            pterms = ", ".join([_t_term_shape_string(n) for n in tupl])
            # print(f"np.einsum('{in_dims}->{out_dims}', {pterms})")
            exp_list.append(f"oe.contract_expression('{in_dims}->{out_dims}', {pterms}),\n")

    return exp_list


def _write_optimized_vemx_paths_function(max_order):
    """Return strings to write all the `oe.contract_expression` calls.
    Unfortunately the code got a lot messier when I had to add in the truncation if statements.
    It should get a rework/factorization at some point
    """
    assert max_order <= 6, "Only implemented up to 6th order"

    string = (
        f"\ndef compute_optimized_vemx_paths(A, N, truncation):\n"
        f'{tab}"""Calculate optimized paths for the einsum calls up to `highest_order`."""\n'
        "\n"
        f"{tab}order_1_list, order_2_list, order_3_list = [], [], []\n"
        f"{tab}order_4_list, order_5_list, order_6_list = [], [], []\n"
        "\n"
    )
    return ""

    optimized_orders = list(range(2, max_order+1))

    for order in optimized_orders:

        # generate all the elements in the `order_{order}_list`
        partitions = generate_un_linked_disconnected_partitions_of_n(order)
        optimized_paths = _contracted_expressions(partitions, order)

        # we always calculate Singles which requires (W^1, W^2, W^3)
        if order == 2:
            string += (
                f"{tab}order_2_list.extend([\n"
                f"{tab}{tab}{optimized_paths[0]}"
                f"{tab}])\n\n"
            )
        elif order == 3:
            string += (
                f"{tab}if truncation.doubles:\n"
                f"{tab}{tab}order_3_list.extend([\n"
                f"{tab}{tab}{tab}{optimized_paths[0]}"
                f"{tab}{tab}{tab}{optimized_paths[1]}"
                f"{tab}{tab}])\n"
                "\n"
                f"{tab}order_3_list.extend([\n"
                f"{tab}{tab}{optimized_paths[2]}"
                f"{tab}])\n\n"
            )


        # for Doubles and higher truncations we want to add multiple`if` statement
        # since we don't need to calculate these optimal paths
        elif order >= 4:
            # the string representation (doubles, triples, quadruples... etc)
            contribution_name = lambda n: taylor_series_order_tag[n].lower()


            string += f"{tab}if truncation.{contribution_name(order-2)}:\n"

            string += (
                f"{tab}{tab}if truncation.{contribution_name(order-1)}:\n"
                f"{tab}{tab}{tab}order_{order}_list.extend([\n"
                f"{tab}{tab}{tab}{tab}{optimized_paths[0]}"
                f"{tab}{tab}{tab}{tab}{optimized_paths[1]}"
                f"{tab}{tab}{tab}])\n\n"
            )

            # we need to a big long string, and also remove the first two opt paths
            optimized_paths = "".join([
                s.replace(f"oe", f"{tab}{tab}{tab}oe") for s in optimized_paths[2:]
            ])

            string += (
                # f"{tab}if truncation.{contribution_name(order-2)}:\n"
                f"{tab}{tab}order_{order}_list.extend([\n"
                f"{optimized_paths}"
                f"{tab}{tab}])\n"
                f"{tab}else:\n"
                f"""{tab}{tab}log.warning("Didn't calculate optimized paths of the W^{order} operator ")\n\n"""
            )


    return_list = ', '.join(['order_1_list',] + [f'order_{order}_list' for order in optimized_orders])

    string += f"{tab}return [{return_list}]\n"

    return string


def _write_optimized_vecc_paths_function(max_order):
    """Return strings to write all the `oe.contract_expression` calls.
    Unfortunately the code got a lot messier when I had to add in the truncation if statements.
    It should get a rework/factorization at some point
    """
    assert max_order <= 6, "Only implemented up to 6th order"

    string = (
        f"\ndef compute_optimized_vecc_paths(A, N, truncation):\n"
        f'{tab}"""Calculate optimized paths for the VECC einsum calls up to `highest_order`."""\n'
        "\n"
        f"{tab}order_1_list, order_2_list, order_3_list = [], [], []\n"
        f"{tab}order_4_list, order_5_list, order_6_list = [], [], []\n"
        "\n"
    )

    # there are no VECC contributions for order < 4
    optimized_vecc_orders = list(range(4, max_order+1))

    for order in optimized_vecc_orders:

        # generate all the elements in the `order_{order}_list`
        partitions = generate_un_linked_disconnected_partitions_of_n(order)
        print('B', partitions)
        optimized_paths = _contracted_expressions(partitions, order)


        # we want to add multiple`if` statements to only calculate paths we need

        # the string representation (doubles, triples, quadruples... etc)
        contribution_name = lambda n: taylor_series_order_tag[n].lower()

        string += f"{tab}if truncation.{contribution_name(order-2)}:\n"

        print(order, order-1, order-2, contribution_name(order-2), contribution_name(order-1))
        print(taylor_series_order_tag)

        string += (
            f"{tab}{tab}if truncation.{contribution_name(order-1)}:\n"
            f"{tab}{tab}{tab}order_{order}_list.extend([\n"
            f"{optimized_paths}"
            f"{tab}{tab}{tab}])\n\n"
        )

        # we need to a big long string, and also remove the first two opt paths
        optimized_paths = "".join([
            s.replace(f"oe", f"{tab}{tab}{tab}oe") for s in optimized_paths[2:]
        ])

        string += (
            # f"{tab}if truncation.{contribution_name(order-2)}:\n"
            f"{tab}{tab}order_{order}_list.extend([\n"
            f"{optimized_paths}"
            f"{tab}{tab}])\n"
            f"{tab}else:\n"
            f"""{tab}{tab}log.warning("Didn't calculate optimized VECC paths of the W^{order} operator")\n\n"""
        )


    return_list = ', '.join(['order_1_list',] + [f'order_{order}_list' for order in optimized_vecc_orders])

    string += f"{tab}return [{return_list}]\n"

    return string


# ------------------------------------------------------- #
def generate_w_operators_string(max_order, s1=75, s2=28):
    """Return a string containing the python code to generate w operators up to (and including) `max_order`.
    Requires the following header: `"import numpy as np\nfrom math import factorial"`.
    """
    spacing_line = "# " + "-"*s1 + " #\n"
    named_line = lambda name, width: "# " + "-"*width + f" {name} " + "-"*width + " #\n"

    # ------------------------------------------------------------------------------------------- #
    # header for default functions (as opposed to the optimized functions)
    string = spacing_line + named_line("DEFAULT FUNCTIONS", s2) + spacing_line
    # ----------------------------------------------------------------------- #
    # header for VECI/CC (MIXED) functions
    string += '\n' + named_line("VECI/CC CONTRIBUTIONS", s2)
    # generate the VECI/CC (MIXED) contributions to W
    string += "".join([_generate_vemx_contributions(order=order) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # header for VECC contributions
    string += '\n' + named_line("VECC CONTRIBUTIONS", s2)
    # generate the VECC contributions to W
    string += "".join([_generate_vecc_contributions(order=order) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # header for w operator functions
    string += '\n' + named_line("W OPERATOR FUNCTIONS", s2)
    # generate the w operator function
    string += "".join([_write_w_function_strings(order=order) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # write the master function
    string += _write_master_w_compute_function(max_order) + "\n"
    # ------------------------------------------------------------------------------------------- #
    #
    # ------------------------------------------------------------------------------------------- #
    # header for optimized functions
    string += spacing_line + named_line("OPTIMIZED FUNCTIONS", s2-1) + spacing_line
    # ----------------------------------------------------------------------- #
    # header for VECI/CC (MIXED) functions
    string += '\n' + named_line("VECI/CC CONTRIBUTIONS", s2)
    # generate the VECI/CC (MIXED) contributions to W
    string += "".join([_generate_vemx_contributions(order=order, opt_einsum=True) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # header for VECC contributions
    string += '\n' + named_line("VECC CONTRIBUTIONS", s2)
    # generate the VECC contributions to W
    string += "".join([_generate_vecc_contributions(order=order, opt_einsum=True) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # header for w operator functions
    string += '\n' + named_line("W OPERATOR FUNCTIONS", s2)
    # generate the w operator function
    string += "".join([_write_w_function_strings(order=order, opt_einsum=True) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # write the master function
    string += _write_master_w_compute_function(max_order, opt_einsum=True) + "\n"
    # ------------------------------------------------------------------------------------------- #
    #
    # ------------------------------------------------------------------------------------------- #
    # header for optimized paths function
    string += '\n' + named_line("OPTIMIZED PATHS FUNCTION", s2)
    # write the code for generating optimized paths for VECI/CC (mixed) contributions
    string += _write_optimized_vemx_paths_function(max_order) + '\n'
    # write the code for generating optimized paths for VECC contributions
    string += _write_optimized_vecc_paths_function(max_order) + '\n'
    # ------------------------------------------------------------------------------------------- #
    return string


def generate_w_operator_equations_file(max_w_order, path="./w_operator_equations.py"):
    """Generates and saves to a file the code to calculate the w operator equations for the CC approach."""

    # start with the import statements
    file_data = (
        "# system imports\n"
        "from math import factorial\n"
        "\n"
        "# third party imports\n"
        "import numpy as np\n"
        "import opt_einsum as oe\n"
        "\n"
        "# local imports\n"
        "from .symmetrize import symmetrize_tensor\n"
        "from ..log_conf import log\n"
        "\n"
    )

    # write the functions to calculate the W operators
    file_data += generate_w_operators_string(max_order=max_w_order)

    # save data
    with open(path, 'w') as fp:
        fp.write(file_data)

    return

# ----------------------------------------------------------------------------------------------- #
# ---------------------  GENERATING T AMPLITUDE DERIVATIVE EQUATIONS  --------------------------- #
# ----------------------------------------------------------------------------------------------- #


dt_terms = [
    None,
    t_term_namedtuple("dt_i", 1, "(A, A, N)"),
    t_term_namedtuple("dt_ij", 2, "(A, A, N, N)"),
    t_term_namedtuple("dt_ijk", 3, "(A, A, N, N, N)"),
    t_term_namedtuple("dt_ijkl", 4, "(A, A, N, N, N, N)"),
    t_term_namedtuple("dt_ijklm", 5, "(A, A, N, N, N, N, N)"),
    t_term_namedtuple("dt_ijklmn", 6, "(A, A, N, N, N, N, N, N)"),
    t_term_namedtuple("dt_ijklmno", 7, "(A, A, N, N, N, N, N, N, N)"),
    t_term_namedtuple("dt_ijklmnop", 8, "(A, A, N, N, N, N, N, N, N, N)"),
]


# ----------------------------------------------------------------------------------------------- #
def _generate_disconnected_einsum_operands_list(dt_index, tupl):
    """Generate a string representing the list of operands for an einsum call.
    the names of the arguments to the einsum call are stored in the lists `t_terms`, `dt_terms`
    and the objects are accessed using each integer in the tuple (2, 1) -> ("t_ij", "dt_i")
    """
    term_list = []
    for t_index, num in enumerate(tupl):
        # print(f"{dt_index=} {t_index=} {num=} {tupl=}")
        if dt_index == t_index:
            term_list.append(dt_terms[num].string)
        else:
            term_list.append(t_terms[num].string)

    einsum_term_list = ", ".join(term_list)
    return einsum_term_list


def _generate_disconnected_einsum_function_call_list(partition, order):
    """Returns a list of strings. Each string represents a call to np.einsum.
    The list is filled relative to the specific `partition` that is being calculated.
    here we want to do all permutations like so:
    2, 1, 1
    1, 2, 1
    1, 1, 2
    but in addition we want to permute which of them is dt term
    and because dt and t don't commute we need to count all perms
    d2, 1, 1
    1, d2, 1
    1, 1, d2
    --------
    d1, 1, 2
    d1, 2, 1
    --------
    1, d1, 2
    2, d1, 1
    --------
    1, 2, d1
    2, 1, d1
    """

    # the unique permutations of the `partition` of the integer `order`
    combinations = unique_permutations(partition)
    # the einsum indices for the surface dimensions
    surf_index = _generate_surface_index(partition)
    # the einsum indices for the normal mode dimensions this is different than the w generators
    mode_index = _generate_mode_index(partition, order)

    return_list = []  # store return values here

    # `combinations` is a list of tuples such as (2, 2, 1, ) or (5,)
    for i, tupl in enumerate(combinations):
        # the input dimensions are two character strings representing the surface dimensions
        # plus 1 or more characters representing the normal mode dimensions
        in_dims = ", ".join([surf_index[a]+mode_index[i][a] for a in range(len(surf_index)) ])
        # the output dimension is the same for all einsum calls for a given `partition` argument
        out_dims = f"ab{tag_str[0:order]}"

        # the names of the arguments to the einsum call are stored in the lists `t_terms`, `dt_terms`
        # and the objects are accessed using each integer in the tuple (2, 1) -> ("t_ij", "dt_i")
        for dt_index in range(len(partition)):
            operands = _generate_disconnected_einsum_operands_list(dt_index, tupl)
            return_list.append(f"np.einsum('{in_dims}->{out_dims}', {operands})")


    return return_list


def _generate_optimized_disconnected_einsum_function_call_list(partition, order, iterator_name='optimized_einsum'):
    """Returns a list of strings. Each string represents a call to np.einsum.
    The list is filled relative to the specific `partition` that is being calculated.

    This function, `_construct_linked_disconnected_definition`, and `_construct_un_linked_disconnected_definition`
    need to agree on the name of the iterator containing the optimized paths. At the moment it is named `optimized_einsum`.
    """

    # the unique permutations of the `partition` of the integer `order`
    combinations = unique_permutations(partition)

    return_list = []  # store return values here

    # `combinations` is a list of tuples such as (2, 2, 1, ) or (5,)
    for i, tupl in enumerate(combinations):
        # the names of the arguments to the einsum call are stored in the lists `t_terms`, `dt_terms`
        # and the objects are accessed using each integer in the tuple (2, 1) -> ("t_ij", "dt_i")
        for dt_index in range(len(partition)):
            operands = _generate_disconnected_einsum_operands_list(dt_index, tupl)
            return_list.append(f"next({iterator_name})({operands})")

    return return_list


# ----------------------------------------------------------------------------------------------- #
def _construct_linked_disconnected_definition(return_string, order, opt_einsum=False, iterator_name='optimized_einsum'):
    """Return the string containing the python code to prepare the function definition and unpack the t_args.
    This function and `_optimized_linked_disconnected_einsum_list` need to agree on the name of the
    iterator containing the optimized paths. At the moment it is named `optimized_einsum`.
    """

    # this line of python code labels the terms in the `t_args` tuple
    t_arg_string = ", ".join([t_terms[n].string for n in range(1, order)] + ["*unusedargs"])
    # this line of python code labels the terms in the `dt_args` tuple
    dt_arg_string = ", ".join([dt_terms[n].string for n in range(1, order)] + ["*unusedargs"])

    if not opt_einsum:
        return_string += (
            f"\ndef _order_{order}_linked_disconnected_terms(A, N, trunc, t_args, dt_args):\n"
            f'{tab}"""Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.\n'
            f'{tab}This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)\n'
            f'{tab}But not terms (5), (3, 2), (2, 2, 1)\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args` and 'dt_args'\n"
            f"{tab}{t_arg_string} = t_args\n"
            f"{tab}{dt_arg_string} = dt_args\n"
        )
    else:
        return_string += (
            f"\ndef _order_{order}_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):\n"
            f'{tab}"""Calculate all uniquely linked disconnected terms generated from the wave operator ansatz.\n'
            f'{tab}This means for order 5 we include terms such as (4, 1), (3, 1, 1), (2, 1, 1, 1)\n'
            f'{tab}But not terms (5), (3, 2), (2, 2, 1), (1, 1, 1, 1, 1)\n'
            f'{tab}Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args` and 'dt_args'\n"
            f"{tab}{t_arg_string} = t_args\n"
            f"{tab}{dt_arg_string} = dt_args\n"
        )
        return_string += (
            f"{tab}# make an iterable out of the `opt_path_list`\n"
            f"{tab}{iterator_name} = iter(opt_path_list)\n"
        )

    return return_string


def _write_linked_disconnected_strings(order, opt_einsum=False):
    """Output the python code which generates all CI/CC terms to subtract from dt."""

    # may want to switch to numerical arguments for einsum
    assert order <= 6, "Can only handle up to 6th order due to `einsum_surface_tags` limit"
    return_string = "" # we store the output in this string

    # prepare the function definition, unpacking of arguments, and initialization of the W array
    return_string += _construct_linked_disconnected_definition(return_string, order, opt_einsum)

    # special exception for order less than 2
    if order < 2:
        lst = return_string.split('"""')
        exception = ( "\n"
            f"{tab}raise Exception(\n"
            f'{tab}{tab}"the first possible purely VECI/CC term is (1, 1) or (dt_i * t_i)"\n'
            f'{tab}{tab}"which requires a residual of at least 2nd order"\n'
            f"{tab})\n"
        )
        return '\"\"\"'.join([lst[0], "Exists for error checking.", exception])

    # name of the return array
    array = "linked_disconnected_terms"

    # initialize the return array
    return_string += (
        f"{tab}# Creating the {num_tag[order]} order return array\n"
        f"{tab}{array} = np.zeros(({', '.join(['A','A',] + ['N',]*order)}), dtype=complex)\n"
    )

    for partition in generate_linked_disconnected_partitions_of_n(order):

        # Label the order of this partition's contribution
        return_string += f"{tab}# the {partition} term\n"

        # make the prefactor
        prefactor = _generate_w_operator_prefactor(partition)

        # compile a list of the einsum calls
        if opt_einsum:
            einsum_list = _generate_optimized_disconnected_einsum_function_call_list(partition, order)
        else:
            einsum_list = _generate_disconnected_einsum_function_call_list(partition, order)

        # join them together
        sum_of_einsums = f" +\n{tab}{tab}".join(einsum_list)

        # place the einsums between the parentheses
        return_string += (
            f"{tab}{array} += {prefactor} * (\n"
            f"{tab}{tab}{sum_of_einsums}\n"
            f"{tab})\n"
        )

    return_string +=  f"\n{tab}return {array}\n"

    return return_string


# ----------------------------------------------------------------------------------------------- #
def _construct_un_linked_disconnected_definition(return_string, order, opt_einsum=False, iterator_name='optimized_einsum'):
    """Return the string containing the python code to prepare the function definition and unpack the t_args.
    This function and `_optimized_un_linked_disconnected_einsum_list` need to agree on the name of the
    iterator containing the optimized paths. At the moment it is named `optimized_einsum`.
    """

    # this line of python code labels the terms in the `t_args` tuple
    t_arg_string = ", ".join([t_terms[n].string for n in range(1, order)] + ["*unusedargs"])
    # this line of python code labels the terms in the `dt_args` tuple
    dt_arg_string = ", ".join([dt_terms[n].string for n in range(1, order)] + ["*unusedargs"])

    if not opt_einsum:
        return_string += (
            f"\ndef _order_{order}_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args):\n"
            f'{tab}"""Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.\n'
            f'{tab}This means for order 5 we include terms such as (3, 2), (2, 2, 1)\n'
            f'{tab}But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args` and 'dt_args'\n"
            f"{tab}{t_arg_string} = t_args\n"
            f"{tab}{dt_arg_string} = dt_args\n"
        )
    else:
        return_string += (
            f"\ndef _order_{order}_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list):\n"
            f'{tab}"""Calculate all uniquely un-linked disconnected terms generated from the wave operator ansatz.\n'
            f'{tab}This means for order 5 we include terms such as (3, 2), (2, 2, 1)\n'
            f'{tab}But not terms (5), (4, 1), (3, 1, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)\n'
            f'{tab}Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `t_args` and 'dt_args'\n"
            f"{tab}{t_arg_string} = t_args\n"
            f"{tab}{dt_arg_string} = dt_args\n"
        )
        return_string += (
            f"{tab}# make an iterable out of the `opt_path_list`\n"
            f"{tab}{iterator_name} = iter(opt_path_list)\n"
        )

    return return_string


def _write_un_linked_disconnected_strings(order, opt_einsum=False):
    """Output the python code which generates all CC terms to subtract from dt."""

    # may want to switch to numerical arguments for einsum
    assert order <= 6, "Can only handle up to 6th order due to `einsum_surface_tags` limit"
    return_string = "" # we store the output in this string

    # prepare the function definition, unpacking of arguments, and initialization of the W array
    return_string += _construct_un_linked_disconnected_definition(return_string, order, opt_einsum)

    # special exception for order less than 2
    if order < 4:
        lst = return_string.split('"""')
        exception = ( "\n"
            f"{tab}raise Exception(\n"
            f'{tab}{tab}"the first possible purely VECC term is (2, 2) or (dt_ij * t_ij)"\n'
            f'{tab}{tab}"which requires a residual of at least 4th order"\n'
            f"{tab})\n"
        )
        return '\"\"\"'.join([lst[0], "Exists for error checking.", exception])

    # name of the return array
    array = "un_linked_disconnected_terms"

    # initialize the return array
    return_string += (
        f"{tab}# Creating the {num_tag[order]} order return array\n"
        f"{tab}{array} = np.zeros(({', '.join(['A','A',] + ['N',]*order)}), dtype=complex)\n"
    )

    for partition in generate_un_linked_disconnected_partitions_of_n(order):

        # Label the order of this partition's contribution
        return_string += f"{tab}{tab}# the {partition} term\n"

        # make the prefactor
        prefactor = _generate_w_operator_prefactor(partition)

        # compile a list of the einsum calls
        if opt_einsum:
            einsum_list = _generate_optimized_disconnected_einsum_function_call_list(partition, order)
        else:
            einsum_list = _generate_disconnected_einsum_function_call_list(partition, order)

        # join them together
        sum_of_einsums = f" +\n{tab}{tab}".join(einsum_list)

        # place the einsums between the parentheses
        return_string += (
            f"{tab}{array} += {prefactor} * (\n"
            f"{tab}{tab}{sum_of_einsums}\n"
            f"{tab})\n"
        )

    return_string +=  f"\n{tab}return {array}\n"

    return return_string


# ----------------------------------------------------------------------------------------------- #
def _construct_dt_amplitude_definition(return_string, order, opt_einsum=False, iterator_name='optimized_einsum'):
    """Return the string containing the python code to prepare the function definition,
    unpack the t_args, and initialize the W array.

    This function and `_construct_dt_amplitude_definition` need to agree on the name of the
    iterator containing the optimized paths. At the moment it is named `optimized_einsum`.
    """

    # this line of python code labels the terms in the `t_args` tuple
    w_arg_string = ", ".join([w_dict[n] for n in range(1, order+1)] + ["*unusedargs"])

    if not opt_einsum:
        return_string += (
            f"\ndef _calculate_order_{order}_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):\n"
            f'{tab}"""Calculate the derivative of the {order} t-amplitude for use in the calculation of the residuals."""\n'
            f"{tab}# unpack the `w_args`\n"
            f"{tab}{w_arg_string} = w_args\n"
        )
    else:
        return_string += (
            f"\ndef _calculate_order_{order}_dt_amplitude_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_epsilon, opt_path_list):\n"
            f'{tab}"""Calculate the derivative of the {order} t-amplitude for use in the calculation of the residuals.\n'
            f'{tab}Uses optimized summation paths generate using `contract_expression` from the `opt_einsum` library.\n'
            f'{tab}"""\n'
            f"{tab}# unpack the `w_args`\n"
            f"{tab}{w_arg_string} = w_args\n"
        )
        if order >= 2:
            return_string += (
                f"{tab}# make an iterable out of the `opt_path_list`\n"
                f"{tab}{iterator_name} = iter(opt_path_list)\n"
            )

    return return_string


def _write_dt_amplitude_strings(order, opt_einsum=False):
    """ output the string of the function `compute_CI_from_CC` which generates all the W's """

    # may want to switch to numerical arguments for einsum
    assert order <= 6, "Can only handle up to 6th order due to `einsum_surface_tags` limit"
    return_string = "" # we store the output in this string

    # prepare the function definition, unpacking of arguments, and initialization of the W array
    return_string += _construct_dt_amplitude_definition(return_string, order, opt_einsum)

    mode_subscripts = f"{tag_str[0:order]}"

    # initialize the W array
    return_string += (
        f"{tab}# Calculate the {num_tag[order]} order residual\n"
        f"{tab}residual = residual_equations.calculate_order_{order}_residual(A, N, trunc, h_args, w_args)\n"
    )


    # generate the appropriate strings `vemx_operations` and `vecc_operations
    # which are used in a f-string about 50 lines below
    if not opt_einsum:
        prefactor = f"1/factorial({order})"

        epsilon_term = f"{prefactor} * np.einsum('ac{mode_subscripts},cb->ab{mode_subscripts}', w_{mode_subscripts}, epsilon)"

        if order < 2:
            vemx_operations = f"{tab}{tab}pass  # no linked disconnected terms for order < 2\n"
            vecc_operations = f"{tab}{tab}pass  # no un-linked disconnected terms for order < 4\n"
        else:
            vemx_operations = (
                f"{tab}{tab}residual -= _order_{order}_linked_disconnected_terms(A, N, trunc, t_args, dt_args)\n"
            )
            if order < 4:
                vecc_operations = (
                    f"{tab}{tab}residual -= _order_{order}_linked_disconnected_terms(A, N, trunc, t_args, dt_args)\n"
                    f"{tab}{tab}pass  # no un-linked disconnected terms for order < 4\n"
                )
            else:
                vecc_operations = (
                    f"{tab}{tab}residual -= _order_{order}_linked_disconnected_terms(A, N, trunc, t_args, dt_args)\n"
                    f"{tab}{tab}residual -= _order_{order}_un_linked_disconnected_terms(A, N, trunc, t_args, dt_args)\n"
                )
    else:
        prefactor = f"1/factorial({order})"

        epsilon_term = f"{prefactor} * opt_epsilon(w_{mode_subscripts}, epsilon)"

        if order < 2:
            vemx_operations = f"{tab}{tab}pass  # no linked disconnected terms for order < 2\n"
            vecc_operations = f"{tab}{tab}pass  # no un-linked disconnected terms for order < 4\n"
        else:
            vemx_operations = (
                f"{tab}{tab}residual -= _order_{order}_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)\n"
            )
            if order < 4:
                vecc_operations = (
                    f"{tab}{tab}residual -= _order_{order}_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)\n"
                    f"{tab}{tab}pass  # no un-linked disconnected terms for order < 4\n"
                )
            else:
                vecc_operations = (
                    f"{tab}{tab}residual -= _order_{order}_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)\n"
                    f"{tab}{tab}residual -= _order_{order}_un_linked_disconnected_terms_optimized(A, N, trunc, t_args, dt_args, opt_path_list)\n"
                )


    return_string += (
        # subtract epsilon term
        f"{tab}# subtract the epsilon term (which is R_0)\n"
        f"{tab}residual -= {epsilon_term}\n"
        "\n"
        # disconnected terms
        f"{tab}# subtract the disconnected terms\n"
        f"{tab}if ansatz.VECI:\n"
        f"{tab}{tab}pass  # veci does not include any disconnected terms\n"
        f"{tab}elif ansatz.VE_MIXED:\n"
        f"{vemx_operations}"
        f"{tab}elif ansatz.VECC:\n"
        f"{vecc_operations}\n"
    )

    return_string += (
        f"{tab}# Symmetrize the residual operator\n"
        f"{tab}dt_{mode_subscripts} = symmetrize_tensor(N, residual, order={order})\n"
        f"{tab}return dt_{mode_subscripts}\n"
    )

    return return_string


def _write_master_dt_amplitude_function(order, opt_einsum=False):
    """ x """
    t_term = t_terms[order].string
    dt_term = dt_terms[order].string
    name = taylor_series_order_tag[order]

    if not opt_einsum:
        string = f'''
            def solve_{name}_equations(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args):
                """Compute the change in the {t_term} term ({name})"""

                if not trunc.at_least_{name}:
                    raise Exception(
                        "It appears that {name} is not true, this cannot be."
                        "Something went terribly wrong!!!"
                    )
                {dt_term} = _calculate_order_{order}_dt_amplitude(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args)
                return {dt_term}
        '''
    else:
        string = f'''
            def solve_{name}_equations_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_path_list):
                """Compute the change in the {t_term} term ({name})"""

                if not trunc.at_least_{name}:
                    raise Exception(
                        "It appears that {name} is not true, this cannot be."
                        "Something went terribly wrong!!!"
                    )
                {dt_term} = _calculate_order_{order}_dt_amplitude_optimized(A, N, ansatz, trunc, epsilon, h_args, t_args, dt_args, w_args, opt_path_list)
                return {dt_term}
        '''


    # remove three indents from string block
    lines = string.splitlines()
    trimmed_string = "\n".join([line[4*3:] for line in lines])

    return trimmed_string


# ----------------------------------------------------------------------------------------------- #
def _write_optimized_dt_amplitude_paths_function(max_order):
    """Return strings to write all the `oe.contract_expression` calls.
    Unfortunately the code got a lot messier when I had to add in the truncation if statements.
    It should get a rework/factorization at some point
    """
    assert max_order < 7, "optimized paths only implemented up to 6th order"

    string = (
        f"\ndef compute_optimized_paths(A, N, truncation):\n"
        f'{tab}"""Calculate optimized paths for the einsum calls up to `highest_order`."""\n'
    )
    # we need three optimized lists:
    # 1 - optimized epsilon call
    # 2 - optimized linked disconnected calls
    # 3 - optimized un-linked disconnected calls

    optimized_orders = list(range(2, max_order+1))

    for order in optimized_orders:
        pass

    string += (
        "\n"
        f"{tab}order_1_list, order_2_list, order_3_list = [], [], []\n"
        f"{tab}order_4_list, order_5_list, order_6_list = [], [], []\n"
        "\n"
    )

    # return_list = ', '.join(['order_1_list',] + [f'order_{order}_list' for order in optimized_orders])
    return_list = None
    string += f"{tab}return [{return_list}]\n"

    return string


def generate_dt_amplitude_string(max_order, s1=75, s2=28):
    """Return a string containing the python code to generate dt up to (and including) `max_order`.
    Requires the following header: `"import numpy as np\nfrom math import factorial"`.
    """
    spacing_line = "# " + "-"*s1 + " #\n"
    named_line = lambda name, width: "# " + "-"*width + f" {name} " + "-"*width + " #\n"

    # ------------------------------------------------------------------------------------------- #
    # header for default functions
    string = spacing_line + named_line("DEFAULT FUNCTIONS", s2) + spacing_line
    # ----------------------------------------------------------------------- #    # header for VECI functions
    string += '\n' + named_line("DISCONNECTED TERMS", s2)
    # generate the linked disconnected function
    string += "".join([_write_linked_disconnected_strings(order=order) for order in range(1, max_order+1)])
    # generate the un-linked disconnected function
    string += "".join([_write_un_linked_disconnected_strings(order=order) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # header for VECI functions
    string += '\n' + named_line("dt AMPLITUDES", s2)
    # generate the default VECI functions
    string += "".join([_write_dt_amplitude_strings(order=order) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # header for master functions
    string += '\n' + named_line("WRAPPER FUNCTIONS", s2)
    # write the master functions
    string += "".join([_write_master_dt_amplitude_function(order) for order in range(1, max_order+1)])

    string += '\n'
    # ------------------------------------------------------------------------------------------- #
    # ------------------------------------------------------------------------------------------- #
    # header for optimized functions
    string += spacing_line + named_line("OPTIMIZED FUNCTIONS", s2-1) + spacing_line
    # ----------------------------------------------------------------------- #    # header for VECI functions
    string += '\n' + named_line("DISCONNECTED TERMS", s2)
    # generate the linked disconnected function
    string += "".join([_write_linked_disconnected_strings(order=order, opt_einsum=True) for order in range(1, max_order+1)])
    # generate the un-linked disconnected function
    string += "".join([_write_un_linked_disconnected_strings(order=order, opt_einsum=True) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # header for VECI functions
    string += '\n' + named_line("dt AMPLITUDES", s2)
    # generate the optimized functions
    string += "".join([_write_dt_amplitude_strings(order=order, opt_einsum=True) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    # header for master functions
    string += '\n' + named_line("WRAPPER FUNCTIONS", s2)
    # write the master functions
    string += "".join([_write_master_dt_amplitude_function(order) for order in range(1, max_order+1)])
    # ----------------------------------------------------------------------- #
    # header for optimized paths function
    string += '\n' + named_line("OPTIMIZED PATHS FUNCTION", s2)
    # write the code for generating the optimized paths
    string += _write_optimized_dt_amplitude_paths_function(max_order)

    string += '\n'
    # ------------------------------------------------------------------------------------------- #
    return string


def generate_dt_amplitude_equations_file(max_w_order, path="./dt_amplitude_equations.py"):
    """Generates and saves to a file the code to calculate the t-amplitude derivative equations for the CC approach."""

    # start with the import statements
    file_data = (
        "# system imports\n"
        "from math import factorial\n"
        "\n"
        "# third party imports\n"
        "import numpy as np\n"
        "import opt_einsum as oe\n"
        "\n"
        "# local imports\n"
        "from ..log_conf import log\n"
        "from .symmetrize import symmetrize_tensor\n"
        "from . import residual_equations\n"
        "\n"
    )

    # write the functions to calculate the derivatives of the t-amplitudes
    file_data += generate_dt_amplitude_string(max_order=max_w_order)

    # save data
    with open(path, 'w') as fp:
        fp.write(file_data)

    return


# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #

if (__name__ == '__main__'):
    # generate_residual_equations_latex(5)

    max_residual_order = 6
    # generate_residual_equations_file(max_residual_order)
    max_w_order = 6
    generate_w_operator_equations_file(max_w_order)
    dt_order = 6
    # generate_dt_amplitude_equations_file(dt_order)

