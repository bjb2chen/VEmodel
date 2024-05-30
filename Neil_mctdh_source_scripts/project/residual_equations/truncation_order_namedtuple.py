"""Rather than using lists,dictionaries,bools or enums using namedtuples makes the code much more concise
we can write things like `self.truncation.at_least_quadratic` and the label of the member explicitly
describes what value it contains making the code more readable and user friendly

The way the `truncation_order_namedtuple` works is that when initializing a `vibronic_hamiltonian` object
one of its members is an instance of `truncation_order_namedtuple`. The arguments we provide to
`create_truncation_order_object` are the 'cc_truncation_order' and `hamiltonian_truncation_order` parameters.
Then we set the appropriate booleans to be true/false and they are used in the following fashion:
```
    if self.trunc.at_least_linear:
        # do some linear stuff
        if self.trunc.at_least_quadratic:
            # do some quadratic stuff

    if self.trunc.singles:
        # do singles stuff
```
We don't have `at_least_*` structure for the CC since BY DEFINITION they are inclusive
(if doubles, then singles by definition)
"""

# system imports
from collections import namedtuple
# third party imports
# local imports


# CC strings (booleans)
cc_string_list = ['singles', 'doubles', 'triples', 'quadruples', 'quintuples', 'sextuples']
# Hamiltonian strings (booleans)
ham_string_list = ['linear', 'quadratic', 'cubic', 'quartic', 'quintic']

field_names = [
    # required initialization arguments (integers)
    'cc_truncation_order', 'hamiltonian_truncation_order',
    # CC fieldnames (booleans)
    *cc_string_list,
    # Hamiltonian fieldnames (booleans)
    *ham_string_list,
    # `at_least_*` fieldnames (booleans)
    *[f"at_least_{word}" for word in ham_string_list],
]

# we only want the first two field names to be required
# so our list of default values needs to be 2 less than the length of `field_names`
defaults = [False]*(len(field_names)-2)


class truncation_order_namedtuple(namedtuple('truncation_order_namedtuple', field_names, defaults=defaults)):
    """ Add some useful functions to the `namedtuple`. """
    def cc_abbreviation(self):
        """Return the abbreviation for the type of CC calculation."""
        # print([word[0] for word in cc_string_list])
        # print([word.capitalize() for word in cc_string_list])
        # print([word.capitalize()[0] for word in cc_string_list[0:self.cc_truncation_order]])
        return "".join([word.capitalize()[0] for word in cc_string_list[0:self.cc_truncation_order]])



def create_truncation_order_object(cc_truncation_order, hamiltonian_truncation_order):
    """A more flexible "initialization" function for `truncation_order_namedtuple`.
    Assigns `True`/`False` to appropriate attributes based on input args.
    User can provide an integer value or string for `cc_truncation_order`.
    The same thing is true of `hamiltonian_truncation_order`.
    We map to an integer value using the two lists (`cc_list`, `ham_list`)
    """

    kwargs = {} # will store the input arguments to `truncation_order_namedtuple` here

    # map the string to an int (if string provided)
    if cc_truncation_order in cc_string_list:
        cc_truncation_order = cc_string_list.index(cc_truncation_order)+1
    # map the string to an int (if string provided)
    if hamiltonian_truncation_order in ham_string_list:
        hamiltonian_truncation_order = ham_string_list.index(hamiltonian_truncation_order)+1
    # set the top level args
    kwargs['cc_truncation_order'] = cc_truncation_order
    kwargs['hamiltonian_truncation_order'] = hamiltonian_truncation_order

    # set the SDTQ attributes
    for i, word in enumerate(cc_string_list):
        if cc_truncation_order-1 >= i:
            kwargs[word] = True
        else:
            kwargs[word] = False

    # set the Hamiltonian order
    kwargs[ham_string_list[hamiltonian_truncation_order-1]] = True

    # set the `at_least_*` booleans
    for i, word in enumerate(ham_string_list):
        if hamiltonian_truncation_order-1 >= i:
            kwargs[f"at_least_{word}"] = True
        else:
            kwargs[f"at_least_{word}"] = False

    return truncation_order_namedtuple(**kwargs)
