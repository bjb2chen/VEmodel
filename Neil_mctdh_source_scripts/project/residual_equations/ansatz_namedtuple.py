"""Rather than using lists,dictionaries,bools or enums using namedtuples makes the code much more concise
we can write things like `self.ansatz.VECC` and the label of the member explicitly
describes what value it contains making the code more readable and user friendly

The way the `ansatz_namedtuple` works is that when initializing a `vibronic_hamiltonian` object
one of its attributes is an instance of `ansatz_namedtuple`.
"""

# system imports
from collections import namedtuple
# third party imports
# local imports

field_names = [
    # fieldnames for the ansatz of the parameterization of the Hamiltonian
    "dS",
    "new_1",
    # fieldnames for the ansatz of the wave operator
    "VECC",
    "VECI",
    "VE_MIXED",
]

defaults = [False]*len(field_names)

ansatz = namedtuple('ansatz_namedtuple', field_names, defaults=defaults)
