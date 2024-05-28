# really simple module to store the choice of boltzman constant
# between all files
#
#
#
import numpy as np

# reference http://ws680.nist.gov/publication/get_pdf.cfm?pub_id=920687
# reference http://ws680.nist.gov/publication/get_pdf.cfm?pub_id=920686

# joules / ev
nist_j_per_ev = np.float64(1.6021766208e-19)

# 1 / mol
avogadro = np.float64(6.022140857e+23)

# ev / kelvin
boltzman = np.float64(1.38064852e-23) / nist_j_per_ev

hbar = 1.0

# convert wavenumbers to electronVolts and back
wavenumber_per_eV = 8065.6


# TODO - the naming of these two functions could probably be improved
def beta(temperature):
    """returns a value of beta for a given temperature in Kelvin"""
    val = 1. / (temperature * boltzman)
    return val


def extract_T_from_beta(beta):
    """returns a temperature (in Kelvin) for a given value of beta"""
    val = 1. / (beta * boltzman)
    return val
