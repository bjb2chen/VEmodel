run_begin = "RUN-SECTION"
run_end = "end-run-section"

basic_propagation = "propagation"
basic_propagation_scheme = "tout = {tout:.2f} tfinal = {tfinal:.2f} tpsi=1.0"

generate_initial_wavefunction = "geninwf"

# this specifies the precision of the wavefunction when written to the output
wavefunction_single_precision = "psi=single"
# wavefunction_double_precision = "psi=double"

# this specifies the 1st order
auto_scheme = "auto=once"
# second order is only needed for filter-diagonalization scheme
# auto_scheme = "auto=twice"

run_section_propagation = "\n".join([
    run_begin,
    "title = wavefunction propagation of {name:s}",
    "name = {name:s}",
    basic_propagation,
    basic_propagation_scheme,
    generate_initial_wavefunction,
    wavefunction_single_precision,
    auto_scheme,
    run_end,
])


# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------


operator_begin = "OPERATOR-SECTION"
operator_end = "end-operator-section"

operator_section = "\n".join([
    operator_begin,
    "opname = {opfile_name:s}",
    operator_end,
])


# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------


spf_begin = "SPF-BASIS-SECTION"
spf_end = "end-spf-basis-section"


def _generate_basic_spf_basis(spf_definitions, multi=False):
    """Generate basic spf basis with same # of spf for all states and modes"""
    return "\n".join([
        spf_begin,
        "multi-set" if multi else "single-set",
        *spf_definitions,
        spf_end,
    ])

def generate_basic_single_set_spf_basis_section(n_BF, N):
    """Generate basic single set with same # of spf for all states and modes"""
    spf_definitions  = [
        f"      m{n:<3d}      =  {n_BF:d}" for n in N.values()
    ]
    return _generate_basic_spf_basis(spf_definitions)

def generate_basic_multi_set_spf_basis_section(n_BF, N, A):
    """Generate basic multi set with same # of spf for all states and modes"""
    spf_definitions  = []
    for n in N.values():
        string = f"      m{n:<3d}      =  1"
        for a in range(A-1):
            string += f", {n_BF:d}"
        spf_definitions.append(string)
    return _generate_basic_spf_basis(spf_definitions, multi=True)


spf_basis_section = "\n".join([
    spf_begin,
    "multi-set",
    "      m1       =  {other_spf:d}, {max_spf:d}, {other_spf:d}",
    "      m2       =  {other_spf:d}, {max_spf:d}, {other_spf:d}",
    "      m3       =  {other_spf:d}, {max_spf:d}, {other_spf:d}",
    spf_end,
])

single_set_spf_basis_section = "\n".join([
    spf_begin,
    "single-set",
    "      m1       =  {n_BF:d}",
    "      m2       =  {n_BF:d}",
    "      m3       =  {n_BF:d}",
    spf_end,
])



# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------


pbs_begin = "PRIMITIVE-BASIS-SECTION"
pbs_end = "end-primitive-basis-section"
electronic_basis = "    el      el     {:d}"

# this is the basis specification where the first element represents...
# the second element is the frequency and the third element is the mass
ho_spec = "0.0   1.0   1.0"
# this specification allows us to define the HO by the first and last grid points
ho_spec_bounded = "xi-xf   0.0   1.0"


def _generate_basic_harmonic_oscillator_pbfs(pbfs_definitions, nof_electronic_states):
    """Generate basic pbfs with same # of H.O. basis functions for all modes"""
    return "\n".join([
        pbs_begin,
        *pbfs_definitions,
        electronic_basis.format(nof_electronic_states),
        pbs_end,
    ])

def generate_basic_harmonic_oscillator_primative_basis_section(n_BF, N, A):
    """Generate PBF's with same # of H.O. basis functions for all modes"""
    pbf_definitions  = [
        f"    m{n:<3d}    HO     {n_BF:d}   {ho_spec:s}" for n in N.values()
    ]
    return _generate_basic_harmonic_oscillator_pbfs(pbf_definitions, A)

# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------


int_begin = "INTEGRATOR-SECTION"
int_end = "end-integrator-section"

# spacing constants
cols = [10, 5, 10, 10]

# mean fields are changed based on the errors of both the MCTDH coefficients AND the single particle functions
# default for propagation and simple relaxation
var = "var"

# mean fields are changed based on the errors of ONLY the single particle functions
# default for block improved relaxation
varphi = "varphi"

cmf_accuracy = "1.0d-05"
cmf_initial_time_interval = 0.5
constant_mean_field = "{cmf_spec:<10} = {R:>5}, {R1:>10}".format(
    cmf_spec=f"CMF/{var}",
    R=cmf_initial_time_interval,
    R1=cmf_accuracy
)

# Bulirsch-Stoer extrapolation integrator
bs_initial_step_size = "2.5d-04"
bs_accuracy = "1.0d-05"
bs_maximal_order = 7
bs_integrator = "{bs_spec:<10} = {I:>5}, {R:>10}, {R1:>10}".format(
    bs_spec='BS/spf',  # could also be `all` or `A`
    I=bs_maximal_order,
    R=bs_accuracy,
    R1=bs_initial_step_size,
)

# SIL integrator
sil_accuracy = "1.0d-05"
sil_maximal_order = 5
sil_integrator = "{sil_spec:<10} = {I:>5}, {R:>10}".format(
    sil_spec='SIL/A',  # could also be `all` or `A`
    I=sil_maximal_order,
    R=sil_accuracy,
)


# Runge-Kutta integrator (order 8, can also replace 8 with 5)
rk_initial_step_size = "0.1"
rk_accuracy = "1.0d-8"
rk_integrator = "{rk_spec:<10} = {R:>10}, {R1:>10}".format(
    rk_spec='RK8/spf',  # could also be `all` or `A`
    R=rk_accuracy,
    R1=rk_initial_step_size,
)

"""
Davidson integrator !!! FOR IMPROVED RELAXATION ONLY !!!
    -   DAV is general (complex Hamiltonians and wavefunctions)
    -   rDav is for real Hamiltonians and real wavefunctions - FASTER
    -   rrDav is rDAV plus more real arithmetic (not as general) - FASTEST
    -   cDav for non-Hermitian Hamiltonians
"""
dav_accuracy = "1.0d-09"
dav_maximal_order = "760"
dav_integrator = "{dav_spec:<10} = {I:>5}, {R:>10}".format(
    dav_spec='rrDAV',
    I=dav_maximal_order,
    R=dav_accuracy,
)

""" The integrators may only be combined in certain ways as described below!!!
For VMF calculations, only the BS, ABM or RK5/8 integrator may be used for the differential equations. Default is ABM. For VMF calculations the integrator must carry the extension S=all (or no extension at all), i.e. there is only one integrator within the VMF scheme. For CMF calculations, the following combinations of integrators are possible: ABM/spf + SIL/A, BS/spf + SIL/A, RKx/spf + SIL/A, BS/all, ABM/all, RKx/all.
"""


propagation_integrator_section = "\n".join([
    int_begin,
    constant_mean_field,
    bs_integrator,
    sil_integrator,
    int_end,
])

# the value of eps_inv is used to regularise the inverse of the reduced density matrices. (Default: 10-8. See eq.(82) review.)
# "energyorb eps_inv=1.0d-09",
# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

int_wf_begin = "INIT_WF-SECTION"
build_begin = "build"
build_end = "end-build"
int_wf_end = "end-init_wf-section"

initial_state_spec = "   init_state={:d}"

# this applies the dipole moment operator onto the wavefunction at every step
operate_spec = "operate={:s}"

# see https://www.pci.uni-heidelberg.de/tc/usr/mctdh/doc/mctdh/input_docu.html#wfbuild
def _generate_basic_wavefunction(basic_HO_wavepacket, nof_electronic_states, operate_string):
    """Generate basic wavefunction with the same harmonic oscillators for each mode"""
    return "\n".join([
        int_wf_begin,
        build_begin,
        #initial_state_spec.format(nof_electronic_states),
        initial_state_spec.format(1),
        basic_HO_wavepacket,
        build_end,
        operate_spec.format(operate_string),
        int_wf_end,
    ])

def generate_basic_harmonic_oscillator_wavefunction_section(N, A, operate_string):
    """Generate basic wavefunction with the same harmonic oscillators for each mode"""

    basic_HO_wavepacket = "\n".join([
        "-----------------------------------------------------------",
        "#  mode   type  center  moment.  freq.    mass",
        "-----------------------------------------------------------",
        *[
            f"    m{n:<3d}    HO     0.0    0.0      1.0     1.0" for n in N.values()
        ],
        "-----------------------------------------------------------",
    ])

    return _generate_basic_wavefunction(basic_HO_wavepacket, A, operate_string)

basic_HO_wavepacket = "\n".join([
    "-----------------------------------------------------------",
    "#  mode   type  center  moment.  freq.    mass",
    "-----------------------------------------------------------",
    "    m1     HO     0.0    0.0      1.0     1.0",
    "    m2     HO     0.0    0.0      1.0     1.0",
    "    m3     HO     0.0    0.0      1.0     1.0",
    "-----------------------------------------------------------",
])

# see https://www.pci.uni-heidelberg.de/tc/usr/mctdh/doc/mctdh/input_docu.html#wfbuild
initial_wavefunction_section = "\n".join([
    int_wf_begin,
    build_begin,
    initial_state_spec.format(2),
    basic_HO_wavepacket,
    build_end,
    int_wf_end,
])

"""  IMPORTANT!!!!!!
this is extremely necessary and the block relaxation program will not run without it
autoblock generates the initial wavefunction and makes sure the vectors are not linearly dependent.
"""
autoblock_spec = "autoblock"

block_initial_wavefunction_section = "\n".join([
    int_wf_begin,
    build_begin,
    initial_state_spec.format(2),
    basic_HO_wavepacket,
    build_end,
    autoblock_spec,
    int_wf_end,
])
