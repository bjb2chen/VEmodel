run_begin = "RUN-SECTION"
run_end = "end-run-section"

basic_relaxation = "relaxation"
basic_relaxation_scheme = "tout = 5.0 tfinal = {tfinal:.2f}"

# this specifies the precision of the wavefunction when written to the output
wavefunction_single_precision = "psi=single"
wavefunction_double_precision = "psi=double"

energy_shift = 0.0
energy_units_wavenumbers = f"rlxunit=cm-1, {energy_shift}"
energy_units_electron_volts = f"rlxunit=eV, {energy_shift}"

run_section_relaxation = "\n".join([
    run_begin,
    "title = ground state relaxation of water",
    "name = h2o_nogs",
    basic_relaxation,
    basic_relaxation_scheme,
    energy_units_wavenumbers,
    wavefunction_double_precision,
    run_end,
])

block_relaxation_olsen = "relaxation=0,olsen precon=50"
# block_relaxation_scheme = "tout = 5.0 tfinal = {tfinal:.2f}"

block_run_section_relaxation = "\n".join([
    run_begin,
    "title = block improved water relaxation",
    "name = h2o_nogs",
    block_relaxation_olsen,
    basic_relaxation_scheme,
    energy_units_wavenumbers,
    wavefunction_double_precision,
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

spf_basis_section = "\n".join([
    spf_begin,
    "multi-set",
    "      v01      =  {other_spf:d}, {max_spf:d}, {other_spf:d}",
    "      v02      =  {other_spf:d}, {max_spf:d}, {other_spf:d}",
    "      v03      =  {other_spf:d}, {max_spf:d}, {other_spf:d}",
    spf_end,
])

single_set_spf_basis_section = "\n".join([
    spf_begin,
    "single-set",
    "      v01      =  {n_BF:d}",
    "      v02      =  {n_BF:d}",
    "      v03      =  {n_BF:d}",
    spf_end,
])

# For block improved relaxation, single-set  must be given. A block improved relaxation requires the packet keyword and a DAV, RDAV or RRDAV "integrator".
block_spf_basis_section = "\n".join([
    spf_begin,
    "packet=10, single-set",
    "      v01      =  {n_BF:d}",
    "      v02      =  {n_BF:d}",
    "      v03      =  {n_BF:d}",
    spf_end,
])


# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------


pbs_begin = "PRIMITIVE-BASIS-SECTION"
pbs_end = "end-primitive-basis-section"
n_electronic_states = 3
electronic_basis = f"    el     el     {n_electronic_states}"

# this is the basis specification where the first element represents...
# the second element is the frequency and the third element is the mass
ho_spec = "0.0   1.0   1.0"
# this specification allows us to define the HO by the first and last grid points
ho_spec_bounded = "xi-xf   0.0   1.0"


primitive_basis_section = "\n".join([
    pbs_begin,
    "    v01    HO     {n_BF:d}   " + ho_spec,
    "    v02    HO     {n_BF:d}   " + ho_spec,
    "    v03    HO     {n_BF:d}   " + ho_spec,
    electronic_basis,
    pbs_end,
])


block_primitive_basis_section = "\n".join([
    pbs_begin,
    "    v01    HO     {n_BF:d}   " + ho_spec,
    "    v02    HO     {n_BF:d}   " + ho_spec,
    "    v03    HO     {n_BF:d}   " + ho_spec,
    electronic_basis,
    pbs_end,
])


# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------


int_begin = "INTEGRATOR-SECTION"
int_end = "end-integrator-section"

# mean fields are changed based on the errors of both the MCTDH coefficients AND the single particle functions
# default for propagation and simple relaxation
var = "var"

# mean fields are changed based on the errors of ONLY the single particle functions
# default for block improved relaxation
varphi = "varphi"

cmf_accuracy = "1.0d-06"
cmf_initial_time_interval = 0.5
constant_mean_field = "CMF/{var_spec:} = {R:}, {R1:}".format(
    var_spec=varphi,
    R=cmf_initial_time_interval,
    R1=cmf_accuracy
)

# Bulirsch-Stoer extrapolation integrator
bs_initial_step_size = "2.5d-06"
bs_accuracy = "1.0d-06"
bs_maximal_order = 7
bs_integrator = "BS/{bs_spec} = {I:}, {R:}, {R1:}".format(
    bs_spec='spf',  # could also be `all` or `A`
    I=bs_maximal_order,
    R=bs_accuracy,
    R1=bs_initial_step_size,
)

# SIL integrator
sil_accuracy = "1.0d-06"
sil_maximal_order = 5
sil_integrator = "SIL/{sil_spec} = {I:}, {R:}".format(
    sil_spec='A',  # could also be `all` or `A`
    I=sil_maximal_order,
    R=sil_accuracy,
)


# Runge-Kutta integrator (order 8, can also replace 8 with 5)
rk_initial_step_size = "0.1"
rk_accuracy = "1.0d-8"
rk_integrator = "RK8/{rk_spec} = {R:}, {R1:}".format(
    rk_spec='spf',  # could also be `all` or `A`
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
dav_integrator = "rrDAV = {I:}, {R:}".format(
    I=dav_maximal_order,
    R=dav_accuracy,
)

""" The integrators may only be combined in certain ways as described below!!!
For VMF calculations, only the BS, ABM or RK5/8 integrator may be used for the differential equations. Default is ABM. For VMF calculations the integrator must carry the extension S=all (or no extension at all), i.e. there is only one integrator within the VMF scheme. For CMF calculations, the following combinations of integrators are possible: ABM/spf + SIL/A, BS/spf + SIL/A, RKx/spf + SIL/A, BS/all, ABM/all, RKx/all.
"""


integrator_section = "\n".join([
    int_begin,
    constant_mean_field,
    rk_integrator,
    sil_integrator,
    int_end,
])


block_integrator_section = "\n".join([
    int_begin,
    constant_mean_field,
    rk_integrator,
    dav_integrator,
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

initial_state = 2
initial_state_spec = f"   init_state={initial_state:d}"

"""  IMPORTANT!!!!!!
this is extremely necessary and the block relaxation program will not run without it
autoblock generates the initial wavefunction and makes sure the vectors are not linearly dependent.
"""
autoblock_spec = "autoblock"

basic_HO_wavepacket = "\n".join([
    "-----------------------------------------------------------",
    "#  mode   type  center  moment.  freq.    mass",
    "-----------------------------------------------------------",
    "    v01    HO     0.0    0.0      1.0     1.0",
    "    v02    HO     0.0    0.0      1.0     1.0",
    "    v03    HO     0.0    0.0      1.0     1.0",
    "-----------------------------------------------------------",
])

# see https://www.pci.uni-heidelberg.de/tc/usr/mctdh/doc/mctdh/input_docu.html#wfbuild
initial_wavefunction_section = "\n".join([
    int_wf_begin,
    build_begin,
    initial_state_spec,
    basic_HO_wavepacket,
    build_end,
    int_wf_end,
])

block_initial_wavefunction_section = "\n".join([
    int_wf_begin,
    build_begin,
    initial_state_spec,
    basic_HO_wavepacket,
    build_end,
    autoblock_spec,
    int_wf_end,
])
