

def convert_qsize_to_realsize(mode_idx, qsize, frequency, constants):
    """ Convert the reduced dimensionless qsize to the actual rsize in sqrt(amu)*Angs unit """
    omega = frequency[mode_idx]  # frequency of a specific mode
    real_size = qsize
    real_size /= pow(constants['atomic_mass_units_to_<me>'], 0.5)
    real_size /= constants['angstrom_to_bohr-radius']
    real_size /= pow(omega * constants['wavemumbers_to_eh'], 0.5)

    # if debug_print: print(mode_idx, omega, real_size, type(real_size))
    return real_size


def convert_qsize_to_realsize(imode, qsize, amu2me, ang2br, wn2eh, debug_print=True, **kwargs):
    """ Convert the reduced dimensionless qsize to the actual rsize in sqrt(amu)*Angs unit """
    constant_dictionarys = kwargs['constants_dictionary']
    amu2me, ang2br, wn2eh = [constant_dictionarys.get(x) for x in ['amu2me', 'ang2br', 'wn2eh']]

    omega = kwargs['freqcm'][imode]
    rsize = qsize / (pow(amu2me, 0.5) * ang2br * pow(omega * wn2eh, 0.5))
    if debug_print: print(imode, omega, rsize, type(rsize))
    return rsize


def _single_mode_displacement(ndim, real_size, coord_dict, debug_print=True, **kwargs):
    """ You should explain what is happening in more detail """

    ref_q = coord_dict['reference_co_ordinate']
    ref_q = coord_dict['reference_co_ordinate']

    coord_dict['plus'] = []
    # Loop over components (and do what?)
    for i_comp in range(1, ndim + 1):

        # list fashion (can lead to issues?)
        coord_dict['plus'].append(ref_q[i_comp] + real_size * normal_mode_array[i_comp, i_mode])

        # dictionary fashion (like how you had previously)
        coord_dict['plus'][i_comp] = ref_q[i_comp] + real_size * normal_mode_array[i_comp, i_mode]


        displacement = real_size * normal_mode_array[i_comp, i_mode]

        coord_dict['q1_p'] = ref_q[i_comp] - displacement
        coord_dict['q1_m'] = ref_q[i_comp] - displacement
        coord_dict['q1_pp'] = ref_q[i_comp] + 2.0 * displacement
        coord_dict['q1_mm'] = ref_q[i_comp] - 2.0 * displacement

        # double check as always (alternative loop)
        for name, constant in (['q1_p', 'q1_m', ..], [1.0, 1.0, 2.0, 2.0]):
            coord_dict2[name] = ref_q[i_comp] + constant * displacement


        distcoord_plus_x2[i_comp] = ref_q[i_comp] + 2.0 * displacement
        distcoord_minus_x2[i_comp] = ref_q[i_comp] - 2.0 * displacement

        if debug_print:
            print(
                imode, i_comp,
                refcoord[i_comp], normal_mode_array[i_comp, i_mode],
                coord_disp_plus, coord_disp_minus,
                distcoord_plus[i_comp], distcoord_minus[i_comp]
            )
    return


def _double_mode_displacement(ndim, debug_print=True, **kwargs):
    """ You should explain what is happening in more detail """

    # Loop over components (and do what?)

    for icomp in range(1, ndim + 1):
        coord_disp_q1p_q2p = distcoord_plus[icomp] + rsizep * nrmmod[icomp, jmode]
        coord_disp_q1p_q2m = distcoord_plus[icomp] - rsizep * nrmmod[icomp, jmode]
        coord_disp_q1m_q2p = distcoord_minus[icomp] + rsizep * nrmmod[icomp, jmode]
        coord_disp_q1m_q2m = distcoord_minus[icomp] - rsizep * nrmmod[icomp, jmode]

        distcoord_pp[icomp] = coord_disp_pp
        distcoord_pm[icomp] = coord_disp_pm
        distcoord_mp[icomp] = coord_disp_mp
        distcoord_mm[icomp] = coord_disp_mm
    return


def _delete_existint_dist_structure_files(filename_list):
    """ Delete the files because... ?

    We use subprocess because ...?
    We choose to not error on failures ... (explain why)
    """
    for dist_file in filename_list:
        try:
            subprocess.run(['rm', '-f', dist_file])
        except Exception as e:
            print(f"Error deleting {dist_file}: {str(e)}")


def _print_single_mode_distorted_structure(natom, single_mode_list, atmlst, chrglst, **kwargs):
    """
        I'll leave this the same for now.
        Personally I might memory-map this since you're looping over it so much.
        I would want to know the performance cost of appending to files.
        Maybe appending is really cheap, whereas if you used 'w'/'r' flags it might be more expensive.
        Something to think about.
    """
    for iatom in range(1, natom + 1):
        with open('dist_structure_plus', 'a') as f_plus, \
                open('dist_structure_minus', 'a') as f_minus, \
                open('dist_structure_plusx2', 'a') as f_plusx2, \
                open('dist_structure_minusx2', 'a') as f_minusx2:

            string = f"{atmlst[iatom]} {chrglst[iatom]} "
            f_plus.write(string)
            f_minus.write(string)
            f_plusx2.write(string)
            f_minusx2.write(string)

            for ixyz in range(1, 4):
                icomp = (iatom - 1) * 3 + ixyz
                f_plus.write(f"{distcoord_plus[icomp]} ")
                f_minus.write(f"{distcoord_minus[icomp]} ")
                f_plusx2.write(f"{distcoord_plus_x2[icomp]} ")
                f_minusx2.write(f"{distcoord_minus_x2[icomp]} ")

            f_plus.write('\n')
            f_minus.write('\n')
            f_plusx2.write('\n')
            f_minusx2.write('\n')

    # actually the better way to do this is to build 4 strings and append them in one-go
    header = f"{atmlst[iatom]} {chrglst[iatom]} "

    for file_name, coord_name in zip(single_mode_list, ['plus', 'minus', 'plus_x2', 'minus_x2']):

        string = header

        for ixyz in range(1, 4):
            icomp = (iatom - 1) * 3 + ixyz
            string += f"{coord_dict[name][icomp]} "

        string += "\n"

        with open(file_name, 'a') as f_plus:
            f_plus.write(string)

    # make sure you check if you need to add '\n' characters to the string
    # I'm not sure if write adds a '\n' by default every time you call it?


def _print_double_mode_distorted_structure(natom, single_mode_list, **kwargs):
    """ Everything I said for the single mode function applies here """
    for iatom in range(1, natom + 1):
        with open(f'dist_structure_pp', 'a') as f_pp, \
                open(f'dist_structure_pm', 'a') as f_pm, \
                open(f'dist_structure_mp', 'a') as f_mp, \
                open(f'dist_structure_mm', 'a') as f_mm:
            f_pp.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
            f_pm.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
            f_mp.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
            f_mm.write(f"{atmlst[iatom]} {chrglst[iatom]} ")
            for ixyz in range(1, 4):
                icomp = (iatom - 1) * 3 + ixyz
                f_pp.write(f"{distcoord_pp[icomp]} ")
                f_pm.write(f"{distcoord_pm[icomp]} ")
                f_mp.write(f"{distcoord_mp[icomp]} ")
                f_mm.write(f"{distcoord_mm[icomp]} ")
            f_pp.write('\n')
            f_pm.write('\n')
            f_mp.write('\n')
            f_mm.write('\n')
    return


def create_single_mode_diabatization_input_files():
    """ Create input files for diabatization calculations """

    # prefix = f'{model_name}_mode{mode_idx}'
    # suffix = f'{qsize}{suffix}.inp'


    for sign in ['+', '-']:

        p_or_m = {'+': 'plus', '-': 'minus'}[sign]

        for suffix in ['', 'x2']:
            file_name = f'{model_name}_mode{imode}_{sign}{qsize}{suffix}'

            shutil.copy('temp.inp', file_name + '.inp')

            input_filename = f'dist_structure_{p_or_m}{suffix}'
            with open(input_filename, 'r') as dist_file:
                data = dist_file.read()

            with open(file_name + '.inp', 'a') as inp_file:
                inp_file.write(data)
                inp_file.write(' $END')

            # Check if the calculation is done already
            slist = ["grep", "grace", file_name + '.out']
            skwargs = {'stdout':subprocess.PIPE, 'stderr':subprocess.PIPE}
            grace1 = subprocess.run(slist, **skwargs)

            # if not os.path.exists(f'{filnam}_mode{imode}_{sign}{qsize}{suffix}.out'):

            if grace1.returncode != 0:
                print(f"Running calculations for {filnam}_mode{imode}_{sign}{qsize}{suffix}")
                try:
                    # subprocess.run(['./subgam.diab', file_name + '.inp', '4', '0', '1'])
                    os.system("sbatch" + " " + my_subgam(file_name + '.inp', 2, 1, 1))

                except Exception as e:
                    print(f"Error running diabatization calculation: {str(e)}")
            else:
                print(f"{filnam}_mode{imode}_{sign}{qsize}{suffix} is done")
    return


def create_double_mode_diabatization_input_files(**kwargs):
    """ create input files for two distinction modes """
    for displacement1 in ['+', '-']:
        for displacement2 in ['+', '-']:
            if displacement1 == '+' and displacement2 == '+':
                suffix1 = 'pp'
            elif displacement1 == '+' and displacement2 == '-':
                suffix1 = 'pm'
            elif displacement1 == '-' and displacement2 == '+':
                suffix1 = 'mp'
            elif displacement1 == '-' and displacement2 == '-':
                suffix1 = 'mm'

            shutil.copy('temp.inp', f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp')
            with open(f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp', 'a') as inp_file:
                with open(f'dist_structure_{suffix1}', 'r') as dist_file:
                    inp_file.write(dist_file.read())
                inp_file.write(' $END ')

            # Check if the calculation is done already
            output_filename = f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.out'
            grace2 = subprocess.run(["grep", "grace", output_filename])

            if grace2.returncode != 0:
                print(f"Running calculations for {output_filename}!")
                try:
                    #subprocess.run(['./subgam.diab', f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp', '4', '0', '1'])
                    os.system("sbatch" + " " + my_subgam(f'{filnam}_mode{imode}_{displacement1}{qsize}_mode{jmode}_{displacement2}{qsize}.inp', 2, 1, 1))
                except Exception as e:
                    print(f"Error running diabatization calculation: {str(e)}")
            else:
                print(f"{output_filename} is already done.")

    return


def diabatization(modes_included, qsize, freqcm, refcoord, debug_print=True, **kwargs):
    """ Preforms the diabatization routine

    This process involves 2 main steps:
        1: first do single mode displacements
        2: then do bi-linear displacements (two distinction modes)

    Note that variables (atmlst, chrglst, filnam) all need to be floats (for some reason) ... (explain some other details).
    <other explanations of small details go here>
    """

    if True:
        # this is a simplistic way to expand the argument list
        modes_included, freqcm, ndim, refcoord, nrmmod, natom, atmlst, chrglst, filnam, qsize, ha2ev, wn2ev, wn2eh, ang2br, amu2me = args
    elif False:
        # you can also split them up
        modes_included, freqcm, ndim, refcoord, nrmmod, natom, atmlst, chrglst, filnam, qsize = args

        # a bunch of similar looking variables can also be stored in kwargs
        ha2ev, wn2ev, wn2eh, ang2br, amu2me = kwargs['constants']

        # this is better though
        constant_dictionary = kwargs['constants_dictionary']
    else:
        ndim, natom, nrmmod, atmlst, chrglst = kwargs['system_parameters']
        filnam = kwargs['model_name']
        # kwargs['constants']  this contain the numbers: (ha2ev, wn2ev, wn2eh, ang2br, amu2me)


    nof_modes = len(modes_included)  # sometimes this makes code more readable?

    if True:
        # define dictionaries for storing results
        distcoord_plus, distcoord_minus, distcoord_plus_x2, distcoord_minus_x2 = {}, {}, {}, {}
        distcoord_pp, distcoord_pm, distcoord_mp, distcoord_mm = {}, {}, {}, {}
    else:
        # you could just use 1 dictionary of dictionaries
        coord_dict = {
            'plus': {}, 'minus': {},
            'plusx2': {}, 'minusx2': {},
            'pp': {}, 'pm': {}, 'mp': {}, 'mm': {},
        }

    single_mode_filename_list = ['dist_structure_plus', 'dist_structure_minus', 'dist_structure_plusx2', 'dist_structure_minusx2']
    double_mode_filename_list = ['dist_structure_pp', 'dist_structure_pm', 'dist_structure_mp', 'dist_structure_mm']

    # Loop over all considered modes and do + and - displacements
    for k_mode in range(1, nof_modes + 1):
        i_mode = modes_included[k_mode]

        # Get co-ordinates in Angstrom units?
        real_size = convert_qsize_to_realsize(i_mode, qsize, debug_print, **kwargs)

        # --- First we displace individual modes ---
        _single_mode_displacement(ndim, real_size, coord_dict, debug_print=True, **kwargs)

        _delete_existint_dist_structure_files(coord_dict, single_mode_filename_list)

        _print_single_mode_distorted_structure(coord_dict, natom, single_mode_filename_list, **kwargs)

        create_single_mode_diabatization_input_files(**kwargs)

        # --- Next we displace pairs of different modes ---


    def _displace_double_mode()
        # Loop over all considered modes and do + and - displacements
        for k_mode in range(1, nof_modes + 1):
            i_mode = modes_included[k_mode]

            # Get co-ordinates in Angstrom units?
            real_size = convert_qsize_to_realsize(i_mode, qsize, debug_print, **kwargs)

            # 2D distortion to get bilinear vibronic coupling
            for l_mode in range(1, k_mode):
                j_mode = modes_included[l_mode]

                # Get co-ordinates in Angstrom units?
                real_size_p = convert_qsize_to_realsize(j_mode, qsize, debug_print, **kwargs)
                if debug_print: print(i_mode, j_mode, real_size, real_size_p)

                _double_mode_displacement(ndim, real_size, debug_print=True, **kwargs)

                _delete_existint_dist_structure_files(double_mode_filename_list)

                _print_double_mode_distorted_structure(natom, double_mode_filename_list, **kwargs)

            create_double_mode_diabatization_input_files(**kwargs)



    sys.exit()
    return