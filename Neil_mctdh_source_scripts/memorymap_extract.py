import os
import mmap

import helper
from helper import StringNotFoundError


def find_byte_begin_and_end(path, memmap, begin_string, end_string_list):
    """ attempts to find the index(in bytes) of the start and end of the important region, which are denoted by the begin_string and end_string
    """
    try:
        if True:
            # find the occurrence of `begin` closest to the end of the file
            memmap.seek(0)  # start looking from the end of the file
            begin = helper.rfind_string_in_file(memmap, path, begin_string)
        else:
            # find the occurrence of `begin` closest to the start of the file
            memmap.seek(0)  # start looking from the begining of the file
            begin = helper.rfind_string_in_file(memmap, path, begin_string)

        for end_string in end_string_list:
            try:
                end = helper.rfind_string_in_file(memmap, path, end_string)
            except StringNotFoundError as e:
                continue
            else:
                break
        else:
            raise Exception("File is improperly formatted, couldn't find end string of {begin:s}")

    except StringNotFoundError as err:
        print(f"Couldn't find {begin_string}")
        raise err

    return begin, end


def extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=1):
    """abstract function for extracting a list of strings for each line in a specific section of the *.op file corresponding to the memmap object"""

    # if you need multiple failback end_string's
    # begin_string = headers[header_index]
    # end_string_list = [headers[header_index + 1:]]

    end_string_list = [end_string, ]
    try:
        begin, end = find_byte_begin_and_end(path, memmap, begin_string, end_string_list)
    except StringNotFoundError as err:
        print(err)
        return list()

    assert end > begin, f"{begin=} < {end=} means you didn't pick up the correct delimiters"

    # go to the beginning of that region
    memmap.seek(begin)

    # skip headers
    helper.readlines(memmap, nof_line_skip)  # you may need to change this

    # get the important bytes
    byteData = memmap.read(end - memmap.tell())

    # decode it
    stringData = byteData.decode(encoding="utf-8")
    # print(stringData); breakpoint()

    # clean it up
    stringData = stringData.strip().replace('=', '').replace(', ev', '')

    # return the list of strings (not including empty strings and comments)
    return [line.split() for line in stringData.splitlines() if line != "" and 'STATE' in line]


def _example_processing_function(path, memmap):

    # if processing block 1
    begin_string, end_string = "block_1_begin", "block_1_end"
    # here is where you call extract string
    lines = extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=1)

    def process_block_1(): return
    b1_out = process_block_1(lines)

    begin_string, end_string = "block_2_begin", "block_2_end"
    # here is where you call extract string
    lines = extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=3)

    def process_block_2(): return
    b2_out = process_block_2(lines)

    stuff = [b1_out, b2_out, ]  # whatever you might need to return?

    return stuff


def extract_from_file(path, processing_function, *args, **kwargs):

    helper.verify_file_exists(path)

    # store the energy offsets, and all the coupling terms
    with open(path, "r+b") as source_file:

        """ Turns out the memory protection arguments for mmap are different
        for windows/unix,  so this is how we handle that
        """
        running_on_windows = bool(os.name == 'nt')
        read_write_prot_kwarg = {}
        if running_on_windows:
            read_write_prot_kwarg["access"] = mmap.ACCESS_READ
        else:
            read_write_prot_kwarg["prot"] = mmap.PROT_READ

        # access the file using memory map for efficiency
        with mmap.mmap(source_file.fileno(), 0, **read_write_prot_kwarg) as mm:
            output = processing_function(path, mm, *args, **kwargs)

    return output


# ---------------------------------------------------------------------------------------
def extract_diabatic_energies_from_refG(path):
    import numpy as np

    raise Exception('Not fully implemented')
    def _extract_energy_from_refG(path, memmap):
        """
        Expects lines of the format (where N = 6):
             --- DIABATIC ENERGIES (DIAGONAL ELEMENT) ---       HARTREE            EV
             STATE #  1'S GMC-PT-LEVEL DIABATIC ENERGY=     -55.780701782       0.000000000
             <...>
             STATE #  6'S GMC-PT-LEVEL DIABATIC ENERGY=     -55.898425032      -3.203412786

             --- DIABATIC COUPLINGS (OFF DIAGONAL ELEMENTS)---  HATREE            EV
             STATE #  1 &  2'S GMC-PT-LEVEL COUPLING  =      -0.018622461      -0.506742966
             <...>
             STATE #  5 &  6'S GMC-PT-LEVEL COUPLING  =      -0.018622461      -0.506742966
        """

        begin_string = " - DM DIABT PT HAMILTONIAN MATRIX ELEMENTS -"
        end_string = "ROTATION TO DIABATIC STATES"
        lines = extract_string_list(path, memmap, begin_string, end_string, nof_line_skip=4)

        print("Ingested:")
        for i, l in enumerate(lines): print(f"Line {i:02d}: ", l)

        hartree_list = [float(l[6]) for l in lines]
        hartree_array = np.array(hartree_list)
        if False: print("Hartree array\n", hartree_array)

        eV_list = [float(l[7]) for l in lines]
        eV_array = np.array(eV_list)
        if False: print("eV array\n", eV_array)

        return eV_array

    """ use memory map to extract data  """
    string = extract_from_file(path, _extract_energy_from_refG)
    return string

# ---------------------------------------------------------------------------------------