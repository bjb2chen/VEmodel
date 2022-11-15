# This code extracts normal mode eigen vectors
import os
import sys
import math


# this me editing the file, can you see this?
def extract_mode(file_name):
    """
    docstring should be here
    """

    # read in the data
    with open(file_name, 'r') as fp:
        data = fp.readlines()

    nof_lines = len(data)

    # loop over the lines and read in the data
    list_of_line_indicies_mapped_to_line_numbers = []

    # loop over each line
    for i in range(nof_lines):

        # split each line along the spaces
        s_list = data[i].split()

        # if the line is not empty
        if len(s_list) > 0:
            if s_list[0] == "IR" and s_list[1] == "INTENSITY:":
                list_of_line_indicies_mapped_to_line_numbers.append([i, len(s_list)])
            if s_list[0] == "TOTAL" and s_list[1] == "NUMBER" and s_list[2] == "OF" and s_list[3] == "ATOMS":
                natom = int(s_list[-1])

    ncoord = 3*natom
    ndim = ncoord
    md = []

    for i in range(ndim):
        md.append([])

    count = 0

    for i in list_of_line_indicies_mapped_to_line_numbers:
        length_of_line = i[1]

        # this is a (-1) times the length of the i'th line which is IR, INTENSITY
        # this is basically 2 - length_of_line
        # this is "backwards" offset
        # this specifies ONLY the last two/three elements
        not_the_first_two = 2 - length_of_line

        #
        for j in range(ncoord):
            for k in range(length_of_line - 2):

                # index
                data_index = i[0]+j+2

                # split a portion of the line into a list (splitting by spaces)
                list_of_substrings = data[data_index].split()
                # index that list by w.t.f. `ni`
                list_of_last_twoorthree_elements = list_of_substrings[not_the_first_two:]
                # then index whatever those 'right-most' elements are
                something = list_of_last_twoorthree_elements[k]
                # store it
                md[k+count].append(something)

        count += i[1]-2

    fout = open("Nmode.out", "w")

    for i in range(ndim):
        for j in range(ncoord):
            fout.write("%12.8f " % float(md[j][i]))
        fout.write("\n")
    fout.close()


# Calling Main Function
if __name__ == "__main__":
    extract_mode("nh3_ccd_mp2_c3v_gh.out")














