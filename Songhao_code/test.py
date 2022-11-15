from distort_structure import *
import numpy as np

def main():
    file_name = "step1_h2o_cct_mp2_c2v_gh.log"
    molecule_name = "h2o"
    # file_name = "nh3_ccd_mp2_c3v_gh.out"
    # initialize the python object
    data = dist_structure(file_name, molecule_name)
    # readin vibrational data from GAMESS output file
    data.read_vibration()
    # scan over each normal mode and generate a 1d grid of distorted structure
    data.dis_structure_1d(ref_file="ref_structure", grid=np.linspace(0,0.8,9), temp_file="temp.inp", diab_info="diab_info")
    # scan over each normal mode and generate a 2d grid of distorted structure
    # data.dis_structure_2d(ref_file="ref_structure", grid=np.linspace(0,0.8,9), temp_file="temp.inp", diab_info="diab_info")
    return


if __name__ == "__main__":
    main()
