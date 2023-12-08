# vmodels

### Repository for vibronic models construction. Developed in collaboration between Nooijen/Zeng group.

### Protocol procedure. Instructions on how to run the code:
#### - You first need to ensure you have GAMESS with Toby's diabatization (gmcpt.src) module installed.
#### - Configure rungms, runG_diab, and subgam.diab scripts to the appropriate $USER.
#### - Create GAMESS scratch directories if necessary (.gamess_ascii_files and gamess_scr)

### The main diabatization function (lines ~300-450) works by first doing linear coupling (+, -, +x2, -x2), then the subroutine loops over bilinear couplings (plusplus, pm, mp, mm)

#### The following are human input to parameters:
#### modes_excluded and qsize.

Video recording on how to run the diabatization code (assuming you have GAMESS with gmcpt module already): https://youtu.be/nLIn-N0oB-0

## Start of procedure:
### Create molecule in MacMolPlt. Dictate and impose symmetry by going to 'Molecule' -> 'Set Point Group'
### Run a geometry optimization MP2 calculation on the molecular input
### Run python3 dist_allmodes_pm.py gamess.out to extract ref_structure coordinates as a template for distorted structure
### Run a single point gmcpt calculation: you must choose active space with appropriate nmofzc, nmodoc, moact, nelact, mstart(1) values and number of states in the $gmcpt group. Append $VEC group using dat2vec script to the input file.
### Use ./cmo2dmo gamess.dat num_orb ini_orb command to pull out $DMO group and make this the third calculation
### Finally, do ./refdet_from_adiab gamess.out to get the $REFDET group and run this fourth calculation
### Can finally compose temp.inp: cat dmo.out refdet.out >> temp.inp, delete initial $data group and set $data title c1
### Run pythons dist_allmodes_pm.py gamess.out again to perform machine-gun diabatization calculation

[Project Poster Reference](https://github.com/bjb2chen/vmodels/files/10171706/SCP2022_bjc_20685630_White.pdf)

Note to self: I have to get good, ignore the math and just focus on the process. The math will fall out once you are making progress.
