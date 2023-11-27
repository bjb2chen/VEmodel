# vmodels

### Repository for vibronic models construction. Deveolped in collaboration between Nooijen/Zeng group.

### Protocol procedure. Instructions on how to run the code:
#### - You first need to ensure you have GAMESS with Toby's diabatization (gmcpt.src) module installed.
#### - Configure rungms, runG_diab, and subgam.diab scripts to the appropriate $USER.
#### - Create GAMESS scratch directories if necessary (.gamess_ascii_files and gamess_scr)
#### - Have an input geometry written in ref_structure
#### - Have temp.inp written

### The main diabatization function (lines ~300-450) works by first doing linear coupling (+, -, +x2, -x2), then the subroutine loops over bilinear couplings (plusplus, pm, mp, mm)

#### The following are human input to parameters:
#### modes_excluded and qsize.

Video recording of diabatization code being run: https://youtu.be/nLIn-N0oB-0

[Project Poster Reference](https://github.com/bjb2chen/vmodels/files/10171706/SCP2022_bjc_20685630_White.pdf)
