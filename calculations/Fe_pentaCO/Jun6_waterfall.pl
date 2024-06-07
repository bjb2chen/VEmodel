set terminal pngcairo size 1200,800
set output 'JSON-models_waterfall_jun6_MCTDH.png'
#set style data line
set nologscale
set xzeroaxis
set xlabel 'Energy[eV]'
set xr [ 0: 13]
set key font 'Helvetica,13'
set yr [ -10: 100]
set samples 1300
set title 'Fe(CO)_5 Vibronic Spectrum' font 'Helvetica,18'
# Plot command with corrected syntax
plot 'Jun6_Z1_fullmodes_constant_200fs_vibronic_constant_tf200.pl' using 1:3 title 'Z1 full modes constant 200fs' with lines dashtype 2, \
     'op_FeCO_17st_noOffDiag_EH_TDM_PBF5_tf100.00_auto_total' using 1:2 title 'MCTDH no off-diagonal EH and TDM constants 100fs' with lines dashtype 3, \
     'mod_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:3 title 'g1 1e-8 TDM full modes constant 10fs' with lines dashtype 9, \
     '27tdm_fullmodes_Z1_10fs_vibronic_constant_tf10.pl' using 1:3 title 'CONVERSIONFACTOR 27 Z1 10fs g1' with lines dashtype 12 linecolor 'orange', \
     '1e-1tdm_fullmodes_Z1_10fs_vibronic_constant_tf10.pl' using 1:3 w lines lc 'black' title '0.1 TDM 10fs g1', \
     '1e-1tdm_fullmodes_Z1_100fs_vibronic_constant_tf100.pl' using 1:3 w lines lc 'gray' title '0.1 TDM 100fs g1', \
     'Z1_lin_20fs_p1T_vibronic_linear_tf20.pl' using 1:3 w lines lc 'pink' title 'Z1 linear 20fs p1T',


     # 'mult2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:($3*25) title 'Z1 multiply by 2 10fs' with lines dashtype 4 linecolor 'orange', 

     # 'div2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:2 title 'g0 Division by 2 on electronic energies 10fs' with lines dashtype 4 linecolor 'pink', \
     # 'div2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:3 title 'g1 Division by 2 on electronic energies 10fs' with lines dashtype 5 linecolor 'blue', \
     # 'div2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:4 title 'g2 Division by 2 on electronic energies 10fs' with lines dashtype 6 linecolor 'red', \
     # 'div2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:5 title 'g5 Division by 2 on electronic energies 10fs' with lines dashtype 7 linecolor 'black',

     # 'div2_Z1_fullmodes_constant_50fs_vibronic_constant_tf50.pl' using 1:2 title 'g0 Division by 2 on electronic energies 50fs' with lines dashtype 8 linecolor 'pink', \
     # 'div2_Z1_fullmodes_constant_50fs_vibronic_constant_tf50.pl' using 1:3 title 'g1 Division by 2 on electronic energies 50fs' with lines dashtype 9 linecolor 'blue', \
     # 'div2_Z1_fullmodes_constant_50fs_vibronic_constant_tf50.pl' using 1:4 title 'g2 Division by 2 on electronic energies 50fs' with lines dashtype 10 linecolor 'red', \
     # 'div2_Z1_fullmodes_constant_50fs_vibronic_constant_tf50.pl' using 1:5 title 'g5 Division by 2 on electronic energies 50fs' with lines dashtype 11 linecolor 'black', \

     # 'mod_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:2 title 'g0 1e-8 TDM full modes constant 10fs' with lines dashtype 8, \
     # 'mod_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:3 title 'g1 1e-8 TDM full modes constant 10fs' with lines dashtype 9, \
     # 'mod_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:4 title 'g2 1e-8 TDM full modes constant 10fs' with lines dashtype 10, \
     # 'mod_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:5 title 'g5 1e-8 TDM full modes constant 10fs' with lines dashtype 11, 


     # 'div2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:2 title 'g0 Division by 2 on electronic energies 10fs' with lines dashtype 4 linecolor 'pink', \
     # 'div2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:3 title 'g1 Division by 2 on electronic energies 10fs' with lines dashtype 5 linecolor 'blue', \
     # 'div2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:4 title 'g2 Division by 2 on electronic energies 10fs' with lines dashtype 6 linecolor 'red', \
     # 'div2_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:5 title 'g5 Division by 2 on electronic energies 10fs' with lines dashtype 7 linecolor 'black',

     # 'mod_Z1_fullmodes_constant_10fs_vibronic_constant_tf10.pl' using 1:3 title '1e-8 TDM full modes constant 10fs' with lines dashtype 8, \
     # 'mod_Z1_fullmodes_constant_50fs_vibronic_constant_tf50.pl' using 1:3 title '1e-8 TDM full modes constant 50fs' with lines dashtype 11, \
     # 'mod_Z1_fullmodes_constant_100fs_vibronic_constant_tf100.pl' using 1:3 title '1e-8 TDM full modes constant 100fs' with lines lw 1.5 dashtype 7, \
     # 'op_FeCO3Q_17st_PBF5_tf100.00_auto_total' using 1:2 title 'MCTDH 1D TDM constants 100fs' with lines dashtype 12 linecolor rgb 'blue', \
     # 'op_FeCO_17st_XYZ_PBF5_tf100.00_auto_total' using 1:2 title 'MCTDH 3D TDM constants 100fs' with lines dashtype 5, \
     # 'op_FeCO_17st_noEHoffdiag_PBF5_tf100.00_auto_total' using 1:2 title 'MCTDH no EH off-diagonal 1D TDM constants 100fs' with lines dashtype 14, \
     # 'op_FeCO_17st_XYZ_noEHoffdiag_PBF5_tf100.00_auto_total' using 1:2 title 'MCTDH no EH off-diagonal 3D TDM constants 100fs' with lines dashtype 15, \

     #'Jun6_Z1_12modes_constant_vibronic_constant_tf200.pl' using 1:3 title 'Z1 screened modes constant 200fs' with lines lw 2 dashtype 3, \
     #'Jun6_Z2_fullmodes_constant_200fs_vibronic_constant_tf200.pl' using 1:3 title 'Z2 full modes constant 200fs' with lines dashtype 4, \
     #'Jun6_Z2_fullmodes_linear_200fs_vibronic_linear_tf200.pl' using 1:3 title 'Z2 full modes linear 200fs' with lines dashtype 5, \
     #'Jun6_Z3_fullmodes_quadratic_10fs_vibronic_quadratic_tf10.pl' using 1:3 title 'Z3 full modes quadratic 10fs' with lines dashtype 6, \
     #'mod_Z1_fullmodes_constant_100fs_vibronic_constant_tf100_TAU80.pl' using 1:($3*-1) title '**TAU=80** Z1 full modes constant 100fs' with lines dashtype 9, \
     #'Jun6_Z3_fullmodes_quadratic_10fs_vibronic_quadratic_tf10_TAU20.pl' using 1:($3*-1) title '**TAU=20** Z3 full modes quadratic 10fs' with lines dashtype 10,