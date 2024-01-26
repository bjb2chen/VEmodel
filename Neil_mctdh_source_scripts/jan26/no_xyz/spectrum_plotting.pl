set terminal png size 1200,800
set output 'output.png'
set style data line
set nologscale
set xzeroaxis
set xlabel 'Energy[eV]'
set xr [ 21: 11]
set yr [ -1.5: 60]
set title 'op H2O3Q 3st spectrum, n-cos = 1, tau: 40.0 1, 10fs'
plot 'auto_h2o_FC_linear_PBF100_tf50.00.txt.pl' using 1:3 lw 2 lc 'black' title 'mctdh g1', 'auto_h2o_FC_linear_PBF100_tf50.00.txt.pl' using 1:4 lw 2 lc 'red' title 'mctdh g2', 'auto_h2o_FC_linear_PBF100_tf50.00.txt.pl' using 1:2 lw 2 lc 'orange' title 'mctdh g0'
 
