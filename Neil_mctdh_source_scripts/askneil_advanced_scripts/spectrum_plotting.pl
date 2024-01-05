set terminal png size 1200,800
set output '/work/ngraymon/projects/mctdh_compare/t_amplitudes_aug2020/ch2o/ch2o_vibronic_linear_SPF1_PBF30_tf10.00/spectrum_ch2o_vibronic_linear_3000_SPF1_PBF30_10fs_40tau.png'
set style data line
set nologscale
set xzeroaxis
set xlabel 'Energy[eV]'
set xr [ 21: 11]
set yr [ -1.5: 60]
set title 'ch_{2}o vibronic linear spectrum, n-cos = 1, tau: 40.0 1, 10fs'
plot             '/work/ngraymon/projects/mctdh_compare/t_amplitudes_aug2020/ch2o/ch2o_vibronic_linear_SPF1_PBF30_tf10.00/mctdh_spectrum_ch2o_vibronic_linear_SPF1_PBF30_tf10.00.pl' using 1:3 lw 2 lc 'black' title 'MCTDH',        