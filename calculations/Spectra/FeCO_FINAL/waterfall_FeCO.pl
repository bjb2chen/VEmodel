set terminal pngcairo size 1200,800
#set style data line
set nologscale; set xzeroaxis; set samples 5000

set xr [ 2: 12]; set yr [ -20: 50]

set title 'Fe(CO)_5 Vibronic Spectrum' font 'Helvetica,18'
set xlabel 'Energy [eV]'

set key font 'Helvetica,13'
set key Left # left align the text
set key left top # move legend to upper left corner

# see https://stackoverflow.com/questions/17564497/gnuplot-break-y-axis-in-two-parts for plotting disconnected y axis
# or https://stackoverflow.com/questions/56005211/break-y-axis-in-three-parts-in-gnuplot

# see http://www.gnuplot.info/docs/loc11700.html for explanation of %b, %d
month_day = "".strftime("%b_%d", time(0))

mctdh_1 = 'op_FeCO27Q_17st_fullmodes_PBF24_tf180.00_auto_total'
mctdh_2 = 'op_FeCO27Q_17st_fullmodes_PBF26_tf180.00_auto_total'
mctdh_3 = 'op_FeCO27Q_17st_fullmodes_PBF30_tf180.00_auto_total'
mctdh_4 = 'op_FeCO27Q_17st_fullmodes_PBF31_tf180.00_auto_total'

screen1 = 'SCREENED_op_FeCO27Q_17st_PBF24_tf180.00_auto_total'
screen2 = 'SCREENED_op_FeCO27Q_17st_PBF25_tf180.00_auto_total'
screen3 = 'SCREENED_op_FeCO27Q_17st_PBF26_tf180.00_auto_total'
screen4 = 'SCREENED_op_FeCO27Q_17st_PBF27_tf180.00_auto_total'


vanilla = 'op_FeCO27Q_17st_fullmodes_PBF5_UNIFORMSPF3'
van_nosoc = 'op_FeCO27Q_17st_fullmodes_PBF5_UNIFORMSPF3_noSOC'

j16_1 = 'FeCO_Z1_H1_SOC_vibronic_linear_tf180.pl'
j16_2 = 'FeCO_Z2_H1_SOC_vibronic_linear_tf180.pl'
j16_3 = 'FeCO_Z3_H1_SOC_vibronic_linear_tf180.pl'
j16_3_tau = 'FeCO_Z3_H1_SOC_vibronic_linear_tf180_300tau.pl'
j16_4 = 'H0CONSTANTSSOC_op_FeCO27Q_17st_fullmodes_PBF1_tf180.00_auto_total'
j16_5 = 'H0CONSTANTSnoSOC_op_FeCO27Q_17st_fullmodes_PBF2_tf180.00_auto_total'

nosoc1 = 'nosoc_Z1_H1_vibronic_linear_tf180.pl'
nosoc2 = 'nosoc_Z2_H1_vibronic_linear_tf180.pl'
nosoc3 = 'nosoc_Z3_H1_vibronic_linear_tf180.pl'
nosoc_tau1 = 'nosoc_Z3_H1_vibronic_linear_tf180_100tau.pl'
nosoc_tau2 = 'nosoc_Z3_H1_vibronic_linear_tf180_300tau.pl'


Z120fs = 'FeCO_Z1_H2_SOC_vibronic_quadratic_tf20.pl'
Z220fs = 'FeCO_Z2_H2_SOC_vibronic_quadratic_tf20.pl'
Z320fs = 'FeCO_Z3_H2_SOC_vibronic_quadratic_tf20.pl'



# ------------------------------------------------------------------------------------------------------------------------------------------------
set output month_day.'_FeCO_waterfall_composite_collated.png'; set yr [ -1: 13]; set xr [ 0: 13]
plot \
mctdh_1  using 1:2 w lines lw 2 lc rgb 'black'          dashtype 1 title 'MCTDH H1 180fs SOC', \
mctdh_2  using 1:2 w lines lw 1 lc rgb 'red'            dashtype 1 title 'MCTDH H1 180fs no SOC', \
mctdh_3  using 1:2 w lines lw 1 lc rgb 'web-green'      dashtype 1 title 'MCTDH H2 180fs SOC', \
mctdh_4  using 1:2 w lines lw 1 lc rgb 'cyan'           dashtype 1 title 'MCTDH H2 180fs no SOC', \
screen1  using 1:2 w lines lw 1 lc rgb 'pink'           dashtype 1 title 'MCTDH H1 180fs SOC screened', \
screen3  using 1:2 w lines lw 1 lc rgb 'brown'          dashtype 1 title 'MCTDH H1 180fs no SOC screened', \
screen2  using 1:2 w lines lw 1 lc rgb 'slategray'      dashtype 1 title 'MCTDH H2 180fs SOC screened', \
screen4  using 1:2 w lines lw 1 lc rgb 'gold'           dashtype 1 title 'MCTDH H2 180fs no SOC screened', \
j16_1    using 1:3 w lines lw 1 lc rgb 'spring-green'   dashtype 1 title 'Z1 H1 180fs SOC', \
j16_2    using 1:3 w lines lw 1 lc rgb 'magenta'        dashtype 1 title 'Z2 H1 180fs SOC', \
j16_3    using 1:3 w lines lw 1 lc rgb 'salmon'         dashtype 1 title 'Z3 H1 180fs SOC', \
nosoc1   using 1:3 w lines lw 1 lc rgb 'dark-plum'      dashtype 1 title 'Z1 H1 180fs no SOC', \
nosoc2   using 1:3 w lines lw 1 lc rgb 'navy'           dashtype 1 title 'Z2 H1 180fs no SOC', \
nosoc3   using 1:3 w lines lw 1 lc rgb 'khaki'          dashtype 1 title 'Z3 H1 180fs no SOC', \

#vanilla  using 1:2 w lines lw 2 lc rgb 'black'          dashtype 1 title 'MCTDH H1 180fs SOC uniform 3SPF/5PBF', \

# ------------------------------------------------------------------------------------------------------------------------------------------------
set output month_day.'_FeCO_vanilla_vs_orig.png'; set yr [ -1: 17]; set xr [ 0: 13]
plot \
mctdh_1  using 1:2 w lines lw 2 lc rgb 'black'          dashtype 1 title 'MCTDH H1 180fs SOC', \
mctdh_2  using 1:2 w lines lw 1 lc rgb 'red'            dashtype 1 title 'MCTDH H1 180fs no SOC', \
vanilla  using 1:2 w lines lw 1 lc rgb 'green'          dashtype 1 title 'MCTDH H1 180fs SOC uniform 3SPF/5PBF', \
van_nosoc using 1:2 w lines lw 1 lc rgb 'web-blue'      dashtype 1 title 'MCTDH H1 180fs no SOC uniform 3SPF/5PBF', \


# ------------------------------------------------------------------------------------------------------------------------------------------------

set output month_day.'_FeCO_VECCtau.png'; set yr [ -1: 17]; set xr [ 0: 13]
plot \
j16_3        using 1:3 w lines lw 2 lc rgb 'black'      dashtype 1 title 'Z3 H1 180fs SOC tau=40', \
j16_3_tau    using 1:3 w lines lw 1 lc rgb 'purple'     dashtype 1 title 'Z3 H1 180fs SOC tau=300', \
nosoc3       using 1:3 w lines lw 1 lc rgb 'green'      dashtype 1 title 'Z3 H1 180fs no SOC tau=40', \
nosoc_tau1   using 1:3 w lines lw 1 lc rgb 'web-blue'   dashtype 1 title 'Z3 H1 180fs no SOC tau=100', \
nosoc_tau2   using 1:3 w lines lw 1 lc rgb 'red'        dashtype 1 title 'Z3 H1 180fs no SOC tau=300', \


# -----------------------------------------------------------------------------------------------------------------------------------------------
set output month_day.'_FeCO_maximal_linear_collated.png'; set yr [ -1: 17]; set xr [ 0: 13]
plot \
vanilla  using 1:2 w lines lw 2 lc rgb 'black'          dashtype 1 title 'MCTDH H1 180fs SOC uniform 3SPF/5PBF', \
mctdh_1  using 1:2 w lines lw 1 lc rgb 'web-blue'       dashtype 1 title 'MCTDH H1 180fs SOC', \
mctdh_3  using 1:2 w lines lw 1 lc rgb 'web-green'      dashtype 1 title 'MCTDH H2 180fs SOC', \
screen1  using 1:2 w lines lw 1 lc rgb 'pink'           dashtype 1 title 'MCTDH H1 180fs SOC screened', \
screen2  using 1:2 w lines lw 1 lc rgb 'slategray'      dashtype 1 title 'MCTDH H2 180fs SOC screened', \
j16_3    using 1:3 w lines lw 1 lc rgb 'salmon'         dashtype 1 title 'Z3 H1 180fs SOC', \




# -----------------------------------------------------------------------------------------------------------------------------------------------

set title 'Fe(CO)_5 Vibronic Spectrum - VECC Linears' font 'Helvetica,18'
set output month_day.'_FeCO_linears_collated.png'; set yr [ -1: 13]; set xr [ 0: 12]
plot \
vanilla using 1:2 w lines lw 2 lc 'black'      dashtype 1 title 'MCTDH H1 180fs SOC uniform 3SPF/5PBF', \
mctdh_1 using 1:2 w lines lw 1 lc 'red'        dashtype 1 title 'MCTDH H1 180fs SOC', \
mctdh_2 using 1:2 w lines lw 1 lc 'web-blue'   dashtype 1 title 'MCTDH H1 180fs no SOC', \
j16_1   using 1:3 w lines lw 1 lc 'pink'       dashtype 1 title 'Z1 H1 180fs', \
j16_2   using 1:3 w lines lw 1 lc 'cyan'       dashtype 1 title 'Z2 H1 180fs', \
j16_3   using 1:3 w lines lw 1 lc 'dark-green' dashtype 1 title 'Z3 H1 180fs', \
nosoc1  using 1:3 w lines lw 1 lc 'gold'       dashtype 1 title 'Z1 H1 180fs no SOC', \
nosoc2  using 1:3 w lines lw 1 lc 'purple'     dashtype 1 title 'Z2 H1 180fs no SOC', \
nosoc3  using 1:3 w lines lw 1 lc 'brown'      dashtype 1 title 'Z3 H1 180fs no SOC', \


# ------------------------------------------------------------------------------------------------------------------------------------------------

set title 'Fe(CO)_5 Vibronic Spectrum - Constants' font 'Helvetica,18'
set output month_day.'_FeCO_constants_collated.png'; set yr [ -1: 50]; set xr [ 0: 12]
plot \
vanilla using 1:2 w lines lw 2 lc 'black'      dashtype 1 title 'MCTDH H1 180fs SOC uniform 3SPF/5PBF', \
mctdh_1 using 1:2 w lines lw 1 lc 'red'        dashtype 1 title 'MCTDH H1 180fs SOC', \
mctdh_2 using 1:2 w lines lw 1 lc 'web-blue'   dashtype 1 title 'MCTDH H1 180fs no SOC', \
j16_1   using 1:3 w lines lw 1 lc 'pink'       dashtype 1 title 'Z1 H1 180fs', \
j16_2   using 1:3 w lines lw 1 lc 'cyan'       dashtype 1 title 'Z2 H1 180fs', \
j16_3   using 1:3 w lines lw 1 lc 'dark-green' dashtype 1 title 'Z3 H1 180fs', \
nosoc1  using 1:3 w lines lw 1 lc 'gold'       dashtype 1 title 'Z1 H1 180fs no SOC', \
nosoc2  using 1:3 w lines lw 1 lc 'purple'     dashtype 1 title 'Z2 H1 180fs no SOC', \
nosoc3  using 1:3 w lines lw 1 lc 'brown'      dashtype 1 title 'Z3 H1 180fs no SOC', \

# ------------------------------------------------------------------------------------------------------------------------------------------------

set title 'Fe(CO)_5 Vibronic Spectrum Comparison - 3 to 4 eV' font 'Helvetica,18'
set output month_day.'_FeCO_lowerpeak_collated.png'; set yr [ -1: 9]; set xr [ 2.8: 4]
plot \
vanilla using 1:2 w lines lw 2 lc 'black'      dashtype 1 title 'MCTDH H1 180fs SOC uniform 3SPF/5PBF', \
mctdh_1 using 1:2 w lines lw 1 lc 'red'        dashtype 1 title 'MCTDH H1 180fs SOC', \
mctdh_2 using 1:2 w lines lw 1 lc 'web-blue'   dashtype 1 title 'MCTDH H1 180fs no SOC', \
j16_1   using 1:3 w lines lw 1 lc 'pink'       dashtype 1 title 'Z1 H1 180fs', \
j16_2   using 1:3 w lines lw 1 lc 'cyan'       dashtype 1 title 'Z2 H1 180fs', \
j16_3   using 1:3 w lines lw 1 lc 'dark-green' dashtype 1 title 'Z3 H1 180fs', \
nosoc1  using 1:3 w lines lw 1 lc 'gold'       dashtype 1 title 'Z1 H1 180fs no SOC', \
nosoc2  using 1:3 w lines lw 1 lc 'purple'     dashtype 1 title 'Z2 H1 180fs no SOC', \
nosoc3  using 1:3 w lines lw 1 lc 'brown'      dashtype 1 title 'Z3 H1 180fs no SOC', \


# ------------------------------------------------------------------------------------------------------------------------------------------------
set title 'Fe(CO)_5 Vibronic Spectrum - Screened 12 Modes' font 'Helvetica,18'
set output month_day.'_FeCO_screened_only_collated.png'; set yr [ -1: 13]; set xr [ 1: 13]
plot \
mctdh_1  using 1:2 w lines lw 2 lc rgb 'black'           dashtype 1 title 'MCTDH H1 180fs SOC', \
screen1  using 1:2 w lines lw 2 lc rgb 'web-blue'        dashtype 1 title 'MCTDH H1 180fs SOC screened', \
screen3  using 1:2 w lines lw 1 lc rgb 'red'             dashtype 1 title 'MCTDH H1 180fs no SOC screened', \
screen2  using 1:2 w lines lw 1 lc rgb 'sea-green'       dashtype 1 title 'MCTDH H2 180fs SOC screened', \
screen4  using 1:2 w lines lw 1 lc rgb 'purple'	         dashtype 1 title 'MCTDH H2 180fs no SOC screened', \


# ------------------------------------------------------------------------------------------------------------------------------------------------
set title 'Fe(CO)_5 Vibronic Spectrum - MCTDH only' font 'Helvetica,18'
set output month_day.'_FeCO_MCTDHonly.png'; set yr [ -1: 17]; set xr [ 0: 13]
plot \
vanilla using 1:2 w lines lw 2 lc 'black'      dashtype 1 title 'MCTDH H1 180fs SOC uniform 3SPF/5PBF', \
mctdh_1 using 1:2 w lines lw 1 lc 'red'        dashtype 1 title 'MCTDH H1 180fs SOC', \
mctdh_2 using 1:2 w lines lw 1 lc 'web-blue'   dashtype 1 title 'MCTDH H1 180fs no SOC', \
j16_1   using 1:3 w lines lw 1 lc 'pink'       dashtype 1 title 'Z1 H1 180fs', \
j16_2   using 1:3 w lines lw 1 lc 'cyan'       dashtype 1 title 'Z2 H1 180fs', \
j16_3   using 1:3 w lines lw 1 lc 'dark-green' dashtype 1 title 'Z3 H1 180fs', \
nosoc1  using 1:3 w lines lw 1 lc 'gold'       dashtype 1 title 'Z1 H1 180fs no SOC', \
nosoc2  using 1:3 w lines lw 1 lc 'purple'     dashtype 1 title 'Z2 H1 180fs no SOC', \
nosoc3  using 1:3 w lines lw 1 lc 'brown'      dashtype 1 title 'Z3 H1 180fs no SOC', \


# ------------------------------------------------------------------------------------------------------------------------------------------------
set title 'Fe(CO)_5 Vibronic Spectrum - 6 to 12 eV' font 'Helvetica,18'
set output month_day.'_FeCO_higherpeak_collated.png'; set yr [ -1: 15]; set xr [ 6: 12]
plot \
vanilla using 1:2 w lines lw 2 lc 'black'      dashtype 1 title 'MCTDH H1 180fs SOC uniform 3SPF/5PBF', \
mctdh_1 using 1:2 w lines lw 1 lc 'red'        dashtype 1 title 'MCTDH H1 180fs SOC', \
mctdh_2 using 1:2 w lines lw 1 lc 'web-blue'   dashtype 1 title 'MCTDH H1 180fs no SOC', \
j16_1   using 1:3 w lines lw 1 lc 'pink'       dashtype 1 title 'Z1 H1 180fs', \
j16_2   using 1:3 w lines lw 1 lc 'cyan'       dashtype 1 title 'Z2 H1 180fs', \
j16_3   using 1:3 w lines lw 1 lc 'dark-green' dashtype 1 title 'Z3 H1 180fs', \
nosoc1  using 1:3 w lines lw 1 lc 'gold'       dashtype 1 title 'Z1 H1 180fs no SOC', \
nosoc2  using 1:3 w lines lw 1 lc 'purple'     dashtype 1 title 'Z2 H1 180fs no SOC', \
nosoc3  using 1:3 w lines lw 1 lc 'brown'      dashtype 1 title 'Z3 H1 180fs no SOC', \


# ------------------------------------------------------------------------------------------------------------------------------------------------
set title 'Fe(CO)_5 Vibronic Spectrum - VECC 20fs' font 'Helvetica,18'
set output month_day.'_FeCO_20fs.png'; set yr [ -1: 17]; set xr [ 0: 13]
plot \
Z120fs  using 1:3 w lines lw 2 lc rgb 'black'          dashtype 1 title 'Z1 H2 20fs SOC', \
Z220fs  using 1:3 w lines lw 1 lc rgb 'web-blue'       dashtype 1 title 'Z2 H2 20fs SOC', \
Z320fs  using 1:3 w lines lw 1 lc rgb 'red'            dashtype 1 title 'Z3 H2 20fs SOC', \


# ------------------------------------------------------------------------------------------------------------------------------------------------
set title 'Fe(CO)_5 Vibronic Spectrum - Linear SOC' font 'Helvetica,18'
set output month_day.'_FeCO_VECC_lin.png'; set yr [ -1: 13]; set xr [ 0: 12]
plot \
j16_1   using 1:3 w lines lw 1 lc 'pink'       dashtype 1 title 'Z1 H1 180fs', \
j16_2   using 1:3 w lines lw 1 lc 'cyan'       dashtype 1 title 'Z2 H1 180fs', \
j16_3   using 1:3 w lines lw 1 lc 'dark-green' dashtype 1 title 'Z3 H1 180fs', \
nosoc1  using 1:3 w lines lw 1 lc 'gold'       dashtype 1 title 'Z1 H1 180fs no SOC', \
nosoc2  using 1:3 w lines lw 1 lc 'purple'     dashtype 1 title 'Z2 H1 180fs no SOC', \
nosoc3  using 1:3 w lines lw 1 lc 'brown'      dashtype 1 title 'Z3 H1 180fs no SOC', \

set title 'Fe(CO)_5 Vibronic Spectrum - Quadratic SOC' font 'Helvetica,18'
set output month_day.'_FeCO_quad_SOC.png'; set yr [ -1: 13]; set xr [ 0: 12]
plot \


# ------------------------------------------------------------------------------------------------------------------------------------------------
set title 'Fe(CO)_5 Vibronic Spectrum - VECC Z1' font 'Helvetica,18'
set output month_day.'_FeCO_Z1.png'; set yr [ -1: 13]; set xr [ 2.8: 4]
plot \
j16_1 using 1:3 w lines lw 1 lc 'pink' dashtype 1 title 'Z1 H1 180fs', \
nosoc1 using 1:3 w lines lw 1 lc 'gold' dashtype 1 title 'Z1 H1 180fs no SOC', \



# ------------------------------------------------------------------------------------------------------------------------------------------------

set title 'Fe(CO)_5 Vibronic Spectrum - VECC Z2' font 'Helvetica,18'
set output month_day.'_FeCO_Z2.png'; set yr [ -1: 13]; set xr [ 2.8: 4]
plot \

# ------------------------------------------------------------------------------------------------------------------------------------------------
set title 'Fe(CO)_5 Vibronic Spectrum - VECC Z3' font 'Helvetica,18'
set output month_day.'_FeCO_Z3.png'; set yr [ -1: 13]; set xr [ 0: 12]
plot \



