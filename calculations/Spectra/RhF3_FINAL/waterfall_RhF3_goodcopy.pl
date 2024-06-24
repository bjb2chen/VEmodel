set terminal pngcairo size 1200,800
#set style data line
set nologscale; set xzeroaxis; set samples 5000

set xr [ 2: 12]; set yr [ -20: 50]

set title 'RhF_3 Vibronic Spectrum' font 'Helvetica,18'
set xlabel 'Energy[eV]'

set key font 'Helvetica,13'
set key Left # left align the text
set key right top # move legend to upper left corner

# see https://stackoverflow.com/questions/17564497/gnuplot-break-y-axis-in-two-parts for plotting disconnected y axis
# or https://stackoverflow.com/questions/56005211/break-y-axis-in-three-parts-in-gnuplot

# see http://www.gnuplot.info/docs/loc11700.html for explanation of %b, %d
month_day = "".strftime("%b_%d", time(0))

###############################################################
# THE FOLLOWING ARE NOW WITH JUN 13 SOC CORRECTIONS!!!
###############################################################

f01 = 'op_RhF36Q_15st_PBF40_tf180.00_auto_total'   # PBF40/180fs/40tau quad SOC
f02 = 'op_RhF36Q_15st_PBF40_tf182.00_auto_total'   # PBF40/180fs/40tau lin SOC
f03 = 'op_RhF36Q_15st_PBF40_tf184.00_auto_total'   # PBF40/180fs/40tau quad no SOC
f04 = 'op_RhF36Q_15st_PBF40_tf186.00_auto_total'   # PBF40/180fs/40tau lin no SOC
f05 = 'op_RhF36Q_15st_PBF5_tf180.00_auto_total'   # vanilla PBF5/SPF3/180fs/40tau quad SOC 
f06 = 'op_RhF36Q_15st_PBF5_tf182.00_auto_total'   # vanilla PBF5/SPF3/180fs/40tau lin SOC 
f07 = 'op_RhF36Q_15st_PBF5_tf184.00_auto_total'   # vanilla PBF5/SPF3/180fs/40tau quad no SOC 
f08 = 'op_RhF36Q_15st_PBF5_tf186.00_auto_total'   # vanilla PBF5/SPF3/180fs/40tau lin no SOC
f09 = 'op_RhF36Q_15st_PBF40_tf1000.00_auto_total'   # PBF40/1000fs/40tau quad SOC
f10 = 'op_RhF36Q_15st_PBF40_tf1002.00_auto_total'   # PBF40/1000fs/40tau lin SOC
f11 = 'op_RhF36Q_15st_PBF40_tf1004.00_auto_total'   # PBF40/1000fs/40tau quad no SOC
f12 = 'op_RhF36Q_15st_PBF40_tf1006.00_auto_total'   # PBF40/1000fs/40tau lin no SOC
f13 = 'op_RhF36Q_15st_PBF5_tf1000.00_auto_total'   # vanilla PBF5/SPF3/1000fs/40tau quad SOC
f14 = 'op_RhF36Q_15st_PBF5_tf1002.00_auto_total'   # vanilla PBF5/SPF3/1000fs/40tau lin SOC
f15 = 'op_RhF36Q_15st_PBF5_tf1004.00_auto_total'   # vanilla PBF5/SPF3/1000fs/40tau quad no SOC
f16 = 'op_RhF36Q_15st_PBF5_tf1006.00_auto_total'   # vanilla PBF5/SPF3/1000fs/40tau lin no SOC
f17 = ''                                           # extra space for now
f18 = ''
f19 = ''
f20 = ''

# VECC runs
f21 = 'RhF3_Z1_H1_SOC_vibronic_linear_tf180.pl'    # SOC runs, they are 40 tau
f22 = 'RhF3_Z2_H1_SOC_vibronic_linear_tf180.pl'    #    "
f23 = 'RhF3_Z3_H1_SOC_vibronic_linear_tf180.pl'
f24 = 'RhF3_Z1_H2_SOC_vibronic_quadratic_tf180.pl'
f25 = 'RhF3_Z2_H2_SOC_vibronic_quadratic_tf180.pl'
f26 = 'RhF3_Z3_H2_SOC_vibronic_quadratic_tf180.pl'
f27 = 'RhF3_Z1_H1_vibronic_linear_tf180.pl'        # non SOC runs, they are 40 tau
f28 = 'RhF3_Z2_H1_vibronic_linear_tf180.pl'        #    "
f29 = 'RhF3_Z3_H1_vibronic_linear_tf180.pl'
f30 = 'RhF3_Z1_H2_vibronic_quadratic_tf180.pl'
f31 = 'RhF3_Z2_H2_vibronic_quadratic_tf180.pl'
f32 = 'RhF3_Z3_H2_vibronic_quadratic_tf180.pl'
f33 = ''                                           # extra space, the VECC long prop are no good.
f34 = 'RhF3_Z1_H1_SOC_vibronic_linear_tf50.pl'     # RhF3 short prop
f35 = 'RhF3_Z1_H1_SOC_vibronic_linear_tf100.pl' 
f36 = ''


#---------------------------------------------------------------------------------------------------

set title 'RhF_3 Vibronic Spectrum Composite - 180 fs ; 40 tau' font 'Helvetica,18'
set output month_day.'_RhF3_waterfall_composite.png'; set yr [ -1: 100]; set xr [ -1.5: 0]
plot \
f02 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H1 SOC', \
f01 using 1:2 w lines lw 0.75 lc 'web-blue'   dashtype 1 title 'MCTDH H2 SOC', \
f04 using 1:2 w lines lw 0.75 lc 'red'        dashtype 1 title 'MCTDH H1', \
f03 using 1:2 w lines lw 0.75 lc 'dark-green' dashtype 1 title 'MCTDH H2', \
f21 using 1:3 w lines lw 0.75 lc 'purple'     dashtype 1 title 'VECC    H1 SOC Z1', \
f22 using 1:3 w lines lw 0.75 lc 'pink'       dashtype 1 title 'VECC    H1 SOC Z2', \
f23 using 1:3 w lines lw 0.75 lc 'brown'      dashtype 1 title 'VECC    H1 SOC Z3', \
f24 using 1:3 w lines lw 0.75 lc 'steelblue'  dashtype 1 title 'VECC    H2 SOC Z1', \
f25 using 1:3 w lines lw 0.75 lc 'coral'      dashtype 1 title 'VECC    H2 SOC Z2', \
f26 using 1:3 w lines lw 0.75 lc 'seagreen'   dashtype 1 title 'VECC    H2 SOC Z3', \
f27 using 1:3 w lines lw 0.75 lc 'goldrenrod' dashtype 1 title 'VECC    H1     Z1', \
f28 using 1:3 w lines lw 0.75 lc 'navy'       dashtype 1 title 'VECC    H1     Z2', \
f29 using 1:3 w lines lw 0.75 lc 'violet'     dashtype 1 title 'VECC    H1     Z3', \
f30 using 1:3 w lines lw 0.75 lc 'turquoise'  dashtype 1 title 'VECC    H2     Z1', \
f31 using 1:3 w lines lw 0.75 lc 'dark-plum'  dashtype 1 title 'VECC    H2     Z2', \
f32 using 1:3 w lines lw 0.75 lc 'dark-pink'  dashtype 1 title 'VECC    H2     Z3', \

#---------------------------------------------------------------------------------------------------

set title 'RhF_3 Vibronic Spectrum - MCTDH SOC Effect Comparison' font 'Helvetica,18'
set output month_day.'_RhF3_MCTDH_SOC.png'; set yr [ -1: 100]; set xr [ -1.0: -0.4]
plot \
f02 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H1 SOC', \
f01 using 1:2 w lines lw 1.00 lc 'web-blue'   dashtype 1 title 'MCTDH H2 SOC', \
f04 using 1:2 w lines lw 1.00 lc 'violet'     dashtype 1 title 'MCTDH H1', \
f03 using 1:2 w lines lw 1.00 lc 'orange-red' dashtype 1 title 'MCTDH H2', \

#---------------------------------------------------------------------------------------------------

set title 'RhF_3 Vibronic Spectrum - VECC SOC Effect Comparison' font 'Helvetica,18'
set output month_day.'_RhF3_VECC_SOC.png'; set yr [ -1: 100]; set xr [ -1.0: -0.4]
plot \
f23 using 1:3 w lines lw 2.25 lc 'black'      dashtype 1 title 'VECC    H1 SOC Z3', \
f26 using 1:3 w lines lw 1.00 lc 'web-blue'   dashtype 1 title 'VECC    H2 SOC Z3', \
f29 using 1:3 w lines lw 1.00 lc 'violet'     dashtype 1 title 'VECC    H1     Z3', \
f32 using 1:3 w lines lw 1.00 lc 'orange-red' dashtype 1 title 'VECC    H2     Z3', \


#---------------------------------------------------------------------------------------------------
set title 'RhF_3 Vibronic Spectrum - Linear Model Comparison' font 'Helvetica,18'
set output month_day.'_RhF3_Linears.png'; set yr [ -1: 100]; set xr [ -1.0: -0.4]
plot \
f02 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H1 SOC', \
f04 using 1:2 w lines lw 1.00 lc 'web-blue'   dashtype 1 title 'MCTDH H1', \
f23 using 1:3 w lines lw 1.00 lc 'violet'     dashtype 1 title 'VECC    H1 SOC Z3', \
f29 using 1:3 w lines lw 1.00 lc 'orange-red' dashtype 1 title 'VECC    H1     Z3', \


#---------------------------------------------------------------------------------------------------
set title 'RhF_3 Vibronic Spectrum - Quadratic Model Comparison' font 'Helvetica,18'
set output month_day.'_RhF3_Quadratics.png'; set yr [ -1: 100]; set xr [ -1.0: -0.4]
plot \
f01 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H2 SOC', \
f03 using 1:2 w lines lw 1.00 lc 'web-blue'   dashtype 1 title 'MCTDH H2', \
f26 using 1:3 w lines lw 1.00 lc 'violet'     dashtype 1 title 'VECC    H2 SOC Z3', \
f32 using 1:3 w lines lw 1.00 lc 'orange-red' dashtype 1 title 'VECC    H2     Z3', \


#---------------------------------------------------------------------------------------------------

set title 'RhF_3 Vibronic Spectrum Composite - 1000 fs ; 200 tau' font 'Helvetica,18'
set output month_day.'_RhF3_longprop_200tau.png'; set yr [ -1: 100]; set xr [ -1.0: -0.4]
plot \
f10 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H1 SOC PBF40', \
f09 using 1:2 w lines lw 1.00 lc 'web-blue'   dashtype 1 title 'MCTDH H2 SOC PBF40', \
f12 using 1:2 w lines lw 1.00 lc 'violet'     dashtype 1 title 'MCTDH H1 PBF40', \
f11 using 1:2 w lines lw 1.00 lc 'dark-green' dashtype 1 title 'MCTDH H2 PBF40', \
f14 using 1:2 w lines lw 1.00 lc 'red'        dashtype 1 title 'MCTDH H1 SOC PBF5', \
f13 using 1:2 w lines lw 1.00 lc 'turquoise'  dashtype 1 title 'MCTDH H2 SOC PBF5', \
f16 using 1:2 w lines lw 1.00 lc 'dark-plum'  dashtype 1 title 'MCTDH H1 PBF5', \
f15 using 1:2 w lines lw 1.00 lc 'chartreuse' dashtype 1 title 'MCTDH H2 PBF5', \

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

set title 'RhF_3 Vibronic Spectrum Composite - 180 fs PBF Comparison' font 'Helvetica,18'
set output month_day.'_RhF3_PBF_180fs.png'; set yr [ -1: 100]; set xr [ -1.0: -0.4]
plot \
f06 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H1 SOC 180fs PBF5', \
f05 using 1:2 w lines lw 1.00 lc 'web-blue'   dashtype 1 title 'MCTDH H2 SOC 180fs PBF5', \
f02 using 1:2 w lines lw 1.00 lc 'violet'     dashtype 1 title 'MCTDH H1 SOC 180fs PBF40', \
f01 using 1:2 w lines lw 1.00 lc 'orange-red' dashtype 1 title 'MCTDH H2 SOC 180fs PBF40', \

#---------------------------------------------------------------------------------------------------

set title 'RhF_3 Vibronic Spectrum Composite - 1000 fs PBF Comparison' font 'Helvetica,18'
set output month_day.'_RhF3_PBF_1000fs.png'; set yr [ -1: 100]; set xr [ -1.0: -0.4]
plot \
f10 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H1 SOC PBF40', \
f09 using 1:2 w lines lw 1.00 lc 'web-blue'   dashtype 1 title 'MCTDH H2 SOC PBF40', \
f14 using 1:2 w lines lw 1.00 lc 'violet'     dashtype 1 title 'MCTDH H1 SOC PBF5', \
f13 using 1:2 w lines lw 1.00 lc 'orange-red' dashtype 1 title 'MCTDH H2 SOC PBF5', \

#---------------------------------------------------------------------------------------------------

set title 'RhF_3 Vibronic Spectrum Composite - 180 fs PBF Comparison' font 'Helvetica,18'
set output month_day.'_RhF3_50_100_180_1000fs.png'; set yr [ -1: 100]; set xr [ -1.0: -0.4]
plot \
f10 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H1 SOC 1000fs', \
f21 using 1:3 w lines lw 1.00 lc 'web-blue'   dashtype 1 title 'VECC    H1 SOC 180fs Z1', \
f34 using 1:3 w lines lw 1.00 lc 'violet'     dashtype 1 title 'VECC    H1 SOC 100fs Z1', \
f35 using 1:3 w lines lw 1.00 lc 'orange-red' dashtype 1 title 'VECC    H1 SOC 50fs Z1', \

# 'black'    
# 'yellow'   
# 'web-green'
# 'pink'     
# 'web-blue' 
# 'red'      
# 'dark-green
# 'gold'     
# 'chartreuse
# 'orange-red
# 'light-cyan
# 'purple'   
# 'cyan'     
# 'brown'    
# 'steelblue'
# 'coral'    
# 'seagreen' 
# 'goldrenrod
# 'navy'     
# 'violet'   
# 'turquoise'
# 'dark-plum'
# 'dark-pink'

# set title 'RhF_3 Vibronic Spectrum Composite - 180 fs ; 40 tau' font 'Helvetica,18'
# set output month_day.'_RhF3_waterfall_composite.png'; set yr [ -1: 100]; set xr [ -1.5: 0]
# plot \
# f02 using 1:2 w lines lw 2.25 lc 'black'      dashtype 1 title 'MCTDH H1 SOC', \
# f01 using 1:2 w lines lw 0.75 lc 'web-blue'   dashtype 1 title 'MCTDH H2 SOC', \
# f04 using 1:2 w lines lw 0.75 lc 'red'        dashtype 1 title 'MCTDH H1', \
# f03 using 1:2 w lines lw 0.75 lc 'dark-green' dashtype 1 title 'MCTDH H2', \
# f21 using 1:3 w lines lw 0.75 lc 'purple'     dashtype 1 title 'VECC    H1 SOC Z1', \
# f22 using 1:3 w lines lw 0.75 lc 'pink'       dashtype 1 title 'VECC    H1 SOC Z2', \
# f23 using 1:3 w lines lw 0.75 lc 'brown'      dashtype 1 title 'VECC    H1 SOC Z3', \
# f24 using 1:3 w lines lw 0.75 lc 'steelblue'  dashtype 1 title 'VECC    H2 SOC Z1', \
# f25 using 1:3 w lines lw 0.75 lc 'coral'      dashtype 1 title 'VECC    H2 SOC Z2', \
# f26 using 1:3 w lines lw 0.75 lc 'seagreen'   dashtype 1 title 'VECC    H2 SOC Z3', \
# f27 using 1:3 w lines lw 0.75 lc 'goldrenrod' dashtype 1 title 'VECC    H1     Z1', \
# f28 using 1:3 w lines lw 0.75 lc 'navy'       dashtype 1 title 'VECC    H1     Z2', \
# f29 using 1:3 w lines lw 0.75 lc 'violet'     dashtype 1 title 'VECC    H1     Z3', \
# f30 using 1:3 w lines lw 0.75 lc 'turquoise'  dashtype 1 title 'VECC    H2     Z1', \
# f31 using 1:3 w lines lw 0.75 lc 'dark-plum'  dashtype 1 title 'VECC    H2     Z2', \
# f32 using 1:3 w lines lw 0.75 lc 'dark-pink'  dashtype 1 title 'VECC    H2     Z3', \