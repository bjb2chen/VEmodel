set terminal pngcairo size 1200,800
#set style data line
set nologscale; set xzeroaxis; set samples 5000

set xr [ 2: 12]; set yr [ -20: 50]

set title 'CoF_3 Vibronic Spectrum' font 'Helvetica,18'
set xlabel 'Energy[eV]'

set key font 'Helvetica,13'
set key Left # left align the text
set key left top # move legend to upper left corner

# see https://stackoverflow.com/questions/17564497/gnuplot-break-y-axis-in-two-parts for plotting disconnected y axis
# or https://stackoverflow.com/questions/56005211/break-y-axis-in-three-parts-in-gnuplot

# see http://www.gnuplot.info/docs/loc11700.html for explanation of %b, %d
month_day = "".strftime("%b_%d", time(0))

###############################################################
# THE FOLLOWING ARE NOW WITH JUN 13 SOC CORRECTIONS!!!
###############################################################

f01 = 'op_CoF36Q_5st_PBF40_tf180.00_auto_total'    # PBF40/180fs/40tau quad SOC
f02 = 'op_CoF36Q_5st_PBF40_tf182.00_auto_total'    # PBF40/180fs/40tau lin SOC
f03 = 'op_CoF36Q_5st_PBF40_tf184.00_auto_total'    # PBF40/180fs/40tau quad no SOC
f04 = 'op_CoF36Q_5st_PBF40_tf186.00_auto_total'    # PBF40/180fs/40tau lin no SOC
f05 = 'op_CoF36Q_5st_PBF5_tf186.00_auto_total'     # vanilla PBF5/SPF3/180fs/40tau quad SOC 
f06 = 'op_CoF36Q_5st_PBF5_tf188.00_auto_total'     # vanilla PBF5/SPF3/180fs/40tau lin SOC 
f07 = 'op_CoF36Q_5st_PBF5_tf190.00_auto_total'     # vanilla PBF5/SPF3/180fs/40tau quad no SOC 
f08 = 'op_CoF36Q_5st_PBF5_tf192.00_auto_total'     # vanilla PBF5/SPF3/180fs/40tau lin no SOC
f09 = 'op_CoF36Q_5st_PBF40_tf1000.00_auto_total'   # PBF40/1000fs/40tau quad SOC
f10 = 'op_CoF36Q_5st_PBF40_tf1002.00_auto_total'   # PBF40/1000fs/40tau lin SOC
f11 = 'op_CoF36Q_5st_PBF40_tf1004.00_auto_total'   # PBF40/1000fs/40tau quad no SOC
f12 = 'op_CoF36Q_5st_PBF40_tf1006.00_auto_total'   # PBF40/1000fs/40tau lin no SOC
f13 = 'op_CoF36Q_5st_PBF5_tf1000.00_auto_total'    # vanilla PBF5/SPF3/1000fs/40tau quad SOC
f14 = 'op_CoF36Q_5st_PBF5_tf1002.00_auto_total'    # vanilla PBF5/SPF3/1000fs/40tau lin SOC
f15 = 'op_CoF36Q_5st_PBF5_tf1004.00_auto_total'    # vanilla PBF5/SPF3/1000fs/40tau quad no SOC
f16 = 'op_CoF36Q_5st_PBF5_tf1006.00_auto_total'    # vanilla PBF5/SPF3/1000fs/40tau lin no SOC
f17 = ''                                           # extra space for now
f18 = ''
f19 = ''
f20 = ''

# VECC runs
f21 = ''
f22 = ''
f23 = ''
f24 = ''
f25 = ''
f26 = ''
f27 = ''
f28 = ''
f29 = ''
f30 = ''
f31 = ''
f32 = ''


#---------------------------------------------------------------------------------------------------

set title 'CoF_3 Vibronic Spectrum Composite; 180fs, 40 tau' font 'Helvetica,18'
set output month_day.'_CoF3_waterfall_composite.png'; set yr [ -1: 110]; set xr [ -1.25: -0.2]
plot \
j8path14 using 1:2 w lines lw 0.75 lc 'pink'       dashtype 1 title 'MCTDH H2 180fs SOC', \
j8path15 using 1:2 w lines lw 0.75 lc 'web-blue'   dashtype 1 title 'MCTDH H1 180fs SOC', \
j8path16 using 1:2 w lines lw 0.75 lc 'red'        dashtype 1 title 'MCTDH H2 180fs no SOC', \
j8path17 using 1:2 w lines lw 0.75 lc 'dark-green' dashtype 1 title 'MCTDH H1 180fs no SOC', \
vanilla1 using 1:2 w lines lw 0.75 lc 'gold'       dashtype 1 title 'MCTDH H2 180fs SOC PBF5/SPF3', \
vanilla2 using 1:2 w lines lw 0.75 lc 'chartreuse' dashtype 1 title 'MCTDH H1 180fs SOC PBF5/SPF3', \
vanilla3 using 1:2 w lines lw 0.75 lc 'orange-red' dashtype 1 title 'MCTDH H2 180fs no SOC PBF5/SPF3', \
vanilla4 using 1:2 w lines lw 0.75 lc 'light-cyan' dashtype 1 title 'MCTDH H1 180fs no SOC PBF5/SPF3', \
j8path18 using 1:3 w lines lw 0.75 lc 'purple'     dashtype 1 title 'Z1 H1 180fs SOC', \
j8path19 using 1:3 w lines lw 0.75 lc 'cyan'       dashtype 1 title 'Z2 H1 180fs SOC', \
j8path20 using 1:3 w lines lw 0.75 lc 'brown'      dashtype 1 title 'Z3 H1 180fs SOC', \
j8path21 using 1:3 w lines lw 0.75 lc 'steelblue'  dashtype 1 title 'Z1 H2 180fs SOC', \
j8path22 using 1:3 w lines lw 0.75 lc 'coral'      dashtype 1 title 'Z2 H2 180fs SOC', \
j8path23 using 1:3 w lines lw 0.75 lc 'seagreen'   dashtype 1 title 'Z3 H2 180fs SOC', \
j8path26 using 1:3 w lines lw 0.75 lc 'goldrenrod' dashtype 1 title 'Z1 H1 180fs no SOC', \
j8path27 using 1:3 w lines lw 0.75 lc 'navy'       dashtype 1 title 'Z2 H1 180fs no SOC', \
j8path28 using 1:3 w lines lw 0.75 lc 'violet'     dashtype 1 title 'Z3 H1 180fs no SOC', \
j8path29 using 1:3 w lines lw 0.75 lc 'turquoise'  dashtype 1 title 'Z1 H2 180fs no SOC', \
j8path30 using 1:3 w lines lw 0.75 lc 'dark-plum'  dashtype 1 title 'Z2 H2 180fs no SOC', \
j8path31 using 1:3 w lines lw 0.75 lc 'dark-pink'  dashtype 1 title 'Z3 H2 180fs no SOC', \


#---------------------------------------------------------------------------------------------------

set title 'CoF_3 Vibronic Spectrum Composite' font 'Helvetica,18'
set output month_day.'_CoF3_longprop.png'; set yr [ -1: 110]; set xr [ -1.25: -0.2]
plot \
j8path25 using 1:2 w lines lw 2.00 lc 'black'      dashtype 1 title 'MCTDH H1 1000fs SOC', \
j8path24 using 1:2 w lines lw 0.75 lc 'web-blue'   dashtype 1 title 'MCTDH H2 1000fs SOC 5BF', \
j8path25 using 1:2 w lines lw 0.75 lc 'pink'       dashtype 1 title 'MCTDH H1 1000fs SOC 5BF', \

#---------------------------------------------------------------------------------------------------


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