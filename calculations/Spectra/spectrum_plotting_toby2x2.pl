#~/LOCAL/mctdh/mctdh86.4/bin/binary/x86_64/autospec86 -o RhF3_Z3_H0_vibronic_constant_tf180.pl -f ACF_ABS_CC_RhF3_Z3_H0_vibronic_constant_tf180_normalized.txt -p 13000 -2 2 ev 40 1 
set terminal png size 1200,800
set output '/home/bjb2chen/gamess/vibronics/RhF3/SOC_15st/jun30_hermitian_testing/tz_testing/E2V1_realonly.png'
set nologscale
set xzeroaxis
set xlabel 'Energy[eV]'
set xr [ 0: 2]
set yr [ -3: 40]
set title 'E2V1'

set style line 1 lc rgb 'red' lw 1.50

set arrow from 0.0499999999997624, graph 0 to 0.0499999999997624, graph 1 nohead  ls 1
set arrow from 0.1499999999999924, graph 0 to 0.1499999999999924, graph 1 nohead  ls 1
set arrow from 0.2499999999999564, graph 0 to 0.2499999999999564, graph 1 nohead  ls 1
set arrow from 0.3499999999999905, graph 0 to 0.3499999999999905, graph 1 nohead  ls 1
set arrow from 0.4499999999988656, graph 0 to 0.4499999999988656, graph 1 nohead  ls 1
set arrow from 0.5499999999992733, graph 0 to 0.5499999999992733, graph 1 nohead  ls 1
set arrow from 0.6499999999990503, graph 0 to 0.6499999999990503, graph 1 nohead  ls 1
set arrow from 0.7500000000000202, graph 0 to 0.7500000000000202, graph 1 nohead  ls 1

set arrow from 0.8499999999999811, graph 0 to 0.8499999999999811, graph 1 nohead  ls 1
set arrow from 0.9499999999999649, graph 0 to 0.9499999999999649, graph 1 nohead  ls 1
set arrow from 1.0499999999999632, graph 0 to 1.0499999999999632, graph 1 nohead  ls 1
set arrow from 1.1499999999997785, graph 0 to 1.1499999999997785, graph 1 nohead  ls 1
set arrow from 1.2499999999993394, graph 0 to 1.2499999999993394, graph 1 nohead  ls 1
set arrow from 1.3499999999999863, graph 0 to 1.3499999999999863, graph 1 nohead  ls 1
set arrow from 1.4499999999999871, graph 0 to 1.4499999999999871, graph 1 nohead  ls 1
set arrow from 1.5499999999999829, graph 0 to 1.5499999999999829, graph 1 nohead  ls 1


set arrow from 1.6500000000000012, graph 0 to 1.6500000000000012, graph 1 nohead  ls 1
set arrow from 1.7500000000000011, graph 0 to 1.7500000000000011, graph 1 nohead  ls 1
set arrow from 1.8499999999998753, graph 0 to 1.8499999999998753, graph 1 nohead  ls 1
set arrow from 1.9499999999999416, graph 0 to 1.9499999999999416, graph 1 nohead  ls 1
set arrow from 2.0499999999999536, graph 0 to 2.0499999999999536, graph 1 nohead  ls 1
set arrow from 2.1499999999999484, graph 0 to 2.1499999999999484, graph 1 nohead  ls 1
set arrow from 2.2499999999999512, graph 0 to 2.2499999999999512, graph 1 nohead  ls 1
set arrow from 2.3499999999997847, graph 0 to 2.3499999999997847, graph 1 nohead  ls 1



set arrow from 2.4499999999999851, graph 0 to 2.4499999999999851, graph 1 nohead  ls 1
set arrow from 2.5500000000000074, graph 0 to 2.5500000000000074, graph 1 nohead  ls 1
set arrow from 2.6499999999999342, graph 0 to 2.6499999999999342, graph 1 nohead  ls 1
set arrow from 2.7500000000000040, graph 0 to 2.7500000000000040, graph 1 nohead  ls 1
set arrow from 2.8499999999999490, graph 0 to 2.8499999999999490, graph 1 nohead  ls 1
set arrow from 2.9499999999998923, graph 0 to 2.9499999999998923, graph 1 nohead  ls 1
set arrow from 3.0500000000000247, graph 0 to 3.0500000000000247, graph 1 nohead  ls 1
set arrow from 3.1499999999998995, graph 0 to 3.1499999999998995, graph 1 nohead  ls 1


set arrow from 3.2499999999999769, graph 0 to 3.2499999999999769, graph 1 nohead  ls 1
set arrow from 3.3500000000000130, graph 0 to 3.3500000000000130, graph 1 nohead  ls 1
set arrow from 3.4499999999979964, graph 0 to 3.4499999999979964, graph 1 nohead  ls 1
set arrow from 3.5499999999999909, graph 0 to 3.5499999999999909, graph 1 nohead  ls 1
set arrow from 3.6500000000000088, graph 0 to 3.6500000000000088, graph 1 nohead  ls 1
set arrow from 3.7499999999999702, graph 0 to 3.7499999999999702, graph 1 nohead  ls 1
set arrow from 3.8499999999999877, graph 0 to 3.8499999999999877, graph 1 nohead  ls 1
set arrow from 3.9499999999999873, graph 0 to 3.9499999999999873, graph 1 nohead  ls 1


f01 = 'E2V1_auto.pl'
f02 = 'E2V1_realonly_auto.pl'
f03 = 'E2V1_auto_10000fs.pl'
f04 = 'E2V1_auto_realonly_10000fs.pl'
f05 = 'E2V1_initstate1_auto.pl'



#plot '/home/bjb2chen/gamess/vibronics/RhF3/SOC_15st/jun30_hermitian_testing/tz_testing/1st/op_NH36Q_3st_PBF30_tf250.00_auto_total' using 1:2 with lines ls 1 lc 'black' title 'MCTDH (vibronic)'
plot \
f01 using 1:3 with lines lw 2 dashtype 1 lc 'black'    title 'E2V1', \
f03 using 1:3 with lines lw 2 dashtype 2 lc 'gray'     title 'E2V1 10000fs', \
f02 using 1:3 with lines lw 2 dashtype 1 lc 'web-blue' title 'E2V1 Real', \
f04 using 1:3 with lines lw 2 dashtype 2 lc 'pink'     title 'E2V1 Real 10000fs', \
f05 using 1:3 with lines lw 2 dashtype 2 lc 'green'    title 'E2V1 Init\_state=1 10000fs 0.01tout', \