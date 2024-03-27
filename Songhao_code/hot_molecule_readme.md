CH2O_VECC_cProfile_data  CH2O_VECC_cProfile_data.txt  CH2O_VECC_thermal_data_TFCC_RK.csv  hot_molecule_test.py	run.sh
bjb2chen@nlogn:~/VECC/hot-molecule$ nano CH2O_VECC_cProfile_data.txt
bjb2chen@nlogn:~/VECC/hot-molecule$ nano CH2O_VECC_cProfile_data
bjb2chen@nlogn:~/VECC/hot-molecule$ nano CH2O_VECC_thermal_data_TFCC_RK.csv
bjb2chen@nlogn:~/VECC/hot-molecule$ python3
Python 3.8.10 (default, Nov 14 2022, 12:59:47) 
[GCC 9.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import pandas as pd
>>> df = pd.readcsv('CH2O_VECC_thermal_data_TFCC_RK.csv')
KeyboardInterrupt
>>> df = pd.read_csv("CH2O_VECC_thermal_data_TFCC_RK.csv")
>>> df(head)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'head' is not defined
>>> df.head()
    temperature  internal energy  partition function  free energy  chemical entropy  q variance 1  ...  state 1 pop   state 2 pop   state 3 pop   state 4 pop   state 5 pop  state 6 pop
0  10000.000000         5.210002         2905.063511    -6.871643          0.001208      5.813187  ...     0.000010  2.370183e-07  3.233681e-08  2.304837e-08  2.039895e-08     0.999989
1   9678.355970         5.044974         2383.917606    -6.485727          0.001191      5.627149  ...     0.000007  1.391580e-07  1.846959e-08  1.466359e-08  1.230753e-08     0.999992
2   9183.684963         4.791364         1735.429973    -5.902977          0.001164      5.341097  ...     0.000004  5.765410e-08  7.323448e-09  7.013389e-09  5.480833e-09     0.999996
3   8872.028525         4.631714         1407.987958    -5.542795          0.001147      5.160922  ...     0.000003  3.147571e-08  3.881087e-09  4.336230e-09  3.231253e-09     0.999997
4   8580.830484         4.482646         1150.177032    -5.211321          0.001130      4.992610  ...     0.000002  1.723270e-08  2.068213e-09  2.747979e-09  1.964257e-09     0.999998

[5 rows x 29 columns]
