from antennalib import AntennaCoupling
import matplotlib.pyplot as plt
import numpy as np
JJ7 = [4*1e3, None, 0, 1000*150]   #[R, L, C, A]
# JJ7 = [28*1e3, None, 0, 150*150]   #[R, L, C, A]
JQ1 = [14*1e3, None, 0, 300*150]
JQ2 = [14*1e3, None, 0, 300*150]

fileJ7 = "Z_syracuse_injector.txt"
fileQ1 = "Z_syracuse_xmon_5.1GHz.txt"
# fileQ1 = "Z_syracuse_xmon_5.1GHz_withRO_ChargeCoupler.txt"
# fileQ2 = "Z_syracuse_xmon_5.1GHz.txt"
# fileQ2_C = "Z_syracuse_xmon_5.1GHz_withRO_ChargeCoupler.txt"
# fileQ2_C = "xmon_coupler_chg-line.txt"

J7 = AntennaCoupling()
J7.import_data(fileJ7, JJ7)
f = J7.Antenna["f"]
ecJ7 = J7.e_c_dB
# J7.plot()

Q1 = AntennaCoupling()
Q1.import_data(fileQ1, JQ1)
f_Q1 = Q1.Antenna["f"]
ecQ1 = Q1.e_c_dB
Z_ReQ1 = Q1.Antenna["Z_Re"]
#
# Q2 = AntennaCoupling()
# Q2.import_data(fileQ2, JQ2)
# ecQ2 = Q2.e_c_dB
# Z_ReQ2 = Q2.Antenna["Z_Re"]
#
# Q2_C = AntennaCoupling()
# Q2_C.import_data(fileQ2_C, JQ2)
# ecQ2_C = Q2_C.e_c_dB
# Z_ReQ2_C = Q2_C.Antenna["Z_Re"]
# f_C = Q2_C.Antenna["f"]

Q1_PSD = np.array([
    [0.0, 100.56], [0.005, 95.14], [0.01, 93.09], [0.015, 105.55], [0.02, 290.07],
    [0.025, 464.17], [0.0275, 389.05], [0.03, 580.15], [0.0325, 1755.69],
    [0.0349, 2195.63], [0.0374, 2354.61], [0.04, 1303.42], [0.041, 1706.05],
    [0.042, 1083.28], [0.043, 1217.15], [0.044, 1030.59], [0.045, 1729.76],
    [0.046, 1727.71], [0.047, 1845.05], [0.048, 3773.84], [0.049, 2620.48],
    [0.05, 4296.21], [0.051, 4929.87], [0.052, 4990.22], [0.053, 6261.39],
    [0.054, 4160.39], [0.055, 6191.23], [0.056, 6714.55], [0.057, 5577.27],
    [0.058, 6579.94], [0.059,6267.85], [0.06, 4146.44], [0.061, 5740.52],
    [0.062, 4549.68], [0.063, 3588.3], [0.064, 1955.06], [0.065, 1598.07],
    [0.066, 1266.78], [0.067, 1964.89], [0.068, 2945.25], [0.069, 3255.49],
    [0.07, 1989.67], [0.07, 1989.67], [0.0725, 1769.13], [0.075, 943.87],
    [0.0775, 866.63], [0.08, 497.36], [0.0825, 411.43], [0.085, 397.27],
    [0.0875, 414.55], [0.09, 447.09], [0.0925, 384.19], [0.095, 443.21],
    [0.0975, 449.53], [0.1, 511.58], [0.1025, 568.87], [0.105, 656.25],
    [0.1075, 756.29], [0.11, 853.78], [0.1125, 907.46], [0.115, 1048.06],
    [0.1175, 1110.32], [0.12, 1217.71], [0.1225, 1255.77], [0.125, 1353.19],
    [0.1275, 1457.7], [0.13, 1673.28], [0.1325, 1746.96], [0.135, 1884.08],
    [0.1375, 1996.61], [0.14, 2121.16], [0.1425, 2342.79], [0.145, 3016.99],
    [0.1475, 3711.24], [0.15, 3845.02], [0.155, 4165.38], [0.16, 4976.58],
    [0.165, 5837.55], [0.17, 5307.72], [0.175, 6209.49], [0.18, 7348.67],
    [0.185, 8213.85], [0.19, 10537.9], [0.195, 9946.68], [0.2, 9379.72],
    [0.205, 9908.34], [0.21, 12343.38], [0.215, 12469.27], [0.22, 10444.79],
    [0.225, 17334.14], [0.23, 14653.09], [0.235, 27709.83], [0.24, 28095.84],
    [0.245, 39425.95]
])

ecJ7 = ecJ7
k = 1

plt.plot(f, ecQ1+ecJ7*k, color='b', label='Q1')
# plt.plot(f, ecQ1+ecJ7*0, color='b', label='Q1 bare')
# plt.plot(f, ecQ2+ecJ7*k, color='r', label='Q2')
# plt.plot(f_C, ecQ2_C+ecJ7*k, color='k', label='Q2_C')
# plt.plot(f, ecQ3+ecJ7*k, color='k', label='Q3')
# plt.plot(f, ecQ4+ecJ7*k, color='y', label='Q4')
plt.plot(f, ecJ7, color='y', label='J7')
plt.xlabel('Freq')
plt.ylabel('Coupling efficiency')
plt.legend()
plt.show()

# plt.plot(f, Z_ReQ1, color='b', label='Q1')
# plt.plot(f, Z_ReQ2, color='r', label='Q2')
# plt.plot(f, Z_ReQ3, color='k', label='Q3')
# plt.plot(f, Z_ReQ4, color='y', label='Q4')
# # plt.plot(f, ecJ7, color='y', label='J7')
# plt.xlabel('Freq')
# plt.ylabel('Coupling efficiency')
# plt.legend()
# plt.show()
