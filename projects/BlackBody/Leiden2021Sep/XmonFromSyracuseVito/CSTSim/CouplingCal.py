from antennalib import AntennaCoupling
import matplotlib.pyplot as plt
import numpy as np
JJ2 = [4*1e3, None, 0, 1000*150]   #[R, L, C, A]
# JJ2 = [28*1e3, None, 0, 150*150]   #[R, L, C, A]
JQ1 = [14*1e3, 18.3*1e-9, 0, 300*150]
JQ2 = [14*1e3, 14.6*1e-9, 0, 300*150]

fileJ2 = "Z_syracuse_injector.txt"
fileQ1 = "Z_syracuse_xmon_4.9GHz.txt"
fileQ2 = "Z_syracuse_xmon_5.1GHz.txt"

J2 = AntennaCoupling()
J2.import_data(fileJ2, JJ2)
f = J2.Antenna["f"]
ecJ2 = J2.e_c_dB
# J2.plot()

Q1 = AntennaCoupling()
Q1.import_data(fileQ1, JQ1)
ecQ1 = Q1.e_c_dB
Z_ReQ1 = Q1.Antenna["Z_Re"]
#
Q2 = AntennaCoupling()
Q2.import_data(fileQ2, JQ2)
ecQ2 = Q2.e_c_dB
Z_ReQ2 = Q2.Antenna["Z_Re"]


# Q3 = AntennaCoupling()
# Q3.import_data(fileQ3, JQ3)
# Q3.plot()
# print(Q3.Junction["C"])
# print(Q3.Junction["L"])
# C = Q3.Junction["C"]
# L = Q3.Junction["L"]
# omega = 1/np.sqrt((C*L))
# print(omega/(6.28*1e9))
# Q3.plot_Q3Mode()

ecJ2 = ecJ2
k = 1

plt.plot(f, ecQ1+ecJ2*k, color='b', label='Q1')
plt.plot(f, ecQ1+ecJ2*0, color='b', label='Q1 bare')
# plt.plot(f, ecQ2+ecJ2*k, color='r', label='Q2')
# plt.plot(f, ecQ3+ecJ2*k, color='k', label='Q3')
# plt.plot(f, ecQ4+ecJ2*k, color='y', label='Q4')
plt.plot(f, ecJ2, color='y', label='J7')
plt.xlabel('Freq')
plt.ylabel('Coupling efficiency')
plt.legend()
plt.show()

# plt.plot(f, Z_ReQ1, color='b', label='Q1')
# plt.plot(f, Z_ReQ2, color='r', label='Q2')
# plt.plot(f, Z_ReQ3, color='k', label='Q3')
# plt.plot(f, Z_ReQ4, color='y', label='Q4')
# # plt.plot(f, ecJ2, color='y', label='J2')
# plt.xlabel('Freq')
# plt.ylabel('Coupling efficiency')
# plt.legend()
# plt.show()
