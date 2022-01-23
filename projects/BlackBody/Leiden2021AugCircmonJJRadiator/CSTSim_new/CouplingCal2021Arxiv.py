from antennalib import AntennaCoupling
import matplotlib.pyplot as plt
import numpy as np

e_eff = 6  # limit (1, 6.5), the voltage can also be built in to have a larger range
C_eff = 75 * 1e-21  # Commonly used (50-100) 75, 88

JJ2 = [4.8*1e3, None, 0, 320*123*2, "Radiator"]   #[R, L, C, A]
JQ1 = [16.6*1e3, 18.3*1e-9, 0, 193.8*121.8*2, "Receiver"]   #
JQ2 = [13.2*1e3, 14.6*1e-9, 0, 350*126.6, "Receiver"]   #
JQ3 = [19*1e3, 21*1e-9, 0, 310*126, "Receiver"]   #    some issue with Q3
JQ4 = [15*1e3, 19.9*1e-9, 0, 184.4*122.5*2, "Receiver"]   #    some issue with Q3

fileJ2 = "J2.txt"
# fileQ1 = "Q1.txt"
fileQ1 = "Q1_leads.txt"
fileQ2 = "Q2.txt"
fileQ3 = "Q3.txt"
fileQ4 = "Q4.txt"

J2 = AntennaCoupling()
J2.import_data(fileJ2, JJ2, C_eff=C_eff)
f = J2.Antenna["f"]
ecJ2 = J2.Antenna["e_c_dB"]
eJ2 = J2.Antenna["e_c"]
pgJ2 = J2.Radiator["Gamma_rad"]
x_qpJ2 = J2.Radiator["X_QP"]
refJ2 = J2.Al_Wall["Ref"]
PhotonFlux = J2.Al_Wall["PhotonFlux"]

# Q1 = AntennaCoupling()
# Q1.import_data(fileQ1, JQ1)
# # Q1.plot()
# ecQ1 = Q1.e_c_dB
# Z_ReQ1 = Q1.Antenna["Z_Re"]
#
# Q2 = AntennaCoupling()
# Q2.import_data(fileQ2, JQ2)
# ecQ2 = Q2.e_c_dB
# Z_ReQ2 = Q2.Antenna["Z_Re"]
#
# Q3 = AntennaCoupling()
# Q3.import_data(fileQ3, JQ3)
# ecQ3 = Q3.e_c_dB
# Z_ReQ3 = Q3.Antenna["Z_Re"]
#
# Q4 = AntennaCoupling()
# Q4.import_data(fileQ4, JQ4)
# ecQ4 = Q4.e_c_dB
# Z_ReQ4 = Q4.Antenna["Z_Re"]
#
# # Q3 = AntennaCoupling()
# # Q3.import_data(fileQ3, JQ3)
# # Q3.plot()
# # print(Q3.Junction["C"])
# # print(Q3.Junction["L"])
# # C = Q3.Junction["C"]
# # L = Q3.Junction["L"]
# # omega = 1/np.sqrt((C*L))
# # print(omega/(6.28*1e9))
# # Q3.plot_Q3Mode()
#

f_scale = np.sqrt(e_eff / 6.0)
f = f / 1e9
f = f * f_scale

plt.plot(f, eJ2)
plt.show()

# f_Q1 = f_Q1 / 1e9
# f_Q1 = f_Q1 * f_scale
# f_Q2 = f_Q2 / 1e9
# f_Q2 = f_Q2 * f_scale

# ecJ2 = ecJ2*1.0
# k = 1
#
# plt.plot(f, ecQ1+ecJ2*k, color='b', label='Q1')
# plt.plot(f, ecQ2+ecJ2*k, color='r', label='Q2')
# plt.plot(f, ecQ3+ecJ2*k, color='k', label='Q3')
# plt.plot(f, ecQ4+ecJ2*k, color='y', label='Q4')
# # plt.plot(f, ecJ2, color='y', label='J2')
# plt.xlabel('Freq')
# plt.ylabel('Coupling efficiency')
# plt.legend()
# plt.show()

# plt.plot(f, Z_ReQ1, color='b', label='Q1')
# plt.plot(f, Z_ReQ2, color='r', label='Q2')
# plt.plot(f, Z_ReQ3, color='k', label='Q3')
# plt.plot(f, Z_ReQ4, color='y', label='Q4')
# # plt.plot(f, ecJ2, color='y', label='J2')
# plt.xlabel('Freq')
# plt.ylabel('Coupling efficiency')
# plt.legend()
# plt.show()
