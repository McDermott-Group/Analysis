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

Q1 = AntennaCoupling()
Q1.import_data(fileQ1, JQ1, C_eff=C_eff)
f_Q1 = Q1.Antenna["f"]
eQ1 = Q1.Antenna["e_c"]

Q2 = AntennaCoupling()
Q2.import_data(fileQ2, JQ2, C_eff=C_eff)
f_Q2 = Q2.Antenna["f"]
ecQ2 = Q2.Antenna["e_c_dB"]
eQ2 = Q2.Antenna["e_c"]
Area = Q2.Receiver["Area"]

Q4 = AntennaCoupling()
Q4.import_data(fileQ4, JQ4, C_eff=C_eff)
f_Q4 = Q4.Antenna["f"]
ecQ4 = Q4.Antenna["e_c_dB"]
eQ4 = Q4.Antenna["e_c"]
# Area = Q4.Receiver["Area"]



f_scale = np.sqrt(e_eff / 6.0)
f = f / 1e9
f = f * f_scale
f_Q1 = f_Q1 / 1e9
f_Q1 = f_Q1 * f_scale
f_Q2 = f_Q2 / 1e9
f_Q2 = f_Q2 * f_scale
f_Q4 = f_Q4 / 1e9
f_Q4 = f_Q4 * f_scale

# plt.plot(f, eJ2)
# plt.show()

# ecJ2 = ecJ2*1.0
# k = 1

# fig, axs = plt.subplots(2, sharex='col', figsize=(12, 8))
l_i = 50
r_i = 500
plt.plot(f, eJ2, color="black", marker="o", markersize=4,
               label='Radiator')
plt.plot(f_Q1, eQ1, color="blue", marker="o", markersize=4, label='Q1')
plt.plot(f_Q2, eQ2, color="red", marker="o", markersize=4, label='Q2')
plt.plot(f_Q4, eQ4, color="green", marker="o", markersize=4, label='Q4')
plt.ylabel('Coupling Efficiency', color="black", fontsize=10)
plt.yscale('log')
plt.xlim([l_i, r_i])
plt.ylim([5e-3, 2e-1])
plt.legend(loc=4)
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
