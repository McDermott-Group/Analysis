from antennalib import AntennaCoupling, EffectiveArea
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

C_eff = 150 * 1e-21
# C_eff = 75 * 1e-21

JQ2 = [16.6 * 1e3, None, 0, 360 * 150, "Receiver"]  #
filePower_Absorbed = "EArea_Xmon_(syracuse)_power_absorbed.txt"
# filePower_Absorbed = "Xmon_(syracuse)_power_absorbed_4.05fF.txt"
fileQ2 = "EArea_Xmon_(syracuse)_impedance.txt"
# fileQ2 = "xmon_full-chip_Q2.txt"

Q2 = AntennaCoupling()
Q2.import_data(fileQ2, JQ2, C_eff=C_eff)
f_Q2 = Q2.Antenna["f"]
ecQ2 = Q2.e_c_dB
eQ2 = Q2.e_c
plt.plot(eQ2)
# plt.show()

EA = EffectiveArea()
EA.import_data(filePower_Absorbed)
f = EA.f
P = EA.P_absorbed*1e11
plt.plot(f, P)
plt.yscale('log')
plt.show()