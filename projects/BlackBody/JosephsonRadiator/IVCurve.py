import numpy as np
import matplotlib.pyplot as plt

Rn = 5e3
R0 = 50

V_JJ = np.linspace(0, 1e-3, 100)
V_src = np.linspace(0, 1e-3, 100)

I_JJ = V_JJ/Rn
I_src_200mDAC = 0.2*202e-6-V_src/R0
I_src_100mDAC = 0.1*202e-6-V_src/R0
I_src_70mDAC = 0.07*202e-6-V_src/R0
I_src_50mDAC = 0.05*202e-6-V_src/R0
I_src_45mDAC = 0.045*202e-6-V_src/R0
I_src_40mDAC = 0.04*202e-6-V_src/R0
I_src_30mDAC = 0.03*202e-6-V_src/R0
I_src_10mDAC = 0.01*202e-6-V_src/R0

plt.plot(V_JJ*1e3, I_JJ*1e6, label='JJ')
# plt.plot(V_src*1e3, I_src_200mDAC*1e6, label='Input_200mDAC')
# plt.plot(V_src*1e3, I_src_100mDAC*1e6, label='Input_100mDAC')
plt.plot(V_src*1e3, I_src_70mDAC*1e6, label='Input_70mDAC')
# plt.plot(V_src*1e3, I_src_50mDAC*1e6, label='Input_50mDAC')
plt.plot(V_src*1e3, I_src_45mDAC*1e6, label='Input_45mDAC')
plt.plot(V_src*1e3, I_src_10mDAC*1e6, label='Input_10mDAC')
plt.xlabel('Voltage(mV)')
plt.ylabel('Current(uA)')
plt.grid()
plt.legend()
plt.show()