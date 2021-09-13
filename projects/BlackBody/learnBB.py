from antennalib import BB_Radiation, getBBTensity
import numpy as np
import matplotlib.pyplot as plt

# BB = BB_Radiation()
# # BB.getDensity_freq()
# # BB.plotDensity_freq()
#
# BB.getDensity_temp()
# BB.plotDensity_temp()

f = np.linspace(1e6, 1e12, 1000)
T1 = 70e-3
T2 = 100e-3
T3 = 300e-3
T4 = 500e-3
T5 = 1000e-3

# print(I)
plt.loglog(f, getBBTensity(f, T1), label='70mK')
plt.loglog(f, getBBTensity(f, T2), label='100mK')
plt.loglog(f, getBBTensity(f, T3), label='300mK')
plt.loglog(f, getBBTensity(f, T4), label='500mK')
plt.loglog(f, getBBTensity(f, T5), label='1000mK')
plt.xlabel('Freq (Hz)')
plt.ylabel('Relative intensity')
plt.legend()
plt.grid()
plt.show()