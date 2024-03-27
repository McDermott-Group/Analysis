import numpy as np
import matplotlib.pyplot as plt


f_SFQ = 2.264 # GHz
# [pulse_length(ns), T_Ramsey(ns),  T1()]
data_2d = np.asarray([
    [0, 989.2], [50, 986.6], [100, 985.9], [150, 983.8], [200, 962.5],
    [250, 964.7], [300, 961.7], [350, 962.8], [400, 961.2]
])
SFQ_pulse = data_2d[:, 0]
T_Ramsey = data_2d[:, 1]
Freq_shift = -1000/T_Ramsey
Phase_slips = SFQ_pulse * f_SFQ * 3
# print('SFQ_pulse=', SFQ_pulse)
# print('T_Ramsey=', T_Ramsey)
# print('Freq_shift=', Freq_shift)

# plt.plot(SFQ_pulse, Freq_shift, label='freq_shift')
# plt.xlabel('SFQ pulse length (ns)')
plt.plot(Phase_slips, Freq_shift-Freq_shift[0], label='freq_shift')
plt.xlabel('Phase_slips')
plt.ylabel('Freq Shift (MHz)')
plt.grid()
plt.legend()
plt.show()
