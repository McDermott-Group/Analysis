from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np

### Before find out the DAC offset
# f1 = 235
# f2 = 500
# f4 = 345
### Before find out the DAC offset
f1 = 26.3*4.8
f2 = 83.5*4.8
f3 = 16.5*4.8
f4 = 51*4.8

r1in = 90
r1out = 182
r2in = 50
r2out = 52.9
r3in = 108
r3out = 500
r4in = 70
r4out = 90

f = np.array([f1, f2, f3, f4])
one_over_rin = 1.0 / np.array([r1in, r2in, r3in, r4in])
one_over_rout = 1.0 / np.array([r1out, r2out, r3out, r4out])
one_over_r_avg = 1.0 / (
            np.array([r1in, r2in, r3in, r4in]) + np.array([r1out, r2out, r3out, r4out]))

f_norm = f / np.mean(f)
one_over_rin_norm = one_over_rin / np.mean(one_over_rin)
one_over_rout_norm = one_over_rout / np.mean(one_over_rout)
one_over_r_avg_norm = one_over_r_avg / np.mean(one_over_r_avg)

plt.plot(f_norm, label='Freq')
plt.plot(one_over_rin_norm, '--', label='1/r_in')
plt.plot(one_over_rout_norm, '--', label='1/r_out')
plt.plot(one_over_r_avg_norm, label='1/r_avg')
# plt.xlabel('J2 Freq (THz)')
# plt.ylabel('P1')
plt.grid()
plt.legend()
plt.show()

# plt.plot(f_norm/one_over_r_avg_norm, label='Freq*r_avg')
# # plt.xlabel('J2 Freq (THz)')
# # plt.ylabel('P1')
# plt.grid()
# plt.legend()
# plt.show()
