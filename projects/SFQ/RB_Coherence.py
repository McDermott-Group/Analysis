"""
Z:\mcdermott-group\data\sfq\MCM_NIST\LIU\MCM13\2022May30RB\RB_Coherence
IRB, use Identity gate as interleaved gate and sweep the identity gate time
"""

import numpy as np
import matplotlib.pyplot as plt

# data = np.array([
#     [10, 99.916, 0.134], [20, 99.525, 0.273], [30, 99.636, 0.208],
#     [40, 99.292, 0.251], [50, 99.375, 0.182], [60, 99.174, 0.191],
#     [70, 99.035, 0.216], [80, 98.623, 0.573]])
# data = np.array([
#     [10, 99.916, 0.134], [20, 99.525, 0.273], [30, 99.636, 0.208],
#     [40, 99.292, 0.251], [50, 99.375, 0.182], [60, 99.174, 0.191],
#     [70, 99.035, 0.216]])

data = np.array([
    [10, 99.867, 0.084], [20, 99.556, 0.099], [30, 99.552, 0.121],
    [40, 99.364, 0.178], [50, 99.427, 0.156], [60, 99.292, 0.140],
    [70, 99.118, 0.116], [80, 98.823, 0.185]])

T1_base = 26.0 # us
Twhite_base= 2 * T1_base # us
Texp_base = 3/(1/T1_base+1/Twhite_base)
print('Texp_base=', Texp_base)
k_base = 1/(Texp_base*1000) * 100 # from Zijun Chen's thesis

r_incoherent = 0.964    # 0.964% per clifford gate from purity
r_RB = 1.19    # 0.964% per clifford gate from RB/IRB

time = data[:, 0]
r = 100 - data[:, 1]
se = data[:, 2]
[[k, b], cov] = np.polyfit(time, r, 1, cov=True)
std = np.sqrt(np.diag(cov))

plt.errorbar(time, r, yerr=se, fmt="o", color='k', capsize=3.0, ms=5.0)
plt.plot(time, time * k + b, 'k--', label='Idle gate error= {:.3f}  $\pm$ {:.3f} % (per 10 ns)'
         .format(k * 10, std[0] * 10))

# plt.plot(time, time * r_RB/91, 'y--', label='IRB error= {:.3f} % (per 10 ns)'
#          .format(r_RB/91 * 10))

plt.plot(time, time * r_incoherent/91, 'r--', label='Purity incoherent error= {:.3f} % (per 10 ns)'
         .format(r_incoherent/91 * 10))

plt.plot(time, time * k_base, 'b--', label='Base incoherent error= {:.3f} % (per 10 ns)'
         .format(k_base * 10))
# plt.errorbar(79, 100-99.151, yerr=0.169, fmt="s", color='r', capsize=3.0, ms=8, label='X')
# plt.errorbar(79, 100-99.125, yerr=0.179, fmt="v", color='b', capsize=3.0, ms=8, label='Y')
# plt.errorbar(39, 100-99.503, yerr=0.160, fmt="*", color='r', capsize=3.0, ms=12, label='X/2')
# plt.errorbar(39, 100-99.383, yerr=0.128, fmt="D", color='b', capsize=3.0, ms=8, label='-X/2')
# plt.errorbar(39, 100-99.377, yerr=0.180, fmt="p", color='r', capsize=3.0, ms=10, label='Y/2')
# plt.errorbar(39, 100-99.338, yerr=0.169, fmt="h", color='b', capsize=3.0, ms=10, label='-Y/2')

plt.plot(79, 100-99.151, marker="s", color='r', ms=7, label='X')
plt.plot(39, 100-99.503, marker="*", color='r', ms=10, label='X/2')
plt.plot(39, 100-99.383, marker="D", color='r', ms=6, label='-X/2')
plt.plot(79, 100-99.125, marker="v", color='b', ms=7, label='Y')
plt.plot(39, 100-99.377, marker="p", color='b', ms=8, label='Y/2')
plt.plot(39, 100-99.338, marker="h", color='b', ms=8, label='-Y/2')

plt.title('Identity Gate Error vs Gate Time')

plt.legend()
plt.xlabel('Identity Gate Time (ns)')
plt.ylabel('Identity Gate Error (%)')
plt.show()
