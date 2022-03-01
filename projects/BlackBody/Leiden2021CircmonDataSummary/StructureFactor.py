import matplotlib.pyplot as plt
from antennalib import AntennaCoupling, UpAndParity, getPhotonRate
import numpy as np
from scipy import interpolate
from scipy.integrate import quad


S_param_file = "S_minusplus.txt"

label_font = 25
tick_font = 22
legend_font = 10

Delta = 46.0
f_range = np.array([100, 200, 300, 400, 500, 600, 700])
new_tick_locations = f_range/Delta

S_param = np.loadtxt(S_param_file)
S_x = S_param[:, 0]
S_f = S_param[:, 0]*Delta
S_Minus = S_param[:, 1]
S_Plus = S_param[:, 2]
# S_MinusOverPlus = np.divide(S_param[:, 1], S_param[:, 2])

fig = plt.figure(figsize=(9, 6))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

ax1.plot(S_x[1:], S_Plus[1:],  linewidth=4, color='orange', label='$S_{+}(\hbar\omega_{\\nu}/\Delta)$')
ax1.plot(S_x[1:], S_Minus[1:],  linewidth=4, color='purple', label='$S_{-}(\hbar\omega_{\\nu}/\Delta)$')
ax1.set_xlim([2, 16])
ax1.set_ylim([0, 16])
ax1.set_xlabel('$\hbar\omega_{\\nu}/\Delta$', fontsize=label_font)
ax1.set_ylabel('$S_{\pm}(\hbar\omega_{\\nu}/\Delta$)', fontsize=label_font)
# ax1.legend(loc=2, prop={'size': 24}, frameon=False)
ax1.legend(loc=2, fontsize='medium', prop={'size': 24}, frameon=False)
ax1.tick_params(labelsize=tick_font, direction='in')

# new_tick_locations = np.array([2, 5, 10, 11, 13])

ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(f_range)
ax2.set_xlabel('$\omega_{\\nu}/2\pi$ (GHz)', fontsize=label_font)
ax2.tick_params(labelsize=tick_font, direction='in')

plt.tight_layout()
path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
plt.savefig(path + '\StructureFactor.pdf', format='pdf', dpi=1200)
# plt.savefig(path + '\StructureFactor.pdf', format='pdf', bbox_inches='tight', dpi=1200)
plt.show()



