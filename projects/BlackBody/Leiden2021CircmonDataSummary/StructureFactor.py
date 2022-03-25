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

UpOverParity = np.zeros(len(S_Minus))

for i in range(len(UpOverParity)):
    ratio = 1.0/(1.0+np.sqrt(8*27.5)*(S_Minus[i]/S_Plus[i]))
    UpOverParity[i] = ratio

# fig, axs = plt.subplots(3, sharex='col', figsize=(7.5, 8),
#                             gridspec_kw={'hspace': 0.15})
fig, axs = plt.subplots(2, sharex='col', figsize=(9, 9), gridspec_kw={'hspace': 0.1})
# ax1 = fig.add_subplot(111)
# ax2 = ax1.twiny()

axs[0].plot(S_x[1:], S_Plus[1:],  linewidth=4, color='orange', label='  ')
axs[0].plot(S_x[1:], S_Minus[1:],  linewidth=4, color='purple', label='  ')
axs[0].set_xlim([2, 16])
axs[0].set_ylim([0, 16])

# axs[0].set_ylabel('$S_{\pm}(\hbar\omega_{\\nu}/\Delta$)', fontsize=label_font)
# ax1.legend(loc=2, prop={'size': 24}, frameon=False)
axs[0].legend(loc=2, fontsize='medium', prop={'size': 24}, frameon=False)
axs[0].tick_params(labelsize=tick_font, direction='in')

axs[0].tick_params(axis="x", direction="in", which='both')
axs[0].tick_params(axis="y", direction="in", which='both')
axs[0].tick_params(axis="x", width=1, length=6, which='both')
axs[0].tick_params(axis="y", width=1, length=6, which='both')

axs[1].plot(S_x[1:], UpOverParity[1:],  linewidth=4, color='black')
axs[1].set_xlim([2, 14])
axs[1].set_ylim([-0.01, 1.01])
# axs[1].set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
axs[1].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
# axs[1].set_yscale('log')
# axs[1].set_ylabel('$\Gamma_{\uparrow}/\Gamma_{\mathrm{p}}$', fontsize=label_font)
axs[1].tick_params(labelsize=tick_font, direction='in')
# axs[1].set_xlabel('$\hbar\omega_{\\nu}/\Delta$', fontsize=label_font)
axs[1].tick_params(axis="x", direction="in", which='both')
axs[1].tick_params(axis="y", direction="in", which='both')
axs[1].tick_params(axis="x", width=1, length=6, which='both')
axs[1].tick_params(axis="y", width=1, length=6, which='both')

# new_tick_locations = np.array([2, 5, 10, 11, 13])

# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticks(new_tick_locations)
# ax2.set_xticklabels(f_range)
# ax2.set_xlabel('$\omega_{\\nu}/2\pi$ (GHz)', fontsize=label_font)
# ax2.tick_params(labelsize=tick_font, direction='in')

plt.tight_layout()
path = 'Z:\mcdermott-group\data\Antenna\PaperWriting\Figs\FiguresFromPythonandOthersForIllustrator'
plt.savefig(path + '\StructureFactor.pdf', format='pdf', bbox_inches='tight', dpi=1200)
plt.show()



