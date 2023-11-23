from antennalib import AntennaCoupling, getNoiseBandwidth, getTbb, getPhotonRate
from scipy import interpolate
from scipy.special import ellipk,ellipe,ellipkm1
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.font_manager as font_manager
from matplotlib.ticker import ScalarFormatter, NullFormatter, FormatStrFormatter
import os
import matplotlib
#matplotlib.use("tkAgg")

def theta(x):
  return np.heaviside(x,0)

def j1(x):
  return -1.0*ellipk(x**2)/ellipk(0) if x < 1 else (-1.0/x)*ellipk(1/x**2)/ellipk(0)

def j2(x):
  return (-1.0/x)*ellipk((x**2-1)/x**2)/ellipk(0) if x > 1 else 0

def qp(x):
  return (2*x*ellipe((x**2-1)/x**2)-(1/x)*ellipk((x**2-1)/x**2))/ellipk(0) if x > 1 else 0

fig = plt.figure(figsize=(20, 12))
axs = fig.gca()

label_font = 44
tick_font = 44
legend_font = 15

legend_font = font_manager.FontProperties(
  # weight='bold',
  style='normal', size=legend_font)

xs = np.arange(0,8,0.00001)

axs.plot([2*x for x in xs],[-1*j1(x) for x in xs], color="red", linestyle='-',  linewidth=8)
axs.plot([2*x for x in xs], [-1*j2(x) for x in xs], color="blue", linestyle='-',  linewidth=8)
axs.plot([2*x for x in xs], [qp(x) for x in xs], color="black", linestyle='-',  linewidth=8)

axs.set_xlim([0, 8])
axs.set_ylim([0, 4])
axs.set_yticks(np.arange(1, 4.1, 1))

# axs[1].legend(loc=2, prop=legend_font, frameon=False)
axs.set_ylabel("Current ($I_0$)", color="black", fontsize=label_font)
axs.set_xlabel("Voltage ($\Delta/e$)", color="black", fontsize=label_font)
axs.tick_params(labelsize=tick_font)
#
def f_to_mV(f):
    return f / 484

def mV_to_f(mV):
    return 484 * mV

def f_to_delta(f):
    return f / 46

def delta_to_f(delta):
    return 46 * delta

def fj_to_delta(f):
    return f / 92

def delta_to_fj(delta):
    return 92 * delta

double_axis = False

axs.axvspan(368,1100, facecolor='0.2', alpha=0.1)
if double_axis:
    secax = axs.secondary_xaxis('top', functions=(delta_to_fj, fj_to_delta))
    secax.tick_params(labelsize=tick_font)
    secax.xaxis.set_major_formatter(FormatStrFormatter('%2d'))
    secax.set_xticks(np.arange(0, 800, 92))
    secax.set_xlabel("Transmitter Josephson Frequency (GHz)", fontsize=label_font, labelpad=20)

    secax.tick_params(axis="x", direction="in", which='both')
    secax.tick_params(axis="y", direction="in", which='both')

    secax.tick_params(axis="x", width=1, length=6, which='both')
    secax.tick_params(axis="y", width=1, length=3, which='minor')
    secax.tick_params(axis="y", width=1, length=6, which='major')

axs.tick_params(axis="x", direction="in", which='both')
axs.tick_params(axis="y", direction="in", which='both')

axs.tick_params(axis="x", width=1, length=6, which='both')
axs.tick_params(axis="y", width=1, length=3, which='minor')
axs.tick_params(axis="y", width=1, length=6, which='major')


plt.savefig("S_JJ_Current_Contributions.png", format="png")
plt.show()


double_axis = True

fig = plt.figure(figsize=(20, 12))
axs = fig.gca()
axs.plot([2*x for x in xs],[(j1(x)**2+j2(x)**2)**0.5 for x in xs], color="purple", linestyle='-',  linewidth=8)

axs.axvspan(368,1100, facecolor='0.2', alpha=0.1)

if double_axis:
    secax = axs.secondary_xaxis('top', functions=(delta_to_fj, fj_to_delta))
    secax.tick_params(labelsize=tick_font)
    secax.xaxis.set_major_formatter(FormatStrFormatter('%2d'))
    secax.set_xticks(np.arange(100, 850, 100))
    secax.set_xlabel("Transmitter Josephson Frequency (GHz)", fontsize=label_font, labelpad=20)
    secax.tick_params(axis="x", direction="in", which='both')
    secax.tick_params(axis="y", direction="in", which='both')

    secax.tick_params(axis="x", width=1, length=6, which='both')
    secax.tick_params(axis="y", width=1, length=3, which='minor')
    secax.tick_params(axis="y", width=1, length=6, which='major')

axs.tick_params(axis="x", direction="in", which='both')
axs.tick_params(axis="y", direction="in", which='both')

axs.tick_params(axis="x", width=1, length=6, which='both')
axs.tick_params(axis="y", width=1, length=3, which='minor')
axs.tick_params(axis="y", width=1, length=6, which='major')


axs.set_xlim([0, 8])
axs.set_ylim([0, 4])
axs.set_yticks(np.arange(1, 4.1, 1))
axs.set_ylabel("Current ($I_0$)", color="black", fontsize=label_font)
axs.set_xlabel("Voltage ($\Delta/e$)", color="black", fontsize=label_font)
axs.tick_params(labelsize=tick_font)

plt.savefig("S_JosephsonCurrents.png", format="png")
plt.show()
