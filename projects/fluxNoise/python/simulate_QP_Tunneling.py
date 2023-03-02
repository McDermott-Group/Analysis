import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import periodogram
import noiselib
import importlib
importlib.reload(noiselib)
import matplotlib.pyplot as plt

N = 5*8192 # length of line to calc self cpsd
n = 10   # number of files to average
charge_dispersion = 3e6 #pk-pk Hz
offset = 0e6
idle_time = 166e-9
cycle_time = 100e-6
T_qp = 10e-3 # QP tunneling rate

# general methods:
#   calc rand ng
#   choose odd or even parity
#   calc detuning due to ng
#   calc rotation in stationary frame
#   project

fs = 1./cycle_time
ng = 0.5 * np.ones((n,N)) # np.random.rand(n, N)
dP = np.random.rand(n, N) < cycle_time/T_qp
parity = 2 * np.mod(np.cumsum(dP,-1),2) - 1
detuning = offset + charge_dispersion/2. * parity * np.abs(np.sin(np.pi*ng)) # this is just using a sin instead of the actual fn for parity bands
rotation = 2. * np.pi * detuning * idle_time
M = np.random.rand(n, N) < (0.5 + 0.3*np.sin(rotation))
print('n=', n)
print('N=', N)
print(len(M))

# avg_cpsd, f = noiselib.partition_and_avg_psd(M, fs)
f, avg_cpsd = periodogram(M, fs=fs, return_onesided=True, axis=1)
avg_cpsd = np.mean(avg_cpsd, axis=0)
window_avg_psd = noiselib.window_averaging(avg_cpsd)

# Fit
def y(x, a, gamma):
    # return a/(4*np.pi**2*x**2 + gamma**2)
    return gamma*a/(np.pi**2*x**2 + gamma**2)
    
window_avg_psd = window_avg_psd[~np.isnan(f)]
f = f[~np.isnan(f)]
f = f[~np.isnan(window_avg_psd)]
window_avg_psd = window_avg_psd[~np.isnan(window_avg_psd)]
popt, pcov = curve_fit(y, f, window_avg_psd, bounds=(0,np.inf))
print(popt[-1], 1./popt[-1], np.sum(dP)/(N*n*cycle_time))

# Plot
fig, ax = plt.subplots()
ax.plot(f, np.abs(window_avg_psd))#, 'DisplayName', [samples(1:end-1),qubit(1:end-1)])
ax.plot(f, y(f,*popt), '--')
ax.set_title('1/f Averaged PSD')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('S_\eta (\eta^2/Hz)')
plt.draw()
plt.pause(0.05)
# plt.show()
