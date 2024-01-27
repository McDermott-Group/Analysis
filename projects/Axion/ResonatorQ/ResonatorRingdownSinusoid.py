import numpy as np
from scipy.optimize import curve_fit
from dataChest import *
import matplotlib.pyplot as plt
from pathlib import Path

def damped_sinusoid(t, A, l, omega, phi):
    return A * np.exp(-l * t) * (np.cos(omega * t + phi))

d = dataChest(r"Axion/Axion4E/DaveAndAbigail/MKID 3 Slot D1/20240126/radiator_ringdown_delay_dense/HDF5Data")
delays = []
Qs = []
for dataSet in d.ls()[0]:
    try:
        d.openDataset(''+Path(dataSet).stem)
        print(''+Path(dataSet).stem)
        data = d.getData(variablesList=["Time","Amplitude t"])
        t = data[:,0]
        mag = data[:,1] - np.mean(data[:,1])
    except:
        continue

    try:
        cutoff = len(t)
        popt, pcov = curve_fit(damped_sinusoid, t[:cutoff]/1e3, mag[:cutoff], p0=[0.001, -1/40.0, 2 * np.pi * 1.25, 0],maxfev=18000)
        Q = 2*np.pi*5.5e3*(1 / popt[1] / 2)
        if Q > 5e5:
            Qs.append(Q)
            delays.append(d.getParameter('J1 Bias to RO')*20.0)
    except:
        pass
    plt.figure()
    plt.plot(t/1e3, mag)
    plt.plot(t/1e3, damped_sinusoid(t/1e3, *popt), label="Ringdown time = {:.2f} us".format(1 / popt[1] / 2))
    plt.xlim(0,5)
    plt.ylabel("Magnitude")
    plt.xlabel("Time (us)")
    plt.legend()
    plt.show()

Qis=np.array([1/(1/Q + 1/1.985e6) for Q in Qs])
plt.figure()

plt.plot(np.convolve(delays, np.ones(3)/3, mode='valid'),np.convolve(Qis, np.ones(3)/3, mode='valid'))
plt.xlim(1000,3000)
plt.ylim(0.675e6,0.725e6)
plt.ylabel('$Q_i$')
plt.xlabel('Delay ($\mu$s)')
plt.show()