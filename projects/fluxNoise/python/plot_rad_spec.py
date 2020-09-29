import numpy as np
import matplotlib.pyplot as plt
import glob
import noiselib

path = 'Z:/mcdermott-group/data/fluxNoise2/sim_data/radiation spectroscopy/'

spec_files = glob.glob(path+'*.csv')

fig, ax = plt.subplots(1,1)

for f in spec_files:
    data = np.loadtxt(f, skiprows=25, delimiter=',')
    label = f.split('\\')[-1]
    ax.plot( data[:,1], data[:,2], label=label)

ax.set_yscale('log')
ax.set_ylabel('Counts')
ax.set_xlabel('Energy (keV)')
noiselib.legend()
plt.draw()
plt.pause(0.05)