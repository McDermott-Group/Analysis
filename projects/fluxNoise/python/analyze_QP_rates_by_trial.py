import numpy as np
import matplotlib.pyplot as plt
import noiselib
import importlib
importlib.reload(noiselib)
from noiselib import movingmean
import QPTunneling
importlib.reload(QPTunneling)
from QPTunneling import *
import ChargeOffset
importlib.reload(ChargeOffset)
from ChargeOffset import *
import TwoMeasDataFile
importlib.reload(TwoMeasDataFile)
from TwoMeasDataFile import *
from dataChest import dataChest
from random import randrange
from random import randrange


offset_path = 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\Parameter'
CO = ChargeOffset()
CO.add_dataset(offset_path + '\cvo2209zfd_correlation.hdf5')
path = 'Z:\\mcdermott-group\\data\\fluxNoise2\\DR1 - 2019-12-17\\CorrFar\\Q1Q2Corr\\General\\03-27-20\\Charge_resetting_QP\MATLABData\\Charge_resetting_QP_{:03d}.mat'
flip_rates = np.array([])
QPT = QPTunneling()
for i in range(0,588,2): #588
    DF = TwoMeasDataFile(path.format(i))
    # DF.apply_infidelity_correction(9)
    flips = np.abs(np.diff(DF.o_meas, axis=1))
    flips = movingmean(flips,500)
    QPT.add_data(flips)
    flip_rates = np.append(flip_rates, np.mean(flips, axis=1) )

QPT.plot_psd()

jumps,sigma = CO.get_jump_sizes()
large_jumps = np.argwhere(jumps['Q1'] > 3*sigma['Q1'])
flip_rates_at_jumps = flip_rates[ np.repeat(large_jumps,10) + list(range(10))*len(large_jumps) ]
indicies = np.arange(len(flip_rates))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('')
ax.set_xlabel('Time')
# ax.set_ylabel('Flip Rate')
ax.plot(np.abs(jumps['Q1']), label='Charge Jump Size')
ax.plot(0.1*indicies, flip_rates, '.-', markevery=10, label='Flip Rate')
ax.legend()
plt.draw()
plt.pause(0.05)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Change in flip rate')
ax.set_xlabel('Flip Rate')
ax.set_ylabel('Counts')
ax.hist(np.diff(flip_rates), 100, label='All')
ax.hist(np.diff(flip_rates_at_jumps), 100, label='Jumps > 3*sigma')
ax.legend()
plt.draw()
plt.pause(0.05)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Flip Rate')
ax.set_xlabel('Flip Rate')
ax.set_ylabel('Counts')
ax.hist(flip_rates, 100, label='All')
ax.hist(flip_rates_at_jumps, 100, label='Jumps > 3*sigma')
ax.legend()
plt.draw()
plt.pause(0.05)
