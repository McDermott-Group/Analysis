import noiselib
import importlib
importlib.reload(noiselib)
import numpy as np
import matplotlib.pyplot as plt

files = noiselib.matpaths('Q2', '05-07-20', 'Calibration_split', list(range(22,279+1)))
o = None

for i,f in enumerate(files):
    data = noiselib.loadmat(f)
    line = np.array(data['Single_Shot_Occupation'])
    if i == 0:
        o = np.zeros((len(files),len(line)))
    o[i] = line

fig, ax = plt.subplots()
im = ax.imshow(o)
ax.set_aspect(1./5.)
fig.colorbar()
ax.set_xlabel('Charge Offset Index (-0.25V - 0.25V)')
ax.set_ylabel('File')
plt.draw()
plt.pause(0.05)

dwells = [[],[],[],[],[],[],[],[],[],[],[]]
for i,f in enumerate(files[135+22:150+22]):
    data = noiselib.loadmat(f)
    sso = np.array(data['Single_Shot_Occupations'])
    for j in range(len(sso)):
        d0,d1 = noiselib.bin_dwell_times(sso[j])
        dwells[j] += d0

# fig, ax = plt.subplots()
# ax.hist(dwells, 30)
# ax.set_xlabel('Dwell Time in 0 State Before Excitation [Reps=500us]')
# ax.set_ylabel('Counts')
# ax.set_yscale('log')
# plt.draw()
# plt.pause(0.05)

fig, ax = plt.subplots()
h, bins = np.histogram(dwells[0], 30)
x = (bins[:-1] + bins[1:])/2
c = 'brrbbbrrrbb'
for i,d in enumerate(dwells):
    h,b = np.histogram(d, bins)
    ax.plot(x, h, c[i])
ax.set_xlabel('Dwell Time in 0 State Before Excitation [Reps=500us]')
ax.set_ylabel('Counts')
ax.set_yscale('log')
plt.draw()
plt.pause(0.05)