import numpy as np
from dataChest import *
import matplotlib.pyplot as plt
import scipy.optimize as sc

# d = dataChest(r"SW_Circmon\2023-03-15 - HQAN\SW\SW_Circmon\03-20-23\T1")
# contents = d.ls()
# files = contents[0]
# f = files[-2]
# print(f)
# d.openDataset(f)

# data = d.getData(variablesList=['Idle Gate Time', 'Amplitude'])
# # print(data)
# amplitude = data[:,1]
# idle = data[:,0]*1e-3
##################################
dc = dataChest(r"SW_Circmon\2023-03-15 - HQAN\SW\03-20-23\SWAP")
contents = dc.ls()
files = contents[0]
f = files[0]
print(f)
dc.openDataset(f)
varsList = dc.getVariables()
# print(varsList)
data = dc.getData(variablesList=['Frequency', 'Idle Gate Time', 'Amplitude'])

freqs = data[:,0]
idle_times = data[:,1]*1e-3
amplitudes = data[:,2]
T1s = []
usable_freqs = []
swap_freqs = np.arange(3.670, 5.835, 0.001)


def T1(x, a, b, c):
    return a*np.exp(-x/b) + c

for freq in swap_freqs:
    vert_cut = np.asarray(np.where(freqs==freq))
    first_index = vert_cut.min()
    last_index = vert_cut.max()

    idle = np.asarray(idle_times[first_index:last_index])
    amplitude = np.asarray(amplitudes[first_index:last_index])

    popt, pcov = sc.curve_fit(T1, idle, amplitude, [0.0007, 10, 0.0048], maxfev=50000)

    perr = np.mean(np.sqrt(np.diag(pcov)))

    # print(popt[1])
    # print(perr)
    if 0 < perr < 5:
        T1s.append(popt[1])
        usable_freqs.append(freq)
        # print(freq)
    else:
        pass
        

######################################
fig, ax = plt.subplots()
plt.plot(usable_freqs, T1s)
# plt.plot(idle, amplitude, ".")
# plt.plot(idle, popt[0]*np.exp(-idle/popt[1]) + popt[2], label="T1 = {:.2f} us".format(popt[1]))
plt.title("T1 at different frequencies")
plt.xlabel("Frequency (GHz)")
plt.ylabel("T1 (us)")
# plt.legend()
plt.show()