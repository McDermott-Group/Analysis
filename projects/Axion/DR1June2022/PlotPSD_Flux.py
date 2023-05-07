from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

qb_id = 3
rad_id = 2

base_path = ['Axion', '2022-06-27 - DR1']
user = 'DCH'
device_name = 'Axion2B'
date = '07-23-22'
experiment_base_name = 'Test_PSD_Q{:d}_J{:d}B_{:d}uA'.format(qb_id,rad_id,-4000)

expt_path = base_path + [user] + [device_name] + [date] + [experiment_base_name.replace(" ", "_")]

path = os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path))

time_per_iteration = 100000
# n500 -  1235 (confusing); n1000 - 1327; n1300 - 1262; n1400 - 1287; n1500 - 1412; n1600 - 1311; n1700 - 1318; n2000 - 1335; n2500 - 1261; n3000 - 1175; n3500 - 1084
# 100 - 1030; 1000 - 1014; 1800 - 1139; 2000 - 1201; 2100 - 1193; 2200 - 1221; 3000 - 1322; 4000 - 1234 ...I think this might just be a trivial heating thing
#other rad  ; 1000 - 1054; 1800 - 1062; 2000 - 1060; 3000 - 1082; 4000 - 1060
# This is super confusing; it does NOT seem like we are seeing only heating, but also we don't really see things come back down.
# Everything looks great in the negative direction and I think I just too too big of steps in the positive direction. Going to bias at -1500 and measure
def getPSD(state,N_max,bias):
    psd = [None]*N_max
    for i in range(N_max):
        freqs, psd[i] = signal.periodogram(state[i],fs=1e9/time_per_iteration,return_onesided=True)
    ave_psd = np.average(psd,axis=0)
    fs = 1e9/time_per_iteration
    def fit_PSD_target_function(f, T_parity, F_map):
        return 2*(4 * F_map ** 2 / T_parity) / ((2 / T_parity) ** 2 + (2 * np.pi * f) ** 2) + 2*(1 - F_map ** 2) / fs


    initial_guess_fidelity = np.arange(0.1,1,0.05)
    initial_guess_parity = np.logspace(-4,-1,5)#[1/10,1/30,1/100,1/300,1/600,1/1000,1/1500,1/2000,1/2500,1/3000]
    # initial_guess = [0.7]
    covariance = float('inf')
    lowest_r_squared  =  float('inf')
    for igf in initial_guess_fidelity:
        for igp in initial_guess_parity:
            params_curr, params_covariance_curr = curve_fit(
                fit_PSD_target_function, freqs[1:], ave_psd[1:], p0=[igp, igf], bounds=([0.0001,0],[0.1,1]),method='trf')

            residuals = ave_psd[1:]-fit_PSD_target_function(freqs[1:],*params_curr)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((ave_psd[1:]-np.mean(ave_psd[1:]))**2)
            r_squared = 1 - (ss_res / ss_tot)
            if r_squared < lowest_r_squared:
                lowest_r_squared = r_squared
                params = params_curr
    #
    T_parity = params[0]
    F_map = params[1]
    fit_psd = fit_PSD_target_function(freqs, T_parity, F_map)
    plt.figure(figsize=(5, 4))
    plt.loglog(freqs[1:], ave_psd[1:])
    plt.loglog(freqs[1:], fit_psd[1:])
    plt.title('Bias = {0:2.0f} mV | Parity Rate = {1:.0f} Hz | Fidelity = {2:2.0f}%'.format(1000*bias,1/T_parity,100*F_map))
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('PSD (1/Hz)')
    plt.ylim([4e-5, 4e-3])
    plt.tight_layout()
    plt.show()
    return 1/T_parity,F_map

    # return T_parity

biases = []
parity_rate =[]
i=0
# dataChest.cd(path)
for dataSet in os.listdir(path):
    try:
        d = dataChest(os.path.join(path),dataSet)
        d.cd(expt_path)
        d.openDataset(dataSet)
        print(d.getDatasetName())
    except:
        continue
    try:
        bias = d.getParameter('J3 Voltage Bias')[0]
        data = d.getData()
        data = data.transpose()
        N_max = int(np.amax(data[0]) + 1) # Number of Iterations
        M_max = int(np.amax(data[1]) + 1) # Number of Reps
        a = []
        for i in range(N_max):
            b = []
            result = list(np.where(data[0].astype(int) == i))
            for j in result:
                b.append((data[2][j]))
            a.append(b[0])
    except:
        continue
    try:
        rate, fidelity = getPSD(a,N_max,bias=bias)
        if rate < 10000 or fidelity < 0.9:
            biases.append(bias)
            parity_rate.append(rate)
    except:
        print('Fitting Failed')
        continue

print(biases)
print(parity_rate)

plt.title('Parity Rate vs Radiator Bias')
# plt.plot(4*times, np.abs(I + Q * 1j))
plt.semilogy([b for b in biases], parity_rate,marker='.',linestyle='None')#(b-0.062)*(490)*(0.02)*484
plt.xlabel('Radiator Frequency (GHz)')
plt.ylabel('Parity Rate (Hz)')
plt.pause(0.1)


