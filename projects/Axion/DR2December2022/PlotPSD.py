from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use('TkAgg')

qb_id = 4
rad_id = 2

base_path = ['Axion', '2023-01-23 - DR2']
user = 'DCH'
device_name = 'Axion4A'
dates = ['01-27-23']
experiment_base_name ='PSD_Q{:d}_J{:d}B'.format(qb_id,rad_id)

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]

paths = [os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path)) for expt_path in expt_paths]

time_per_iteration = 100000

#
def getPSD(state,N_max,bias):
    print(np.average(state))
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
    if False:
        print('Here')
        print(bias)
        plt.figure(figsize=(5, 4))
        plt.loglog(freqs[1:], ave_psd[1:])
        plt.loglog(freqs[1:], fit_psd[1:])
        plt.title('Bias = {0:2.0f} mV | Parity Rate = {1:.0f} Hz | Fidelity = {2:2.0f}%'.format(1000*bias,1/T_parity,100*F_map))
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('PSD (1/Hz)')
        #plt.ylim([4e-5, 4e-3])
        ply.ylim([1.6e-4,2.8e-4])
        plt.tight_layout()
        plt.show()
    return 1/T_parity,F_map

    # return T_parity

biases = []
parity_rate =[]
i=0
# dataChest.cd(path)
for i in range(len(paths)):
    path = paths[i]
    expt_path = expt_paths[i]
    for dataSet in os.listdir(path):
        print(dataSet)
        try:
            d = dataChest(os.path.join(path),dataSet)
            d.cd(expt_path)
            d.openDataset(dataSet)
        except Exception as e:
            print(e)
            print('Error retreiving data')
            continue
        try:
            bias = d.getParameter('J3 Voltage Bias')[0]
            # if np.abs(bias/0.49)>=0.178:#bias/0.49 < -0.20:
            #     continue
            # if bias <=-0.06699999999999992 or bias in np.arange(-0.30,-.345,0.01):#in np.arange(-.23,.275,0.01):
            #     continue
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
            print('Error retreiving data')
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

plt.figure()
plt.title('Parity Rate vs Radiator Bias')
# plt.plot(4*times, np.abs(I + Q * 1j))
plt.semilogy([b for b in biases], parity_rate,marker='.',linestyle='None')#(b-0.062)*(490)*(0.02)*484
plt.xlabel('Radiator Bias (V from DAC)')
plt.ylabel('Parity Rate (Hz)')
plt.pause(0.1)


