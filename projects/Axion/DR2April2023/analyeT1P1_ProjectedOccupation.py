from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from axionlib import Parity,T1
import matplotlib
matplotlib.use('TkAgg')


base_path = ['Axion', '2023-04-05 - DR2']
user = 'DCH'
device_name = 'Axion4A Blackbody'
dates = ['04-19-23','04-20-23','04-21-23']
experiment_base_name ='T1_BB_Projected_Occupation'

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]
paths = [os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path)) for expt_path in expt_paths]

for qb_id in [1,2]:
    t1s=[]
    p1s=[]
    temps=[]
    gammas=[]
    for i in range(len(paths)):
        path=paths[i]
        expt_path=expt_paths[i]
        for dataSet in os.listdir(path):
            try:
                d = dataChest(expt_path)
                d.openDataset(dataSet,modify=True)
                if d.getParameter('Qubit ID') != qb_id:
                    continue
                print(dataSet)
                t1 = T1(expt_path, dataSet,dependent_variable='Projected Occupation')
                fit_params = t1.fit(save=False)
                temps.append(d.getParameter('Blackbody Temperature')[0])
                p1s.append(fit_params[2])
                gammas.append((-1e6*fit_params[1])*fit_params[2]/(1-fit_params[2]))
                print(-1/fit_params[1])#This is T1
            except Exception as e:
                print(e)
    try:
        plot_gammas = []
        plot_temps = []
        for i in range(len(gammas)):
            if 1/gammas[i] > 1e-2:
                continue
            plot_temps.append(temps[i])
            plot_gammas.append(gammas[i])
        plt.semilogy([1000/temp for temp in plot_temps],[1/g for g in plot_gammas],linestyle='',marker='*',markersize=12)
    except:
        pass

plt.xlabel('1/T (1/K)')
plt.ylabel('1 Hz/$\\Gamma_\\uparrow$')
# plt.ylim(1e-1,4e-1)
plt.xlim(2,15)
plt.show()


