from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from axionlib import Parity,T1
import matplotlib
matplotlib.use('TkAgg')

qb_id = 4
rad_id = 1

base_path = ['Axion', '2023-01-23 - DR2']
user = 'DCH'
device_name = 'Axion4A'
dates = ['01-25-23']
experiment_base_name ='T1'

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]
paths = [os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path)) for expt_path in expt_paths]


parity_rates=[]
fidelities=[]
biases=[]
for i in range(len(paths)):
    path=paths[i]
    expt_path=expt_paths[i]
    for dataSet in os.listdir(path):
        print(dataSet)
        t1 = T1(expt_path, dataSet)
        t1.plot(fit=True)


