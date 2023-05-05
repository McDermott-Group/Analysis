from axionlib import *
from dataChest import *
import os

base_path = ['Axion', '2023-04-05 - DR2']
user = 'DCH'
device_name = 'Axion4A JJ Radiator'
dates = ['04-07-23']
experiment_base_name ='P1_vs_JB1'

expt_paths =  [base_path + [user] + [device_name] + [d] + [experiment_base_name.replace(" ", "_")] for d in dates]
paths = [os.path.join(*([r'Z:\mcdermott-group\data'] + expt_path)) for expt_path in expt_paths]

qb_data = {i:[] for i in range(1,7)}
for path,expt_path in zip(paths,expt_paths):
    for dataSet in os.listdir(path):
        P1(expt_path,dataSet).plot_versus_radiator_frequency(save=True)