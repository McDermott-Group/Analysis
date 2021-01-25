"""
Under development by Chuan-Hong Liu
reference to Alex's adapative spectroscopy
TO DO:
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from noiselib import loadmat

from dataChest import dataChest
from dataChestUtils import FilePicker

class P1_Detuning(object):

    def __init__(self):
        self.one_state_data = np.array([])
        self.detuning_data = np.array([])
        self.time_data = np.asarray([])
        self.average_n2 = 1000  # this is the rolling avg data
        self.data_set = ''
        self.file_name = ''
        self.files = []
        self.label = ''

    def one_state_data_retrieve(self, data_path, file_name, label=None):
        self.label = label
        self.data_set = file_name

        self.file_name = file_name
        d = dataChest(data_path)
        d.openDataset(file_name)
        self.files = d.getParameter("files")
        variable_info = d.getVariables()
        independent_variable_info = variable_info[0]
        dependent_variable_info = variable_info[1]

        # print('dependent_variable_info=', dependent_variable_info)

        occupation_data = d.getData(
            variablesList=[dependent_variable_info[-1][0]])
        occupation_data = np.asarray(occupation_data).flatten()

        detuning_data = d.getData(
            variablesList=[dependent_variable_info[-2][0]])
        detuning_data = 10 * np.asarray(detuning_data).flatten()

        time_data = d.getData(variablesList=[independent_variable_info[0][0]])
        time_data = np.asarray(time_data).flatten()

        self.one_state_data = occupation_data
        self.detuning_data = detuning_data
        self.time_data = time_data
        # print(len(occupation_data))
        # print(len(detuning_data))

    def get_rolling_data_pre(self, rolling_avg):
        r_avg = rolling_avg
        data = self.one_state_data
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]

data_path = ['GapEngineer', 'Nb_GND_Dev06_Trap', 'Leiden_2020Jul', 'P1PSD',
             'LIU', 'DataAnalysis_test', 'P1_Parity_Interleave', 'HDF5Data']

file_Q4_NoPoison = 'dgg1555cla_2020Dec28_P1_Parity_Interleave_Reset_Interleave_P1_Neg100.hdf5'
# file_Q4_NoPoison = 'dgm1557hyi_2021Jan04_P1_Parity_Interleave_Reset_Interleave_P1_Neg100.hdf5'
# file_Q4_NoPoison = 'dgm1558gkv_2021Jan04_P1_Parity_Interleave_Reset_Interleave_P1_Neg16.hdf5'
# file_Q4_NoPoison = 'dgm1932vkn_2021Jan05_P1_Parity_Interleave_Reset_Interleave_P1_Neg16.hdf5'
# file_Q4_NoPoison = 'dgm1933hem_2021Jan05_P1_Parity_Interleave_Reset_Interleave_P1_Neg13.hdf5'
# file_Q4_NoPoison = 'dgm1933vdb_2021Jan05_P1_Parity_Interleave_Reset_Interleave_P1_Neg10.hdf5'
Q4_NoPoison = P1_Detuning()
Q4_NoPoison.one_state_data_retrieve(data_path, file_Q4_NoPoison, label='NoPoison')
# print('len')

P1 = Q4_NoPoison.one_state_data
df = Q4_NoPoison.detuning_data
fig = plt.figure()
plt.scatter(P1, df)
plt.xlabel('P1')
plt.ylabel('df (MHz)')
plt.title('P1 df Correlation (No Poison)')
plt.legend(loc=2)
plt.grid()
plt.show()

