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

class OneStateAnalyze(object):

    def __init__(self):
        self.one_state_data = np.array([])
        self.time_data = np.asarray([])
        self.average_n2 = 1000  # this is the rolling avg data
        self.data_set = ''
        self.file_name = ''
        self.files = []

    def one_state_data_retrieve(self, data_path, file_name):
        self.data_set = file_name

        # file_info = FilePicker()
        # data_path = file_info.relativePath
        # file_name = file_info.selectedFileName
        self.file_name = file_name
        d = dataChest(data_path)
        d.openDataset(file_name)
        self.files = d.getParameter("files")
        variable_info = d.getVariables()
        independent_variable_info = variable_info[0]
        dependent_variable_info = variable_info[1]

        # print('dependent_variable_info=', dependent_variable_info)

        occupation_data_pre = d.getData(
            variablesList=[dependent_variable_info[0][0]])
        occupation_data_pre = np.asarray(occupation_data_pre).flatten()

        time_data = d.getData(variablesList=[independent_variable_info[0][0]])
        time_data = np.asarray(time_data).flatten()

        self.one_state_data = occupation_data_pre
        self.time_data = time_data

    def get_rolling_data_pre(self, rolling_avg):
        r_avg = rolling_avg
        data = self.one_state_data
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]




data_path = ['GapEngineer', 'Nb_GND_Dev06_Trap', 'Leiden_2020Jul', 'P1PSD',
             'LIU', 'DataAnalysis_test', 'P1_Parity_Interleave', 'HDF5Data']

file_Q6_NoPoison = 'dee1649ynm_2020Nov04_P1_Parity_Interleave_Interleave_P1_Neg100.hdf5'
Q6_NoPoison = OneStateAnalyze()
Q6_NoPoison.one_state_data_retrieve(data_path, file_Q6_NoPoison)

file_Q6_Neg20dBm = 'dee1659qjg_2020Nov04_P1_Parity_Interleave_Interleave_P1_Neg20.hdf5'
Q6_Neg20dBm = OneStateAnalyze()
Q6_Neg20dBm.one_state_data_retrieve(data_path, file_Q6_Neg20dBm)

file_Q6_Neg16dBm = 'dee1708hev_2020Nov05_P1_Parity_Interleave_Interleave_P1_Neg16.hdf5'
Q6_Neg16dBm = OneStateAnalyze()
Q6_Neg16dBm.one_state_data_retrieve(data_path, file_Q6_Neg16dBm)

file_Q6_Neg13dBm = 'dee1711qjz_2020Nov05_P1_Parity_Interleave_Interleave_P1_Neg13.hdf5'
Q6_Neg13dBm = OneStateAnalyze()
Q6_Neg13dBm.one_state_data_retrieve(data_path, file_Q6_Neg13dBm)


file_Q6_Neg10dBm = 'dee1727frd_2020Nov05_P1_Parity_Interleave_Interleave_P1_Neg10.hdf5'
Q6_Neg10dBm = OneStateAnalyze()
Q6_Neg10dBm.one_state_data_retrieve(data_path, file_Q6_Neg10dBm)

# python to pandas
r_avg = 10  # rolling average
t = Q6_NoPoison.time_data[:1800]

f, ax = plt.subplots(5, 2)

ax[0][0].plot(t, Q6_NoPoison.get_rolling_data_pre(r_avg)[0][:1800], label='No Poison')
ax[0][0].set_title('Q6 P1 Time Trace')

ax[1][0].plot(t, Q6_Neg20dBm.get_rolling_data_pre(r_avg)[0][:1800], label='Neg 20 dBm')
ax[2][0].plot(t, Q6_Neg16dBm.get_rolling_data_pre(r_avg)[0][:1800], label='Neg 16 dBm')
ax[3][0].plot(t, Q6_Neg13dBm.get_rolling_data_pre(r_avg)[0][:1800], label='Neg 13 dBm')
ax[4][0].plot(t, Q6_Neg10dBm.get_rolling_data_pre(r_avg)[0][:1800], label='Neg 10 dBm')

for i in range(5):
    ax[i][0].legend(loc=1)
    ax[i][0].grid()
    ax[i][0].set_ylim(0, 0.5)
    ax[i][0].set_ylabel('P1')
    ax[i][0].set_xlabel('file number (5 file = 3.7 min)')

bin_list = np.arange(0.1, 0.4, 0.002)

weights = np.ones_like(Q6_NoPoison.one_state_data) / len(Q6_NoPoison.one_state_data)
ax[0][1].hist(Q6_NoPoison.one_state_data, bins=bin_list, weights=weights, label='No Poison')
ax[0][1].set_title('P1 Histogram Scaled Probability')

weights = np.ones_like(Q6_Neg20dBm.one_state_data) / len(Q6_Neg20dBm.one_state_data)
ax[1][1].hist(Q6_Neg20dBm.one_state_data, bins=bin_list, weights=weights, label='Neg 20 dBm')

weights = np.ones_like(Q6_Neg16dBm.one_state_data) / len(Q6_Neg16dBm.one_state_data)
ax[2][1].hist(Q6_Neg16dBm.one_state_data, bins=bin_list, weights=weights, label='Neg 16 dBm')

weights = np.ones_like(Q6_Neg13dBm.one_state_data) / len(Q6_Neg13dBm.one_state_data)
ax[3][1].hist(Q6_Neg13dBm.one_state_data, bins=bin_list, weights=weights, label='Neg 13 dBm')

weights = np.ones_like(Q6_Neg10dBm.one_state_data) / len(Q6_Neg10dBm.one_state_data)
ax[4][1].hist(Q6_Neg10dBm.one_state_data, bins=bin_list, weights=weights, label='Neg 10 dBm')

for i in range(5):
    ax[i][1].legend(loc=1)
    ax[i][1].grid()
    ax[i][1].set_ylim(0, 0.1)
    ax[i][1].set_xlabel('P1')

# f.tight_layout()
plt.show()




