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
        self.label = ''

    def one_state_data_retrieve(self, data_path, file_name, label=None):
        self.label = label
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

        occupation_data = d.getData(
            variablesList=[dependent_variable_info[-1][0]])
        occupation_data = np.asarray(occupation_data).flatten()

        time_data = d.getData(variablesList=[independent_variable_info[0][0]])
        time_data = np.asarray(time_data).flatten()

        self.one_state_data = occupation_data
        self.time_data = time_data

    def get_rolling_data_pre(self, rolling_avg):
        r_avg = rolling_avg
        data = self.one_state_data
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]

data_path = ['GapEngineer', 'Nb_GND_Dev06_Trap', 'Leiden_2021Jan', 'P1PSD',
             'LIU', 'Q6_withQ5Poison_DataAnalysis_2021Feb06to08', 'P1_Parity_Interleave', 'HDF5Data']

file_Q6_NoPoison = 'dhv1637zpl_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg100.hdf5'
Q6_NoPoison = OneStateAnalyze()
Q6_NoPoison.one_state_data_retrieve(data_path, file_Q6_NoPoison, label='NoPoison')

file_Q6_Neg20dBm = 'dhv1642sww_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg20.hdf5'
Q6_Neg20dBm = OneStateAnalyze()
Q6_Neg20dBm.one_state_data_retrieve(data_path, file_Q6_Neg20dBm, label='Neg20dBm')

file_Q6_Neg19dBm = 'dhv1645bsw_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg19.hdf5'
Q6_Neg19dBm = OneStateAnalyze()
Q6_Neg19dBm.one_state_data_retrieve(data_path, file_Q6_Neg19dBm, label='Neg19dBm')

file_Q6_Neg18dBm = 'dhv1652foy_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg18.hdf5'
Q6_Neg18dBm = OneStateAnalyze()
Q6_Neg18dBm.one_state_data_retrieve(data_path, file_Q6_Neg18dBm, label='Neg18dBm')

file_Q6_Neg17dBm = 'dhv1653por_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg17.hdf5'
Q6_Neg17dBm = OneStateAnalyze()
Q6_Neg17dBm.one_state_data_retrieve(data_path, file_Q6_Neg17dBm, label='Neg17dBm')

file_Q6_Neg16dBm = 'dhv1653wxf_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg16.hdf5'
Q6_Neg16dBm = OneStateAnalyze()
Q6_Neg16dBm.one_state_data_retrieve(data_path, file_Q6_Neg16dBm, label='Neg16dBm')

file_Q6_Neg15dBm = 'dhv1654dpl_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg15.hdf5'
Q6_Neg15dBm = OneStateAnalyze()
Q6_Neg15dBm.one_state_data_retrieve(data_path, file_Q6_Neg15dBm, label='Neg15dBm')

file_Q6_Neg14dBm = 'dhv1654hju_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg14.hdf5'
Q6_Neg14dBm = OneStateAnalyze()
Q6_Neg14dBm.one_state_data_retrieve(data_path, file_Q6_Neg14dBm, label='Neg14dBm')

file_Q6_Neg13dBm = 'dhv1654fmw_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg13.hdf5'
Q6_Neg13dBm = OneStateAnalyze()
Q6_Neg13dBm.one_state_data_retrieve(data_path, file_Q6_Neg13dBm, label='Neg13dBm')

file_Q6_Neg12dBm = 'dhv1712oxw_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg12.hdf5'
Q6_Neg12dBm = OneStateAnalyze()
Q6_Neg12dBm.one_state_data_retrieve(data_path, file_Q6_Neg12dBm, label='Neg12dBm')

file_Q6_Neg11dBm = 'dhv1713jwt_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg11.hdf5'
Q6_Neg11dBm = OneStateAnalyze()
Q6_Neg11dBm.one_state_data_retrieve(data_path, file_Q6_Neg11dBm, label='Neg11dBm')

file_Q6_Neg10dBm = 'dhv1713yyu_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg10.hdf5'
Q6_Neg10dBm = OneStateAnalyze()
Q6_Neg10dBm.one_state_data_retrieve(data_path, file_Q6_Neg10dBm, label='Neg10dBm')

file_Q6_Neg9dBm = 'dhv1714ezv_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg9.hdf5'
Q6_Neg9dBm = OneStateAnalyze()
Q6_Neg9dBm.one_state_data_retrieve(data_path, file_Q6_Neg9dBm, label='Neg9dBm')

file_Q6_Neg8dBm = 'dhv1714ija_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg8.hdf5'
Q6_Neg8dBm = OneStateAnalyze()
Q6_Neg8dBm.one_state_data_retrieve(data_path, file_Q6_Neg8dBm, label='Neg8dBm')

file_Q6_Neg7dBm = 'dhv1714jtw_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg7.hdf5'
Q6_Neg7dBm = OneStateAnalyze()
Q6_Neg7dBm.one_state_data_retrieve(data_path, file_Q6_Neg7dBm, label='Neg7dBm')

file_Q6_Neg6dBm = 'dhv1714niq_2021Feb06_P1_Parity_Interleave_Interleave_P1_Neg6.hdf5'
Q6_Neg6dBm = OneStateAnalyze()
Q6_Neg6dBm.one_state_data_retrieve(data_path, file_Q6_Neg6dBm, label='Neg6dBm')

# # python to pandas
r_avg = 10  # rolling average
file_length = 400 * 5
t = Q6_NoPoison.time_data[:file_length]

# Q6_list = [Q6_NoPoison, Q6_Neg20dBm, Q6_Neg19dBm, Q6_Neg18dBm,
#            Q6_Neg17dBm, Q6_Neg16dBm, Q6_Neg15dBm, Q6_Neg14dBm]

# Q6_list = [Q6_Neg13dBm, Q6_Neg12dBm, Q6_Neg11dBm, Q6_Neg10dBm,
#            Q6_Neg9dBm, Q6_Neg8dBm, Q6_Neg7dBm, Q6_Neg6dBm]

# Q6_list = [Q6_NoPoison, Q6_Neg20dBm, Q6_Neg13dBm, Q6_Neg10dBm, Q6_Neg7dBm]

Q6_list = [Q6_NoPoison, Q6_Neg20dBm, Q6_Neg19dBm, Q6_Neg18dBm, Q6_Neg17dBm, Q6_Neg16dBm,
           Q6_Neg15dBm, Q6_Neg14dBm, Q6_Neg13dBm, Q6_Neg12dBm, Q6_Neg11dBm,
           Q6_Neg10dBm, Q6_Neg9dBm, Q6_Neg8dBm, Q6_Neg7dBm, Q6_Neg6dBm]

len_list = len(Q6_list)
p1_max = 0.4
f, ax = plt.subplots(len_list, 2)
ax[0][0].set_title('Q6 P1 Time Trace')
for i, Q6 in enumerate(Q6_list):
    ax[i][0].plot(t, Q6.get_rolling_data_pre(r_avg)[0][:file_length], label=Q6.label)
    ax[i][0].legend(loc=1)
    ax[i][0].grid()
    ax[i][0].set_ylim(0, p1_max)
    ax[i][0].set_ylabel('P1')
    ax[i][0].set_xlabel('file number (5 file = 3.7 min)')

### Start the histogram

bin_list = np.arange(0, p1_max, 0.002)

ax[0][1].set_title('P1 Histogram Scaled Probability')

for i, Q6 in enumerate(Q6_list):
    weights = np.ones_like(Q6.one_state_data) / len(Q6.one_state_data)
    ax[i][1].hist(Q6.one_state_data, bins=bin_list, weights=weights,
                  label=Q6.label)
    ax[i][1].legend(loc=1)
    ax[i][1].grid()
    ax[i][1].set_ylim(0, 0.3)
    ax[i][1].set_xlabel('P1')

plt.show()




