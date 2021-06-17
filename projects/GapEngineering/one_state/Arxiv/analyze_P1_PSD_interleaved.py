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
        self.parity_avg_data = np.array([])
        self.parity_jump_data = np.array([])
        self.time_data = np.asarray([])
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
        # print('dependent_variable_info[0]=', dependent_variable_info[0])

        occupation_data = d.getData(
            variablesList=[dependent_variable_info[0][0]])
        occupation_data = np.asarray(occupation_data).flatten()

        parity_avg_data = d.getData(
            variablesList=[dependent_variable_info[1][0]])
        parity_avg_data = np.asarray(parity_avg_data).flatten()

        parity_jump_data = d.getData(
            variablesList=[dependent_variable_info[2][0]])
        parity_jump_data = np.asarray(parity_jump_data).flatten()

        time_data = d.getData(variablesList=[independent_variable_info[0][0]])
        time_data = np.asarray(time_data).flatten()

        self.one_state_data = occupation_data
        self.parity_avg_data = parity_avg_data
        self.parity_jump_data = parity_jump_data
        self.time_data = time_data

    def get_rolling_one_state_data(self, rolling_avg):
        r_avg = rolling_avg
        data = self.one_state_data
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]

    def get_rolling_parity_avg_data(self, rolling_avg):
        r_avg = rolling_avg
        data = self.parity_avg_data
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]

    def get_rolling_parity_jump_data(self, rolling_avg):
        r_avg = rolling_avg
        data = self.parity_jump_data
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]

data_path = ['GapEngineer', 'Nb_GND_Dev06_Trap', 'Leiden_2020Jul', 'Debug',
             'LIU', 'DataAnalysis', 'one_state', 'HDF5Data']

file_I1 = 'Q4_2020Oct17_P1_Parity_Interleave_Interleave_P1.hdf5'

os_I1 = OneStateAnalyze()
os_I1.one_state_data_retrieve(data_path, file_I1)
# python to pandas
r_avg = 50  # rolling average
t = os_I1.time_data
t = t*22000/(1*t[-1])
files = os_I1.files
n_file = 10 # file number step

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
# ax2 = ax1.twiny()

# ax1.plot(t, os_I1.one_state_data, label='Dirty Clean Monitor Off')
ax1.plot(t, os_I1.get_rolling_one_state_data(r_avg)[0], label='P1')
ax1.plot(t, os_I1.get_rolling_parity_avg_data(r_avg)[0], label='Parity Avg')
ax1.plot(t, os_I1.get_rolling_parity_jump_data(r_avg)[0], label='Parity Jump Rate Normalized')
ax1.legend()
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('Occupation')
ax1.set_title('Rolling Mean')
# ax1.set_title('Rolling Mean -10dBm Poison ~Al Gap')
# ax1.set_ylim(0, 0.2)

"""Set second x axis so the file number is obvious to track"""
# new_tick_locations = t[::n_file*len(t)/(len(files))]
# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticks(new_tick_locations)
# ax2.set_xticklabels(files[::n_file])
# ax2.set_xlabel('File Number')
# ax2.tick_params(axis="x", direction="in", pad=-22)

# fig = plt.figure(2)
#
# ax1 = fig.add_subplot(111)
# ax1.plot(t, [-4.5*(np.log(p))**(-1) for p in os_X1.get_rolling_data_filtered(r_avg)[0]], label='T1 Filtered Off')
# ax1.plot(t, [-4.5*(np.log(p))**(-1) for p in os_X2.get_rolling_data_filtered(r_avg)[0]], label='T1 Filtered -10 dBm')
# ax1.legend()
# ax1.set_xlabel('Time (sec)')
# ax1.set_ylabel('T1 (us)')
# ax1.set_title('T1 Conversion')

plt.show()






