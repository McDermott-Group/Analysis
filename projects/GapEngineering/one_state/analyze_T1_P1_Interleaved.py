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
        self.one_state_data_pre = np.array([])
        self.one_state_data_filtered = np.array([])
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
        # print('dependent_variable_info[0]=', dependent_variable_info[0])

        occupation_data_pre = d.getData(
            variablesList=[dependent_variable_info[0][0]])
        occupation_data_pre = np.asarray(occupation_data_pre).flatten()

        occupation_data = d.getData(
            variablesList=[dependent_variable_info[1][0]])
        occupation_data = np.asarray(occupation_data).flatten()

        occupation_data_filtered = d.getData(
            variablesList=[dependent_variable_info[2][0]])
        occupation_data_filtered = np.asarray(occupation_data_filtered).flatten()

        time_data = d.getData(variablesList=[independent_variable_info[0][0]])
        time_data = np.asarray(time_data).flatten()

        self.one_state_data_pre = occupation_data_pre
        self.one_state_data = occupation_data
        self.one_state_data_filtered = occupation_data_filtered
        self.time_data = time_data

    def get_rolling_data_pre(self, rolling_avg):
        r_avg = rolling_avg
        data = self.one_state_data_pre
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]

    def get_rolling_data(self, rolling_avg):
        r_avg = rolling_avg
        data = self.one_state_data
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]

    def get_rolling_data_filtered(self, rolling_avg):
        r_avg = rolling_avg
        data = self.one_state_data_filtered
        d = pd.Series(data)
        return [d.rolling(r_avg).mean(), d.rolling(r_avg).std()]



data_path = ['GapEngineer', 'Nb_GND_Dev06_Trap', 'Leiden_2020Jul', 'Debug',
             'LIU', 'DataAnalysis', 'one_state', 'HDF5Data']

# file_I = 'Q4_2020Oct12_one_state_Interleave_T1_Off.hdf5'
# file_X = 'Q4_2020Oct12_one_state_Interleave_T1_Off.hdf5'
# file_I = 'Q4_2020Oct12_one_state_Interleave_One_State_Off_NoFiltered.hdf5'

file_X1 = 'Q4_one_state_Interleave_T1_Off_2020Oct12.hdf5'
file_I1 = 'Q4_one_state_Interleave_One_State_Off_2020Oct12.hdf5'
file_X2 = 'Q4_2020Oct12_one_state_Interleave_T1_Neg10.hdf5'
file_I2 = 'Q4_2020Oct12_one_state_Interleave_One_State_Neg10.hdf5'

os_X1 = OneStateAnalyze()
os_X1.one_state_data_retrieve(data_path, file_X1)
os_X2 = OneStateAnalyze()
os_X2.one_state_data_retrieve(data_path, file_X2)
os_I1 = OneStateAnalyze()
os_I1.one_state_data_retrieve(data_path, file_I1)
os_I2 = OneStateAnalyze()
os_I2.one_state_data_retrieve(data_path, file_I2)

# python to pandas
r_avg = 100  # rolling average
t = os_X1.time_data
# t = t*13140/(1*t[-1])
t = t*16500/(1*t[-1])
# files = os_X.files
n_file = 10 # file number step

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
# ax2 = ax1.twiny()

# ax1.plot(t, os_I.get_rolling_data_pre(r_avg)[0], label='Clean Dirty Monitor')
ax1.plot(t, os_X1.get_rolling_data_filtered(r_avg)[0], label='T1 Filtered Off')
ax1.plot(t, os_X2.get_rolling_data_filtered(r_avg)[0], label='T1 Filtered -10 dBm')
ax1.plot(t, os_I1.get_rolling_data_pre(r_avg)[0], label='Dirty Clean Monitor Off')
ax1.plot(t, os_I2.get_rolling_data_pre(r_avg)[0], label='Dirty Clean Monitor -10 dBm')
ax1.legend()
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('Occupation')
ax1.set_title('Rolling Mean Poison')
# ax1.set_title('Rolling Mean -10dBm Poison ~Al Gap')
ax1.set_ylim(0, 1.0)

"""Set second x axis so the file number is obvious to track"""
# new_tick_locations = t[::n_file*len(t)/(len(files))]
# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticks(new_tick_locations)
# ax2.set_xticklabels(files[::n_file])
# ax2.set_xlabel('File Number')
# ax2.tick_params(axis="x",direction="in", pad=-22)

fig = plt.figure(2)

ax1 = fig.add_subplot(111)
ax1.plot(t, [-4.5*(np.log(p))**(-1) for p in os_X1.get_rolling_data_filtered(r_avg)[0]], label='T1 Filtered Off')
ax1.plot(t, [-4.5*(np.log(p))**(-1) for p in os_X2.get_rolling_data_filtered(r_avg)[0]], label='T1 Filtered -10 dBm')
ax1.legend()
ax1.set_xlabel('Time (sec)')
ax1.set_ylabel('T1 (us)')
ax1.set_title('T1 Conversion')

plt.show()

# one point is 50 * 2000 * 10 us = 1 second

#
# fig = plt.figure(2)
# bin_list = np.arange(0, 0.4, 0.001)
# avg = np.mean(os_I.one_state_data, axis=1)
# plt.hist(os_I.one_state_data_pre, bins=bin_list)
# plt.hist(os_I.get_rolling_data_pre(r_avg)[0], bins=bin_list)
# plt.hist(a, bins=bin_list)
# plt.title('Hist Poison Off')
# plt.title('Histm Poison -10dBm Poison')
plt.show()

# o_data_p.plot(style='k--')

# print (type(o_data_p.rolling(100).mean()))






