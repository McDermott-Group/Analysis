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

ExptInfo = {
    ## processed data save information
    'Device Name': 'DataAnalysis_test',
    'User': 'LIU',
    'Base Path': r'Z:\mcdermott-group\data\GapEngineer\Nb_GND_Dev06_Trap\Leiden_2020Jul\Debug',
    'Experiment Name': 'one_state',
    'Poison Method': 'Cavity bare resonance',
    'Poison Resonator': 'R5',
    'Measurement Qubit': 'Q6',

    ## matlab data import info:
    'path': 'Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap'
            '/Leiden_2020Jul/Debug/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}',
    'expt_name': 'One_State_Poison_Neg8dBm',
    'date': '09-30-20',
    'files': np.arange(0, 24, 1),
}

class OneState(object):
    def __init__(self, t=10e-6):
        self.t = t
        self.one_state_avg = np.array([])
        self.time = []
        self.average_n1 = 1000  # this is the average number from raw matlab data to datachest
        self.expt_info = {}  # this stores the data source and path, etc

    def add_data_from_matlab(self, filenames,
                             data_type='Single_Shot_Occupations'):
        """
        import the data from Matlab and do first level average
        :param filenames: a list of file names
        :param data_type:
        :return: array of single_shot_occupations
        """
        n = self.average_n1
        for file in filenames:
            data = loadmat(file)
            one_state_list = np.array(data[data_type])
            for i_os in one_state_list:
                self._get_one_state_avg(os=i_os, n=n)
        t_avg = n * self.t
        self.time = (np.arange(0, len(self.one_state_avg), 1)) * t_avg

    def _get_one_state_avg(self, os, n):
        one_state_array = np.asarray(os)
        avg = np.mean(one_state_array.reshape(-1, n), axis=1)
        self.one_state_avg = np.append(self.one_state_avg, avg)

    def averaged_data_to_dataChest(self):
        d = create_datachest_object(self.expt_info)

        time_axis = self.time * 10 ** 3

        d.createDataset(ExptInfo['Experiment Name'],
                        [("time", [len(time_axis)], "float64", "ms")],
                        [("Occupation", [len(time_axis)], "float64", "")]
                        )

        d.addParameter("X Lable", "Time")
        d.addParameter("Y Lable", "Occupation")

        for key in self.expt_info:
            d.addParameter(key, self.expt_info[key])
            
        d.addData([[time_axis, self.one_state_avg]])

def create_datachest_object(expt_info):
    path_information_dict = expt_info
    experiment_path = path_information_dict['Base Path'].split("\\")[3:]
    experiment_path.append(path_information_dict['User'])
    experiment_path.append(path_information_dict['Device Name'])
    experiment_path.append(path_information_dict['Experiment Name'])
    experiment_path.append('HDF5Data')
    d = dataChest(experiment_path)
    return d

def run_matlab_data_to_datachest(ExptInfo):
    os = OneState()

    for key in ExptInfo:
        os.expt_info[key] = ExptInfo[key]

    path = ExptInfo['path']
    expt_name = ExptInfo['expt_name']
    date = ExptInfo['date']
    files = ExptInfo['files']
    filenames = [path.format(date, expt_name, expt_name) + '_{:03d}.mat'.format(i) for i in files]

    os.add_data_from_matlab(filenames)
    os.averaged_data_to_dataChest()

def one_state_data_retrieve():
    """
    get time and occupation data back
    :return:
    """
    # file_info = FilePicker()
    # data_path = file_info.relativePath
    # file_name = file_info.selectedFileName
    file_name = ''
    data_path = ['GapEngineer', 'Nb_GND_Dev06_Trap', 'Leiden_2020Jul', 'Debug',
                 'LIU', 'DataAnalysis', 'one_state', 'HDF5Data']
    # data_set = 'Q4Neg6dBm'
    data_set = 'Q4Neg8dBm'
    # data_set = 'Q4Off'
    # data_set = 'Q6Neg6dBm'
    # data_set = 'Q6Off'
    if data_set == 'Q6Off':
        file_name = 'Q6_PoisonOff_0_2499.hdf5'
    elif data_set == 'Q6Neg6dBm':
        file_name = 'Q6_PoisonNeg6dBm_0_2499.hdf5'
    elif data_set == 'Q4Neg6dBm':
        file_name = 'Q4_PoisonNeg6dBm_0_4370.hdf5'
    elif data_set == 'Q4Neg8dBm':
        file_name = 'Q4_PoisonNeg8dBm_0_2960.hdf5'
    elif data_set == 'Q4Off':
        file_name = 'Q4_PoisonOff_0_3779.hdf5'
    d = dataChest(data_path)
    d.openDataset(file_name)
    variable_info = d.getVariables()
    independent_variable_info = variable_info[0]
    dependent_variable_info = variable_info[1]

    occupation_data = d.getData(variablesList=[dependent_variable_info[0][0]])
    occupation_data = np.asarray(occupation_data).flatten()

    time_data = d.getData(variablesList=[independent_variable_info[0][0]])
    time_data = np.asarray(time_data).flatten()

    return [time_data, occupation_data, data_set]


run_matlab_data_to_datachest(ExptInfo)

# t_data, o_data, name = one_state_data_retrieve()[0], one_state_data_retrieve()[1], one_state_data_retrieve()[2]
# t_data, o_data = np.full(100000, 0.5, dtype=np.float), np.full(100000, 0.5, dtype=np.float)

# plt.plot(t_data, o_data)
# n = 100 # before this, each point is 10*1000 data, 10ms, in lab time, that is 0.1728 second per point
# avg_n = np.mean(o_data.reshape(-1, n), axis=1)

# python to pandas
# o_data_p = pd.Series(o_data)
# o_data_rolling_mean = o_data_p.rolling(n).mean()
# o_data_rolling_std = o_data_p.rolling(n).std()

# pandas to numpy
# o_data_avg_np = o_data_avg.to_numpy

# print type(o_data_avg_np)
# fig = plt.figure(1)
# plt.plot(o_data_rolling_mean)
# plt.plot(o_data_rolling_std + 1)
# plt.ylim(0, 1.4)
# plt.title(name+' Rolling Mean and std')
#
# fig = plt.figure(2)
# bin_list = np.arange(0, 1, 0.02)
# plt.hist(o_data_rolling_mean, bins=bin_list)
# plt.title(name+' Rolling Mean histgram')
# plt.show()

# o_data_p.plot(style='k--')

# print (type(o_data_p.rolling(100).mean()))






