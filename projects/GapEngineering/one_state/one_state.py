"""
Under development by Chuan-Hong Liu
reference to Alex's adapative spectroscopy
TO DO:
1. To see if I can retrieve the data
2. make the data saving with more detailed information, the original data location
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from noiselib import loadmat

from dataChest import dataChest
from dataChestUtils import FilePicker

ExptInfo = {
    'Device Name': 'DataAnalysis',
    'User': 'LIU',
    'Base Path': r'Z:\mcdermott-group\data\GapEngineer\Nb_GND_Dev06_Trap\Leiden_2020Jul\Debug',
    'Experiment Name': 'one_state',
    'Poison Method': 'Cavity bare resonance',
    'Poison Resonator': 'R5',
    'Measurement Qubit': 'Q6',
    'Comments': 'hard',

    ## data import info:
    'path': 'Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap'
            '/Leiden_2020Jul/Debug/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}',
    'exp_name': 'One_State_Poison_Neg80dBm8GHz',
    'date': '09-29-20',
    'files': np.arange(0, 3780, 1),
}


class one_state(object):
    def __init__(self, t=10e-6):
        self.t = t
        self.one_state_avg = np.array([])
        self.time = []

    def add_data_from_matlab(self, filenames,
                             data_type='Single_Shot_Occupations'):
        """

        :param filenames: a list of file names
        :param data_type:
        :return: array of single_shot_occupations
        """
        n = 1000
        for file in filenames:
            data = loadmat(file)
            one_state_list = np.array(data[data_type])
            for os in one_state_list:
                self._get_one_state_avg(os=os, n=n)
        t_avg = n * self.t
        self.time = (np.arange(0, len(self.one_state_avg), 1)) * t_avg

    def _get_one_state_avg(self, os, n):
        one_state_array = np.asarray(os)
        avg = np.mean(one_state_array.reshape(-1, n), axis=1)
        self.one_state_avg = np.append(self.one_state_avg, avg)

def create_datachest_object(Exp):
    path_information_dict = Exp
    experiment_path = path_information_dict['Base Path'].split("\\")[3:]
    experiment_path.append(path_information_dict['User'])
    experiment_path.append(path_information_dict['Device Name'])
    experiment_path.append(path_information_dict['Experiment Name'])
    experiment_path.append('HDF5Data')
    d = dataChest(experiment_path)
    return d

def python_to_dataChest(Exp, one_state):
    d = create_datachest_object(Exp)
    timeAxis = (one_state.time)*10**3
    d.createDataset(ExptInfo['Experiment Name'],
                    [("time", [len(timeAxis)], "float64", "ms")],
                    [("Occupation", [len(timeAxis)], "float64", "")]
                    )
    d.addParameter("X Lable", "Time")
    d.addParameter("Y Lable", "Occupation")
    d.addParameter("Poison Method", Exp["Poison Method"])
    d.addParameter("Poison Resonator", Exp["Poison Resonator"])
    d.addParameter("Measurement Qubit", Exp["Measurement Qubit"])

    d.addParameter("Import data path", Exp["path"])
    d.addParameter("exp_name", Exp["exp_name"])
    d.addParameter("date", Exp["date"])
    d.addParameter("files", Exp["files"])

    d.addData([[timeAxis, one_state.one_state_avg]])

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


## import data
# path = ExptInfo['path']
# exp_name = ExptInfo['exp_name']
# date = ExptInfo['date']
# files = ExptInfo['files']
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os = one_state()
# os.add_data_from_matlab(filenames)
# python_to_dataChest(ExptInfo, os)

t_data, o_data, name = one_state_data_retrieve()[0], one_state_data_retrieve()[1], one_state_data_retrieve()[2]
# t_data, o_data = np.full(100000, 0.5, dtype=np.float), np.full(100000, 0.5, dtype=np.float)

# plt.plot(t_data, o_data)
n = 100 # before this, each point is 10*1000 data, 10ms, in lab time, that is 0.1728 second per point
avg_n = np.mean(o_data.reshape(-1, n), axis=1)

# python to pandas
o_data_p = pd.Series(o_data)
o_data_rolling_mean = o_data_p.rolling(n).mean()
o_data_rolling_std = o_data_p.rolling(n).std()

# pandas to numpy
# o_data_avg_np = o_data_avg.to_numpy

# print type(o_data_avg_np)
fig = plt.figure(1)
plt.plot(o_data_rolling_mean)
plt.plot(o_data_rolling_std + 1)
plt.ylim(0, 1.4)
plt.title(name+' Rolling Mean and std')

fig = plt.figure(2)
bin_list = np.arange(0, 1, 0.02)
plt.hist(o_data_rolling_mean, bins=bin_list)
plt.title(name+' Rolling Mean histgram')
plt.show()


# o_data_p.plot(style='k--')

# print (type(o_data_p.rolling(100).mean()))






