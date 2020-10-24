"""
Under development by Chuan-Hong Liu
reference to Alex's adapative spectroscopy
TO DO:
"""
import pandas as pd
import numpy as np
import datetime
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
    'Poison Method': 'No Poison',
    # 'Poison Method': 'Cavity bare resonance',
    'Poison Resonator': 'R5',
    'Measurement Qubit': 'Q4',

    ## matlab data import info:
    'path': 'Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap'
            '/Leiden_2020Jul/Debug/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}',
    # 'expt_name': 'One_State_X_Total50us_RO2us_Idle3us',
    # 'expt_name': 'Interleave_One_State_Off',
    # 'expt_name': 'Interleave_T1_Off',
    # 'expt_name': 'Interleave_T1_Neg10',
    # 'expt_name': 'Interleave_One_State_Neg10',
    # 'expt_name': 'Interleave_I_Off',
    # 'expt_name': 'Interleave_X_Off',
    'expt_name': 'Interleave_I_Off_1',
    'date': '10-14-20',
    'files': np.arange(0, 500, 1),
}


class OneState(object):
    def __init__(self, t=10e-6):
        self.t = t
        self.one_state_avg_pre = np.array([])
        self.one_state_avg = np.array([])
        self.one_state_avg_filtered = np.array([])
        self.time = []
        self.average_n1 = 10000  # this is the average number from raw matlab data to datachest
        self.expt_info = {}  # this stores the data source and path, etc

    def add_data_from_matlab(self, filenames,
                             data_type='Single_Shot_Occupations',
                             data_type_pre='Single_Shot_Occupations_Pre',
                             data_type_filtered='Single_Shot_Occupation_Filtered'):
        """
        import the data from Matlab and do first level average
        :param filenames: a list of file names
        :param data_type:
        :return: array of single_shot_occupations
        """
        n = self.average_n1
        for f in filenames:
            data = loadmat(f)
            one_state_filtered_list = np.array(data[data_type_filtered])
            one_state_pre_list = np.array(data[data_type_pre])
            one_state_list = np.array(data[data_type])
            for i_os_filtered in one_state_filtered_list:
                avg_filtered = self._get_one_state_avg(os=i_os_filtered, n=1)
                self.one_state_avg_filtered = np.append(
                    self.one_state_avg_filtered, avg_filtered)
            for i_os_pre in one_state_pre_list:
                avg_pre = self._get_one_state_avg(os=i_os_pre, n=n)
                self.one_state_avg_pre = np.append(self.one_state_avg_pre,
                                                   avg_pre)
            for i_os in one_state_list:
                avg = self._get_one_state_avg(os=i_os, n=n)
                self.one_state_avg = np.append(self.one_state_avg, avg)
        t_avg = n * self.t
        self.time = (np.arange(0, len(self.one_state_avg), 1)) * t_avg

    def _get_one_state_avg(self, os, n):
        one_state_array = np.asarray(os)
        avg = np.mean(one_state_array.reshape(-1, n), axis=1)
        return avg

    def averaged_data_to_dataChest(self):
        d = create_datachest_object(self.expt_info)
        time_axis = self.time * 10 ** 3
        date = datetime.strptime(ExptInfo['date'], '%m-%d-%y')
        date = date.strftime('%Y%b%d')
        # d.createDataset(date+'_'+ExptInfo['Experiment Name']+'_'+ExptInfo['expt_name'],
        #                 [("time", [len(time_axis)], "float64", "ms")],
        #                 [("Occupation Pre", [len(time_axis)], "float64", ""),
        #                  ("Occupation", [len(time_axis)], "float64", ""),
        #                  ("Occupation Filtered", [len(time_axis)], "float64",
        #                   "")]
        #                 )
        d.createDataset(date+'_'+ExptInfo['Experiment Name']+'_'+ExptInfo['expt_name'],
                        [("time", [len(time_axis)], "float64", "ms")],
                        [("Occupation Pre", [len(time_axis)], "float64", ""),
                         ("Occupation", [len(time_axis)], "float64", "")]
                        )
        d.addParameter("X Lable", "Time")
        d.addParameter("Y Lable", "Occupation")

        for key in self.expt_info:
            d.addParameter(key, self.expt_info[key])

        # d.addData([[time_axis, self.one_state_avg_pre, self.one_state_avg,
                    # self.one_state_avg_filtered]])

        d.addData([[time_axis, self.one_state_avg_pre, self.one_state_avg]])

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
    filenames = [
        path.format(date, expt_name, expt_name) + '_{:03d}.mat'.format(i) for i
        in files]

    os.add_data_from_matlab(filenames)
    os.averaged_data_to_dataChest()


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

    def one_state_data_retrieve(self, data_path, data_set):
        self.data_set = data_set

        if data_set == 'Q6Off':
            file_name = 'Q6_PoisonOff_0_2499.hdf5'
        elif data_set == 'Q6Neg10dBm':
            file_name = 'Q6_PoisonNeg10dBm_0_2499.hdf5'
        elif data_set == 'Q6Neg8dBm':
            file_name = 'Q6_PoisonNeg8dBm_0_2499.hdf5'
        elif data_set == 'Q6Neg7dBm':
            file_name = 'Q6_PoisonNeg7dBm_0_2499.hdf5'
        elif data_set == 'Q6Neg6dBm':
            file_name = 'Q6_PoisonNeg6dBm_0_2499.hdf5'
        elif data_set == 'Q6Neg4dBm':
            file_name = 'Q6_PoisonNeg4dBm_0_2499.hdf5'

        elif data_set == 'Q4Neg6dBm':
            file_name = 'Q4_PoisonNeg6dBm_0_4370.hdf5'
        elif data_set == 'Q4Neg8dBm':
            file_name = 'Q4_PoisonNeg8dBm_0_2960.hdf5'
        elif data_set == 'Q4Off':
            file_name = 'Q4_PoisonOff_0_3779.hdf5'

        elif data_set == 'Interleave_T1':
            file_name = 'Q4_Off_Interleave_T1_Cal_0_700.hdf5'

        elif data_set == 'Interleave_One_State':    #13140 seconds
            file_name = 'Q4_Off_Interleave_One_State_Cal_0_700.hdf5'

        elif data_set == 'Interleave_T1_Oct_11':
            file_name = 'Q4_Off_Interleave_T1_Cal_382_1116.hdf5'

        elif data_set == 'Interleave_One_State_Oct_11':
            file_name = 'Q4_Off_Interleave_One_State_Cal_382_1116.hdf5'

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


run_matlab_data_to_datachest(ExptInfo)
