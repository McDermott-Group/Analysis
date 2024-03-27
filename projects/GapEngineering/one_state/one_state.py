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
from noiselib import loadmat, loadmat_ExptVars
from dataChest import dataChest
from QPTunneling_LIU import QPTunneling_Liu, plotMultiFittedPSD, QPTunneling_Wilen
from dataChestUtils import FilePicker
# import Markov_Python2.analyze_QPTunneling_pomegranate as pome
from Markov_Python2.analyze_QPTunneling_pomegranate import *

ExptInfo = {
    ## processed data save information
    'Device Name': 'Q4_DataAnalysis',
    'User': 'LIU',
    # 'Base Path': r'Z:\mcdermott-group\data\GapEngineer\Nb_GND_Dev06_Trap\Leiden_2021Jan\P1PSD',
    'Base Path': r'Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2',
    'Experiment Name': 'P1_Parity_Interleave',
    'Poison Method': 'BB 305mK',
    # 'Poison Resonator': 'BB',
    'Measurement Qubit': 'Q4',

    ## matlab data import info:
    # 'path': 'Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap'
    #         '/Leiden_2021Jan/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}',
    'path': 'Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}',
    'expt_name_p1': 'Interleave_P1_Neg305',
    'expt_name_parity_switch': 'Interleave_PSD_Neg305',
    'Comment': '14 us T1',
    'date': '04-16-21',
    # 'date': '2021Feb06',
    'files': np.arange(0, 200, 1),
}


class OneState(object):
    def __init__(self, t=50e-6):
        self.t = t
        self.p1_sso_avg_pre = np.array([])
        self.p1_sso_avg = np.array([])
        self.p1_weighted = np.array([])
        self.parity_trace_avg = np.array([])
        self.parity_time_trace_list = []    # HMM Check
        self.parity_trace_HMM_list = []    # HMM Check
        self.parity_jump = np.array([])
        self.parity_jump_HMM = np.array([])
        self.stripped_charge_bias = np.array([])
        self.freq_detuning = np.array([])
        self.charge_parity_sso_avg = np.array([])
        self.charge_parity_sso_avg_pre = np.array([])
        self.time = []
        self.file_axis = []
        self.average_n1 = 10000  # this is the average number from raw matlab data to datachest
        self.expt_info = {}  # this stores the data source and path, etc

    def add_p1_data_from_matlab(self, filenames,
                             data_type='Single_Shot_Occupations',
                                data_type_p1_weighted='Weighted_Occupation',
                                data_type_stripped_charge='Stripped_Charge_Bias'):
        """
        import the data from Matlab and do first level average
        :param filenames: a list of file names
        :param data_type:
        :return: array of single_shot_occupations
        """
        n = self.average_n1

        f_l = filenames[0]
        sample_rate = loadmat_ExptVars(f_l)['Total_Time']
        sample_rate = sample_rate * 10**(-6)
        self.t = sample_rate

        for f in filenames:
            data = loadmat(f)
            one_state_list = np.array(data[data_type])
            weighted_occupation_list = np.array(data[data_type_p1_weighted])
            stripped_charge_bias_list = np.array(data[data_type_stripped_charge])
            for i_os in one_state_list:
                avg = self._get_one_state_avg(os=i_os, n=n)
                self.p1_sso_avg = np.append(self.p1_sso_avg, avg)

            for charge_bias in stripped_charge_bias_list:
                self.stripped_charge_bias = np.append((self.stripped_charge_bias), charge_bias)

            for WO in weighted_occupation_list:
                self.p1_weighted = np.append((self.p1_weighted), WO)

        t_avg = n * self.t
        self.time = (np.arange(0, len(self.p1_sso_avg), 1)) * t_avg

        self.file_axis = self._get_interpolate_file_axis()

        # print('self.stripped_charge_bias=', self.stripped_charge_bias)

        freq_detuning = [0.1*np.abs(1.25*np.cos(2*np.pi*(-x0))) for x0 in self.stripped_charge_bias]
        # stripped_charge_bias = self.stripped_charge_bias
        # freq_detuning = self.stripped_charge_bias
        # for i in range(len(stripped_charge_bias)-1):
        #     x0 = stripped_charge_bias[i]
        #     x1 = stripped_charge_bias[i+1]
        #     freq_detuning[i] = 0.1*np.abs(1.25*np.cos(2*np.pi*(x1-x0)))
        self.freq_detuning = freq_detuning


    def _get_interpolate_file_axis(self):
        files = self.expt_info['files']
        time_axis = self.time
        n_inter = len(time_axis)/len(files)
        file_axis = []
        # print('n_inter=', n_inter)
        for f in files:
            for n in range(n_inter):
                file_axis.append(1.0*f+1.0*n/n_inter)
        # print(file_axis)
        return file_axis


    def add_parity_data_from_matlab(self, filenames,
                                    data_type_parity_trace='Charge_Parity_Trace',
                                    data_type_parity_trace_jump_count='Charge_Parity_Jump_Counts'):
        """
        import the data from Matlab and do first level average
        :param filenames: a list of file names
        :param data_type:
        :return: array of single_shot_occupations
        """
        n = self.average_n1
        for i, f in enumerate(filenames):

            data = loadmat(f)
            parity_trace_list = np.array(data[data_type_parity_trace])
            parity_jump_count_list = np.array(data[data_type_parity_trace_jump_count])
            for i_os in parity_trace_list:

                """plot"""
                # print('i_os=', i_os[:100])
                # i_os is an array, i_os = [-1.0, -1.0, 1.0, 1.0, ...]
                avg = self._get_one_state_avg(os=i_os, n=n)
                self.parity_trace_avg = np.append(self.parity_trace_avg, avg)
                ### Add parity data to HMM
                i_os_HMM = self._apply_HMM(i_os)
                jump_count = transitions_count(i_os_HMM)

                self.parity_trace_HMM_list.append(i_os_HMM)
                self.parity_jump_HMM = np.append(self.parity_jump_HMM, jump_count)

            for i_os in parity_jump_count_list:
                avg = self._get_one_state_avg(os=i_os, n=1)
                self.parity_jump = np.append((self.parity_jump), avg)


    def _get_one_state_avg(self, os, n):
        one_state_array = np.asarray(os)
        avg = np.mean(one_state_array.reshape(-1, n), axis=1)
        return avg

    def _apply_HMM(self, os, apply=False):
        random_seed = 2
        os_0and1 = self._parity_offset_convert(os)
        if apply:
            recovered_signal_BW_Result = observed_to_recovered_signal_BW(os_0and1, seed=random_seed)
            os_HMM = recovered_signal_BW_Result[0]
            print(recovered_signal_BW_Result[1:])
        else:
            os_HMM = os_0and1
        return os_HMM

    def _parity_offset_convert(self, os):
        """
        Convert the -1, 1 parity value to 0, and 1 for the HMM code already written
        :param os:
        :return:
        """
        os_0and1 = [int(i * 0.5 + 0.5) for i in os]
        return os_0and1

    def averaged_data_to_dataChest(self):
        d = create_datachest_object(self.expt_info)
        time_axis = self.time * 10 ** 3
        file_axis = self.file_axis
        try:
            date = datetime.strptime(ExptInfo['date'], '%m-%d-%y')
            date = date.strftime('%Y%b%d')
        except ValueError:
            date = ExptInfo['date']

        # try:
        #     fit_data = fit_charge(DACArray, data, e_period)
        # except RuntimeError:
        #     print('Rerun charge stabilize Starts')
        #     fit_data = self.charge_stabilize()

        d.createDataset(date+'_'+ExptInfo['Experiment Name']+'_'+ExptInfo['expt_name_p1'],
                        [("file_number", [len(time_axis)], "float64", "")],
                        [("P1 SSO Avg", [len(time_axis)], "float64", ""),
                         ("Parity Trace Avg", [len(time_axis)], "float64", ""),
                         ("Parity Jump Counts", [len(time_axis)], "float64", ""),
                         ("Parity Jump Counts HMM", [len(time_axis)], "float64", ""),
                         ("Stripped Charge Bias", [len(time_axis)], "float64", "10 e"),
                         ("Freq Detuning", [len(time_axis)], "float64", "10 MHz"),
                         ("P1 Weighted", [len(time_axis)], "float64", "")])
        d.addParameter("X Lable", "Time")
        d.addParameter("Y Lable", "Occupation")

        for key in self.expt_info:
            d.addParameter(key, self.expt_info[key])

        # d.addData([[time_axis, self.p1_sso_avg_pre, self.p1_sso_avg, self.parity_trace_avg,
        #             self.charge_parity_sso_avg, self.charge_parity_sso_avg_pre]])
        # d.addData([[time_axis, self.p1_sso_avg, self.parity_trace_avg, (self.parity_jump)*0.0001]])
        # d.addData([[time_axis, self.p1_sso_avg, self.parity_trace_avg,
        #             (self.parity_jump)*0.0001, self.parity_jump_HMM*0.001]])
        d.addData([[self.file_axis, self.p1_sso_avg, self.parity_trace_avg,
                    (self.parity_jump)*0.0001, self.parity_jump_HMM*0.001,
                    self.stripped_charge_bias * 0.1, self.freq_detuning, self.p1_weighted]])
                   # self.stripped_charge_bias])
        # d.addData([[time_axis, self.p1_sso_avg]])

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
    date = ExptInfo['date']
    files = ExptInfo['files']
    expt_name_p1 = ExptInfo['expt_name_p1']
    filenames_p1 = [
        path.format(date, expt_name_p1, expt_name_p1) + '_{:03d}.mat'.format(i) for i
        in files]
    os.add_p1_data_from_matlab(filenames_p1)

    expt_name_parity_switch = ExptInfo['expt_name_parity_switch']
    filenames_parity_switch = [
        path.format(date, expt_name_parity_switch, expt_name_parity_switch) + '_{:03d}.mat'.format(i) for i
        in files]
    os.add_parity_data_from_matlab(filenames_parity_switch)

    os.averaged_data_to_dataChest()

def run_matlab_data_to_datachest_HMM(ExptInfo):
    os = OneState()
    for key in ExptInfo:
        os.expt_info[key] = ExptInfo[key]

    path = ExptInfo['path']
    date = ExptInfo['date']
    files = ExptInfo['files']
    expt_name_p1 = ExptInfo['expt_name_p1']
    filenames_p1 = [
        path.format(date, expt_name_p1,
                    expt_name_p1) + '_{:03d}.mat'.format(i) for i
        in files]
    os.add_p1_data_from_matlab(filenames_p1)

    expt_name_parity_switch = ExptInfo['expt_name_parity_switch']
    filenames_parity_switch = [
        path.format(date, expt_name_parity_switch,
                    expt_name_parity_switch) + '_{:03d}.mat'.format(i) for
        i
        in files]
    os.add_parity_data_from_matlab(filenames_parity_switch)

    os.averaged_data_to_dataChest()

    ### HMM debug starts

    # Recovered_Signal_BW = os.parity_trace_HMM_list[0]
    # fig = plt.figure(figsize=(12, 4))
    # plt.plot(asarray(Recovered_Signal_BW), 'o-', label=r"{} Recovered Transitions".format(transitions_count(Recovered_Signal_BW)))
    # plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
    # plt.show()
    ### HMM debug ends

run_matlab_data_to_datachest_HMM(ExptInfo)

def run_matlab_data_to_datachest_plotParity(ExptInfo):
    os = OneState()
    for key in ExptInfo:
        os.expt_info[key] = ExptInfo[key]

    path = ExptInfo['path']
    date = ExptInfo['date']
    files = ExptInfo['files']
    expt_name_p1 = ExptInfo['expt_name_p1']
    filenames_p1 = [
        path.format(date, expt_name_p1,
                    expt_name_p1) + '_{:03d}.mat'.format(i) for i
        in files]
    os.add_p1_data_from_matlab(filenames_p1)

    expt_name_parity_switch = ExptInfo['expt_name_parity_switch']
    filenames_parity_switch = [
        path.format(date, expt_name_parity_switch,
                    expt_name_parity_switch) + '_{:03d}.mat'.format(i) for
        i
        in files]
    os.add_parity_data_from_matlab(filenames_parity_switch)

run_matlab_data_to_datachest_plotParity(ExptInfo)

