"""
Under development by Chuan-Hong Liu
reference to Alex's adapative spectroscopy
TO DO:
"""
import numpy as np
import datetime
import matplotlib.pyplot as plt
from datetime import datetime
from noiselib import loadmat, loadmat_ExptVars
from dataChest import dataChest
from Markov_Python2.analyze_QPTunneling_pomegranate import *
import matplotlib.pyplot as plt

ExptInfo = {
    ## processed data save information
    'Device Name': 'DataAnalysis_test',
    'User': 'LIU',
    'Base Path': r'Z:\mcdermott-group\data\GapEngineer\Nb_GND_Dev06_Trap\Leiden_2020Jul\P1PSD',
    # 'Base Path': r'Z:\mcdermott-group\data\GapEngineer\Nb_GND_Dev06_Trap\Leiden_2020Jul\Debug',
    'Experiment Name': 'P1_Parity_Interleave',
    'Poison Method': 'Bare Cavity',
    'Poison Resonator': 'None',
    'Measurement Qubit': 'Q4',

    ## matlab data import info:
    # 'path': 'Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap'
    #         '/Leiden_2020Jul/P1PSD/LIU/Q6_withQ5Poison/{}/{}/MATLABData/{}',
    'path': 'Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap'
            '/Leiden_2020Jul/P1PSD/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}',
    # 'path': 'Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap'
    #         '/Leiden_2020Jul/P1PSD/LIU/Q3_withQ2Poison/{}/{}/MATLABData/{}',
    # 'expt_name_p1': 'Tenus_Interleave_P1_Neg100',
    # 'date': '01-07-21',
    # 'files': np.arange(40, 60, 1), # Clean
    # 'files': np.arange(65, 75, 1),  # Dirty
    # 'expt_name_p1': 'Tenus_Interleave_P1_Neg13',
    # 'date': '01-07-21',
    # 'files': np.arange(0, 5, 1), # Clean

    # 'expt_name_p1': 'Interleave_P1_Neg13',
    # 'date': '12-31-20',
    # 'files': np.arange(50, 60, 1),  # Clean
    # 'files': np.arange(65, 75, 1),  # Dirty

    # 'expt_name_p1': 'Interleave_P1_Neg7',
    # 'date': '01-01-21',
    # # 'files': np.arange(23, 30, 1),  # Clean
    # 'files': np.arange(50, 60, 1),  # Dirty

    # 'expt_name_p1': 'Interleave_P1_Neg6',
    # 'date': '01-01-21',
    # 'files': np.arange(130, 140, 1), # Clean
    # 'files': np.arange(0, 10, 1),  # Dirty
    # 'expt_name_parity_switch': 'Interleave_PSD_Neg10',
    # 'Comment': '10 us T1, No Actual Reset',

    # 'expt_name_p1': 'Twentyus_Interleave_P1_Neg100',
    # 'date': '01-08-21',
    # 'files': np.arange(0, 10, 1), # Clean
    # 'files': np.arange(10, 20, 1),  # Dirty

    # 'expt_name_p1': 'Twentyus_Interleave_P1_Neg13',
    # 'date': '01-08-21',
    # 'files': np.arange(2, 10, 1), # Clean
    # 'files': np.arange(10, 20, 1),  # Dirty

    # 'expt_name_p1': 'Twentyus_Interleave_P1_Neg7',
    # 'date': '01-08-21',
    # # 'files': np.arange(13, 19, 1), # Clean
    # 'files': np.arange(20, 30, 1),  # Dirty

    # 'expt_name_p1': 'Interleave_P1_Neg100',
    # 'date': '12-30-20',
    # 'files': np.arange(0, 2, 1), # Clean
    # 'files': np.arange(20, 30, 1),  # Dirty

    # 'expt_name_p1': 'Fiveus_Interleave_P1_Neg100',
    # 'date': '01-09-21',
    # 'files': np.arange(0, 5, 1), # Clean
    # 'files': np.arange(40, 58, 1), # Clean
    # 'files': np.arange(15, 23, 1),  # Dirty

    # 'expt_name_p1': 'Fiveus_Interleave_P1_Neg16',
    # 'date': '01-10-21',
    # 'files': np.arange(231, 237, 1), # Clean
    # 'files': np.arange(220, 230, 1),  # Dirty

    'expt_name_p1': 'Fiveus_Interleave_P1_Neg13',
    'date': '01-10-21',
    # 'files': np.arange(235, 240, 1), # Clean
    'files': np.arange(10, 20, 1),  # Dirty

    # 'expt_name_p1': 'Fiveus_Interleave_P1_Neg7',
    # 'date': '01-10-21',
    # 'files': np.arange(210, 220, 1),  # Clean
    # 'files': np.arange(20, 30, 1),  # Dirty
    # 'files': np.arange(178, 180, 1),  # Dirty
}


class OneState(object):
    def __init__(self, t=50e-6):
        self.t = t
        self.time = []
        self.file_axis = []
        self.average_n1 = 10000  # this is the average number from raw matlab data to datachest
        self.expt_info = {}  # this stores the data source and path, etc
        self.t_list = np.arange(0, 10000, 1)
        self.g_list = np.zeros(10000)
        self.e_list = np.zeros(10000)
        self.g_list_norm = np.zeros(10000)
        self.e_list_norm = np.zeros(10000)
        self.g_tpt_list = np.zeros(10000)
        self.e_tpt_list = np.zeros(10000)
        self.g_tpt_pred_list = np.zeros(10000)
        self.e_tpt_pred_list = np.zeros(10000)
        self.g_count_sum = 1
        self.e_count_sum = 1
        self.g_tau_mean = -1
        self.e_tau_mean = -1

    def add_p1_data_from_matlab(self, filenames,
                                data_type='Single_Shot_Occupations'):
        """
        import the data from Matlab and do first level average
        :param filenames: a list of file names
        :param data_type:
        :return: array of single_shot_occupations
        """
        n = self.average_n1

        f_l = filenames[0]
        sample_rate = loadmat_ExptVars(f_l)['Total_Time']
        sample_rate = sample_rate
        self.t = sample_rate    # in unites of us
        self.t_list = sample_rate * self.t_list
        # print('sample_rate=', sample_rate)

        for f in filenames:
            data = loadmat(f)
            one_state_list = np.array(data[data_type])
            for i_os in one_state_list:
                # i_os = generate_hidden_signal(length=10000, charge_burst_time=4000, p_QP=[0.1, 0.1])
                # print('i_os=', i_os[:10])
                i_os_int = [int(i) for i in i_os]
                print('i_os_int=', i_os_int[:10])
                # print('Here', type(i_os_int))
                i_os_HMM = observed_to_recovered_signal_BW(i_os_int)[0]
                print('i_os_HMM=', i_os_HMM[:10])
                self._get_dwell_count(i_os_HMM)


    def _get_dwell_count(self, os):
        curr_state = os[0]
        dwell = 0
        for i, s in enumerate(os):
            if s == curr_state:
                dwell += 1
            else:
                if curr_state == 0:
                    self.g_list[dwell] += 1
                else:
                    self.e_list[dwell] += 1
                curr_state = s
                dwell = 1
            if i == 10000 - 1:
                if curr_state == 0:
                    self.g_list[dwell] += 1
                else:
                    self.e_list[dwell] += 1

        # return g_list, e_list

    def _get_count_sum(self):
        self.g_count_sum = np.sum(self.g_list)
        self.e_count_sum = np.sum(self.e_list)

    def _normalize_hist(self):
        g_list, e_list = self.g_list[:], self.e_list[:]
        g_sum, e_sum = self.g_count_sum, self.e_count_sum
        self.g_list_norm = [float(i) /g_sum for i in g_list]
        self.e_list_norm = [float(i) /e_sum for i in e_list]

    def _get_tpt_list(self):
        g_tpt_list, e_tpt_list = self.g_list_norm[:], self.e_list_norm[:]
        # print ('g_tpt_list=', g_tpt_list[:10])
        # g_tpt_list, e_tpt_list = self.g_list, self.e_list
        t = self.t  # unit time in units of us
        for i, p in enumerate(g_tpt_list):
            g_tpt_list[i] = i * p * t
        for i, p in enumerate(e_tpt_list):
            e_tpt_list[i] = i * p * t
        # print ('g_tpt_list=', g_tpt_list[:10])
        self.g_tpt_list = g_tpt_list
        self.e_tpt_list = e_tpt_list

    def _get_tau_mean(self):
        self.g_tau_mean = np.sum(self.g_tpt_list)
        self.e_tau_mean = np.sum(self.e_tpt_list)

    def _get_tpt_prediction_list(self):
        g_tau = self.g_tau_mean
        e_tau = self.e_tau_mean
        t = self.t  # unit time in units of us

        g_tpt_pred_list = self.g_tpt_pred_list[:]
        e_tpt_pred_list = self.e_tpt_pred_list[:]

        for i in range(len(g_tpt_pred_list)):
            # tau = i * t
            g_tpt_pred_list[i] = t*(i*t/g_tau)*(np.exp(-i*t/g_tau))   # why another t
            e_tpt_pred_list[i] = t*(i*t/e_tau)*(np.exp(-i*t/e_tau))   # why another t

        self.g_tpt_pred_list = g_tpt_pred_list
        self.e_tpt_pred_list = e_tpt_pred_list

    def get_parameters(self):
        self._get_count_sum()
        self._normalize_hist()
        # print('g_count_list=', self.g_list_norm[:10])
        self._get_tpt_list()
        self._get_tau_mean()
        self._get_tpt_prediction_list()


    def plot_data(self):
        t_list = self.t_list
        g_list_norm, e_list_norm = self.g_list_norm[:], self.e_list_norm[:]
        g_list, e_list = self.g_tpt_list[:], self.e_tpt_list[:]
        g_p_list, e_p_list = self.g_tpt_pred_list[:], self.e_tpt_pred_list[:]

        # print('np.sum(g_list_norm)=', g_list_norm[:10])
        # print('np.sum(e_list_norm)=', e_list_norm[:10])
        # print('np.sum(g_list)=', np.sum(g_list))
        # print('np.sum(e_list)=', np.sum(e_list))
        # print('np.sum(g_p_list)=', np.sum(g_p_list))
        # print('np.sum(e_p_list)=', np.sum(e_p_list))
        l_max = 2000
        # plt.semilogy(t_list[1:l_max], g_list_norm[1:l_max], 'o', label='0 State [tau={:.1f} us]'.format(self.g_tau_mean))
        # plt.semilogy(t_list[1:l_max], e_list_norm[1:l_max], 'v-', label='1 State [tau={:.1f} us]'.format(self.e_tau_mean))
        plt.semilogx(t_list[1:l_max], g_list[1:l_max], 'o', label='0 State')
        plt.semilogx(t_list[1:l_max], g_p_list[1:l_max], '--', label='0 State Predicted [tau={:.1f} us]'.format(self.g_tau_mean))
        plt.semilogx(t_list[1:l_max], e_list[1:l_max], 'v', label='1 State')
        plt.semilogx(t_list[1:l_max], e_p_list[1:l_max], '-', label='1 State Predicted [tau={:.1f} us]'.format(self.e_tau_mean))
        plt.xlabel('Dwell Time (us)')
        plt.ylabel('Scaled Count/Probability')
        # plt.title('Normalized Distribution Neg13dBm Jan07 Dirty')
        # plt.title('Normalized Distribution Neg13dBm Dec31 Dirty')
        # plt.title('Normalized Distribution Neg7dBm Jan01 Dirty')
        # plt.title('Normalized Distribution Neg13dBm Jan08 Dirty (20us)')
        # plt.title('Normalized Distribution Neg7dBm Jan08 Dirty (20us)')
        # plt.title('Normalized Distribution Neg100dBm Dec30 Clean (50us)')
        # plt.title('Normalized Distribution Neg100dBm Jan08 Clean (20us)')
        # plt.title('Normalized Distribution Neg100dBm Jan09 Clean (5us)')
        # plt.title('Normalized Distribution Neg100dBm Jan09 Dirty (5us)')
        # plt.title('Normalized Distribution Neg16dBm Jan10 Dirty (5us)')
        plt.title('Normalized Distribution Neg13dBm Jan10 Dirty (5us)')
        # plt.title('Normalized Distribution Neg7dBm Jan10 Clean (5us)')
        plt.legend(loc=2)
        plt.grid()
        plt.show()


def run_matlab_data_to_datachest(ExptInfo):
    os = OneState()
    for key in ExptInfo:
        os.expt_info[key] = ExptInfo[key]

    path = ExptInfo['path']
    date = ExptInfo['date']
    files = ExptInfo['files']
    expt_name_p1 = ExptInfo['expt_name_p1']
    filenames_p1 = [
        path.format(date, expt_name_p1, expt_name_p1) + '_{:03d}.mat'.format(i)
        for i in files]
    os.add_p1_data_from_matlab(filenames_p1)
    os.get_parameters()
    os.plot_data()

run_matlab_data_to_datachest(ExptInfo)

