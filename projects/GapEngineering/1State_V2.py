from noiselib import loadmat
import numpy as np
import matplotlib.pyplot as plt

class one_state(object):
    ###TO DO:
    #1. Import data
    #2. dummy average => moving average
    def __init__(self, t = 10e-6):
        self.t = t
        # self.one_state_array = []
        self.one_state_avg = np.array([])
        self.time = []
        # self.name = name
    def add_data_from_matlab(self, filenames, data_type='Single_Shot_Occupations'):
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
                # self.one_state_array.append(os)
                self._get_one_state_avg(os=os, n=n)
        t_avg = n * self.t
        self.time = (np.arange(0, len(self.one_state_avg), 1))*t_avg
                
    def _get_one_state_avg(self, os, n):
        one_state_array = np.asarray(os)
        avg = np.mean(one_state_array.reshape(-1, n), axis=1)
        self.one_state_avg = np.append(self.one_state_avg, avg)
        
    def plot_one_state_avg_vs_time(self):
        plt.xlabel('time (ms)')
        plt.ylabel('avg_one_state')
        plt.ylim(-0.1, 1.1)
        plt.plot((self.time)*10**3, self.one_state_avg)
        # plt.title(self.name)

date = '09-28-20'
path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/Debug/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')

exp_name = ('One_State_Poison_Neg7dBm')
date = '09-29-20'
files = np.arange(0, 500, 1)
filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
os_test_Neg7dBm = one_state()
os_test_Neg7dBm.add_data_from_matlab(filenames)

# exp_name = ('One_State_Poison_Neg80dBm8GHz')
# date = '09-29-20'
# files = np.arange(0, 36, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_Neg80dBm = one_state(name=str(Poison))
# os_test_Neg80dBm.add_data_from_matlab(filenames)

# fig = plt.figure(1)
# plt.xlabel('time (ms)')
# plt.ylabel('avg_one_state')
# plt.ylim(-0.1, 1.1)
# plt.plot((os_test_Neg6dBm.time)*10**3, os_test_Neg6dBm.one_state_avg)
# plt.title('Poison -6 dBm')

fig = plt.figure(2)
os_test_Neg7dBm.plot_one_state_avg_vs_time()
plt.title('Poison -7 dBm')
plt.show()



