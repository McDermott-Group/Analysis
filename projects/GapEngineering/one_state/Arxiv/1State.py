from noiselib import loadmat
import numpy as np
import matplotlib.pyplot as plt

class one_state(object):
    ###TO DO:
    #1. Import data
    #2. dummy average => moving average
    def __init__(self, t = 10e-6):
        self.t = t
        self.one_state_array = []
        self.one_state_avg = []
        self.time = []
    def add_data_from_matlab(self, filenames, data_type='Single_Shot_Occupations'):
    # def add_data_from_matlab(self, filenames, data_type='Amplitudes'):
    # def add_data_from_matlab(self, filenames, data_type='Is'):
        """

        :param filenames: a list of file names
        :param data_type:
        :return: array of single_shot_occupations
        """
        for file in filenames:
            data = loadmat(file)
            one_state_list = np.array(data[data_type])
            for os in one_state_list:
                self.one_state_array.append(os)
            # print type(self.one_state_array)
                

    def get_one_state_avg(self, n=1000):
        one_state_array = np.asarray(self.one_state_array)
        self.one_state_avg = np.mean(one_state_array.reshape(-1, n), axis=1)
        t_avg = n * self.t
        self.time = (np.arange(0, len(self.one_state_avg), 1))*t_avg
        # print len(self.one_state_avg)


    def plot_PSD(self, fit = False):
        self.get_PSD()
        fig = plt.figure()
        if fit:
            self.fit_PSD()
            plt.loglog(self.f_data, self.fit_PSD_target_function(self.f_data, self.params[0], self.params[1]), label='Fitted')
        plt.loglog(self.f_data, self.Spp_avg, label='Test')
        axes = plt.gca()
        axes.set_ylim([1e-6, 1e-2])
        # axes.set_ylim([1e-6, 1e-1])
        plt.xlabel('frequency [Hz]')
        plt.ylabel('PSD [1/Hz]')
        plt.show()


date = '09-28-20'
path = ('Z:/mcdermott-group/data/GapEngineer/Nb_GND_Dev06_Trap/Leiden_2020Jul/Debug/LIU/Q4_withQ5Poison/{}/{}/MATLABData/{}')

# exp_name = ('One_State_Poison_Off')
# files = np.arange(49, 149, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_Off = one_state()
# os_test_Off.add_data_from_matlab(filenames)
# os_test_Off.get_one_state_avg()

# exp_name = ('One_State_Poison_Neg4dBm')
# date = '09-27-20'
# files = np.arange(0, 50, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_Neg4dBm = one_state()
# os_test_Neg4dBm.add_data_from_matlab(filenames)
# os_test_Neg4dBm.get_one_state_avg()

# exp_name = ('One_State_Poison_Neg2dBm')
# files = np.arange(0, 500, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_Neg2dBm = one_state()
# os_test_Neg2dBm.add_data_from_matlab(filenames)
# os_test_Neg2dBm.get_one_state_avg()

# exp_name = ('One_State_Poison_0dBm')
# files = np.arange(50, 550, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_0dBm = one_state()
# os_test_0dBm.add_data_from_matlab(filenames)
# os_test_0dBm.get_one_state_avg()

# exp_name = ('One_State_Poison_1dBm')
# files = np.arange(50, 550, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_1dBm = one_state()
# os_test_1dBm.add_data_from_matlab(filenames)
# os_test_1dBm.get_one_state_avg()

# exp_name = ('One_State_Poison_2dBm')
# files = np.arange(50, 550, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_2dBm = one_state()
# os_test_2dBm.add_data_from_matlab(filenames)
# os_test_2dBm.get_one_state_avg()

# exp_name = ('One_State_Poison_3dBm')
# files = np.arange(50, 550, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_0dBm = one_state()
# os_test_0dBm.add_data_from_matlab(filenames)
# os_test_0dBm.get_one_state_avg()

# exp_name = ('One_State_Poison_4dBm')
# files = np.arange(50, 150, 1)
# filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
# os_test_4dBm = one_state()
# os_test_4dBm.add_data_from_matlab(filenames)
# os_test_4dBm.get_one_state_avg()

exp_name = ('One_State_Poison_Neg6dBm')
date = '09-28-20'
files = np.arange(0, 100, 1)
filenames = [path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in files]
os_test_Neg6dBm = one_state()
os_test_Neg6dBm.add_data_from_matlab(filenames)
os_test_Neg6dBm.get_one_state_avg()

# fig = plt.figure(1)
# plt.xlabel('time (ms)')
# plt.ylabel('avg_one_state')
# plt.plot((os_test_Off.time)*10**3, os_test_Off.one_state_avg)
# plt.title('Poison Off')

# fig = plt.figure(2)
# plt.xlabel('time (ms)')
# plt.ylabel('avg_one_state')
# plt.ylim(-0.1, 1.1)
# plt.plot((os_test_Neg4dBm.time)*10**3, os_test_Neg4dBm.one_state_avg)
# plt.title('Poison -4 dBm')

# fig = plt.figure(3)
# plt.xlabel('time (ms)')
# plt.ylabel('avg_one_state')
# plt.plot((os_test_Neg2dBm.time)*10**3, os_test_Neg2dBm.one_state_avg)
# plt.title('Poison -2 dBm')

# fig = plt.figure(4)
# plt.xlabel('time (ms)')
# plt.ylabel('avg_one_state')
# plt.plot((os_test_0dBm.time)*10**3, os_test_0dBm.one_state_avg)
# plt.title('Poison 0 dBm')

# fig = plt.figure(5)
# plt.xlabel('time (ms)')
# plt.ylabel('avg_one_state')
# plt.plot((os_test_1dBm.time)*10**3, os_test_1dBm.one_state_avg)
# plt.title('Poison 1 dBm')

# fig = plt.figure(6)
# plt.xlabel('time (ms)')
# plt.ylabel('avg_one_state')
# plt.plot((os_test_2dBm.time)*10**3, os_test_2dBm.one_state_avg)
# plt.title('Poison 2 dBm')

# fig = plt.figure(7)
# plt.xlabel('time (ms)')
# plt.ylabel('avg_one_state')
# plt.plot((os_test_4dBm.time)*10**3, os_test_4dBm.one_state_avg)
# plt.title('Poison 4 dBm')

fig = plt.figure(1)
plt.xlabel('time (ms)')
plt.ylabel('avg_one_state')
# plt.ylim(0.8*10**(-4), 2.0*10**(-4))
plt.ylim(-0.1, 1.1)
plt.plot((os_test_Neg6dBm.time)*10**3, os_test_Neg6dBm.one_state_avg)
plt.title('Poison -6 dBm')

plt.show()



