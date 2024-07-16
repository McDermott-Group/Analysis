from dataChest import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
import os
from itertools import chain
from datetime import datetime
from dateStamp import *
# matplotlib.use("TkAgg")
mpl.rcParams['mathtext.default'] = 'regular'  # this makes all TeX symbols non-italicized

DATACHEST_ROOT = 'Z:\\mcdermott-group\\data\\'
class FitFunctions:
    '''
    for fitting data
    '''
    def exp(t, a, b, c):
        return a * np.exp(-b*t) + c

class T1(object):

    def __init__(self, path, file, dependent_variable='Single Shot Occupation'):

        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        independents = variables[0][0]
        dependents = variables[1]

        data = d.getData(variablesList=[independents[0], dependent_variable]) # want to grab idle time and both dependent variables

        self.path = path
        self.file = file
        self.qb_id = d.getParameter("Qubit ID")
        self.idle_gate_times = data[:,0]/1000  # in microseconds
        self.dependent_variable_name = dependent_variable
        self.dependent_variable_values = data[:,1]

    def plot(self, fit=True, save=False, save_name=None):

        fig = plt.figure(figsize=(8, 6))
        plt.plot(self.idle_gate_times, self.dependent_variable_values, marker='.', linestyle='None')
        plt.xlabel('Idle Gate Time ($\mu$s)', fontsize=25)
        plt.ylabel(self.dependent_variable_name, fontsize=25)
        plt.ion()

        if fit:
            fit_params = self.fit()
            tau0 = 1/fit_params[1] # T1 in us
            plt.plot(self.idle_gate_times, FitFunctions.exp(self.idle_gate_times, *fit_params),
                     color='red',
                     linewidth=2)
            title = 'Qb {0} T1 = {1:.2f} $\mu$s'.format(self.qb_id, tau0)
        else:
            title = 'Qb {} T1'.format(self.qb_id)

        if save:
            plt.title(title, fontsize=25)
            if save_name is None:
                save_name = self.file.replace('.','_')
            plt.ioff()
            if not os.path.exists(DATACHEST_ROOT + self.path + '\\Figures\\'):
                os.makedirs(DATACHEST_ROOT + self.path + '\\Figures\\')
                print(DATACHEST_ROOT + self.path + '\\Figures\\')
            plt.savefig(DATACHEST_ROOT + self.path + '\\Figures\\' + save_name)
            plt.close(fig)
        else:
            plt.title(title, fontsize=25)
            plt.ion()
            plt.show(block=True)


    def fit(self):
        popt, pcov = curve_fit(FitFunctions.exp, self.idle_gate_times, self.dependent_variable_values, p0=[1.0, 0.01, 0.1])
        return popt

class Swap(object):
    '''
    For plotting, fitting, and saving SWAP experiments. Can also add new dataset to dataChest
    of (frequency, T1) data after fitting.
    '''
    def __init__(self, path, file, dependent_variable='SSO'):

        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        independents = variables[0]
        dependents = variables[1]

        data = d.getData(variablesList=[independents[0][0], independents[1][0], dependent_variable])

        self.path = path
        self.file = file
        self.qb_id = d.getParameter("Qubit ID")
        self.dependent_variable_name = dependent_variable
        self.dependent_variable_values = data[:,2]
        self.frequencies = np.unique(data[:,0]) # in GHz
        self.idle_gate_times = np.unique(data[:,1] * 1e-3) # in us

        # cutting off highest 10 MHz because LO bleedthrough was messing up fits; should make frequency exceptions a method
        # self.frequencies = self.frequencies[:-9]
        # self.dependent_variable_values = self.dependent_variable_values[:-9*len(self.idle_gate_times)]

    def plot(self, save=False, save_name=None):
        '''
        plots 2d color plot with frequencies on the x-axis, idle times on the y-axis,
        and SSO on the z-axis.
        '''

        # vertical cut of SWAP scan
        print(self.dependent_variable_values)
        line_cuts = np.array(np.split(self.dependent_variable_values, len(self.frequencies)))

        fig = plt.figure(figsize=(12, 8))

        cmap = plt.get_cmap('rainbow')
        plt.set_cmap(cmap)
        cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap), orientation='vertical')
        cbar.ax.tick_params(labelsize=20)

        plt.pcolormesh(self.frequencies, self.idle_gate_times, line_cuts.transpose())
        plt.title("Q{} SWAP Spectroscopy".format(self.qb_id), fontsize=25)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel("Qubit Frequency (GHz)", fontsize=25)
        plt.ylabel("Idle Time ($\mu$s)", fontsize=25)

        if save:
            if save_name is None:
                save_name =self.file.replace('.','_')
            plt.ioff()
            if not os.path.exists(DATACHEST_ROOT + self.path + '\\Figures\\'):
                os.makedirs(DATACHEST_ROOT + self.path + '\\Figures\\')
                print(DATACHEST_ROOT + self.path + '\\Figures\\')
            plt.savefig(DATACHEST_ROOT + self.path + '\\Figures\\' + save_name)
            plt.close(fig)
        else:
            plt.ion()
            plt.show(block=True)


    def fit(self, add_fitted_dataset=False):
        '''
        fits each vertical line cut of the SWAP scan and extracts T1. returns list of tuples
        (frequency, T1). if save, creates new dataset in folder "fitted_data" inside data folder
        '''

        line_cuts = np.array(np.split(self.dependent_variable_values, len(self.frequencies)))

        fit_data = []
        for i, cut in enumerate(line_cuts):
            popt, pcov = curve_fit(FitFunctions.exp, self.idle_gate_times, cut, p0=[0.9, 0.01, 0.1])
            tau0 = 1/popt[1]  # T1 in us
            # if tau0 > 325:
            #     plt.plot(self.idle_gate_times, cut)
            #     plt.plot(self.idle_gate_times, FitFunctions.exp(self.idle_gate_times, *popt))
            #     plt.title("T1 = {}".format(tau0))
            #     plt.show(block=True)
            fit_data.append([self.frequencies[i], tau0])

        mean_T1 = np.average([item[1] for item in fit_data])
        stdev_T1 = np.std(list(zip(*fit_data))[1])
        max_T1 = np.max([item[1] for item in fit_data])
        print('Average T1 = {:.2f} $\mu$s'.format(mean_T1))
        print('Standard deviation of T1 = {:.2f} $\mu$s'.format(stdev_T1))
        print('Maximum T1 = {:.2f} $\mu$s'.format(max_T1))

        if add_fitted_dataset:
            self.add_fitted_dataset(fit_data, mean_T1, stdev_T1)

        return fit_data

    def add_fitted_dataset(self, fit_data, mean_T1, stdev_T1, save_name=None):

        if save_name is None:
            save_name = self.file + "_fitted"
        if not os.path.exists(DATACHEST_ROOT + self.path + "\\fitted_data\\"):
            os.makedirs(DATACHEST_ROOT + self.path + "\\fitted_data\\")

        d = dataChest(self.path + "\\fitted_data")
        d.createDataset(save_name,
                        [("Frequency", [1], "float64", "Hz")],
                        [("T1", [1], "float64", "us")])

        d.addData(fit_data)
        d.addParameter("Qubit ID", self.qb_id)
        d.addParameter("Average T1", mean_T1, "us")
        d.addParameter("Std T1", stdev_T1, "us")

class SwapVsTime(object):
    '''
    for plotting SWAP vs Time data. the folder passed in should contain only datasets of T1 vs frequency, generated using the Swap class.
    '''
    def __init__(self, path):

        d = dataChest(path)

        self.path = path
        self.T1s = []
        for file in os.listdir(DATACHEST_ROOT + path):

            d.openDataset(file)
            variables = d.getVariables()
            independents = variables[0][0]
            dependents = variables[1][0]

            data = d.getData(variablesList=[independents[0],dependents[0]])

            T1s = data[:,1]
            self.frequencies = data[:,0]
            self.T1s.append(T1s)

    def plot(self, save=False, time_step=None, save_name=None):
        '''
        plots the time
        :param save:
        :param time_step:
        :param save_name:
        :return:
        '''
        if time_step is None:  # we will extract the timestep using the difference in timestamps from two consecutive files
            timestamp_01 = os.listdir(DATACHEST_ROOT + self.path + "\..")[0][:10]  # grabs just the timestamp from file name
            timestamp_02 = os.listdir(DATACHEST_ROOT + self.path + "\..")[1][:10]

            utcTime_01 = dateStamp().invertDateStamp(timestamp_01)  # converts timestamp to utc time
            utcTime_02 = dateStamp().invertDateStamp(timestamp_02)

            datetime_01 = datetime.strptime(utcTime_01, '%Y-%m-%dT%H:%M:%S.%f')
            datetime_02 = datetime.strptime(utcTime_02, '%Y-%m-%dT%H:%M:%S.%f')

            time_step = (datetime_02 - datetime_01).total_seconds()

            # times = np.arange(0, time_step * len(self.T1s)/ 3600, time_step/3600)
            times = np.linspace(0, time_step * len(self.T1s) / 3600, len(self.T1s))
        else:
            times = np.arange(0, time_step * len(self.T1s)/ 3600, time_step/3600)

        df = self.frequencies[1] - self.frequencies[0]

        fig = plt.figure(figsize=(8, 6))
        plt.pcolormesh(self.frequencies, times, self.T1s)

        cmap = plt.get_cmap('Spectral')
        plt.set_cmap(cmap)
        cbar = plt.colorbar(orientation='vertical')
        cbar.ax.tick_params(labelsize=20)
        cbar.set_label('T1 ($\mu$s)', fontsize=20, rotation=0, labelpad=50)

        plt.title("Q1 SWAP Spectroscopy vs Time", fontsize=25)
        plt.xticks(np.arange(self.frequencies[0], self.frequencies[-1] + df , 100*df), fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel("Qubit Frequency (GHz)", fontsize=25)
        plt.ylabel("Time (hrs)", fontsize=25)

        if save:
            if save_name is None:
                save_name = self.file.replace('.','_')
            plt.ioff()
            if not os.path.exists(DATACHEST_ROOT + self.path + '\\Figures\\'):
                os.makedirs(DATACHEST_ROOT + self.path + '\\Figures\\')
                print(DATACHEST_ROOT + self.path + '\\Figures\\')
            plt.savefig(DATACHEST_ROOT + self.path + '\\Figures\\' + save_name)
            plt.close(fig)
        else:
            plt.ion()
            plt.show(block=True)

    def histogram(self, num_bins):
        T1s = [sublist.tolist() for sublist in self.T1s]  # each sublist is n T1 values, one for each frequency bin
        T1s_split = np.array_split(T1s, num_bins)  # this splits the 2D list into chunks
        # use itertools' chain method here as built=in sum is O(n^2)
        T1s_split = [list(chain.from_iterable(chunk)) for chunk in T1s_split]  # this combines all of the sublists within each chunk into one large list

        return T1s_split, self.T1s


class QubitSpecVsTime(object):
    '''
    class for analyzing repeated qubit spectroscopy, run interleaved with repeated swap. for each file, this code will find the minimum in amplitude and make
    a new dataset containing chronological data of min_amplitude vs qubit frequency
    '''

    def __init__(self, path):

        self.path = path
        self.min_frequencies = []

        d = dataChest(path)

        for file in sorted(os.listdir(DATACHEST_ROOT + self.path)):
            d.openDataset(file)
            variables = d.getVariables()
            indep_var = variables[0][0]
            data = d.getData(variablesList=[indep_var[0], 'Amplitude'])
            frequencies = data[:, 0]
            amplitude = data[:, 1]
            self.min_frequencies.append(frequencies[np.argmin(amplitude)])
            print(frequencies[np.argmin(amplitude)])


        save_name = "extracted_minimums"
        if not os.path.exists(DATACHEST_ROOT + self.path + "\\fitted_data\\"):
            os.makedirs(DATACHEST_ROOT + self.path + "\\fitted_data\\")

        newData = dataChest(self.path + "\\fitted_data")
        newData.createDataset(save_name,
                        [("Frequencies", [1], "float64", "Hz")],
                              [("Amplitudes", [1], "float64", "V")])

        data = [[freq, 1.0] for freq in self.min_frequencies]
        print(data)
        newData.addData(data)

    def plot(self, save=False):
        pass
