from datetime import datetime
from scipy import signal
from dataChest import *
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import platform
import matplotlib
matplotlib.use("TkAgg")

DATACHEST_ROOT = '/Volumes/smb/mcdermott-group/data/' if 'macOS' in platform.platform() else 'Z:\\mcdermott-group\\data\\'

class FitFunc:
    """Class contains only static fit functions for use with curve_fit; never any need to instantiate.

    """
    def __init__(self):
        pass

    @staticmethod
    def exp(t, a, b, c):
        return a * np.exp(b * t) + c

class SaveHelper:
    """Helper class contains only static methods; never any need to instantiate.

    """
    def __init__(self):
        pass

class DataChestHelper:
    """Helper class contains only static methods; never any need to instantiate.

    """
    def __init__(self):
        pass


class T1(object):
    """For fitting, plotting, and analyzing T1 experiments.

    """
    def __init__(self, path, file, dependent_variable='Single Shot Occupation'):
        """Instantiate with given file for analysis.

        :param path: path containing data files.
        :param file: hdf5 file to open in data chest.
        :param dependent_variable: dependent variable for fitting, plotting, and analyzing.
        """

        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        independents = variables[0]
        if 'Idle Gate Time' not in independents:
            raise Exception('For T1 analysis, independent variable must be idle gate time.')
        dependents = variables[1]
        data = d.getData(variablesList=[independents[0], dependent_variable])
        self._path = path
        self._file = file
        self._qubit_id = d.getParameter('Qubit ID')
        self._idle_gate_times = data[:,0]/1000
        self._dependent_variable_name = dependent_variable
        self._dependent_variable_values = data[:,1]

    def plot(self, fit=True, save=False, save_path=None, save_name=None):
        """Plot the T1 data.

        :param fit: fit the data and plot fit with data?
        :param save: save the plot (and the extracted T1 and P1, if fit)?
        :param save_path: save path for figures.  If None, defaults to '../Figures/'.
        :param save_name: save name for figure.  If None, defaults to hdf5 file name.
        :return:
        """
        #need to edit
        fig = plt.figure()
        plt.title('Q{:1d}: T1'.format(self._qubit_id))
        plt.plot(self._idle_gate_times, self._dependent_variable_values, marker='.', linestyle='None')
        plt.xlabel('Idle Gate Time (us)')
        plt.ylabel(self._dependent_variable_name)
        if fit:
            fit_parameters = self.fit(save)
            t1=(-1/fit_parameters[1])
            if t1 > 4:
                plt.title('Q{:1d}: T1 = {:2d}us'.format(self._qubit_id,int(t1)))
            else:
                plt.title('Q{:1d}: T1 = {:.2f}us'.format(self._qubit_id, t1))
            plt.plot(self._idle_gate_times, FitFunc.exp(self._idle_gate_times,*fit_parameters))
        if save:
            if save_path is None:
                save_path = os.path.join(*self._path[:-1])
            if save_name is None:
                save_name =self._file.replace('.','_')
            plt.ioff()
            if not os.path.exists(DATACHEST_ROOT + save_path + '\\Figures\\'):
                os.makedirs(DATACHEST_ROOT + save_path + '\\Figures\\')
            plt.savefig(DATACHEST_ROOT + save_path + '\\Figures\\' + save_name)
            plt.close(fig)
        else:
            plt.ion()
            plt.show()

    def fit(self, save):
        """Fit the data

        :param save: update data chest with T1 and P1?
        :return: [T1, P1]
        """
        try:
            popt, pcov = curve_fit(FitFunc.exp, self._idle_gate_times, self._dependent_variable_values,p0=[1,(-1/15),0])
        except:
            #The fit above tends to fail in the BB experiments when T1 is small and P1 is big. This is a quick fix:
            popt, pcov = curve_fit(FitFunc.exp, self._idle_gate_times,
                                   self._dependent_variable_values, p0=[0.5, (-1 / 1), 0.5])
        if save:
            self._update_dataChest(int(-100 / popt[1]) / 100, round(popt[2], 4))
        return popt

    def _update_dataChest(self, t1, p1):
        """Updates data chest with the fit values for T1 and P1.

        :param t1: T1 to save to data chest
        :param p1: P1 to save to data chest
        :return:
        """
        d = dataChest(self._path)
        d.openDataset(self._file, modify = True)
        d.addParameter('Fit T1', t1, 'us', overwrite = True)
        d.addParameter('Fit P1', p1, overwrite=True)
        d.addParameter('Fit Date Stamp', datetime.now().strftime("%Y-%m-%d"), overwrite = True)

class T2(object):
    """For fitting, plotting, and analyzing T2* experiments.

    """
    def __init__(self, path, file, dependent_variable='Single Shot Occupation'):
        """Instantiate with given file for analysis.

        :param path: path containing data files.
        :param file: hdf5 file to open in data chest.
        :param dependent_variable: dependent variable for fitting, plotting, and analyzing.
        """
        pass

class P1(object):
    """For fitting, plotting, and analyzing all types of P1 experiments

    """
    def __init__(self, path, file, dependent_variable='Single Shot Occupation'):
        """Instantiate with given file for analysis.

        :param path: path containing data files.
        :param file: hdf5 file to open in data chest.
        :param dependent_variable: dependent variable for fitting, plotting, and analyzing.
        """
        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        independents = variables[0]
        data = d.getData(variablesList=[independents[0], dependent_variable])
        self._file = file
        self._path = path
        self._biases = data[:,0]
        self._dependent_variable_values = data[:,1]
        self._dependent_variable_name = dependent_variable

    def plot_versus_radiator_frequency(self, bias_offset=0, bias_conversion=2*10*484, save=False, save_path=None, save_name=None):
        """Plots P1 data vs JJ Radiator Frequency

        :param bias_offset: Voltage offset to substract from OPX output voltage
        :param bias_conversion: Factor to take into account line attenutation to get on chip voltage, then 484GHz/mV, to convert V to f
        :param save: Save the plot?
        :param save_path: Save path. If None, defaults to '../Figures'
        :param save_name: Save name. If None, defaults to hdf5 file name.
        :return: None.
        """

        fig = plt.figure()
        plt.title('P1 vs Radiator Bias')
        plt.semilogy([(b-bias_offset)*bias_conversion for b in self._biases], self._dependent_variable_values, marker='.', linestyle='None')
        plt.xlabel('Radiator Bias (GHz)')
        plt.ylabel(self._dependent_variable_name)
        if save:
            if save_path is None:
                save_path = os.path.join(*self._path[:-1])
            if save_name is None:
                save_name =self._file.replace('.','_')
            plt.ioff()
            if not os.path.exists(DATACHEST_ROOT + save_path + '\\Figures\\'):
                os.makedirs(DATACHEST_ROOT + save_path + '\\Figures\\')
            plt.savefig(DATACHEST_ROOT + save_path + '\\Figures\\' + save_name)
            plt.close(fig)
        else:
            plt.ion()
            plt.show()

    def plot_versus_radiator_flux(self, bias_offset=0, bias_conversion=1, save=False, save_path=None, save_name=None):
        """Plots P1 data at fixed JJ Radiator Frequency versus JJ Radiator Flux Bias

        :param bias_offset: Bias offset
        :param bias_conversion: Conversion factor between OPX voltage and current on-chip, based on line attenuation.
        :param save: Save the plot?
        :param save_path: Save path. If None, defaults to '../Figures'
        :param save_name: Save name. If None, defaults to hdf5 file name.
        :return: None.
        """
        #need to edit
        fig = plt.figure()
        plt.title('P1 vs Radiator Flux')
        plt.semilogy([b*bias_conversion+bias_offset for b in self._biases], self._dependent_variable_values, marker='.', linestyle='None')
        plt.xlabel('Radiator Flux Bias (mA)')
        plt.ylabel(self._dependent_variable_name)
        if save:
            if save_path is None:
                save_path = os.path.dirname(os.path.abspath(__file__))
            if save_name is None:
                save_name = self._file.replace('.','_')
            plt.ioff()
            plt.savefig(save_path+ '\\Figures\\' + save_name)
            plt.close(fig)
        else:
            plt.ion()
            plt.show()

class Parity(object):
    """For fitting, plotting, and analyzing parity experiments

    """
    def __init__(self, path, file, rad_id=None, parity_guess=100):
        """Instantiate with given file for analysis.

        :param path: path containing data files.
        :param file: hdf5 file to open in data chest.
        :param rad_id: the radiator ID, for old files where the value was not saved to data chest.
        :param parity_guess: a guess for the parity switching rate, in Hz.  Used to assist in fitting.
        """
        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        data = d.getData()
        data = data.transpose()
        self._iterations = int(np.amax(data[0]) + 1)  # Number of Iterations
        self._reps = int(np.amax(data[1]) + 1)  # Number of Reps
        self._data = []
        self._fit_psd = None
        self._parity_guess = parity_guess
        for i in range(self._iterations):
            b = []
            result = list(np.where(data[0].astype(int) == i))
            for j in result:
                b.append((data[2][j]))
            self._data.append(b[0])
        self._path = path
        self._file = file
        self._qubit_id = d.getParameter('Qubit ID')
        self._rad_id = rad_id if rad_id is not None else d.getParameter('Radiator ID')
        self._bias = d.getParameter('JB{:1d} Voltage Bias'.format(self._rad_id))[0]
        self._repetition_rate = 1e9/d.getParameter('Time Per Iteration')[0]
        self._psd()

    def _psd(self):
        """Estimate the PSD from the data using signal.periodogram

        :return: None
        """
        state = self._data
        n_max = self._iterations

        psd = [None] * n_max
        for i in range(n_max):
            freqs, psd[i] = signal.periodogram(state[i], fs=self._repetition_rate, return_onesided=True)
        ave_psd = np.average(psd, axis=0)
        self._psd = ave_psd[1:]
        self._f = freqs[1:]

    def psd(self):
        """f, PSD computed from the time record

        :return: [f, PSD]
        """
        return self._f, self._psd

    def fit_psd(self):
        """f, fit PSD

        :return: f, fit PSD
        """
        if self._fit_psd is None:
            self.fit()
        return self._f, self._fit_psd

    def fit(self, save=False):
        """fit to the PSD data

        :param save: Update data chest with parity rate and fidelity from fit?
        :return: parity switching rate, fidelity
        """
        fs=self._repetition_rate
        def fit_PSD_target_function(f, T_parity, F_map):
            return 2 * (4 * F_map ** 2 / T_parity) / ((2 / T_parity) ** 2 + (2 * np.pi * f) ** 2) + 2 * (
                        1 - F_map ** 2) / fs

        initial_guess_fidelity = np.arange(0, 1, 0.05)
        lower_parity_guess = int(np.floor(np.log10(self._parity_guess/100)))
        initial_guess_parity = np.logspace(lower_parity_guess, lower_parity_guess+6, 7)
        covariance = float('inf')
        best_r_squared = float('inf')
        freqs,ave_psd = self.psd()
        for igf in initial_guess_fidelity:
            for igp in initial_guess_parity:
                params_curr, params_covariance_curr = curve_fit(
                    fit_PSD_target_function, freqs[5:], ave_psd[5:], p0=[igp, igf], bounds=([10**lower_parity_guess, 0], [10**(lower_parity_guess+6), 1]),
                    method='trf')

                residuals = ave_psd - fit_PSD_target_function(freqs, *params_curr)
                ss_res = np.sum(residuals ** 2)
                ss_tot = np.sum((ave_psd - np.mean(ave_psd)) ** 2)
                r_squared = 1-(ss_res / ss_tot) #r^2 is really 1-this but this is what I want to use in the next line....
                if np.abs(1-r_squared) < np.abs(1-best_r_squared):
                    best_r_squared = r_squared
                    params = params_curr
        T_parity = params[0]
        F_map = params[1]
        fit_psd = fit_PSD_target_function(freqs, T_parity, F_map)

        self._fit_psd = fit_psd
        self._parity_rate = 1/T_parity
        self._fit_fidelity = F_map
        if save:
            self.update_dataChest()

        return 1 / T_parity, F_map

    def plot(self, fit=True, save=False, save_path=None, save_name=None, title=None):
        """Plot the PSD data and, optionally, the fit to the data

        :param fit: Plot the fit? (And do the fitting if not yet done.)
        :param save: Save the plot and update data chest with fit? If False, plots in interactive mode.
        :param save_path: Path to save the plot. If None, defaults to '../Figures/'
        :param save_name: Save name for the plot. If None, defaults to the hdf5 file name
        :param title: human-readable plot title.  Currently not implemented, but can add in future if useful.
        :return: None
        """
        fig=plt.figure()
        f,psd = self.psd()
        plt.loglog(f,psd)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('PSD (1/Hz)')
        if fit:
            f, fit_psd = self.fit_psd()
            plt.loglog(f,fit_psd)
            title = 'Q{:1d} | J{:1d} Bias = {:2.0f} mV | Parity Rate = {:3.0f} Hz | Fidelity = {:2.0f}%'.format(self._qubit_id,self._rad_id,1000*self._bias,self._parity_rate,100*self._fit_fidelity) if title is None else title
        else:
            title = 'Q{:1d} | J{:1d} Bias = {:2.0f} mV'.format(self._qubit_id,self._rad_id,1000*self._bias) if title is None else title
        plt.title(title)
        if save:
            self.update_dataChest()
            if save_path is None:
                save_path = os.path.join(*self._path[:-1])
            if save_name is None:
                save_name = self._file.replace('.', '_')
            plt.ioff()
            if not os.path.exists(DATACHEST_ROOT+save_path+'/Figures/'):
                os.makedirs(DATACHEST_ROOT+save_path+'/Figures/')
            print((DATACHEST_ROOT+save_path+ '/Figures/' + save_name))
            plt.savefig(DATACHEST_ROOT+save_path+ '/Figures/' + save_name)
            plt.close(fig)
        else:
            plt.ion()
            plt.show()

    def update_dataChest(self):
        """Update DataChest with fit information.

        :return: None
        """
        d = dataChest(self._path)
        d.openDataset(self._file, modify = True)
        d.addParameter('Fit Parity Rate', self._parity_rate, 'Hz', overwrite = True)
        d.addParameter('Fit Fidelity', self._fit_fidelity, overwrite=True)
        d.addParameter('Fit Date Stamp', datetime.now().strftime("%Y-%m-%d"), overwrite = True)

class IQBlobs(object):
    """For plotting IQ blobs

    """
    _Igs = []
    _Qgs = []
    _Ies = []
    _Qes = []
    _qubit_id = 0

    def __init__(self, path, files):
        """Instantiate with files containing IQ data to plot.

        :param path: path containing data files.
        :param files: list hdf5 files with Ig, Ie, Qg, and Qe data to plot.
        """
        d = dataChest(path)
        self._path = path
        for file in files:
            d.openDataset(file)
            self._qubit_id = d.getParameter('Qubit ID')
            variables = d.getVariables()
            data = d.getData().transpose()
            dependents = [variables[i][0][0].decode() for i in range(2)]
            for i,var in enumerate(dependents):
                if var == 'Ig':
                    self._Igs = data[i]
                elif var == 'Ie':
                    self._Ies= data[i]
                elif var == 'Qg':
                    self._Qgs = data[i]
                elif var == 'Qe':
                    self._Qes = data[i]
                else:
                    raise Exception('Problem with Dependent Variables. Check that this is the right experiment.')

    def plot(self, save=False, save_path=None, save_name=None):
        """Plot the IQ Blobs

        :param save: Save the plot? If False, plots in interactive mode.
        :param save_path: Path for saving the plot. If None, defaults to '../Figures'.
        :param save_name: Save name for plot. If None, defaults to hdf5 file name.
        :return:
        """
        fig = plt.figure()
        plt.title('Q{:1d}: IQ Data'.format(self._qubit_id))
        plt.xlabel('I (a.u.)')
        plt.ylabel('Q (a.u.)')
        plt.plot(self._Igs,self._Qgs,linestyle="None",marker='o',markersize='1')
        plt.plot(self._Ies,self._Qes,linestyle="None",marker='o',markersize='1')
        if save:
            if save_path is None:
                save_path = os.path.join(*self._path[:-1])
            if save_name is None:
                save_name = 'Q{:1d} IQ Data.png'.format(self._qubit_id)
            plt.ioff()
            if not os.path.exists(DATACHEST_ROOT + save_path + '\\Figures\\'):
                os.makedirs(DATACHEST_ROOT + save_path + '\\Figures\\')
            plt.savefig(DATACHEST_ROOT + save_path + '\\Figures\\' + save_name)
            plt.close(fig)
        else:
            plt.ion()
            plt.show()

class TwoDimensionalRamsey(object):
    """For plotting 2D Ramsey Scans

    """
    def __init__(self, path, file, dependent_variable="Single Shot Occupation", sequence_name="X/2-Idle-Y/2"):
        """Instantiate with given file for plotting.

        :param path: path containing data files.
        :param file: hdf5 file to open in data chest.
        :param dependent_variable: dependent variable for fitting, plotting, and analyzing.
        :param sequence_name: human-readable name of Ramsey sequence, for plotting, etc.
        """
        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        data = d.getData()
        independents = variables[0]
        dependents = variables[1]
        self._path = path
        self._file = file
        self._sequence_name = sequence_name
        self._qubit_id = d.getParameter('Qubit ID')
        transposed_data = np.transpose(data)
        self._independent_1 = independents[0]
        self._independent_2 = independents[1]
        self._independent_1_values = np.sort(np.unique(transposed_data[0]))
        self._independent_2_values = np.sort(np.unique(transposed_data[1]))
        len_1 = len(self._independent_1_values)
        len_2 = len(self._independent_2_values)
        self._dependent_variable_values = np.empty([len_1,len_2])
        data = data[np.lexsort((data[:,1],data[:,0]))]

        idx = 0
        for i,dependent in enumerate(dependents):
            if dependent[0] == dependent_variable.encode():
                idx = i + 2
                self._dependent_variable = dependent
                break
        for i in range(len_1):
            for j in range(len_2):
                self._dependent_variable_values[i][j] = data[i*len_2+j][idx]
        self._dependent_variable_values = self._dependent_variable_values.transpose()

    def plot(self, save=False, save_path=None, save_name=None):
        """Plot the 2D Ramsey Scan

        :param save: Save the plot? If False, displays in interactive mode. (Not yet implemented).
        :param save_path: Save path. If None, defaults to '../Figures'. (Not yet implemented).
        :param save_name: Save name. If None, defaults to hdf5 file name. (Not yet implemented).
        :return:
        """
        fig = plt.figure()
        plt.xlabel('{0} ({1})'.format(self._independent_1[0].decode(),self._independent_1[3].decode()))
        plt.ylabel('{0} ({1})'.format(self._independent_2[0].decode(),self._independent_2[3].decode()))
        if self._dependent_variable[3].decode() != '':
            plt.title('{0} \n {1} ({2})'.format(self._sequence_name,self._dependent_variable[0].decode(),self._dependent_variable[3].decode()))
        else:
            plt.title('{0} \n {1}'.format(self._sequence_name,self._dependent_variable[0].decode()))
        c = plt.pcolor(np.log(self._independent_1_values),self._independent_2_values,self._dependent_variable_values)
        fig.colorbar(c)
        plt.show()