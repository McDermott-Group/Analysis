from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib

DATACHEST_ROOT = 'Z:\\mcdermott-group\\data\\'

class T1(object):

    def __init__(self, path, file, dependent_variable='Single Shot Occupation'):
        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        data = d.getData()
        independents = variables[0]
        dependents = variables[1]
        self._path = path
        self._file = file
        self._qubit_id = d.getParameter('Qubit ID')
        self._idle_gate_times = data[:,0]/1000
        self._dependent_variable_values = []
        self._dependent_variable_name = dependent_variable
        #there absolutely has to be a better way of doing this, but:
        for i in range(len(dependents)):
            if dependents[i][0].decode() == dependent_variable:
                self._dependent_variable_values = data[:, i + len(independents)]
                break

    def plot(self, fit=True, save=False, save_path=None, save_name=None):
        #need to edit
        fig = plt.figure()
        plt.title('Q{:1d}: T1'.format(self._qubit_id))
        plt.plot(self._idle_gate_times, self._dependent_variable_values, marker='.', linestyle='None')  # (b-0.062)*(490)*(0.02)*484
        plt.xlabel('Idle Gate Time (us)')
        plt.ylabel(self._dependent_variable_name)
        if fit:
            fit_parameters = self.fit(save)
            plt.title('Q{:1d}: T1 = {:2d}us'.format(self._qubit_id,int(-1/fit_parameters[1])))
            plt.plot(self._idle_gate_times, [fit_parameters[0]*np.exp(fit_parameters[1]*t)+fit_parameters[2] for t in self._idle_gate_times])
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

    def fit(self, save):
        popt, pcov = curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, self._idle_gate_times, self._dependent_variable_values,p0=[1,(-1/15),0])
        if save:
            self.update_dataChest(int(-100/popt[1])/100)
        return popt

    def update_dataChest(self, value):
        d = dataChest(self._path)
        d.openDataset(self._file, modify = True)
        d.addParameter('Fit T1', value, 'us', overwrite = True)
        d.addParameter('Fit Date Stamp', datetime.now().strftime("%Y-%m-%d"), overwrite = True)

class T2(object):

    def __init__(self, path, file, dependent_variable='Single Shot Occupation'):
        pass

class P1(object):

    def __init__(self, path, file, dependent_variable='Single Shot Occupation'):
        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        data = d.getData()
        independents = variables[0]
        dependents = variables[1]
        #there absolutely has to be a better way of doing this, but:
        self._file = file
        self._biases = data[:,0]
        self._dependent_variable_values = []
        self._dependent_variable_name = dependent_variable
        for i in range(len(dependents)):
            if dependents[i][0].decode() == dependent_variable:
                self._dependent_variable_values = data[:, i + len(independents)]
                break

    def plot_versus_radiator_frequency(self, bias_offset=0, bias_conversion=2*10*484, save=False, save_path=None, save_name=None):
        #need to edit
        fig = plt.figure()
        plt.title('P1 vs Radiator Bias')
        plt.semilogy([(b-bias_offset)*bias_conversion for b in self._biases], self._dependent_variable_values, marker='.', linestyle='None')  # (b-0.062)*(490)*(0.02)*484
        plt.xlabel('Radiator Bias (GHz)')
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

    def plot_versus_radiator_flux(self, bias_conversion=2*10*484, save=False, save_path=None, save_name=None):
        #need to edit
        fig = plt.figure()
        plt.title('P1 vs Radiator Flux')
        plt.semilogy([(b)*bias_conversion for b in self._biases], self._dependent_variable_values, marker='.', linestyle='None')  # (b-0.062)*(490)*(0.02)*484
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

    def __init__(self, path, file, rad_id=None):
        d = dataChest(path)
        d.openDataset(file)
        variables = d.getVariables()
        data = d.getData()
        data = data.transpose()
        self._iterations = int(np.amax(data[0]) + 1)  # Number of Iterations
        self._reps = int(np.amax(data[1]) + 1)  # Number of Reps
        self._data = []
        self._fit_psd = None
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
        self._bias = d.getParameter('JB{:1d} Voltage Bias'.format(rad_id))[0]
        try:
            self._repetition_rate = 1e9/d.getParameter('Time Per Iteration')[0]
        except:
            self._repetition_rate = 1e9/100000
        self._psd()

    def _psd(self):
        state = self._data
        n_max = self._iterations

        psd = [None] * n_max
        for i in range(n_max):
            freqs, psd[i] = signal.periodogram(state[i], fs=self._repetition_rate, return_onesided=True)
        ave_psd = np.average(psd, axis=0)
        self._psd = ave_psd[1:]
        self._f = freqs[1:]

    def psd(self):
        return self._f, self._psd

    def fit_psd(self):
        if self._fit_psd is None:
            self.fit()
        return self._f, self._fit_psd

    def fit(self, save=False):
        fs=self._repetition_rate
        def fit_PSD_target_function(f, T_parity, F_map):
            return 2 * (4 * F_map ** 2 / T_parity) / ((2 / T_parity) ** 2 + (2 * np.pi * f) ** 2) + 2 * (
                        1 - F_map ** 2) / fs

        initial_guess_fidelity = np.arange(0.1, 1, 0.05)
        initial_guess_parity = np.logspace(-4, -1, 5)
        covariance = float('inf')
        best_r_squared = float('inf')
        freqs,ave_psd = self.psd()
        for igf in initial_guess_fidelity:
            for igp in initial_guess_parity:
                params_curr, params_covariance_curr = curve_fit(
                    fit_PSD_target_function, freqs, ave_psd, p0=[igp, igf], bounds=([0.0001, 0], [0.1, 1]),
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
            if not os.path.exists(DATACHEST_ROOT+save_path+'\\Figures\\'):
                os.makedirs(DATACHEST_ROOT+save_path+'\\Figures\\')
            plt.savefig(DATACHEST_ROOT+save_path+ '\\Figures\\' + save_name)
            plt.close(fig)
        else:
            plt.ion()
            plt.show()


    def update_dataChest(self):
        d = dataChest(self._path)
        d.openDataset(self._file, modify = True)
        d.addParameter('Fit Parity Rate', self._parity_rate, 'Hz', overwrite = True)
        d.addParameter('Fit Fidelity', self._fit_fidelity, overwrite=True)
        d.addParameter('Fit Date Stamp', datetime.now().strftime("%Y-%m-%d"), overwrite = True)
        d.addParameter('Time Per Iteration', int(1e9/self._repetition_rate), overwrite=True)