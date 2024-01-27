from datetime import datetime
import numpy as np
from numpy import diff
from scipy import signal
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
from dataChest import *
import platform

DATACHEST_ROOT = '/Volumes/smb/mcdermott-group/data/' if 'macOS' in platform.platform() else 'Z:\\mcdermott-group\\data\\'

class JJ(object):
    """"
        Class for analyzing JJ IV Sweeps
    """
    def __init__(self, path, files, device_names, gap=190e-6):
        """Instantiate JJ for files associated with given devices.

        :param path: path containing the IV files.
        :param files: hdf5 files to open in data chest.
        :param device_names: human-readable device identifiers (mostly used for plotting).
        :param gap: Delta/e, in Volts
        """
        d = dataChest(path)
        self.I = []
        self.V = []
        self.R = 0
        self.files = []
        self.device_names = []
        for i, file in enumerate(files):
            d.openDataset(file)
            data = d.getData(variablesList=['Current','Voltage'])
            self.Rs = d.getParameter('AC Resistance In [kOhms]') * 1e3
            self.I.append(data[:,0])
            self.V.append(data[:,1])
            self.files.append(files[i])
            self.device_names.append(device_names[i])
        self.path = path
        self.gap = gap


    def plot_I_vs_V(self, save=False, save_path=None, save_name=None, autocenter_mode='supercurrent', remove_jumps = False):
        """Plot the centered IV on a linear scale.

        :param save: specifies whether to save the plot. When false, plot will be displayed in interactive mode.
        :param save_path:  path for saving the plot. When None, defaults to '../Figures/'.
        :param save_name: name for saving the plot. When None, defaults to 'Log_' + the hdf5 file name.
        :param autocenter_mode: 'supercurrent','gap', or 'both'. Specifies whether zero voltage is set by using the
                                vertical feature at V=0, at twice the gap, or both.
        :param remove_jumps: remove the jump from the supercurrent branch to the resistive branch?
        :return: None
        """
        if remove_jumps:
            self._plot_I_vs_V_removing_jumps()
        else:
            self._plot_I_vs_V_raw()

    def plot_log_I_vs_V(self, save=False, save_path=None, save_name=None, autocenter_mode='supercurrent', remove_jumps=False):
        """Plot the centered IV on a semilog(y) scale.

        :param save: specifies whether to save the plot. When false, plot will be displayed in interactive mode.
        :param save_path:  path for saving the plot. When None, defaults to '../Figures/'.
        :param save_name: name for saving the plot. When None, defaults to 'Log_' + the hdf5 file name.
        :param autocenter_mode: 'supercurrent','gap', or 'both'. Specifies whether zero voltage is set by using the
                                vertical feature at V=0, at twice the gap, or both.
        :param remove_jumps: remove the jump from the supercurrent branch to the resistive branch?
        :return: None
        """
        if remove_jumps:
            self._plot_I_vs_V_removing_jumps()
        else:
            self._plot_I_vs_V_raw()
        if remove_jumps:
            self._plot_log_I_vs_V_removing_jumps()
        else:
            self._plot_log_I_vs_V_raw()

    def _plot_log_I_vs_V_raw(self, save=False, save_path=None, save_name=None, autocenter_mode='supercurrent'):
        """Plot the centered IV on a semilog(y) scale.

        :param save: specifies whether to save the plot. When false, plot will be displayed in interactive mode.
        :param save_path:  path for saving the plot. When None, defaults to '../Figures/'.
        :param save_name: name for saving the plot. When None, defaults to 'Log_' + the hdf5 file name.
        :param autocenter_mode: 'supercurrent','gap', or 'both'. Specifies whether zero voltage is set by using the
                                vertical feature at V=0, at twice the gap, or both.
        :return: None
        """
        self.autocenter(mode=autocenter_mode)
        for idx, file in enumerate(self.files):
            fig = plt.figure()
            plt.plot(self.V[idx],[np.log10(i) if i > 0 else np.log10(-1*i) for i in self.I[idx]])
            plt.title(file+ ' ({0})'.format(self.device_names[idx]))
            plt.xlabel('Voltage (V)')
            plt.ylabel('Log |Current| (Log(A))i')
            if save:
                if save_path is None:
                    save_path = os.path.join(*self.path[:-1])
                if save_name is None:
                    save_name = file.replace('.', '_')
                plt.ioff()
                if not os.path.exists(DATACHEST_ROOT + save_path + '//Figures//'):
                    os.makedirs(DATACHEST_ROOT + save_path + '//Figures//')
                plt.savefig(DATACHEST_ROOT + save_path + '//Figures//Log_' + save_name)
                plt.close(fig)
            else:
                plt.ion()
                plt.show()

    def _plot_I_vs_V_raw(self, save=False, save_path=None, save_name=None, autocenter_mode ='supercurrent'):
        """Plot the centered IV on a linear scale.

        :param save: specifies whether to save the plot. When false, plot will be displayed in interactive mode.
        :param save_path:  path for saving the plot. When None, defaults to '../Figures/'.
        :param save_name: name for saving the plot. When None, defaults to the hdf5 file name.
        :param autocenter_mode: 'supercurrent','gap', or 'both'. Specifies whether zero voltage is set by using the
                                vertical feature at V=0, at twice the gap, or both.
        :return: None
        """
        self.autocenter(mode=autocenter_mode)
        for idx, file in enumerate(self.files):
            fig = plt.figure()
            plt.axvline(0,color='black')
            plt.axhline(0,color='black')
            plt.plot(self.V[idx], self.I[idx])
            plt.title(file + ' ({:<})\n R={:2d}$\Omega$'.format(self.device_names[idx],int(self.R)))
            plt.xlabel('Voltage (V)')
            plt.ylabel('Current (A)')
            if save:
                if save_path is None:
                    save_path = os.path.join(*self.path[:-1])
                if save_name is None:
                     save_name = file.replace('.', '_')
                plt.ioff()
                if not os.path.exists(DATACHEST_ROOT + save_path + '//Figures//'):
                    os.makedirs(DATACHEST_ROOT + save_path + '//Figures//')
                plt.savefig(DATACHEST_ROOT + save_path + '//Figures//' + save_name)
                plt.close(fig)
            else:
                plt.ion()
                plt.show()

    def _plot_I_vs_V_removing_jumps(self, save=False, save_path=None, save_name=None, autocenter_mode='supercurrent'):
        self.autocenter(mode=autocenter_mode)
        for idx, file in enumerate(self.files):
            v = self.V[idx]
            i = self.I[idx]
            fig = plt.figure()
            plt.axvline(0)
            # to make arrays of all the local mins and maxes
            v_localmin = []
            v_localmax = []
            for j in range(len(v)-1):
                if v[j] < -1.3*self.gap and v[j-1] > v[j] and v[j] < v[j+1]:
                    v_localmin.append(v[j])
                elif v[j] > 1.3*self.gap and v[j-1] < v[j] and v[j] > v[j+1]:
                    v_localmax.append(v[j])
                else:
                    pass
            # to make data for negative voltage
            for l in range(len(v_localmin)):
                v_neg = []
                i_neg = []
                # fill up v_neg starting from local minimum going until zero for each cycle in the data
                for j in range(len(v)):
                    if v_localmin[l] == v[j]:
                        v_short = v[j:(len(v)-1)]
                        for k in range(len(v_short)):
                            if v_short[k] < 0.000005 and v_short[k] > -0.000005:
                                v_neg = v[j:(j+k)]
                                i_neg = i[j:(j+k)]
                                break
                plt.plot(v_neg, i_neg, color='blue')
            # to make data for positive voltage
            for l in range(len(v_localmax)):
                v_pos = []
                i_pos = []
                for j in range(len(v)):
                    if v_localmax[l] == v[j]:
                        v_short2 = v[j:(len(v)-1)]
                        for k in range(len(v_short2)):
                            if v_short2[k] < 0.000005 and v_short2[k] > -0.000005:
                                v_pos = v[j:(j+k)]
                                i_pos = i[j:(j+k)]
                                break
                plt.plot(v_pos, i_pos, color='blue')
            # to make data for zero voltage vertical
            v_vert = []
            i_vert = []
            for l in range(len(v)):
                if v[l] < 0.000001 and v[l] > -0.000001:
                    v_vert.append(v[l])
                    i_vert.append(i[l])
            plt.plot(v_vert, i_vert, color='blue')

            plt.title(file + ' ({0})'.format(self.device_names[idx]))
            plt.xlabel('Voltage (V)')
            plt.ylabel('Current (A)')
            if save:
                if save_path is None:
                    save_path = os.path.join(*self.path[:-1])
                if save_name is None:
                     save_name = file.replace('.', '_')
                plt.ioff()
                if not os.path.exists(DATACHEST_ROOT + save_path + '//Figures//'):
                    os.makedirs(DATACHEST_ROOT + save_path + '//Figures//')
                plt.savefig(DATACHEST_ROOT + save_path + '//Figures//' + save_name)
                plt.close(fig)
            else:
                plt.ion()
                plt.show()

    def _plot_log_I_vs_V_removing_jumps(self, save=False, save_path=None, save_name=None, autocenter_mode='supercurrent'):
        self.autocenter(mode=autocenter_mode)
        for idx, file in enumerate(self.files):
            v = self.V[idx]
            i = self.I[idx]
            fig = plt.figure()
            plt.axvline(0)
            # to make arrays of all the local mins and maxes
            v_localmin = []
            v_localmax = []
            for j in range(len(v) - 1):
                if v[j] < -1.3 * self.gap and v[j - 1] > v[j] and v[j] < v[j + 1]:
                    v_localmin.append(v[j])
                elif v[j] > 1.3 * self.gap and v[j - 1] < v[j] and v[j] > v[j + 1]:
                    v_localmax.append(v[j])
                else:
                    pass
            # to make data for negative voltage
            for l in range(len(v_localmin)):
                v_neg = []
                i_neg = []
                # fill up v_neg starting from local minimum going until zero for each cycle in the data
                for j in range(len(v)):
                    if v_localmin[l] == v[j]:
                        v_short = v[j:(len(v) - 1)]
                        for k in range(len(v_short)):
                            if v_short[k] < 0.000005 and v_short[k] > -0.000005:
                                v_neg = v[j:(j + k)]
                                i_neg = i[j:(j + k)]
                                break
                plt.plot(v_neg, [np.log10(a) if a > 0 else np.log10(-1*a) for a in i_neg], color='blue')  # color='blue')
            # to make data for positive voltage
            for l in range(len(v_localmax)):
                v_pos = []
                i_pos = []
                for j in range(len(v)):
                    if v_localmax[l] == v[j]:
                        v_short2 = v[j:(len(v) - 1)]
                        for k in range(len(v_short2)):
                            if v_short2[k] < 0.000005 and v_short2[k] > -0.000005:
                                v_pos = v[j:(j + k)]
                                i_pos = i[j:(j + k)]
                                break
                plt.plot(v_pos, [np.log10(a) if a > 0 else np.log10(-1 * a) for a in i_pos], color='blue')  # color='blue')
            # to make data for zero voltage vertical
            v_vert = []
            i_vert = []
            for l in range(len(v)):
                if v[l] < 0.000001 and v[l] > -0.000001:
                    v_vert.append(v[l])
                    i_vert.append(i[l])
            plt.plot(v_vert, [np.log10(a) if a > 0 else np.log10(-1*a) for a in i_vert], color='blue')  # color='blue')
            plt.title(file+ ' ({0})'.format(self.device_names[idx]))
            plt.xlabel('Voltage (V)')
            plt.ylabel('Log |Current| (Log(A))i')
            if save:
                if save_path is None:
                    save_path = os.path.join(*self.path[:-1])
                if save_name is None:
                    save_name = file.replace('.', '_')
                plt.ioff()
                if not os.path.exists(DATACHEST_ROOT + save_path + '//Figures//'):
                    os.makedirs(DATACHEST_ROOT + save_path + '//Figures//')
                plt.savefig(DATACHEST_ROOT + save_path + '//Figures//Log_' + save_name)
                plt.close(fig)
            else:
                plt.ion()
                plt.show()

    def autocenter(self,mode='supercurrent'):
        """Modify the voltage(V) attribute of JJ object to center the supercurrent at V=0.

        :param mode: 'supercurrent','gap', or 'both'. Specifies whether zero voltage is set by using the
                    vertical feature at V=0, at twice the gap, or the average of both.
        :return: None
        """
        for idx, file in enumerate(self.files):
            self.V[idx] = self.V[idx] - np.mean(self.V[idx])
            v = self.V[idx]
            i = self.I[idx]
            i = np.abs(i)
            if mode == 'supercurrent' or mode == 'both':
                offset = self._vert_search(v, i, -0.8 * self.gap, 0.8 * self.gap)
                v = v - offset
            if mode == 'gap' or mode == 'both':
                a = 4.5 if mode == 'gap' else 3.0
                b = 0.1 if mode == 'gap' else 1.0
                gap1 = self._vert_search(v, i, -a * self.gap, -b * self.gap)
                gap2 = self._vert_search(v, i, b * self.gap, a * self.gap)
                avg_factor = 2 if mode == 'gap' else 3
                v = v -(gap1+gap2)/avg_factor

            a,b = self._fit_line(v, i, np.min(v), 0.99 * np.min(v))
            c,d = self._fit_line(v, i, 0.99 * np.max(v), np.max(v))
            self.R = 1/np.min(np.abs([a,c]))
            self.V[idx] = v

    def _vert_search(self, v, i, vmin, vmax, n_roll = 21):
        """Search for vertical features in specified portion of IV curve, used to autocenter.

        :param v: voltage.
        :param i: current.
        :param vmin: minimum voltage in span to search.
        :param vmax: maximum voltage in span to search.
        :param n_roll: number of points in rolling average to smooth curve.
        :return: voltage corresponding to steepest jump.
        """
        v_tmp = []
        i_tmp = []
        for j in range(len(v)):
            if v[j] >= vmin and v[j] <= vmax:
                v_tmp.append(v[j])
                i_tmp.append(i[j])

        v_tmp = np.convolve(v_tmp, np.ones(n_roll) / n_roll, mode='valid')
        i_tmp = np.convolve(i_tmp, np.ones(n_roll) / n_roll, mode='valid')
        didv = np.abs(np.diff(i_tmp) / np.diff(v_tmp))
        return v_tmp[np.argmax(didv)]

    def _fit_line(self, v, i, vmin, vmax, slope = None):
        """Fit a line to specified portion of IV curve

        :param v: voltage.
        :param i: current.
        :param vmin: minimum voltage in span to fit.
        :param vmax: maximum voltage in span to fit.
        :param slope: slope of line, if goal is to find intercept only.  If None, both slope and intercept will be fit.
        :return: [slope, intercept]
        """
        v_tmp = []
        i_tmp = []
        for j in range(len(v)):
            if v[j] >= vmin and v[j] <= vmax:
                v_tmp.append(v[j])
                i_tmp.append(i[j])
        def lin(x, a, b):
            return a*x + b
        def lin_fixed_slope(x, b):
            return slope*x + b

        if slope is None:
            popt, cov = curve_fit(lin, v_tmp, i_tmp)
            return popt
        else:
            popt, cov = curve_fit(lin_fixed_slope, v_tmp, i_tmp)
            return slope,popt[0]

    def _fit_reference_line(self, v, i):
        """Used to fit the calibration data.

        :param v: voltage.
        :param i: current.
        :return: [slope, intercept]
        """
        def lin(x, a, b):
            return a * x + b

        popt, cov = curve_fit(lin, v, i,p0=[2/350.0, 0.0])
        return popt
    def _diode_current_from_output_voltage(self, cal_file='/windowJJs/Calibrations/20231201_cal.csv'):
        """Correct the current (I) attribute of JJ object for the non-linearity of diode box.

        :param cal_file: File containing IV data using the diode box going into a known, linear resistance.
        :return: None
        """

        cal_data = np.loadtxt(DATACHEST_ROOT+'/windowJJs/Calibrations/20231201_cal.csv',delimiter=',')
        v_cal = cal_data[:, 0].T
        i_cal = cal_data[:, 1].T
        posSlope,posInt = self._fit_reference_line(v_cal[v_cal > 1.5], i_cal[v_cal > 1.5])
        negSlope, negInt = self._fit_reference_line(v_cal[v_cal < -1.5], i_cal[v_cal < -1.5])
        i_cal = [negSlope*v+negInt for v in np.arange(-20, np.min(v_cal), 0.1)] + list(i_cal) + [posSlope*v+posInt for v in np.arange(np.max(v_cal),20, 0.1)]
        v_cal = [v for v in np.arange(-20, np.min(v_cal), 0.1)] + list(v_cal) + [v for v in np.arange(np.max(v_cal), 20, 0.1)]

        for idx, file in enumerate(self.files):
            v_out = self.I[idx]*self.Rs
            i_mapped = np.interp(v_out,v_cal,i_cal)
            self.I[idx] = i_mapped