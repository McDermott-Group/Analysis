from datetime import datetime
import numpy as np
from numpy import diff
from scipy import signal
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
from dataChest import dataChest as dc
from dataChest import *
#DATACHEST_ROOT = '/Volumes/smb/mcdermott-group/data/'
DATACHEST_ROOT = r'Z:\\mcdermott-group\data'

class JJ(object):
    def __init__(self, path, files, device_names, gap=190e-6):
        d = dataChest(path)
        self.I = []
        self.V = []
        self.files = []
        self.device_names = []
        for i, file in enumerate(files):
            d.openDataset(file)
            variables = d.getVariables()
            data = d.getData(variablesList=['Current','Voltage'])
            self.SeriesResistance = d.getParameter('AC Resistance In [kOhms]')*1e3
            self.I.append(data[:,0])
            self.V.append(data[:,1])
            self.files.append(files[i])
            self.device_names.append(device_names[i])
        self.path = path
        self.gap = gap

    def plotLogIvsV(self, save=False, save_path=None, save_name=None):
        self.autocenter()
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

    def plotIvsV(self, save=False, save_path=None, save_name=None):
        self.autocenter()
        for idx, file in enumerate(self.files):
            fig = plt.figure()
            plt.axvline(0)
            plt.plot(self.V[idx], self.I[idx])
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

    def plotAppended(self, save=False, save_path=None, save_name=None):
        gap2 = self.Sophiaautocenter()
        for idx, file in enumerate(self.files):
            v = self.V[idx]
            i = self.I[idx]
            fig = plt.figure()
            plt.axvline(0)
            # to make arrays of all the local mins and maxes
            v_localmin = []
            v_localmax = []
            for j in range(len(v)-1):
                if v[j] < -1.3*gap2 and v[j-1] > v[j] and v[j] < v[j+1]:
                    v_localmin.append(v[j])
                elif v[j] > 1.3*gap2 and v[j-1] < v[j] and v[j] > v[j+1]:
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

    def plotAppendedLogIvsV(self, save=False, save_path=None, save_name=None):
        self.autocenter()
        gap2 = self.Sophiaautocenter()
        for idx, file in enumerate(self.files):
            v = self.V[idx]
            i = self.I[idx]
            fig = plt.figure()
            plt.axvline(0)
            # to make arrays of all the local mins and maxes
            v_localmin = []
            v_localmax = []
            for j in range(len(v) - 1):
                if v[j] < -1.3 * gap2 and v[j - 1] > v[j] and v[j] < v[j + 1]:
                    v_localmin.append(v[j])
                elif v[j] > 1.3 * gap2 and v[j - 1] < v[j] and v[j] > v[j + 1]:
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

    def autocenter(self):
        for idx, file in enumerate(self.files):
            self.V[idx] = self.V[idx] - np.mean(self.V[idx])
            v = self.V[idx]
            i = self.I[idx]
            i = np.abs(i)
            # offset = v[np.argmin(i)]
            # v = v - offset
            offset = self.vertSearch(v, i, -0.75*self.gap, 0.75*self.gap)
            v = v - offset
            gap1 = self.vertSearch(v, i, -2.5*self.gap, -1.5*self.gap)
            gap2 = self.vertSearch(v, i, 1.5*self.gap, 2.5*self.gap)
            v = v -(gap1+gap2)/3 #This is the average of the THREE vertical search results
            # plt.figure()
            # plt.axvline(0)
            # plt.axvline(gap1-(gap1+gap2)/2)
            # plt.axvline(gap2-(gap1+gap2)/2)
            # v=v-(gap1+gap2)/2
            # a,b = self.fitLine(v,i,-5*self.gap,-2.5*self.gap)
            # plt.plot([np.min(v),np.max(v)],[a*np.min(v)+b,a*np.max(v)+b])
            # c,d = self.fitLine(v,i,2.5*self.gap,5*self.gap,slope=a)
            # plt.plot([np.min(v),np.max(v)],[c*np.min(v)+d,c*np.max(v)+d])
            # plt.plot(v,i)
            plt.show()
            self.V[idx] = v

    def Sophiaautocenter(self):
        for idx, file in enumerate(self.files):
            self.V[idx] = self.V[idx] - np.mean(self.V[idx])
            v = self.V[idx]
            i = self.I[idx]
            i = np.abs(i)
            # offset = v[np.argmin(i)]
            # v = v - offset
            offset = self.vertSearch(v, i, -0.75*self.gap, 0.75*self.gap)
            v = v - offset
            gap1 = self.vertSearch(v, i, -2.5*self.gap, -1.5*self.gap)
            gap2 = self.vertSearch(v, i, 1.5*self.gap, 2.5*self.gap)
            v = v -(gap1+gap2)/3 #This is the average of the THREE vertical search results
            # plt.figure()
            # plt.axvline(0)
            # plt.axvline(gap1-(gap1+gap2)/2)
            # plt.axvline(gap2-(gap1+gap2)/2)
            # v=v-(gap1+gap2)/2
            # a,b = self.fitLine(v,i,-5*self.gap,-2.5*self.gap)
            # plt.plot([np.min(v),np.max(v)],[a*np.min(v)+b,a*np.max(v)+b])
            # c,d = self.fitLine(v,i,2.5*self.gap,5*self.gap,slope=a)
            # plt.plot([np.min(v),np.max(v)],[c*np.min(v)+d,c*np.max(v)+d])
            # plt.plot(v,i)
            plt.show()
            self.V[idx] = v
            return(gap2)

    def vertSearch(self, v, i, vmin, vmax):
        n_roll = 21
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
    def fitLine(self, v, i, vmin, vmax, slope = None):
        v_tmp = []
        i_tmp = []
        for j in range(len(v)):
            if v[j] >= vmin and v[j] <= vmax:
                v_tmp.append(v[j])
                i_tmp.append(i[j])
        def lin(x, a, b):
            return a*x + b

        def linFixedSlope(x, b):
            return slope*x + b

        if slope is None:
            popt, cov = curve_fit(lin, v_tmp, i_tmp)
            return popt
        else:
            popt, cov = curve_fit(linFixedSlope, v_tmp, i_tmp)
            return slope,popt[0]
    def removeJump(self, v, i):
        v_tmp = []
        i_tmp = []
        positive_sweep = True
        for j in range(len(v)):
            if v[j] >= 0 and v[j] <= 2*self.gap and positive_sweep:
                continue
            if v[j] <= 0 and v[j] >= 2*self.gap and not positive_sweep:
                continue
    def fitReferenceLine(self, v, i):
        def lin(x, a, b):
            return a * x + b

        popt, cov = curve_fit(lin, v, i,p0=[2/350.0, 0.0])
        return popt
    def diodeCurrentFromOutputVoltage(self):
        cal_data = np.loadtxt(DATACHEST_ROOT+'/windowJJs/Calibrations/20231201_cal.csv',delimiter=',')
        v_cal = cal_data[:, 0].T
        i_cal = cal_data[:, 1].T
        posSlope,posInt = self.fitReferenceLine(v_cal[v_cal > 1.5],i_cal[v_cal > 1.5])
        negSlope, negInt = self.fitReferenceLine(v_cal[v_cal < -1.5],i_cal[v_cal < -1.5])
        i_cal = [negSlope*v+negInt for v in np.arange(-20, np.min(v_cal), 0.1)] + list(i_cal) + [posSlope*v+posInt for v in np.arange(np.max(v_cal),20, 0.1)]
        v_cal = [v for v in np.arange(-20, np.min(v_cal), 0.1)] + list(v_cal) + [v for v in np.arange(np.max(v_cal), 20, 0.1)]

        for idx, file in enumerate(self.files):
            v_out = self.I[idx]*self.SeriesResistance
            i_mapped = np.interp(v_out,v_cal,i_cal)
            self.I[idx] = i_mapped