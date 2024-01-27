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
    def __init__(self, path, files, device_names, gap=190e-6):
        d = dataChest(path)
        self.I = []
        self.V = []
        self.R = 0
        self.files = []
        self.device_names = []
        for i,file in enumerate(files):
            d.openDataset(file)
            data = d.getData(variablesList=['Current','Voltage'])
            self.SeriesResistance = d.getParameter('AC Resistance In [kOhms]')*1e3
            self.I.append(data[:,0])
            self.V.append(data[:,1])
            self.files.append(files[i])
            self.device_names.append(device_names[i])
        self.path = path
        self.gap = gap

    def plotLogIvsV(self, save=False, save_path=None, save_name=None,autocenter_mode='supercurrent'):
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

    def plotIvsV(self, save=False, save_path=None, save_name=None,autocenter_mode = 'supercurrent'):
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

    def autocenter(self,mode='supercurrent'):
        for idx, file in enumerate(self.files):
            self.V[idx] = self.V[idx] - np.mean(self.V[idx])
            v = self.V[idx]
            i = self.I[idx]
            i = np.abs(i)
            # offset = v[np.argmin(i)]
            # v = v - offset
            if mode == 'supercurrent' or mode == 'both':
                offset = self.vertSearch(v, i, -1.25*self.gap, 1.25*self.gap)
                v = v - offset

            if mode == 'gap' or mode == 'both':
                a = 4.5 if mode == 'gap' else 3.0
                b = 0.1 if mode == 'gap' else 1.0
                gap1 = self.vertSearch(v, i, -a*self.gap, -b*self.gap)
                gap2 = self.vertSearch(v, i, b*self.gap, a*self.gap)
                avg_factor = 2 if mode == 'gap' else 3
                v = v -(gap1+gap2)/avg_factor #This is the average of the THREE vertical search results
            # plt.figure()
            # plt.axvline(0)
            # plt.axvline(gap1-(gap1+gap2)/2)
            # plt.axvline(gap2-(gap1+gap2)/2)
            # v=v-(gap1+gap2)/2
            a,b = self.fitLine(v,i,np.min(v),0.99*np.min(v))
            c,d = self.fitLine(v,i,0.99*np.max(v),np.max(v))
            self.R = 1/np.min(np.abs([a,c]))
            # plt.plot([np.min(v),np.max(v)],[a*np.min(v)+b,a*np.max(v)+b])
            # c,d = self.fitLine(v,i,2.5*self.gap,5*self.gap,slope=a)
            # plt.plot([np.min(v),np.max(v)],[c*np.min(v)+d,c*np.max(v)+d])
            # plt.plot(v,i)
            # plt.show()
            self.V[idx] = v

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
        cal_data = np.loadtxt('/Volumes/smb/mcdermott-group/data/windowJJs/Calibrations/20231201_cal.csv',delimiter=',')
        v_cal = cal_data[:,0].T
        i_cal = cal_data[:, 1].T
        posSlope,posInt = self.fitReferenceLine(v_cal[v_cal > 1.5],i_cal[v_cal > 1.5])
        negSlope, negInt = self.fitReferenceLine(v_cal[v_cal < -1.5],i_cal[v_cal < -1.5])
        i_cal = [negSlope*v+negInt for v in np.arange(-20, np.min(v_cal), 0.1)] + list(i_cal) + [posSlope*v+posInt for v in np.arange(np.max(v_cal),20, 0.1)]
        v_cal = [v for v in np.arange(-20, np.min(v_cal), 0.1)] + list(v_cal) + [v for v in np.arange(np.max(v_cal), 20, 0.1)]

        for idx, file in enumerate(self.files):
            v_out = self.I[idx]*self.SeriesResistance
            i_mapped = np.interp(v_out,v_cal,i_cal)
            self.I[idx] = i_mapped