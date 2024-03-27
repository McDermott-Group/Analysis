import os
import numpy as np
import noiselib
import importlib
importlib.reload(noiselib)
from dataChest import *
import matplotlib.pyplot as plt
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
from scipy.stats import norm
from dateStamp import dateStamp
import glob

class ChargeOffset(object):

    def __init__(self):
        self.file_tags = []
        self.labels = {}
        self.time = {}
        self.abs_time = {} #actual utc time, not relative to t0 like ^
        self.unwrapped_charge = {}
        self.fit_R2 = {}
    
    def analyze_charge_jump_correlations(self, file_path):
        pass
    
    def add_dataset(self, file_path):
        fpath, fname = os.path.split(file_path)
        ftag = fname.split('_')[0]
        self.file_tags += [ftag]
        dc = dataChest(fpath.split('\\'))
        dc.openDataset(fname)
        
        varsList = dc.getVariables()
        timeVar = varsList[0][0][0]
        offsetVars = [var[0] for var in varsList[1] if 'R2' not in var[0]]
        data = dc.getData()
        data = data.transpose()
        labels = [l[-2:] for l in offsetVars]
        if labels == ['ge']:
            labels = [fpath.split('\\')[3]]
        
        abs_time = data[0]
        time = data[0] - data[0][0]
        period2e = {l: dc.getParameter('2e Period {}'.format(l)) for l in labels}
        unwrapped_charge = { l: noiselib.unwrap_voltage_to_charge(
                                data[2*i+1], period2e[l]/2, 2/period2e[l] )
                                for i,l in enumerate(labels) }
        unwrapped_charge = { l: unwrapped_charge[l] - unwrapped_charge[l][0] 
                                for l in labels } # zero the trace
        fit_R2 = { l: data[2*i+2] for i,l in enumerate(labels) }
        # charge_jump_sizes = {l: unwrapped_charge[l][1:] - unwrapped_charge[l][:-1]
                                # for l in labels}
        # large_jump[ftag]s = {l:charge_jump_sizes[l] > 0.1 for l in labels}
        
        self.labels[ftag] = labels
        self.time[ftag] = time
        self.abs_time[ftag] = abs_time
        self.unwrapped_charge[ftag] = unwrapped_charge
        self.fit_R2[ftag] = fit_R2
        self._last_dataset_added = ftag
        return self
    
    def limit_dataset(self, label, dataset=None, start=None, end=None):
        """limits the data for label by setting the values in [start,end) to nan"""
        if dataset is None:
            dataset = self._last_dataset_added
        self.unwrapped_charge[dataset][label][start:end] = np.nan
        self.fit_R2[dataset][label][start:end] = np.nan
        return self
        
    def _above(self, array, threshold=0.1, abs=True):
        """Equivilent to array > threshold, but deals with nans"""
        array = array.copy()
        if abs:
            jump_size = np.abs(array, array, where=~np.isnan(array))
        else:
            jump_size = array
        is_above = np.greater( jump_size, threshold, jump_size, where=~np.isnan(array) )
        is_above[np.isnan(is_above)] = False
        return is_above.astype(bool)
    
    def get_charge_offset(self, datasets = None):
        if datasets is None:
            datasets = self.file_tags
        time = np.array([0])
        offset = {}
        for ftag in datasets:
            n_old = time.size
            time = np.append( time, time[-1] + self.time[ftag] )
            n_new = self.time[ftag].size
            for l in self.labels[ftag]:
                # fill with nan for length of other traces
                new_charge = self.unwrapped_charge[ftag][l]
                if l not in offset:
                    offset[l] = np.append( np.full(n_old-1, np.nan), [0] )
                lastValIndex = np.argwhere(~np.isnan(offset[l]))[-1][0]
                offset[l] = np.append( offset[l], offset[l][lastValIndex] + new_charge )
            for l in set(offset.keys()) - set(self.labels[ftag]):
                offset[ftag] = np.append( offset[ftag], np.full(n_new, np.nan) )
        return time, offset
    
    def get_jump_sizes(self, datasets = None, plot=False, ax=None, qubits=None):
        time, offset = self.get_charge_offset(datasets)
        if qubits is None:
            qubits = list(offset.keys())
        jumps = {l: offset[l][1:] - offset[l][:-1] for l in list(offset.keys()) if l in qubits}
        
        # find 2*sigma jump size for each label
        # def gaus(x, b, A, mu, sigma):#, A2, mu2, sigma2):
            # return b + A * exp(-(x-mu)**2/(2*sigma**2))# + A2 * exp(-(x-mu2)**2/(2*sigma2**2))
        def gaus(x, b, A, sigma):
            mu = 0
            return b + A * exp(-(x-mu)**2/(2*sigma**2))
        sigma = {}
        for l in list(jumps.keys()):
            if l not in qubits:
                continue
            h, bins = np.histogram(jumps[l][~np.isnan(jumps[l])], bins=500, range=(-0.5,0.5))
            # h = noiselib.movingmean(h, 10)
            x = (bins[1:]+bins[:-1])/2
            popt, pcov = curve_fit(gaus, x, h, p0=[0,h.max(),0.1])
            # popt, pcov = curve_fit(gaus, x, h, p0=[0,h.max(),0,0.1,h.max(),0.005,0.3])
            sigma[l] = popt[-1]
            
            if plot:
                if ax is None:
                    # fig = plt.figure()
                    # axi = fig.add_subplot(111)
                    fig, axi = plt.subplots(1,1)
                else:
                    axi = ax
                axi.set_title('Charge Jumps {}'.format(l))
                axi.set_xlabel('Jump Size {} [e]'.format(l))
                axi.set_ylabel('')
                center = (bins[:-1] + bins[1:])/2
                widths = np.diff(bins)
                # center2 = center**2 * np.sign(center)
                # widths = np.diff(bins**2)
                axi.step(center, h)#, color='C2', width=widths)
                # axi.bar(center, h-gaus(center, 0, *popt[1:]), 
                        # color='C3', width=widths)
                axi.plot(center, gaus(center, 0, *popt[1:]), 'k-' )
                axi.set_yscale('log')
                axi.set_ylim([10e-1, 1.5*h.max()])
                plt.draw()
                plt.pause(0.05)
        
        return jumps, sigma#{'Q1':0.05, 'Q2':0.05, 'Q3':0.05, 'Q4':0.05}#sigma
        
    def plot_charge_offset(self, datasets = None, ax=None):
        time, offset = self.get_charge_offset(datasets)
        jumps, sigma = self.get_jump_sizes(datasets)
        largeJumps = {l: np.append(self._above(jumps[l], 2*sigma[l]), True) | 
                         np.append(True, self._above(jumps[l], 2*sigma[l]))
                        for l in list(jumps.keys())}
        
        if ax is None:
            fig, ax = plt.subplots(1,1)
        ax.set_title('Charge Offset')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Charge Offset [e]')
        for l in list(offset.keys()):
            ax.plot(time, offset[l], label=l)
            ax.plot(time, np.ma.masked_where(~largeJumps[l], offset[l]), label=l)
        ax.legend()
        plt.draw()
        plt.pause(0.05)
        
    def plot_jump_sizes(self, datasets = None, ax=None):
        jumps, sigma = self.get_jump_sizes(datasets)
        
        for l in list(jumps.keys()):
            n_tot = np.sum(~np.isnan(jumps[l]))
            n_big = np.sum(self._above(jumps[l], 2*sigma[l]))
            print('{}: {} / {} = {:.2f}% = {:.2f}*T between jumps'.format(
                    l, n_big, n_tot, 100.*n_big/n_tot, 1.*n_tot/n_big))
        print( {l: '{:.3f}'.format(2*s) for l,s in list(sigma.items())} )
        
        if ax is None:
            fig, ax = plt.subplots(1,1)
        ax.set_title('Jump Sizes')
        ax.set_xlabel('N')
        ax.set_ylabel('Jump Size ($e$)')
        for l in list(jumps.keys()):
            ax.plot(np.abs(jumps[l], jumps[l], where=~np.isnan(jumps[l])), label=l)
        ax.legend()
        plt.draw()
        plt.pause(0.05)
        
    def plot_charge_correlation(self, label1, label2, thresh=None, 
                                      datasets=None, ax=None, plot=True):
        jumps, sigma = self.get_jump_sizes(datasets)
        jumps1 = jumps[label1]
        jumps2 = jumps[label2]
        dt = self.get_time_steps(datasets)
        tA = np.nansum( dt[np.isfinite(jumps1)] )
        tB = np.nansum( dt[np.isfinite(jumps2)] )
        
        if thresh is None:
            thresh = (2*sigma[label1], 2*sigma[label2])
        
        corrJumps = self._above(jumps1, thresh[0]) & self._above(jumps2, thresh[1])
        eitherJumps = self._above(jumps1, thresh[0]) | self._above(jumps2, thresh[1])
        bothMeas = np.isfinite(jumps1) & np.isfinite(jumps2)
        t = np.nanmean(dt[bothMeas])
        d_t = 0#np.nanstd(dt[bothMeas])
        
        nA, nB = np.sum(np.isfinite(jumps1)), np.sum(np.isfinite(jumps2))
        jA = np.sum( self._above( jumps1,thresh[0],abs=True) )
        jB = np.sum( self._above( jumps2,thresh[1],abs=True) )
        pA = 1. * jA / nA
        d_pA = 1. * np.sqrt(jA) / nA
        pB = 1. * jB / nB
        d_pB = 1. * np.sqrt(jB) / nB
        # gammas (rates)
        gA, gB = 1.*jA/tA, 1.*jB/tB
        d_gA, d_gB = gA*np.sqrt(1./jA+(d_t/t)**2), gB*np.sqrt(1./jB+(d_t/t)**2)
        gAB = (gA+gB)/2.
        d_gAB = np.sqrt( d_gA**2 + d_gB**2 )
        
        jumps1[np.isnan(jumps2)] = np.nan
        jumps2[np.isnan(jumps1)] = np.nan
        print( 'Correlated Jumps: {}-{}'.format(label1,label2) )
        print( '    {} / {} = {:.2f}%'.format( 
                np.sum(corrJumps), np.sum(eitherJumps), 
                100.*np.sum(corrJumps)/np.sum(eitherJumps) ))
        # for l,j in [(label1,jumps1),(label2,jumps2)]:
            # print( '    {} / {} = {:.2f}%'.format( 
                    # np.sum(corrJumps), 
                    # np.sum(self._above(j,2*sigma[l])),
                    # 100.*np.sum(corrJumps)/np.sum(self._above(j,2*sigma[l])) ) )
        q1 = np.sum(self._above( jumps1,thresh[0],abs=False) & 
                    self._above( jumps2,thresh[1],abs=False))
        q2 = np.sum(self._above(-jumps1,thresh[0],abs=False) & 
                    self._above( jumps2,thresh[1],abs=False))
        q3 = np.sum(self._above(-jumps1,thresh[0],abs=False) & 
                    self._above(-jumps2,thresh[1],abs=False))
        q4 = np.sum(self._above( jumps1,thresh[0],abs=False) & 
                    self._above(-jumps2,thresh[1],abs=False))
        a1324 = 1.*(q1+q3)/(q2+q4)
        d_a1324 = a1324 * np.sqrt( 1./(q1+q3) + 1./(q2+q4) )
        print('    13/24 = ({}+{})/({}+{}) = {:.2f} \u00B1 {:.3f}'.format(
                    q1, q3, q2, q4, a1324, d_a1324  ))
        
        # nA, nB = np.sum(np.isfinite(jumps1)), np.sum(np.isfinite(jumps2))
        # jA = np.sum( self._above( jumps1,thresh[0],abs=True) )
        # jB = np.sum( self._above( jumps2,thresh[1],abs=True) )
        # pA = 1. * jA / nA
        # d_pA = 1. * np.sqrt(jA) / nA
        # pB = 1. * jB / nB
        # d_pB = 1. * np.sqrt(jB) / nB
        # pA,d_pA = 0.05549, 0.00244  # Q1
        # pB,d_pB = 0.06050, 0.00278  # Q2
        # pA,d_pA = 0.06049, 0.00256  # Q3
        # pB,d_pB = 0.06035, 0.00269  # Q4
        p_obs = 1.*(q1 + q2 + q3 + q4) / np.sum( bothMeas )
        d_p_obs = p_obs * np.sqrt( 1./(q1 + q2 + q3 + q4) + 1./np.sum( bothMeas ) )
        pC = 1.*(p_obs-pA*pB) / (1.+p_obs-pA-pB)
        d_pC_top = np.sqrt( d_p_obs**2 + (pA*pB)**2*((d_pA/pA)**2+(d_pB/pB)**2) )
        d_pC_bot = np.sqrt( d_p_obs**2 + d_pA**2 + d_pB**2 )
        d_pC = pC * np.sqrt( (d_pC_top/(p_obs-pA*pB))**2 + (d_pC_bot/(1+p_obs-pA-pB))**2 )
        pAp = (pA-pC)/(1-pC)
        d_pAp = pAp * np.sqrt( (d_pA**2+d_pC**2)/(pA-pC)**2 + d_pC**2/pC**2 )
        pBp = (pB-pC)/(1-pC)
        d_pBp = pBp * np.sqrt( (d_pB**2+d_pC**2)/(pB-pC)**2 + d_pC**2/pC**2 )
        pCp = 1.*pC/np.mean([pA,pB])
        d_pCp = 1.*pC/np.mean([pA,pB])*np.sqrt( (d_pC/pC)**2 + (d_pA**2+d_pB**2)/(pA+pB)**2 )
        
        # gammas (rates)
        gC = pCp*gAB
        d_gC = gC * np.sqrt( (d_pCp/pCp)**2 + (d_gAB/gAB)**2 )
        
        # print('    pA,pB,p_obs = {:.2f},{:.2f},{:.2f}'.format(pA,pB,p_obs))
        print('    pA = {}/{} = {:.5f} \u00B1 {:.5f}'.format( jA, nA, pA, d_pA ))
        print('    pB = {}/{} = {:.5f} \u00B1 {:.5f}'.format( jB, nB, pB, d_pB ))
        print("    pA' = {:.4f} \u00B1 {:.3f}".format(pAp,d_pAp))
        print("    pB' = {:.4f} \u00B1 {:.3f}".format(pBp,d_pBp))
        print('    p_obs = {}/{} = {:.4f} \u00B1 {:.3f}'.format( (q1 + q2 + q3 + q4),
                                            np.sum( bothMeas ), p_obs, d_p_obs ))
        print('    pC = {:.3f} \u00B1 {:.4f}'.format(pC, d_pC))
        print('    pC/mean(pA,pB) = {:.3f} \u00B1 {:.4f}'.format( pCp, d_pCp ))
        print('    gA,gB = {:.3f} \u00B1 {:.4f} mHz, {:.3f} \u00B1 {:.4f} mHz'.format( 
                           1e3*gA, 1e3*d_gA, 1e3*gB, 1e3*d_gB ))
        print('    gC = {:.3f} \u00B1 {:.4f} mHz'.format( 1e3*gC, 1e3*d_gC ))
        
        if plot:
            if ax is None:
                fig, ax = plt.subplots(1,1)
            ax.set_aspect(1)
            ax.set_xlim(-0.5,0.5)
            ax.set_ylim(-0.5,0.5)
            ax.set_title('Charge Jumps')
            ax.set_xlabel('Jump Size {} [e]'.format(label1))
            ax.set_ylabel('Jump Size {} [e]'.format(label2))
            # ax.scatter(jumps1, jumps2, c=corrJumps, cmap='rainbow', marker='.')
            ax.scatter(jumps1, jumps2, c='k', marker='.', s=4)
            plt.draw()
            plt.pause(0.05)
        return pCp, d_pCp, a1324, d_a1324, gA, d_gA, gB, d_gB
        
    def plot_time_steps(self, datasets=None, ax=None):
        dt = self.get_time_steps(datasets)
        
        if ax is None:
            fig, ax = plt.subplots(1,1)
        ax.set_title('Time Steps')
        ax.set_xlabel('N')
        ax.set_ylabel('Step Size [s]')
        ax.plot(dt)
        plt.draw()
        plt.pause(0.05)
    
    def get_time_steps(self, datasets=None):
        if datasets is None:
            datasets = self.file_tags
        dt = np.array([])
        for ftag in datasets:
            dt = np.append(dt, [np.nan])
            dt = np.append(dt, np.diff(self.time[ftag]))
            
        print('Average Measurement Interval T = {:.2f}'.format(np.nanmean(dt)))
        
        return dt
        
    def get_files_triggered_on_bad_fit(self, dataset, label, path, thresh=0.9):
        """Takes in a specific dataset tag and the path to the folder holding the
        files that were fit to get each individual offset curve."""
        large_err_indicies = np.argwhere(self.fit_R2[dataset][label] < thresh)
        times = self.abs_time[dataset][large_err_indicies]
        times = times.transpose()[0]
        return self._files_closest_to_times_mat(path, times)
        
    
    def _files_closest_to_times_mat(self, paths, times):
        if not isinstance(paths, list):
            paths = [paths]
        all_files = []
        for path in paths:
            all_files += glob.glob(path+'*.mat')
        all_times = [os.path.getctime(f) for f in all_files]
        all = list(zip(all_times, all_files))
        all.sort()
        all_times, all_files = list(zip(*all)) # undoes zip, now sorted
        closest_time_index = np.searchsorted(all_times, times)
        return [all_files[i] for i in closest_time_index]
        
    
    def _files_closest_to_times_hdf5(self, path, times):
        dc = dataChest(path.split('\\'))
        all_files,_ = dc.ls()
        all_files.sort()
        all_tags = [fname[:10] for fname in all_files]
        DS = dateStamp()
        all_times = [DS.utcDateStrToFloat(DS.invertDateStamp(tag)[:26])
                        for tag in all_tags]
        closest_time_index = np.searchsorted(all_times, times)
        return [all_files[i] for i in closest_time_index.transpose()[0]]

