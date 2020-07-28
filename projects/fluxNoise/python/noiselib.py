import scipy.io as spio
import numpy as np
import os
# from Markov_Python2.analyze_QPTunneling_pomegranate import observed_to_recovered_signal
from scipy.signal import periodogram

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    Taken from https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
    '''

    def _check_keys(d):
        '''
        checks if entries in dictionary are mat-objects. If yes
        todict is called to change them to nested dictionaries
        '''
        for key in d:
            if isinstance(d[key], spio.matlab.mio5_params.mat_struct):
                d[key] = _todict(d[key])
        return d

    def _todict(matobj):
        '''
        A recursive function which constructs from matobjects nested dictionaries
        '''
        d = {}
        for strg in matobj._fieldnames:
            elem = matobj.__dict__[strg]
            if isinstance(elem, spio.matlab.mio5_params.mat_struct):
                d[strg] = _todict(elem)
            elif isinstance(elem, np.ndarray):
                d[strg] = _tolist(elem)
            else:
                d[strg] = elem
        return d

    def _tolist(ndarray):
        '''
        A recursive function which constructs lists from cellarrays
        (which are loaded as numpy ndarrays), recursing into the elements
        if they contain matobjects.
        '''
        elem_list = []
        for sub_elem in ndarray:
            if isinstance(sub_elem, spio.matlab.mio5_params.mat_struct):
                elem_list.append(_todict(sub_elem))
            elif isinstance(sub_elem, np.ndarray):
                elem_list.append(_tolist(sub_elem))
            else:
                elem_list.append(sub_elem)
        return elem_list

    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True, 
                                  verify_compressed_data_integrity=False)
    data = _check_keys(data)
    keys = [key for key in data.keys() if key[0] != '_']
    if len(keys) == 1:
        data = data[keys[0]]['Data']
    return data


def unwrap_voltage_to_charge(vs, wrap_voltage, dedv):
    unwrap_vs = vs[:]
    wrap_value = 0
    for i in range(len(vs)):
        unwrap_vs[i] = vs[i] + wrap_value / dedv
        if vs[i] > wrap_voltage:
            wrap_value = wrap_value + 1
        elif vs[i] < -wrap_voltage:
            wrap_value = wrap_value - 1
    unwrap_es = unwrap_vs * dedv

    # take into account aliasing
    delta_q = unwrap_es[1:] - unwrap_es[:-1]
    delta_q = alias(delta_q, 0.5)
    first_point = alias(unwrap_es[0], 0.5)
    return np.append(first_point, first_point + np.cumsum(delta_q))


def alias(x0, bound):
    """Takes a value or array x0 and aliases it into the range (-bound,bound)"""
    return np.mod(bound + x0, 2 * bound) - bound


def window_averaging(psd):
    psd_avg = np.zeros(psd.size)
    for i in range(psd.size):
        filter_fl = int(np.max([0, np.round(i - i / 4.)]))
        filter_fh = int(np.min([psd.size, np.round(i + i / 4.)]))
        psd_avg[i] = np.mean(psd[filter_fl:filter_fh])
    return psd_avg


def partition_data(data):
    """Partition data into even and odd for CPSD on single data set"""
    return data[0::2], data[1::2]


def crosspsd(z1, z2, fs):
    fft1 = np.fft.fft(z1)
    fft2 = np.fft.fft(z2)
    psd = (2. / fs / z1.size) * fft1 * np.conj(fft2)
    psd[0] = psd[0] / 2
    return psd


def partition_and_avg_psd(data, fs):
    """Splits the data into odd and even segments and takes the self cpsd.
    data is input as a matrix, where each line (the second axis) will
    give a cpsd that is then averaged along the first axis."""
    n, N = data.shape
    cpsd_avg = np.zeros(N / 4 + 1)
    for i in range(n):
        seg_even, seg_odd = partition_data(data[i])
        L = min(len(seg_even), len(seg_odd))
        seg_cpsd = crosspsd(seg_even[:L], seg_odd[:L], fs / 2)
        cpsd_avg = cpsd_avg + seg_cpsd[:N / 4 + 1]
    cpsd_avg = cpsd_avg / n
    cpsd_freq = np.arange(0, fs / 4. + 0.0001, fs / N)
    return cpsd_avg, cpsd_freq


def crosscorrelate(x, y, n):
    """Finds the cross correlation between the ones in two binary strings x and
    y of the same length.  n is the number of lags.  ccf is the cross
    correlated result, and lags is the lags=-n:n.  normalized to number of
    ones in x."""

    def ccorr(x, y, n):
        acf = np.zeros(n + 1)
        for i in range(n):
            if np.sum(x[:-1]) == 0:
                acf[i] = 0
            else:
                # we divide by the sum (not length) so that it is normalized to the
                # number of 1's and we get 100% corr if the 1s are correlated
                xx = x[:-i] if i > 0 else x
                acf[i] = np.sum(xx * y[i:]) / np.sum(xx)
        return acf

    ccf = np.zeros(2 * n + 1)
    ccf[n:] = ccorr(x, y, n)
    ccf[:n + 1] = np.flipud(ccorr(y, x, n))
    lags = np.arange(-n, n)
    return ccf, lags


def movingmean(a, n, axis=-1):
    """Should be equivilent to matlabs movmean"""

    def _movingmean(a, n):
        a_new = np.zeros(len(a), dtype=np.float)
        for i in range(len(a)):
            i_0 = max(0, i - int(np.floor(n / 2.)))
            i_end = min(len(a), i + int(np.ceil(n / 2.)))
            a_new[i] = np.mean(a[i_0:i_end])
        return a_new

    return np.apply_along_axis(_movingmean, axis, a, n)


def bin_dwell_times(o):
     last_val = o[0]
     count = 0
     x0 = []
     x1 = []
     for v in o:
         if v == last_val:
             count += 1
         else:
             if last_val == 0:
                 x0 += [count]
             else:
                 x1 += [count]
             count = 1
         last_val = v
     return x0,x1


def path_to_num(path):
    beginning, end = path.rsplit('_', 1)
    num, ext = end.split('.', 1)
    path = beginning + '_{:03d}.' + ext
    return path, int(num)
    
    
def matpaths(Q, date, fileName, fileNums, fluxNoise='fluxNoise2', **kwargs):
    if type(fileNums) == str:
        fileNums = kwargs[fileNums]
    if not isinstance(fileNums, (list, np.array)):
        fileNums = [[fileNums]]
    if not (type(fileNums[0]) in (list, np.array)):
        fileNums = [fileNums]
    if not isinstance(date, list):
        date = [date]
    if len(date) != len(fileNums):
        print len(date), len(fileNums)
        raise Exception('The number of dates must match the number of lists of file numbers')
    dataroot = os.getenv('DATA_ROOT').replace('\\','/')
    file_paths = []
    for j,d in enumerate(date):
        file_paths += [dataroot + '/{}/DR1 - 2019-12-17/CorrFar/{}/General/{}/'
                       '{}/MATLABData/{}_{:03d}.mat'.format(fluxNoise, Q, d, 
                                                          fileName, fileName, i)
                        for i in fileNums[j]]
    return file_paths


def apply_infidelity_correction(o, n_bins=9, thresh=0.5):
    o = o.astype(np.float)
    for trial in o:
        # trial[...] = movingmean(trial, n_bins) > thresh
        a = trial.astype(np.float)
        for i in range(trial.size):
            i_0 = max(0, i - int(np.floor(n_bins / 2.)))
            i_end = min(trial.size, i + int(np.ceil(n_bins / 2.)))
            trial[i] = np.mean(trial[i_0:i_end]) > thresh
    return o


# def apply_infidelity_correction_HMM(o, fidelity=[0.95, 0.8]):
#     o = o.astype(np.float)
#     for trial in o:
#         observed_signal = list(trial)   # ndarray -> list
#         observed_signal = [int(x) for x in observed_signal] # float -> int
#         recovered_signal = observed_to_recovered_signal(observed_signal, readout_fidelity=fidelity)
#         recovered_signal = list(np.float_(recovered_signal)) # int-> float
#         recovered_signal = np.asarray(recovered_signal)    # list-> ndarray
#         for i in range(trial.size):
#             trial[i] = recovered_signal[i]
#     return o
