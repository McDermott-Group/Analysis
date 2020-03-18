import scipy.io as spio
import numpy as np


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
        
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    data = _check_keys(data)
    keys = [key for key in data.keys() if key[0] != '_']
    if len(keys) == 1:
        data = data[keys[0]]['Data']
    return data

def unwrap_voltage_to_charge(vs, wrap_voltage, dedv):
    unwrap_vs = vs[:]
    wrap_value = 0
    for i in range(len(vs)):
        unwrap_vs[i] = vs[i] + wrap_value/dedv
        if vs[i] > wrap_voltage:
            wrap_value = wrap_value + 1
        elif vs[i] < -wrap_voltage:
            wrap_value = wrap_value - 1
    unwrap_es = unwrap_vs * dedv
    
    # take into account aliasing
    delta_q = unwrap_es[1:] - unwrap_es[:-1]
    delta_q = alias(delta_q, 0.5)
    first_point = alias(unwrap_es[0], 0.5)
    return [first_point] + first_point+np.cumsum(delta_q)

def alias(x0, bound):
    """Takes a value or array x0 and aliases it into the range (-bound,bound)"""
    return np.mod(bound+x0, 2*bound) - bound

def window_averaging(psd):
    psd_avg = np.zeros(psd.size)
    for i in range(psd.size):
        filter_fl = int(np.max( [0, np.round(i-i/4.)] ))
        filter_fh = int(np.min( [psd.size, np.round(i+i/4.)] ))
        psd_avg[i] = np.mean( psd[filter_fl:filter_fh] )
    return psd_avg

def partition_data(data):
    """Partition data into even and odd for CPSD on single data set"""
    return data[0::2], data[1::2]

def crosspsd(z1, z2, fs):
    fft1 = np.fft.fft(z1)
    fft2 = np.fft.fft(z2)
    psd = (2./fs/z1.size)*fft1*np.conj(fft2)
    psd[0] = psd[0]/2
    return psd

def partition_and_avg_psd(data, fs):
    """Splits the data into odd and even segments and takes the self cpsd.
    data is input as a matrix, where each line (the second axis) will
    give a cpsd that is then averaged along the first axis."""
    n, N = data.shape
    cpsd_avg = np.zeros(N/4+1)
    for i in range(n):
        seg_even, seg_odd = partition_data(data[i])
        seg_cpsd = crosspsd(seg_even, seg_odd, fs/2)
        cpsd_avg = cpsd_avg + seg_cpsd[:N/4+1]
    cpsd_avg = cpsd_avg/n
    cpsd_freq = np.arange(0, fs/4.+0.0001, fs/N)
    return cpsd_avg, cpsd_freq