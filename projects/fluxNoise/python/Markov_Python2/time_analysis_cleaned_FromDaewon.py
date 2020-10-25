from __future__ import print_function, division
import pomegranate as pome
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import json
import jsonpickle


def generate_hidden_signal(length=8192, charge_burst_time=1000, p_QP=[0.01, 0.01], seed=0):
    np.random.seed(seed)

    p_QP_before = p_QP[0]
    p_QP_after = p_QP[1]
    hidden_signal = [0] * length

    hidden_signal[0] = np.random.choice([0, 1], p=[0.5, 0.5])

    for i in range(1, charge_burst_time):
        hidden_signal[i] = np.random.choice([hidden_signal[i - 1], 1 - hidden_signal[i - 1]],
                                            p=[1 - p_QP_before, p_QP_before])
    for i in range(charge_burst_time, length):
        hidden_signal[i] = np.random.choice([hidden_signal[i - 1], 1 - hidden_signal[i - 1]],
                                            p=[1 - p_QP_after, p_QP_after])

    return hidden_signal




def hidden_to_observed_signal(hidden_signal, readout_fidelity=[0.95, 0.75], seed=0):
    np.random.seed(seed+1)
    observed_signal = [0] * len(hidden_signal)

    # Readout Fidelity
    p_0g = readout_fidelity[0]
    p_1g = 1 - p_0g
    p_1e = readout_fidelity[1]
    p_0e = 1 - p_1e

    for i in range(len(hidden_signal)):
        if hidden_signal[i] == 0: observed_signal[i] = np.random.choice(['0', '1'], p=[p_0g, p_1g])
        if hidden_signal[i] == 1: observed_signal[i] = np.random.choice(['0', '1'], p=[p_0e, p_1e])

    return observed_signal



# here, BW estimates transition / emission probability
def observed_to_recovered_signal_BW(observed_signal, seed=0):
    seq = observed_signal
    flat_seq = seq[0]

    np.random.seed(seed+2)

    # random initial emission probability
    temp1 = np.random.uniform(low=0.5, high=1)  # indicates P[g|g]
    temp2 = np.random.uniform(low=0.5, high=1)  # indicates P[e|e]
    dists = [pome.DiscreteDistribution({'0': temp1, '1': 1-temp1}), pome.DiscreteDistribution({'0': 1-temp2, '1': temp2})]
    
    
    # random initial guess for transition matrix, we can polish if we know rough ranges for this
    trans_mat_guess = np.array(np.random.rand(2,2))

    # estimate values / hidden sequence by BW
    starts = np.array([0.5, 0.5])
    ends = np.array([0.5, 0.5])
    model = pome.HiddenMarkovModel.from_matrix(trans_mat_guess, dists, starts, ends)
    print('len(seq[0])=', asarray(seq).shape)
    model.fit(seq, algorithm='baum-welch', verbose=False)
    transition_count_matrix, emissions_matrix = model.forward_backward(flat_seq)
    est_tran_mat = model.dense_transition_matrix()
    recovered_seq = model.predict(flat_seq)
    
    # obtain estimated values...
    traneg = est_tran_mat[0][1]
    trange = est_tran_mat[1][0]

    empJSON0 = jsonpickle.encode(model.states[0])
    state_dicts0g = json.loads(empJSON0)
    empJSON1 = jsonpickle.encode(model.states[1])
    state_dicts1e = json.loads(empJSON1)

    emis0g = state_dicts0g['py/reduce'][1]['py/tuple'][0]['py/reduce'][1]['py/tuple'][0]['0']
    emis1e = state_dicts1e['py/reduce'][1]['py/tuple'][0]['py/reduce'][1]['py/tuple'][0]['1']

    return np.asarray(recovered_seq), emis0g, emis1e, traneg, trange


# observation length
var_length = 10000

# True transition probability
t_meas = 0.03  # units ms
t_QP = 10.0  # units ms
p_QP = 1.0 - np.exp(-t_meas/t_QP)

# True emission probability
p_0g = 0.95
p_1e = 0.85
readout_fidelity = [p_0g, p_1e]

print('Trans prob =%.5f'%p_QP, 'read_fidelity =', readout_fidelity, 'at length =', var_length)




# control the randomness
random_seed = 0

Observed_Signal_List=[]
Hidden_Signal = generate_hidden_signal(length = var_length, p_QP=[p_QP, p_QP], seed=random_seed)
for i in range(1):
    Observed_Signal = hidden_to_observed_signal(Hidden_Signal, readout_fidelity = readout_fidelity, seed=random_seed)
    Observed_Signal_List.append(Observed_Signal)

# hidden state reconstruction by BW, it also returns estimated transition / emission probability
Recovered_Signal_BW_List = list(observed_to_recovered_signal_BW(Observed_Signal_List, seed=random_seed))
Recovered_Signal_BW = Recovered_Signal_BW_List[:-4]

print('len(BW)=', len(Recovered_Signal_BW_List[0]))


# calculate and print error
emis0g_est = np.asarray(Recovered_Signal_BW_List[-4:-3])
emis1e_est = np.asarray(Recovered_Signal_BW_List[-3:-2])
traneg_est = np.asarray(np.asarray(Recovered_Signal_BW_List[-2:-1]))
trange_est = np.asarray(np.asarray(Recovered_Signal_BW_List[-1]))

error = np.average(np.absolute(np.array(Hidden_Signal) - (np.array(Recovered_Signal_BW))   ))

print('True P[e|g]=%.5f'%p_QP, 'P[g|e]=%.5f'%p_QP)
print('Estimated P[e|g]=%.5f'%traneg_est, 'P[g|e]=%.5f'%trange_est)

print('True P[0|g]=%.5f'%p_0g, 'P[1|e]=%.5f'%p_1e)
print('Estimated P[e|g]=%.5f'%emis0g_est, 'P[g|e]=%.5f'%emis1e_est)



fig = plt.figure(figsize=(12, 4))
plt.plot(asarray(Hidden_Signal)+1.5, 'o-', label=r"Hidden Transitions")
plt.plot(asarray(Observed_Signal), 'o-', label=r"Observed Transitions")
plt.plot(asarray(Recovered_Signal_BW_List[0])-1.5, 'o-', label=r"Recovered Transitions")
plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
plt.show()