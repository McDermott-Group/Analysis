"""
This is to extract QP tunneling rate from limited data points without high readout fidelities.
Date: VBT_lenVBT_lenA April
Author: Chuanhong Liu
"""

import numpy as np
import matplotlib.pyplot as plt

# Simulation data
VBT_len = 100
n_VBT_len = 10 # multiple of VBT_len
t = np.arange(0, VBT_len * n_VBT_len, 1)
state_trace = np.zeros(len(t))
meas_trace = np.zeros(len(t))
recover_trace = np.zeros(len(t))


t_QP_ = 10.0    # units ms
t_meas_ = 0.1     # units ms
p_QP_ = 1.0 - np.exp(-t_meas_/t_QP_)
# Transition matrix
p_ee_ = 1.0 - p_QP_
p_eg_ = p_QP_
p_ge_ = p_QP_
p_gg_ = 1.0 - p_QP_

# Initial Probabilities, also known as the prior probability, calculated from Transition matrix
p_e_ = 0.5
p_g_ = 1.0-p_e_

# Emission Probabilities, or measurement fidelities
p_e1_ = 0.75
p_e0_ = 1.0 - p_e1_
p_g0_ = 0.95
p_g1_ = 1.0 - p_g0_
readout_fidelity = [p_g0_, p_e1_]
state_trace[0] = 1
for i in range(1, len(state_trace)):
    state_trace[i] = np.random.choice([state_trace[i - 1], 1 - state_trace[i - 1]], p=[1 - p_QP_, p_QP_])
for i in range(len(meas_trace)):
    if state_trace[i] == 0:
        meas_trace[i] = np.random.choice([0, 1], p=[p_g0_, p_g1_])
    if state_trace[i] == 1:
        meas_trace[i] = np.random.choice([0, 1], p=[p_e0_, p_e1_])

def meas_trace_recover(meas_trace, readout_fidelity, t_QP):
    """

    :param meas_trace:
        A list of integer 0s ands 1s
    :param readout_fidelity:
        A list of ground and excited state readout fidelity, [f_g, f_e]
    :param t_QP:
        A number which is the assumed QP tunneling rate in units of ms
    :return: recover_trace, type array
        A list of integer 0s ands 1s based Markov Chain and smoothed Viberti algorithm
    Examples
    meas_trace = [0,1,0,1,1,1,1,0,0,1]
    readout_fidelity = [0.95, 0.75]
    t_QP = 10
    >>> meas_trace_recover(meas_trace, readout_fidelity, t_QP)
    [0,1,1,1,1,1,1,1,1,1]

    """

    recover_trace = np.zeros(len(meas_trace))
    # Get the optimized parameters from the readout_fidelity
    # viberti_param = {[0.95, 0.75]: [0.6, 0.75]}

    t_QP_VBT = 10  # _VBT means this variable is for Viberti algorithm
    # t_QP_VBT = t_QP   # _VBT means this variable is for Viberti algorithm
    t_meas_VBT = 0.1
    p_QP_VBT = 1.0 - np.exp(-t_meas_VBT/t_QP_VBT)
    # Transition matrix
    p_ee_VBT = (1.0 - p_QP_VBT)*5
    p_eg_VBT = (p_QP_VBT)*5
    p_ge_VBT = (p_QP_VBT)*5
    p_gg_VBT = (1.0 - p_QP_VBT)*5

    # Initial Probabilities, also known as the prior probability, calculated from Transition matrix
    p_e_VBT = 0.5
    p_g_VBT = 1.0-p_e_VBT

    # Emission Probabilities, this is the human choice and needs to be optimized, the p_g0_VBT is the most important one
    p_e1_VBT = 0.75*5
    p_e0_VBT = (1.0 - p_e1_VBT)*5
    p_g0_VBT = 0.95*(5)
    p_g1_VBT = (1.0 - p_g0_VBT)*5

    n_VBT_len = len(meas_trace)/VBT_len
    for chunck in range(n_VBT_len):
        probabilities = []
        if meas_trace[VBT_len * chunck] == 0:
            probabilities.append((p_e_VBT * p_e1_VBT, p_g_VBT * p_g1_VBT))
        else:
            probabilities.append((p_e_VBT * p_e0_VBT, p_g_VBT * p_g0_VBT))

        for i in range(VBT_len * chunck + 1, VBT_len * chunck + VBT_len):
            prev_e, prev_g = probabilities[-1]
            if meas_trace[i] == 0:
                curr_e = max(prev_e * p_ee_VBT * p_e1_VBT, prev_g * p_ge_VBT * p_e1_VBT)
                curr_g = max(prev_e * p_eg_VBT * p_g1_VBT, prev_g * p_gg_VBT * p_g1_VBT)
                probabilities.append((curr_e, curr_g))
            else:
                curr_e = max(prev_e * p_ee_VBT * p_e0_VBT, prev_g * p_ge_VBT * p_e0_VBT)
                curr_g = max(prev_e * p_eg_VBT * p_g0_VBT, prev_g * p_gg_VBT * p_g0_VBT)
                probabilities.append((curr_e, curr_g))

        for p_index in range(len(probabilities)):
            p = probabilities[p_index]
            if p[0] > p[1]:
                recover_trace[VBT_len * chunck + p_index] = 0
            else:
                recover_trace[VBT_len * chunck + p_index] = 1

    return recover_trace

def state_recover_diff(state_trace, recover_trace):
    """

    :param state_trace: The hidden state behind the measurement, [0,1,1,1...]
    :param recover_trace: The recovered state, [0, 1, 1, 1]
    :return: miscounts per VBT_len points (VBT_len ms)
    """
    state_transitions = 0
    recover_transitions = 0
    for i in range(len(state_trace)-1):
        if state_trace[i] != state_trace[i+1]:
            state_transitions += 1
        if recover_trace[i] != recover_trace[i+1]:
            recover_transitions += 1
    print('state:', state_transitions, 'recover:', recover_transitions)
    return state_transitions, recover_transitions
recover_trace = meas_trace_recover(meas_trace, readout_fidelity, t_QP_)
state_recover_diff(state_trace, recover_trace)


# plot the result
fig = plt.figure(figsize=(VBT_len, 10))
plt.plot(t, state_trace + 1.5, 'o-', label=r"State Trace")
plt.plot(t, meas_trace, 'o-', label=r"Meas Trace")
plt.plot(t, recover_trace - 1.5, 'o-', label=r"Recover Trace")
plt.plot(t, recover_trace-state_trace - 3, 'o-', label=r"Recover-State Trace")
plt.legend(bbox_to_anchor=(0.85, 0.75), loc=2)
plt.show()

#TO DO
# 1. Map readout fidelity to viberti fidelity
# 2. Non-VBT_len multiple length
# 3. Find a good metric to compare the state_transitions and recover_transitions
# 4. Target: find the suitable parameters to make the error < 10%, currently is around 15-25% for close readout, but <10% for 0.65e, 0.5g



