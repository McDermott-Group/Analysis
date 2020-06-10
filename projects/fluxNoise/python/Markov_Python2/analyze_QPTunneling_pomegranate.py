"""Model Pomegranate Author: Jacob Schreiber.
Author Vincent Liu and Sohair Abdullah."""

import numpy as np
from pomegranate import *
from numpy import *
import string
import matplotlib.pyplot as plt


def generate_hidden_signal(length=8192, charge_burst_time=4000, p_QP=[0.01, 0.01]):
    """
    Generate hidden signal from given parameters
    :param length: the length of the hidden signal, eg. 8192
    :param charge_burst_time: the time where we have the charge burst, eg. 2456
    :param p_QP: p_QP = [p_QP_before, p_QP_after], a two element list consists QP tunneling rate before and after the burst
        p_QP = 1.0-np.exp(-t_meas/t_QP), e.g. t_meas = 100 us, t_QP = 10 ms, p_QP=0.01
    :return: the hidden signal [0,1,0,0,1,1,1,...]
    """

    p_QP_before = p_QP[0]
    p_QP_after = p_QP[1]
    hidden_signal = [0] * length

    hidden_signal[0] = random.choice([0, 1], p=[0.5, 0.5])  # Initial parity

    for i in range(1, charge_burst_time):
        hidden_signal[i] = random.choice([hidden_signal[i - 1], 1 - hidden_signal[i - 1]],
                                         p=[1 - p_QP_before, p_QP_before])
    for i in range(charge_burst_time, length):
        hidden_signal[i] = random.choice([hidden_signal[i - 1], 1 - hidden_signal[i - 1]],
                                         p=[1 - p_QP_after, p_QP_after])

    return hidden_signal


def hidden_to_observed_signal(hidden_signal, readout_fidelity=[0.95, 0.75]):
    """
    Convert the hidden signal to the observed signal
    :param hidden_signal: the hidden signal list [0,1,1,1,0,...]
    :param readout_fidelity: readout_fidelity = [p_0g, p_1e]
    :return: observed signal [0,1,1,1,1, ...]
    """

    observed_signal = [0] * len(hidden_signal)

    # Readout Fidelity
    p_0g = readout_fidelity[0]
    p_1g = 1 - p_0g
    p_1e = readout_fidelity[1]
    p_0e = 1 - p_1e

    for i in range(len(hidden_signal)):
        if hidden_signal[i] == 0: observed_signal[i] = random.choice([0, 1], p=[p_0g, p_1g])
        if hidden_signal[i] == 1: observed_signal[i] = random.choice([0, 1], p=[p_0e, p_1e])

    return observed_signal


def observed_to_recovered_signal(observed_signal, readout_fidelity=[0.95, 0.75], p_QP=0.01):
    """

    :param observed_signal: observed signal for parity from given readout fidelity
    :param readout_fidelity: the given readout fidelity will be used as emission probability for Viberti Algorithm
    :param p_QP: QP tunneling rate
    :return: recovered signal
    """
    p_0g_VTB = readout_fidelity[0]
    p_1g_VTB = 1 - p_0g_VTB
    p_1e_VTB = readout_fidelity[1]
    p_0e_VTB = 1 - p_1e_VTB

    p_QP_VTB = p_QP

    observed_signal = str(observed_signal).strip('[]')
    observed_sequence = ''.join(observed_signal)
    observed_sequence = observed_sequence.replace(',', '')
    observed_sequence = "".join(observed_sequence.split())
    d1 = DiscreteDistribution({'0': p_0g_VTB, '1': p_1g_VTB})
    d2 = DiscreteDistribution({'0': p_0e_VTB, '1': p_1e_VTB})

    g = State(d1, name="g")
    e = State(d2, name="e")

    model = HiddenMarkovModel('example')
    model.add_states([g, e])

    model.add_transition(model.start, g, 0.5)
    model.add_transition(model.start, e, 0.5)
    model.add_transition(g, g, 1 - p_QP_VTB)
    model.add_transition(g, e, p_QP_VTB)
    model.add_transition(e, e, 1 - p_QP_VTB)
    model.add_transition(e, g, p_QP_VTB)
    model.add_transition(e, model.end, 0.5)
    model.bake()

    recovered_signal = ", ".join(state.name for i, state in model.viterbi(observed_sequence)[1])

    # Data post-processing
    recovered_signal = recovered_signal.split(',')
    recovered_signal = recovered_signal[1:-1]
    for i in range(len(recovered_signal)):
        if recovered_signal[i] == ' e':
            recovered_signal[i] = 1
        if recovered_signal[i] == ' g':
            recovered_signal[i] = 0

    recovered_signal = asarray(recovered_signal)
    return recovered_signal


def observed_to_recovered_signal_LIU(observed_signal, readout_fidelity=[0.95, 0.75], p_QP=0.01):
    """

    :param observed_signal:
        A list of integer 0s ands 1s
    :param readout_fidelity:
        A list of ground and excited state readout fidelity, [f_g, f_e]
    :param t_QP:
        A number which is the assumed QP tunneling rate in units of ms
    :return: recovered_signal, type array
        A list of integer 0s ands 1s based Markov_Python3 Chain and smoothed Viberti algorithm
    Examples
    observed_signal = [0,1,0,1,1,1,1,0,0,1]
    readout_fidelity = [0.95, 0.75]
    t_QP = 10
    # >>> observed_to_recovered_signal_LIU(observed_signal, readout_fidelity, p_QP)
    [0,1,1,1,1,1,1,1,1,1]

    """

    recovered_signal = np.zeros(len(observed_signal))

    # Transition matrix
    p_ee_VTB = np.log(1.0 - p_QP)
    p_eg_VTB = np.log(p_QP)
    p_ge_VTB = np.log(p_QP)
    p_gg_VTB = np.log(1.0 - p_QP)

    # Initial Probabilities, also known as the prior probability, calculated from Transition matrix
    p_e_VTB = np.log(0.5)
    p_g_VTB = np.log(1.0-p_e_VTB)

    # Emission Probabilities, this is the human choice and needs to be optimized, the p_0g_VTB is the most important one
    p_0g_VTB = np.log(readout_fidelity[0])
    p_1g_VTB = np.log(1.0-readout_fidelity[0])
    p_1e_VTB = np.log(readout_fidelity[1])
    p_0e_VTB = np.log(1.0-readout_fidelity[1])

    print p_0g_VTB, p_1g_VTB, p_0e_VTB, p_1e_VTB

    probabilities = []
    if observed_signal[0] == 0:
        probabilities.append((p_e_VTB + p_1e_VTB, p_g_VTB + p_1g_VTB))
    else:
        probabilities.append((p_e_VTB + p_0e_VTB, p_g_VTB + p_0g_VTB))

    for i in range(len(observed_signal)):
        prev_e, prev_g = probabilities[-1]
        if observed_signal[i] == 0:
            curr_e = max(prev_e + p_ee_VTB + p_1e_VTB, prev_g + p_ge_VTB + p_1e_VTB)
            curr_g = max(prev_e + p_eg_VTB + p_1g_VTB, prev_g + p_gg_VTB + p_1g_VTB)
            probabilities.append((curr_e, curr_g))
        else:
            curr_e = max(prev_e + p_ee_VTB + p_0e_VTB, prev_g + p_ge_VTB + p_0e_VTB)
            curr_g = max(prev_e + p_eg_VTB + p_0g_VTB, prev_g + p_gg_VTB + p_0g_VTB)
            probabilities.append((curr_e, curr_g))

    for p_index in range(len(probabilities)-1):
        p = probabilities[p_index]
        if p[0] < p[1]:
            recovered_signal[p_index] = 0
        else:
            recovered_signal[p_index] = 1

    return recovered_signal


def transitions_count(signal):
    """

    :param signal: signal trace
    :return: the number of 0->1, 1->0 counts
    """
    transitions = 0
    for i in range(len(signal) - 1):
        if signal[i] != signal[i + 1]:
            transitions += 1
    return transitions


"""
Parameters Setup
"""

# For Generating Simulation Data
# t_meas = 0.1  # units ms
# t_QP_before = 10.0  # units ms
# t_QP_after = 10.0  # units ms
# p_QP_before = 1.0 - exp(-t_meas/t_QP_before)
# p_QP_after = 1.0 - exp(-t_meas/t_QP_after)
# p_QP = [p_QP_before, p_QP_after]
# p_0g = 0.95
# p_1e= 0.75
# readout_fidelity = [p_0g, p_1e]
#
# # For VTB
# t_QP_VTB = 10  # units ms
# p_QP_VTB = 1.0 - exp(-t_meas/t_QP_VTB)
# p_0g_VTB = 0.95
# p_1e_VTB = 0.75
# readout_fidelity_VTB = [p_0g_VTB, p_1e_VTB]
#
# Hidden_Signal = generate_hidden_signal(p_QP=[0.01, 0.01])
# Observed_Signal = hidden_to_observed_signal(Hidden_Signal)
# Recovered_Signal = observed_to_recovered_signal(Observed_Signal, p_QP=0.01)

# fig = plt.figure(figsize=(12, 4))
# plt.plot(asarray(Hidden_Signal)+1.5, 'o-', label=r"{} Hidden Transitions".format(transitions_count(Hidden_Signal)))
# plt.plot(asarray(Recovered_Signal), 'o-', label=r"{} Recovered Transitions".format(transitions_count(Recovered_Signal)))
# plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
# plt.show()

"""
TO DO
2. test its robustness against several fidelity and QP values
"""
