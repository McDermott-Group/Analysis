"""Model Pomegranate Author: Jacob Schreiber.
Author Vincent Liu and Sohair Abdullah."""


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
# Hidden_Signal = generate_hidden_signal()
# Observed_Signal = hidden_to_observed_signal(Hidden_Signal)
# Recovered_Signal = observed_to_recovered_signal(Observed_Signal)
# #
# fig = plt.figure(figsize=(12, 4))
# plt.plot(asarray(Hidden_Signal)+1.5, 'o-', label=r"{} Hidden Transitions".format(transitions_count(Hidden_Signal)))
# plt.plot(asarray(Recovered_Signal), 'o-', label=r"{} Recovered Transitions".format(transitions_count(Recovered_Signal)))
# plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
# plt.show()

"""
TO DO
2. test its robustness against several fidelity and QP values
"""
