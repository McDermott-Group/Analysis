from pomegranate import *
from numpy import *
import string
import matplotlib.pyplot as plt

n = 100
Time = arange(0, n, 1)
Hidden_Signal = [0]*n
Observed_Signal = [0]*n

# Readout Fidelity
p_e1 = 0.75
p_e0 = 1-p_e1
p_g0 = 0.95
p_g1 = 1-p_g0

# Measurement setup
t_QP_before = 10.0     # units ms
t_QP_after = 10.0     # units ms
t_meas = 0.1    # units ms
p_QP_before = 1.0 - exp(-t_meas/t_QP_before)
p_QP_after = 1.0 - exp(-t_meas/t_QP_after)

# Numbers for pomegranate which are fixed
p_e1_VBT = p_e1 + 0.0
p_e0_VBT = 1-p_e1_VBT
p_g0_VBT = p_g0 + 0.0
p_g1_VBT = 1-p_g0_VBT

t_QP_VBT = 10    # units ms
p_QP_VBT = 1.0 - exp(-t_meas/t_QP_VBT)

# Generating simulation signal

Hidden_Signal[0] = 1
for i in range(1, n//2):
    Hidden_Signal[i] = random.choice([Hidden_Signal[i - 1], 1 - Hidden_Signal[i - 1]], p=[1 - p_QP_before, p_QP_before])
for i in range(n//2, n):
    Hidden_Signal[i] = random.choice([Hidden_Signal[i - 1], 1 - Hidden_Signal[i - 1]], p=[1 - p_QP_after, p_QP_after])

for i in range(n):
    if Hidden_Signal[i] == 0: Observed_Signal[i] = random.choice([0, 1], p=[p_g0, p_g1])
    if Hidden_Signal[i] == 1: Observed_Signal[i] = random.choice([0, 1], p=[p_e0, p_e1])

Observed_Signal_copy = Observed_Signal
Observed_Signal = str(Observed_Signal).strip('[]')
Observed_Sequence = ''.join(Observed_Signal)
Observed_Sequence = Observed_Sequence.replace(',', '')
Observed_Sequence = Observed_Sequence.translate({ord(c): None for c in string.whitespace})

d1 = DiscreteDistribution({'0': p_g0_VBT, '1': p_g1_VBT})
d2 = DiscreteDistribution({'0': p_e0_VBT, '1': p_e1_VBT})

g = State(d1, name ="g")
e = State(d2, name ="e")

model = HiddenMarkovModel('example')
model.add_states([g, e])

model.add_transition(model.start, g, 0.5)
model.add_transition(model.start, e, 0.5)
model.add_transition(g, g, 1-p_QP_VBT)
model.add_transition(g, e, p_QP_VBT)
model.add_transition(e, e, 1-p_QP_VBT)
model.add_transition(e, g, p_QP_VBT)
model.add_transition(e, model.end, 0.5)
model.bake()

Recovered_Signal = ", ".join(state.name for i, state in model.viterbi(Observed_Sequence)[1])

# Data post-processing
Recovered_Signal = Recovered_Signal.split(',')
Recovered_Signal = Recovered_Signal[1:-1]
for i in range(len(Recovered_Signal)):
    if(Recovered_Signal[i] == ' e'):
        Recovered_Signal[i] = 1
    if(Recovered_Signal[i] == ' g'):
        Recovered_Signal[i] = 0

Hidden_Signal = asarray(Hidden_Signal)
Observed_Signal = asarray(Observed_Signal)
Recovered_Signal = asarray(Recovered_Signal)

def transitions_count(signal):
    """

    :param signal: signal trace
    :return: the number of 0->1, 1->0 counts
    """
    transitions = 0
    for i in range(len(signal)-1):
        if signal[i] != signal[i+1]:
            transitions += 1
    return transitions

fig = plt.figure(figsize=(12, 6))
plt.plot(Time, Hidden_Signal + 1.5, 'o-', label=r"{} Hidden Transitions".format(transitions_count(Hidden_Signal)))
plt.plot(Time, Recovered_Signal, 'o-', label=r"{} Recovered Transitions".format(transitions_count(Recovered_Signal)))
# plt.plot(Time, Observed_Signal_copy - 1.5, 'o-', label=r"Observed")
plt.legend(bbox_to_anchor=(0.85, 0.55), loc=2)
plt.show()