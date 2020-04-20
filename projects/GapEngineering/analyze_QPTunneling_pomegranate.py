from pomegranate import *
from numpy import *
import string
import matplotlib.pyplot as plt

n = 5000
Time = arange(0, n, 1)
Hidden_Signal = [0]*n
Observed_Signal = [0]*n

p_e1 = 0.75
p_e0 = 0.25
p_g0 = 0.95
p_g1 = 0.05

t_QP_ = 10.0    # units ms
t_meas_ = 0.1     # units ms
p_QP_ = 1.0 - exp(-t_meas_/t_QP_)

Hidden_Signal[0] = 1
for i in range(1, len(Hidden_Signal)):
    Hidden_Signal[i] = random.choice([Hidden_Signal[i - 1], 1 - Hidden_Signal[i - 1]], p = [1 - p_QP_, p_QP_])

for i in range(n):
    if Hidden_Signal[i] == 0: Observed_Signal[i] = random.choice([0, 1], p = [p_g0, p_g1])
    if Hidden_Signal[i] == 1: Observed_Signal[i] = random.choice([0, 1], p = [p_e0 , p_e1])


Observed_Signal_copy = Observed_Signal
Observed_Signal = str(Observed_Signal).strip('[]')
Observed_Sequence = ''.join(Observed_Signal)
Observed_Sequence = Observed_Sequence.replace(',', '')
Observed_Sequence = Observed_Sequence.translate({ord(c): None for c in string.whitespace})

d1 = DiscreteDistribution({'0' : p_g0, '1' : p_g1})
d2 = DiscreteDistribution({'0' : p_e0, '1' : p_e1})

g = State(d1, name = "g")
e = State(d2, name = "e")

model = HiddenMarkovModel('example')
model.add_states([g, e])

model.add_transition(model.start, g, 0.5)
model.add_transition(model.start, e, 0.5)
model.add_transition(g, g, 0.99)
model.add_transition(g, e, 0.01)
model.add_transition(e, e, 0.99)
model.add_transition(e, g, 0.01)
model.add_transition(e, model.end, 0.5)
model.bake()


Recovered_Signal = ", ".join(state.name for i, state in model.viterbi(Observed_Sequence)[1])


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

fig = plt.figure(figsize=(12, 6))
plt.plot(Time, Hidden_Signal + 1.5, 'o-' )
plt.plot(Time, Observed_Signal_copy, 'o-')
plt.plot(Time, Recovered_Signal - 1.5, 'o-')
plt.show()
