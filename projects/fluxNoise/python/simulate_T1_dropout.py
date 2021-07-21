import numpy as np
import matplotlib.pyplot as plt
import noiselib
import importlib
importlib.reload(noiselib)
from noiselib import movingmean

T1_0 = 40   # max T1
tau = 5000  # QP relaxation time
N = 3000    # number of reps in a trial
n = 70      # number of trials to average
t_m = 10    # time between X gate and measurement
T_rep = 500 # Period of one rep
dtrigbins = 30 # uncertainty in trigger timing

def T1_decay(dt):
    """dt should always be negative"""
    return 1.-np.exp(-T_rep*(np.arange(0,N)-dt)/tau)

# T1 = T1_0 * np.ones(N)
# T1[1000:] = T1[1000:] * (1.-np.exp(-T_rep*np.linspace(0,1999,2000)/tau))

T1 = T1_0 * np.ones((n,N))
delta_t = np.random.random_sample(n) * dtrigbins
for i,dt in enumerate(delta_t):
    b = int(np.ceil(dt))
    T1[i,1000+b:] = T1[i,1000+b:] * T1_decay(dt - b)[:N-1000-b]

r = np.random.random_sample((n,N))
m = r < np.exp(-t_m/T1)

r = np.random.random_sample((n,N))
m[r < 0.05] = np.logical_not(m[r < 0.05])

plt.figure()
# plt.plot(T1.T)
plt.plot(np.mean(movingmean(m,30),axis=0))
plt.draw()
plt.pause(0.05)