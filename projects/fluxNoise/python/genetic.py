import numpy as np
import matplotlib.pyplot as plt
import noiselib
reload(noiselib)
from noiselib import movingmean

def T1_decay(dt):
    """dt should always be negative"""
    return 1.-np.exp(-T_rep*(np.arange(0,N)-dt)/tau)

def gen_data():
    T1 = T1_0 * np.ones((n,N))
    delta_t = np.random.random_sample(n) * dtrigbins

    for i,dt in enumerate(delta_t):
        b = int(np.ceil(dt))
        T1[i,1000+b:] = T1[i,1000+b:] * T1_decay(dt - b)[:N-1000-b]

    r = np.random.random_sample((n,N))
    m = r < np.exp(-t_m/T1)

    r = np.random.random_sample((n,N))
    m[r < 0.05] = np.logical_not(m[r < 0.05])

    return m

def gen_noise():
    m = np.zeros((n,N))
    mu = 0
    sigma = 0.1

    for r in range(n):
        m[r] = mu + np.random.normal(mu, sigma, N)

    return m

def get_fitness(m, shifts, fitness_old):
    y = np.mean(movingmean(m,30), axis=0)
    idx_min = np.argmin(y)
    #bound_l = int(idx_min + shifts.min())
    #bound_u = int(idx_min + shifts.max())
    bound_l = int(1000 + shifts.min())
    bound_u = int(1000 + shifts.max())

    #fitness = np.sum(y[bound_l:bound_u])
    f = max(y) - y
    fitness = np.trapz(f[bound_l:bound_u])
    #fitness = abs(min(y) - np.median(y))

    if fitness <= fitness_old:
        more_fit = False
        fitness = fitness_old
    else:
        more_fit = True

    return y, fitness, more_fit

def shift_and_pad_array(arr, shift):
    ret = np.empty_like(arr)
    if shift > 0:
        ret[:shift] = np.nan
        ret[shift:] = arr[:-shift]
    elif shift < 0:
        ret[shift:] = np.nan
        ret[:shift] = arr[-shift:]
    else:
        ret[:] = arr

    return ret

def shift_rand_rows(m):
    n_to_shift = 10

    # get rid of repeated indices
    rows_to_shift = np.unique(np.random.randint(0, n-1, n_to_shift))
    n_to_shift = len(rows_to_shift)

    shifts = np.zeros(n)
    shifts[rows_to_shift] = np.random.randint(-30, 30, n_to_shift)

    for r in rows_to_shift:
        m[r] = shift_and_pad_array(m[r], int(shifts[r]))

    return m, shifts

T1_0 = 40      # max T1
tau = 5000     # QP relaxation time
N = 3000       # number of reps in a trial
n = 70         # number of trials to average
t_m = 10       # time between X gate and measurement
T_rep = 500    # Period of one rep
dtrigbins = 30 # uncertainty in trigger timing

fitness = 0
m = gen_data()
#m = gen_noise()
for i in range(10):
    print i
    m, shifts = shift_rand_rows(m)
    y_tmp, fitness, more_fit = get_fitness(m, shifts, fitness)
    if more_fit == True:
        y_best = y_tmp

x = np.arange(0, N)

plt.figure()
plt.plot(x, y_best)
#plt.plot(x, np.median(y) * np.ones(N))

plt.draw()
plt.pause(0.05)

