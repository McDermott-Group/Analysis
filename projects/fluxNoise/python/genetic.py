import numpy as np
import matplotlib.pyplot as plt
import pickle
from noiselib import movingmean
from scipy.optimize import curve_fit

def T1_decay(dt):
    """dt should always be negative"""
    return 1.-np.exp(-T_rep*(np.arange(0,N)-dt)/tau)

def fit_exp_dropoff(y, domain_length, origin):
    f = np.max(y) - y

    f = f[origin:(origin+domain_length)]

    x = np.arange(0, domain_length)
    popt, _ = curve_fit(lambda t,a,b: a*np.exp(-b*t),  x,  f)
    a, b = popt[0], popt[1]

    fit_x = np.zeros(N)
    fit_x[0:origin], fit_x[(origin+domain_length):N] = np.nan, np.nan
    fit_x[origin:(origin+domain_length)] = np.arange(0,domain_length)

    fit_y = np.zeros(N)
    fit_y[0:origin], fit_y[(origin+domain_length):N] = np.nan, np.nan
    fit_y[origin:(origin+domain_length)] = -a*np.exp(-b*fit_x[origin:(origin+domain_length)])

    return fit_x, fit_y

def gen_data(T1_dropoff=True):
    T1 = T1_0 * np.ones((n,N))

    delta_t = np.random.random_sample(n) * dtrigbins

    if T1_dropoff:
        for i,dt in enumerate(delta_t):
            b = int(np.ceil(dt))
            T1[i,1000+b:] = T1[i,1000+b:] * T1_decay(dt - b)[:N-1000-b]

    r = np.random.random_sample((n,N))
    m = r < np.exp(-t_m/T1)

    r = np.random.random_sample((n,N))
    m[r < 0.05] = np.logical_not(m[r < 0.05])

    return m

def natural_selection(m, m_best, fitness_best):
    y = np.mean(movingmean(m, 30), axis=0)

    # invert y axis so larger sum -> higher fitness
    f = max(y) - y

    # sum +/- some amount around T1 @ 1000
    fitness = np.sum(f[(1000 - 5):(1000 + 5)])
    #fitness = np.trapz(f[990:1010])
    #fitness = abs(min(y) - np.median(y))

    # return previous best if mutation does not increase fitness
    if fitness <= fitness_best:
        m = m_best
        fitness = fitness_best

    return m, fitness

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
        #m[r] = shift_and_pad_array(m[r], int(shifts[r]))
        m[r] = np.roll(m[r], int(shifts[r]))

    return m

T1_0 = 40      # max T1
tau = 5000     # QP relaxation time
N = 3000       # number of reps in a trial
n = 70         # number of trials to average
t_m = 10       # time between X gate and measurement
T_rep = 500    # Period of one rep
dtrigbins = 30 # uncertainty in trigger timing


#m_best = gen_data()
#m_best = gen_data(T1_dropoff=False)

#f = open('T1_data_simulated.pckl', 'wb')
#pickle.dump(m_best, f)
#f.close()

f = open('T1_data_simulated.pckl', 'rb')
m_best = pickle.load(f)
f.close()

fitness_best = 0
for i in range(0):
    print i
    m = shift_rand_rows(m_best)
    m_best, fitness_best = natural_selection(m, m_best, fitness_best)

x = np.arange(0, N)
y = np.mean(movingmean(m_best, 30), axis=0)

plt.figure()
plt.plot(x, y)

#fit_domain_length = 500
#origin = 1000
#fit_x, fit_y = fit_exp_dropoff(y, fit_domain_length, origin)
#plt.plot(origin + fit_x, fit_y + max(y[origin:(origin + fit_domain_length)]) )

plt.draw()
plt.pause(0.05)

