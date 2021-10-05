import numpy as np
import matplotlib.pyplot as plt
from noiselib import movingmean
from scipy.optimize import curve_fit
import sys

T1_0 = 40      # max T1
tau = 5000     # QP relaxation time
N = 3000       # number of reps in a trial
n = 70         # number of trials to average
t_m = 10       # time between X gate and measurement
T_rep = 500    # Period of one rep
dtrigbins = 30 # uncertainty in trigger timing

def T1_decay(dt):

    # dt should always be negative
    return 1.-np.exp(-T_rep*(np.arange(0,N)-dt)/tau)

def fit_exp_dropoff(y, dl, o):

    # dl: domain length
    # o: origin

    f = max(y) - y

    f = f[o:(o + dl)]

    x = np.arange(0, dl)
    [a, b], _ = curve_fit(lambda t,a,b: a*np.exp(-b*t),  x,  f)

    fit_x = np.zeros(N)
    fit_x[0:o], fit_x[(o + dl):N] = np.nan, np.nan
    fit_x[o:(o + dl)] = np.arange(0, dl)

    fit_y = np.zeros(N)
    fit_y[0:o], fit_y[(o + dl):N] = np.nan, np.nan
    fit_y[o:(o + dl)] = -a*np.exp(-b*fit_x[o:(o + dl)])

    return fit_x, fit_y

def gen_data(T1_dropoff=True):

    T1 = T1_0 * np.ones((n,N))

    delta_t = (2*dtrigbins) * np.random.random_sample(n) - dtrigbins

    if T1_dropoff:
        for i,dt in enumerate(delta_t):
            b = int(np.ceil(dt))
            T1[i,1000+b:] = T1[i,1000+b:] * T1_decay(dt - b)[:N-1000-b]

    r = np.random.random_sample((n,N))
    m = r < np.exp(-t_m/T1)

    r = np.random.random_sample((n,N))
    m[r < 0.05] = np.logical_not(m[r < 0.05])

    return m, delta_t

def gen_shifts(shifts):

    n_to_shift = 10

    # get rid of repeated indices
    rows_to_shift = np.unique(np.random.randint(0, n-1, n_to_shift))
    n_to_shift = len(rows_to_shift)

    shifts[rows_to_shift] = np.random.randint(-30, 30, n_to_shift)

    return shifts

def mutate_m(m_old, shifts):

    m_new = np.zeros((n,N))
    for r in range(n):
        m_new[r] = shift_and_pad_array(m_old[r], int(shifts[r]))
        #m_new[r] = np.roll(m_old[r], int(shifts[r]))

    return m_new

def natural_selection(m, shifts, shifts_best, fitness_best):

    for p in range(n_offspring):
        y = np.sum(m[p], axis=0)

        # invert y axis so larger sum -> higher fitness
        f = max(y) - y

        # sum +/- some amount around T1 @ 1000
        fitness = np.sum(f[(1000 - 4):(1000 + 4)])

        # store new bests if mutation increases fitness
        if fitness > fitness_best:
            shifts_best = shifts[p]
            fitness_best = fitness

        return shifts_best, fitness_best

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

# fix seed for data generation
np.random.seed(3742373417)
m_orig, delta_t = gen_data()
#m_orig, delta_t = gen_data(T1_dropoff=False)

# reseed for everything else
np.random.seed()

n_gen = 100
n_offspring = 50 # per generation
shifts_best = np.zeros(n)
fitness_best = np.zeros(n_gen + 1)
shifts = np.zeros((n_offspring, n))
m = np.zeros((n_offspring, n, N))

# run fitness test on unshifted data
y = np.sum(m_orig, axis=0)
f = max(y) - y
fitness_best[0] = np.sum(f[(1000 - 4):(1000 + 4)])

for i in range(n_gen):
    print '{:03.1f}%'.format(float(i + 1)/float(n_gen) * 100)
    for j in range(n_offspring):

        shifts[j] = gen_shifts(shifts_best)
        m[j] = mutate_m(m_orig, shifts[j])

    shifts_best, fitness_best[i+1] = natural_selection(m, shifts, shifts_best, fitness_best[i])
    sys.stdout.write('\x1b[1A')
    sys.stdout.write('\x1b[2K')

x = np.arange(0, N)
m_best = mutate_m(m_orig, shifts_best)
y = np.mean(movingmean(m_best, 1), axis=0)

f, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [2, 1, 1]})
a0.plot(x, np.mean(movingmean(m_orig, 1), axis=0))
a0.plot(x, y)
a0.set_xlabel(r'Time $(\mu s)$')
a1.plot(delta_t)
a1.plot(shifts_best)
a1.set_xlabel('Time series index')
a1.set_ylabel(r'Shift $(\mu s)$')
a2.plot(fitness_best/max(fitness_best))
a2.set_xlabel('Generation')
a2.set_ylabel('Fitness')

#fit_domain_length = 500
#origin = 1000
#fit_x, fit_y = fit_exp_dropoff(y, fit_domain_length, origin)
#ax0.plot(origin + fit_x, fit_y + max(y[origin:(origin + fit_domain_length)]) )

#plt.draw()
#plt.pause(0.05)
plt.tight_layout()
plt.show()

