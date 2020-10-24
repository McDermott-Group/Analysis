import numpy as np
import matplotlib.pyplot as plt

def gen_hidden(length=8192, charge_burst_time=4000, p_QP=[0.01, 0.01]):

    hid = np.zeros(length, dtype=int)

    # initial parity
    hid[0] = np.random.choice([0, 1], p=[0.5, 0.5])

    for i in range(1, charge_burst_time):
        hid[i] = np.random.choice(
                [hid[i - 1], 1 - hid[i - 1]],
                p=[1 - p_QP[0], p_QP[0]]
                )

    for i in range(charge_burst_time, length):
        hid[i] = np.random.choice(
                [hid[i - 1], 1 - hid[i - 1]],
                p=[1 - p_QP[1], p_QP[1]]
                )

    return hid




def hidden_to_observed(hidden_signal, readout_fidelity=[0.95, 0.75]):

    obs = np.zeros(len(hidden_signal), dtype=int)

    p_0g = readout_fidelity[0]
    p_1g = 1 - p_0g
    p_1e = readout_fidelity[1]
    p_0e = 1 - p_1e

    obs[hidden_signal == 0] = np.random.choice([0, 1],
            size=sum(hidden_signal == 0), p=[p_0g, p_1g]
            )

    obs[hidden_signal == 1] = np.random.choice([0, 1],
            size=sum(hidden_signal == 1), p=[p_0e, p_1e]
            )

    return obs




def states_to_observables(state_seq, emission_states):

    # map hidden states to recovered states
    # emission_states: list element corresponds to state of the same index
    #                  e.g. [0, 0, 1, 1] maps states 0 and 1 to observable
    #                  0, and states 2 and 3 to observable 1.

    obs_seq = np.zeros(len(state_seq), dtype=int)
    for i,s in enumerate(emission_states):
        obs_seq[state_seq == i] = s*np.ones(sum(state_seq == i))

    return obs_seq




def forward(V, a, b, init_dist):
    alpha = np.zeros((V.shape[0], a.shape[0]))
    alpha[0, :] = init_dist * b[:, V[0]]

    scale = np.ones(V.shape[0])

    t = 1
    while t < V.shape[0]:

        for j in range(a.shape[0]):
            alpha[t, j] = scale[t] * alpha[t - 1, :].dot(a[:, j]) * b[j, V[t]]

        # go back 1 step in loop and scale inversely with probabilities once
        # values start getting very small to avoid underflows.  these scalars
        # are accounted for when a and b are normalized.
        if (alpha[t, :] < 1e-20).any():
            scale[t-1] = 1/sum(alpha[t-1, :])
            t -= 1
        else:
            t += 1

    return alpha, scale




def backward(V, a, b, scale):
    beta = np.zeros((V.shape[0], a.shape[0]))

    # beta(T) = 1
    beta[V.shape[0] - 1] = np.ones((a.shape[0]))

    # loop backwards from T-2 to 0
    for t in range(V.shape[0] - 2, -1, -1):
        for j in range(a.shape[0]):
            # beta needs to be scaled _exactly_ as alpha is.
            beta[t, j] = scale[t] * (beta[t + 1, :] * b[:, V[t + 1]]).dot(a[j, :])

    return beta




def baum_welch(V, a, b, init_dist, conv_thres=1e-9, max_iter=100):

    M = a.shape[0]
    if V.ndim == 1:
        V = np.array([V])
    R, T = V.shape

    it = 1
    old_a = np.ones_like(a)
    old_b = np.ones_like(b)

    # until convergence defined by threshold
    # or max # iterations has been reached
    while ( (abs(np.log(a/old_a)) > conv_thres).all() \
            and ( abs(np.log(b/old_b)) > conv_thres).all() ) \
            and it < max_iter:

        xi = np.zeros((R, M, M, T - 1))
        alpha = np.zeros((R, T, M))
        beta = np.zeros((R, T, M))
        # xi[r, i, j, t] = P(X_t = i, X_t+1 = j | r, theta)
        # gamma[r, i, t] = P(X_t = i | r, theta)
        for r in range(R):
            alpha[r], scale = forward(V[r], a, b, init_dist)
            beta[r] = backward(V[r], a, b, scale)

            for t in range(T - 1):
                den = np.dot(
                        np.dot(alpha[r, t, :].T, a) * b[:, V[r, t + 1]].T,
                        beta[r, t + 1, :]
                        )
                for i in range(M):
                    num = alpha[r, t, i] * a[i, :] * \
                            b[:, V[r, t + 1]].T * beta[r, t + 1, :].T
                    xi[r, i, :, t] = num / den

        gamma = np.sum(xi, axis=2)

        # store for conv test
        old_a = a
        old_b = b

        a = np.sum(np.sum(xi, 0), 2) / \
                np.sum(np.sum(gamma, 0), 1).reshape((-1, 1))

        # append T'th element in gamma
        tmp = gamma
        gamma = np.zeros((R, M, T))
        for r in range(R):
            gamma[r, :, :] = np.hstack(
                    (tmp[r], np.sum(xi[r, :, :, T - 2], axis=0).reshape((-1, 1)))
                    )

        den = np.sum(np.sum(gamma, 0), 1)
        for l in range(b.shape[1]):
            g = 0
            for r in range(R):
                g += (V[r] == l) * gamma[r, :, :]
            b[:, l] = np.sum(g, 1)

        b = b / den.reshape((-1, 1))

        it += 1

    if it != max_iter:
        print("BW converged in", it, "iterations")
    else:
        print("BW did not converge in", max_iter, "iterations")

    return a, b




def viterbi(V, a, b, init_dist):
    T = V.shape[0]
    M = a.shape[0]

    omega = np.zeros((T, M))
    omega[0, :] = np.log(init_dist * b[:, V[0]])

    prev = np.zeros((T - 1, M))

    for t in range(1, T):
        for j in range(M):
            # same as fwd probability
            p = omega[t - 1] + np.log(a[:, j]) + np.log(b[j, V[t]])

            # most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(p)

            # probability of the most probable state (2)
            omega[t, j] = np.max(p)

    S = np.zeros(T)

    # most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])

    S[0] = last_state

    rev_idx = 1
    for i in range(T - 2, -1, -1):
        S[rev_idx] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        rev_idx += 1

    # path is flipped from backtracking
    S = np.flip(S, axis=0)

    return S




t_meas = 0.03 # ms
t_QP = 10.0   # ms
p_QP = 1.0 - np.exp(-t_meas/t_QP)

p_0g = 0.90
p_1e = 0.80

np.random.seed(1454146880) # fix for data generation (in [0, 2^32-1])
length = 8192
n_trials = 10

hid_arr = np.zeros((n_trials, length), dtype=int)
obs_arr = np.zeros((n_trials, length), dtype=int)
for ii in range(n_trials):
    hid_arr[ii, :] = gen_hidden(length, charge_burst_time=4000, p_QP=[p_QP, p_QP])
    obs_arr[ii, :] = hidden_to_observed(hid_arr[ii,:], readout_fidelity=[p_0g, p_1e])


np.random.seed()

# transitions, normalized s.t. each row sums to 1
a = np.random.random((2, 2))
a = a / np.sum(a, axis=1).reshape((-1, 1))

# emissions, normalized s.t. each row sums to 1
b = np.random.random((2, 2))
b = b / np.sum(b, axis=1).reshape((-1, 1))

init_dist = np.array((0.5, 0.5))
init_dist /= np.sum(init_dist)

a, b = baum_welch(obs_arr, a, b, init_dist, conv_thres=1e-9, max_iter=1000)

# adjust manually if ICs give flipped emission probabilities
# i.e. if BW converges toward 'opposite' local optimum
if b[0,1] > b[0,0]:
    b = np.flip(b, axis=0)

print("a = \n{}\nb = \n{}".format(a, b))


seq_arr = np.zeros((n_trials, length), dtype=int)
rec_arr = np.zeros((n_trials, length), dtype=int)
err_arr = np.zeros((n_trials, length), dtype=int)
n_err_arr = np.zeros(n_trials, dtype=int)
for ii in range(n_trials):
    seq_arr[ii,:] = viterbi(obs_arr[ii,:], a, b, init_dist)
    rec_arr[ii,:] = states_to_observables(seq_arr[ii,:], emission_states=[0, 1])
    err_arr[ii,:] = hid_arr[ii,:] - rec_arr[ii,:]
    n_err_arr[ii] = sum(err_arr[ii,:] != 0)


#trial = np.random.randint(0, n_trials)
f, ax = plt.subplots(4, 1)
for trial in range(n_trials):

    ax[0].plot(hid_arr[trial], 'o-')
    ax[0].set_title('Hidden')

    ax[1].plot(obs_arr[trial], 'o-')
    ax[1].set_title('Observed')

    ax[2].plot(rec_arr[trial], 'o-')
    ax[2].set_title('Recovered')

    ax[3].plot(err_arr[trial], 'o')
    ax[3].set_title('Error (Hidden - Recovered) (n = {})'.format(n_err_arr[trial]))

    f.tight_layout()
    f.savefig('hmm_{}.png'.format(trial))
    for ii in range(len(ax)):
        ax[ii].clear()
    #plt.show()

