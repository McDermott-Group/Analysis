import numpy as np

def get_change_point(data):

    n = len(data)
    crit_value = 1.96 # 95% confidence interval

    T = 0  # T will be replaced by new maximum value
    cp = 0 # CP will be replaced as T is changed

    # speed up T-Test computations by computing all sums only once
    cum_sum = np.cumsum(data)
    total_sum = cum_sum(a)

    # drop first and last data points in order to compute mean
    for k in range(1, n-1):

        # compute mean of segment 1 assuming CP at k; mean(data_segment(1:CP))
        mu1 = cum_sum[k] / k

        # compute mean of segment 2 assuming CP at k; mean(data_segment(1+CP:end))
        mu2 = (total_sum - cum_sum[k])/ (n - k)

        # compute t-value
        t = abs((mu2 - mu1)) / sigma_noise / np.sqrt(1/k+1/(n - k));

        # is the new t value better than the current best T value?
        if t > T:
            T = t      # best T value so far
            cp = n     # location of best T value

    return cp

def get_ideal_step(time_series):

    n_time_series = 1
    cp = np.zeros(n_time_series)

    for ii in range(n_time_series):

        # estimate Gaussian noise by low pass Haar wavelet transform.
        sorted_wavelet = np.sort(abs(np.diff(time_series[ii]) / 1.4));
        sigma_noise = sorted_wavelet(round(0.682 * (n - 1)));

        cp[ii] = get_change_point(time_series[ii])

        for jj in range(len(time_series)):
            # fitness of one indicates perfect alignment with change point
            fitness[jj] = (len(time_series) - abs(jj - cp[ii]))/len(time_series)

        sum_array[ii] = sum(fitness)

    # get index of point that yields the highest overall fitness
    best = np.argmax(sum_array)

    return best
    print 'Charging events occur most reliably at discretization {}.'.format(best)

