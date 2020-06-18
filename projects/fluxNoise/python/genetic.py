import numpy as np
import ruptures as rpt

def get_change_point(data, sigma_noise):

    n = len(data)
    crit_value = 1.96 # 95% confidence interval

    T = 0  # T will be replaced by new maximum value
    cp = 0 # CP will be replaced as T is changed

    cum_sum = np.cumsum(data)
    total_sum = sum(data)

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

def get_ideal_origin(all_time_series):

    # all_times_series: array w/ one time series per row

    n_time_series = len(all_time_series)
    cp = np.zeros(n_time_series)

    # allocate assuming all time series are of the same length
    n_points = len(all_time_series[0])
    fitness = np.zeros(n_time_series, n_points)
    fit_sum_arr = np.zeros(n_points)

    for ii in range(n_time_series):

        time_series = all_time_series[ii]

        # estimate Gaussian noise by low pass Haar wavelet transform
        #sorted_wavelet = np.sort(abs(np.diff(time_series) / 1.4))
        #sigma_noise = sorted_wavelet[round(0.682 * (n_points - 1)) - 1]

        # change point detection
        algo = rpt.Dynp(model='l1').fit(time_series)
        cp[ii] = algo.predict(n_bkps=1)

        for jj in range(n_points):
            # walk time series to find ideal point
            # fitness of one indicates perfect alignment with change point
            fitness[ii][jj] = (n_points - abs(jj - cp[ii]))/n_points

    # sum all rows (fitnesses per point) together
    fit_sum_arr = sum(fit)

    # get index of point that yields the highest overall fitness
    best = np.argmax(fit_sum_arr)

    return best

