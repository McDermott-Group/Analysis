# 6/16/23
# Gabe Spahn
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, Model, report_fit, fit_report
from dataChest import *
from labrad import units

device = ''
# date = '08-04-23'
tstart = 0.2  # us, relative to end of pulse
tstop = 1600  # us
n_periods = 20  # number of IF periods per window of rolling average
tau_guess = 100  # us, guess for characteristic decay time
ind = -45
plot_all = True
save = False

# path = r'{}\{}\{}\{}'.format('cavityResonators\V1\Gabe',device,date,'Ringdown_3d_cavity')
# if save:
#     savepath = r'{}\fits'.format(path)
# else:
#     savepath = None


def get_raw_data(path, ind=-2, tstart=5, tstop=150, tof_adjust=-10, paramlist=[]):
    # Import data
    d = dataChest(path)
    # print(d.ls())
    rd_dataset_name = d.ls()[0][ind]
    d.openDataset(rd_dataset_name)
    return_params = []
    for param in paramlist:
        num = d.getParameter(param)
        print(param + ' is: {}'.format(num))
        return_params.append(num)
    ru_dataset_name = d.ls()[0][ind - 1]
    print("Reading out from dataset " + str(ru_dataset_name))
    d.openDataset(ru_dataset_name)
    t = np.asarray(d.getData()[:, dc_get_var_ind(d, 'Time', True)])  # us
    # I = np.asarray(d.getData()[:, 1]) #V
    # Q = np.asarray(d.getData()[:, 2]) #V
    mag = np.asarray(d.getData()[:, dc_get_var_ind(d, 'Amplitude')])  # V
    tstep = t[1] - t[0]  # us
    if tof_adjust <= 0:
        indstart = 0
    else:
        indstart = int(tof_adjust / tstep)
    indstop = int((pulse_len * 1e-3 + tof_adjust) / tstep)

    return [t[indstart:indstop], mag[indstart:indstop], return_params]


def make_path(device, date, basepath='cavityResonators\V1\Gabe', meas_type='Ringup_and_ringdown_3d_cavity'):
    return r'{}\{}\{}\{}'.format('cavityResonators\V1\Gabe', device, date, 'Ringup_and_ringdown_3d_cavity')


def ringdown_fit_bare(t, I, Q, f0=4.8E9, tstart=5, tstop=150, n_roll_avg=20, plot_all=False):
    tstep = t[1] - t[0]  # us
    indstart = int(tstart / tstep)
    indstop = int(tstop / tstep)

    # Center magnitude around x axis, convert to power
    I -= np.mean(I)
    Q -= np.mean(Q)
    mag = np.abs(I + 1j * Q)
    y = 0.02 * (mag[indstart:indstop] ** 2)
    x = t[indstart:indstop]

    if plot_all:
        fig, axs = plt.subplots(1, 3)
        plt.subplot(131)
        plt.plot(t, mag)
        plt.axvline(tstart, color='g', label='fit start')
        plt.axvline(tstop, color='r', label='fit stop')
        plt.legend()
        plt.xlabel("Time since pulse end (us)")
        plt.ylabel("Magnitude w/ arbitrary offset (V)")

    if plot_all:
        plt.subplot(132)
        plt.plot(x, y)
        plt.ylabel("Power (W)")
        plt.xlabel("Time (us)")
    # Try rolling average to get a normal line, then
    # fit to exponential decay with some positive offset as noise floor
    N = n_roll_avg
    y_smooth = np.convolve(y, np.ones(N) / N, mode='valid')
    x_smooth = list(range(len(y_smooth))) * tstep

    def ringdown(params, xs):
        A = params['A'].value
        k = params['k'].value
        y0 = params['y0'].value

        return A * np.exp(-xs * k) + y0

    def ringdown_obj(params):
        return np.abs(y_smooth - ringdown(params, x_smooth))

    # Fit phase vs w
    params = Parameters()
    params.add('A', value=y_smooth[0], min=0.01 * y_smooth[0], max=10 * y_smooth[0])
    params.add('k', value=1 / tau_guess, min=0.002, max=5)
    params.add('y0', value=y_smooth[-1], min=0, max=y_smooth[0])

    out = minimize(ringdown_obj, params, method='leastsq', nan_policy='omit')

    A = out.params['A'].value
    k = out.params['k'].value
    k_err = out.params['k'].stderr
    y0 = out.params['y0'].value

    report_fit(out.params, min_correl=0.5)
    Qtot = 2 * np.pi * f0 / (k * 1E6)  # have to convert linewidth from 1/us
    Qtot_err = Qtot * k_err / k
    print("Qtot = %.3e +/- %.1e" % (Qtot, Qtot_err))

    if plot_all:
        plt.subplot(133)
        plt.plot(x_smooth, y_smooth, '.b', label="Data")
        plt.plot(x_smooth, ringdown(out.params, x_smooth), '-r', label=('Fit: Qtot = %.3e +/- %.1e' % (Qtot, Qtot_err)))
        plt.legend()
        plt.title("N = " + str(N))
        plt.ylabel("Time-averaged power (W)")
        plt.xlabel("Time (us)")
        plt.subplots_adjust(left=0.05, right=0.95)  # get rid of white space
        figManager = plt.get_current_fig_manager()
        # figManager.window.state('zoomed') #maximizes plot size, for convenience and for saving, for TkAgg
        figManager.window.showMaximized()  # maximizes plot size, for convenience and for saving, for Qt
        plt.show()

    result = [Qtot, Qtot_err, pulse_len / 1000000]
    return result


def ringdown_fit_single_old(path, tstart=5, tstop=150, n_periods=1, tau_guess=100, ind=-1, plot_all=True,
                            savepath=None):
    # Import data
    d = dataChest(path)
    print(d.ls())
    dataset_name = d.ls()[0][ind]
    print("Reading out from dataset " + str(dataset_name))
    d.openDataset(dataset_name)
    f0 = d.getParameter("Resonator frequency")[0]
    print(f"Resonator frequency {f0:.7E} GHz")
    LO = d.getParameter("Resonator LO")[0]
    temp_bool = "Temperature" in d.getParameterList()
    pulse_len = d.getParameter("Pulse Length")[0]
    print(f"Pulse length {pulse_len / 1e6} ms")
    if temp_bool:
        temp = d.getParameter("Temperature")[0]
    # n_avgs = d.getParameter("Averages")
    n_avgs = 1
    data = d.getData()
    t = np.asarray(d.getData()[:, 0])  # us
    I = np.asarray(d.getData()[:, 1])  # V
    Q = np.asarray(d.getData()[:, 2])  # V
    # mag = np.asarray(d.getData()[:, 3]) #V
    # mag  = np.asarray(d.getData()[:, 2])
    tstep = t[1] - t[0]  # us
    indstart = int(tstart / tstep)
    indstop = int(tstop / tstep)

    # plot initial data
    if plot_all:
        fig, axs = plt.subplots(1, 3)
        if n_avgs == 1:
            if temp_bool:
                fig.suptitle('{} {}, {}mK ({})'.format(date, device, temp, dataset_name))
            else:
                fig.suptitle('{} {} ({})'.format(date, device, dataset_name))
        else:
            if temp_bool:
                fig.suptitle('{} {}, {}mK {} avgs ({})'.format(date, device, temp, n_avgs, dataset_name))
            else:
                fig.suptitle('{} {}, {} avgs ({})'.format(date, device, temp, n_avgs, dataset_name))

        plt.subplot(131)
        plt.plot(t, I, label='I')
        plt.plot(t, Q, label='Q')
        plt.axvline(tstart, color='g', label='fit start')
        plt.axvline(tstop, color='r', label='fit stop')
        plt.legend()
        plt.xlabel("Time since pulse end (us)")
        plt.ylabel("Response w/ arbitrary offset (V)")

    # Center magnitude around x axis, convert to power
    I -= np.mean(I)
    Q -= np.mean(Q)
    mag = np.abs(I + 1j * Q)
    y = 0.02 * (mag[indstart:indstop] ** 2)
    x = t[indstart:indstop]
    if plot_all:
        plt.subplot(132)
        plt.plot(x, y)
        plt.ylabel("Power (W)")
        plt.xlabel("Time (us)")

    # Try rolling average to get a normal line, then
    # fit to exponential decay with some positive offset as noise floor
    N = int(n_periods / (abs(f0 - LO) * 1E-6 * tstep))
    y_smooth = np.convolve(y, np.ones(N) / N, mode='valid')
    x_smooth = list(range(len(y_smooth))) * tstep

    def ringdown(params, xs):
        A = params['A'].value
        k = params['k'].value
        y0 = params['y0'].value

        return A * np.exp(-xs * k) + y0

    def ringdown_obj(params):
        return np.abs(y_smooth - ringdown(params, x_smooth))

    # Fit phase vs w
    params = Parameters()
    params.add('A', value=y_smooth[0], min=0.01 * y_smooth[0], max=10 * y_smooth[0])
    params.add('k', value=1 / tau_guess, min=1E-6, max=5)
    params.add('y0', value=y_smooth[-1], min=0, max=y_smooth[0])

    out = minimize(ringdown_obj, params, method='leastsq', nan_policy='omit')

    A = out.params['A'].value
    k = out.params['k'].value
    k_err = out.params['k'].stderr
    y0 = out.params['y0'].value

    report_fit(out.params, min_correl=0.5)
    Qtot = 2 * np.pi * f0 / (k * 1E6)  # have to convert linewidth from 1/us
    Qtot_err = Qtot * k_err / k
    print("Qtot = %.3e +/- %.1e" % (Qtot, Qtot_err))
    if plot_all:
        plt.subplot(133)
        plt.plot(x_smooth, y_smooth, '.b', label="Data")
        plt.plot(x_smooth, ringdown(out.params, x_smooth), '-r', label=('Fit: Qtot = %.3e +/- %.1e' % (Qtot, Qtot_err)))
        plt.legend()
        plt.title("N = " + str(N))
        plt.ylabel("Time-averaged power (W)")
        plt.xlabel("Time (us)")
        plt.subplots_adjust(left=0.05, right=0.95)  # get rid of white space
        figManager = plt.get_current_fig_manager()
        # figManager.window.state('zoomed') #maximizes plot size, for convenience and for saving, for TkAgg
        figManager.window.showMaximized()  # maximizes plot size, for convenience and for saving, for Qt
        plt.show()
    if savepath:
        plt.pause(1)  # to allow maximization to effect the savefig
        # plt.savefig(savepath+'test')

    result = [Qtot, Qtot_err, pulse_len / 1000000]
    if temp_bool:
        result.append(temp)
    if n_avgs > 1:
        result.append(n_avgs)
    print(result)
    return result


def ringdown_fit_single(path, tstart=5, tstop=150, n_periods=1, tau_guess=100, ind=-1, plot_all=True, savepath=None):
    # Import data
    d = dataChest(path)
    print(d.ls())
    dataset_name = d.ls()[0][ind]
    print("Reading out from dataset " + str(dataset_name))
    d.openDataset(dataset_name)
    f0 = d.getParameter("Resonator frequency")[0]
    print(f"Resonator frequency {f0 * 1E-9:.7E} GHz")
    LO = d.getParameter("Resonator LO")[0]
    temp_bool = "Temperature" in d.getParameterList()
    pulse_len = d.getParameter("Readout Length")[0]
    print(f"Pulse length {pulse_len / 1e6} ms")
    if temp_bool:
        temp = d.getParameter("Temperature")[0]
    # n_avgs = d.getParameter("Averages")
    n_avgs = 1
    data = d.getData()
    t = np.asarray(data[:, dc_get_var_ind(d, 'Time', True)])  # us
    I = np.asarray(data[:, dc_get_var_ind(d, 'I')])  # V
    Q = np.asarray(data[:, dc_get_var_ind(d, 'Q')])  # V
    # mag = np.asarray(d.getData()[:, 3]) #V
    # mag  = np.asarray(d.getData()[:, 2])
    tstep = t[1] - t[0]  # us
    indstart = int(tstart / tstep)
    indstop = int(tstop / tstep)

    # plot initial data
    if plot_all:
        fig, axs = plt.subplots(1, 3)
        if n_avgs == 1:
            if temp_bool:
                fig.suptitle('{} {}, {}mK ({})'.format(date, device, temp, dataset_name))
            else:
                fig.suptitle('{} {} ({})'.format(date, device, dataset_name))
        else:
            if temp_bool:
                fig.suptitle('{} {}, {}mK {} avgs ({})'.format(date, device, temp, n_avgs, dataset_name))
            else:
                fig.suptitle('{} {}, {} avgs ({})'.format(date, device, temp, n_avgs, dataset_name))

        plt.subplot(131)
        plt.plot(t, I, label='I')
        plt.plot(t, Q, label='Q')
        plt.axvline(tstart, color='g', label='fit start')
        plt.axvline(tstop, color='r', label='fit stop')
        plt.legend()
        plt.xlabel("Time since pulse end (us)")
        plt.ylabel("Response w/ arbitrary offset (V)")

    # Center magnitude around x axis, convert to power
    I -= np.mean(I)
    Q -= np.mean(Q)
    mag = np.abs(I + 1j * Q)
    y = 0.02 * (mag[indstart:indstop] ** 2)
    x = t[indstart:indstop]
    if plot_all:
        plt.subplot(132)
        plt.plot(x, y)
        plt.ylabel("Power (W)")
        plt.xlabel("Time (us)")
    '''
    #Try rolling average to get a normal line, then
    #fit to exponential decay with some positive offset as noise floor
    N = int(n_periods/(abs(f0-LO)*1E-6*tstep))
    y_smooth = np.convolve(y, np.ones(N)/N, mode = 'valid')
    x_smooth = range(len(y_smooth))*tstep
    '''

    def ringdown(params, xs):
        A = params['A'].value
        k = params['k'].value
        y0 = params['y0'].value
        phi0 = params['phi0'].value
        df = params['df'].value

        return A * np.exp(-xs * k) * (np.sin(phi0 + (np.abs(f - LO) + df) * xs / (2 * np.pi)) ** 2) + y0

    def ringdown_obj(params):
        return np.abs(y - ringdown(params, x))

    # Fit phase vs w
    params = Parameters()
    params.add('A', value=np.max(y), min=0.01 * np.max(y), max=10 * np.max(y))
    params.add('k', value=1 / tau_guess, min=1E-4, max=5)
    params.add('y0', value=y[-1], min=0, max=y[0])
    params.add('phi0', value=0, min=-np.pi, max=np.pi)
    params.add('df', value=1, min=0, max=1E3)

    out = minimize(ringdown_obj, params, method='leastsq', nan_policy='omit')

    A = out.params['A'].value
    k = out.params['k'].value
    k_err = out.params['k'].stderr
    y0 = out.params['y0'].value

    report_fit(out.params, min_correl=0.5)
    Qtot = 2 * np.pi * f0 / (k * 1E6)  # have to convert linewidth from 1/us
    Qtot_err = Qtot * k_err / k
    print("Qtot = %.3e +/- %.1e" % (Qtot, Qtot_err))
    if plot_all:
        plt.subplot(133)
        plt.plot(x_smooth, y_smooth, '.b', label="Data")
        plt.plot(x_smooth, ringdown(out.params, x_smooth), '-r', label=('Fit: Qtot = %.3e +/- %.1e' % (Qtot, Qtot_err)))
        plt.legend()
        plt.title("N = " + str(N))
        plt.ylabel("Time-averaged power (W)")
        plt.xlabel("Time (us)")
        plt.subplots_adjust(left=0.05, right=0.95)  # get rid of white space
        figManager = plt.get_current_fig_manager()
        # figManager.window.state('zoomed') #maximizes plot size, for convenience and for saving, for TkAgg
        figManager.window.showMaximized()  # maximizes plot size, for convenience and for saving, for Qt
        plt.show()
    # if savepath:

    # plt.pause(1) #to allow maximization to effect the savefig
    # plt.savefig(savepath+'test')

    result = [Qtot, Qtot_err, pulse_len / 1000000]
    if temp_bool:
        result.append(temp)
    if n_avgs > 1:
        result.append(n_avgs)
    print(result)
    return result


def dc_get_var_ind(open_dataset, var_name, indep=False):
    varlist = open_dataset.getVariables()[int(not (indep))]
    for i, val in enumerate(varlist):
        if val[0] == var_name:
            return i + int(not (indep)) * len(open_dataset.getVariables()[0])


def ringup_fit_single(path, tof_adjust=-10, n_periods=100, Qc_guess=50e6, Qi_guess=10E6, ind=-2, plot_all=True,
                      savepath=None):
    # Import data
    date = get_date_from_path(path)
    print(path)
    d = dataChest(path)
    print(d.ls())
    ru_dataset_name = d.ls()[0][ind]
    d.openDataset(ru_dataset_name)
    print("Reading out from dataset " + str(ru_dataset_name))
    f0 = d.getParameter("Resonator frequency")[0]  # Hz
    print(f"Resonator frequency {f0 * 1E-9:.7f} GHz")
    LO = d.getParameter("Resonator LO")[0]  # Hz
    temp_bool = "Temperature" in d.getParameterList()
    pulse_len = d.getParameter("Pulse Length")[0]  # ns
    print(f"Pulse length {pulse_len / 1e6} ms")
    if temp_bool:
        temp = d.getParameter("Temperature")[0]
    try:
        n_avgs = d.getParameter("Averages")
    except:
        n_avgs = 1
    # ru_dataset_name = d.ls()[0][ind-1]
    # print("Reading out from dataset "+ str(ru_dataset_name))
    # d.openDataset(ru_dataset_name)
    data = d.getData()
    t = np.asarray(data[:, dc_get_var_ind(d, 'Time', True)])  # us
    I = np.asarray(data[:, dc_get_var_ind(d, 'I')])  # V
    Q = np.asarray(data[:, dc_get_var_ind(d, 'Q')])  # V

    tstep = t[101] - t[100]  # us
    if tof_adjust <= 0:
        indstart = 0
    else:
        indstart = int(tof_adjust / tstep)
    indstop = int((pulse_len * 1e-3 + tof_adjust) / tstep)
    # plot initial data
    if plot_all:
        fig, axs = plt.subplots(2, 2)
        if n_avgs == 1:
            if temp_bool:
                fig.suptitle('{} {}, {}mK ({})'.format(date, device, temp, ru_dataset_name))
            else:
                fig.suptitle('{} {} ({})'.format(date, device, ru_dataset_name))
        else:
            if temp_bool:
                fig.suptitle('{} {}, {}mK {} avgs ({})'.format(date, device, temp, n_avgs, ru_dataset_name))
            else:
                fig.suptitle('{} {}, {} avgs ({})'.format(date, device, temp, n_avgs, ru_dataset_name))

        plt.subplot(221)
        plt.plot(t, I, '.b', label='I')
        plt.plot(t, Q, '.b', label='Q')
        plt.axvline(tof_adjust, color='g', label='fit start')
        plt.axvline(pulse_len * 1e-3 + tof_adjust, color='r', label='fit stop')
        plt.legend(loc=1)
        plt.xlabel("Time (us)")
        plt.ylabel("Magnitude w/ arbitrary offset (V)")

    # Center magnitude around x axis, convert to power
    I -= np.mean(I)
    Q -= np.mean(Q)
    mag = np.abs(I + 1j * Q)
    y = 0.02 * (mag[indstart:indstop] ** 2)
    x = t[indstart:indstop]
    if plot_all:
        plt.subplot(222)
        plt.plot(x, y)
        plt.ylabel("Power (W)")
        plt.xlabel("Time (us)")

    # Try rolling average to get a normal line, then
    # fit to exponential decay with some positive offset as noise floor
    N = int(n_periods / (abs(f0 - LO) * 1E-6 * tstep))
    y_smooth = np.convolve(y, np.ones(N) / N, mode='valid')
    x_smooth = list(range(len(y_smooth))) * tstep

    def ringup(params, xs):
        Qc = params['Qc'].value
        Qi = params['Qi'].value
        # Pf = params['Pf'].value
        d = params['d'].value
        f = f0
        Pf = y_smooth.max()

        prefactor = Pf / ((2 * Qc * d * Qi) ** 2 + (f * (Qc + Qi)) ** 2)
        line1 = (2 * Qc * Qi * d) ** 2 + 4 * np.exp(-2 * np.pi * f * xs * (Qc + Qi) / (Qc * Qi)) * ((Qi * f) ** 2) + (
                    f * (Qc - Qi)) ** 2
        lines23 = 4 * Qi * f * np.exp(-1 * np.pi * f * xs * (Qc + Qi) / (Qc * Qi)) * (
                    (Qc - Qi) * f * np.cos(2 * np.pi * d * xs) - 2 * Qc * Qi * d * np.sin(2 * np.pi * d * xs))
        return prefactor * (line1 + lines23)

    def ringup_obj(params):
        return np.abs(y_smooth - ringup(params, x_smooth * 1e-6))

    # Fit phase vs w
    params = Parameters()
    # params.add('Pf', value= y_smooth.max(), min=0.1 * y_smooth.max(), max=2*y_smooth.max())
    params.add('Qc', value=Qc_guess, min=5E6, max=150E6)
    params.add('Qi', value=Qi_guess, min=2E6, max=100E6)
    params.add('d', value=0, min=-500, max=500)
    print("starting fit")
    out = minimize(ringup_obj, params, method='leastsq', nan_policy='omit')
    print("done with fit")
    Qi = out.params['Qi'].value
    Qc = out.params['Qc'].value
    Qi_err = out.params['Qi'].stderr
    Qc_err = out.params['Qc'].stderr
    d = out.params['d'].value
    # Pf = out.params['Pf'].value

    print(Qi_err)
    print(Qc_err)
    print(type(Qi_err))
    print(type(Qc_err))

    print(fit_report(out))
    Qtot = 1 / (1 / Qc + 1 / Qi)  # have to convert linewidth from 1/us

    if plot_all:
        plt.subplot(223)
        plt.plot(x_smooth, y_smooth, '.b', label="Data")
        plt.plot(x_smooth, ringup(out.params, x_smooth * 1e-6), '-r', label=('Fit: Qtot = %.3e +/- %.1e' % (Qtot, 5)))
        plt.legend(loc=1)
        plt.title("N = " + str(N))
        plt.ylabel("Time-averaged power (W)")
        plt.xlabel("Time (us)")
        plt.subplots_adjust(left=0.05, right=0.95)  # get rid of white space
        figManager = plt.get_current_fig_manager()
        # figManager.window.state('zoomed') #maximizes plot size, for convenience and for saving, for TkAgg
        figManager.window.showMaximized()  # maximizes plot size, for convenience and for saving, for Qt
        plt.show()
    # if savepath:
    # plt.pause(1) #to allow maximization to effect the savefig
    # plt.savefig(savepath+'test')
    Qtot_err = np.sqrt((Qi_err / (Qi ** 2)) ** 2 + (Qc_err / (Qc ** 2)) ** 2) / (Qtot ** 2)
    print("d = %.3f Hz" % d)
    print("Qi = %.3e +/- %.1e, Qc = %.3e +/- %.1e" % (Qi, Qi_err, Qc, Qc_err))
    print("Qtot = %.3e +/- %.1e" % (Qtot, Qtot_err))
    result = [Qtot, Qtot_err, pulse_len / 1e6, Qi, Qc]
    if temp_bool:
        result.append(temp)
    if n_avgs > 1:
        result.append(n_avgs)
    print(result)
    return result


if __name__ == "__main__":
    ringdown_fit_single_old(path, tstart, tstop, n_periods, tau_guess, ind, plot_all, savepath)
    # ringdown_fit_single(path, n_periods = 100, ind = ind, plot_all = True)
