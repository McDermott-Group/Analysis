"""

"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, OneStateCleanDirty,QPTunneling_Harrison
import matplotlib.pyplot as plt
import numpy as np
#See if this works -- curve fitting with parameter uncertainties
import scipy.odr
import scipy.stats
from scipy.optimize import curve_fit


"""Q4"""
QP_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/{}/CW20180514A_Ox2/{}/{}/MATLABData/{}')
cr=1
n=25
igf=False
excludedFiles=[]
# Poisoning thru continuous drive @R3
# First file -
# user = 'DaveHarrison'
# date = '07-08-21'
# temp = '76mK'
# qubit = 'Q4'

# Second file -
# user = 'LIU'
# date = '07-07-21'
# temp = '76mK'
# qubit ='Q4'

# Radiator experiments
# qubit ='Q2'
# user = 'LIU'
# date = 'PSD2021Jun28'
# temp=='451mK' #NOGO
# temp = '501mK' #NOGO
# cr=0.05
# n=10
# igf=False
# ep=5
# n=10
# ep=500
# temp='404mK' #Could be tuned
# ep = 10
# temp = '372mK'
# temp = '331mK'
# ep = 2
# excludedFiles=range(0,10)
# temp = '310mK'#^^^ exclude files ^^^

# For the following--need to concatenate several records
# cr=5
# n=10
# temp = '284mK'
# temp = '251mK'
# temp = '199mK'
# temp = '102mK'
# temp = '76mK'

qubit ='Q2'
user = 'LIU'
date = 'PSDTemp2021Jul31'
ep=10
cr=0.5
n=10
# temp = '199mK'
# temp = '252mK'
# temp = '275mK'
# temp = '300mK'
# temp = '325mK'
# temp = '352mK'
# cr=1
# n=10
# temp = '375mK'
# temp='396mK'
# cr=0.2
# temp='428mK'

# qubit ='Q4'
# user = 'LIU'
# date = 'PSD2021Jun28'
# # temp='501mK' #NOGO
# temp = '404mK'
# cr=0.05
# n=10
# igf=True
# ep=1
# ep=100
# n=5
# temp= '451mK'  #Could be tuned

# ep =10
# temp = '372mK'
# ep =5
# excludedFiles = range(10,20)
# temp = '331mK' # ^^^ exclude files ^^^
# temp = '310mK'
# temp = '284mK'
# excludedFiles = range(30,40)
# temp = '251mK' # ^^^ exclude files ^^^
# temp = '199mK'
# excludedFiles = range(40,50)
# temp = '102mK'# ^^^ exclude files ^^^
# temp = '76mK'

# qubit ='Q4'
# user = 'LIU'
# date = 'PSDTemp2021Jul31'
# ep=10
# cr=0.2
# n=10
# temp = '199mK'
# temp = '252mK'
# temp = '275mK'
# temp = '300mK'
# temp = '325mK'
# cr=1
# n=10
# temp = '375mK'
# temp='396mK'
# ep=5
# temp='428mK'

# qubit ='Q3'
# user = 'LIU'
# date = 'PSD2021Jun28'
# ep=2
# temp = '284mK' #Breaking everything up into 1 part instead of 20 gives something reasonable-ish but I do not trust it.
# temp = '199mK' #Same --giving up on these.

qubit ='Q1'
user = 'LIU'
date = 'PSD2021Jun28'
# temp='501mK' #NOGO
# temp='451mK'  #This looks OK when set up like this. Fits still questionable, though
# cr=0.05
# n=10
# igf=True
# ep=1
# cr=0.05
# n=10
# igf=True
# ep=1
# temp='404mK' #This looks OK when set up like this. Fits still questionable, though
# n=25
# ep = 50
# igf=True
#--------
#I *THINK* the problem here has to do with the low freq 1/f noise messing up the area conservation
#Do not reqiure fidelity mapping model
# temp = '372mK'
# temp = '331mK'
#--------
ep=1
n=20
# excludedFiles=range(90,100)
# temp = '310mK' #^^^ exclude files ^^^ Data and fits good
# ep = 10
# n=20
# excludedFiles=range(0,20)
# temp = '284mK' # ^^^ exclude files ^^^ Data and fits good
# ep=10
# n=25
# temp = '251mK'  #Data and fits good
# temp = '199mK' #Data and fits good
# temp = '102mK' #Data and fits good
# temp = '76mK' #Data and fits good

# qubit ='Q1'
# user = 'LIU'
# date = 'PSDTemp2021Jul31'
# ep=0
# cr=0.2
# n=10
# temp = '199mK'
# temp = '252mK'
# temp = '275mK'
# temp = '300mK'
# temp = '325mK'
# temp = '352mK'
# ep=10
# cr=0.2
# n=10
# temp = '375mK'
# temp='396mK'
# ep=5
# temp='428mK'



experiment_name_PSD = (qubit+'_PSD_'+temp)


def fit_PSD_target_function(f, f_parity, F_map):
    return (4 * F_map ** 2 * f_parity) / ((2 * f_parity) ** 2 + (2 * np.pi * f) ** 2) + (1 - F_map ** 2) *(50e-6)#0.0015

def fit_wrapper_for_odr(beta, f): # parameter order for odr
    return fit_PSD_target_function(f, *beta)

def fit(psd,f,f_guess,F_guess):
    psd = psd[~np.isnan(f)]
    f = f[~np.isnan(f)]
    f = f[~np.isnan(psd)]
    psd = psd[~np.isnan(psd)]

    parameters, cov = curve_fit(fit_PSD_target_function, f, psd, bounds=[(0.001, F_guess/2), (1e6, 2*F_guess)], p0=[f_guess, F_guess],method='trf',sigma=f**1)

    model = scipy.odr.odrpack.Model(fit_wrapper_for_odr)
    data = scipy.odr.odrpack.Data(f, psd)
    myodr = scipy.odr.odrpack.ODR(data, model, beta0=parameters, maxit=0)
    myodr.set_job(fit_type=2)
    parameterStatistics = myodr.run()
    df_e = len(f) - len(parameters)  # degrees of freedom, error
    cov_beta = parameterStatistics.cov_beta  # parameter covariance matrix from ODR
    sd_beta = parameterStatistics.sd_beta * parameterStatistics.sd_beta
    ci = []
    t_df = scipy.stats.t.ppf(0.975, df_e)
    ci = []
    for i in range(len(parameters)):
        ci.append([parameters[i] - t_df * parameterStatistics.sd_beta[i],
                   parameters[i] + t_df * parameterStatistics.sd_beta[i]])

    # tstat_beta = parameters / parameterStatistics.sd_beta  # coeff t-statistics
    # pstat_beta = (1.0 - scipy.stats.t.cdf(np.abs(tstat_beta), df_e)) * 2.0  # coef. p-values
    #
    # for i in range(len(parameters)):
    #     print('parameter:', parameters[i])
    #     print('   conf interval:', ci[i][0], ci[i][1])
    #     print('   tstat:', tstat_beta[i])
    #     print('   pstat:', pstat_beta[i])
    #     print()
    return [[parameters[0], ci[0][0], ci[0][1]],[parameters[1], ci[1][0], ci[1][1]]]





# kneeFreqs = []
#
# delta=4
# for fileNo in np.arange(0,100,delta):
#     PSD_file = [QP_path.format(user, date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in range(fileNo,fileNo+delta)]
#     QPT = QPTunneling_Wilen()#start fresh each iteration
#     QPT.add_datasets(PSD_file)
#
#     psd, f = QPT.get_psd(window_averaging=True)
#     psd_fit, f_fit = QPT.get_fit() #need to actually LOOK at these fits and see if they make sense
#     print fit(psd,f,1/QPT.params[0],QPT.params[1])
#     print(QPT.params[0])
#     plt.loglog(f,psd,'--')
#     plt.loglog(f_fit, psd_fit, '-')
#     plt.show()
#     # raw_input("Press Enter to continue...")
#     if fileNo>0:
#         kneeFreqs.append(1/QPT.params[0])


kneeFreqs = []

QPT = QPTunneling_Harrison()

#This is to exclude files where the charge stabilization has failed.
for i in range(50):
    if i not in excludedFiles:
        PSD_file = [QP_path.format(user, date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i)]
        QPT.add_datasets(PSD_file)

fitBool=True
fidelity = []
psd, f1 = QPT.get_psd(number=n,window_averaging=True,concatenate_records=0.1)
fit,f2 = QPT.get_fit(excluded_points=1,ignore_fidelity=False)
for i in range(0,len(psd)):
    plt.loglog(f1, psd[i], '--')
    plt.loglog(f2, fit[i], '-')
    plt.show()

print(len(psd))
print(n)

# if fitBool:
#     psd_fit, f_fit = QPT.get_fit(excluded_points=ep,ignore_fidelity=igf)
#     fidelity=QPT.fidelity
#
# for i in range(n):
#     if fitBool:
#       #  plt.title('Loop #{:03d} f_parity={:.2e}'.format(i,1/QPT.T_parity[i]))
#         pass
#     else:
#       #  plt.title('Loop #{:03d}'.format(i))
#         pass
#     plt.loglog(f,psd[i],'--')
#     if fitBool:
#         plt.loglog(f_fit, psd_fit[i], '-')
#     plt.show()
# if fitBool:
#     f_parity=np.reciprocal(QPT.T_parity)
#     print('{:.2e}'.format(np.mean(f_parity)))
#     print('{:.2e}'.format(np.std(f_parity)/np.sqrt(len(f_parity))))
#     print('{:.2e}'.format(np.mean(fidelity)))

