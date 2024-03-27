"""
2021Jul
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison, \
    OneStateCleanDirty, plotFittedPSD_Harrison
import matplotlib.pyplot as plt
import numpy as np
import sys

"""Q2"""
QP_path = ('Z:/mcdermott-group/data/testProject/Keysight/DCH/NA/{}/{}/MATLABData/{}')
# Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\08-19-21
# date = '08-19-21'
# date = 'JJRadiatorQPT_2021Aug19'
date = '02-14-22'
# QBs = ['Q1','Q2','Q4']
QBs=['Q1']
#J2Biaslist = np.arange(35000,80001,500)
#J2Biaslist=list(np.arange(684,909,4))+list(np.arange(916,1160,4))#list(np.arange(382,559,4))
start_file = 0
num_files =25
J2Biaslist = [0]
fparity = {}
uparity = {}
fidelity={}
# for J2Bias in J2Biaslist:
for QB_id in QBs:
    fparity[QB_id]=[]
    uparity[QB_id]=[]
    fidelity[QB_id]=[]
    if QB_id=='Q1':
        ylim=[10 ** (-5), 2*10 ** (-4)]
    if QB_id=='Q2':
        ylim=[5 * 10 ** (-5), 10 ** (-2)]
    if QB_id=='Q4':
        ylim=[10 ** (-5), 10 ** (-3)]
    for J2Bias in J2Biaslist:
        experiment_name_PSD = ('PSD_'+str(J2Bias)+'mVJB2_'+QB_id)
        PSD_file_Number = np.arange(start_file, start_file + num_files, 1)
        if(J2Bias < 909):
            date = '02-14-22'
        else:
            date = '02-14-22'
        PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
        QPT_Q = QPTunneling_Harrison(name='{} with J2 = {} mV, {}GHz'.
                                  format(QB_id, str(J2Bias), str(0.48*2*(J2Bias-30))))
        QPT_Q.add_datasets(PSD_file)

        cr = 1;
        ep = 1;


        avg_fidelity, avg_parity, parity_uncertainty =plotFittedPSD_Harrison(QPT_Q, save=True, name='00 {} with J2 ={} mV {}GHz'.
                                  format(QB_id, str(J2Bias), str(0.48*2*(J2Bias-30))),excluded_points=ep,concatenate_records=cr,ylim=[10 ** (-5), 10 ** (-3)])
        # print(QB_id)
        # print(J2Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
        fparity[QB_id].append(avg_parity)
        uparity[QB_id].append(parity_uncertainty)
        fidelity[QB_id].append(avg_fidelity)

plt.figure(figsize=(8, 6))
plt.title('Parity Rate vs J2Bias')
plt.xlabel('J2 Radiator Frequency (GHz)')
plt.ylabel('Parity Rate (Hz)')
plt.yscale('log')
for QB_id in QBs:
    filtered_freq = []
    filtered_parity = []
    for i,bias in enumerate(J2Biaslist):
        if fidelity[QB_id][i] > 0.2 and uparity[QB_id][i]/fparity[QB_id][i]<0.15 and fparity[QB_id][i]>100 and fparity[QB_id][i]<4000:
            filtered_freq.append(4.8*2*(bias-30)/10)
            filtered_parity.append(fparity[QB_id][i])
    plt.plot(filtered_freq,filtered_parity)
    np.set_printoptions(threshold=sys.maxsize)
    print(filtered_freq)
    print(filtered_parity)
    plt.ylim(10**2,10**4)
    plt.show()
# plt.savefig('ParityRateVsBias.png')