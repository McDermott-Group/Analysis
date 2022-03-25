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
date = '02-19-22'
# QBs = ['Q1','Q2','Q4']
QBs=['Q3']
#J3Biaslist = np.arange(35000,80001,500)
J3Biaslist=list(np.arange(700,-301,-20))+list(np.arange(570,-301,-20))
start_file = 0
num_files =10
# J3Biaslist = [40, 42, 44]
fparity = {}
uparity = {}
fidelity={}
# for J3Bias in J3Biaslist:
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
    for J3Bias in J3Biaslist:
        experiment_name_PSD1 = ('TestPSD_'+str(J3Bias)+'mDAC_JB3_'+QB_id+'_Narrow')
        experiment_name_PSD2 = ('TestPSD_' + str(J3Bias) + 'mDAC_JB3_' + QB_id + '_Narrow')
        PSD_file_Number = np.arange(start_file, start_file + num_files, 1)
        PSD_file = [QP_path.format(date, experiment_name_PSD1, experiment_name_PSD2) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
        QPT_Q = QPTunneling_Harrison(name='{} with J3 = {} mDAC, {}GHz'.
                                  format(QB_id, str(J3Bias), str(0.48*2*(J3Bias-30))))
        QPT_Q.add_datasets(PSD_file)

        cr = 1;
        ep = 1;


        avg_fidelity, avg_parity, parity_uncertainty =plotFittedPSD_Harrison(QPT_Q, save=True, name='00 {} with J3 ={} mDAC {}GHz'.
                                  format(QB_id, str(J3Bias), str(0.48*2*(J3Bias-30))),excluded_points=ep,concatenate_records=cr,ylim=[10 ** (-5), 10 ** (-3)])
        # print(QB_id)
        # print (J3Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
        fparity[QB_id].append(avg_parity)
        uparity[QB_id].append(parity_uncertainty)
        fidelity[QB_id].append(avg_fidelity)

plt.figure(figsize=(8, 6))
plt.title('Parity Rate vs J3Bias')
plt.xlabel('J3 Radiator Frequency (GHz)')
plt.ylabel('Parity Rate (Hz)')
plt.yscale('log')
for QB_id in QBs:
    filtered_freq = []
    filtered_parity = []
    for i,bias in enumerate(J3Biaslist):
        if fidelity[QB_id][i] > 0.2 and uparity[QB_id][i]/fparity[QB_id][i]<0.15 and fparity[QB_id][i]>100 and fparity[QB_id][i]<4000 and fparity[QB_id][i] !=0:
            filtered_freq.append(4.8*2*(bias)/10)
            filtered_parity.append(fparity[QB_id][i])
    plt.plot(filtered_freq,filtered_parity)
    np.set_printoptions(threshold=sys.maxsize)
    print(filtered_freq)
    print(filtered_parity)
    plt.ylim(10**2,10**4)
    plt.show()
# plt.savefig('ParityRateVsBias.png')