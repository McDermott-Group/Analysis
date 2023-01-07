"""
2021Jul
"""
from QPTunneling_LIU import QPTunneling_Wilen, plotMultiFittedPSD, QPTunneling_Liu, QPTunneling_Harrison, \
    OneStateCleanDirty, plotFittedPSD_Harrison
from bisect import bisect
import matplotlib.pyplot as plt
import numpy as np


QP_path = ('Z:/mcdermott-group/data/Antenna/SUXmon/LIU/VitoChip1/{}/{}/MATLABData/{}')
# Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\08-19-21
# date = '08-19-21'
# date = 'JJRadiatorQPT_2021Aug19'
date = '2021Nov09_Q12_PSD_J1Radiator'
QBs = ['Q1','Q2']
fparity = {}
uparity = {}
fidelity={}
# for J1Bias in J1Biaslist:

def get_bias_list(QB_id):
    return [bias for bias in (list(np.arange(0,100,10))+list(np.arange(105,695,5)))]

def get_cr(freq):
    #parity_string_array is length 125, so want to be able to divide that evenly.
    cr_list = [0.04,0.2,1,5,25]
    if freq <400:#prev 100
        return 25
    elif freq < 800:
        return 5
    elif freq < 2000:
        return 1
    elif freq < 4000:
        return 0.2
    return 0.02

    # if 10**2/freq > 1:
    #     return 25
    # print cr_list[bisect(cr_list,10**2/freq)]
    # return cr_list[bisect(cr_list,10**2/freq)]

for QB_id in QBs:
    fparity[QB_id]=[]
    uparity[QB_id]=[]
    fidelity[QB_id]=[]
    # if QB_id=='Q1':
    #     ylim=[10 ** (-5), 2*10 ** (-2)]
    # if QB_id=='Q2':
    #     ylim=[10 ** (-5), 10 ** (-2)]
    # if QB_id=='Q12':
    #     ylim=[10 ** (-5), 10 ** (-2)]
    for J1Bias in get_bias_list(QB_id):
        experiment_name_PSD = (QB_id+'_PSD_'+str(J1Bias)+'mVJ1')
        print(experiment_name_PSD)
        PSD_file_Number = np.arange(0, 25, 1)
        PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
        QPT_Q = QPTunneling_Harrison(name='{} with J1 = {} mV, {}GHz'.
                                  format(QB_id, str(J1Bias), str(J1Bias*4.8)))
        QPT_Q.add_datasets(PSD_file)

        cr = 1
        ep = 8

        avg_fidelity, avg_parity, parity_uncertainty =plotFittedPSD_Harrison(QPT_Q, save=False, name='Zoomed Optimized {} with J1 ={} mV {}GHz'.
                                  format(QB_id, str(J1Bias), str(J1Bias*484)),excluded_points=ep, concatenate_records=cr)#, ylim=ylim)

        cr=get_cr(avg_parity)
        ep=2

        avg_fidelity, avg_parity, parity_uncertainty = plotFittedPSD_Harrison(QPT_Q, save=True,
                                                                              name='Zoomed Optimized {} with J1 ={} mV {}GHz'.
                                                                              format(QB_id, str(J1Bias),
                                                                                     str(J1Bias * 484)),
                                                                              excluded_points=ep,
                                                                              concatenate_records=cr)#, ylim=ylim)

        fparity[QB_id].append(avg_parity)
        uparity[QB_id].append(parity_uncertainty)
        fidelity[QB_id].append(avg_fidelity)

plt.figure(figsize=(8, 6))
plt.title('Parity Rate vs J1Bias')
plt.xlabel('J1Bias (mV)')
plt.ylabel('Parity Rate (Hz)')
for QB_id in QBs:
    plt.plot(get_bias_list(QB_id),fparity[QB_id])
    print('{}_J1Bias='.format(QB_id),get_bias_list(QB_id))
    print('{}_J1ParityRate='.format(QB_id), fparity[QB_id])
    print('{}_J1ParityUncertainty='.format(QB_id), uparity[QB_id])
    print('{}_J1Fidelity='.format(QB_id), fidelity[QB_id])
plt.show()
# plt.savefig('ParityRateVsBias.png')