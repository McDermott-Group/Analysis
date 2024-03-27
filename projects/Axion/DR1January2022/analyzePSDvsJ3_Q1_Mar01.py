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
date = '03-03-22'
# QBs = ['Q1','Q2','Q4']
QBs=['Q1']
#J3Biaslist = np.arange(35000,80001,500)
J3Biaslist= list(np.arange(-310,691,20))
start_file = 0
num_files = 10
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
        #Really JB3--mislabeled in experiment name
        experiment_name_PSD = ('PSD_{}mDAC_JB3_Q{}'.format(str(J3Bias),str(1)))
        PSD_file_Number = np.arange(start_file, start_file + num_files, 1)
        PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
        QPT_Q = QPTunneling_Harrison(name='{} with J3 = {} mDAC'.
                                  format(QB_id, str(J3Bias)))
        QPT_Q.add_datasets(PSD_file)

        cr = 0.1;
        ep = 1;


        avg_fidelity, avg_parity, parity_uncertainty =plotFittedPSD_Harrison(QPT_Q, save=True, name='{} with J3 ={} mDAC'.
                                  format(QB_id, str(J3Bias)),excluded_points=ep,concatenate_records=cr,ylim=[10 ** (-5), 10 ** (-3)])
        # print(QB_id)
        # print(J3Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
        fparity[QB_id].append(avg_parity)
        uparity[QB_id].append(parity_uncertainty)
        fidelity[QB_id].append(avg_fidelity)

plt.figure(figsize=(8, 6))
plt.title('Parity Rate vs J3Bias')
plt.xlabel('J3 Bias (mDAC)')
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
    plt.ylim(10**3,2*10**3)
    plt.show()
# plt.savefig('ParityRateVsBias.png')

#[-249.6, -230.4, -211.2, -192.0, -172.8, -153.6, -134.4, -115.2, -96.0, -76.8, -57.6, -38.4, -19.2, 0.0, 19.2, 38.4, 57.6, 76.8, 115.2, 134.4, 153.6, 172.8, 192.0, 211.2, 230.4, 249.6, 268.8]
#[3031.133819894432, 1357.3852132987386, 1448.6700846546637, 1352.7893429066328, 1301.4830579320878, 1296.7120006507155, 1283.4792120965346, 1378.2413441824951, 1434.3583760515965, 1346.8702366894859, 1290.6398156164864, 1339.1554429311686, 1247.9195180199345, 1286.5555832407074, 1346.9099357553719, 1345.0749539886262, 1380.8262730519205, 1313.8228185442708, 1113.3172267273862, 1160.355879012291, 1097.734285280708, 1092.9749132974277, 1285.9483616257749, 1059.5915627872641, 1459.7321828246531, 1416.2866629969835, 1199.468694006187]
#[-220.8, -201.6, -182.4, -163.2, -144.0, -124.8, -105.6, -86.4, -67.2, -48.0, -28.8, 9.6, 86.4, 105.6, 124.8, 144.0, 182.4, 220.8, 297.6, 316.8, 336.0, 355.2, 374.4, 451.2, 470.4, 489.6, 508.8, 528.0, 547.2, 566.4, 604.8, 624.0, 643.2, 662.4]
#[1387.929252192186, 1340.481785642099, 1273.5988547906336, 1276.5478520894958, 1305.8706176815, 1254.1065194945386, 1613.393063861564, 1539.6945712564805, 1335.5120130041705, 1272.248515913974, 1293.6025556727473, 1235.4781424651246, 1201.496705579338, 1212.2721932061008, 1112.8001832084494, 1065.3811704047928, 1202.8076573650978, 1345.5131219274117, 1143.2300901894096, 1256.8976796209513, 1300.0659795877987, 1300.6690722004928, 1455.4181343383655, 1308.459499385584, 1134.022376910409, 1359.3572673645604, 1595.2184926577122, 1463.1677268140484, 1536.3994874449704, 1699.196086335938, 1631.8633950774329, 1884.2199039974669, 1773.0094260491878, 2006.0544577792386]
