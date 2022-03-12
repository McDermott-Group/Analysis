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
date = '03-04-22'
# QBs = ['Q1','Q2','Q4']
QBs=['Q1']
#J3Biaslist = np.arange(35000,80001,500)
J3Biaslist1= []#list(np.arange(-300, 441, 20))
J3Biaslist2= list(np.arange(450,690,20))+list(np.arange(-305,55,20))#list(np.arange(480, 681, 20))+ list(np.arange(-310, 11, 20))
start_file = 0
num_files = 50
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
    for J3Bias in J3Biaslist1:
        #Really JB3--mislabeled in experiment name
        experiment_name_PSD = ('PSD_{}mDAC_JB3_Q{}_HighRO'.format(str(J3Bias),str(1)))
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
        # print (J3Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
        fparity[QB_id].append(avg_parity)
        uparity[QB_id].append(parity_uncertainty)
        fidelity[QB_id].append(avg_fidelity)

    for J3Bias in J3Biaslist2:
        date='03-05-22'
        #Really JB3--mislabeled in experiment name
        experiment_name_PSD = ('PSD_{}mDAC_JB3_Q{}_HighRO'.format(str(J3Bias),str(1)))
        PSD_file_Number = np.arange(start_file, start_file + num_files, 1)
        PSD_file = [QP_path.format(date, experiment_name_PSD, experiment_name_PSD) + '_{:03d}.mat'.format(i) for i in PSD_file_Number]
        QPT_Q = QPTunneling_Harrison(name='{} with J3 = {} mDAC HighRO'.
                                  format(QB_id, str(J3Bias)))
        QPT_Q.add_datasets(PSD_file)

        cr = 0.1;
        ep = 1;


        avg_fidelity, avg_parity, parity_uncertainty =plotFittedPSD_Harrison(QPT_Q, save=True, name='{} with J3 ={} mDAC HighRO'.
                                  format(QB_id, str(J3Bias)),excluded_points=ep,concatenate_records=cr,ylim=[10 ** (-5), 10 ** (-3)])
        # print(QB_id)
        # print (J3Bias, '[{:.1f}]'.format(QPT_Q.params[0]))
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
    for i,bias in enumerate(J3Biaslist1+J3Biaslist2):
        if fidelity[QB_id][i] > 0.2 and uparity[QB_id][i]/fparity[QB_id][i]<0.15 and fparity[QB_id][i]>100 and fparity[QB_id][i]<4000 and fparity[QB_id][i] !=0:
            filtered_freq.append(bias)
            filtered_parity.append(fparity[QB_id][i])
    plt.plot(filtered_freq,filtered_parity)
    np.set_printoptions(threshold=sys.maxsize)
    print(filtered_freq)
    print(filtered_parity)
    plt.ylim(10**3,2*10**3)
    plt.show()

#[-300, -280, -260, -240, -220, -200, -180, -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, -310, -290, -270, -250, -230, -210, -190, -170, -150, -130, -110, -90, -70, -50, -30, -10, 10]
#[1736.925704875256, 1652.1843071952192, 1569.16269759059, 1505.6833617333837, 1446.419551769339, 1436.9502444342097, 1403.3547039652708, 1333.5526261097511, 1332.8428555673117, 1386.60184615875, 1534.9658783761906, 1438.5249286831695, 1352.2854422277042, 1280.2125793492833, 1250.0767446556113, 1248.5489662487528, 1336.9051012324167, 1312.8205885125549, 1452.9709155109026, 1388.513753445553, 1252.5006786510335, 1123.93083683362, 1132.2668364138622, 1107.8159196548684, 1132.142291132424, 1165.8419091619867, 1281.6281572594035, 1492.9486238181275, 1326.6900951383097, 1328.7034624653977, 1320.507993033514, 1236.739277156397, 1255.3488188505798, 1283.1546913242917, 1372.8563710053045, 1517.8609504021365, 1502.316102869668, 1319.7792758249082, 1313.1997117691603, 1431.0021170627776, 1435.544492925471, 1502.8598419516397, 1513.5524518670345, 1565.1670603501445, 1755.719927806276, 1787.8366748817407, 1891.7350174705175, 1975.3763686692294, 2200.2794161211914, 1756.8052706084131, 1687.7385253237892, 1600.9697519077658, 1532.6479247924344, 1505.9766544925692, 1428.6279823720972, 1378.8138407858455, 1345.8957442214194, 1299.2072982624009, 1345.7487240217304, 1508.5294839119488, 1584.7832457325112, 1411.0844678809167, 1312.3438971074315, 1283.3417332533745, 1263.3246362111604, 1312.1136125296669]



