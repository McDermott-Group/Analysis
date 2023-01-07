from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')

# Z:\mcdermott-group\data\Axion\2022-09-22 - DR2\DCH\Axion3A\09-29-22\PSD_Q3_J2B
biases_data =[-0.3, -0.295, -0.29, -0.285, -0.27999999999999997, -0.27499999999999997, -0.26999999999999996, -0.26499999999999996, -0.25999999999999995, -0.25499999999999995, -0.24999999999999994, -0.24499999999999994, -0.23999999999999994, -0.23499999999999993, -0.22999999999999993, -0.22499999999999992, -0.21999999999999992, -0.2149999999999999, -0.2099999999999999, -0.2049999999999999, -0.1999999999999999, -0.1949999999999999, -0.1899999999999999, -0.1849999999999999, -0.17999999999999988, -0.17499999999999988, -0.16999999999999987, -0.16499999999999987, -0.15999999999999986, -0.15499999999999986, -0.14999999999999986, -0.14499999999999985, -0.13999999999999985, -0.13499999999999984, -0.12999999999999984, -0.12499999999999983, -0.11999999999999983, -0.11499999999999982, -0.10999999999999982, -0.10499999999999982, -0.09999999999999981, -0.0949999999999998, -0.0899999999999998, -0.0849999999999998, -0.0799999999999998, -0.07499999999999979, -0.06999999999999978, -0.06499999999999978, -0.059999999999999776, -0.05499999999999977, -0.04999999999999977, -0.04499999999999976, -0.03999999999999976, -0.034999999999999754, -0.02999999999999975, -0.024999999999999745, -0.01999999999999974, -0.014999999999999736, -0.009999999999999731, -0.004999999999999727, 2.7755575615628914e-16, 0.005000000000000282, 0.010000000000000286, 0.015000000000000291, 0.020000000000000295, 0.0250000000000003, 0.030000000000000304, 0.03500000000000031, 0.04000000000000031, 0.04500000000000032, 0.05000000000000032, 0.055000000000000326, 0.06000000000000033, 0.06500000000000034, 0.07000000000000034, 0.07500000000000034, 0.08000000000000035, 0.08500000000000035, 0.09000000000000036, 0.09500000000000036, 0.10000000000000037, 0.10500000000000037, 0.11000000000000038, 0.11500000000000038, 0.12000000000000038, 0.1250000000000004, 0.1300000000000004, 0.1350000000000004, 0.1400000000000004, 0.1450000000000004, 0.1500000000000004, 0.15500000000000042, 0.16000000000000042, 0.16500000000000042, 0.17000000000000043, 0.17500000000000043, 0.18000000000000044, 0.18500000000000044, 0.19000000000000045, 0.19500000000000045, 0.20000000000000046, 0.20500000000000046, 0.21000000000000046, 0.21500000000000047, 0.22000000000000047, 0.22500000000000048, 0.23000000000000048, 0.2350000000000005, 0.2400000000000005, 0.2450000000000005, 0.2500000000000005, 0.2550000000000005, 0.2600000000000005, 0.2650000000000005, 0.2700000000000005, 0.2750000000000005, 0.2800000000000005, 0.28500000000000053, 0.29000000000000054, 0.29500000000000054, 0.30000000000000054]
parity_rate_data =[33.0268476068771, 3804.595323648049, 3707.1646870903946, 19.365173054119047, 19.397302147030896, 3403.241999539749, 3293.5345065715724, 3332.448611477359, 3344.9599639975013, 36.62889656772163, 2928.3120598922537, 2912.1330443168363, 33.16783369013689, 3029.8866167461756, 2892.5200765718077, 2883.464249735162, 2653.1000671112984, 2695.3395277685995, 2567.6584720232345, 2446.552778879517, 2291.7701048531044, 2257.9242998885265, 1838.1382281767612, 1799.5283691767497, 1853.2938291520738, 1799.4426642285473, 1682.839035526108, 1657.1064830139428, 1301.4591312294335, 1225.2306733546511, 1151.6627433039548, 1097.241210667024, 1006.3564813343008, 977.8279249976375, 848.7694199723323, 628.1406430923062, 547.3612290363957, 490.03616922692873, 454.5495655222479, 480.4696854627902, 585.952662086949, 686.8485961222626, 872.8017823796807, 1069.6649956667684, 1128.9582172800108, 979.4271331892515, 848.7982986182833, 754.7048590245752, 759.7551145707127, 877.8578333234274, 1063.9634210547606, 1371.658617692911, 1205.5152827757333, 969.0636606636518, 805.1971996912454, 742.0675921633012, 750.0002410597579, 644.8504994457983, 711.700973067462, 658.2372898017547, 585.8863154437788, 359.5088532830582, 280.63189197458485, 281.57341648318976, 278.9469299958782, 224.18402472918117, 263.8739380099864, 221.8897969845285, 249.30599829357905, 286.18934522241466, 278.6799083287889, 255.41353681135234, 261.83083097748727, 308.6750877288074, 367.3485088616239, 372.55264855058664, 644.6019733128367, 701.8065043874678, 724.7134768526228, 776.0070157807079, 835.871313185051, 913.7184970287606, 1038.3464222179903, 1367.0583857562558, 1041.6938014340237, 788.5896208042135, 712.2134406106974, 761.904041109732, 801.4191862492045, 900.1511341187075, 1125.3588857734642, 1014.9865800289626, 791.1944770929969, 681.1726038845375, 561.4032383005333, 538.9002612917874, 542.9581472037697, 586.3227014409117, 541.8730696386116, 678.7268936614348, 944.8322072957125, 1035.4798816983498, 986.4878310552995, 1160.029030043856, 1325.7920272242363, 1218.4853840554433, 1373.4158024089088, 1640.7790976276665, 1618.3606479988898, 1835.3913221553653, 1731.8508375380036, 1726.6490135821293, 1896.4755122537642, 2234.3924480887995, 2399.362211601743, 2306.7739273523594, 2446.296070860676, 2570.6811549257436, 2701.2383855875387, 2968.500214526717, 26.551276916381635]

# Z:\mcdermott-group\data\Axion\2022-09-22 - DR2\DCH\Axion3A\09-29-22\PSD_Q3_J2B_Retightened
biases_data1 =[-0.3, -0.295, -0.29, -0.285, -0.27999999999999997, -0.27499999999999997, -0.26999999999999996, -0.26499999999999996, -0.25999999999999995, -0.25499999999999995, -0.24999999999999994, -0.24499999999999994, -0.23999999999999994, -0.23499999999999993, -0.22999999999999993, -0.22499999999999992, -0.21999999999999992, -0.2149999999999999, -0.2099999999999999, -0.2049999999999999, -0.1999999999999999, -0.1949999999999999, -0.1899999999999999, -0.1849999999999999, -0.17999999999999988, -0.17499999999999988, -0.16999999999999987, -0.16499999999999987, -0.15999999999999986, -0.15499999999999986, -0.14999999999999986, -0.14499999999999985, -0.13999999999999985, -0.13499999999999984, -0.12999999999999984, -0.12499999999999983, -0.11999999999999983, -0.11499999999999982, -0.10999999999999982, -0.10499999999999982, -0.09999999999999981, -0.0949999999999998, -0.0899999999999998, -0.0849999999999998, -0.0799999999999998, -0.07499999999999979, -0.06999999999999978, -0.06499999999999978, -0.059999999999999776, -0.05499999999999977, -0.04999999999999977, -0.04499999999999976, -0.03999999999999976, -0.034999999999999754, -0.02999999999999975, -0.024999999999999745, -0.01999999999999974, -0.014999999999999736, -0.009999999999999731, -0.004999999999999727, 2.7755575615628914e-16, 0.005000000000000282, 0.010000000000000286, 0.015000000000000291, 0.020000000000000295, 0.0250000000000003, 0.030000000000000304, 0.03500000000000031, 0.04000000000000031, 0.04500000000000032, 0.05000000000000032, 0.055000000000000326, 0.06000000000000033, 0.06500000000000034, 0.07000000000000034, 0.07500000000000034, 0.08000000000000035, 0.08500000000000035, 0.09000000000000036, 0.09500000000000036, 0.10000000000000037, 0.10500000000000037, 0.11000000000000038, 0.11500000000000038, 0.12000000000000038, 0.1250000000000004, 0.1300000000000004, 0.1350000000000004, 0.1400000000000004, 0.1450000000000004, 0.1500000000000004, 0.15500000000000042, 0.16000000000000042, 0.16500000000000042, 0.17000000000000043, 0.17500000000000043, 0.18000000000000044, 0.18500000000000044, 0.19000000000000045, 0.19500000000000045, 0.20000000000000046, 0.20500000000000046, 0.21000000000000046, 0.21500000000000047, 0.22000000000000047, 0.22500000000000048, 0.23000000000000048, 0.2350000000000005, 0.2400000000000005, 0.2450000000000005, 0.2500000000000005, 0.2550000000000005, 0.2600000000000005, 0.2650000000000005, 0.2700000000000005, 0.2750000000000005, 0.2800000000000005, 0.28500000000000053, 0.29000000000000054, 0.29500000000000054, 0.30000000000000054]
parity_rate_data1 = [34.4645025401347, 40.52850710748375, 38.33124490311449, 31.61707126834793, 29.370404211570573, 3243.1006900091547, 3356.7004670448323, 25.227754787612827, 29.285861947081703, 3147.2624045120424, 25.989471607356542, 3141.6226196906064, 32.107368323788535, 3182.490671436952, 23.320984935339382, 53.58135170372689, 2956.5639140986264, 3162.7214544970234, 2977.614496723117, 2733.0460405472204, 29.836661200868992, 46.02710067983333, 35.63199141243225, 37.94435876444115, 43.60728327187585, 2303.909671969774, 2335.1783689362164, 37.41237358175952, 244.3873436412472, 1678.425483735762, 1592.607867736693, 1914.9447112965108, 1829.1820658215122, 1538.6662078788274, 1512.716141038692, 1356.716769833392, 1018.8415218783516, 879.5013291431165, 1072.222563365139, 924.8291913931255, 952.1631425714805, 1039.9048460777306, 1207.5838943734366, 1370.0317328233662, 1327.3755481673634, 1362.3406280767558, 1118.5396417754548, 1179.619527516466, 1065.3911031458854, 1264.4284662247749, 1484.4319208039067, 1843.8665502535946, 1551.1586516672762, 1331.0984679675373, 1628.5157303495807, 1632.9228484262687, 1282.8405040831963, 1225.6635576848712, 1296.4017040393835, 1313.3302281349845, 1033.595935784898, 856.0790232838619, 912.8898233902504, 898.6091160237672, 705.9632727477058, 614.0140266443442, 775.5042888868717, 810.5665543144465, 963.5397743670318, 710.2249093246352, 867.9678095029708, 1001.3583486832212, 700.1633244305517, 664.7747642715253, 1006.5081702619409, 814.3526150048681, 951.8325749652516, 948.158702663741, 941.9596221702759, 1015.7915376129669, 1124.62856765482, 1434.8165118836496, 1723.907229206477, 1431.3398554495127, 1231.0159947407722, 1000.226196960157, 1121.775740896926, 1219.0629533880322, 1170.5355844257006, 1060.838486234636, 1388.5501780321429, 1528.9350429757328, 1320.538195987509, 922.1366941832885, 796.5040509422189, 680.689340007343, 878.090894343147, 897.0809371967308, 890.3795490578901, 963.8452836286162, 1210.6235787772541, 1339.2569174276869, 1368.2900917531995, 1649.4661840201052, 1534.4190246694968, 1577.3242426725953, 1650.7920603808468, 1960.9390494080997, 2107.489718320197, 1984.9288633648086, 2034.9018446466753, 2029.6251969244406, 2200.0309296494315, 2341.1013704336174, 2449.5650210842473, 2693.4203640980563, 2869.7318728300706, 2895.5448450120457, 2758.7796625334104, 2981.9635926898513, 2929.9464772363585]

# Z:\mcdermott-group\data\Axion\2022-09-22 - DR2\DCH\Axion3A\09-30-22\PSD_Q3_J2B_OA
biases_data2 = [0.07350000000000015, -0.029399999999999944, 0.00980000000000009, 0.08820000000000015, 0.07840000000000015, 0.0196000000000001, -0.098, 0.004900000000000086, -0.04409999999999996, -0.0931, -0.03429999999999995, -0.06859999999999998, 0.09800000000000017, -0.03919999999999995, -0.07349999999999998, 0.04900000000000013, 0.09310000000000017, 0.014700000000000095, 0.029400000000000107, -0.06369999999999998, -0.08329999999999999, 0.03920000000000012, 0.05880000000000013, 0.06370000000000013, -0.02449999999999994, -0.004899999999999923, -0.07839999999999998, -0.05389999999999996, -0.014699999999999932, -0.019599999999999937, 8.1601392309949e-17, 0.08330000000000015, 0.03430000000000011, -0.05879999999999997, 0.04410000000000012, -0.009799999999999927, -0.0882, 0.0245000000000001, -0.04899999999999996, 0.06860000000000015, 0.05390000000000013]
parity_rate_data2 =[1057.9831043660172, 1168.613156574155, 616.7674760828385, 1918.3012413365855, 1295.6259926300352, 731.4345499526738, 1996.7500579732591, 659.6423382922535, 1390.1811527767502, 1573.1050453799228, 1291.5655899035262, 974.9146118100405, 1956.2740477429804, 1567.2827310263363, 1026.1357642207322, 1257.534365343513, 1743.4383455831642, 668.2023336133393, 1034.324317901774, 1272.2449921965765, 1610.1098263472918, 1440.5145263730774, 1632.3327050360854, 1103.393985460492, 973.57897572413, 664.7393108119918, 1252.5930929049805, 1295.6270533046272, 709.6557434070465, 867.7035510603439, 673.530917697228, 1773.0318421362192, 1263.726189837023, 1704.3232128250647, 1006.0027615333212, 633.3430670419692, 1615.7422899118465, 827.2631163014008, 1111.4004294293047, 1100.8917226522785, 1402.0408965183067]

#Z:\mcdermott-group\data\Axion\2022-09-22 - DR2\DCH\Axion3A\10-10-22\PSD_Q3_Baseline
#[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#[300.70431344643583, 223.01631468457734, 224.16291971493425, 209.57515804559102, 206.78899863191626, 237.73031274534225, 229.61751354241846, 253.07091604029694, 240.49387820955624, 240.90326693229989, 234.62506624742855, 231.41098933134143, 255.46585696930705, 243.64647466782444, 238.35086155618728, 263.5713155687315, 216.3051599431929, 251.98317381648405, 194.18500306823307, 230.88936302451415, 225.58556053047568, 218.23378493936949, 288.7812187135538, 209.49737669214494, 305.1408846020174, 209.51002539812697, 234.77832736393066, 227.03587829752703, 239.54880991489475, 246.1166135070409, 228.84441448279907, 228.7072473372517, 221.3312631854884, 222.4018575084531, 231.41589414198918, 233.49462051476996, 195.6063637516116, 228.03925429838404, 226.83445621896038, 206.8331790356066, 276.03777270876344, 223.6873501795623, 270.69236641777024, 274.0288210723051, 254.88063916166524, 271.6480652430632, 229.4652873987789, 245.19582955408697, 277.176196141596, 290.71717295762886, 219.8214349756942, 209.30466553523195, 285.5391869126666, 217.20122252846033, 221.54981293709469, 233.58065949464367, 213.02420557299217, 224.76495907471318, 220.18646865091517, 245.08053272583913, 214.76843232974812, 229.85744341034487, 235.648025273138, 259.90995620298554, 232.89446362166566, 232.22162596279452, 215.57092514961457, 216.83055541877297, 238.07196161025337, 240.76512270032092, 241.19946725036942, 219.6284589059276, 210.84977350297834, 264.5264632144029, 259.9811000963801, 220.46031060744747, 222.86276819453892, 215.6048353773156]


biases2=[]
parity_rate2=[]
for i in range(len(parity_rate_data2)):
    if parity_rate_data2[i] < 5000 and parity_rate_data2[i] > 200:
        biases2.append(biases_data2[i])
        parity_rate2.append(parity_rate_data2[i])

biases1=[]
parity_rate1=[]
for i in range(len(parity_rate_data1)):
    if parity_rate_data1[i] < 5000 and parity_rate_data1[i] > 200:
        biases1.append(biases_data1[i])
        parity_rate1.append(parity_rate_data1[i])

biases=[]
parity_rate=[]
for i in range(len(parity_rate_data)):
    if parity_rate_data[i] < 5000 and parity_rate_data[i] > 200:
        biases.append(biases_data[i])
        parity_rate.append(parity_rate_data[i])


plt.title('Parity Rate vs Radiator Bias')
# plt.plot(4*times, np.abs(I + Q * 1j))
plt.semilogy([np.abs(b-0.035) for b in biases1], parity_rate1,marker='o',linestyle='None')
plt.semilogy([np.abs(b-0.035) for b in biases], parity_rate,marker='+',linestyle='None')
plt.semilogy([np.abs(b/0.49) for b in biases2], parity_rate2,marker='+',linestyle='None')

plt.grid(which='both')
plt.ylim(2e2,4.5e3)
#plt.xlim(-600,600)
# plt.xlabel('Radiator Bias (mV)')
plt.xlabel('Radiator Bias (mV)')
plt.ylabel('Parity Rate (Hz)')
plt.pause(0.1)

