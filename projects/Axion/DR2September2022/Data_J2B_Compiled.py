from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')

parity_rate_data = {}
biases_data = {}
for qb in ['Q1J2','Q3J2','Q4J2','Q5J2']:
    parity_rate_data[qb] = []
    biases_data[qb] = []

## With antenna JJ
# Q1
biases_data['Q1J2'].extend([-0.06467999999999995, 0.04214000000000014, -0.08035999999999996, -0.029399999999999916, -0.05879999999999994, -0.07545999999999996, 0.0009800000000001097, -0.09015999999999998, 0.05586000000000016, -0.09799999999999998, 0.05194000000000015, -0.10387999999999999, -0.12152, -0.04213999999999993, -0.10681999999999998, -0.06957999999999995, 0.005880000000000114, 0.06076000000000016, -0.025479999999999912, -0.10877999999999999, -0.03919999999999992, -0.005879999999999897, 0.07056000000000016, -0.050959999999999936, -0.05193999999999994, -0.1225, 0.007840000000000116, -0.06761999999999994, 0.06958000000000017, 0.014700000000000121, 0.026460000000000133, 0.07350000000000018, -0.04311999999999993, 0.049000000000000155, 0.05684000000000016, -0.08623999999999997, -0.08133999999999997, -0.03429999999999992, -0.03331999999999992, 0.0029400000000001114, 0.05978000000000016, -0.09603999999999997, 0.020580000000000126, 0.08624000000000018, 0.039200000000000144, -0.000979999999999892, -0.08329999999999996, 0.03234000000000014, 0.07644000000000017, -0.0107799999999999, -0.026459999999999914, 1.0880185641326534e-16, -0.017639999999999906, -0.09407999999999997, -0.03527999999999992, 0.05488000000000016, -0.07643999999999995, -0.10583999999999999, -0.04409999999999993, 0.04704000000000015, -0.10779999999999998, -0.12054, -0.008819999999999899, -0.09995999999999998, -0.006859999999999897, -0.07937999999999996, 0.021560000000000128, -0.06565999999999995, -0.003919999999999895, -0.08917999999999997, 0.03626000000000014, -0.020579999999999907, -0.08427999999999997, 0.016660000000000123, 0.08918000000000019, -0.023519999999999913, -0.040179999999999924, 0.07252000000000017, 0.04606000000000015, 0.07742000000000017, 0.03528000000000014, 0.07938000000000017, 0.01568000000000012, 0.03332000000000014, 0.04116000000000015, -0.05585999999999994, 0.06762000000000017, 0.02548000000000013, 0.04508000000000015, -0.02253999999999991, -0.06663999999999995, -0.09701999999999998, 0.08428000000000019, -0.04507999999999993, 0.05880000000000016, -0.004899999999999896, -0.08525999999999996, 0.07448000000000017, -0.036259999999999924, 0.07546000000000018, 0.06370000000000016, -0.030379999999999917, -0.11367999999999999, -0.018619999999999907, 0.02450000000000013, -0.07839999999999996, -0.05683999999999994, -0.03233999999999992, -0.028419999999999918, -0.014699999999999904, 0.034300000000000136, 0.08820000000000018, -0.054879999999999936, -0.016659999999999904, 0.08722000000000019, 0.052920000000000154, 0.08232000000000018, -0.001959999999999893, 0.02744000000000013, 0.006860000000000115, 0.023520000000000128, 0.013720000000000121, -0.09897999999999997, -0.01959999999999991, 0.028420000000000133, -0.009799999999999899, -0.11857999999999999, 0.029400000000000134, -0.049979999999999934, 0.06468000000000017, 0.040180000000000146, 0.07154000000000017, 0.017640000000000124, 0.06860000000000017, 0.04998000000000015, -0.05389999999999994, 0.08330000000000018, 0.04802000000000015, 0.010780000000000118, 0.019600000000000124, -0.10093999999999997, 0.07840000000000018, -0.027439999999999916, -0.09211999999999997, -0.07251999999999996, -0.007839999999999897, 0.03724000000000014, 0.02254000000000013, 0.01176000000000012, -0.07153999999999995, -0.07349999999999995, -0.11956, -0.07741999999999996, 0.03136000000000014, -0.05781999999999994, -0.11465999999999998, 0.003920000000000113, -0.09309999999999997, 0.06174000000000016, -0.060759999999999946, -0.046059999999999934, 0.06664000000000017, 0.08526000000000018, 0.030380000000000136, -0.07447999999999996, -0.09505999999999998, 0.008820000000000116, 0.044100000000000146, 0.05096000000000015, -0.11073999999999999, 0.057820000000000156, -0.08721999999999996, 0.09016000000000018, -0.04703999999999993, -0.08231999999999996, -0.0029399999999998937, -0.06859999999999995, -0.05291999999999994, -0.10485999999999998, -0.03821999999999993, -0.11171999999999999, 0.009800000000000118, -0.02155999999999991, 0.018620000000000126, 0.08134000000000018, 0.03822000000000014, -0.09113999999999997, -0.04899999999999993, -0.024499999999999914, -0.10975999999999998, -0.012739999999999902, 0.01274000000000012, -0.0117599999999999, 0.08036000000000018, -0.11563999999999999, -0.10289999999999998, 0.06566000000000016, -0.041159999999999926, -0.037239999999999926, -0.013719999999999902, -0.015679999999999906, -0.06369999999999995, -0.07055999999999996, 0.0049000000000001135, -0.1127, -0.10191999999999998, -0.059779999999999944, -0.08819999999999997, -0.06271999999999994, -0.031359999999999916, -0.1225, -0.11661999999999999, -0.1176, 0.053900000000000156, 0.043120000000000144, 0.0019600000000001105, -0.04801999999999993, -0.06173999999999995, 0.06272000000000016])
parity_rate_data['Q1J2'].extend([938.5556449858478, 1758.5033608655074, 1093.0344469189674, 1730.9935228970428, 1630.5021623511145, 856.2658893944944, 839.2457953265761, 1411.437114577877, 1638.2883539536253, 1664.4692502059265, 2260.321629984574, 1867.3331486019956, 3247.009094356949, 2016.110835442547, 2773.2297576419637, 903.4959234734224, 964.3162423529557, 1297.2641966550343, 890.2385701101097, 3185.171061553501, 1907.865395949799, 641.6247682382782, 796.497527576919, 891.6560029667644, 925.2492524461347, 2239.911130118741, 823.8003018871647, 839.0132732160185, 840.1604754743407, 701.5685951095347, 895.9497353866697, 942.5860386971954, 1526.8371012919938, 1207.1860696732556, 1394.6684918820977, 1445.2683887757498, 1156.4625296656902, 1093.122722196025, 2876.4149782419595, 807.5202883027669, 2084.287066656266, 1679.1927230992499, 1141.7486488355048, 1282.6297068193835, 3068.8322868195055, 546.6253691164197, 1451.977510668609, 1327.1811507125292, 1149.7268796523788, 674.5059116831927, 889.602472778845, 468.87342079851265, 698.8612825891645, 1600.6651581243875, 1157.0281788340953, 1247.0731641501384, 884.2418087147346, 2649.103404836825, 1232.947753702982, 1023.8659742765427, 2766.712517698653, 4198.328081513512, 665.745662299386, 1803.614751862763, 641.140643997905, 1083.4054472750836, 1014.8085514834622, 872.6051347738752, 672.3876470362512, 1316.7187224587976, 1316.4120997958653, 1304.3398262439425, 1510.3203664966622, 844.3635047858558, 1365.324484882487, 1411.5888847052418, 2130.7344812759156, 880.1258292826329, 765.9323015417169, 949.0288261549933, 1831.3841188003444, 1095.1832410047653, 787.7685536375136, 1477.7404350575687, 1407.2822347113374, 1281.130253055326, 881.7254764623344, 894.873315580371, 958.4020618530379, 1190.292152734327, 859.6144475437807, 1875.4975116469507, 1512.0003638814276, 952.2448034902319, 3228.9245892163003, 640.2224670900114, 1446.853789410875, 1052.976103817752, 2442.5312240339067, 1084.7179056570628, 1094.2753276436092, 1910.2157339432877, 1872.7677137461046, 759.5348367459899, 895.2514255973966, 960.7292968271372, 1376.3726643322775, 2205.4595135505315, 913.7551880129332, 715.6588471463426, 3813.1060974587776, 1414.6687459504913, 1249.5648051649252, 652.4593065217907, 1240.1755113069626, 1813.1524415185809, 1760.5074389390732, 622.6161214070562, 896.4168637123813, 1177.2349139835922, 1036.6746179623872, 602.2167774904927, 1728.0416892935013, 1223.985349904342, 944.6348248331883, 681.34727225151, 2218.8287928029613, 1075.3616910327441, 906.0909653137318, 1239.6603427081952, 2153.5453861250276, 863.6767156514187, 991.3517209931446, 850.2702287634781, 1468.7546420507338, 1140.3693008041573, 1191.2875857456477, 1239.602288077468, 750.4680745730466, 1123.8961277867502, 1714.8403112849699, 962.0849995917699, 862.341336021298, 1475.204513818187, 767.3178046249692, 668.6024104181081, 2790.450217630967, 988.7840101542608, 671.816250702069, 882.8042330514111, 920.9315069319122, 3964.405842823174, 938.0828800983156, 1215.9610760186572, 1408.4901374373821, 1942.5674130573545, 506.94226316495735, 1440.6312656072246, 1253.085063563647, 1315.4258516242035, 1000.2230336872507, 935.4392921371896, 1419.8516057895918, 1144.2381054255518, 807.9628999130516, 1737.1567376462426, 759.38189220445, 1749.0975897789874, 1866.4547905908225, 3507.4034009178154, 1386.9292468254187, 1585.6291719759258, 1270.8124864591537, 969.804829878776, 1354.556042200875, 668.8445274161818, 877.1860974506208, 1238.714488157356, 2173.4928111869895, 1540.5942756780369, 3323.170360425849, 735.7054883196053, 1125.5650641137986, 1122.7499405242131, 1538.394764212286, 2512.3054855234645, 1396.523486419362, 845.5215898111505, 957.0313688875416, 3300.8359285243628, 641.6915381393518, 694.4767213464037, 602.4738311125427, 1146.4118756541118, 2051.6668028043846, 1786.0015993220995, 1057.6290142991866, 2200.171204956437, 1109.6440714625971, 743.1539857409566, 940.8155604292488, 991.3476665052452, 954.3732871503205, 653.3173302452354, 2668.0748903697986, 1757.5702783429401, 1586.2098957865835, 1780.1191803143092, 1298.568630399104, 2091.752824832879, 3048.851358714067, 3383.3979928587514, 2743.389116266757, 1176.556954244013, 1843.8388168379217, 522.3124207823796, 925.1896506139426, 1472.6257692475872, 1887.6788889402974])
# Q4
biases_data['Q4J2'].extend([-0.09603999999999997, -0.07741999999999996, 0.04998000000000007, 0.07056000000000016, -0.013719999999999902, 0.06860000000000017, -0.11661999999999999, -0.11563999999999999, -0.11073999999999999, -0.04311999999999993, -0.09995999999999998, 0.012740000000000038, -0.20579999999999996, 0.02450000000000005, 0.08036000000000018, -0.12054, 0.1078000000000002, -0.06173999999999995, 0.12250000000000014, -0.026459999999999997, -0.08525999999999996, -0.06565999999999995, 0.022540000000000046, 0.020580000000000046, -0.08721999999999996, -0.049979999999999934, 0.0970200000000001, 0.05096000000000015, -0.08819999999999997, -0.05193999999999994, -0.006859999999999979, -0.07153999999999995, -0.08231999999999996, 0.02744000000000013, -0.09015999999999998, -0.03821999999999993, -0.05879999999999994, -0.2646, -0.2841999999999999, -0.01861999999999999, 0.00294000000000003, 0.03822000000000006, 0.05390000000000007, 0.021560000000000128, -0.027439999999999916, -0.08427999999999997, -0.023519999999999913, 0.11270000000000012, -0.05389999999999994, 0.05684000000000016, 0.06762000000000008, 0.06272000000000016, -0.004899999999999977, 0.06566000000000008, -0.11857999999999999, -0.08329999999999996, -0.08133999999999997, -0.07643999999999995, 0.07154000000000009, 0.09506000000000012, 0.0960400000000002, 0.01568000000000012, -0.05291999999999994, 0.06958000000000009, -0.046059999999999934, -0.040179999999999924, 0.07938000000000009, -0.012739999999999984, -0.18619999999999995, 0.11172000000000021, 0.03528000000000014, -0.03527999999999992, -0.16659999999999991, 0.019600000000000124, -0.017639999999999906, 0.09310000000000011, -0.04213999999999993, 0.11074000000000013, -0.09897999999999997, 0.0891800000000001, -0.0117599999999999, -0.245, 0.08134000000000009, 0.11858000000000013, -0.034300000000000004, 0.11662000000000013, 0.10290000000000012, -0.08035999999999996, 0.09212000000000019, 0.11564000000000021, -0.06467999999999995, 0.02548000000000013, -0.07055999999999996, -0.050959999999999936, -0.0009799999999999737, 0.017640000000000124, -0.10583999999999999, 0.007840000000000116, 0.1176000000000002, -0.3626, -0.022539999999999994, 0.049000000000000155, 0.10682000000000012, -0.10387999999999999, -0.02155999999999991, 0.08624000000000018, 0.052920000000000154, -0.09799999999999998, -0.1469999999999999, -0.09701999999999998, 0.04802000000000007, -0.028419999999999997, 0.06174000000000008, 0.11956000000000021, -0.03234, 0.05194000000000007, 0.004900000000000031, 0.0833000000000001, -0.07545999999999996, 0.05488000000000016, -0.060759999999999946, -0.04801999999999993, 0.07252000000000017, 0.09898000000000011, 0.039200000000000144, 0.06664000000000017, 0.043120000000000144, -0.05585999999999994, -0.09309999999999997, 0.055860000000000076, 1.0880185641326534e-16, -0.10877999999999999, 0.010780000000000036, -0.06761999999999994, -0.07349999999999995, -0.03919999999999992, 0.05782000000000008, -0.01959999999999991, 0.005880000000000114, -0.1127, -0.10093999999999997, -0.015679999999999906, 0.046060000000000066, -0.06369999999999995, -0.005879999999999897, -0.05781999999999994, 0.07350000000000009, 0.0774200000000001, -0.11956, -0.029399999999999916, 0.1058400000000002, 0.008820000000000034, -0.07937999999999996, -0.03038, 0.10094000000000011, 0.059780000000000076, 0.0754600000000001, 0.03724000000000014, -0.1176, 0.023520000000000128, 0.0009800000000000282, -0.10975999999999998, -0.10191999999999998, 0.10192000000000019, 0.07448000000000017, 0.032340000000000056, 0.0999600000000002, -0.08623999999999997, 0.0852600000000001, 0.11368000000000021, -0.02057999999999999, 0.12152000000000021, -0.2254, 0.1038800000000002, 0.02646000000000005, -0.024499999999999994, -0.007839999999999897, -0.06663999999999995, -0.07447999999999996, -0.07839999999999996, 0.09016000000000018, -0.009799999999999899, -0.003919999999999895, -0.04507999999999993, 0.08428000000000019, 0.006860000000000033, -0.11367999999999999, 0.06076000000000016, -0.014699999999999986, 0.016660000000000043, 0.04410000000000006, -0.10485999999999998, 0.1097600000000002, -0.04899999999999993, 0.11466000000000012, 0.08232000000000018, -0.04703999999999993, 0.05880000000000016, -0.12347999999999999, -0.07251999999999996, -0.06859999999999995, -0.03331999999999992, -0.016659999999999987, -0.001959999999999893, 0.0872200000000001, 0.06370000000000008, 0.07840000000000018, 0.04018000000000006, -0.041159999999999926, -0.002939999999999975, 0.013720000000000121, -0.00881999999999998, 0.10878000000000013, -0.010779999999999982, -0.09113999999999997, -0.10681999999999998, 0.009800000000000118, -0.05683999999999994, -0.059779999999999944, -0.32339999999999997, 0.09114000000000011, 0.036260000000000056, 0.0980000000000002, 0.01470000000000004, 0.10486000000000012, -0.10289999999999998, -0.037239999999999926, -0.06271999999999994, -0.10779999999999998, 0.0019600000000001105, -0.34299999999999997, 0.12054000000000013, -0.09505999999999998, -0.12152, -0.025479999999999912, -0.09407999999999997, 0.04116000000000015, -0.11171999999999999, -0.036259999999999924, -0.031359999999999916, 0.03136000000000014, 0.018620000000000043, -0.1225, 0.06468000000000017, -0.08917999999999997, -0.06957999999999995, 0.003920000000000113, 0.03332000000000014, -0.09211999999999997, -0.10779999999999988, 0.042140000000000066, 0.04704000000000015, -0.1273999999999999, 0.09408000000000019, 0.01176000000000012, 0.04508000000000015, -0.11465999999999998, 0.08820000000000018, -0.30379999999999996, 0.030380000000000053, -0.054879999999999936, 0.07644000000000017, 0.029400000000000134, 0.028420000000000053, 0.03430000000000006, -0.04409999999999993])
parity_rate_data['Q4J2'].extend([1199.6664565123556, 433.64697848299915, 522.2045129937384, 331.9159802287976, 114.65126126862377, 338.813063565473, 2175.95795863804, 2165.899348496834, 1637.021548168034, 829.6392729033487, 1261.6626382645707, 91.99978866498451, 27.502607411282167, 383.4007450255458, 744.3216430597158, 2232.106468349297, 1552.094871648386, 755.2914224806196, 2153.193339688825, 471.24379612506795, 727.6412525186121, 460.9547194429838, 397.75946413906536, 398.0145733106211, 844.1964397632894, 511.72240120291826, 1245.190739206556, 521.6179221221801, 850.6974324603818, 543.0490763084381, 71.09400377536574, 329.1308428655885, 682.0454115143978, 493.0912957076494, 856.7336638089519, 989.3287705186996, 990.6678187817083, 24.363720073318213, 61.66505557182625, 318.39340585781554, 73.90969467638854, 1146.983909149044, 694.9265645554823, 451.2994511162882, 482.65404211627174, 704.3496365796422, 337.8468669296772, 1985.9608521486477, 632.6168242664452, 924.2964324567515, 350.17754033361183, 634.5553225671524, 69.55291336086766, 418.85750358686323, 2183.114372417507, 655.7407015916694, 689.8638203374625, 406.0684115907181, 330.3910964796653, 1139.761369587299, 1241.471651606855, 194.77717242452033, 577.1279610808689, 325.85420766298375, 554.090137708125, 983.4101366181566, 628.4510552029213, 87.11285858947676, 3417.8899170365107, 1981.661979404863, 755.3604680376219, 664.619533608338, 3567.2930791980184, 265.3847422914757, 289.59243025315556, 965.2301998094254, 952.4664339604368, 1889.9694855005212, 1243.6197064261198, 839.0708932684654, 89.01094855489738, 31.079983658172285, 670.7326673963195, 2299.475509671417, 637.7051738680449, 2219.119085519484, 1428.1975566497294, 665.9117490849386, 926.1896845872633, 2147.360993243212, 497.75635345577734, 482.8930984331809, 332.2227329367967, 502.7821880773372, 67.20558286463968, 230.4092320076787, 1493.4500293501635, 75.88278521475614, 2322.2650010204707, 26.28301507409139, 468.137601088983, 515.6520920905034, 1517.841661556814, 1439.017117039653, 415.84111675043295, 814.8004360037388, 617.7998298665921, 1202.2490565418518, 2838.3502962190273, 1242.906820106023, 536.9446419848053, 473.0882019623335, 687.0246921668323, 2328.600378234576, 664.1554350431462, 578.9286057925385, 75.1750296242799, 697.3638990513339, 400.7519448961872, 793.772449667653, 819.1857226855535, 518.072079434302, 355.70112901880304, 1314.0487457414888, 1066.4147446947998, 386.52936538931925, 676.0334962840633, 800.0770932658379, 880.8499972893112, 866.6451320648707, 66.95500223766874, 1504.0872505878451, 89.70386795731572, 374.3274453811288, 358.2119436443163, 1287.820017605301, 1031.2360051595429, 407.0682905769141, 68.43765274096317, 2049.725396627705, 1273.6792950390131, 202.89066125181083, 551.5005362811864, 559.5058886653753, 75.32107350937184, 977.7347034038105, 375.7469013335493, 473.0483719776412, 2281.0527735712126, 563.301107307344, 1522.4291495795371, 77.4229533764707, 537.72214601089, 594.1341654057546, 1328.5967514279428, 842.0146846956656, 389.1470855720121, 986.9618366814957, 2215.022792124647, 360.33641309614535, 72.7561566252549, 1543.1130865712962, 1272.3127594243567, 1362.4383994127827, 405.31162922706073, 660.3013292448294, 1317.5178197769033, 779.2593613474883, 796.1754340128589, 2024.6084943026142, 493.9952047327824, 2326.00445923998, 27.094704639306865, 1463.4468589037226, 484.4717114380383, 348.13945145534257, 77.6279794553426, 408.4250933877005, 394.8405351809085, 475.545791812874, 843.7493495446536, 73.99054619100073, 69.52685721074805, 613.6100253566499, 753.9879326190475, 74.91398134506075, 2040.1408590833082, 799.6843080252171, 140.85864087263323, 213.70111523225407, 626.9811809337241, 1461.1922221545199, 1571.4745826684011, 505.46060476490595, 2048.2264889431267, 709.977150911144, 531.2968523954995, 989.4112834245022, 2271.1765569412432, 339.66947213899505, 329.71357036312264, 671.2175588486527, 230.15462903384505, 66.9115824337707, 863.563816599716, 532.2614005225603, 504.5411318818797, 963.512729516543, 1000.9029946333061, 64.95947896482046, 125.35064097280203, 76.51406885650762, 1499.7759692525913, 74.67376748763202, 892.9580987587715, 1488.4337523682013, 70.58461422296539, 906.1973960670819, 966.3335058093552, 57.99173400318439, 906.6876460774124, 831.9550585037668, 1219.5163051650452, 141.55352220007072, 1479.3840046235837, 1432.1121768473101, 877.5409255585334, 637.6909222837797, 1538.135708024238, 65.73109626163534, 27.173665034663312, 2348.7676593897245, 1065.3801715480286, 2264.629325864625, 400.6469974420539, 972.8093978493301, 1037.01500988296, 1898.5380453353205, 785.7327993426735, 647.6793513063578, 689.6332444981834, 276.05883066321616, 2076.1472165970795, 482.29762176394985, 833.8652692075099, 330.21532559286214, 64.62129553052894, 625.3560913321313, 914.2732229582215, 1479.7590931695192, 874.6781498470522, 521.1002738831592, 2188.5156540174667, 1056.2969760270691, 85.94266556912626, 562.0870638972115, 2024.9794928602303, 854.4052699506793, 53.388892684038346, 622.5027803014399, 733.8108723271167, 418.749377116645, 600.4114555076497, 548.4426311705606, 663.5603764352204, 703.3715143102123])

## Reference (without antenna JJ)
# Q3
biases_data['Q3J2'].extend([0.49*(x-0.035) for x in [-0.3, -0.295, -0.29, -0.285, -0.27999999999999997, -0.27499999999999997, -0.26999999999999996, -0.26499999999999996, -0.25999999999999995, -0.25499999999999995, -0.24999999999999994, -0.24499999999999994, -0.23999999999999994, -0.23499999999999993, -0.22999999999999993, -0.22499999999999992, -0.21999999999999992, -0.2149999999999999, -0.2099999999999999, -0.2049999999999999, -0.1999999999999999, -0.1949999999999999, -0.1899999999999999, -0.1849999999999999, -0.17999999999999988, -0.17499999999999988, -0.16999999999999987, -0.16499999999999987, -0.15999999999999986, -0.15499999999999986, -0.14999999999999986, -0.14499999999999985, -0.13999999999999985, -0.13499999999999984, -0.12999999999999984, -0.12499999999999983, -0.11999999999999983, -0.11499999999999982, -0.10999999999999982, -0.10499999999999982, -0.09999999999999981, -0.0949999999999998, -0.0899999999999998, -0.0849999999999998, -0.0799999999999998, -0.07499999999999979, -0.06999999999999978, -0.06499999999999978, -0.059999999999999776, -0.05499999999999977, -0.04999999999999977, -0.04499999999999976, -0.03999999999999976, -0.034999999999999754, -0.02999999999999975, -0.024999999999999745, -0.01999999999999974, -0.014999999999999736, -0.009999999999999731, -0.004999999999999727, 2.7755575615628914e-16, 0.005000000000000282, 0.010000000000000286, 0.015000000000000291, 0.020000000000000295, 0.0250000000000003, 0.030000000000000304, 0.03500000000000031, 0.04000000000000031, 0.04500000000000032, 0.05000000000000032, 0.055000000000000326, 0.06000000000000033, 0.06500000000000034, 0.07000000000000034, 0.07500000000000034, 0.08000000000000035, 0.08500000000000035, 0.09000000000000036, 0.09500000000000036, 0.10000000000000037, 0.10500000000000037, 0.11000000000000038, 0.11500000000000038, 0.12000000000000038, 0.1250000000000004, 0.1300000000000004, 0.1350000000000004, 0.1400000000000004, 0.1450000000000004, 0.1500000000000004, 0.15500000000000042, 0.16000000000000042, 0.16500000000000042, 0.17000000000000043, 0.17500000000000043, 0.18000000000000044, 0.18500000000000044, 0.19000000000000045, 0.19500000000000045, 0.20000000000000046, 0.20500000000000046, 0.21000000000000046, 0.21500000000000047, 0.22000000000000047, 0.22500000000000048, 0.23000000000000048, 0.2350000000000005, 0.2400000000000005, 0.2450000000000005, 0.2500000000000005, 0.2550000000000005, 0.2600000000000005, 0.2650000000000005, 0.2700000000000005, 0.2750000000000005, 0.2800000000000005, 0.28500000000000053, 0.29000000000000054, 0.29500000000000054, 0.30000000000000054]])
parity_rate_data['Q3J2'].extend([34.4645025401347, 40.52850710748375, 38.33124490311449, 31.61707126834793, 29.370404211570573, 3243.1006900091547, 3356.7004670448323, 25.227754787612827, 29.285861947081703, 3147.2624045120424, 25.989471607356542, 3141.6226196906064, 32.107368323788535, 3182.490671436952, 23.320984935339382, 53.58135170372689, 2956.5639140986264, 3162.7214544970234, 2977.614496723117, 2733.0460405472204, 29.836661200868992, 46.02710067983333, 35.63199141243225, 37.94435876444115, 43.60728327187585, 2303.909671969774, 2335.1783689362164, 37.41237358175952, 244.3873436412472, 1678.425483735762, 1592.607867736693, 1914.9447112965108, 1829.1820658215122, 1538.6662078788274, 1512.716141038692, 1356.716769833392, 1018.8415218783516, 879.5013291431165, 1072.222563365139, 924.8291913931255, 952.1631425714805, 1039.9048460777306, 1207.5838943734366, 1370.0317328233662, 1327.3755481673634, 1362.3406280767558, 1118.5396417754548, 1179.619527516466, 1065.3911031458854, 1264.4284662247749, 1484.4319208039067, 1843.8665502535946, 1551.1586516672762, 1331.0984679675373, 1628.5157303495807, 1632.9228484262687, 1282.8405040831963, 1225.6635576848712, 1296.4017040393835, 1313.3302281349845, 1033.595935784898, 856.0790232838619, 912.8898233902504, 898.6091160237672, 705.9632727477058, 614.0140266443442, 775.5042888868717, 810.5665543144465, 963.5397743670318, 710.2249093246352, 867.9678095029708, 1001.3583486832212, 700.1633244305517, 664.7747642715253, 1006.5081702619409, 814.3526150048681, 951.8325749652516, 948.158702663741, 941.9596221702759, 1015.7915376129669, 1124.62856765482, 1434.8165118836496, 1723.907229206477, 1431.3398554495127, 1231.0159947407722, 1000.226196960157, 1121.775740896926, 1219.0629533880322, 1170.5355844257006, 1060.838486234636, 1388.5501780321429, 1528.9350429757328, 1320.538195987509, 922.1366941832885, 796.5040509422189, 680.689340007343, 878.090894343147, 897.0809371967308, 890.3795490578901, 963.8452836286162, 1210.6235787772541, 1339.2569174276869, 1368.2900917531995, 1649.4661840201052, 1534.4190246694968, 1577.3242426725953, 1650.7920603808468, 1960.9390494080997, 2107.489718320197, 1984.9288633648086, 2034.9018446466753, 2029.6251969244406, 2200.0309296494315, 2341.1013704336174, 2449.5650210842473, 2693.4203640980563, 2869.7318728300706, 2895.5448450120457, 2758.7796625334104, 2981.9635926898513, 2929.9464772363585])
# Q5
biases_data['Q5J2'].extend([-0.11073999999999999, -0.08329999999999996, -0.07545999999999996, -0.06369999999999995, -0.06761999999999994, -0.10093999999999997, -0.06957999999999995, -0.08721999999999996, -0.08917999999999997, -0.10681999999999998, -0.1225, -0.07741999999999996, -0.08133999999999997, -0.09505999999999998, -0.11661999999999999, -0.10877999999999999, -0.12054, -0.07937999999999996, -0.10289999999999998, -0.07349999999999995, -0.09701999999999998, -0.09113999999999997, -0.11465999999999998, -0.1127, -0.11857999999999999, -0.08525999999999996, -0.10485999999999998, -0.09897999999999997, -0.06565999999999995, -0.07153999999999995, -0.09309999999999997, 0.07154000000000003, -0.09113999999999997, 0.09114, -0.026459999999999914, 0.02842, 0.0029400000000001114, -0.30379999999999996, -0.06957999999999995, 0.10289999999999999, 0.09898, 0.08918, -0.12054, -0.016659999999999904, 0.07546, -0.08133999999999997, -0.05193999999999994, -0.05389999999999994, 0.05586, 0.03626000000000001, 0.018619999999999998, 0.06566000000000002, -0.006859999999999897, 0.03234, 0.034300000000000004, -0.09897999999999997, 0.10878000000000002, -0.05781999999999994, -0.0107799999999999, -0.05585999999999994, -0.04801999999999993, -0.046059999999999934, -0.07545999999999996, -0.10779999999999977, 0.05586000000000002, 0.1225, 0.02450000000000013, -0.018619999999999907, 0.06566, 0.020580000000000126, 0.00294, -0.16659999999999983, -0.016659999999999904, 0.01274000000000012, 0.02646, 0.04606000000000001, -0.06369999999999995, -0.11857999999999999, 0.10486000000000001, 0.08722000000000005, 0.10682000000000001, 0.11661999999999999, -0.08525999999999996, -0.30379999999999996, 0.008820000000000116, -0.06173999999999995, -0.024499999999999914, -0.1225, 0.038220000000000004, 0.016659999999999998, 0.04410000000000001, -0.028419999999999918, -0.04899999999999972, -0.008819999999999899, -0.08819999999999975, -0.059779999999999944, -0.00294, 0.018620000000000126, -0.000979999999999892, -0.1273999999999998, -0.02253999999999991, -0.10877999999999999, 0.08330000000000004, -0.012739999999999902, -0.2645999999999999, -0.11465999999999998, 0.05978000000000003, 0.12054, -0.10485999999999998, 0.08525999999999999, -0.08133999999999997, 0.0931, 0.061739999999999996, 0.00098, -0.1273999999999998, -0.10093999999999997, 0.006860000000000115, -0.10681999999999998, -0.03821999999999993, 0.06762000000000003, 0.08721999999999999, 0.09702000000000001, -0.030379999999999917, -0.04213999999999993, 0.07742, 0.048020000000000014, -0.2841999999999999, -0.008819999999999899, -0.05585999999999994, -0.00098, -0.09897999999999997, 0.00686, -0.11661999999999999, 0.022539999999999998, 0.05978000000000001, 0.06370000000000003, 0.06957999999999999, 0.09702000000000001, -0.20579999999999984, -0.0293999999999997, -0.11465999999999998, -0.34299999999999997, 0.07546000000000004, 0.03626000000000001, -0.10485999999999998, -0.08819999999999975, -0.06173999999999995, -0.07937999999999996, -0.08329999999999996, -0.08721999999999996, -0.024499999999999914, 0.10486, -0.18619999999999984, -0.03429999999999992, -0.049979999999999934, -0.09701999999999998, -0.03233999999999992, -0.014699999999999904, 0.04018000000000001, 0.0049000000000001135, 0.10682, 0.05194000000000001, -0.07153999999999995, 0.028420000000000133, 0.11074, -0.0107799999999999, -0.1469999999999998, -0.3626, 0.07350000000000004, 0.07938000000000005, 0.057820000000000024, -0.11661999999999999, 0.061740000000000024, 0.11465999999999998, -0.018619999999999907, 0.07938, -0.03429999999999992, -0.11073999999999999, -0.05389999999999994, 0.06958000000000003, 0.044100000000000014, 0.026460000000000133, -0.11073999999999999, 0.024499999999999997, 0.07153999999999999, -0.06859999999999973, -0.04899999999999972, -0.07349999999999995, -0.09113999999999997, -0.036259999999999924, -0.09505999999999998, 0.0735, -0.009799999999999682, -0.07741999999999996, -0.07545999999999996, 0.05194000000000002, -0.09309999999999997, 0.01078, -0.10779999999999977, -0.3626, 0.0009800000000001097, 0.09114, -0.08721999999999996, 0.08134000000000004, -0.10289999999999998, -0.049979999999999934, -0.11857999999999999, -0.04801999999999993, 0.08329999999999999, -0.10093999999999997, -0.06565999999999995, 0.03234, -0.08917999999999997, -0.08917999999999997, -0.0029399999999998937, -0.020579999999999907, -0.006859999999999897, -0.24499999999999988, -0.014699999999999904, -0.09505999999999998, -0.028419999999999918, -0.04409999999999993, 0.04802, 0.0049, -0.18619999999999984, 0.0147, -0.059779999999999944, 0.02254000000000013, -0.026459999999999914, -0.1225, 0.016660000000000123, -0.07937999999999996, -0.036259999999999924, 0.030380000000000004, 0.012740000000000001, 0.1127, -0.04409999999999993, 0.049980000000000004, 0.0637, -0.0293999999999997, 0.08526000000000004, 0.06761999999999999, -0.07349999999999995, 0.034300000000000004, -0.32339999999999997, -0.04213999999999993, -0.24499999999999988, -0.030379999999999917, -0.10681999999999998, 0.12054, -0.009799999999999682, 0.10878, -0.09701999999999998, -0.20579999999999984, 0.1225, -0.03233999999999992, -0.004899999999999896, 0.11857999999999999, -0.046059999999999934, 0.11857999999999999, 0.04998000000000002, -0.020579999999999907, 0.046060000000000004, -0.09309999999999997, -0.040179999999999924, 0.053899999999999997, -0.06369999999999995, -0.07741999999999996, -0.012739999999999902, 0.02058, -0.1127, 0.1029, -0.06565999999999995, -0.05781999999999994, 0.10094, -0.040179999999999924, -0.34299999999999997, -0.22539999999999988, 0.05390000000000002, 0.0931, 0.008820000000000001, 0.09506, -0.05193999999999994, -0.06761999999999994, -0.2841999999999999, 0.07742000000000004, -0.10877999999999999, -0.1469999999999998, 0.04018, -0.16659999999999983, 0.09897999999999998, -0.07153999999999995, 0.010780000000000118, 0.03822, 0.042140000000000004, 0.08134, 0.04214000000000001, 0.05782, 0.09506, -0.06761999999999994, 0.014700000000000121, -0.06859999999999973, -0.22539999999999988, -0.0049, -0.12054, -0.06957999999999995, -0.10289999999999998, 0.10093999999999999, -0.08525999999999996, -0.32339999999999997, -0.02253999999999991, -0.03821999999999993, -0.1127, -0.08329999999999996, -0.2645999999999999])
parity_rate_data['Q5J2'].extend([1296.0949896723541, 603.1753014167397, 289.0274813941571, 556.3412919256104, 339.89811472789137, 1171.625479731675, 295.29978406962465, 566.047143858423, 647.4672198767306, 1207.0765718738419, 2021.8370838659812, 331.59399393154405, 579.3597605459657, 851.7678120205107, 1617.7274043200503, 1195.5395233395943, 1835.4351894444412, 463.05222828864265, 1111.3348913224008, 287.36058997567574, 916.7226982617294, 698.9474523828826, 1617.7471356438891, 1550.7972535443248, 1785.3360208154772, 636.6840349677367, 1124.8064499969662, 1113.8343802066422, 448.15708540677093, 274.97275406391276, 764.7730860907068, 301.2064158680899, 728.7479578153573, 834.0425163078781, 272.42677940775064, 333.3819216610802, 134.65274313812233, 40.337949849102614, 330.026117044851, 1178.9286657899556, 1271.652444203869, 669.9184846512665, 1904.598000947346, 200.95577137882506, 319.37273671991596, 590.4009014160372, 457.0736062726403, 714.1858082228993, 820.5300693101277, 805.919045376194, 248.03132881695484, 428.9621964593179, 103.34163725263357, 464.8329089121441, 599.4788504997451, 1089.0141434318612, 1449.9429442006779, 1122.9348043287757, 119.83551971182024, 742.359140714928, 368.63322440563303, 393.4743491840681, 334.03106012293057, 1283.5054991887941, 864.1058547486439, 1925.0173978889393, 295.3527539850211, 345.7432022777067, 414.8329160702682, 342.6790059886254, 107.81104670746495, 3614.2592150605924, 205.47542842849148, 146.31714955004892, 268.24719706439237, 352.94495280060926, 608.5043299458871, 1637.9460045879312, 1410.03318974038, 604.0986689387545, 1404.068185591649, 1779.622825597508, 679.7616021140406, 50.672086687315, 140.130834312075, 766.8648279547497, 251.71802647554046, 1867.3274923945405, 1013.9309692762138, 168.68784930578633, 480.51752799336094, 306.5298334829083, 350.59149904815087, 111.6462766437163, 612.0466253458856, 911.0990659595719, 103.65249715279198, 272.51181089401206, 117.33181056433973, 2204.0084163596916, 294.14222254789644, 1195.0585645353362, 654.7608048596778, 108.46082018938648, 348.6687912311336, 1439.558615056488, 878.1600122842966, 1972.289981506595, 1081.3060861752558, 694.9562274640627, 582.5501025296247, 879.0213165214406, 760.4003661705012, 115.82523118030143, 2307.5252344711016, 1076.9466150625724, 145.29677811278134, 1209.0243470359733, 925.9981131916562, 351.55113116496767, 637.9681038687057, 1074.2005073360892, 392.8330290261832, 884.5133830875992, 417.60212421233666, 379.6426316231477, 42.54765416833615, 104.77493786654023, 741.313008499574, 120.45983324873899, 1102.0445831858788, 95.4269281166595, 1720.044585430511, 293.9238671926485, 904.4977840114482, 594.3960069812181, 278.0687943264725, 1078.1741780711684, 23.365396704908463, 340.89127931726625, 1671.2259166875608, 173.02479490057806, 287.62254498178675, 755.9870983495601, 1064.0983673235285, 650.1637025424568, 740.0671589578631, 475.10953701893254, 639.2617185344732, 583.0425545069279, 265.33071564637834, 1253.1206072053194, 4024.4720616143386, 558.2825599652275, 432.1124364494782, 866.6606623828337, 443.0153739734168, 141.99903729811533, 946.2588865885381, 144.27979358775266, 1214.0186748380024, 489.4364881030731, 316.8724739765569, 374.8164217831096, 1605.4712523712114, 101.7811238605513, 3208.0847110115624, 126.34984436842608, 308.9411895429868, 454.7799505646329, 1123.9538097864686, 1614.429364789214, 747.0333959121103, 1716.5125330314816, 319.8342380342225, 475.2549146268194, 595.6229431911241, 1299.8590824869095, 724.9244965681239, 308.3615714249165, 498.5971336368923, 290.497777896631, 1263.476165245244, 244.90656039286023, 296.7387886817229, 313.01101999558233, 372.23856732253273, 308.05436442630526, 725.6219272053899, 723.6868142527308, 853.9750530369082, 316.73752797077697, 94.77964444260537, 352.9129402928623, 318.714576624736, 500.95258908770717, 733.9043807316965, 40.71689272691627, 1310.697898907902, 157.23265825543817, 128.46462024402032, 779.8187247245851, 580.7580984245628, 662.8987355817446, 1085.5093773681226, 398.48189123340063, 1688.6102010364407, 354.05169624037677, 662.3951536966389, 1078.044480410303, 479.64504897223657, 465.6055305352598, 634.3842873328123, 620.7521395468124, 105.702850691213, 409.95598927053453, 103.34491393707764, 217.74726439773133, 142.57101770469097, 845.9771133180672, 320.0432342637166, 587.6481283690565, 361.91024869898985, 95.59975117670635, 20.25298210849727, 119.84376352822187, 917.5171523706364, 336.043094600031, 288.64396299794504, 1808.0532316340925, 212.82601428769584, 472.61615649708835, 674.8239811820572, 428.63869473789197, 107.04609944176126, 1806.504690288807, 546.4425414017022, 351.6071793988336, 583.3176584952068, 349.27534111446283, 671.3958588741406, 320.06769583324035, 286.8167280286189, 581.1125942190923, 155.85233262987362, 876.2786396964239, 508.41629052637677, 377.1596514212859, 1180.2943314400718, 2092.5262257628488, 84.51059827742581, 1321.4019054117473, 891.501094606805, 24.637373742239834, 2115.2213632682488, 458.14347669796933, 104.28545779191946, 1657.7406333012611, 381.6945052924948, 1922.4424208696187, 423.3789225602916, 414.6869743263313, 352.64507542924747, 748.7061309597651, 897.174724956409, 694.5779593235202, 613.6474598990468, 360.24806257910836, 110.572503198508, 308.3764022758823, 1557.0409972678322, 1299.7746010448573, 479.4659995110507, 1106.0960069139192, 1173.366665068961, 18.080228021923848, 492.80404530512345, 42.17254139475663, 742.9143041675944, 870.5516328312038, 90.53944844131158, 963.2432774867682, 436.8147618780767, 351.39965156905555, 977.6527575521253, 416.963229310713, 1184.2487412901448, 41.72699365515707, 992.9203257216993, 18.889331737548687, 1194.7693602592512, 291.7785630071074, 137.24946890948442, 968.2380246056596, 764.3645197025756, 653.6815422387093, 818.0846197953199, 1186.75023516431, 912.0046487965619, 371.41101384785856, 178.71904645722347, 320.7605424152988, 28.988367423101074, 95.43184273813809, 1770.8423756509926, 336.6110754148159, 1077.4857106828283, 1170.7893770821581, 697.9030333287466, 178.01241453440792, 306.6335875450627, 933.5488164456918, 1592.4958099295657, 597.211431074419, 47.19189711760283])


color = iter(cm.rainbow(np.linspace(0,1,2)))
for qb in ['Q5J2','Q4J2']:##['Q1J2','Q3J2','Q4J2','Q5J2']:
    biases=[]
    parity_rate=[]
    for i in range(len(parity_rate_data[qb])):
        if parity_rate_data[qb][i] < 5000 and parity_rate_data[qb][i] > 50:
            biases.append(biases_data[qb][i])
            parity_rate.append(parity_rate_data[qb][i])
    plt.title('Parity Rate vs Radiator Bias')
    # plt.plot(4*times, np.abs(I + Q * 1j))
    print(qb)
    x=400*(400/271)
    plt.semilogy([x*np.abs(b)/0.058 for b in biases], parity_rate, marker='o' if qb == 'Q5J2' else '+', color=next(color), linestyle='None', label='Qubit+Reference' if qb == 'Q5J2' else 'Qubit+Detector')
    #plt.semilogy([np.abs(b) for b in biases], parity_rate,marker='+',color=next(color),linestyle='None',label=qb)

    plt.grid(which='both')
    # plt.ylim(0.4e3,8e3)
    plt.xlim(0,x*0.12/0.058)
    plt.legend()
    # plt.xlabel('Radiator Bias (mV)')
    plt.xlabel('Radiator Frequency (GHz)')
    plt.ylabel('Parity Rate (Hz)')
    plt.pause(0.1)

