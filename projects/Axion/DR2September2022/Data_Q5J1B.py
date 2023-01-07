from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')

# Z:\mcdermott-group\data\Axion\2022-09-22 - DR2\DCH\Axion3A\10-07-22\PSD_Q5_J2B
biases_data =[0.07154000000000017, -0.02253999999999991, -0.06565999999999995, 0.08330000000000018, -0.10779999999999977, -0.025479999999999912, 0.06370000000000016, -0.007839999999999897, 0.018620000000000126, 0.04214000000000014, -0.08231999999999996, -0.03331999999999992, -0.10583999999999999, -0.09897999999999997, 0.09702000000000019, 0.054880000000000026, -0.06859999999999995, 0.040180000000000146, -0.0107799999999999, -0.34299999999999997, -0.12152, -0.027439999999999916, 0.05292000000000002, 0.008820000000000116, -0.05683999999999994, 0.09408000000000005, -0.07545999999999996, -0.023519999999999913, -0.05291999999999994, -0.10681999999999998, -0.08819999999999997, 0.08232000000000005, -0.03821999999999993, 0.04312000000000002, 0.08722000000000019, 0.12152000000000007, -0.009799999999999682, 0.030380000000000136, 0.06566000000000016, 0.08428000000000005, 0.1029000000000002, -0.041159999999999926, -0.07643999999999995, 0.02450000000000013, -0.08035999999999996, -0.08525999999999996, -0.12347999999999999, -0.037239999999999926, -0.1469999999999998, -0.11465999999999998, -0.07349999999999995, 0.03822000000000014, -0.0293999999999997, -0.08133999999999997, 0.10780000000000006, -0.05585999999999994, -0.015679999999999906, -0.10779999999999998, -0.09505999999999998, -0.040179999999999924, -0.08329999999999996, -0.32339999999999997, -0.02155999999999991, 0.05684000000000003, -0.054879999999999936, -0.07839999999999996, -0.0029399999999998937, -0.07251999999999996, -0.10877999999999999, 0.026460000000000133, 0.1068200000000002, -0.04507999999999993, 0.044100000000000146, 0.06174000000000016, -0.04409999999999993, -0.11073999999999999, 0.034300000000000136, -0.013719999999999902, 0.10388000000000007, -0.029399999999999916, 0.04802000000000015, -0.06467999999999995, -0.059779999999999944, -0.018619999999999907, 0.01176000000000012, 0.04704000000000002, 0.07840000000000004, 0.1146600000000002, 0.009800000000000118, -0.008819999999999899, 0.09996000000000006, 0.07742000000000017, -0.10975999999999998, -0.06663999999999995, -0.000979999999999892, -0.09799999999999998, -0.06369999999999995, 0.04900000000000002, -0.1127, 0.09212000000000005, 0.017640000000000124, -0.09309999999999997, -0.001959999999999893, 0.12250000000000022, 0.02156, 0.03626000000000014, -0.16659999999999983, -0.003919999999999895, -0.10191999999999998, 0.06664000000000003, -0.10093999999999997, 0.10878000000000021, -0.03233999999999992, 0.06860000000000004, 0.12054000000000022, -0.05389999999999994, 0.050960000000000026, 0.031360000000000006, 0.027440000000000003, 0.08624000000000005, -0.05193999999999994, 0.04998000000000015, 0.07546000000000018, -0.11563999999999999, 0.007840000000000116, -0.2841999999999999, -0.22539999999999988, 0.005880000000000114, -0.30379999999999996, 0.10976000000000007, -0.06957999999999995, 0.025480000000000003, -0.07447999999999996, -0.017639999999999906, -0.03919999999999992, 0.07938000000000017, -0.11661999999999999, 0.06076000000000003, -0.016659999999999904, -0.06271999999999994, 0.09506000000000019, 0.0049000000000001135, 0.1107400000000002, -0.10289999999999998, 0.09016000000000005, 0.0911400000000002, 0.05194000000000015, 0.019600000000000124, 0.08526000000000018, -0.11367999999999999, -0.07153999999999995, 0.1009400000000002, -0.11857999999999999, 1.0880185641326534e-16, -0.03429999999999992, -0.004899999999999896, -0.09211999999999997, -0.11956, -0.07055999999999996, 0.09898000000000019, 0.09604000000000006, -0.12054, 0.03724000000000001, -0.11171999999999999, -0.06761999999999994, -0.028419999999999918, -0.09113999999999997, -0.09015999999999998, 0.10584000000000006, -0.012739999999999902, 0.02254000000000013, 0.06272000000000003, 0.0009800000000001097, -0.036259999999999924, -0.3626, 0.020580000000000126, -0.06859999999999973, -0.060759999999999946, -0.050959999999999936, -0.08721999999999996, -0.1225, -0.08917999999999997, 0.07644000000000004, 0.01568000000000012, 0.05586000000000016, 0.05880000000000003, 0.06762000000000017, 0.11956000000000008, -0.10387999999999999, -0.09407999999999997, -0.0117599999999999, -0.04703999999999993, 0.08036000000000004, 0.1048600000000002, 0.07056000000000004, 0.08134000000000018, -0.020579999999999907, -0.05781999999999994, -0.09603999999999997, 0.006860000000000115, 0.0931000000000002, 0.08918000000000019, 0.003920000000000113, -0.03527999999999992, 0.041160000000000016, 0.057820000000000156, -0.030379999999999917, -0.08427999999999997, -0.046059999999999934, 0.07252000000000004, -0.04213999999999993, 0.07448000000000003, -0.07741999999999996, -0.2645999999999999, 0.11172000000000007, -0.04311999999999993, 0.053900000000000156, 0.028420000000000133, -0.09995999999999998, 0.04606000000000015, 0.11662000000000021, -0.20579999999999984, 0.013720000000000121, 0.03528000000000001, -0.04899999999999972, 0.11564000000000008, -0.026459999999999914, 0.08820000000000006, 0.06958000000000017, -0.014699999999999904, 0.010780000000000118, 0.11368000000000007, 0.10192000000000007, -0.08819999999999975, 0.03234000000000014, 0.02352, 0.07350000000000018, 0.014700000000000121, -0.049979999999999934, -0.024499999999999914, 0.1127000000000002, -0.01959999999999991, 0.016660000000000123, 0.03332000000000001, -0.031359999999999916, 0.01274000000000012, -0.006859999999999897, -0.05879999999999994, 0.0019600000000001105, -0.10485999999999998, -0.04801999999999993, -0.1273999999999998, -0.18619999999999984, 0.029400000000000006, -0.07937999999999996, -0.09701999999999998, 0.11760000000000008, -0.1176, -0.24499999999999988, 0.09800000000000006, -0.04899999999999993, -0.009799999999999899, 0.03920000000000001, 0.045080000000000016, 0.05978000000000016, -0.005879999999999897, -0.06173999999999995, 0.0029400000000001114, -0.08623999999999997, 0.06468000000000003, 0.11858000000000021]
parity_rate_data =[2491.061992559313, 825.0720254172667, 686.5081976178908, 1252.0115786046783, 1842.6683282309475, 593.995590884203, 777.9347563161738, 113.14409767225861, 1874.9715884025854, 610.1126552455887, 908.444379616252, 347.1738850700219, 1932.8379413565485, 1459.616191281773, 1662.4304615884996, 536.1440354166211, 1657.0365430317095, 410.9846750717367, 116.03282317677534, 75.28092889293598, 2511.4669048152305, 582.4563737830117, 768.6696269786835, 123.82965133672575, 435.3172077541508, 1545.8270723593903, 2729.419490260484, 697.8882474431684, 1051.4764531900735, 1808.5513866244373, 1141.635137012393, 981.9606733272758, 421.95300237283413, 717.113941178461, 1696.3507298270322, 2696.394575723373, 122.56122441567089, 325.70739132447045, 1308.0960239546266, 1221.916813080704, 1825.760728067476, 547.6786250954287, 2686.948734645015, 511.3719142125426, 3832.1259206252753, 1040.7055641462402, 2625.460589224721, 374.7829336345078, 3655.9358351107676, 1891.4103430213986, 13.819657374116481, 427.7475850398675, 617.8914382607381, 3932.114575518189, 2043.2195464916833, 456.8063724535228, 159.67130105744204, 1644.4801108762915, 1559.4823143403637, 493.2815684463864, 1044.1708752060983, 19.192288393410916, 634.8789628866598, 418.56415979190393, 662.0988843080753, 3233.590653003041, 108.9589618277132, 2504.540168028832, 2163.989393909654, 551.422354008018, 2085.6377672622602, 857.2702468881353, 883.6238411289188, 529.7783698685583, 697.7709224071895, 1867.962540950567, 299.11904369560096, 118.58785132919913, 1895.7495474695204, 625.0554998814129, 1563.5662161355988, 714.815874997299, 367.34546079987007, 241.93577182763653, 121.46086174529209, 1532.6205807643878, 3190.6304555718216, 2383.1016382434345, 113.23275507597165, 114.66292873526666, 1723.2051948203218, 3306.35987747654, 2022.8544867578198, 1025.7460223130472, 107.75685787136227, 1359.2649315124609, 501.8940680413171, 1807.1401003974302, 1944.6500209321562, 1550.63781209681, 269.06358611111307, 1200.391130838633, 128.8296650344565, 2824.600978672975, 525.243220497888, 354.1881760251612, 4417.13930189272, 111.06404958566314, 1580.1261369785666, 1585.218025082401, 1644.6821866108291, 1895.9322959928852, 370.8507754461587, 2187.535681880896, 2305.0625489951367, 781.4962546177168, 1288.7108784768395, 274.6440433848496, 417.8360925302556, 1419.7583055618747, 1011.9845941701321, 1011.1630896121496, 3218.753357291856, 2188.112050681457, 116.11768111093859, 29.268443050576188, 23.764927554710713, 108.9656129171312, 21.89984236596645, 2126.6965044011604, 1659.152974970206, 494.3585938676125, 2993.529354288621, 235.19580347988116, 482.4987298881584, 3904.5076235892434, 2291.649111709259, 538.8693128666338, 193.5502060064451, 509.6829072231813, 1248.9485606444434, 117.74943219379134, 2058.486676852946, 1704.317833100232, 1296.1757389161558, 1197.4504556198476, 950.2631932138702, 551.9507837900098, 943.8963997050745, 2027.7597742948435, 2414.859448680575, 1728.85882421814, 56.64643929439163, 122.38036550114049, 294.8599283916196, 102.17551636014646, 1264.027344115126, 2661.0013959877147, 1768.986282384366, 1560.1306074423362, 1482.0879191247668, 2450.176969873062, 443.46077443865937, 2044.5884159501227, 1023.1547305108933, 651.634010500106, 1478.288130646904, 1333.7843261180994, 1955.9982555343527, 114.67449696518143, 640.1166370565655, 685.1694254358576, 106.22565254513579, 317.67253581281443, 20.827109806408277, 576.966204818759, 1457.3284286993376, 393.42818284991404, 1379.3598613743306, 953.8720811264136, 2603.1429972907404, 1438.1117233266007, 2937.9511845918723, 208.87484981429088, 450.9614214268553, 401.6290626464447, 1595.151935574528, 2663.295561308439, 1802.839959622373, 1329.7104936498133, 114.90035093164828, 1055.5479343719953, 1799.096017488334, 1868.5934608146797, 2779.2198093175557, 1183.176159150945, 353.2962761339148, 406.60322500250425, 1523.5416163647799, 113.00079026184693, 1670.02427369949, 1474.6386577565504, 102.96139146464374, 344.37403780117114, 477.655614516401, 442.6357149042778, 476.0999479732537, 843.3654243924066, 953.762717351596, 2894.3242477970643, 488.6320590714686, 2767.376747574757, 2881.2857441739257, 19.868035077842638, 2273.381821585018, 575.4200545835909, 642.0338475765983, 409.9840114948347, 1554.4907780981107, 1241.6948863098385, 2301.539670832867, 20.601719719325583, 142.3963044601741, 311.4313297914518, 1313.7655986360342, 2450.2557222134965, 561.9508524297755, 1477.6179251188767, 2349.814289378256, 135.52924313129088, 121.16318583602992, 2381.9178171246767, 1916.7346316296062, 1201.8842524861736, 273.27858145267015, 1306.3615250668186, 2688.774061767672, 217.64451526608678, 1336.0776747659304, 677.7350128526642, 2163.0799451378607, 336.1217947510936, 275.59651713271546, 276.4207216773793, 450.377664256645, 140.87011966951644, 104.88334798738542, 441.87748013989085, 115.8644476702147, 1798.2965265480505, 1049.9183157343457, 2629.5861683863623, 21.213889649695258, 347.87895436307434, 3276.942802926226, 1236.5695176851093, 2423.673566179544, 2339.6120890027387, 24.793489926340428, 1513.0883728299802, 1301.0846140676724, 109.92242043227473, 453.0295813106455, 1097.9541767485573, 494.91028212100224, 105.23753738844434, 385.718713936573, 99.71365711054898, 1244.673657403614, 1089.8450294545232, 2695.0553624004046]

# biases_data1  = [-0.06467999999999995, 0.04214000000000014, -0.08035999999999996, -0.029399999999999916, -0.05879999999999994, -0.07545999999999996, 0.0009800000000001097, -0.09015999999999998, 0.05586000000000016, -0.09799999999999998, 0.05194000000000015, -0.10387999999999999, -0.12152, -0.04213999999999993, -0.10681999999999998, -0.06957999999999995, 0.005880000000000114, 0.06076000000000016, -0.025479999999999912, -0.10877999999999999, -0.03919999999999992, -0.005879999999999897, 0.07056000000000016, -0.050959999999999936, -0.05193999999999994, -0.1225, 0.007840000000000116, -0.06761999999999994, 0.06958000000000017, 0.014700000000000121, 0.026460000000000133, 0.07350000000000018, -0.04311999999999993, 0.049000000000000155, 0.05684000000000016, -0.08623999999999997, -0.08133999999999997, -0.03429999999999992, -0.03331999999999992, 0.0029400000000001114, 0.05978000000000016, -0.09603999999999997, 0.020580000000000126, 0.08624000000000018, 0.039200000000000144, -0.000979999999999892, -0.08329999999999996, 0.03234000000000014, 0.07644000000000017, -0.0107799999999999, -0.026459999999999914, 1.0880185641326534e-16, -0.017639999999999906, -0.09407999999999997, -0.03527999999999992, 0.05488000000000016, -0.07643999999999995, -0.10583999999999999, -0.04409999999999993, 0.04704000000000015, -0.10779999999999998, -0.12054, -0.008819999999999899, -0.09995999999999998, -0.006859999999999897, -0.07937999999999996, 0.021560000000000128, -0.06565999999999995, -0.003919999999999895, -0.08917999999999997, 0.03626000000000014, -0.020579999999999907, -0.08427999999999997, 0.016660000000000123, 0.08918000000000019, -0.023519999999999913, -0.040179999999999924, 0.07252000000000017, 0.04606000000000015, 0.07742000000000017, 0.03528000000000014, 0.07938000000000017, 0.01568000000000012, 0.03332000000000014, 0.04116000000000015, -0.05585999999999994, 0.06762000000000017, 0.02548000000000013, 0.04508000000000015, -0.02253999999999991, -0.06663999999999995, -0.09701999999999998, 0.08428000000000019, -0.04507999999999993, 0.05880000000000016, -0.004899999999999896, -0.08525999999999996, 0.07448000000000017, -0.036259999999999924, 0.07546000000000018, 0.06370000000000016, -0.030379999999999917, -0.11367999999999999, -0.018619999999999907, 0.02450000000000013, -0.07839999999999996, -0.05683999999999994, -0.03233999999999992, -0.028419999999999918, -0.014699999999999904, 0.034300000000000136, 0.08820000000000018, -0.054879999999999936, -0.016659999999999904, 0.08722000000000019, 0.052920000000000154, 0.08232000000000018, -0.001959999999999893, 0.02744000000000013, 0.006860000000000115, 0.023520000000000128, 0.013720000000000121, -0.09897999999999997, -0.01959999999999991, 0.028420000000000133, -0.009799999999999899, -0.11857999999999999, 0.029400000000000134, -0.049979999999999934, 0.06468000000000017, 0.040180000000000146, 0.07154000000000017, 0.017640000000000124, 0.06860000000000017, 0.04998000000000015, -0.05389999999999994, 0.08330000000000018, 0.04802000000000015, 0.010780000000000118, 0.019600000000000124, -0.10093999999999997, 0.07840000000000018, -0.027439999999999916, -0.09211999999999997, -0.07251999999999996, -0.007839999999999897, 0.03724000000000014, 0.02254000000000013, 0.01176000000000012, -0.07153999999999995, -0.07349999999999995, -0.11956, -0.07741999999999996, 0.03136000000000014, -0.05781999999999994, -0.11465999999999998, 0.003920000000000113, -0.09309999999999997, 0.06174000000000016, -0.060759999999999946, -0.046059999999999934, 0.06664000000000017, 0.08526000000000018, 0.030380000000000136, -0.07447999999999996, -0.09505999999999998, 0.008820000000000116, 0.044100000000000146, 0.05096000000000015, -0.11073999999999999, 0.057820000000000156, -0.08721999999999996, 0.09016000000000018, -0.04703999999999993, -0.08231999999999996, -0.0029399999999998937, -0.06859999999999995, -0.05291999999999994, -0.10485999999999998, -0.03821999999999993, -0.11171999999999999, 0.009800000000000118, -0.02155999999999991, 0.018620000000000126, 0.08134000000000018, 0.03822000000000014, -0.09113999999999997, -0.04899999999999993, -0.024499999999999914, -0.10975999999999998, -0.012739999999999902, 0.01274000000000012, -0.0117599999999999, 0.08036000000000018, -0.11563999999999999, -0.10289999999999998, 0.06566000000000016, -0.041159999999999926, -0.037239999999999926, -0.013719999999999902, -0.015679999999999906, -0.06369999999999995, -0.07055999999999996, 0.0049000000000001135, -0.1127, -0.10191999999999998, -0.059779999999999944, -0.08819999999999997, -0.06271999999999994, -0.031359999999999916, -0.1225, -0.11661999999999999, -0.1176, 0.053900000000000156, 0.043120000000000144, 0.0019600000000001105, -0.04801999999999993, -0.06173999999999995, 0.06272000000000016]
# parity_rate_data1 = [938.5556449858478, 1758.5033608655074, 1093.0344469189674, 1730.9935228970428, 1630.5021623511145, 856.2658893944944, 839.2457953265761, 1411.437114577877, 1638.2883539536253, 1664.4692502059265, 2260.321629984574, 1867.3331486019956, 3247.009094356949, 2016.110835442547, 2773.2297576419637, 903.4959234734224, 964.3162423529557, 1297.2641966550343, 890.2385701101097, 3185.171061553501, 1907.865395949799, 641.6247682382782, 796.497527576919, 891.6560029667644, 925.2492524461347, 2239.911130118741, 823.8003018871647, 839.0132732160185, 840.1604754743407, 701.5685951095347, 895.9497353866697, 942.5860386971954, 1526.8371012919938, 1207.1860696732556, 1394.6684918820977, 1445.2683887757498, 1156.4625296656902, 1093.122722196025, 2876.4149782419595, 807.5202883027669, 2084.287066656266, 1679.1927230992499, 1141.7486488355048, 1282.6297068193835, 3068.8322868195055, 546.6253691164197, 1451.977510668609, 1327.1811507125292, 1149.7268796523788, 674.5059116831927, 889.602472778845, 468.87342079851265, 698.8612825891645, 1600.6651581243875, 1157.0281788340953, 1247.0731641501384, 884.2418087147346, 2649.103404836825, 1232.947753702982, 1023.8659742765427, 2766.712517698653, 4198.328081513512, 665.745662299386, 1803.614751862763, 641.140643997905, 1083.4054472750836, 1014.8085514834622, 872.6051347738752, 672.3876470362512, 1316.7187224587976, 1316.4120997958653, 1304.3398262439425, 1510.3203664966622, 844.3635047858558, 1365.324484882487, 1411.5888847052418, 2130.7344812759156, 880.1258292826329, 765.9323015417169, 949.0288261549933, 1831.3841188003444, 1095.1832410047653, 787.7685536375136, 1477.7404350575687, 1407.2822347113374, 1281.130253055326, 881.7254764623344, 894.873315580371, 958.4020618530379, 1190.292152734327, 859.6144475437807, 1875.4975116469507, 1512.0003638814276, 952.2448034902319, 3228.9245892163003, 640.2224670900114, 1446.853789410875, 1052.976103817752, 2442.5312240339067, 1084.7179056570628, 1094.2753276436092, 1910.2157339432877, 1872.7677137461046, 759.5348367459899, 895.2514255973966, 960.7292968271372, 1376.3726643322775, 2205.4595135505315, 913.7551880129332, 715.6588471463426, 3813.1060974587776, 1414.6687459504913, 1249.5648051649252, 652.4593065217907, 1240.1755113069626, 1813.1524415185809, 1760.5074389390732, 622.6161214070562, 896.4168637123813, 1177.2349139835922, 1036.6746179623872, 602.2167774904927, 1728.0416892935013, 1223.985349904342, 944.6348248331883, 681.34727225151, 2218.8287928029613, 1075.3616910327441, 906.0909653137318, 1239.6603427081952, 2153.5453861250276, 863.6767156514187, 991.3517209931446, 850.2702287634781, 1468.7546420507338, 1140.3693008041573, 1191.2875857456477, 1239.602288077468, 750.4680745730466, 1123.8961277867502, 1714.8403112849699, 962.0849995917699, 862.341336021298, 1475.204513818187, 767.3178046249692, 668.6024104181081, 2790.450217630967, 988.7840101542608, 671.816250702069, 882.8042330514111, 920.9315069319122, 3964.405842823174, 938.0828800983156, 1215.9610760186572, 1408.4901374373821, 1942.5674130573545, 506.94226316495735, 1440.6312656072246, 1253.085063563647, 1315.4258516242035, 1000.2230336872507, 935.4392921371896, 1419.8516057895918, 1144.2381054255518, 807.9628999130516, 1737.1567376462426, 759.38189220445, 1749.0975897789874, 1866.4547905908225, 3507.4034009178154, 1386.9292468254187, 1585.6291719759258, 1270.8124864591537, 969.804829878776, 1354.556042200875, 668.8445274161818, 877.1860974506208, 1238.714488157356, 2173.4928111869895, 1540.5942756780369, 3323.170360425849, 735.7054883196053, 1125.5650641137986, 1122.7499405242131, 1538.394764212286, 2512.3054855234645, 1396.523486419362, 845.5215898111505, 957.0313688875416, 3300.8359285243628, 641.6915381393518, 694.4767213464037, 602.4738311125427, 1146.4118756541118, 2051.6668028043846, 1786.0015993220995, 1057.6290142991866, 2200.171204956437, 1109.6440714625971, 743.1539857409566, 940.8155604292488, 991.3476665052452, 954.3732871503205, 653.3173302452354, 2668.0748903697986, 1757.5702783429401, 1586.2098957865835, 1780.1191803143092, 1298.568630399104, 2091.752824832879, 3048.851358714067, 3383.3979928587514, 2743.389116266757, 1176.556954244013, 1843.8388168379217, 522.3124207823796, 925.1896506139426, 1472.6257692475872, 1887.6788889402974]



biases=[]
parity_rate=[]
for i in range(len(parity_rate_data)):
    if parity_rate_data[i] < 5000 and parity_rate_data[i] > 50:
        biases.append(biases_data[i])
        parity_rate.append(parity_rate_data[i])

# biases1 = []
# parity_rate1 = []
# for i in range(len(parity_rate_data1)):
#     if parity_rate_data1[i] < 5000 and parity_rate_data1[i] > 200:
#         biases1.append(biases_data1[i])
#         parity_rate1.append(parity_rate_data1[i])


plt.title('Parity Rate vs Radiator Bias')
# plt.plot(4*times, np.abs(I + Q * 1j))
# plt.semilogy([np.abs(b/0.49) for b in biases1], parity_rate1,marker='o',linestyle='None')
plt.semilogy([np.abs(b/0.49) for b in biases], parity_rate,marker='+',linestyle='None')

plt.grid(which='both')
plt.ylim(0.5e2,4.5e3)
plt.xlim(0,0.4)
# plt.xlabel('Radiator Bias (mV)')
plt.xlabel('Radiator Bias (mV)')
plt.ylabel('Parity Rate (Hz)')
plt.pause(0.1)

