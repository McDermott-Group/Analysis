from datetime import datetime
from scipy import signal
from dataChest import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('TkAgg')

# Z:\mcdermott-group\data\Axion\2022-09-22 - DR2\DCH\Axion3A\10-09-22\PSD_Q4_J1B
biases_data = [0.07643999999999991, -0.08917999999999997, -0.016659999999999904, -0.05683999999999994, -0.028419999999999918, 0.009799999999999965, 0.016660000000000123, -0.06369999999999995, -0.03233999999999992, 0.1068200000000002, 0.013719999999999968, -0.10093999999999997, -0.012739999999999902, 0.09702000000000019, -0.10387999999999999, -0.049979999999999934, 0.04899999999999993, -0.06859999999999995, 0.037239999999999954, 0.054879999999999936, -0.07839999999999996, -0.09995999999999998, 0.11367999999999988, 0.04606000000000015, -0.07545999999999996, 0.0029400000000001114, 0.04802000000000015, -0.05585999999999994, 0.0803599999999999, -0.004899999999999896, -0.01764000000000001, -0.09211999999999997, -0.020579999999999907, 0.018620000000000126, -0.06957999999999995, -0.08819999999999975, 0.06958000000000017, 0.06075999999999993, 0.07055999999999993, -0.10583999999999999, 0.005879999999999975, 0.010780000000000118, -0.08231999999999996, -0.1127, -0.006859999999999897, 0.08231999999999991, -0.1225, 0.0911400000000002, -0.11956, -0.0107799999999999, -0.014699999999999904, -0.05781999999999994, -0.041159999999999926, 0.003919999999999973, 0.014700000000000121, 0.1146600000000002, 0.040180000000000146, -0.06565999999999995, 0.07938000000000017, -0.013720000000000012, 0.03919999999999995, -0.09309999999999997, 0.08623999999999991, 0.08134000000000018, -0.018619999999999907, -0.32339999999999997, 0.019599999999999958, -0.11073999999999999, -0.07447999999999996, -0.054879999999999936, 0.05879999999999993, -0.1176, -0.3626, -0.06761999999999994, -0.05879999999999994, 0.03822000000000014, 0.06271999999999991, 0.10975999999999987, -0.050959999999999936, 0.04998000000000015, -0.02253999999999991, -0.08819999999999997, 0.03234000000000014, 0.09898000000000019, 0.06467999999999992, -0.04507999999999993, 0.12054000000000022, 0.08427999999999991, 0.050959999999999936, 0.10779999999999988, 0.07546000000000018, -0.07055999999999996, -0.08427999999999997, 0.02547999999999995, -0.09603999999999997, -0.07349999999999995, 0.01567999999999996, -0.029399999999999916, -0.05193999999999994, 0.07742000000000017, 0.06859999999999992, 0.1029000000000002, 0.11858000000000021, 0.07154000000000017, 0.1019199999999999, 0.008820000000000116, 0.09407999999999989, 0.0921199999999999, 0.08918000000000019, -0.05291999999999994, 0.011759999999999967, -0.08721999999999996, -0.20579999999999984, 0.06663999999999992, -0.11171999999999999, 0.11759999999999989, -0.10681999999999998, -0.10191999999999998, -0.06663999999999995, 0.09603999999999989, -0.07937999999999996, 0.11171999999999987, -0.01568000000000001, 0.06762000000000017, -0.09701999999999998, -0.019600000000000006, 0.05978000000000016, -0.059779999999999944, 0.04703999999999993, -0.09407999999999997, -0.026459999999999914, 0.0999599999999999, 0.00783999999999997, 0.06370000000000016, -0.16659999999999983, -0.0029399999999998937, -0.008819999999999899, -0.07251999999999996, 0.023519999999999954, 0.03626000000000014, 0.05586000000000016, 0.07251999999999992, 0.0009800000000001097, 0.057820000000000156, 0.05194000000000015, 0.03527999999999995, -0.000979999999999892, -0.30379999999999996, -0.08133999999999997, -0.08035999999999996, -0.10975999999999998, -0.037239999999999926, -0.03919999999999992, -0.001960000000000022, -0.060759999999999946, -0.009799999999999682, -0.024499999999999914, -0.11661999999999999, 0.08526000000000018, -0.07643999999999995, 0.020580000000000126, -0.11857999999999999, 0.11563999999999988, -0.03527999999999992, 0.1038799999999999, 0.11662000000000021, -0.08329999999999996, 0.07447999999999992, -0.005880000000000019, 0.12151999999999986, 0.045079999999999947, 0.02254000000000013, 0.10878000000000021, -0.036259999999999924, -0.04703999999999993, 0.031359999999999943, 0.0931000000000002, 0.03331999999999995, -0.12152, 0.05291999999999994, 0.02155999999999996, -0.04311999999999993, 0.04311999999999994, -0.11465999999999998, -0.12054, -0.046059999999999934, -0.04409999999999993, 0.1058399999999999, -0.0293999999999997, 0.06566000000000016, 0.01274000000000012, -0.03429999999999992, 0.029399999999999954, -0.18619999999999984, -0.09015999999999998, -0.04213999999999993, 0.006860000000000115, -0.04899999999999993, 0.09799999999999989, -0.06859999999999973, -0.10877999999999999, 0.02450000000000013, -0.12347999999999999, -0.2841999999999999, -0.025480000000000003, -0.34299999999999997, -0.02744, 0.028420000000000133, -0.07741999999999996, 0.1107400000000002, -0.07153999999999995, 0.034300000000000136, 0.027439999999999954, -0.09799999999999998, -0.04801999999999993, -0.09897999999999997, 0.08330000000000018, -0.023520000000000003, -0.10289999999999998, 0.1127000000000002, -0.06467999999999995, -0.11563999999999999, 0.12250000000000022, 0.1048600000000002, 0.11955999999999989, 0.1009400000000002, -2.380040609040179e-17, -0.040179999999999924, -0.06271999999999994, -0.009800000000000015, 0.0901599999999999, 0.0783999999999999, 0.07350000000000018, 0.08722000000000019, -0.1469999999999998, -0.04899999999999972, 0.026460000000000133, 0.053900000000000156, -0.03331999999999992, 0.056839999999999925, -0.2645999999999999, -0.08623999999999997, 0.0019599999999999743, 0.04214000000000014, 0.06174000000000016, -0.11367999999999999, 0.0049000000000001135, -0.10779999999999977, -0.06173999999999995, -0.021560000000000006, 0.09506000000000019, -0.030379999999999917, -0.031359999999999916, -0.03821999999999993, 0.0881999999999999, -0.05389999999999994, -0.24499999999999988, -0.08525999999999996, -0.22539999999999988, 0.030380000000000136, -0.003920000000000021, -0.09505999999999998, -0.011760000000000013, -0.10485999999999998, -0.007840000000000017, 0.017639999999999958, 0.04115999999999994, 0.044100000000000146, -0.10779999999999998, -0.09113999999999997, -0.1273999999999998]
parity_rate_data = [2127.128353530031, 2776.5871155574428, 351.88846928819265, 445.60960736495224, 681.4329923543885, 74.51581621689672, 312.44542636412075, 796.6263072175615, 582.394513471299, 2594.971380359349, 153.82305015746894, 2038.9108537349293, 116.51236736872599, 2147.4846653386157, 2266.337320706093, 1081.9771811841345, 1213.5072498091288, 1422.4422602808045, 494.2256061723046, 503.283095715458, 2142.871892892029, 2203.958592949726, 3154.3392562061745, 1062.232845399451, 1978.2406379254785, 128.37753663238252, 1386.352357824658, 522.6653286677264, 2091.473746072216, 132.19156562354908, 212.43045818576016, 2436.6790925422674, 938.2514839743665, 743.7377442690894, 1594.6279644078165, 2550.408000338136, 1935.5989625603186, 697.0928845912734, 2017.8261548437351, 2494.0470738881313, 72.51721687489947, 148.88710403436895, 1784.518909168744, 2969.8117306368367, 124.03193412072564, 1941.9790722635557, 43.61820175245118, 1963.772802707976, 3150.318075315444, 128.36177643970683, 222.04932575817682, 462.2023987451458, 697.7953244302969, 68.03430359100032, 337.94677859288083, 3159.051143245562, 542.877123128546, 1010.500962290259, 2358.8842499802518, 76.9939203471446, 510.02287187819223, 2237.525220606304, 2763.9400061508545, 1895.1737457649438, 330.71861882498723, 27.051489946381807, 648.431986237476, 2606.439419381767, 2070.0435540666917, 564.0377263478249, 528.1776912982953, 3081.886707291999, 53.323569081264026, 1299.8203676878065, 449.23882011756655, 534.1558126896157, 916.0625507167882, 2885.3329873763123, 988.5788877984924, 1112.6572080383755, 982.8020845179794, 2466.271823569172, 527.859231293913, 2219.6367286624804, 1235.899265344452, 743.84731836324, 3053.485327086674, 1199.5250838670581, 973.4144137477363, 2883.1844247596964, 2124.6228843093018, 1775.655629580513, 1769.2417391350116, 567.0162384887926, 2122.3312898667473, 2021.3094161029182, 194.1115598765784, 723.4649299694777, 867.3513549379489, 2158.997376532809, 1819.0064046499165, 2456.6509823810698, 3082.5192954235963, 2015.556315675711, 2346.3673311079665, 150.6870883500097, 2115.515472468224, 1902.8715030025917, 2341.463995471879, 793.4958093333748, 90.9027811903757, 1533.3569272378118, 26.271828855049986, 1493.1043320518775, 2954.905857912861, 3048.676270759837, 2746.2681971259594, 2181.1702876600298, 1160.8056236693699, 2330.0035888832067, 2136.7248490047446, 2950.6277661476247, 145.36301077775786, 1664.0430550669325, 2052.1863933800896, 454.797339154485, 654.2128624513975, 512.7686239132374, 1261.3413578145796, 1745.0858750999, 781.5160407963798, 2281.666799569447, 82.23033461917869, 1081.662410007663, 50.51158306649808, 118.16766054451229, 113.61232910551111, 1961.6970224628883, 732.6664591121369, 550.6046594516071, 487.63957142402916, 2064.403554106313, 118.18608622663054, 537.0902517365595, 832.7217191891475, 508.6940162663585, 114.00228562476849, 136.2987783434918, 2473.632726560675, 2068.360961507429, 2700.7283109604846, 585.5170293065842, 652.7416717273233, 59.142815756053835, 502.6990111430161, 84.76902739667716, 1068.3151434002161, 3253.0901662186284, 1991.2752704068237, 2063.625471252556, 780.3004707863581, 3141.1480678923485, 3070.750971529463, 536.7444483634155, 2521.916563699552, 3089.224672366231, 2284.187935330902, 2045.6446197990367, 70.15192866078031, 3111.037708450054, 913.207158786539, 1026.1455365534082, 2847.0133567483867, 614.7702433342992, 791.279293694854, 475.52466709188855, 2070.4340364377117, 481.86021293951706, 3178.832942049146, 643.0472482095503, 889.2820554638234, 714.3188834720172, 672.4122120275648, 2898.1792127062363, 3214.231443651788, 801.740542910314, 745.8302064862268, 2674.3443494264034, 697.5082315578911, 1382.4802229479462, 265.4283639673508, 553.0495328318945, 432.89789942461323, 31.99732157337011, 2103.222075976526, 746.3766704593824, 148.04227743553704, 970.5549196584947, 2211.10145377708, 1449.651789654271, 2613.1450963501293, 623.9105266713699, 3015.3257573393475, 27.159195359785706, 717.3513126623545, 51.06159439478452, 684.9169256503094, 519.6109461578247, 2152.9349268343803, 2960.5278176862157, 1857.298936527472, 538.3281355953493, 532.2203229900975, 2556.1017632851745, 919.4739230036695, 2138.454564662494, 1639.1689422780112, 1100.5127668959258, 2339.9885093552334, 2999.67266957163, 873.5490813576118, 3062.565410394787, 3276.100195084653, 2473.3179175429195, 3047.233922688248, 2248.63324820394, 65.04055702337799, 698.085083843244, 670.5563301442761, 91.52947691500782, 2473.7858802260953, 2150.8565312207775, 1999.6432980297793, 2553.1839137673833, 46.017248567975486, 969.2313608202123, 553.2785012231647, 537.7092867234222, 561.5696041919958, 492.63155792316763, 24.996984705360664, 1299.412998691859, 67.02621390050827, 649.5509786566446, 823.0341390668102, 3064.2101251041236, 130.83980526644052, 2594.555356855172, 615.7920104709773, 845.552305190267, 2327.2770743403853, 652.5151952069187, 549.7480621538485, 704.2465131879503, 2047.8973308802263, 708.5943768261344, 23.343838633728176, 1797.4633086192339, 21.864267922210697, 529.8185820720275, 63.08887839011054, 2107.59825344184, 82.26626527297199, 2343.1027926212414, 84.21614528452729, 409.04002263779097, 533.8997701349306, 783.5495715157648, 2742.7924259113347, 2304.111033373177, 3329.1644629296165]

biases=[]
parity_rate=[]
for i in range(len(parity_rate_data)):
    if parity_rate_data[i] < 5000 and parity_rate_data[i] > 55:
        biases.append(biases_data[i])
        parity_rate.append(parity_rate_data[i])

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

