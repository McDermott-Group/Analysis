"""### To analyze the P1 vs Duration
### Z:\mcdermott-group\data\Antenna\SUXmon2\Liu\VitoChip2\10-27-21

"""
import noiselib
import matplotlib.pyplot as plt
import numpy as np

Q2_22mDAC = np.array([
 [0.0, 0.0122, 0.00287], [10.0, 0.02519, 0.00337], [20.0, 0.02875, 0.00487],
 [30.0, 0.03476, 0.00477], [40.0, 0.03667, 0.00533], [50.0, 0.0361, 0.00177],
 [60.0, 0.03242, 0.0047], [70.0, 0.03687, 0.00368], [80.0, 0.03684, 0.0061],
 [90.0, 0.03855, 0.00798], [100.0, 0.03933, 0.00504], [110.0, 0.03784, 0.00478],
 [120.0, 0.03684, 0.00649], [130.0, 0.03532, 0.00645], [140.0, 0.03853, 0.00912],
 [150.0, 0.0379, 0.008], [160.0, 0.03647, 0.00322], [170.0, 0.03666, 0.00495],
 [180.0, 0.04022, 0.00487], [190.0, 0.03891, 0.00352], [200.0, 0.03794, 0.00536],
 [210.0, 0.04196, 0.00871], [220.0, 0.03469, 0.00482], [230.0, 0.0402, 0.00397],
 [240.0, 0.03813, 0.00607], [250.0, 0.03668, 0.00579], [260.0, 0.03569, 0.00393],
 [270.0, 0.04148, 0.00792], [280.0, 0.03921, 0.00663], [290.0, 0.03875, 0.00587],
 [300.0, 0.03825, 0.00394], [310.0, 0.03973, 0.00388], [320.0, 0.04137, 0.00559],
 [330.0, 0.03862, 0.00585], [340.0, 0.03858, 0.005], [350.0, 0.04005, 0.00601],
 [360.0, 0.03808, 0.00549], [370.0, 0.03906, 0.00636], [380.0, 0.03813, 0.00588],
 [390.0, 0.03675, 0.00668]
])

Q2_25mDAC = np.array([
 [0.0, 0.01439, 0.00469], [10.0, 0.01639, 0.00602], [20.0, 0.01847, 0.00588],
 [30.0, 0.02055, 0.00753], [40.0, 0.01865, 0.00478], [50.0, 0.02261, 0.00849],
 [60.0, 0.01922, 0.00378], [70.0, 0.02308, 0.00603], [80.0, 0.02086, 0.00509],
 [90.0, 0.01982, 0.00363], [100.0, 0.02046, 0.00575], [110.0, 0.01909, 0.00248],
 [120.0, 0.02019, 0.00793], [130.0, 0.0195, 0.0053], [140.0, 0.02188, 0.00366],
 [150.0, 0.02228, 0.00344], [160.0, 0.02008, 0.00358], [170.0, 0.02059, 0.0055],
 [180.0, 0.02039, 0.00298], [190.0, 0.01884, 0.00425], [200.0, 0.02158, 0.00743],
 [210.0, 0.02238, 0.00384], [220.0, 0.02114, 0.007], [230.0, 0.02233, 0.00192],
 [240.0, 0.02334, 0.00598], [250.0, 0.02233, 0.00843], [260.0, 0.02053, 0.00432],
 [270.0, 0.02037, 0.00377], [280.0, 0.02005, 0.00561], [290.0, 0.02067, 0.00362],
 [300.0, 0.01915, 0.00551], [310.0, 0.02224, 0.00899], [320.0, 0.02225, 0.00575],
 [330.0, 0.01912, 0.00475], [340.0, 0.02248, 0.0053], [350.0, 0.02266, 0.00847],
 [360.0, 0.02106, 0.00649], [370.0, 0.02147, 0.00252], [380.0, 0.02128, 0.0056],
 [390.0, 0.02424, 0.00767]
])

Q2_32mDAC = np.array([
 [0.0, 0.01384, 0.00329], [10.0, 0.03286, 0.00656], [20.0, 0.04454, 0.00498],
 [30.0, 0.04587, 0.00646], [40.0, 0.05221, 0.00551], [50.0, 0.05113, 0.00768],
 [60.0, 0.05032, 0.00584], [70.0, 0.05224, 0.00565], [80.0, 0.05083, 0.00696],
 [90.0, 0.05293, 0.00622], [100.0, 0.05458, 0.00641], [110.0, 0.05671, 0.00736],
 [120.0, 0.04962, 0.00528], [130.0, 0.05484, 0.00922], [140.0, 0.04851, 0.00542],
 [150.0, 0.05541, 0.00419], [160.0, 0.05223, 0.0063], [170.0, 0.04959, 0.00511],
 [180.0, 0.05611, 0.00544], [190.0, 0.05189, 0.00493], [200.0, 0.04955, 0.00551],
 [210.0, 0.05184, 0.00478], [220.0, 0.05379, 0.0065], [230.0, 0.05419, 0.00734],
 [240.0, 0.05213, 0.00676], [250.0, 0.05505, 0.00705], [260.0, 0.05168, 0.00228],
 [270.0, 0.05039, 0.00562], [280.0, 0.05242, 0.0047], [290.0, 0.05114, 0.00638],
 [300.0, 0.05362, 0.00746], [310.0, 0.05449, 0.00734], [320.0, 0.05132, 0.00672],
 [330.0, 0.05716, 0.00584], [340.0, 0.05751, 0.00577], [350.0, 0.05913, 0.00351],
 [360.0, 0.05466, 0.0074], [370.0, 0.05373, 0.00574], [380.0, 0.05461, 0.00652],
 [390.0, 0.05499, 0.00927]
])

Q2_56mDAC = np.array([
 [0.0, 0.01441, 0.00269], [10.0, 0.11557, 0.00637], [20.0, 0.0947, 0.00895],
 [30.0, 0.07446, 0.00721], [40.0, 0.05819, 0.00524], [50.0, 0.0499, 0.0056],
 [60.0, 0.04295, 0.00463], [70.0, 0.04286, 0.00542], [80.0, 0.04138, 0.0049],
 [90.0, 0.04291, 0.00595], [100.0, 0.04312, 0.00649], [110.0, 0.0398, 0.00463],
 [120.0, 0.03926, 0.00645], [130.0, 0.03812, 0.00441], [140.0, 0.03788, 0.00516],
 [150.0, 0.04128, 0.00521], [160.0, 0.04028, 0.00426], [170.0, 0.0394, 0.00524],
 [180.0, 0.04156, 0.0055], [190.0, 0.04218, 0.00555], [200.0, 0.03907, 0.00357],
 [210.0, 0.03857, 0.00275], [220.0, 0.04322, 0.00526], [230.0, 0.03948, 0.00398],
 [240.0, 0.04275, 0.00575], [250.0, 0.03769, 0.00596], [260.0, 0.04173, 0.00564],
 [270.0, 0.04004, 0.00388], [280.0, 0.04331, 0.00449], [290.0, 0.04033, 0.00445],
 [300.0, 0.04413, 0.00659], [310.0, 0.04304, 0.00312], [320.0, 0.04613, 0.00592],
 [330.0, 0.04159, 0.00648], [340.0, 0.04025, 0.00677], [350.0, 0.04104, 0.00377],
 [360.0, 0.0416, 0.00335], [370.0, 0.04406, 0.00374], [380.0, 0.0426, 0.00394],
 [390.0, 0.04658, 0.00628]
])

Q2_300mDAC = np.array([
 [0.0, 0.01197, 0.00419], [10.0, 0.08682, 0.00616], [20.0, 0.0321, 0.00394],
 [30.0, 0.02488, 0.00387], [40.0, 0.02203, 0.00279], [50.0, 0.02162, 0.00465],
 [60.0, 0.02248, 0.00398], [70.0, 0.02189, 0.00223], [80.0, 0.02281, 0.0036],
 [90.0, 0.02278, 0.00275], [100.0, 0.02396, 0.00378], [110.0, 0.02393, 0.00212],
 [120.0, 0.02409, 0.00448], [130.0, 0.02512, 0.00405], [140.0, 0.02446, 0.00351],
 [150.0, 0.02546, 0.00343], [160.0, 0.02634, 0.00462], [170.0, 0.02707, 0.00445],
 [180.0, 0.02708, 0.00329], [190.0, 0.02583, 0.00205], [200.0, 0.02902, 0.00326],
 [210.0, 0.02955, 0.0046], [220.0, 0.02977, 0.00329], [230.0, 0.02974, 0.00334],
 [240.0, 0.03001, 0.00526], [250.0, 0.0283, 0.00165], [260.0, 0.03218, 0.00234],
 [270.0, 0.02929, 0.0048], [280.0, 0.03182, 0.00496], [290.0, 0.03196, 0.00358],
 [300.0, 0.03201, 0.00392], [310.0, 0.03208, 0.00453], [320.0, 0.03134, 0.00443],
 [330.0, 0.03247, 0.00261], [340.0, 0.03198, 0.00326], [350.0, 0.0309, 0.00275],
 [360.0, 0.0323, 0.00527], [370.0, 0.03312, 0.00298], [380.0, 0.03594, 0.00541],
 [390.0, 0.03708, 0.003]
])

Q2_450mDAC = np.array([
 [0.0, 0.01162, 0.00442], [10.0, 0.06866, 0.00687], [20.0, 0.03876, 0.00293],
 [30.0, 0.03693, 0.00477], [40.0, 0.03787, 0.00284], [50.0, 0.03597, 0.00502],
 [60.0, 0.03111, 0.00668], [70.0, 0.0368, 0.00393], [80.0, 0.04193, 0.00282],
 [90.0, 0.04071, 0.00527], [100.0, 0.04117, 0.00386], [110.0, 0.04032, 0.00482],
 [120.0, 0.03998, 0.00567], [130.0, 0.04339, 0.00629], [140.0, 0.04392, 0.00579],
 [150.0, 0.04514, 0.00336], [160.0, 0.05, 0.00691], [170.0, 0.04974, 0.00558],
 [180.0, 0.04985, 0.00312], [190.0, 0.05286, 0.00677], [200.0, 0.05304, 0.00819],
 [210.0, 0.05551, 0.00876], [220.0, 0.05138, 0.00392], [230.0, 0.05528, 0.00592],
 [240.0, 0.05514, 0.00651], [250.0, 0.06116, 0.00698], [260.0, 0.0597, 0.00343],
 [270.0, 0.05491, 0.00411], [280.0, 0.06042, 0.00582], [290.0, 0.05828, 0.00525],
 [300.0, 0.06087, 0.00362], [310.0, 0.06329, 0.00419], [320.0, 0.06016, 0.00347],
 [330.0, 0.05801, 0.00379], [340.0, 0.06425, 0.0044], [350.0, 0.07018, 0.00526],
 [360.0, 0.06346, 0.00519], [370.0, 0.06703, 0.00473], [380.0, 0.06777, 0.00338],
 [390.0, 0.07009, 0.00404]
])

Q4_300mDAC = np.array([
 [0.0, 0.02426, 0.00768], [10.0, 0.15948, 0.00764], [20.0, 0.05608, 0.00427],
 [30.0, 0.03679, 0.00306], [40.0, 0.02922, 0.0042], [50.0, 0.02827, 0.0027],
 [60.0, 0.02823, 0.00412], [70.0, 0.02571, 0.00227], [80.0, 0.02869, 0.00272],
 [90.0, 0.02704, 0.00218], [100.0, 0.02803, 0.00521], [110.0, 0.02984, 0.00501],
 [120.0, 0.02979, 0.00361], [130.0, 0.02955, 0.00242], [140.0, 0.02986, 0.00534],
 [150.0, 0.02748, 0.00362], [160.0, 0.02986, 0.00471], [170.0, 0.03151, 0.00394],
 [180.0, 0.03261, 0.00328], [190.0, 0.03298, 0.00444], [200.0, 0.03346,0.00341],
 [210.0, 0.03296, 0.00294], [220.0, 0.03655, 0.00429], [230.0, 0.0349, 0.00379],
 [240.0, 0.03185, 0.00226], [250.0, 0.03358, 0.00407], [260.0, 0.03534, 0.00348],
 [270.0, 0.0335, 0.00266], [280.0, 0.03489, 0.00371], [290.0, 0.03681, 0.00377],
 [300.0, 0.03674, 0.00244],[310.0, 0.03532, 0.00412], [320.0, 0.03615, 0.00321],
 [330.0, 0.03726, 0.00427], [340.0, 0.03789, 0.00507], [350.0, 0.0386, 0.00255],
 [360.0, 0.03687, 0.00431], [370.0, 0.0391, 0.00428], [380.0, 0.03889, 0.00372],
 [390.0, 0.03978, 0.00519]
])

Q4_450mDAC = np.array([
 [0.0, 0.02221, 0.00484], [10.0, 0.10304, 0.00915], [20.0, 0.06585, 0.00546],
 [30.0, 0.05722, 0.00514], [40.0, 0.05533, 0.00533], [50.0, 0.05409, 0.0042],
 [60.0, 0.0537, 0.00749], [70.0, 0.05759, 0.00662], [80.0, 0.05588, 0.00587],
 [90.0, 0.05491, 0.00582], [100.0, 0.06382, 0.00692], [110.0, 0.06643, 0.00949],
 [120.0, 0.06656, 0.00592], [130.0, 0.07058, 0.0074], [140.0, 0.06584, 0.00672],
 [150.0, 0.07191, 0.00575], [160.0, 0.07204, 0.00685], [170.0, 0.07025, 0.00743],
 [180.0, 0.0745, 0.00706], [190.0, 0.07929, 0.00485], [200.0, 0.07594, 0.00982],
 [210.0, 0.08108, 0.00706], [220.0, 0.07855, 0.0116], [230.0, 0.08127, 0.00979],
 [240.0, 0.08388, 0.00658], [250.0, 0.07817, 0.01158], [260.0, 0.08407, 0.00459],
 [270.0, 0.08281, 0.00691], [280.0, 0.08707, 0.00831], [290.0, 0.08556, 0.00665],
 [300.0, 0.08833, 0.00917], [310.0, 0.0866, 0.00777], [320.0, 0.08832, 0.00497],
 [330.0, 0.09332, 0.00959], [340.0, 0.09004, 0.00692], [350.0, 0.09167, 0.01008],
 [360.0, 0.09227, 0.00446], [370.0, 0.08834, 0.00709], [380.0, 0.09731, 0.00634],
 [390.0, 0.0951, 0.00896]
])

### J7 Bias Conversion

plt.errorbar(Q2_22mDAC[:, 0], Q2_22mDAC[:, 1], yerr=Q2_22mDAC[:, 2]/np.sqrt(10), label='Q2_22mDAC')
plt.errorbar(Q2_25mDAC[:, 0], Q2_25mDAC[:, 1], yerr=Q2_25mDAC[:, 2]/np.sqrt(10), label='Q2_25mDAC')
plt.errorbar(Q2_32mDAC[:, 0], Q2_32mDAC[:, 1], yerr=Q2_32mDAC[:, 2]/np.sqrt(10), label='Q2_32mDAC')
plt.errorbar(Q2_56mDAC[:, 0], Q2_56mDAC[:, 1], yerr=Q2_56mDAC[:, 2]/np.sqrt(10), label='Q2_56mDAC')
# plt.errorbar(Q2_300mDAC[:, 0], Q2_300mDAC[:, 1], yerr=Q2_300mDAC[:, 2]/np.sqrt(10), label='Q2_300mDAC')
# plt.errorbar(Q2_450mDAC[:, 0], Q2_450mDAC[:, 1], yerr=Q2_450mDAC[:, 2]/np.sqrt(10), label='Q2_450mDAC')
# plt.errorbar(Q4_300mDAC[:, 0], Q4_300mDAC[:, 1], yerr=Q4_300mDAC[:, 2]/np.sqrt(10), linestyle='--', label='Q4_300mDAC')
# plt.errorbar(Q4_450mDAC[:, 0], Q4_450mDAC[:, 1], yerr=Q4_450mDAC[:, 2]/np.sqrt(10), linestyle='--', label='Q4_450mDAC')

# plt.axvline(x=DAC_Al * f, color='k', linestyle='--', linewidth=4, label='JJ Al Gap')

# plt.xlabel('Radiator Josephson Frequency (GHz)')
plt.xlabel('Josephson Bias time (us)')
plt.ylabel('P1')
plt.yscale('log')
plt.grid()
plt.legend(loc=1)
# plt.xlim([0, 1500])
# plt.ylim([10, 100000])
plt.show()


