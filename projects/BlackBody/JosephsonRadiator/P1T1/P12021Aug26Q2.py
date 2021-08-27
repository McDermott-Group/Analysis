### To analyze the P1 vs J2 Bias data
### src: Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorP1_2021Aug25_HighDensity\

from antennalib import P1_JSweep
import matplotlib.pyplot as plt
import numpy as np
import copy

file_path = ('Z:/mcdermott-group/data/BlackBody/Circmon/LIU/CW20180514A_Ox2/{}/{}/MATLABData/{}')

date = 'JJRadiatorP1_2021Aug26_HighDensity'

file_Number = np.arange(0, 200, 1)


experiment_name_Q2 = ('P1_J2_Q2')
file_list_Q2 = [file_path.format(date, experiment_name_Q2, experiment_name_Q2) + '_{:03d}.mat'.format(i) for i in file_Number]
P1_J2_Q2 = P1_JSweep()
P1_J2_Q2.add_data_from_matlab(file_list_Q2)

### J2 Bias offset
# J2_Offset = 19.2
J2_Bias_Q2 = copy.deepcopy(P1_J2_Q2.J2_Bias)
P1_J2_Q2 = copy.deepcopy(P1_J2_Q2.occ_1D_avg)
J2_Bias_Q2 = 1000*J2_Bias_Q2
print('J2_Bias_Q2=', J2_Bias_Q2)
print('P1_J2_Q2=', P1_J2_Q2)

plt.plot(J2_Bias_Q2, P1_J2_Q2, label='Q2')
plt.xlabel('J2 Bias (mDAC)')
plt.ylabel('P1')
# plt.yscale('log')
plt.grid()
plt.legend()
plt.show()

# plt.plot(J2_Bias_Q1*4.8, P1_J2_Q1, label='Q1')
# plt.plot(J2_Bias_Q4*4.8, P1_J2_Q4, label='Q2')
# plt.xlabel('J2 Bias (GHz)')
# plt.ylabel('P1')
# plt.grid()
# plt.legend()
# plt.show()
# ('J2_Bias_Q2=', array([70. , 70.5, 71. , 71.5, 72. , 72.5, 73. , 73.5, 74. , 74.5, 75. ,
#        75.5, 76. , 76.5, 77. , 77.5, 78. , 78.5, 79. , 79.5, 80. , 80.5,
#        81. , 81.5, 82. , 82.5, 83. , 83.5, 84. , 84.5, 85. , 85.5, 86. ,
#        86.5, 87. , 87.5, 88. , 88.5, 89. , 89.5, 90. , 90.5, 91. , 91.5,
#        92. , 92.5, 93. , 93.5, 94. , 94.5, 95. , 95.5, 96. , 96.5, 97. ,
#        97.5, 98. , 98.5, 99. , 99.5]))
# ('P1_J2_Q2=', array([0.03404607, 0.03425766, 0.03458928, 0.03405381, 0.03437726,
#        0.03440869, 0.03430476, 0.03512444, 0.03452862, 0.03493009,
#        0.03499414, 0.03537695, 0.0353496 , 0.03557982, 0.03466431,
#        0.03529033, 0.03562535, 0.03611968, 0.03565287, 0.0359003 ,
#        0.03616611, 0.03683911, 0.0371199 , 0.03710679, 0.0376104 ,
#        0.03746472, 0.03762441, 0.03795804, 0.03784616, 0.03796717,
#        0.03724041, 0.03694722, 0.03676988, 0.03708419, 0.03651043,
#        0.0371083 , 0.03673545, 0.03637769, 0.0364179 , 0.03564328,
#        0.03590894, 0.03551412, 0.03451617, 0.03570741, 0.0348187 ,
#        0.03502556, 0.03441208, 0.03516457, 0.03432769, 0.03412752,
#        0.03393104, 0.03442792, 0.03492402, 0.03357021, 0.03400459,
#        0.03412774, 0.03393825, 0.03388221, 0.03387263, 0.03409479]))

