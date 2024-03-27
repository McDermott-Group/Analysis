from .SFQlib import RB_ORBIT
import numpy as np


file_path = ('Z:/mcdermott-group/data/sfq/MCM_NIST/LIU/MCM13/{}/{}/MATLABData/{}')


# date = '06-24-22'
date = '2022Jun25Q1Over2'


# exp_name = 'ORBIT1_RF'
# data_type0 = 'SFQ_Drive_Frequency'

# exp_name = 'ORBIT4_Pulses_Over2'
exp_name = 'ORBIT4_Pulses'
# exp_name = 'ORBIT4_Pulses_YPhase'
# exp_name = 'ORBIT4_Pulses_YPhaseNeg'
# exp_name = 'ORBIT4_Pulses_YPhaseFine'
# data_type0 = 'SFQ_Y_Gate_Phase_Offset'
# data_type0 = 'SFQ_Pi_Over_2_Duration'
data_type0 = 'SFQ_Buffer_Between_Pulses'
# data_type0 = 'SFQ_Pulse_Duration'

# exp_name = 'ORBIT2_Mixer'
# data_type0 = 'SFQ_Pi_Amplitude'
# data_type0 = 'SFQ_I_Offset'
# data_type0 = 'SFQ_Q_Offset'

# exp_name = 'ORBIT2_MixerGain'
# data_type0 = 'SFQ_Gain_Skewness'

# exp_name = 'ORBIT2_MixerPhase'
# exp_name = 'ORBIT3_IBias_Offset'
# data_type0 = 'SFQ1_Offset'

# file_Number = [0, 1]
file_Number = np.arange(0, 8, 1)
RB_file = [file_path.format(date, exp_name, exp_name) + '_{:03d}.mat'.format(i) for i in file_Number]
RB_data = RB_ORBIT()
RB_data.add_data_from_matlab(RB_file, data_type0=data_type0)
RB_data.plot()
