import matplotlib.pyplot as plt
import numpy as np

power = np.arange(-20, -5, 1)
f_Q4_Clean = [1/0.00081588, 1/0.0008277, 1/0.00089667, 1/0.000797989, 1/0.0009013, 1/0.000706, 1/0.0007448,
              1/0.000773, 1/0.000701, 1/0.000768, 1/0.000893, 1/0.0007419, 1/0.00088785, 1/0.0008675, 1/0.000864]
f_Q4_Dirty = [1/0.0008065, 1/0.0007252, 1/0.00087836, 1/0.000856551, 1/0.0008549, 1/0.0008078, 1/0.0007774,
              1/0.0006835, 1/0.0006353, 1/0.000768, 1/0.0007711, 1/0.0006484, 1/0.000738, 1/0.0007476, 1/0.0005892]
f_Q4_Base_Clean = [1/0.00095836]
f_Q4_Base_Dirty = [1/0.00090757]
# '+', 'o', 'v', 's'
fig = plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(power, f_Q4_Clean, '-+', label='Q4 Poison Clean')
plt.plot(power, f_Q4_Dirty, '-o', label='Q4 Poison Dirty')
# plt.plot(power, f_Q4_Base_Clean, '--v', label='Q4 No Poison Clean')
# plt.plot(power, f_Q4_Base_Dirty, '--s', label='Q4 No Poison Dirty')
plt.xlabel('Relative Poison Power (dBm)')
plt.ylabel('PSD Knee Freq (Hz)')
plt.title('QP Tunneling Rate vs Poison Power')
plt.legend(loc=2)
plt.grid()
plt.show()