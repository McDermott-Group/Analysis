import matplotlib.pyplot as plt

power = [-20, -16, -13, -10]
f_Q6_Clean = [1000/0.4855, 1000/0.4705, 1000/0.47385, 1000/0.49975]
f_Q6_Dirty = [1000/0.46729, 1000/0.45107, 1000/0.44786,  1000/0.45201]
f_Q6_Base_Clean = [1000/0.49529, 1000/0.49529, 1000/0.49529, 1000/0.49529]
f_Q6_Base_Dirty = [1000/0.44417, 1000/0.44417, 1000/0.44417, 1000/0.44417]
fig = plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(power, f_Q6_Clean, '-+', label='Q6 Poison Clean')
plt.plot(power, f_Q6_Dirty, '-o', label='Q6 Poison Dirty')
plt.plot(power, f_Q6_Base_Clean, '--v', label='Q6 No Poison Clean')
plt.plot(power, f_Q6_Base_Dirty, '--s', label='Q6 No Poison Dirty')
plt.xlabel('Relative Poison Power (dBm)')
plt.ylabel('PSD Knee Freq (Hz)')
plt.title('QP Tunneling Rate vs Poison Power')
plt.legend(loc=2)
plt.grid()
plt.show()


# format is [Clean_P1, Dirty_P1, Clean_PSD(sec), Dirty_PSD(sec)]
Neg100 = [0.01088, 0.163234, 0.00326769, 0.0030713]
Neg20 = [0.11223, 0.163234, 0.003305995, 0.0030117]
Neg19 = [0.01205445, 0.13125, 0.0033016, 0.0031228]
Neg18 = [0.0122843, 0.153764665, 0.00336, 0.0030478]
Neg17 = [0.0126, 0.1511, 0.00338, 0.00295]
Neg16 = [0.014697, 0.17147, 0.00334986, 0.00313996]
Neg15 = [0.0135014, 0.20299277, 0.00337222, 0.0030034]
Neg14 = [0.01409983, 0.135112, 0.0033582, 0.0030353]
Neg13 = [0.014199, 0.1678699, 0.0033385, 0.002923369]
Neg12 = [0.017434, 0.13180, 0.003337, 0.003108686]
Neg11 = [0.016272, 0.20571894, 0.0032894, 0.0029882]
Neg10 = [0.01881, 0.16950, 0.0031829, 0.0029656]
Neg9 = [0.02114436, 0.1204981, 0.00309277, 0.0029963]
Neg8 = [0.0227389, 0.1171387, 0.00316297, 0.0029226]
Neg7 = [0.0234465, 0.1220778, 0.00310084, 0.002934545]
Neg6 = [0.023717, 0.1345021, 0.0031196697, 0.00289809]