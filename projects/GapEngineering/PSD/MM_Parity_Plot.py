import numpy as np
import matplotlib.pyplot as plt

parity = np.array(
    [1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., 1., -1.,
     -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.,
     -1., -1., -1., -1., -1., 1., -1., 1., 1., 1., 1., 1., 1.,
     1., 1., -1., 1., -1., 1., 1., 1., 1., 1., 1., -1., 1.,
     1., -1., 1., -1., 1., 1., -1., -1., 1., -1., -1., -1., -1.,
     1., 1., -1., -1., -1., -1., -1., 1., -1., -1., -1., -1., -1.,
     -1., -1., -1., -1., -1., -1., 1., -1., -1., -1., -1., -1., -1.,
     -1., -1., 1., -1., -1., -1., -1., 1., -1.])
l = 0
r = 100
parity = parity[l:r]
t = np.linspace(0, len(parity), len(parity)) * 0.05
# x = np.linspace(0, 8*np.pi, 1000)
# s = 0.3*np.sin(x)
xticks = np.linspace(0, 5, 11)
# print(len(t))

fig = plt.figure(figsize=(8, 3))
# fig = plt.figure()
plt.plot(t, parity, 'o')
# plt.plot(x, s, 'o')
plt.xlabel('time (ms)', fontsize=20)
plt.ylabel('Parity', fontsize=20)
plt.title('Parity Trace', fontsize=24)
plt.yticks([-1, 1], fontsize=16)
plt.xticks(xticks, fontsize=16)
plt.tight_layout()
plt.show()
