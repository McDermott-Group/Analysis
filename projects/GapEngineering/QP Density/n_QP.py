from noiselib import loadmat
import numpy as np
import matplotlib.pyplot as plt

rest = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
n_QP = [1.493, 1.71, 1.718, 1.448, 1.193, 1.671, 1.254, 0.7727, 0.6728, 1.254, 0.287, 0.6307, 0.3953, 0.0787, 0.14]
    
fig = plt.figure()
plt.plot(rest, n_QP)
axes = plt.gca()
plt.xlabel('rest time (us)')
plt.ylabel('n_QP')
# plt.legend(bbox_to_anchor=(0.75, 0.58), loc=2)
plt.show()