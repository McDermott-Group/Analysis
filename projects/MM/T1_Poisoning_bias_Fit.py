from dataChest import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#Z:\mcdermott-group\data\Micromachining\2022-09-22 - DR2\MAS\MM_1\11-10-22\T1_Poisoning_vs_bias_vary_idle
path = ['Micromachining', '2023-02-28 - DR2', 'MAS','MM_SC_1.1','03-01-23', 'T1_Poisoning_vs_bias_vary_idle']
dc = dataChest(path)
contents = dc.ls()
files = contents[0]
f = files[0]
print(f)
dc.openDataset(f)
varsList = dc.getVariables()
# print(varsList)
data = dc.getData(variablesList=['Single Shot Occupation', 'Idle Time', 'Bias'])
states = data[:,0]
# print(data)
idle = data[:,1]
bias = data[:,2]

#find the dimension of one of the axis and then find the other dimension
count = 0
for i in range(len(bias)):
    if(bias[i] != bias[i-1]):
        count = count + 1

dim_bias = int(count)
dim_idle = int(len(states)/dim_bias)

#take the first couple entries of the axis that doesn't repeat immediately
idle_sep = idle[:dim_idle]

bias_sep = []
#for every 61'st element (the idle dimension) we want to add a bias to our list
for i in range(dim_bias):
    bias_sep.append(bias[dim_idle*i])

#Reshape the states vector into the relevant matrix
state_sep = np.reshape(states, (dim_bias, dim_idle))
# print(state_sep[:][0])

#fit
def func(x, a, b, c):
    return a * np.exp(-b * x) + c


# print(len(idle_sep), len(state_sep[:][0]))

gamma = []
for j in range(dim_bias):
    # plt.plot(idle_sep, state_sep[:][j])
    # plt.show
    #fit to f(t)=A*exp(-B*t)+C
    popt, pcov = curve_fit(func, idle_sep, state_sep[:][j], [0.7, 1/60000, 0.2])
    gamma.append(popt[1])

# print(gamma)
gamma = np.array(gamma)*1000

plt.plot(bias_sep, gamma)
# plt.semilogy(bias_sep, gamma)
plt.xlabel("bias [mV]")
plt.ylabel(r'$\Gamma_1$ [1\us]')
plt.show()