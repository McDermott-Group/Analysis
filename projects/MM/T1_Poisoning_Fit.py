from dataChest import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#Z:\mcdermott-group\data\Micromachining\2022-09-22 - DR2\MAS\MM_1\11-10-22\T1_Poisoning_vary_idle
path = ['Micromachining', '2022-09-22 - DR2', 'MAS','MM_1','11-14-22', 'T1_Poisoning_vary_idle']
dc = dataChest(path)
contents = dc.ls()
files = contents[0]
f = files[3]
print(f)
dc.openDataset(f)
varsList = dc.getVariables()
data = dc.getData(variablesList=['Single Shot Occupation', 'Idle Time', 'Delay Between Bias and X'])
states = data[:,0]
# print(data)
idle = data[:,1]
delay = data[:,2]

#find the dimension of one of the axis and then find the other dimension
count = 0
for i in range(len(delay)):
    if(delay[i] != delay[i-1]):
        count = count + 1

dim_delay = int(count)
dim_idle = int(len(states)/dim_delay)

#take the first couple entries of the axis that doesn't repeat immediately
idle_sep = idle[:dim_idle]

delay_sep = []
#for every 61'st element (the idle dimension) we want to add a delay to our list
for i in range(dim_delay):
    delay_sep.append(delay[dim_idle*i])

#Reshape the states vector into the relevant matrix
state_sep = np.reshape(states, (dim_delay, dim_idle))
# print(state_sep[:][0])

#fit
def func(x, a, b, c):
    return a * np.exp(-b * x) + c


# print(len(idle_sep), len(state_sep[:][0]))

gamma = []
for j in range(dim_delay):
    # plt.plot(idle_sep, state_sep[:][j])
    # plt.show
    #fit to f(t)=A*exp(-B*t)+C
    popt, pcov = curve_fit(func, idle_sep, state_sep[:][j], [0.7, 1/20000, 0.2])
    gamma.append(popt[1])

gamma = np.array(gamma)*1000
delay_sep = np.array(delay_sep)/1000

plt.plot(delay_sep, gamma)
plt.xlabel("delay [us]")
plt.ylabel(r'$\Gamma_1$ [1/us]')
plt.show()