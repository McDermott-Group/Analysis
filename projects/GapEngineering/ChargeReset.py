# First to just follow an example
from scipy import optimize
from numpy import pi, sin, cos
import numpy as np
import matplotlib.pyplot as plt




# np.random.seed(0)
n = 1000
x_off = 0.2
noise_factor = 0.0
x_data = np.linspace(-0.5, 0.5, num=n)
y_data = 0.5 + 0.5 * cos(3.5 * pi * sin(2 * pi * (x_data-x_off)))
y_data_noise = 1.0 * y_data + noise_factor * np.random.normal(size=n)

def test_func(x, a):
    return 0.5 + 0.5 * cos(3.5 * pi * sin(2 * pi * (x - a)))

params, params_covariance = optimize.curve_fit(test_func, x_data, y_data_noise, p0=[0.1], method='trf')

plt.figure(figsize=(6, 4))
# plt.scatter(x_data, y_data, label='Data')
plt.scatter(x_data, y_data_noise, label='Data_Noise')
plt.plot(x_data, test_func(x_data, params[0]), label='Fitted function')
plt.legend(loc='best')
print(params)
plt.show()
