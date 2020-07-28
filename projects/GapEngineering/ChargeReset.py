from scipy import optimize
from numpy import pi, sin, cos
import numpy as np
import matplotlib.pyplot as plt


n = 50
x_off = round(np.random.uniform(-0.5, 0.5), 3)
noise_factor = 0.1
dispersion = 1
x_data = np.linspace(-0.5, 0.5, num=n)
y_data = 0.4 + 0.35 * cos(dispersion * pi * sin(2 * pi * (x_data-x_off)))
y_data_noise = 1.0 * y_data + noise_factor * np.random.normal(size=n)

initial_guess = np.linspace(-0.5, 0.5, num=6)
covariance = float('inf')
# print initial_guess
def test_func(x, a):
    return 0.5 + 0.5 * cos(dispersion * pi * sin(2 * pi * (x - a)))

for ig in initial_guess:
    params_curr, params_covariance_curr = optimize.curve_fit(test_func, x_data, y_data_noise, p0=[ig], method='trf')
    if params_covariance_curr[0] < covariance:
        params = params_curr
        params_covariance = params_covariance_curr
        covariance = params_covariance_curr[0]


print('x_off = ', round(x_off-0.5,3), round(x_off, 3), round(x_off+0.5, 3))
print('fitted_off = ', round(params[0], 3))
plt.figure(figsize=(6, 4))
plt.scatter(x_data, y_data_noise, label='Data_Noise')
plt.plot(x_data, test_func(x_data, params[0]), label='Fitted function')
plt.legend(loc='best')
plt.show()
