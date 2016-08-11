function voltage = freq2volt(frequency)
% freq2volt(FREQUENCY) convert frequency to voltage correspoinding to 
% charge of a single electron. FREQUENCY should be given in GHz, 
% the returned result is voltage given in volts.

h = 6.62606957e-34; % m^2 * kg / s
e =  1.6021766208e-19; % C 

frequency = 1e9 * frequency; % Hz
voltage = h * frequency / e;

end