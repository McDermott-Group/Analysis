function voltage = freq2volt(frequency)
% freq2volt(frequency) converge voltage to voltage correspoinding to 
% single electron charge. FREQUENCY should be given in GHz, 
% the output voltage is in volts.

h = 6.62606957e-34; %m^2 * kg / s
e =  1.6021766208e-19; % C 

frequency = 1e9 * frequency; % Hz
voltage = h * frequency / e;

end

