function photon_fluence = photonFluence(power, attenuation, frequency,...
    measurement_time)
% photonFluence(power, attenuation, frequency, measurement_time) calculates
% the photon fluence. POWER should be given in dBm, ATTENUATION — in dB, 
% FREQUENCY — in GHz, and MEASUREMENT_TIME — in nanosec.

h = 6.62606957e-34; % m^2 * kg / s

power = 1e-3 * 10^(power/10) * 10^(-attenuation/10); % W
frequency = 1e9 * frequency; % Hz
measurement_time = 1e-9 * measurement_time; % sec

photon_fluence = power * measurement_time / (h * frequency);

end

