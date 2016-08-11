function photons = temperature2photons(frequency, temperature)
% temperature2photons(FREQUENCY, TEMPERATURE) compute number of photos
% correpsonding to a specific temperature. FREQUENCY should be given
% in GHz, TEMPERATURE - in K. The returned result is the number of photons.

h = 6.62606957e-34; % m^2 * kg / s
k_B = 1.38064852e-23; % m^2 s^-2 / K

frequency = 1e9 * frequency; % Hz
photons = k_B * temperature / (h * frequency);

end