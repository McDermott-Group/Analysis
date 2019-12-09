% Simulate QP Tunneling curves

N = 5*8192; % length of line to calc self cpsd
n = 100; % number of files to average
charge_dispersion = 3e6; %pk-pk Hz
offset = 3e6;
idle_time = 200e-9;
cycle_time = 100e-6;
T_qp = 0.01e-3; % QP tunneling rate

Fs = 1/cycle_time;
ng = rand(n, N);
dP = rand(n, N) < 1/T_qp*cycle_time;
parity = mod(cumsum(dP,2),2);
% parity = ones(n, N); % &&& this needs to be random, changing on the ms timescale
detuning = offset + charge_dispersion/2*parity.*abs(sin(pi*ng)); % this is just using a sin instead of the actual fn for parity bands
rotation = 2*pi*detuning*idle_time;
M = rand(n, N) < 0.5+0.5*sin(rotation);

[avg_cpsd, psd_freq] = noiselib.partition_and_avg_psd(M, Fs);
window_avg_psd = noiselib.window_averaging(avg_cpsd');

%Plot
figure(111);hold on
title('1/f Averaged PSD')
plot(psd_freq,abs(window_avg_psd), 'DisplayName', [samples(1:end-1),qubit(1:end-1)])
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Frequency (Hz)')
ylabel('S_\eta (\eta^2/Hz)')
grid on

% calc rand ng
% choose odd or even parity
% calc detuning due to ng
% calc rotation in stationary frame
% project