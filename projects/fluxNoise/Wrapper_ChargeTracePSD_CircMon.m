function [measurement_time, expons, ampls, combined_vs] = Wrapper_ChargeTracePSD_CircMon(charge_ext, qubit, start, stop, combined_vs)

ext_q1 = strcat(['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-06-10\Circ1\Q', num2str(qubit), '\General\Parameter\']);
dedv_list = [0.3844, 0.411, nan, 0.5754];
wrapping_voltage_list = [0.3844, 0.411, nan, 0.5754];
% combined_vs = {[0], [0], [0], [0]};

% Data set 1
dedv = 2/dedv_list(qubit);%e/V  %#0.3844 # 0.411 # # 0.5754
wrapping_voltage = wrapping_voltage_list(qubit)/2;%V
unwrap_mult = 1;%If we unwrap at 5 V, then we are dropping 2e.  dedv converts
%voltage to 1e, so we need a multiplier for the unwrapping.
% charge_ext = 'ckw0658lgh_parameters.hdf5'; % Q1
% charge_ext = 'cky2155mll_parameters.hdf5'; % Q1
% charge_ext = 'ckv0334jpp_parameters.hdf5'; % Q2
% charge_ext = 'ckx2109xfx_parameters.hdf5'; % Q4
% charge_ext = 'clb0821qqz_parameters.hdf5'; % Q4
charge_file = loadMeasurementData(strcat([ext_q1,charge_ext,'_parameters.hdf5']));
if (stop==0)
    time = charge_file.Time(start:end);
    vs = charge_file.Offset_Voltage(start:end);
else
    time = charge_file.Time(start:stop);
    vs = charge_file.Offset_Voltage(start:stop);
end
measurement_time = (time(end) - time(1))/length(time);
[Aqpsd1,F1,es1,t1] = Wrapper_AnalyzeChargeTracePSD(vs,wrapping_voltage,unwrap_mult,...
                dedv,measurement_time);

temp_vs = combined_vs{qubit};
combined_vs{qubit} = [temp_vs, temp_vs(end) + es1'];

%Plot time trace
figure(111);hold on
title('Combined Time traces')
plot(measurement_time*(1:length(combined_vs{qubit})), combined_vs{qubit})
xlabel('Time (s)')
ylabel('Offset charge (e)')

% Find fill breaks
T = time(2:end) - time(1:end-1);
breaks = find(T>5*measurement_time)
nBreaksForFill = (time(breaks+1)-time(breaks))/measurement_time
            
%Plot time trace for combined measurements
figure(101);hold on
title('Time traces')
plot(t1,es1)
xlabel('Time (s)')
ylabel('Offset charge (e)')
            
%Plot PSDs
figure(102);hold on
title('PSDs')
plot(F1,Aqpsd1)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Frequency (Hz)')
ylabel('S_q (e^2/Hz)')
grid on

%Find All Exponents\
expons = zeros(1,1);
ampls = zeros(1,1);
[fr,~] = fit_psd(log(F1),log(Aqpsd1));
expons(1) = fr.b;
ampls(1) = fr.a;

%Plot histogram
[delta_1] = calc_delta(combined_vs{qubit},1);
delta_1 = mod(0.5 + delta_1, 1) - 0.5;
figure(130+qubit);h=histogram(delta_1,50);
xlabel('Jump size (e)')
ylabel('N')
set(gca,'yscale','log')

end

function [delta] = calc_delta(es,stepsize)
delta = zeros(length(es)-stepsize,1);
for i = 1:length(delta)
    delta(i) = es(i+stepsize) - es(i);
end
end

function [es_jump, es_smooth] = filterJumps(es,delta)
thresh = 0.07;%0.001;
es_jump = zeros(length(es),1) + es(1);
es_smooth = zeros(length(es),1) + es(1);
for i = 1:length(delta)
    if abs(delta(i)) < thresh
        es_smooth(i+1) = es_smooth(i) + delta(i);
        es_jump(i+1) = es_jump(i);
    else
        es_jump(i+1) = es_jump(i) + delta(i);
        es_smooth(i+1) = es_smooth(i);
    end
end
end

function [psd_cut,freqs] = crosspsd(z1,z2,ts)
Fs = length(ts)/(ts(end) - ts(1));
NumP = length(z1);
freqs = 0:Fs/NumP:(Fs/2);
%Cross PSD
NumP = length(z1);
fft_seq1 = fft(z1);
fft_seq2 = conj(fft(z2));
psd = (1/(Fs*NumP))*(fft_seq1.*fft_seq2);
psd(2:(NumP/2)) = 2*psd(2:(NumP/2));
%End PSD
psd_cut = psd(1:(NumP/2+1));
end

function [apsd] = window_averaging(psd)
% apsd = psd;
%Averaging with f window
apsd = zeros(size(psd));
for i = 1:size(psd,1)
    filter_fl = max(1,round(i - i/4));
    filter_fh = min(size(psd,1),round(i + i/4));
    apsd(i) = mean(psd(filter_fl:filter_fh));
end
end

function [fitresult, gof] = fit_psd(f,psd)
[xData, yData] = prepareCurveData(f, psd');

ft = fittype( 'a + b*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-7 -1.5];

[fitresult, gof] = fit( xData, yData, ft, opts );

% figure;plot( fitresult, xData, yData );
end