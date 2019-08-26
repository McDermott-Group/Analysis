function [Aqpsd,F,es,t] = Wrapper_AnalyzeChargeTracePSD(vs,wrapping_voltage,unwrap_mult,...
                dedv,measurement_time)

%Needed additions:
% Fit for 1/f spectrum amplitude and exponent

Fs = 1/measurement_time;%samples per second
NumP = length(vs);
F = 0:Fs/NumP:(Fs/2);

[es] = unwrap_charge(vs, wrapping_voltage, dedv, unwrap_mult);
t = (0:(length(vs)-1))*(measurement_time);

%%%%% Add in noise to make QB2 look like QB1
% for i = 1:length(es)
%     sig_add = randn(1) * sqrt(3)*1e-2;
%     es(i) = es(i) + sig_add;
% end
%%%%%

[Qpsd] = crosspsd(es,es,Fs/2);
[Aqpsd] = window_averaging(Qpsd);
Aqpsd = Aqpsd(1:(NumP/2+1));

% figure(100);hold on;plot(t,es)
end

function [unwrap_es] = unwrap_charge(vs, wrap_voltage, dedv, unwrapMult)

unwrap_vs = vs;
wrap_value = 0;
for i = 1:length(vs)
    unwrap_vs(i) = vs(i) + unwrapMult*wrap_value/dedv;
    if vs(i) > wrap_voltage% + 7.5
        wrap_value = wrap_value + 1;
    end
    if vs(i) < -wrap_voltage% + 7.5
        wrap_value = wrap_value - 1;
    end
end
unwrap_es = unwrap_vs * dedv;
end

function [psd] = crosspsd(z1,z2,Fs)
%Cross PSD
NumP = length(z1);
fft_seq1 = fft(z1);
fft_seq2 = conj(fft(z2));
psd = (1/(Fs*NumP))*(fft_seq1.*fft_seq2);
psd(2:(NumP/2)) = 2*psd(2:(NumP/2));
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