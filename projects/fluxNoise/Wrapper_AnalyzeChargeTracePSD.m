function [Aqpsd,F,es,t] = Wrapper_AnalyzeChargeTracePSD(vs,wrapping_voltage,...
    unwrap_mult, dedv,measurement_time)

%Needed additions:
% Fit for 1/f spectrum amplitude and exponent

Fs = 1/measurement_time;%samples per second
NumP = length(vs);
F = 0:Fs/NumP:(Fs/2);

[es] = noiselib.unwrap_charge(vs, wrapping_voltage, dedv, unwrap_mult);
t = (0:(length(vs)-1))*(measurement_time);

%%%%% Add in noise to make QB2 look like QB1
% for i = 1:length(es)
%     sig_add = randn(1) * sqrt(3)*1e-2;
%     es(i) = es(i) + sig_add;
% end
%%%%%

[Qpsd] = noiselib.crosspsd(es,es,Fs/2);
[Aqpsd] = noiselib.window_averaging(Qpsd);
Aqpsd = Aqpsd(1:(NumP/2+1));

% figure(100);hold on;plot(t,es)
end