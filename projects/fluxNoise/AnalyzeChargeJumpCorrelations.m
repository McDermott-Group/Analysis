data = loadMeasurementData;

period_2e_QA = data.x2e_Period_Q3;
period_2e_QB = data.x2e_Period_Q4;
vs_QA = data.Offset_Voltage_Q3;
vs_QB = data.Offset_Voltage_Q4;

% figure; plot(t, vs_QA(1:end-1), t, vs_QB(1:end-1));
% xlabel('time [s]'); ylabel('voltage [V]'); title('Measured (Wrapped) Voltage');

t = data.Time(1:end-1) - data.Time(1);
T = (t(end)-t(1))/length(t);
unwrapped_qs_QA = noiselib.unwrap_charge(vs_QA, period_2e_QA/2, 2/period_2e_QA, 1);
unwrapped_qs_QB = noiselib.unwrap_charge(vs_QB, period_2e_QB/2, 2/period_2e_QB, 1);

figure; plot(t, [unwrapped_qs_QA(1:end-1), unwrapped_qs_QB(1:end-1)]);
xlabel('time [s]'); ylabel('charge [e]'); title('Unwrapped Charge');

delta_QA = (unwrapped_qs_QA(2:end) - unwrapped_qs_QA(1:end-1));
delta_QA = mod(0.5 + delta_QA, 1) - 0.5;
delta_QB = (unwrapped_qs_QB(2:end) - unwrapped_qs_QB(1:end-1));
delta_QB = mod(0.5 + delta_QB, 1) - 0.5;

figure; plot(t, abs(delta_QA), t, abs(delta_QB));
xlabel('time [s]'); ylabel('\Delta Q [e]'); title('Charge Jump Size');

figure; plot(t, unwrapped_qs_QA(1:end-1), t, unwrapped_qs_QB(1:end-1), t, cumsum(delta_QA), t, cumsum(delta_QB));

jumps_QA = abs(delta_QA) > 0.1;
jumps_QB = abs(delta_QB) > 0.1;
figure;
plot(t, jumps_QA, ...
     t, jumps_QB, ...
     t, jumps_QA + jumps_QB, ...
     t, jumps_QA + [[0];jumps_QB(1:end-1)])
 n_jumps_QA = sum(jumps_QA)
 n_jumps_QB = sum(jumps_QB)
 nPerT_QA = n_jumps_QA/t(end)
 nPerT_QB = n_jumps_QB/t(end)
 correlations = sum(jumps_QA + jumps_QB > 1.5)
 correlationsOffset = sum(jumps_QA + [[0];jumps_QB(1:end-1)] > 1.5)
 expectedCorrelations = (nPerT_QA*T)*(nPerT_QB*T)*t(end)/T
 
 eitherjumps = (jumps_QA | jumps_QB);
 sum( eitherjumps.*(delta_QA.*delta_QB)./sum(eitherjumps) )
 sum( eitherjumps.*abs(delta_QA.*delta_QB)./sum(eitherjumps) )
 
correlated = jumps_QA + jumps_QB > 1.5;
figure; hold on; 
plot(~correlated.*delta_QA, ~correlated.*delta_QB, '.'); 
plot(correlated.*delta_QA, correlated.*delta_QB, '.');