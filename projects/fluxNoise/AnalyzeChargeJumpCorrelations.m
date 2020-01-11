data = loadMeasurementData;

% period_2e_QA = data.x2e_Period_Q3;
% period_2e_QB = data.x2e_Period_Q4;
period_2e_QA = 0.4126;
period_2e_QB = 0.418;
vs_QA = data.Offset_Voltage_Q3;
vs_QB = data.Offset_Voltage_Q4;

% figure; plot(t, vs_QA(1:end-1), t, vs_QB(1:end-1));
% xlabel('time [s]'); ylabel('voltage [V]'); title('Measured (Wrapped) Voltage');

t = data.Time - data.Time(1);
T = t(end)/length(t);
unwrapped_qs_QA = noiselib.unwrap_charge(vs_QA, period_2e_QA/2, 2/period_2e_QA, 1);
unwrapped_qs_QB = noiselib.unwrap_charge(vs_QB, period_2e_QB/2, 2/period_2e_QB, 1);

delta_QA = (unwrapped_qs_QA(2:end) - unwrapped_qs_QA(1:end-1));
delta_QA = noiselib.alias(delta_QA, 0.5);
delta_QB = (unwrapped_qs_QB(2:end) - unwrapped_qs_QB(1:end-1));
delta_QB = noiselib.alias(delta_QB, 0.5);
unwrapped_qs_QA = [[unwrapped_qs_QA(1)]; unwrapped_qs_QA(1)+cumsum(delta_QA)];
unwrapped_qs_QB = [[unwrapped_qs_QB(1)]; unwrapped_qs_QB(1)+cumsum(delta_QB)];

figure; plot(t, [unwrapped_qs_QA(1:end-1), unwrapped_qs_QB(1:end-1)]);
xlabel('time [s]'); ylabel('charge [e]'); title('Unwrapped Charge');

figure; plot(t, abs(delta_QA), t, abs(delta_QB));
xlabel('time [s]'); ylabel('\Deltan [e]'); title('Charge Jump Size');

jumps_QA = abs(delta_QA) > 0.1;
jumps_QB = abs(delta_QB) > 0.1;
correlated = jumps_QA + jumps_QB > 1.5;
correlatedOffset = jumps_QA + [[0];jumps_QB(1:end-1)] > 1.5;

% figure;
% plot(t, jumps_QA, ...
%      t, jumps_QB, ...
%      t, jumps_QA + jumps_QB, ...
%      t, jumps_QA + [[0];jumps_QB(1:end-1)])
 
 n_jumps_QA = sum(jumps_QA)
 n_jumps_QB = sum(jumps_QB)
 nPerT_QA = n_jumps_QA/t(end)
 nPerT_QB = n_jumps_QB/t(end)
 correlations = sum(correlated)
 correlationsOffset = sum(correlatedOffset)
 expectedCorrelations = (nPerT_QA*T)*(nPerT_QB*T)*t(end)/T
 
 eitherjumps = (jumps_QA | jumps_QB);
 sum( eitherjumps.*(delta_QA.*delta_QB)./sum(eitherjumps) )
 sum( eitherjumps.*abs(delta_QA.*delta_QB)./sum(eitherjumps) )
 
figure; hold on; 
plot(~correlated.*delta_QA, ~correlated.*delta_QB, '.', 'MarkerSize', 10, 'DisplayName', 'Uncorrelated'); 
plot(correlated.*delta_QA, correlated.*delta_QB, '.', 'MarkerSize', 10, 'DisplayName', 'Correlated');
xlabel('\Deltan_A [e]'); ylabel('\Deltan_B [e]'); title('Charge Jump Size'); legend;
