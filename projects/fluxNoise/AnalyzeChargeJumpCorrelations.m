data = loadMeasurementData;

%period_2e_QA = data.x2e_Period_Q3;
%period_2e_QB = data.x2e_Period_Q4;
period_2e_QA = 0.4126;
period_2e_QB = 0.418;
vs_QA = data.Offset_Voltage_Q3;
vs_QB = data.Offset_Voltage_Q4;

% figure; plot(t, vs_QA(1:end-1), t, vs_QB(1:end-1));
% xlabel('time [s]'); ylabel('voltage [V]'); title('Measured (Wrapped) Voltage');

t = data.Time - data.Time(1);
T = t(end)/length(t);
unwrapped_qs_QA = noiselib.unwrap_voltage_to_charge(vs_QA, period_2e_QA/2, 2/period_2e_QA);
unwrapped_qs_QB = noiselib.unwrap_voltage_to_charge(vs_QB, period_2e_QB/2, 2/period_2e_QB);

delta_QA = (unwrapped_qs_QA(2:end) - unwrapped_qs_QA(1:end-1));
delta_QB = (unwrapped_qs_QB(2:end) - unwrapped_qs_QB(1:end-1));

% figure; plot(t, [unwrapped_qs_QA, unwrapped_qs_QB]);
% xlabel('time [s]'); ylabel('charge [e]'); title('Unwrapped Charge');
% 
% figure; plot(t(1:end-1), abs(delta_QA), t(1:end-1), abs(delta_QB));
% xlabel('time [s]'); ylabel('\Deltan [e]'); title('Charge Jump Size');

jumps_QA = abs(delta_QA) > 0.1;
jumps_QB = abs(delta_QB) > 0.1;
correlated = jumps_QA + jumps_QB > 1.5;
correlatedOffset = jumps_QB + [jumps_QA(2:end); [0]] > 1.5;

figure; hold on;
plot(t, unwrapped_qs_QA, 'Color', [0 0.4470 0.7410], 'DisplayName', 'A')
plot(t, unwrapped_qs_QB, 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'B');
% all jumps
jump_pointsA = [jumps_QA; [0]] | [[0]; jumps_QA];
jump_pointsB = [jumps_QB; [0]] | [[0]; jumps_QB];
qs_jumpA = unwrapped_qs_QA.*jump_pointsA;
qs_jumpB = unwrapped_qs_QB.*jump_pointsB;
qs_jumpA(~jump_pointsA) = NaN;
qs_jumpB(~jump_pointsB) = NaN;
plot(t, qs_jumpA, 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
plot(t, qs_jumpB, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2);
% normal correlation
jump_points = [correlated; [0]] | [[0]; correlated];
qs_corrA = unwrapped_qs_QA.*jump_points;
qs_corrB = unwrapped_qs_QB.*jump_points;
qs_corrA(~jump_points) = NaN;
qs_corrB(~jump_points) = NaN;
plot(t, qs_corrA, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2)
plot(t, qs_corrB, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
% offset correlation
jump_points = [correlatedOffset; [0]] | [[0]; correlatedOffset];
qs_corrA = unwrapped_qs_QA.*[[0]; jump_points(1:end-1)];
qs_corrB = unwrapped_qs_QB.*jump_points;
qs_corrA(~[[0]; jump_points(1:end-1)]) = NaN;
qs_corrB(~jump_points) = NaN;
plot(t, qs_corrA, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
plot(t, qs_corrB, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);

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
 
figure(1131); hold on; 
corrOffsetSelectB = correlatedOffset; 
corrOffsetSelectA = [[0]; correlatedOffset(1:end-1)]; 
noCorrSelect = ones(size(correlated)); noCorrSelect(correlated|corrOffsetSelectA|corrOffsetSelectB) = 0;
corrA = correlated.*delta_QA; corrB = correlated.*delta_QB;
noCorrA = noCorrSelect.*delta_QA; noCorrB = noCorrSelect.*delta_QB;
corrOffsetA = correlatedOffset.*[delta_QA(2:end); [0]]; corrOffsetB = correlatedOffset.*delta_QB;
corrA(~corrA) = NaN; corrB(~corrB) = NaN;
noCorrA(~noCorrA) = NaN; noCorrB(~noCorrB) = NaN;
corrOffsetA(~corrOffsetA) = NaN; corrOffsetB(~corrOffsetB) = NaN;
plot(noCorrA, noCorrB, '.', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [0 0.4470 0.7410], 'DisplayName', 'Uncorrelated'); 
plot(corrA, corrB, '.', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'DisplayName', 'Correlated');
plot(corrOffsetA, corrOffsetB, '.', 'MarkerSize', 10, ...
    'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'DisplayName', 'Correlated (Offset)');
xlabel('\Deltan_A [e]'); ylabel('\Deltan_B [e]'); title('Charge Jump Size'); 
legend; axis square; grid on; box on;
