data = loadMeasurementData;
period_2e_Q2 = 0.4137;
period_2e_Q3 = 0.2526;
t = data.Time(1:end-1) - data.Time(1);
T = (t(end)-t(1))/length(t);
delta_Q2 = (data.Offset_Voltage_Q2(2:end) - data.Offset_Voltage_Q2(1:end-1))*2/period_2e_Q2;
delta_Q3 = (data.Offset_Voltage_Q3(2:end) - data.Offset_Voltage_Q3(1:end-1))*2/period_2e_Q3;
jumps_Q2 = abs(delta_Q2) > 0.1;
jumps_Q3 = abs(delta_Q3) > 0.1;
figure;
plot(t, jumps_Q2, ...
     t, jumps_Q3, ...
     t, jumps_Q2 + jumps_Q3, ...
     t, jumps_Q2 + [[0];jumps_Q3(1:end-1)])
 n_jumps_Q2 = sum(jumps_Q2)
 n_jumps_Q3 = sum(jumps_Q3)
 nPerT_Q2 = n_jumps_Q2/t(end)
 nPerT_Q3 = n_jumps_Q3/t(end)
 correlations = sum(jumps_Q2 + jumps_Q3 > 1.5)
 correlationsOffset = sum(jumps_Q2 + [[0];jumps_Q3(1:end-1)] > 1.5)
 expectedCorrelations = nPerT_Q2*nPerT_Q3*(T^2)*t(end)/T