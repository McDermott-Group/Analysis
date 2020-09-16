fig_path = 'Z:\mcdermott-group\users\ChrisWilen\FluxNoise\figs\';


path = 'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Q4Corr\General\08-30-20\Charge_resetting\MATLABData\';
% 23-24, 28-29, 34-35, 64-65, 71-73
data1 = loadMeasurementData([path 'Charge_resetting_024.mat']);
data2 = loadMeasurementData([path 'Charge_resetting_023.mat']);
o1 = data1.Single_Shot_Occupation_RO2_SB1;
o2 = data2.Single_Shot_Occupation_RO2_SB1;
o1 = flip(o1); o2 = flip(o2);
b = -0.5:0.01:0.49;
bp = -0.5:0.01:0.49;
break_index = 48;
guess_offset = mean(o1);
guess_vis = max(o1) - min(o1);
[f1, gof1] = fit(b',o1','offset+vis/2*cos(pi*cos(pi*(x-noff)/nper))', ...
               'Start',[-0.1,0.3,guess_offset,guess_vis] );
[f2, gof2] = fit(b',o2','offset+vis/2*cos(pi*cos(pi*(x-noff)/nper))', ...
               'Start',[0.1,0.3,guess_offset,guess_vis], ...
               'Exclude',b<b(break_index));
f1.offset = f1.offset - 0.02;
f2p = f1; f2p.noff = f2.noff;
[gof1.rsquare, gof2.rsquare]
figure; hold on;
p1 = plot(b,o1,'.');
p2 = plot(b,o2,'.');
p3 = plot(bp', f1(bp), 'Color', [0 0.4470 0.7410]); % b
p4 = plot(bp(break_index:end)', f2p(bp(break_index:end)), 'Color', [0.8500 0.3250 0.0980]); % r
p5 = plot(bp(1:break_index)', f2p(bp(1:break_index)), ':', 'Color', [0.8500 0.3250 0.0980]); % r
p3(1).LineWidth = 1.5;
p4(1).LineWidth = 1.5;
p5(1).LineWidth = 1.5;
p1(1).Marker = '.';
p2(1).Marker = '.';
p1(1).MarkerSize = 10;
p2(1).MarkerSize = 10;
xlabel('Applied Offset Voltage (V)')
ylabel('P_{01}')
ylim([0,0.9]);
box on;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [75/10, 50/10]);
set(gcf,'units','centimeters','position',[3,3,75/10,50/10]);
set(gca,'FontSize',10);
saveas(gcf, [fig_path 'ramsey_fit.pdf'])


path = 'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1\General\07-19-20\Qubit_spectroscopy\MATLABData\Qubit_spectroscopy_001.mat';
data = loadMeasurementData(path);
data.Qubit_Frequency = (data.Qubit_Frequency - 4.56415)*1000;
data.Single_Shot_Occupation_SB1 = flip(data.Single_Shot_Occupation_SB1);
b = data.Qubit_Slow_Bias_1;
plotMeasData('Single Shot Occupation SB1', data);
title('');
xlabel('Applied Offset Voltage (V)');
ylabel('f-f_{01} (MHz)');
colorbar off;
offset = 0;  vis = 4.1;  noff = 0.08;  nper = 0.2711;
% o = offset + vis/2*cos((pi*b-noff)/nper);
o = offset + vis/2*cos((pi*b-noff)/nper);
hold on; plot(b, o, b, -o)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [75/10, 50/10]);
set(gcf,'units','centimeters','position',[3,3,75/10,50/10]);
set(gca,'FontSize',10);
saveas(gcf, [fig_path 'qubit_spec.pdf'])