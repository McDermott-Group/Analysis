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
plot(bp', f1(bp), 'b')
plot(bp(break_index:end)', f2p(bp(break_index:end)), 'r')
plot(bp(1:break_index)', f2p(bp(1:break_index)), 'r:')
p1(1).LineWidth = 1;
p2(1).LineWidth = 1;
p1(1).Marker = '.';
p2(1).Marker = '.';
p1(1).MarkerSize = 6;
p2(1).MarkerSize = 6;
xlabel('Applied Offset Voltage (V)')
ylabel('Occupation')
ylim([0,0.9]);
box on;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [75/10, 50/10]);
set(gcf,'units','centimeters','position',[3,3,75/10,50/10])
set(gca,'FontSize',10)

path = 'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1\General\07-19-20\Qubit_spectroscopy\MATLABData\Qubit_spectroscopy_001.mat';
data = loadMeasurementData(path);
data.Qubit_Frequency = (data.Qubit_Frequency - 4.56415)*1000;
plotMeasData('Single Shot Occupation SB1', data);
title('');
xlabel('Applied Offset Voltage (V)');
ylabel('f-f_{01} (MHz)');
colorbar off;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [75/10, 50/10]);
set(gcf,'units','centimeters','position',[3,3,75/10,50/10])
set(gca,'FontSize',10)