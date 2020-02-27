% [filenames, pathnames, status] = selectMeasurementDataFile;
% if ~status
%     return
% end
% 
% if ~iscell(filenames)
%     filenames = {filenames};
% end
% if ~iscell(pathnames)
%     pathnames = {pathnames};
% end

% data = loadMeasurementData;
% dvde = data.x2e_Period/2;
% es = noiselib.unwrap_voltage_to_charge(data.Offset_Voltage, dvde, 1/dvde);
% figure; plot([es,data.Offset_Voltage_R2])


index = 289;
start_index = 1748;   % cuk2242kfj 2/26 1748
% start_index = 741;    % cuk1058ibk 2/26 741
% start_index = -524;   % cuk0242cvg 2/25-2/26 30 -524
% start_index = 1263;   % cuc1833nzp
% start_index = 1672;   % cud0112txh

% set(gcf, 'Position', get(0, 'Screensize'));

pathname = 'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\02-26-20\Charge_resetting\MATLABData';
filenames = {};
for i=1:6
    filenames{i} = ['Charge_resetting_',num2str(start_index+index-3+i,'%03d'),'.mat'];
end

n = length(filenames);

% Read the data file, convert the variable names, and specify the units.
try
    data = cell(0);
    for k = 1:n
        file = fullfile(pathname, filenames{k});
        data{k} = processMeasurementData(importMeasurementData(file));
    end
catch
    error('The selected files contain unmatched data.')
end

mino = 1;
maxo = 0;
% figure; hold on;
for i = 1:length(filenames)
    mino = min([mino, data{i}.Single_Shot_Occupation_SB2]);
    maxo = max([maxo, data{i}.Single_Shot_Occupation_SB2]);
%     plot(data{i}.Qubit_Slow_Bias_2, data{i}.Single_Shot_Occupation_SB1)
end


figure;
set(gcf, 'Position', get(0, 'Screensize'));
for i = 1:length(filenames)
    j = floor((i-1)/ceil(n/2));
    subplot(4,ceil(n/2),j*ceil(n/2)+i);
    imagesc(movmean(data{i}.Single_Shot_Occupations_SB2',20)');
    caxis([0 0.4]);
%     colorbar
    title(num2str(index-3+i));
    subplot(4,ceil(n/2),(j+1)*ceil(n/2)+i);
    plot(data{i}.Single_Shot_Occupation_SB2);
    title(num2str(index-3+i));
    ylim([mino maxo]);
    yyaxis right;
    plot(data{i}.Single_Shot_Occupation_SB1);
    yyaxis left;
end


% figure(1122); hold on;
% imagesc(autocorrelate(data{3}.Single_Shot_Occupations_SB2))


function [acf] = autocorrelate(x)

acf = zeros(201,1);
for j = 1:201
    i = j-1;
    acf(j) = sum(x(1:end-i).*x(1+i:end))/(length(x)-i);
end

end