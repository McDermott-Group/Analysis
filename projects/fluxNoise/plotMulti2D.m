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

index = 142;
% subtract_index = 1262; % cuc1833nzp
subtract_index = 1672; % cud0112txh
pathname = 'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\02-18-20\Charge_resetting\MATLABData';
pathnames = {};
filenames = {};
for i=1:6
    pathnames{i} = pathname;
    filenames{i} = ['Charge_resetting_',num2str(subtract_index+index-3+i,'%03d'),'.mat'];
end

n = length(filenames);

% Read the data file, convert the variable names, and specify the units.
try
    data = cell(0);
    for k = 1:n
        file = fullfile(pathnames{k}, filenames{k});
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



% x=40;
figure;
for i = 1:length(filenames)
    j = floor((i-1)/ceil(n/2));
%     figure(111+x);
%     subplot(2,ceil(n/2),i);
    subplot(4,ceil(n/2),j*ceil(n/2)+i);
    imagesc(movmean(data{i}.Single_Shot_Occupations_SB2',20));
    caxis([mino maxo]);
    title(num2str(index-3+i));
%     figure(112+x);
%     subplot(2,ceil(n/2),i);
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