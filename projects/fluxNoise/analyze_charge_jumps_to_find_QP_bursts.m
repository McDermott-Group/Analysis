charge_jump_indicies = [53,59,142,214,226,287,303,514]%,780,788,803];  % index corresponding to a file number (offset below)
trial_indicies       = [6,3,3,4,7,7,2,2]%,2,2,2];    % index of trial within file where jump seems to occur

% subtract_index = 1262; % cuc1833nzp
subtract_index = 1672; % cud0112txh
% pathname = 'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\02-18-20\Charge_resetting\MATLABData';
pathname = '/Volumes/smb/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/Q1Q2Corr/General/02-18-20/Charge_resetting/MATLABData';

n = length(charge_jump_indicies);

acf = zeros(1,2*500+1);
acf_baseline = zeros(1,2*500+1);

for i = 1:n
    
    filename = ['Charge_resetting_',num2str(subtract_index+charge_jump_indicies(i),'%03d'),'.mat'];
    % Read the data file, convert the variable names, and specify the units.
    try
        data = cell(0);
        for k = 1:n
            file = fullfile(pathname, filename);
            data = processMeasurementData(importMeasurementData(file));
        end
    catch
        error('The selected files contain unmatched data.')
    end
    
%     figure(111); hold on;
%     plot(data.Single_Shot_Occupation_SB1)
    
    o = data.Single_Shot_Occupations_SB2';
    ti = trial_indicies(i);
    acf = acf + 1/n * noiselib.crosscorrelate( o(:,ti), o(:,ti), 500 );
    for j = [1:ti-1 ti+1:10]
        acf_baseline = acf_baseline + 1/n/9 * noiselib.crosscorrelate( o(:,j), o(:,j), 500 );
    end
    
end

figure; hold on;
plot(-500:500, acf_baseline)
plot(-500:500, acf)