

start_index = 741;   % cuk1058ibk 2/26 741
date = '02-26-20';
charge_jump_indicies = [19,80,114,121,139,178,190,205,214,315,320,343,344,378,379,399,400,429,464,497,507,508,562,586,591,632,633,641,665,742,745,803,815,818,832,849,874,896,935,941];
trial_indicies       = [8, 7, 2,  0,  0,  0,  9,  4,  8,  0,  8,  9,  0,  0,  7,  0,  7,  ];
rep_indicies         = [1454];

% start_index = 30;   % cuk0242cvg 2/25-2/26 30 -524
% date = '02-25-20';
% start_index = -524;   % cuk0242cvg 2/25-2/26 30 -524
% date = '02-26-20';
% charge_jump_indicies = [];
% trial_indicies       = [];
% rep_indicies         = [];

% start_index = ??;   % cuc1833nzp
% date = '02-??-20';
% charge_jump_indicies = [];
% trial_indicies       = [];
% rep_indicies         = [];

% start_index = 1672;   % cud0112txh 2/18 1672
% date = '02-18-20';
% charge_jump_indicies = [];
% trial_indicies       = [];
% rep_indicies         = [];

% charge_jump_indicies = [53,59,142,214,226,287,303,514]%,780,788,803];  % index corresponding to a file number (offset below)
% trial_indicies       = [6,3,3,4,7,7,2,2]%,2,2,2];    % index of trial within file where jump seems to occur

pathname = ['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\' date '\Charge_resetting\MATLABData'];
% pathname = '/Volumes/smb/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/Q1Q2Corr/General/02-18-20/Charge_resetting/MATLABData';

n = length(charge_jump_indicies);

acf = zeros(1,2*50+1);
acf_baseline = zeros(1,2*50+1);

for i = 1:n
    
    filename = ['Charge_resetting_',num2str(start_index+charge_jump_indicies(i),'%03d'),'.mat'];
    % Read the data file, convert the variable names, and specify the units.
    try
        data = cell(0);
        for k = 1:n
            file = fullfile(pathname, filename);
            data = loadMeasurementData(file);
        end
    catch
        error('The selected files contain unmatched data.')
    end
    
% %     (1) Plot Ramsey Curves to find trial_indicies
    figure(111); hold on;
    plot(data.Single_Shot_Occupation_SB1, 'DisplayName',num2str(charge_jump_indicies(i)))
    
    if i <= length(trial_indicies)
        ti = trial_indicies(i);
        if ti > 0
            o = data.Single_Shot_Occupations_SB2;
            oi = o(ti,:);

    %         (2) Plot Trial to find rep_indicies
            figure(112); hold on;
            area(movmean(oi',20), 'DisplayName',num2str(charge_jump_indicies(i)))
            area(oi', 'DisplayName',num2str(charge_jump_indicies(i)))

% %             (3) Plot autocorr, etc around that specific rep
%             if i <= length(rep_indicies)
%                 ri = rep_indicies(i);
%                 oi_at_rep = oi(ri-100:ri+100);
%                 acf = acf + 1/n * noiselib.crosscorrelate( oi_at_rep, oi_at_rep, 50 );
%                 for j = [1:ti-1 ti+1:10]
%                     acf_baseline = acf_baseline + 1/n/9 * noiselib.crosscorrelate( o(j,:), o(j,:), 50 );
%                 end
%             end

        end
    end
end

figure; hold on;
plot(-50:50, acf_baseline)
plot(-50:50, acf)