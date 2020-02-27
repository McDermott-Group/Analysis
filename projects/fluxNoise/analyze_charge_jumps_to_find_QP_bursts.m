% jump_indicies()'
% return

start_index = 1748;   % cuk2242kfj 2/26 741
date = '02-26-20';
charge_jump_indicies = [49,119,136,157,164,213,216,266,273,289,315,359,420,424,433,470,496,515,552,578,628,635,641,645,647,669,714,764];
trial_indicies       = [0, 0,  6,  4,  0,  0,  0,  0,  0,  2,  7,  2,  2,  0,  0,  0,  2,  6,  0,  6,  4,  0,  8,  5,  0,  0,  0,  6];
rep_indicies         = [0,0,3314,1805,0,0,0,0,0,1231,1442,1152,3363,0,0,0,3001,872,0,1455,1746,0,171,2366,0,0,0,380];

% start_index = 741;   % cuk1058ibk 2/26 741
% date = '02-26-20';
% charge_jump_indicies = [19,80,114,121,139,178,190,205,214,315,320,343,344,378,379,399,400,429,464,497,507,508,562,586,591,632,633,641,665,742,745,803,815,818,832,849,874,896,935,941];
% trial_indicies       = [8, 7, 2,  0,  0,  0,  9,  4,  8,  0,  8,  9,  0,  0,  7,  0,  7,  0,  8,  0,  0,  7,  8,  5,  7,  0,  7,  0,  3,  7,  3,  0,  0,  6,  3,  0,  0,  0,  0,  10];
% rep_indicies         = [1454,1965,1547,0,0,0,1907,1851,2933,0,753,629,0,0,2567,0,665,0,407,0,0,3907,1386,1255,3524,0,1814,0,1726,3382,844,0,0,2679,2625,0,0,0,0,];

% start_index = 30;   % cuk0242cvg 2/25-2/26 30 -524
% date = '02-25-20';
% start_index = -524;   % cuk0242cvg 2/25-2/26 30 -524
% date = '02-26-20';
% charge_jump_indicies = [8,45,46,114,274,292,294,317,381,412,472,474,552,632,647,776,784,788,800,824,832,860,875,880,886,906,1005,1058,1063,1076,1079,1157,1170];
% trial_indicies       = [];
% rep_indicies         = [];

% start_index = 1263;   % cuc1833nzp 2/18 1263
% date = '02-18-20';
% charge_jump_indicies = [5  13  18  19  25  75  121  163  177  192  229  269  315  353];
% trial_indicies       = [2  9   0   0   0   0   0    4    0    9    0    0    10   4];
% rep_indicies         = [];

% start_index = 1672;   % cud0112txh 2/18 1672
% date = '02-18-20';
% charge_jump_indicies = [3,8,26,46,53,59,94,95,110,118,119,128,142,186,204,214,220,226,287,289,303,371,403,426,461,472,514,540,542,546,590,733,780,788,792,802,803,839,865,876];
% trial_indicies       = [];
% rep_indicies         = [];

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
    
% % %     (1) Plot Ramsey Curves to find trial_indicies
%     figure(111); hold on;
%     plot(data.Single_Shot_Occupation_SB1, 'DisplayName',num2str(charge_jump_indicies(i)))
    
    if i <= length(trial_indicies)
        ti = trial_indicies(i);
        if ti > 0
            o = data.Single_Shot_Occupations_SB2;
            oi = o(ti,:);

%     %         (2) Plot Trial to find rep_indicies
%             figure(113); hold on;
%             area(movmean(oi',20), 'DisplayName',num2str(charge_jump_indicies(i)))
%             area(oi', 'DisplayName',num2str(charge_jump_indicies(i)))

%             (3) Plot autocorr, etc around that specific rep
            if i <= length(rep_indicies)
                ri = rep_indicies(i);
                oi_at_rep = oi( max(1,ri-100):min(4000,ri+100) );
                acf = acf + 1/n * noiselib.crosscorrelate( oi_at_rep, oi_at_rep, 50 );
                for j = [1:ti-1 ti+1:10]
                    acf_baseline = acf_baseline + 1/n/9 * noiselib.crosscorrelate( o(j,:), o(j,:), 50 );
                end
            end

        end
    end
end

figure; hold on;
plot(-50:50, acf_baseline)
plot(-50:50, acf)


function [jump_idxs] = jump_indicies()

data = loadMeasurementData;
dvde = data.x2e_Period/2;
es = noiselib.unwrap_voltage_to_charge(data.Offset_Voltage, dvde, 1/dvde);

jumps = es(2:end) - es(1:end-1);
jump_idxs = find( abs(jumps) > 0.1 );

figure; hold on;
plot([es,data.Offset_Voltage_R2])
jump_points = es;
jump_points(~(abs(jumps) > 0.1)) = nan;
plot(jump_points, 'o')
title(data.Filename);

end