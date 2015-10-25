function findBiasRange
%findBiasRange Find bias voltage range for biasing a JPM.

% Select a file to plot.
[filename, pathname, status] = selectMeasurementDataFile(1);
if ~status
    return
end

% Read probability data file, convert the variable names, and define
% the units.
data = processMeasurementData(importMeasurementData(fullfile(pathname, filename)));

prob = {'JPM_A_Switching_Probability',...
        'JPM_B_Switching_Probability',...
        'Switching_Probability'};
     
time = {'JPM_A_Detection_Time',...
        'JPM_B_Detection_Time',...
        'Detection_Time'};

bias = {'JPM_A_Bias_Voltage',...
        'JPM_B_Bias_Voltage',...
        'Bias_Voltage'};

name = {'JPM A',...
        'JPM B',...
        'JPM'};
    
for k = 1:length(prob)
    if isfield(data, prob{k}) && isfield(data, time{k}) &&...
        strcmp(data.rels.(prob{k}), bias{k}) && strcmp(data.rels.(time{k}), bias{k}) 
        P = data.(prob{k});
        t = data.(time{k});
        if strcmp(data.units.(bias{k}), 'V')
            BV = data.(bias{k});
        elseif strcmp(data.units(bias{k}), 'mV')
            BV = data.(bias{k}) / 1000;
        else
            error(['''', bias{k}, ''' units ''', data.units.(bias{k}), ''' are not yet supported.'])
        end
        t = t(P == 1);
        BV = BV(P == 1);
        t_threshold = .75 * max(t(t < 1000));
        BV = BV(t > t_threshold & t < 1253);
        if ~isempty(BV)
            if min(BV) < max(BV)
                disp(['Suitable bias voltage range for ', name{k}, ' is from ', num2str(min(BV)),...
                    ' V to ', num2str(max(BV)), ' V.'])
            else
                disp(['Suitable bias voltage for ', name{k}, ' is ', num2str(min(BV)), ' V.'])  
            end
        else
            disp(['Suitable bias voltage range could not be determined for ', name{k}, '.'])
        end

    end
end
