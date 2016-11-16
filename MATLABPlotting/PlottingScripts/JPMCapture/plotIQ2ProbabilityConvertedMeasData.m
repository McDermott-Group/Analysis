function plotIQ2ProbabilityConvertedMeasData
%plotIQ2ProbabilityConvertedMeasData Convert the quadrature data to the
%measurement probability using the caalibration files.

% Select the calibration files.
[filenames, pathnames, status] = selectMeasurementDataFile(2,...
{'Select the file corresponding to the phase particle in one well...',...
 'Select the file corresponding to the phase particle in the other well...'});

if ~status
    return
end

data = loadMeasurementData;
if isempty(fields(data))
    return
end
[~, filename, ext] = fileparts(data.Filename);

filename1 = fullfile(pathnames{1}, filenames{1});
filename2 = fullfile(pathnames{2}, filenames{2});
calib1 = processMeasurementData(importMeasurementData(filename1));
calib2 = processMeasurementData(importMeasurementData(filename2));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

if isfield(calib1, 'Is') && isfield(calib1, 'Qs')
    I1 = mean(calib1.Is(:));
    Q1 = mean(calib1.Qs(:));
elseif isfield(calib1, 'I') && isfield(calib1, 'Q')
    I1 = mean(calib1.I(:));
    Q1 = mean(calib1.Q(:));
else
    error('No quadrature data in the first calibration file.')
end

if isfield(calib2, 'Is') && isfield(calib2, 'Qs')
    I2 = mean(calib2.Is(:));
    Q2 = mean(calib2.Qs(:));
elseif isfield(calib2, 'I') && isfield(calib2, 'Q')
    I2 = mean(calib2.I(:));
    Q2 = mean(calib2.Q(:));
else
    error('No quadrature data in the second calibration file.')
end

if isfield(data, 'Is') && isfield(data, 'Qs')
    I = data.Is;
    Q = data.Qs;
else
    error('The selected datafile does not contain single shot quadrature data.')
end

for k = 1:length(data.rels.Is)
    if strcmp(strrep(data.rels.Is{k}, '_', ' '), 'Repetition Index')
        rep_index = k;
    end
end

distance1 = hypot(I - I1, Q - Q1);
distance2 = hypot(I - I2, Q - Q2);

result = distance2 <= distance1;

probability = sum(result, rep_index) ./ size(result, rep_index);
prob_err = 1.95 * sqrt(probability .* (1 - probability) ./...
        size(result, rep_index));

prob_name = 'Switching_Probability';
data.(prob_name) = probability;
data.error.(prob_name) = prob_err;
data.units.(prob_name) = '';
rels = data.rels.Is;
rels(rep_index) = [];
data.rels.(prob_name) = rels;
data.dep{length(data.dep) + 1} = prob_name;
data.plotting.(prob_name).full_name = strrep(prob_name, '_', ' ');
data.plotting.(prob_name).plot_title =...
    {['Data: ', strrep(filename, '_', '\_'), ext,...
     ' [', data.Timestamp, ']'],...
     ['Calibration A: ', strrep(filenames{1}, '_', '\_'),...
     ' [', calib1.Timestamp, ']'],...
     ['Calibration B: ', strrep(filenames{2}, '_', '\_'),...
     ' [', calib2.Timestamp, ']']};
data.plotting.(prob_name).plot_filename =...
    fullfile(plts_path, [filename, '_', prob_name]);

plotDataVar(data, prob_name);
saveMeasData(data, [filename, '_', prob_name])
end