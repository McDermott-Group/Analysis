function plotCalGain(data_variable)
%plotCalGain  Plot the calibrated gain through an amplifier.
%   plotCalGain(data_variable) Plots the difference in S21 (S43) given two
%   data sets. If DATA_VARIABLE is specified, only the difference for this data
%   variable is plotted.

% Select files to compute the difference.
[filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the first data file to subtract from...',...
     'Select the second data file to be subtracted...'});
if ~status
    return
end

% Read the data files, convert the variable names, and specify the units.
file1 = fullfile(pathnames{1}, filenames{1});
file2 = fullfile(pathnames{2}, filenames{2});
data1 = processMeasurementData(importMeasurementData(file1));
data2 = processMeasurementData(importMeasurementData(file2));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

[~, base_filename1] = fileparts(filenames{1});
[~, base_filename2] = fileparts(filenames{2});

if exist('data_variable', 'var')
    % Check that the data variable exists (compute it if necessary).
    [~, data_variable] = checkDataVar(data1, data_variable);
    [~, data_variable] = checkDataVar(data2, data_variable);
    range = {data_variable};
else
    range = data1.dep;
end

for data_index = 1:length(range)
    dep_name = range{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev')) ||...
             ~isempty(strfind(dep_name, '_Error'))
        continue
    end
    dep_vals1 = data1.(dep_name);
    if ~isfield(data2, dep_name)
        disp(dep_name)
        error('The selected files do not match.')
    else
        dep_vals2 = data2.(dep_name);
    end
    dep_rels1 = data1.rels.(dep_name);
    dep_rels2 = data2.rels.(dep_name);
     
    if isempty(dep_rels1) || isempty(dep_rels2)
        continue
    end
    
%     if length(dep_rels1) ~= length(dep_rels2)
%         error('The selected files do not match.')
%     end
%     for k = 1:length(dep_rels1)
%         if ~strcmp(dep_rels1{k}, dep_rels2{k}) ||...
%                any(dep_rels1{k} ~= dep_rels2{k})
%             error('The selected files do not match.')
%         end
%     end
    
    data = data1;
    diff_name = ['Difference_in_', dep_name];

    %interpS21 = interp1(data2(data2.indep{1}),data2(dep_name),data1.(data1.indep{1}));
    interpS21 = interp1(data2.('RF_Frequency'),data2.('S43'),data1.('RF_Frequency'));
    %data.(diff_name) = dep_vals1 - dep_vals2;
    [nBias, ~] = size(data1.('S43'));
    
    % An offset introduced to deal with the calibration line not having the
    % same attenuation nor the two bias tees (each with insertion loss of
    % order ~ dB).
    offset_dB = 2 * 10 + 2;
    for bias_idx = 1:nBias
        size(data1.('S43')(bias_idx,:));
        size(interpS21(:,:)');
        data.(diff_name)(bias_idx,:) = data1.('S43')(bias_idx,:)' - interpS21(:,:) + offset_dB;
    end
        
    %data.(diff_name) = data1.('S43') - interpS21;

    if isfield(data1, 'error') && isfield(data1.error, dep_name) &&...
       isfield(data2, 'error') && isfield(data2.error, dep_name)
        data.error.(diff_name) = sqrt(data1.error.(dep_name).^2 +...
                                      data2.error.(dep_name).^2);
    end
    data.units.(diff_name) = data1.units.(dep_name);
    data.rels.(diff_name) = dep_rels1;
    data.dep{length(data.dep) + 1} = diff_name;
    data.plotting.(diff_name).full_name = strrep(diff_name, '_', ' ');
    data.plotting.(diff_name).plot_title =...
        {[strrep(diff_name, '_', ' '), getUnits(data1, dep_name),...
        ' between Two Datasets:'],...
         ['   ', strrep(filenames{1}, '_', '\_'),...
         ' [', data1.Timestamp, ']'],...
         [' - ', strrep(filenames{2}, '_', '\_'),...
         ' [', data2.Timestamp, ']']};
    data.plotting.(diff_name).plot_filename =...
        fullfile(plts_path, [base_filename1, '-', base_filename2, '_',...
        diff_name]);
     
    plotDataVar(data, diff_name);
    plotDataVar(data1, 'Voltage');
end