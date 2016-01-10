function plotNMeasContrast(n)
%plotNMeasContrast(n) Plot estimated contrast for multiple measurements.
%   plotNMeasContrast(n) plot estimated contrast for multiple measurements,
%   assuming n number of measurements.

% Select files for contrast estimation.
[filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the first data file containing P[bright]...',...
     'Select the second data file containing P[dark]...'});
if ~status
    return
end

if ~exist('n', 'var')
    n = 1;
end

% Read probability data file, convert the variable names, and define
% the units.
file1 = fullfile(pathnames{1}, filenames{1});
file2 = fullfile(pathnames{2}, filenames{2});
data1 = processMeasurementData(importMeasurementData(file1));
data2 = processMeasurementData(importMeasurementData(file2));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

[~, base_filename1] = fileparts(filenames{1});
[~, base_filename2] = fileparts(filenames{2});

for data_index = 1:length(data1.dep)
    dep_name = data1.dep{data_index};
    if ~strcmp(dep_name, 'Switching_Probability') &&...
            ~strcmp(dep_name, 'JPM_A_Switching_Probability') &&...
            ~strcmp(dep_name, 'JPM_B_Switching_Probability')
        continue
    end
    dep_vals1 = data1.(dep_name);
    if ~isfield(data2, dep_name)
        error('The selected files do not match.')
    else
        dep_vals2 = data2.(dep_name);
    end
    dep_rels1 = data1.rels.(dep_name);
    dep_rels2 = data2.rels.(dep_name);
     
    if isempty(dep_rels1) || isempty (dep_rels2)
        continue
    end
    
    if length(dep_rels1) ~= length(dep_rels2)
        error('The selected files do not match.')
    end
    for k = 1:length(dep_rels1)
        if ~strcmp(dep_rels1{k}, dep_rels2{k}) ||...
               length(dep_rels1{k}) ~= length(dep_rels2{k}) ||...
               any(dep_rels1{k} ~= dep_rels2{k})
            error('The selected files do not match.')
        end
    end
    
    contrast = (1 - dep_vals2).^n - (1 - dep_vals1).^n;
    
    data = data1;
    contrast_name = 'Contrast';
    data.(contrast_name) = contrast;
    
    data.units.(contrast_name) = '';
    data.rels.(contrast_name) = dep_rels1;
    data.dep{length(data.dep) + 1} = contrast_name;
    data.plotting.(contrast_name).full_name =...
        'Contrast';
    data.plotting.(contrast_name).plot_title =...
        {['Contrast for ', num2str(n), ' Measurements Estimated ',...
        'from Two Datasets:'],...
        [' P(bright): ', strrep(filenames{1}, '_', '\_'),...
        ' [', data1.Timestamp, ']'],...
        [' P(dark):   ', strrep(filenames{2}, '_', '\_'),...
        ' [', data2.Timestamp, ']']};
    data.plotting.(contrast_name).plot_filename =...
            fullfile(plts_path, [base_filename1, '_',...
            base_filename2, '_contrast', num2str(n)]);
    
    plotDataVar(data, contrast_name)
end