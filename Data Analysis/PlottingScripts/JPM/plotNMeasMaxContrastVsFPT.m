function plotNMeasMaxContrastVsFPT(N)
%plotNMeasMaxContrastVsFPT Plot estimated maximum contrast assuming
%multiple measurments.
%   plotNMeasMaxContrastVsFPT(N) plots estimated maximum contrast assuming
%   up to N measurments.

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

[~, filename1] = fileparts(filenames{1});
[~, filename2] = fileparts(filenames{2});

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
    
    if length(dep_rels1) == 1 &&...
            ~isempty(strfind(dep_rels1{1}, 'Time'))
        t_name = dep_rels1{1};
        t = data1.(t_name);
    elseif length(dep_rels1) == 2 &&...
            ~isempty(strfind(dep_rels1{1}, 'Time'))
        t_name = dep_rels1{1};
        t = data1.(t_name);
        t = t * ones(1, length(data1.(dep_rels1{2})));
        t = t(:);
    elseif length(dep_rels1) == 2 &&...
            ~isempty(strfind(dep_rels1{2}, 'Time'))
        t_name = dep_rels1{2};
        t = data1.(t_name)';
        t = ones(length(data1.(dep_rels1{1})), 1) * t;
        t = t(:);
    else
        return
    end

    time = nan(1, N);
    max_contrast = nan(1, N);
    for n = 1:N
        contrast = (1 - dep_vals2).^n - (1 - dep_vals1).^n;
        max_contrast(n) = max(contrast(:));
        time(n) = n * t(contrast == max_contrast(n));
    end
    
    labels = cellstr(num2str((1:N)'));
    xunits = getUnits(data1, t_name);

    createFigure;
    scatter(time, max_contrast, 'filled')
    text(time + (max(time) - min(time)) / 100,...
        max_contrast + (max(max_contrast) - min(max_contrast)) / 50,...
        labels)
    grid on
    axis tight
    xlabel([strrep(t_name, '_', ' '), xunits], 'FontSize', 14);
    ylabel('Maximum Contrast', 'FontSize', 14);
    title({['Maximum Contrast for Repeated Measurements Estimated ',...
        'from Two Datasets:'],...
        [' P(bright): ', strrep(filenames{1}, '_', '\_'),...
        ' [', data1.Timestamp, ']'],...
        [' P(dark):   ', strrep(filenames{2}, '_', '\_'),...
        ' [', data2.Timestamp, ']']}, 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename1, '-', filename2,...
        '_max_contrast']));
end