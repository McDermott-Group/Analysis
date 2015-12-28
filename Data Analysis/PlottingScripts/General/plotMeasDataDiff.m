function plotMeasDataDiff
%plotMeasDataDiff   Plot the difference between two data sets.

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

for data_index = 1:length(data1.dep)
    dep_name = data1.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev')) ||...
             ~isempty(strfind(dep_name, '_Error'))
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
     
    if isempty(dep_rels1) || isempty(dep_rels2)
        continue
    end
    
    if length(dep_rels1) ~= length(dep_rels2)
        error('The selected files do not match.')
    end
    for k = 1:length(dep_rels1)
        if ~strcmp(dep_rels1{k}, dep_rels2{k}) ||...
               any(dep_rels1{k} ~= dep_rels2{k})
            error('The selected files do not match.')
        end
    end
    
    data = data1;
    diff_name = ['Difference_in_', dep_name];
    data.(diff_name) = dep_vals1 - dep_vals2;
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
         ['   ', filenames{1}, ' [', data1.Timestamp, ']'],...
         [' - ', filenames{2}, ' [', data2.Timestamp, ']']};
    data.plotting.(diff_name).plot_filename =...
        fullfile(plts_path, [base_filename1, '-', base_filename2, '_',...
        diff_name]);
     
    plotDataVar(data, diff_name);
end