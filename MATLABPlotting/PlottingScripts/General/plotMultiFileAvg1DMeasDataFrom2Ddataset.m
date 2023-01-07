function plotMultiFileAvg1DMeasDataFrom2Ddataset(data_variable)
% plotMultiFileAvg2DMeasData(
% DATA_VARIABLE) Plots a mean of 2D data file along first variable
% that you swept in those files. You can use this file to plot average of
% 1D plots

CTC = 0;
% CutsToConsider = [];
% CutsToConsider = [1, 8];
CutsToConsider = [1, 2, 5, 6, 7, 8, 9, 10];
% CutsToConsider = [1, 2, 3, 4, 6, 7, 8, 9, 10];
% Select files to plot.
[filenames, pathnames, status] = selectMeasurementDataFile;
if ~status
    return
end

if ~iscell(filenames)
    filenames = {filenames};
end
if ~iscell(pathnames)
    pathnames = {pathnames};
end


% Read the data file, convert the variable names, and specify the units.
try
    data = cell(0);
    for k = 1:length(filenames)
        file = fullfile(pathnames{k}, filenames{k});
        data{k} = processMeasurementData(importMeasurementData(file));
    end
catch
    error('The selected files contain unmatched data.')
end

if ~exist('data_variable', 'var')
    dep_vars = selectDepDataVars(data);
elseif ~iscell(data_variable)
    dep_vars = {data_variable};
else
    dep_vars = data_variable;
end

for k = 1:length(dep_vars)
    data_variable = dep_vars{k};
   

    for q = 1:length(filenames)
        % Check that the data variable exists (compute it if necessary).
        [data{q}, data_variable] = checkDataVar(data{q}, data_variable);
        dep_vals = data{q}.(data_variable);
        dep_rels = data{q}.rels.(data_variable);
 
        if ~isfield(data{q}, data_variable)
            error(['Could not find data variable ''', data_variable, ''' in ',...
                'at least one of the files.'])
        elseif isempty(dep_rels)
                error(['Independent (sweep) variables for data variable ''',...
                      strrep(data_variable, '_', ' '),...
                      ''' are not specified.'])
        end
        if CTC == 0
            CTC
            avg_data = mean(dep_vals, 1);
        else
            CTC
            avg_data1 = dep_vals(CutsToConsider,:);
            avg_data = mean(avg_data1, 1);
        end

    processed_data_var = ['Average_', data_variable];
    data2plot = data{q};
    data2plot.(processed_data_var) = avg_data;
    data2plot.units.(processed_data_var) = data{q}.units.(data_variable);
    data2plot.rels.(processed_data_var) = {dep_rels{2}};
    data2plot.dep{length(data{q}.dep) + 1} = processed_data_var;
    data2plot.plotting.(processed_data_var).plot_title =...
            {[strrep(filenames{q}, '_', '\_'), ' - ',...
              strrep(filenames{end}, '_', '\_')],...
             ['[', data{1}.Timestamp, ' - ',...
              data{end}.Timestamp, ']']};

    plotDataVar(data2plot, processed_data_var);
[~, filename, ~] = fileparts(filenames{q});
    saveMeasData(data2plot, [filename, '_avg'])

    end
end