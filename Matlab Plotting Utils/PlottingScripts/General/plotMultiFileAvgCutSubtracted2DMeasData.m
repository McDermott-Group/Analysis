function plotMultiFileAvgCutSubtracted2DMeasData(normalization_direction,...
    data_variable)
% plotMultiFileAvgCutSubtracted2DMeasData(NORMALIZATION_DIRECTION,
% DATA_VARIABLE) Plot a mean of several 2D data files subtracting
% the median cuts from each slice.
%   plotMultiFileAvgCutSubtracted2DMeasData(DATA_VARIABLE,
%   NORMALIZATION_DIRECTION) plots mean of several cut-substracted data
%   sets. Before the averaging, the median of all slices is computed and
%   removed from the data.

if exist('normalization_direction', 'var') &&...
        ~isempty(strfind(normalization_direction, 'x'))
    normalization_direction = 'along_x';
else
    normalization_direction = 'along_y';
end

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
                      strrep(data_variable, '_', ' '), ''' are not specified.'])
        end

        % Plot 2D data.
        if length(data{q}.rels.(data_variable)) == 2
            if strcmp(normalization_direction, 'along_y')
                dep_vals = dep_vals - (ones(size(dep_vals, 1), 1) *...
                    median(dep_vals));
            elseif strcmp(normalization_direction, 'along_x')
                dep_vals = dep_vals - (median(dep_vals, 2) *...
                    ones(1, size(dep_vals, 2)));
            end
            if q == 1
                avg_data = dep_vals;
            else
                if any(data{1}.(dep_rels{1}) ~= data{q}.(dep_rels{1}))
                    if any(flip(data{1}.(dep_rels{1})) == data{q}.(dep_rels{1}))
                        dep_vals = flip(dep_vals, 1);
                    end
                end
                if any(data{1}.(dep_rels{2}) ~= data{q}.(dep_rels{2}))
                    if any(flip(data{1}.(dep_rels{2})) == data{q}.(dep_rels{2}))
                        dep_vals = flip(dep_vals, 2);
                    end
                end
                avg_data = avg_data + dep_vals;
            end
        end
    end
    processed_data_var = ['CutSubtracted_', data_variable];
    data2plot = data{1};
    data2plot.(processed_data_var) = avg_data / length(data);
    data2plot.units.(processed_data_var) = data{1}.units.(data_variable);
    data2plot.rels.(processed_data_var) = data{1}.rels.(data_variable);
    data2plot.dep{length(data{1}.dep) + 1} = processed_data_var;
    data2plot.plotting.(processed_data_var).plot_title =...
                    {strrep(dep_vars{k}, '_', ' '),...
            [strrep(filenames{1}, '_', '\_'), ' - ',...
            strrep(filenames{end}, '_', '\_')],...
            ['[', data{1}.Timestamp, ' - ',...
             data{end}.Timestamp, ']']};
    data2plot.plotting.(processed_data_var).full_name =...
        ['_', normalization_direction];
    data2plot.plotting.(processed_data_var).extra_filename =...
        ['_', normalization_direction];

    plotDataVar(data2plot, processed_data_var);
    
    saveMeasData(data2plot, [data_variable, '_avgCutSubtr'])

end
end