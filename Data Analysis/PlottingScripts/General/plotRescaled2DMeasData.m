function data = plotRescaled2DMeasData(normalization_direction,...
    data_variable)
%plotRescaled2DMeasData(NORMALIZATION_DIRECTION, DATA_VARIABLE) Plot
%a line-by-line rescaled 2D data.
%   DATA = plotRescaled2DMeasData(DATA_VARIABLE, NORMALIZATION_DIRECTION)
%   plots line-by-line rescaled data for vairable DATA_VARIABLE along
%   NORMALIZATION_DIRECTION. NORMALIZATION_DIRECTION should be either
%   'along_x' or 'along_y'. The function returns structure DATA containing
%   the normalized data.

if exist('normalization_direction', 'var') &&...
        ~isempty(strfind(normalization_direction, 'x'))
    normalization_direction = 'along_x';
else
    normalization_direction = 'along_y';
end

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

if ~exist('data_variable', 'var')
    data_variable = selectDepDataVars(data, true);
    if isempty(data_variable)
        return
    end
    data_variable = data_variable{1};
end

% Check that the data variable exists (compute it if necessary).
[data, data_variable] = checkDataVar(data, data_variable);

dep_vals = data.(data_variable);
dep_rels = data.rels.(data_variable);
    
if isempty(dep_rels)
    error(['Independent (sweep) variables for data variable ''',...
          strrep(data_variable, '_', ' '), ''' are not specified.'])
end

% Plot 2D data.
if length(dep_rels) == 2
    if strcmp(normalization_direction, 'along_y')
        row_min = min(dep_vals, [], 2) * ones(1, size(dep_vals, 2));
        row_max = max(dep_vals, [], 2) * ones(1, size(dep_vals, 2));
        dep_vals = (dep_vals - row_min) ./ (row_max - row_min);
    elseif strcmp(normalization_direction, 'along_x')
        col_min = ones(size(dep_vals, 1), 1) * min(dep_vals);
        col_max = ones(size(dep_vals, 1), 1) * max(dep_vals);
        dep_vals = (dep_vals - col_min) ./ (col_max - col_min);   
    end
    processed_data_var = ['Rescaled_', data_variable];
    data.(processed_data_var) = dep_vals;
    data.units.(processed_data_var) = '';
    data.rels.(processed_data_var) = data.rels.(data_variable);
    data.dep{length(data.dep) + 1} = processed_data_var;
    data.plotting.(processed_data_var).full_name =...
        ['Line-by-Line Rescaled ', strrep(data_variable, '_', ' ')];
    data.plotting.(processed_data_var).extra_filename =...
        ['_', normalization_direction];

    plotDataVar(data, processed_data_var);
end